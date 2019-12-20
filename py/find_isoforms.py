#!/usr/bin/env python3
import os
import argparse
from gurobipy import *

MIN_READS = 5
MAX_ROUND = 20

def parse_args():
    def str2bool(v):
        if isinstance(v, bool):
           return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-d",
                        "--data-matrix",
                        type=str,
                        required=True,
                        help="Path to DATA file of reads by exons matrix")
    parser.add_argument("-k",
                        "--isoform-count",
                        type=int,
                        default=3,
                        help="Number of isoform clusters")
    parser.add_argument("-et",
                        "--exons-tsv",
                        type=str,
                        default=None,
                        help="Path to exons TSV file")
    parser.add_argument("-ug",
                        "--unaligned-gaps",
                        type=str,
                        default=None,
                        help="Path to unligned gap tuples per read")
    parser.add_argument("-uo",
                        "--unaligned-offset",
                        type=int,
                        default=20,
                        help="Slack +- value for exons and the unaligned gaps")
    parser.add_argument("-e",
                        "--epsilon",
                        type=float,
                        default=0.2,
                        help="Epsilon value for how much can unaligned gaps can cover")
    parser.add_argument("-irp",
                        "--incomp-read-pairs",
                        type=str,
                        default='',
                        help="Path to file of read pairs that can not originate from the same isoform")
    parser.add_argument("-garb",
                        "--garbage-isoform",
                        type=str,
                        nargs='?',
                        const=True,
                        default=False,
                        help="Include in the model a garbage isoform collecting unassigned reads"),
    parser.add_argument("-names",
                        "--valid-read-names",
                        type=str,
                        required=True,
                        help="Valid read ids to do ilp"),
    parser.add_argument("-ptinfo",
                        "--poly-tail-status",
                        type=str,
                        required=True,
                        help="A file which tells if polyA tail is cut or no, first col is read name second cold is 1 if cut 0 is uncut"),
    parser.add_argument("-rg",
                        "--recycle-garbage",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help="Recycle garbage reads. To have effect -garb needs to be True."),
    parser.add_argument("-oi",
                        "--order-isoforms",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help="Force to label isoforms and order them increasingly")
    parser.add_argument("-to",
                        "--timeout",
                        type=int,
                        default=15,
                        help="Gurobi time-out in minutes. Default: 15")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=8,
                        help="Number of threads to use")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return(args)

def read_matrix(data_matrix, M):
    result = dict()
    for rid,line in enumerate(open(data_matrix)):
        line = line.rstrip()
        row  = [int(x) for x in line]
        assert len(row)==M, 'Each row in data matrix must be equal to M ({})'.format(M)
        if 1 in row:
            result[rid] = row
    return(result)

def read_incomp_read_pairs(pairs_file, rids):
    # Assumption: read pairs are represented by pairs of integers indicating the index in the data matrix of the incompatible reads
    result = list()
    for line in open(pairs_file):
        line = line.rstrip().split()
        assert len(line)==2, 'Incompatible read pairs are pairs of reads'
        r1,r2  = [int(x) for x in line]
        assert 0<=r1<len(rids) and 0<=r2<len(rids) and r1!=r1, 'Incompatible read pairs are pairs of different read indices'
        if not r1 in rids or not r2 in rids:
            continue
        result.append(row)
    return(result)

def read_exon_lengths(exons_file):
    result = [int(r.split('\t')[1])-int(r.split('\t')[0]) for r in open(exons_file).readlines()]
    for e in result:
        assert e>0, 'All exons must be positive'
    return(result)

def read_unaligned_gaps(gaps_file,nb_exons,rids):
    result = dict()
    for rid,line in enumerate(open(gaps_file)):
        if not rid in rids:
            continue
        line = line.rstrip()
        gaps = list()
        if len(line) > 0:
            for gap in line.split('\t'):
                gap = [int(x) for x in gap.split('-')]
                assert len(gap)==3, 'Each gap must be a 3-tuple'
                assert 0 <= gap[0] < nb_exons, 'Exon id ({eid}) of read id ({rid}) are between 0 and {M}'.format(eid=gap[0],rid=rid, M=nb_exons)
                assert 0 <= gap[1] < nb_exons, 'Exon id ({eid}) of read id ({rid}) are between 0 and {M}'.format(eid=gap[1],rid=rid, M=nb_exons)
                assert gap[2]>=0 , 'Unaligned length ({}) must be non-negative'.format(gap[2])
                gaps.append(gap)
        result[rid] = gaps
    return(result)

def find_segment_read(M,i):
    min_i = -1
    for j in range(len(M[i])):
        if min_i==-1 and M[i][j]==1:
            min_i = j
        if M[i][j]==1:
            max_i = j
    return((min_i,max_i))

def garbage_cost(C):
    # Sum of 1, i.e. exons that could be corrected, minus 0.5 to make sure that assigning anyway
    # to an isoform with maximum correction score is not at least as good as assigning to the garbage isoform
    return(max(sum(C.values())-0.5,1))
# TO DO: think about better ways to define the cost to assign to the garbagte isoform


def run_ILP(MATRIX, RIDS, GAP_L, EXON_L, K, polyA_trimmed, EPSILON, OFFSET, INCOMP_RIDS, garbage_isoform, order_isoforms, timeout, threads, out_prefix):
    os.makedirs(out_prefix[:out_prefix.rfind('/')], exist_ok=True)
    log_file = open('{}.log'.format(out_prefix),'w+')

    # Variables directly based on the input ------------------------------------------------------
    # I[i,j] = 1 if MATRIX[i][j]==1 and 0 if MATRIX[i][j]==0 or 2
    # C[i,j] = 1 if exon j is between the first and last exons (inclusively) covered by read i and is not in read i but can be turned into a 1
    M = len(EXON_L)
    MAX_ISOFORM_LG = sum(EXON_L)
    I = dict()
    C = dict()
    UB_NB_CORRECTIONS = 0 # Upper bound on the number of 0s that can be corrected into 1
    for i in RIDS:
        I[i] = dict()
        for j in range(0,M):
            I[i][j] = MATRIX[i][j]%2
        C[i] = dict()
        (min_i,max_i) = find_segment_read(I,i)
        if polyA_trimmed[i]:
            max_i = M-1
        for j in range(0,M):
            if min_i <= j <= max_i and MATRIX[i][j]==0:
                C[i][j]   = 1
                UB_NB_CORRECTIONS += 1
            else:
                C[i][j]   = 0
    log_file.write('# Maximum number of 0\'s that can be corrected into 1: {X}\n'.format(X=UB_NB_CORRECTIONS))

    # List of incompatible pairs of reads
    INCOMP_READ_PAIRS_AUX1 = list()
    # Read i1 and i2 have a 1 and a 2 in a given position j
    for idx,i1 in enumerate(RIDS):
        for i2 in RIDS[idx+1:]:
            incomp = False
            for j in range(M):
                if MATRIX[i1][j]*MATRIX[i2][j]==2:
                    incomp=True
                    break
            if incomp:
                INCOMP_READ_PAIRS_AUX1.append((i1,i2))
    log_file.write('# Number of incompatible read pairs due to a (1,2) pattern: {X}\n'.format(X=INCOMP_READ_PAIRS_AUX1))

    # Reading the list of provided incompatible pairs
    INCOMP_READ_PAIRS_AUX2 = list()
    for (i1,i2) in INCOMP_RIDS:
        if i1>i2:
            i1,i2=i2,i1
        INCOMP_READ_PAIRS_AUX2.append((i1,i2))
    log_file.write('# Number of incompatible read pairs provided as input: {X}\n'.format(X=INCOMP_READ_PAIRS_AUX2))
    # Fusing both lists
    INCOMP_READ_PAIRS = list(set(INCOMP_READ_PAIRS_AUX1+INCOMP_READ_PAIRS_AUX2))
    log_file.write('# Total number of incompatible read pairs: {X}\n'.format(X=INCOMP_READ_PAIRS))

    # Assigning a cost to the assignment of reads to the garbage isoform if any
    if garbage_isoform:
        GARBAGE_COST = {}
        ISOFORM_INDEX_START = 1 # The garbage isoform is isoform 0
        for i in RIDS:
            GARBAGE_COST[i] = garbage_cost(C[i])
    else:
        ISOFORM_INDEX_START = 0

    # ILP model ------------------------------------------------------
    ILP_ISOFORMS = Model('isoforms_v4_31072019')
    ILP_ISOFORMS.setParam(GRB.Param.Threads, threads)

    # Decision variables
    # R2I[i,k] = 1 if read i assigned to isoform j
    R2I    = {}
    R2I_C1 = {} # Constraint enforcing that each read is assigned to exactly one isoform
    # Could be relaxed if we want to account for possibly wrong reads, i.e. reads not from the given gene
    for i in RIDS:
        R2I[i] = {}
        for k in range(0,K):
            R2I[i][k] = ILP_ISOFORMS.addVar(
                vtype=GRB.BINARY,
                name='R2I[{i}][{k}]'.format(i=i,k=k)
            )
        R2I_C1[i] = ILP_ISOFORMS.addLConstr(
            lhs   = quicksum(R2I[i][k] for k in range(0,K)),
            sense = GRB.EQUAL,
            rhs   = 1,
            name  = 'R2I_C1[{i}]'.format(i=i)
        )

    # Implied variable: canonical exons presence in isoforms
    # E2I[j,k]    = 1 if canonical exon j is in isoform k
    # E2IR[j,k,i] = 1 if read i assigned to isoform k AND exon j covered by read i
    # Auxiliary variable
    # E2IR[j,k,i] = R2I[i,k] AND I[i,j]
    # E2I[j,k]    = max over  all reads i of E2IR[j,k,i]
    E2I     = {}
    E2I_C1  = {}
    E2IR    = {}
    E2IR_C1 = {}
    for j in range(0,M):
        E2I[j]     = {}
        E2I_C1[j]  = {}
        E2IR[j]    = {}
        E2IR_C1[j] = {}
        for k in range(ISOFORM_INDEX_START,K): # We do not assign exons for the garbage isoform if any
            E2I[j][k]     = ILP_ISOFORMS.addVar(
                vtype = GRB.BINARY,
                name  = 'E2I[{j}][{k}]'.format(j=j,k=k)
            )
            E2IR[j][k]    = {}
            E2IR_C1[j][k] = {}
            for  i in RIDS:
                E2IR[j][k][i]    = ILP_ISOFORMS.addVar(
                    vtype = GRB.BINARY,
                    name  = 'E2IR[{j}][{k}][{i}]'.format(j =j,k =k,i =i)
                )
                E2IR_C1[j][k][i] = ILP_ISOFORMS.addLConstr(
                    lhs   = E2IR[j][k][i],
                    sense = GRB.EQUAL,
                    rhs   = R2I[i][k]*I[i][j],
                    name  = 'E2IR_C1[{j}][{k}][{i}]'.format(j=j,k=k,i=i)
                )
            E2I_C1[j][k] = ILP_ISOFORMS.addGenConstrMax(
                resvar   = E2I[j][k],
                vars     = [E2IR[j][k][i] for i in RIDS],
                constant = 0.0,
                name     = 'E2I_C1[{j}][{k}]'.format(j=j,k=k)
            )
        if garbage_isoform:
            # No exon is assignd to the garbage isoform
            E2I[j][0]    = ILP_ISOFORMS.addVar(
                vtype = GRB.BINARY,
                name  = 'E2I[{j}][{k}]'.format(j=j,k=0)
            )
            E2I_C1[j][0] = ILP_ISOFORMS.addLConstr(
                lhs       = E2I[j][0],
                sense     = GRB.EQUAL,
                rhs       = 0,
                name      = 'E2I_C1[{j}][{k}]'.format(j=j,k=0)
            )

    # Adding constraints for unaligned gaps
    # If read i is assigned to isoform k, and GAP_L[i] contains (j1,j2,l), and
    # the sum of the lengths of exons in isoform k between exons j1 and j2 is L
    # the (1-EPSILON)L <= l <= (1+EPSILON)L
    GAPI    = {} # GAPI[(j1,j2,k)] = sum of the length of the exons between exons j1 and j2 (exclusively) in isoform k
    GAPI_C1 = {} # Constraint fixing the value of GAPI
    GAPR_C1 = {} # Constraint ensuring that the unaligned gap is not too short for every isoform and gap
    GAPR_C2 = {} # Constraint ensuring that the unaligned gap is not too long for every isoform and gap
    for i in RIDS:
        for (j1,j2,l) in GAP_L[i]:
            for k in range(ISOFORM_INDEX_START,K): # No such constraint on the garbage isoform if any
                if not (j1,j2,k) in GAPI:
                    GAPI[(j1,j2,k)]      = ILP_ISOFORMS.addVar(
                        vtype = GRB.INTEGER,
                        name  = 'GAPI[({j1},{j2},{k})]'.format(j1=j1,j2=j2,k=k)
                    )
                    GAPI_C1[(j1,j2,k)]   = ILP_ISOFORMS.addLConstr(
                        lhs   = GAPI[(j1,j2,k)],
                        sense = GRB.EQUAL,
                        rhs   = quicksum(E2I[j][k]*EXON_L[j] for j in range(j1+1,j2)),
                        name  = 'GAPI_C1[({j1},{j2},{k})]'.format(j1=j1,j2=j2,k=k)
                    )
                GAPR_C1[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr(
                    lhs   = (1.0-EPSILON)*GAPI[(j1,j2,k)]-OFFSET-((1-R2I[i][k])*MAX_ISOFORM_LG),
                    sense = GRB.LESS_EQUAL,
                    rhs   = l,
                    name  = 'GAPR_C1[({i},{j1},{j2},{k})]'.format(i=i,j1=j1,j2=j2,k=k)
                )
                GAPR_C2[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr(
                    lhs   = (1.0+EPSILON)*GAPI[(j1,j2,k)]+OFFSET+((1-R2I[i][k])*MAX_ISOFORM_LG),
                    sense = GRB.GREATER_EQUAL,
                    rhs   = l,
                    name  = 'GAPR_C2[({i},{j1},{j2},{k})]'.format(i=i,j1=j1,j2=j2,k=k)
                )
    # Adding constraints for incompatible read pairs
    INCOMP_READ_PAIRS_C1 = {}
    for (i1,i2) in INCOMP_READ_PAIRS:
        for k in range(ISOFORM_INDEX_START,K): # Again, no such constraint on the garbage isoform if any
            INCOMP_READ_PAIRS_C1[(i1,i2,k)] = ILP_ISOFORMS.addLConstr(
                lhs   = R2I[i1][k]+R2I[i2][k],
                sense = GRB.LESS_EQUAL,
                rhs   = 1,
                name  = 'INCOMP_READ_PAIRS_C1[({i1},{i2},{k})]'.format(i1 =i1,i2 =i2,k =k)
            )

    # [OPTIONAL] Labeling non-garbage isoforms by their exon content and forcing them to occur in increasing label order
    if order_isoforms==1:
        LABEL_I    = {}
        LABEL_I_C1 = {}
        LABEL_I_C2 = {}
        for k in range(ISOFORM_INDEX_START,K):
            LABEL_I[k]    = ILP_ISOFORMS.addVar(
                vtype = GRB.INTEGER,
                name  = 'LABEL_I[{k}]'.format(k=k)
            )
            LABEL_I_C1[k] = ILP_ISOFORMS.addLConstr(
                lhs   = LABEL_I[k],
                sense = GRB.EQUAL,
                rhs   = quicksum(E2I[j][k]*(2**j) for j in range(0,M)),
                name  = 'LABEL_I_C1[{k}]'.format(k =k)
            )
            if k > ISOFORM_INDEX_START:
                LABEL_I_C2[k] = ILP_ISOFORMS.addLConstr(
                    lhs   = LABEL_I[k],
                    sense = GRB.LESS_EQUAL,
                    rhs   = LABEL_I[k-1]-0.1,
                    name  = 'LABEL_I_C2[{k}]'.format(k=k)
                )

    # Objective function
    # OBJ[i][j][k] = 1 if read i assigned to isoform k and exon j in isoform k, 0 otherwise (k non-garbage isoform)
    OBJ    = {}
    OBJ_C1 = {}
    OBJ_SUM = LinExpr(0.0)
    for i in RIDS:
        OBJ[i]    = {}
        OBJ_C1[i] = {}
        for j in range(0,M):
            if C[i][j]==1: # 1 if exon j not in read i but can be added to it
                OBJ[i][j]    = {}
                OBJ_C1[i][j] = {}
                for k in range(ISOFORM_INDEX_START,K):
                    OBJ[i][j][k]    = ILP_ISOFORMS.addVar(
                        vtype = GRB.BINARY,
                        name  = 'OBJ[{i}][{j}][{k}]'.format(i=i,j=j,k=k)
                    )
                    OBJ_C1[i][j][k] = ILP_ISOFORMS.addGenConstrAnd(
                        resvar = OBJ[i][j][k],
                        vars   = [R2I[i][k],E2I[j][k]],
                        name   = 'OBJ_C1[{i}][{j}][{k}]'.format(i=i,j=j,k=k)
                    )
                    OBJ_SUM.addTerms(
                        coeffs = 1.0,
                        vars   = OBJ[i][j][k]
                    )
    # We add the chosen cost for each isoform assigned to the garbage isoform if any
    if garbage_isoform:
        for i in RIDS:
            OBJ_SUM.addTerms(
                coeffs = 1.0,#*GARBAGE_COST[i],
                vars   = R2I[i][0]
            )

    ILP_ISOFORMS.setObjective(
        expr  = OBJ_SUM,
        sense = GRB.MINIMIZE
    )

    ILP_ISOFORMS.write('{}.lp'.format(out_prefix))

    # Optimization
    #ILP_ISOFORMS.Params.PoolSearchMode=2
    #ILP_ISOFORMS.Params.PoolSolutions=5
    ILP_ISOFORMS.setParam('LogFile', '{}.glog'.format(out_prefix))
    ILP_ISOFORMS.setParam('TimeLimit', timeout*60)
    ILP_ISOFORMS.optimize()

    ILP_ISOFORMS_STATUS = ILP_ISOFORMS.Status

    iid_to_isoform = dict()
    iid_to_rids = dict()
    rid_to_corrections = dict()
    rid_to_iid = dict()
    print('STATUS: {}'.format(ILP_ISOFORMS.Status))
    if ILP_ISOFORMS_STATUS == GRB.Status.INF_OR_UNBD or \
           ILP_ISOFORMS_STATUS == GRB.Status.INFEASIBLE or \
           ILP_ISOFORMS_STATUS == GRB.Status.UNBOUNDED:
        status = 'NO_SOLUTION'
        log_file.write(' The model can not be solved because it is infeasible or unbounded\n')
        ILP_ISOFORMS.computeIIS()
        ILP_ISOFORMS.write('{}.ilp'.format(out_prefix))
    elif ILP_ISOFORMS_STATUS == GRB.Status.OPTIMAL or ILP_ISOFORMS_STATUS == GRB.Status.SUBOPTIMAL or ILP_ISOFORMS_STATUS == GRB.Status.TIME_LIMIT:
        if ILP_ISOFORMS_STATUS == GRB.Status.SUBOPTIMAL:
            status = 'SUBOPTIMAL'
            log_file.write('The model was stopped with status {}'.format(ILP_ISOFORMS_STATUS))
        elif ILP_ISOFORMS_STATUS == GRB.Status.TIME_LIMIT:
            status = 'TIME_LIMIT'
            log_file.write('The model was stopped with status {}'.format(ILP_ISOFORMS_STATUS))
        else:
            status = 'OPTIMAL'
        # Writing the optimal solution to disk
        solution_file = open('{}.sol'.format(out_prefix), 'w+')
        for v in ILP_ISOFORMS.getVars():
            solution_file.write('{}\t{}\n'.format(v.VarName, v.X))
        solution_file.close()
        # Isoform id to isoform structure
        if garbage_isoform:
            iid_to_isoform[0] = [0 for _ in range(0,M)]
        for k in range(ISOFORM_INDEX_START,K):
            iid_to_isoform[k] = [int(E2I[j][k].getAttr(GRB.Attr.X)>0.9) for j in range(0,M)]
        # Isoform id to read ids set
        for i in RIDS:
            rid_to_iid[i] = set()
            for k in range(0,K):
                if R2I[i][k].getAttr(GRB.Attr.X) > 0.9:
                    rid_to_iid[i].add(k)
            assert len(rid_to_iid[i])==1, 'Read {} has been assigned to {} isoforms!'.format(i,rid_to_iid[i])
            rid_to_iid[i] = next(iter(rid_to_iid[i]))
        for k in range(0,K):
            iid_to_rids[k] = set()
        for i in RIDS:
            k = rid_to_iid[i]
            iid_to_rids[k].add(i)
        # Read id to its exon corrections
        if garbage_isoform:
            for i in iid_to_rids[0]:
                rid_to_corrections[i] = [str(MATRIX[i][j]) for j in range(M)]
                for j in range(M):
                    if C[i][j] == 1:
                        rid_to_corrections[i][j] = 'X'
        for k in range(ISOFORM_INDEX_START,K):
            for i in iid_to_rids[k]:
                rid_to_corrections[i] = [str(MATRIX[i][j]) for j in range(M)]
                for j in range(0,M):
                    if C[i][j] == 1 and OBJ[i][j][k].getAttr(GRB.Attr.X) > 0.9:
                        rid_to_corrections[i][j] = 'X'

        assert sum([len(x) for x in iid_to_rids.values()])==len(RIDS), 'Some reads have been assigned to more than one isoform!'
    log_file.close()
    return(status,iid_to_isoform,iid_to_rids,rid_to_corrections)

def output_isoform_clusters(read_name_path, clusters, matrix, isoforms, unaligned_gaps, rid_corrected_exons, exons_lengths, outpath):
    rid_to_name = [name.rstrip() for name in open(read_name_path, 'r')]
    out_file = open(outpath, 'w+')
    for k,rids in clusters.items():
        output = list()
        output.append('isoform_{:03d}'.format(k))
        output.append(''.join([str(e) for e in isoforms[k]]))
        output.extend([str(l*isoforms[k][eid]) for eid,l in enumerate(exons_lengths)])
        out_file.write('#{}\n'.format('\t'.join(output)))
        for i in rids:
            output = list()
            output.append('isoform_{:03d}'.format(k))
            output.append(str(rid_to_name[i]))
            output.append(str(i))
            output.append(''.join([str(e) for e in matrix[i]]))
            exon_strs = [str(x) for x in rid_corrected_exons[i]]
            for (j1,j2,l) in unaligned_gaps[i]:
                exon_strs[j1]+='({})'.format(l)
            output.extend(exon_strs)
            out_file.write('{}\n'.format('\t'.join(output)))
    out_file.close()

def make_polyAcut_info_lookup_table(names_file,cuts_file):

    name2index = {}
    lookup = {}
    with open(names_file, 'r') as nh , open(cuts_file, 'r') as ch:
        for index, line in enumerate(nh):
            name2index[line.rstrip()] = index

        for line in ch:
            name, cut_info = fields = line.rstrip().split("\t")
            if name in name2index:
                lookup[name2index[name]] = cut_info

    return lookup

def main():
    args = parse_args()

    exon_lens         = read_exon_lengths(exons_file=args.exons_tsv)
    M = len(exon_lens)
    matrix            = read_matrix(args.data_matrix, M =M)
    rids              = sorted(matrix.keys())
    N = len(rids)
    rid_to_unaln_gaps = read_unaligned_gaps(gaps_file=args.unaligned_gaps, nb_exons=M, rids=matrix.keys())
    incomp_rids       = list()
    if args.incomp_read_pairs != '':
        incomp_rids   = read_incomp_read_pairs(pairs_file=args.incomp_read_pairs, rids=matrix.keys())

    polyA_trimmed = make_polyAcut_info_lookup_table(args.valid_read_names, args.poly_tail_status)

    round = -1
    iid_to_isoform     = dict()
    iid_to_rids        = dict()
    rid_to_corrections = dict()
    while True:
        round += 1
        print('Running {}-th round with {} reads...'.format(round, len(rids)))
        status,iid_to_isoform_cur,iid_to_rids_cur,rid_to_corrections_cur = run_ILP(
            polyA_trimmed   = polyA_trimmed,
            RIDS            = rids,
            MATRIX          = matrix,
            EXON_L          = exon_lens,
            GAP_L           = rid_to_unaln_gaps,
            INCOMP_RIDS     = args.incomp_read_pairs,
            K               = args.isoform_count,
            EPSILON         = args.epsilon,
            OFFSET          = args.unaligned_offset,
            garbage_isoform = args.garbage_isoform,
            order_isoforms  = args.order_isoforms,
            timeout         = args.timeout,
            threads         = args.threads,
            out_prefix      = '{}.gurobi_logs/round.{}'.format(args.out_prefix, round)
        )
        if status != 'OPTIMAL':
            break
        for iid_cur in iid_to_isoform_cur.keys():
            if args.recycle_garbage and iid_cur == 0:
                continue
            iid = len(iid_to_isoform)
            assert not iid in iid_to_isoform, 'Isoform ids should be sequential'
            iid_to_isoform[iid]     = list(iid_to_isoform_cur[iid_cur])
            iid_to_rids[iid]        = list(iid_to_rids_cur[iid_cur])
            for rid in iid_to_rids[iid]:
                rid_to_corrections[rid] = rid_to_corrections_cur[rid]
        if args.recycle_garbage == False or round >= MAX_ROUND:
            break
        rids = sorted(iid_to_rids_cur[0])
        if len(rids) < MIN_READS:
            break
    output_isoform_clusters(
        read_name_path      = args.valid_read_names,
        clusters            = iid_to_rids,
        isoforms            = iid_to_isoform,
        matrix              = matrix,
        rid_corrected_exons = rid_to_corrections,
        unaligned_gaps      = rid_to_unaln_gaps,
        exons_lengths       = exon_lens,
        outpath             = '{}.tsv'.format(args.out_prefix)
    )


if __name__ == "__main__":
    main()
