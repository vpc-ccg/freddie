#!/usr/bin/env python3
import os
import argparse
from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    recycle_models = ['constant', 'exons', 'introns', 'relative']
    parser.add_argument("-tr",
                        "--transcripts",
                        type=str,
                        required=True,
                        help="Path to TSV file of transcripts")
    parser.add_argument("-d",
                        "--data",
                        type=str,
                        required=True,
                        help="Path to DATA file of reads by segments matrix")
    parser.add_argument("-s",
                        "--segments",
                        type=str,
                        required=True,
                        help="Path to segments TXT file")
    parser.add_argument("-rm",
                        "--recycle-model",
                        type=str,
                        default='constant',
                        help="Model type: {}. Default: {}".format(', '.join(recycle_models), 'constant'))
    parser.add_argument("-g",
                        "--gaps",
                        type=str,
                        required=True,
                        help="Path to unligned gap tuples per read including tail information")
    parser.add_argument("-go",
                        "--gap-offset",
                        type=int,
                        default=20,
                        help="Slack +- value for exons and the unaligned gaps")
    parser.add_argument("-e",
                        "--epsilon",
                        type=float,
                        default=0.2,
                        help="Epsilon value for how much can unaligned gaps can cover")
    parser.add_argument("-names",
                        type=str,
                        required=True,
                        help="Valid read ids to do ilp"),
    parser.add_argument("-mr",
                        "--max-rounds",
                        type=int,
                        default=20,
                        help="Maximum number of ILP rounds. Default {}".format(20))
    parser.add_argument("-is",
                        "--min-isoform-size",
                        type=int,
                        default=5,
                        help="Minimum isoform size. Default {}".format(5))
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
    assert args.recycle_model in recycle_models
    return(args)

def get_data(data, segs, reads):
    for rid,line in enumerate(open(data)):
        line = line.rstrip()
        reads[rid]['data']  = [int(x) for x in line]
        assert len(reads[rid]['data'])==len(segs), 'NR ({}):  Length if each row in data segment data must be equal to M ({})'.format(rid,len(segs))
        assert all(x in [0,1,2] for x in reads[rid]['data']), 'NR ({}):  {} is not a valid entry in segment data'.format(rid,x)

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

def get_segments(segments):
    segs = [int(l.rstrip()) for l in open(segments)]
    segs = [(s,e,e-s) for s,e in zip(segs[:-1],segs[1:])]
    assert all(l>0 and e>s>=0 for s,e,l in segs), 'Something is wrong with segments_txt file'
    return segs

def get_gaps(gaps, segs, reads):
    for read,line in enumerate(open(gaps)):
        line = line.rstrip()
        read = reads[read]
        read['gaps']=dict()
        read['tail']=dict(
            SA  = (0,0),
            ST  = (0,0),
            EA  = (0,0),
            ET  = (0,0),
            SSC = 0,
            ESC = 0,
        )
        if len(line)==0:
            continue
        for gap in line.split('\t'):
            if any(gap.startswith('{}:'.format(x)) for x in read['tail'].keys()):
                gap=gap.split(':')
                assert len(gap)==2
                sc = gap[0]
                size = int(gap[1])
                assert size>=0
                read['tail'][sc]=size
            elif any(gap.startswith('{}_'.format(x)) for x in read['tail'].keys()):
                gap=gap.split(':')
                assert len(gap)==2
                sc = gap[0].split('_')[0]
                p_size = int(gap[0].split('_')[1])
                size = int(gap[1])
                read['tail'][sc]=(p_size,size)
            else:
                interval,size = gap.split(':')
                interval = tuple([int(x) for x in interval.split('-')])
                assert len(interval)==2
                size=int(size)
                assert 0 <= interval[0] < len(segs), 'Exon id ({eid}) of read id ({rid}) are between 0 and {M}'.format(eid=interval[0],rid=rid, M=len(segs))
                assert 0 <= interval[1] < len(segs), 'Exon id ({eid}) of read id ({rid}) are between 0 and {M}'.format(eid=interval[1],rid=rid, M=len(segs))
                assert size>=0 , 'Unaligned length ({}) must be non-negative'.format(gap[2])
                read['gaps'][interval]=size

def find_segment_read(M,i):
    min_i = -1
    max_i = len(M[i])-1
    for j in range(len(M[i])):
        if min_i==-1 and M[i][j]==1:
            min_i = j
        if M[i][j]==1:
            max_i = j
    return((min_i,max_i))

def garbage_cost_introns(C):
    # Sum of 1, i.e. exons that could be corrected, minus 0.5 to make sure that assigning anyway
    # to an isoform with maximum correction score is not at least as good as assigning to the garbage isoform
    return(max(sum(C.values())-0.5,1))
def garbage_cost_exons(I):
    # Sum of exons present in the read
    return(max(sum(I.values())-0.5,1))
# TO DO: think about better ways to define the cost to assign to the garbagte isoform

def preprocess_ilp(reads, segs, ilp_settings):
    print('Preproessing ILP with {} reads and the following settings:\n{}'.format(len(reads), ilp_settings))
    N = len(reads)
    M = len(segs)
    I = dict()
    C = dict()
    for i in range(N):
        I[i] = dict()
        for j in range(0,M):
            I[i][j] = reads[i]['data'][j]%2
        C[i] = dict()
        (min_i,max_i) = find_segment_read(I,i)
        if ilp_settings['strand']=='-' and reads[i]['tail']['ST'][0]>10:
            reads[i]['gaps'][(-1,min_i)] = reads[i]['tail']['ST'][1]
            min_i = 0
        elif ilp_settings['strand']=='+' and reads[i]['tail']['EA'][0]>10:
            reads[i]['gaps'][(max_i,M)] = reads[i]['tail']['EA'][1]
            max_i = M-1
        for j in range(0,M):
            if min_i <= j <= max_i and reads[i]['data'][j]==0:
                C[i][j]   = 1
            else:
                C[i][j]   = 0
    # List of incompatible pairs of reads
    INCOMP_READ_PAIRS_AUX1 = list()
    # Read i1 and i2 have a 1 and a 2 in a given position j
    for i1 in range(N):
        for i2 in range(i1+1,N):
            incomp = False
            for j in range(M):
                if reads[i1]['data'][j]*reads[i2]['data'][j]==2:
                    incomp=True
                    break
            if incomp:
                INCOMP_READ_PAIRS_AUX1.append((i1,i2))
    print('# Number of incompatible read pairs due to a (1,2) pattern: {X}\n'.format(X=len(INCOMP_READ_PAIRS_AUX1)))

    # If r1 and r2 don't have any 1 in common, then they are incompatible
    INCOMP_READ_PAIRS_AUX2 = list()
    for i1 in range(N):
        for i2 in range(i1+1,N):
            incomp = True
            for j in range(M):
                if reads[i1]['data'][j]*reads[i2]['data'][j]==1:
                    incomp=False
                    break
            if incomp:
                INCOMP_READ_PAIRS_AUX2.append((i1,i2))
    print('# Number of incompatible read pairs due to not having any common exons: {X}\n'.format(X=len(INCOMP_READ_PAIRS_AUX2)))

    INCOMP_READ_PAIRS = list(set(INCOMP_READ_PAIRS_AUX1)|set(INCOMP_READ_PAIRS_AUX2))
    print('# Total number of incompatible read pairs: {X}\n'.format(X=len(INCOMP_READ_PAIRS)))

    # Assigning a cost to the assignment of reads to the garbage isoform
    GARBAGE_COST = {}
    for i in range(N):
        if ilp_settings['recycle_model'] == 'exons':
            GARBAGE_COST[i] = garbage_cost_exons(I=I[i])
        elif ilp_settings['recycle_model'] == 'introns':
            GARBAGE_COST[i] = garbage_cost_introns(C=C[i])
        elif ilp_settings['recycle_model'] == 'constant':
            GARBAGE_COST[i] = 1
    ilp_settings['I'] = I
    ilp_settings['C'] = C
    ilp_settings['INCOMP_READ_PAIRS'] = INCOMP_READ_PAIRS
    ilp_settings['GARBAGE_COST'] = GARBAGE_COST

def run_ILP(reads, remaining_rids, segs, ilp_settings, out_prefix):
    os.makedirs(out_prefix[:out_prefix.rfind('/')], exist_ok=True)
    # Variables directly based on the input ------------------------------------
    # I[i,j] = 1 if reads[i]['data'][j]==1 and 0 if reads[i]['data'][j]==0 or 2
    # C[i,j] = 1 if exon j is between the first and last exons (inclusively)
    #   covered by read i and is not in read i but can be turned into a 1
    ISOFORM_INDEX_START = 1
    M = len(segs)
    MAX_ISOFORM_LG = sum(seg[2] for seg in segs)
    I                 = ilp_settings['I']
    C                 = ilp_settings['C']
    INCOMP_READ_PAIRS = ilp_settings['INCOMP_READ_PAIRS']
    GARBAGE_COST      = ilp_settings['GARBAGE_COST']

    # ILP model ------------------------------------------------------
    ILP_ISOFORMS = Model('isoforms_v6_20200418')
    ILP_ISOFORMS.setParam(GRB.Param.Threads, ilp_settings['threads'])

    # Decision variables
    # R2I[i,k] = 1 if read i assigned to isoform k
    R2I    = {}
    R2I_C1 = {} # Constraint enforcing that each read is assigned to exactly one isoform
    for i in remaining_rids:
        R2I[i] = {}
        for k in range(ilp_settings['K']):
            R2I[i][k] = ILP_ISOFORMS.addVar(
                vtype=GRB.BINARY,
                name='R2I[{i}][{k}]'.format(i=i,k=k)
            )
        R2I_C1[i] = ILP_ISOFORMS.addLConstr(
            lhs   = quicksum(R2I[i][k] for k in range(0,ilp_settings['K'])),
            sense = GRB.EQUAL,
            rhs   = 1,
            name  = 'R2I_C1[{i}]'.format(i=i)
        )

    # Implied variable: canonical exons presence in isoforms
    # E2I[j,k]     = 1 if canonical exon j is in isoform k
    # E2I_min[j,k] = 1 if canonical exon j is in isoform k and is shared by all reads of that isoform
    # E2IR[j,k,i]  = 1 if read i assigned to isoform k AND exon j covered by read i
    # Auxiliary variable
    # E2IR[j,k,i]  = R2I[i,k] AND I[i,j]
    # E2I[j,k]     = max over  all reads i of E2IR[j,k,i]
    # E2I_min[j,k] = min over  all reads i of E2IR[j,k,i]
    E2I     = {}
    E2I_C1  = {}
    E2I_min = {}
    E2I_min_C1  = {}
    E2IR    = {}
    E2IR_C1 = {}
    for j in range(0,M):
        E2I[j]     = {}
        E2I_C1[j]  = {}
        E2I_min[j]     = {}
        E2I_min_C1[j]  = {}
        E2IR[j]    = {}
        E2IR_C1[j] = {}
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
        # We start assigning exons from the first isoform
        for k in range(ISOFORM_INDEX_START, ilp_settings['K']):
            E2I[j][k]     = ILP_ISOFORMS.addVar(
                vtype = GRB.BINARY,
                name  = 'E2I[{j}][{k}]'.format(j=j,k=k)
            )
            E2I_min[j][k]     = ILP_ISOFORMS.addVar(
                vtype = GRB.BINARY,
                name  = 'E2I_min[{j}][{k}]'.format(j=j,k=k)
            )
            E2IR[j][k]    = {}
            E2IR_C1[j][k] = {}
            for  i in remaining_rids:
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
                vars     = [E2IR[j][k][i] for i in remaining_rids],
                constant = 0.0,
                name     = 'E2I_C1[{j}][{k}]'.format(j=j,k=k)
            )
            E2I_min_C1[j][k] = ILP_ISOFORMS.addGenConstrMin(
                resvar   = E2I_min[j][k],
                vars     = [E2IR[j][k][i] for i in remaining_rids],
                constant = 0.0,
                name     = 'E2I_min_C1[{j}][{k}]'.format(j=j,k=k)
            )

    # Adding constraints for unaligned gaps
    # If read i is assigned to isoform k, and reads[i]['gaps'] contains ((j1,j2),l), and
    # the sum of the lengths of exons in isoform k between exons j1 and j2 is L
    # then (1-EPSILON)L <= l <= (1+EPSILON)L
    GAPI    = {} # GAPI[(j1,j2,k)] = sum of the length of the exons between exons j1 and j2 (inclusively) in isoform k
    GAPI_C1 = {} # Constraint fixing the value of GAPI
    GAPR_C1 = {} # Constraint ensuring that the unaligned gap is not too short for every isoform and gap
    GAPR_C2 = {} # Constraint ensuring that the unaligned gap is not too long for every isoform and gap
    for i in remaining_rids:
        for ((j1,j2),l) in reads[i]['gaps'].items():
            for k in range(ISOFORM_INDEX_START,ilp_settings['K']): # No such constraint on the garbage isoform if any
                if not (j1,j2,k) in GAPI:
                    GAPI[(j1,j2,k)]      = ILP_ISOFORMS.addVar(
                        vtype = GRB.INTEGER,
                        name  = 'GAPI[({j1},{j2},{k})]'.format(j1=j1,j2=j2,k=k)
                    )
                    GAPI_C1[(j1,j2,k)]   = ILP_ISOFORMS.addLConstr(
                        lhs   = GAPI[(j1,j2,k)],
                        sense = GRB.EQUAL,
                        rhs   = quicksum(E2I[j][k]*segs[j][2] for j in range(j1+1,j2)),
                        name  = 'GAPI_C1[({j1},{j2},{k})]'.format(j1=j1,j2=j2,k=k)
                    )
                GAPR_C1[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr(
                    lhs   = (1.0-ilp_settings['epsilon'])*GAPI[(j1,j2,k)]-ilp_settings['offset']-((1-R2I[i][k])*MAX_ISOFORM_LG),
                    sense = GRB.LESS_EQUAL,
                    rhs   = l,
                    name  = 'GAPR_C1[({i},{j1},{j2},{k})]'.format(i=i,j1=j1,j2=j2,k=k)
                )
                GAPR_C2[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr(
                    lhs   = (1.0+ilp_settings['epsilon'])*GAPI[(j1,j2,k)]+ilp_settings['offset']+((1-R2I[i][k])*MAX_ISOFORM_LG),
                    sense = GRB.GREATER_EQUAL,
                    rhs   = l,
                    name  = 'GAPR_C2[({i},{j1},{j2},{k})]'.format(i=i,j1=j1,j2=j2,k=k)
                )
    # Adding constraints for incompatible read pairs
    INCOMP_READ_PAIRS_C1 = {}
    for (i1,i2) in INCOMP_READ_PAIRS:
        if not (i1 in remaining_rids and i2 in remaining_rids):
            continue
        for k in range(ISOFORM_INDEX_START,ilp_settings['K']): # Again, no such constraint on the garbage isoform if any
            INCOMP_READ_PAIRS_C1[(i1,i2,k)] = ILP_ISOFORMS.addLConstr(
                lhs   = R2I[i1][k]+R2I[i2][k],
                sense = GRB.LESS_EQUAL,
                rhs   = 1,
                name  = 'INCOMP_READ_PAIRS_C1[({i1},{i2},{k})]'.format(i1 =i1,i2 =i2,k =k)
            )

    # [OPTIONAL] Labeling non-garbage isoforms by their exon content and forcing them to occur in increasing label order
    # LABEL_I    = {}
    # LABEL_I_C1 = {}
    # LABEL_I_C2 = {}
    # for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
    #     LABEL_I[k]    = ILP_ISOFORMS.addVar(
    #         vtype = GRB.INTEGER,
    #         name  = 'LABEL_I[{k}]'.format(k=k)
    #     )
    #     LABEL_I_C1[k] = ILP_ISOFORMS.addLConstr(
    #         lhs   = LABEL_I[k],
    #         sense = GRB.EQUAL,
    #         rhs   = quicksum(E2I[j][k]*(2**j) for j in range(0,M)),
    #         name  = 'LABEL_I_C1[{k}]'.format(k =k)
    #     )
    #     if k > ISOFORM_INDEX_START:
    #         LABEL_I_C2[k] = ILP_ISOFORMS.addLConstr(
    #             lhs   = LABEL_I[k],
    #             sense = GRB.LESS_EQUAL,
    #             rhs   = LABEL_I[k-1]-0.1,
    #             name  = 'LABEL_I_C2[{k}]'.format(k=k)
    #         )

    # Objective function
    # For i,j,k such that i ∈ remaining_rids, C[i,j]=1 (read has a zero that can be
    #   corrected), and E2I[j,k]=1 (isoform k has exon j), OBJ[i][j][k] = 1
    OBJ    = {}
    OBJ_C1 = {}
    OBJ_SUM = LinExpr(0.0)
    for i in remaining_rids:
        OBJ[i]    = {}
        OBJ_C1[i] = {}
        for j in range(0,M):
            if C[i][j]==1: # 1 if exon j not in read i but can be added to it
                OBJ[i][j]    = {}
                OBJ_C1[i][j] = {}
                for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
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
    GAR_OBJ   = {}
    GAR_OBJ_C = {}
    for i in remaining_rids:
        if ilp_settings['recycle_model'] in ['constant', 'exons', 'introns']:
            OBJ_SUM.addTerms(
                coeffs = 1.0*GARBAGE_COST[i],
                vars   = R2I[i][0]
            )
        elif RE_MOD == 'relative':
            GAR_OBJ[i]   = {}
            GAR_OBJ_C[i] = {}
            for j in range(0,M):
                GAR_OBJ[i][j]   = {}
                GAR_OBJ_C[i][j] = {}
                for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
                    if I[i][j] == 1:
                        GAR_OBJ[i][j][k]    = ILP_ISOFORMS.addVar(
                            vtype = GRB.BINARY,
                            name  = 'GAR_OBJ[{i}][{j}][{k}]'.format(i=i,j=j,k=k)
                        )
                        GAR_OBJ_C[i][j][k] = ILP_ISOFORMS.addGenConstrAnd(
                            resvar = GAR_OBJ[i][j][k],
                            vars   = [R2I[i][0],E2I_min[j][k]],
                            name   = 'GAR_OBJ_C[{i}][{j}][{k}]'.format(i=i,j=j,k=k)
                        )
                        OBJ_SUM.addTerms(
                            coeffs = 1.0,
                            vars   = GAR_OBJ[i][j][k]
                        )
                    elif I[i][j] == 0 and C[i][j] == 1:
                        pass

    ILP_ISOFORMS.setObjective(
        expr  = OBJ_SUM,
        sense = GRB.MINIMIZE
    )

    ILP_ISOFORMS.write('{}.lp'.format(out_prefix))

    # Optimization
    #ILP_ISOFORMS.Params.PoolSearchMode=2
    #ILP_ISOFORMS.Params.PoolSolutions=5
    ILP_ISOFORMS.setParam('LogFile', '{}.glog'.format(out_prefix))
    ILP_ISOFORMS.setParam('TimeLimit', ilp_settings['timeout']*60)
    ILP_ISOFORMS.optimize()

    ILP_ISOFORMS_STATUS = ILP_ISOFORMS.Status

    isoforms = {k:dict() for k in range(ISOFORM_INDEX_START,ilp_settings['K'])}
    print('STATUS: {}'.format(ILP_ISOFORMS_STATUS))
    if ILP_ISOFORMS_STATUS == GRB.Status.INF_OR_UNBD or \
           ILP_ISOFORMS_STATUS == GRB.Status.INFEASIBLE or \
           ILP_ISOFORMS_STATUS == GRB.Status.UNBOUNDED:
        status = 'NO_SOLUTION'
        ILP_ISOFORMS.computeIIS()
        ILP_ISOFORMS.write('{}.ilp'.format(out_prefix))
    elif ILP_ISOFORMS_STATUS == GRB.Status.OPTIMAL or \
            ILP_ISOFORMS_STATUS == GRB.Status.SUBOPTIMAL or \
            ILP_ISOFORMS_STATUS == GRB.Status.TIME_LIMIT:
        if ILP_ISOFORMS_STATUS == GRB.Status.SUBOPTIMAL:
            status = 'SUBOPTIMAL'
        elif ILP_ISOFORMS_STATUS == GRB.Status.TIME_LIMIT:
            status = 'TIME_LIMIT'
        else:
            status = 'OPTIMAL'
        # Writing the optimal solution to disk
        solution_file = open('{}.sol'.format(out_prefix), 'w+')
        for v in ILP_ISOFORMS.getVars():
            solution_file.write('{}\t{}\n'.format(v.VarName, v.X))
        solution_file.close()
        # Isoform id to isoform structure
        for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
            isoforms[k]['exons'] = [int(E2I[j][k].getAttr(GRB.Attr.X)>0.9) for j in range(0,M)]
            isoforms[k]['rid_to_corrections'] = dict()
        # Isoform id to read ids set
        for i in remaining_rids:
            isoform_id = -1
            for k in range(0,ilp_settings['K']):
                if R2I[i][k].getAttr(GRB.Attr.X) > 0.9:
                    assert isoform_id == -1,'Read {} has been assigned to multiple isoforms!'.format(i)
                    isoform_id = k
            assert isoform_id!=-1,'Read {} has not been assigned to any isoform!'.format(i)
            if isoform_id == 0:
                continue
            isoforms[isoform_id]['rid_to_corrections'][i]=-1
        # Read id to its exon corrections
        for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
            for i in isoforms[k]['rid_to_corrections'].keys():
                isoforms[k]['rid_to_corrections'][i] = [str(reads[i]['data'][j]) for j in range(M)]
                for j in range(0,M):
                    if C[i][j] == 1 and OBJ[i][j][k].getAttr(GRB.Attr.X) > 0.9:
                        isoforms[k]['rid_to_corrections'][i][j] = 'X'
    return status,isoforms

def output_isoforms(isoforms, reads, garbage_rids, segs, outpath):
    out_file = open(outpath, 'w+')
    for idx,isoform in enumerate(isoforms):
        output = list()
        isoform_name = 'isoform_{:03d}'.format(idx)
        output.append(isoform_name)
        output.append('.')
        output.append('.')
        output.append(''.join([str(e) for e in isoform['exons']]))
        output.extend([str(seg[2]*isoform['exons'][eid]) for eid,seg in enumerate(segs)])
        out_file.write('#{}\n'.format('\t'.join(output)))
        for i,corrections in isoform['rid_to_corrections'].items():
            output = list()
            output.append(isoform_name)
            output.append(reads[i]['name'])
            output.append(str(i))
            output.append(''.join([str(x) for x in corrections]))
            exon_strs = [str(x) for x in corrections]
            for (j1,j2),l in reads[i]['gaps'].items():
                exon_strs[j1]+='({})'.format(l)
            output.extend(exon_strs)
            for k,v in sorted(reads[i]['tail'].items()):
                output.append('{}:{}'.format(k,v))
            out_file.write('{}\n'.format('\t'.join(output)))
    for i in garbage_rids:
        output = list()
        output.append('isoform_GARB')
        output.append(reads[i]['name'])
        output.append(str(i))
        output.append(''.join([str(x) for x in reads[i]['data']]))
        exon_strs = [str(x) for x in reads[i]['data']]
        for (j1,j2),l in reads[i]['gaps'].items():
            exon_strs[j1]+='({})'.format(l)
        output.extend(exon_strs)
        for k,v in sorted(reads[i]['tail'].items()):
            output.append('{}:{}'.format(k,v))
        out_file.write('{}\n'.format('\t'.join(output)))
    out_file.close()

def main():
    args = parse_args()
    segs = get_segments(segments=args.segments)
    reads = [dict(name=name.rstrip()) for name in open(args.names)]
    if len(reads)==0:
        print('No reads for this gene!')
        out_file = open('{}.tsv'.format(args.out_prefix),'w+')
        out_file.close()
        exit()
    get_data(data=args.data, segs=segs, reads=reads)
    get_gaps(gaps=args.gaps, segs=segs, reads=reads)

    strand = list({l.rstrip().split('\t')[2] for l in open(args.transcripts)})
    assert len(strand)==1
    strand = strand[0]

    print('Matrix size: {}x{}'.format(len(reads),len(segs)))

    remaining_rids = set(rid for rid in range(len(reads)))
    for read in reads:
        read['isoform']=-1

    isoforms           = list()
    ilp_settings=dict(
        recycle_model   = args.recycle_model,
        strand  = strand,
        K       = 2,
        epsilon = args.epsilon,
        offset  = args.gap_offset,
        timeout = args.timeout,
        threads = args.threads,
    )
    preprocess_ilp(reads, segs, ilp_settings)
    print('ILP params: {}'.format(ilp_settings.keys()))
    for round in range(args.max_rounds):
        print('==========\nRunning {}-th round with {} reads...'.format(round, len(remaining_rids)))
        if len(remaining_rids) < args.min_isoform_size:
            break
        status,round_isoforms = run_ILP(
            reads          = reads,
            remaining_rids = remaining_rids,
            segs           = segs,
            ilp_settings   = ilp_settings,
            out_prefix     = '{}.gurobi_logs/round.{}'.format(args.out_prefix, round)
        )
        if status != 'OPTIMAL' or max([0]+[len(i['rid_to_corrections']) for i in round_isoforms.values()]) < args.min_isoform_size:
            break
        for k,isoform in round_isoforms.items():
            print('Isoform {} size: {}'.format(k, len(isoform['rid_to_corrections'])))
            if len(isoform['rid_to_corrections']) < args.min_isoform_size:
                continue
            isoforms.append(isoform)
            for rid,corrections in isoform['rid_to_corrections'].items():
                assert rid in remaining_rids
                remaining_rids.remove(rid)
                reads[rid]['corrections']=corrections
                reads[rid]['isoform']=len(isoforms)-1
        print('------->')
        print('Remaining reads: {}\n'.format(len(remaining_rids)))
        print('<-------')

    output_isoforms(
        isoforms     = isoforms,
        reads        = reads,
        garbage_rids = remaining_rids,
        segs         = segs,
        outpath      = '{}.tsv'.format(args.out_prefix)
    )


if __name__ == "__main__":
    main()
