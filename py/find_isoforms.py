import os
import argparse
from gurobipy import *

def parse_args():
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
    parser.add_argument("-e",
                        "--epsilon",
                        type=float,
                        default=0.2,
                        help="Epsilon value for how much can unaligned gaps can cover")
    parser.add_argument("-irp",
                        "--incomp-read-pairs",
                        type=str,
                        default='./incomp_read_pairs',
                        help="Path to file of read pairs that can not originate from the same isoform")
    parser.add_argument("-oi",
                        "--order-isoforms",
                        type=bool,
                        default=False,
                        help="Force to label isoforms and order them increasingly")
    parser.add_argument("-to",
                        "--timeout",
                        type=int,
                        default=15,
                        help="Gurobi time-out")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=8,
                        help="Number of threads to use")
    parser.add_argument("-op",
                        "--out_prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def read_matrix(data_matrix):
    result = list()
    for line in open(data_matrix):
        line = line.rstrip()
        row  = [int(x) for x in line]
        assert (1 in row), 'Some row of the data matrix has only 0s'
        result.append(row)
    return(result)

def read_incomp_read_pairs(pairs_file,nb_reads):
    # Assumption: read pairs are represented by pairs of integers indicating the index in the data matrix of the incompatible reads
    result = list()
    for line in open(pairs_file):
        line = line.rstrip().split()
        assert len(line)==2, 'Incompatible read pairs are pairs of reads'
        row  = [int(x) for x in line]
        assert row[0]>=0 and row[0]<nb_reads and row[1]>=0 and row[1]<nb_reads and row[0]!=row[1], 'Incompatible read pairs are pairs of different read index'
        result.append(row)
    return(result)

def read_exon_length(exons_file,nb_exons):
    result = [int(r.split('\t')[1])-int(r.split('\t')[0]) for r in open(exons_file).readlines()]
    assert len(result) == nb_exons, 'Number of lines should be equal to number of exons ({})'.format(nb_exons)
    for e in result:
        assert e>0, 'All exons must be positive'
    return result

# TO WRITE
def read_unaligned_gaps(gaps_file,nb_reads,nb_exons):
    result = list()
    for rid,line in enumerate(open(gaps_file)):
        line = line.rstrip()
        if len(line) == 0:
            continue
        for gap in line.split('\t'):
            gap = [int(x) for x in gap.split('-')]
            assert len(gap)==3, 'Each gap must be a 3-tuple'
            assert gap[0]>=0 and gap[0]<nb_exons, 'Exon id ({}) are between 0 and {}'.format(gap[0],nb_exons)
            assert gap[1]>=0 and gap[1]<nb_exons, 'Exon id ({}) are between 0 and {}'.format(gap[1],nb_exons)
            assert gap[2]>=0 , 'Unaligned length ({}) must be non-negative {}'.format(gap[2])
            result.append([rid]+gap)
    assert rid+1==nb_reads, 'Number of lines should match number of reads'
    return(result)

def find_segment_read(M,i):
    min_i = -1
    for j in range(len(M[i])):
        if min_i==-1 and M[i][j]==1:
            min_i = j
        if M[i][j]==1:
            max_i = j
    return((min_i,max_i))

def run_ILP(data_matrix, exons_length, unaligned_gaps, K, epsilon, incomp_read_pairs, timeout, order_isoforms, threads, out_file, log_file_name):
    INPUT_1 = read_matrix(data_matrix) # 0/1/2 matrix
    N = len(INPUT_1)
    M = len(INPUT_1[0])
    # log_file.write('# Input Matrix has '+str(N)+' rows(reads) and '+str(M)+' columns(canonical exons)\n')

    INPUT_2 = read_exon_length(exons_length,M)
    print INPUT_2
    # Format: INPUT_2[j] = length of exon j
    MAX_ISOFORM_LG = sum(INPUT_2)

    INPUT_3 = read_unaligned_gaps(unaligned_gaps,N,M)
    # Format: INPUT_3 = list of (i,j1,j2,l) where l = length of unaligned gap of read i between exons j1 and j2

    if not os.path.isfile(incomp_read_pairs):
        # log_file.write('# Empty list of incompatible read pairs\n')
        INPUT_4 = list()
    else:
        INPUT_4 = read_incomp_read_pairs(incomp_read_pairs,N)

    # Variables directly based on the input
    # I[i,j] = 1 if INPUT_1[i][j]==1 and 0 if INPUT_1[i][j]==0 or 2
    # C[i,j] = 1 if exon j is between the first and last exons (included) covered by read i and is not in read i but can be turned into a 1
    I = {}
    C = {}
    UB_NB_CORRECTIONS = 0 # Upper bound on the number of 0s that can be corrected into 1
    for i in range(0,N):
        I[i] = {}
        for j in range(0,M):
            I[i][j] = INPUT_1[i][j]%2
        C[i] = {}
        (min_i,max_i) = find_segment_read(I,i)
        for j in range(0,M):
            if j>=min_i and j<=max_i and INPUT_1[i][j]==0:
                C[i][j]   = 1
                UB_NB_CORRECTIONS += 1
            else:
                C[i][j]   = 0
    # log_file.write('# Maximum number of 0 that can be corrected into 1: '+str(UB_NB_CORRECTIONS)+'\n')

    # List of incompatible pairs of reads
    INCOMP_READ_PAIRS_AUX1 = list()
    # Read i1 and i2 have a 1 and a 2 in a given position j
    for i1 in range(N):
        for i2 in range(i1,N):
            incomp = False
            for j in range(M):
                if INPUT_1[i1][j]*INPUT_1[i2][j]==2:
                    incomp=True
                    break
            if incomp:
                INCOMP_READ_PAIRS_AUX1.append((i1,i2))
    # log_file.write('# Number of incompatible read pairs due to a (1,2) pattern: '+str(len(INCOMP_READ_PAIRS_AUX1))+'\n')
    # From the list of provided incompatible pairs
    INCOMP_READ_PAIRS_AUX2 = list()
    for (i1,i2) in INPUT_4:
        if i1>i2:
            i1,i2=i2,i1
        INCOMP_READ_PAIRS_AUX2.append((i1,i2))
    # log_file.write('# Number of incompatible read pairs provided as input: '+str(len(INCOMP_READ_PAIRS_AUX2))+'\n')
    # Fusing both lists
    INCOMP_READ_PAIRS = list(set(INCOMP_READ_PAIRS_AUX1+INCOMP_READ_PAIRS_AUX2))
    # log_file.write('# Toral number of incompatible read pairs: '+str(len(INCOMP_READ_PAIRS))+'\n')


    # ILP model
    ILP_ISOFORMS = Model('isoforms_v3_25072019')
    ILP_ISOFORMS.setParam(GRB.Param.Threads, threads)

    # Decision variables
    # R2I[i,k] = 1 if read i assigned to isoform j
    R2I    = {}
    R2I_C1 = {} # Constraint enforcing that each read is assigned to exactly one isoform
    # Could be relaxed if we want to account for possibly wrong reads, i.e. reads not from the given gene
    for i in range(0,N):
        R2I[i] = {}
        for k in range(0,K):
            R2I[i][k] = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='R2I['+str(i)+']['+str(k)+']')
        R2I_C1[i] = ILP_ISOFORMS.addLConstr(quicksum(R2I[i][k] for k in range(0,K)),GRB.EQUAL,1,'R2I_C1['+str(i)+']')

    # Implied variable: canonical exons presence in isoforms
    # E2I[j,k]    = 1 if canonical exon j is in isoform k
    # Auxiliary variable
    # E2IR[j,k,i] = R2I[i,k] AND I[i,j] (read i assigned to isoform k and exon j covered by read i)
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
        for k in range(0,K):
            E2I[j][k]     = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='E2I['+str(j)+']['+str(k)+']')
            E2IR[j][k]    = {}
            E2IR_C1[j][k] = {}
            for  i in range(0,N):
                E2IR[j][k][i]    = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='E2IR['+str(j)+']['+str(k)+']['+str(i)+']')
                E2IR_C1[j][k][i] = ILP_ISOFORMS.addLConstr(E2IR[j][k][i],GRB.EQUAL,R2I[i][k]*I[i][j],'E2IR_C1['+str(j)+']['+str(k)+']['+str(i)+']')
            E2I_C1[j][k] = ILP_ISOFORMS.addGenConstrMax(E2I[j][k],[E2IR[j][k][i] for i in range(0,N)],0.0,'E2I_C1['+str(j)+']['+str(k)+']')

    # Adding constraints for unaligned gaps
    # If read i is assigned to isoform k, and INPUT_3 contains (i,j1,j2,l), and
    # the sum of the lengths of exons in isoform k between exons j1 and j2 is L
    # the (1-epsilon)L <= l <= (1+epsilon)L
    GAPI    = {} # GAPI[(j1,j2,k)] = sum of the length of the exons between exons j1 and j2 (exlusively) in isoform k
    GAPI_C1 = {} # Constraint fixing the value of GAPI
    GAPR_C1 = {} # Constraint ensuring that the unaligned gap is not too short for every isoform and gap
    GAPR_C2 = {} # Constraint ensuring that the unaligned gap is not too long for every isoform and gap
    print MAX_ISOFORM_LG
    print epsilon
    for (i,j1,j2,l) in INPUT_3:
        for k in range(0,K):
            GAPI[(j1,j2,k)]      = ILP_ISOFORMS.addVar(vtype=GRB.INTEGER,name='GAPI['+str(j1)+','+str(j2)+','+str(k)+')]')
            GAPI_C1[(j1,j2,k)]   = ILP_ISOFORMS.addLConstr(GAPI[(j1,j2,k)],GRB.EQUAL,quicksum(E2I[j][k]*INPUT_2[j] for j in range(j1+1,j2)),'GAPI_C1['+str(j1)+','+str(j2)+','+str(k)+')]')
            GAPR_C1[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr((1.0-epsilon)*GAPI[(j1,j2,k)]-((1-R2I[i][k])*MAX_ISOFORM_LG),GRB.LESS_EQUAL,l,'GAPR_C1['+str(i)+','+str(j1)+','+str(j2)+','+str(k)+')]')
            GAPR_C2[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr((1.0+epsilon)*GAPI[(j1,j2,k)]+((1-R2I[i][k])*MAX_ISOFORM_LG),GRB.GREATER_EQUAL,l,'GAPR_C2['+str(i)+','+str(j1)+','+str(j2)+','+str(k)+')]')
            # GAPR_C1[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr((1.0-epsilon)*(GAPI[(j1,j2,k)]-(1)*MAX_ISOFORM_LG),GRB.LESS_EQUAL,l,'GAPR_C1['+str(i)+','+str(j1)+','+str(j2)+','+str(k)+')]')
            # GAPR_C2[(i,j1,j2,k)] = ILP_ISOFORMS.addLConstr((1.0+epsilon)*(GAPI[(j1,j2,k)]+(1)*MAX_ISOFORM_LG),GRB.GREATER_EQUAL,l,'GAPR_C2['+str(i)+','+str(j1)+','+str(j2)+','+str(k)+')]')

    # Adding constraints for incompatible read pairs
    INCOMP_READ_PAIRS_C1 = {}
    for (i1,i2) in INCOMP_READ_PAIRS:
        for k in range(K):
            INCOMP_READ_PAIRS_C1[(i1,i2,k)] = ILP_ISOFORMS.addLConstr(R2I[i1][k]+R2I[i2][k],GRB.LESS_EQUAL,1,'INCOMP_READ_PAIRS_C1[('+str(i1)+','+str(i2)+','+str(k)+')]')

    # [OPTIONAL] Labeling strains by their exon content and forcing them to occur in increasing label order
    if order_isoforms==1:
        LABEL_I    = {}
        LABEL_I_C1 = {}
        LABEL_I_C2 = {}
        for k in range(0,K):
            LABEL_I[k]    = ILP_ISOFORMS.addVar(vtype=GRB.INTEGER,name='LABEL_I['+str(k)+']')
            LABEL_I_C1[k] = ILP_ISOFORMS.addLConstr(LABEL_I[k],GRB.EQUAL,quicksum(E2I[j][k]*(2**j) for j in range(0,M)),'LABEL_I_C1['+str(k)+']')
            if k > 0:
                LABEL_I_C2[k] = ILP_ISOFORMS.addLConstr(LABEL_I[k],GRB.LESS_EQUAL,LABEL_I[k-1]-0.1,'LABEL_I_C2['+str(k)+']')

    # Objective function
    # OBJ[i][j][k] = 1 if read i assigned to isoform k and exon j in isoform k, 0 otherwise
    OBJ    = {}
    OBJ_C1 = {}
    for i in range(0,N):
        OBJ[i]    = {}
        OBJ_C1[i] = {}
        for j in range(0,M):
            OBJ[i][j]    = {}
            OBJ_C1[i][j] = {}
            for k in range(0,K):
                OBJ[i][j][k]    = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='OBJ['+str(i)+']['+str(j)+']['+str(k)+']')
                OBJ_C1[i][j][k] = ILP_ISOFORMS.addGenConstrAnd(OBJ[i][j][k],[R2I[i][k],E2I[j][k]],'OBJ_C1['+str(i)+']['+str(j)+']['+str(k)+']')

    # Objective function
    OBJ_SUM = LinExpr(0.0)
    for i in range(0,N):
       for j in range(0,M):
           if (1-I[i][j])*C[i][j]==1: # 1 if exon j not in read i but can be added to it
               for k in range(0,K):
                   OBJ_SUM.addTerms(1.0, OBJ[i][j][k])
    ILP_ISOFORMS.setObjective(OBJ_SUM,GRB.MINIMIZE)

    ILP_ISOFORMS.write('TMP.lp')
    # Optimization
    #ILP_ISOFORMS.Params.PoolSearchMode=2
    #ILP_ISOFORMS.Params.PoolSolutions=5
    ILP_ISOFORMS.setParam('LogFile', log_file_name)
    ILP_ISOFORMS.setParam('TimeLimit', timeout*60)
    ILP_ISOFORMS.optimize()

    # Solution
    VAL_R2I = {}
    for i in range(0,N):
        VAL_R2I[i] = {}
        for k in range(0,K):
            VAL_R2I[i][k] = int(R2I[i][k].getAttr(GRB.Attr.X))
    VAL_E2I = {}
    for k in range(0,K):
        # VAL_E2I[k] = [int(E2I[j][k].getAttr(GRB.Attr.X)) for j in range(0,M)]
        # out_file.write(str(k)+':'+ str(VAL_E2I[k]) + '\n')
        for i in range(0,N):
            if VAL_R2I[i][k] == 1:
                out_file.write(str(i)+'\t'+str(k))
                out_file.write('\n')

def main():
    args = parse_args()
    print args
    out_file = open('{}.tsv'.format(args.out_prefix), 'w+')
    run_ILP(data_matrix       = args.data_matrix,
            exons_length      = args.exons_tsv,
            unaligned_gaps    = args.unaligned_gaps,
            K                 = args.isoform_count,
            epsilon           = args.epsilon,
            incomp_read_pairs = args.incomp_read_pairs,
            order_isoforms    = args.order_isoforms,
            timeout           = args.timeout,
            threads           = args.threads,
            out_file          = out_file,
            log_file_name     = args.out_prefix+'.glog')
    out_file.close()


if __name__ == "__main__":
    main()
