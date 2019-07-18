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

def find_segment_read(M,i):
    if not 1 in M:
        return(0,len(M))
    min_i = -1
    for j in range(len(M[i])):
        if min_i==-1 and M[i][j]==1:
            min_i = j
        if M[i][j]==1:
            max_i = j
    return((min_i,max_i))

def run_ILP(data_matrix, K, threads, out_file):
    INPUT_1 = read_matrix(data_matrix=data_matrix)
    N = len(INPUT_1)
    M = len(INPUT_1[0])

    ILP_ISOFORMS = Model('isoforms_v1')
    ILP_ISOFORMS.setParam(GRB.Param.Threads, threads)

    # Decision variables
    # R2I[i,k] = 1 if read i assigned to isoform j
    R2I    = {}
    R2I_C1 = {} # Constraint enforcing that each read is assigned to exactly one isoform
    for i in range(0,N):
        R2I[i] = {}
        for k in range(0,K):
            R2I[i][k] = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='R2I['+str(i)+']['+str(k)+']')
        R2I_C1[i] = ILP_ISOFORMS.addLConstr(quicksum(R2I[i][k] for k in range(0,K)),GRB.EQUAL,1,'R2I_C1['+str(i)+']')

    # Input variables
    # I[i,j]     = INPUT_1[i,j]
    # I_INV[i,j] = 1-INPUT_1[i,j]
    # C[i,j]     = 1 if exon j is between the first and last exons (included) covered by read i
    I        = {}
    I_C1     = {}
    I_INV    = {}
    I_INV_C1 = {}
    C        = {}
    C_C1    = {}
    for i in range(0,N):
        I[i]        = {}
        I_C1[i]     = {}
        I_INV[i]    = {}
        I_INV_C1[i] = {}
        C[i]        = {}
        C_C1[i]     = {}
        (min_i,max_i) = find_segment_read(INPUT_1,i)
        max_i=M
        for j in range(0,M):
            I[i][j]    = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='I['+str(i)+']['+str(j)+']')
            I_C1[i][j] = ILP_ISOFORMS.addLConstr(I[i][j],GRB.EQUAL,INPUT_1[i][j],'I_C1['+str(i)+']['+str(j)+']')
            C[i][j]    = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='C['+str(i)+']['+str(j)+']')
            if j>=min_i and j<=max_i:
                C_C1[i][j] = ILP_ISOFORMS.addLConstr(C[i][j],GRB.EQUAL,1,'C_C1['+str(i)+']['+str(j)+']')
            else:
                C_C1[i][j] = ILP_ISOFORMS.addLConstr(C[i][j],GRB.EQUAL,0,'C_C1['+str(i)+']['+str(j)+']')
            I_INV[i][j]    = ILP_ISOFORMS.addVar(vtype=GRB.BINARY,name='I_INV['+str(i)+']['+str(j)+']')
            I_INV_C1[i][j] = ILP_ISOFORMS.addLConstr(I_INV[i][j],GRB.EQUAL,1-I[i][j],'I_INV_C1['+str(i)+']['+str(j)+']')

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
                E2IR_C1[j][k][i] = ILP_ISOFORMS.addGenConstrAnd(E2IR[j][k][i],[R2I[i][k],I[i][j]],'E2IR_C1['+str(j)+']['+str(k)+']['+str(i)+']')
            E2I_C1[j][k] = ILP_ISOFORMS.addGenConstrMax(E2I[j][k],[E2IR[j][k][i] for i in range(0,N)],0.0,'E2I_C1['+str(j)+']['+str(k)+']')

    # Objective function
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
                OBJ_C1[i][j][k] = ILP_ISOFORMS.addGenConstrAnd(OBJ[i][j][k],[R2I[i][k],E2I[j][k],I_INV[i][j],C[i][j]],'OBJ_C1['+str(i)+']['+str(j)+']['+str(k)+']')

    # Objective function
    OBJ_SUM = LinExpr(0.0)
    for i in range(0,N):
       for j in range(0,M):
            for k in range(0,K):
                OBJ_SUM.addTerms(1.0, OBJ[i][j][k])
    ILP_ISOFORMS.setObjective(OBJ_SUM,GRB.MINIMIZE)

    # Optimization
    ILP_ISOFORMS.optimize()

    # Solution
    # VAL_I = {}
    # VAL_C = {}
    # for i in range(0,N):
    #     VAL_I[i] = {}
    #     VAL_C[i] = {}
    #     VAL_I[i] = [int(I[i][j].getAttr(GRB.Attr.X)) for j in range(0,M)]
    #     VAL_C[i] = [int(C[i][j].getAttr(GRB.Attr.X)) for j in range(0,M)]
    # out_file.write('Input')
    # for i in range(0,N):
    #     out_file.write(str(i)+'\t'+str(VAL_I[i])+'\t'+str(VAL_C[i]))

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

def read_matrix(data_matrix):
    result = list()
    for line in open(data_matrix):
        line = line.rstrip()
        result.append([int(x) for x in line])
    return result

def main():
    args = parse_args()

    out_file = open('{}.tsv'.format(args.out_prefix), 'w+')
    run_ILP(data_matrix=args.data_matrix, K=args.isoform_count, threads=args.threads, out_file=out_file)
    out_file.close()



if __name__ == "__main__":
    main()
