#!/usr/bin/env python3
import os
from multiprocessing import Pool

import argparse
import re

from networkx.algorithms import components
from networkx import Graph
from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    recycle_models = ['constant', 'exons', 'introns', 'relative']
    parser.add_argument("-s",
                        "--segment-tsv",
                        type=str,
                        required=True,
                        help="Path to Freddie segment TSV file of the reads")
    parser.add_argument("-rm",
                        "--recycle-model",
                        type=str,
                        default='constant',
                        help="Model type: {}. Default: {}".format(', '.join(recycle_models), 'constant'))
    parser.add_argument("-go",
                        "--gap-offset",
                        type=int,
                        default=20,
                        help="Slack +- value for exons and the unaligned gaps. Default: 20")
    parser.add_argument("-e",
                        "--epsilon",
                        type=float,
                        default=0.2,
                        help="Epsilon percent value for how much can unaligned gaps can cover. Default: 0.2")
    parser.add_argument("-mr",
                        "--max-rounds",
                        type=int,
                        default=30,
                        help="Maximum number of ILP rounds. Default {}".format(30))
    parser.add_argument("-is",
                        "--min-isoform-size",
                        type=int,
                        default=3,
                        help="Minimum isoform size in terms of number supporting reads. Default {}".format(3))
    parser.add_argument("-to",
                        "--timeout",
                        type=int,
                        default=4,
                        help="Gurobi time-out in minutes. Default: 4")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads to use")
    parser.add_argument("-l",
                        "--logs-dir",
                        type=str,
                        default=None,
                        help="Directory path where logs will be outputted. Default: No log")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_cluster.tsv',
                        help="Path to output file. Default: freddie_cluster.tsv")
    args = parser.parse_args()

    assert args.recycle_model in recycle_models
    assert args.gap_offset       >= 0
    assert args.epsilon          >= 0
    assert args.timeout          >  0
    assert args.threads          >  0
    assert args.min_isoform_size >= 0
    assert args.max_rounds       >= 0

    return args

def read_segment(segment_tsv):
    tints = dict()
    tint_prog   = re.compile(r'#%(chr_re)s\t%(cid_re)s\t%(positions_re)s\n$' % {
        'chr_re'       : '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
        'cid_re'       : '(?P<cid>[0-9]+)',
        'positions_re' : '(?P<positions>([0-9]+)(,([0-9]+))*)',
    })
    internal_gap_re = '(\d+)-(\d+):(\d+),'
    softclip_gap_re = '([ES]SC):(\d+),'
    poly_gap_re     = '([ES][AT])_(\d+):(\d+),'
    read_prog = re.compile(r'%(rid_re)s\t%(name_re)s\t%(chr_re)s\t%(strand_re)s\t%(cid_re)s\t%(data_re)s\t%(gaps)s\n$' % {
        'rid_re'    : '(?P<rid>[0-9]+)',
        'name_re'   : '(?P<name>[!-?A-~]{1,254})',
        'chr_re'    : '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
        'strand_re' : '(?P<strand>[+-])',
        'cid_re'    : '(?P<cid>[0-9]+)',
        'data_re'   : '(?P<data>[012]+)',
        'gaps'      : r'(?P<gaps>(%(i)s|%(s)s|%(p)s)*)' % {'i':internal_gap_re, 's':softclip_gap_re, 'p':poly_gap_re},
    })
    internal_gap_prog = re.compile(internal_gap_re)
    softclip_gap_prog = re.compile(softclip_gap_re)
    poly_gap_prog    = re.compile(poly_gap_re)
    for line in open(segment_tsv):
        if line[0]=='#':
            re_dict = tint_prog.match(line).groupdict()
            tint = dict(
                id    = int(re_dict['cid']),
                chr   = re_dict['chr'],
                segs  = [int(x) for x in re_dict['positions'].split(',')],
                reads = list(),
            )
            assert all(a<b for a,b in zip(tint['segs'][:-1],tint['segs'][1:])),(tint['segs'])
            tint['segs'] = [(s,e,e-s) for s,e in zip(tint['segs'][:-1],tint['segs'][1:])]
            assert not tint['id'] in tints, 'Transcriptional interval with id {} is repeated!'.format(tint['id'])
            tints[tint['id']] = tint
        else:
            re_dict = read_prog.match(line).groupdict()
            read = dict(
                id         = int(re_dict['rid']),
                name       = re_dict['name'],
                chr        = re_dict['chr'],
                strand     = re_dict['strand'],
                tint       = int(re_dict['cid']),
                data       = [int(d) for d in re_dict['data']],
                gaps       = {(int(g[0]),int(g[1])) : int(g[2])             for g in internal_gap_prog.findall(re_dict['gaps'])},
                softclip   = {s[0]                  : int(s[1])             for s in softclip_gap_prog.findall(re_dict['gaps'])},
                poly_tail  = {p[0]                  : (int(p[1]),int(p[2])) for p in poly_gap_prog.findall(    re_dict['gaps'])},
            )
            tind_id = read['tint']
            tints[tind_id]['reads'].append(read)
            assert len(read['data'])==len(tints[tind_id]['segs']), (read['data'],tints[tind_id]['segs'])
            assert read['chr']==tints[tind_id]['chr']
            assert all(0<=g[0]<g[1]<len(read['data']) for g in read['gaps'].keys())
    return tints

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

def partition_reads(tint):
    reads = tint['reads']
    I  = tint['ilp_data']['I']
    FL = tint['ilp_data']['FL']
    tint['partitions'] = list()

    rids = sorted(I.keys())
    unique_data = dict()
    edges = list()
    for i in rids:
        d = (tuple(I[i]),(FL[i][0], FL[i][1], reads[i]['poly_tail_category']))
        if d in unique_data:
            unique_data[d].append(i)
        else:
            unique_data[d] = [i]
    unique_data = list(unique_data.items())
    N = len(unique_data)
    for i in range(N):
        for j in range(i+1, N):
            d1,(f1,l1,t1) = unique_data[i][0]
            d2,(f2,l2,t2) = unique_data[j][0]
            f = max(f1,f2)
            l = min(l1,l2)
            o = l-f+1
            d = sum(x!=y      for x,y in zip(d1[f:l+1],d2[f:l+1]))
            w = sum(x==y==1 for x,y in zip(d1[f:l+1],d2[f:l+1]))
            if t1!='N' and t2!='N' and t1!=t2:
                continue
            if w < 1:
                continue
            if (o > 3 and d < 3) or (1<=o<=3 and d==0):
                edges.append((i,j))
    G = Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)
    while True:
        edges_to_remove = list()
        for i,j in G.edges:
            n1 = set(G.neighbors(i))
            n2 = set(G.neighbors(j))
            if len(n1)==1 or len(n2)==1 or len(n1 & n2) > 0:
                continue
            edges_to_remove.append((i,j))
        G.remove_edges_from(edges_to_remove)
        if len(edges_to_remove) == 0:
            break
    for c in components.connected_components(G):
        rids   = list()
        incomp = list()
        c = list(c)
        for idx,i in enumerate(c):
            rids.extend(unique_data[i][1])
            for j in c[idx+1:]:
                i,j = min(i,j),max(i,j)
                assert i<j
                if G.has_edge(i,j):
                    continue
                for rid_1 in unique_data[i][1]:
                    for rid_2 in unique_data[j][1]:
                        incomp.append((rid_1,rid_2))
        tint['partitions'].append((rids,incomp))

def preprocess_ilp(tint, ilp_settings):
    print('Preproessing ILP with {} reads and the following settings:\n{}'.format(len(tint['reads']), ilp_settings))
    reads = tint['reads']
    N = len(reads)
    M = len(tint['segs'])
    I = dict()
    C = dict()
    FL = dict()
    tail_s_rids = list()
    tail_e_rids = list()

    for i,read in enumerate(reads):
        I[i] = [0 for _ in range(M)]
        for j in range(0,M):
            I[i][j] = read['data'][j]%2

        C[i] = [0 for _ in range(M)]
        (min_i,max_i) = find_segment_read(I,i)
        read['poly_tail_category'] = 'N'
        if len(read['poly_tail'])==1:
            tail_key = next(iter(read['poly_tail']))
            tail_val = read['poly_tail'][tail_key]
            if tail_key in ['SA','ST'] and tail_val[0] > 10:
                tail_s_rids.append(i)
                read['poly_tail_category'] = 'S'
                read['gaps'][(-1,min_i)] = tail_val[1]
                min_i = 0
            elif tail_key in ['EA','ET'] and tail_val[0] > 10:
                read['poly_tail_category'] = 'E'
                tail_e_rids.append(i)
                read['gaps'][(max_i,M)] = tail_val[1]
                max_i = M-1
        FL[i] = (min_i,max_i)
        for j in range(0,M):
            if min_i <= j <= max_i and reads[i]['data'][j]==0:
                C[i][j]   = 1
            else:
                C[i][j]   = 0
    # Assigning a cost to the assignment of reads to the garbage isoform
    garbage_cost = {}
    for i in range(N):
        if ilp_settings['recycle_model'] == 'exons':
            garbage_cost[i] = garbage_cost_exons(I=I[i])
        elif ilp_settings['recycle_model'] == 'introns':
            garbage_cost[i] = garbage_cost_introns(C=C[i])
        elif ilp_settings['recycle_model'] == 'constant':
            garbage_cost[i] = 3
    tint['ilp_data']=dict(
        FL=FL,
        I=I,
        C=C,
        garbage_cost=garbage_cost,
    )

def informative_segs(tint, remaining_rids):
    M = len(tint['segs'])
    I = tint['ilp_data']['I']
    seg_content = [set() for _ in range(M)]
    informative = [True for _ in range(M)]
    for j in range(M):
        for i in remaining_rids:
            seg_content[j].add(I[i][j])
            if seg_content[j] == {0,1}:
                break
    for j in range(1,M-1):
        if len(seg_content[j])==1 and (seg_content[j-1]==seg_content[j]==seg_content[j+1]):
            informative[j]=False
    # for j in range(M):
    #     print(informative[j], seg_content[j])
    return informative

def run_ilp(tint, remaining_rids, incomp_rids, ilp_settings, log_prefix):
    # Variables directly based on the input ------------------------------------
    # I[i,j] = 1 if reads[i]['data'][j]==1 and 0 if reads[i]['data'][j]==0 or 2
    # C[i,j] = 1 if exon j is between the first and last exons (inclusively)
    #   covered by read i and is not in read i but can be turned into a 1
    ISOFORM_INDEX_START = 1
    M = len(tint['segs'])
    MAX_ISOFORM_LG = sum(seg[2] for seg in tint['segs'])
    I                 = tint['ilp_data']['I']
    C                 = tint['ilp_data']['C']
    INCOMP_READ_PAIRS = incomp_rids
    GARBAGE_COST      = tint['ilp_data']['garbage_cost']
    informative = informative_segs(tint, remaining_rids)

    # ILP model ------------------------------------------------------
    ILP_ISOFORMS = Model('isoforms_v7_20200608')
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
        if not informative[j]:
            continue
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
        for ((j1,j2),l) in tint['reads'][i]['gaps'].items():
            for k in range(ISOFORM_INDEX_START,ilp_settings['K']): # No such constraint on the garbage isoform if any
                if not (j1,j2,k) in GAPI:
                    assert informative[j1%M]
                    assert informative[j2%M]
                    assert not any(informative[j+1:j2])
                    GAPI[(j1,j2,k)]      = ILP_ISOFORMS.addVar(
                        vtype = GRB.INTEGER,
                        name  = 'GAPI[({j1},{j2},{k})]'.format(j1=j1,j2=j2,k=k)
                    )
                    GAPI_C1[(j1,j2,k)]   = ILP_ISOFORMS.addLConstr(
                        lhs   = GAPI[(j1,j2,k)],
                        sense = GRB.EQUAL,
                        rhs   = quicksum(E2I[j][k]*tint['segs'][j][2] for j in range(j1+1,j2) if informative[j]),
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
            if not informative[j]:
                continue
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
                    OBJ_SUM.addTerms(1.0, OBJ[i][j][k])
                    #     coeffs = 1.0,
                    #     vars   = OBJ[i][j][k]
                    # )
    # We add the chosen cost for each isoform assigned to the garbage isoform if any
    GAR_OBJ   = {}
    GAR_OBJ_C = {}
    for i in remaining_rids:
        if ilp_settings['recycle_model'] in ['constant', 'exons', 'introns']:
            OBJ_SUM.addTerms(1.0*GARBAGE_COST[i], R2I[i][0])
        elif ilp_settings['recycle_model'] == 'relative':
            GAR_OBJ[i]   = {}
            GAR_OBJ_C[i] = {}
            for j in range(0,M):
                if not informative[j]:
                    continue
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
                        OBJ_SUM.addTerms(1.0, GAR_OBJ[i][j][k])
                    elif I[i][j] == 0 and C[i][j] == 1:
                        pass

    ILP_ISOFORMS.setObjective(
        expr  = OBJ_SUM,
        sense = GRB.MINIMIZE
    )


    # Optimization
    #ILP_ISOFORMS.Params.PoolSearchMode=2
    #ILP_ISOFORMS.Params.PoolSolutions=5
    ILP_ISOFORMS.setParam('TuneOutput', 1)
    if not log_prefix == None:
        ILP_ISOFORMS.setParam('LogFile', '{}.glog'.format(log_prefix))
        ILP_ISOFORMS.write('{}.lp'.format(log_prefix))
    ILP_ISOFORMS.setParam('TimeLimit', ilp_settings['timeout']*60)
    ILP_ISOFORMS.optimize()

    ILP_ISOFORMS_STATUS = ILP_ISOFORMS.Status

    isoforms = {k:dict() for k in range(ISOFORM_INDEX_START,ilp_settings['K'])}
    print('STATUS: {}'.format(ILP_ISOFORMS_STATUS))
    if ILP_ISOFORMS_STATUS != GRB.Status.OPTIMAL:
        status = 'NO_SOLUTION'
    else:
        status = 'OPTIMAL'
        # Writing the optimal solution to disk
        if not log_prefix == None:
            solution_file = open('{}.sol'.format(log_prefix), 'w+')
            for v in ILP_ISOFORMS.getVars():
                solution_file.write('{}\t{}\n'.format(v.VarName, v.X))
            solution_file.close()
        # Isoform id to isoform structure
        for k in range(ISOFORM_INDEX_START,ilp_settings['K']):
            isoforms[k]['exons'] = list()
            for j in range(0,M):
                if informative[j]:
                    isoforms[k]['exons'].append(int(E2I[j][k].getAttr(GRB.Attr.X)>0.9))
                else:
                    isoforms[k]['exons'].append(I[next(iter(remaining_rids))][j])
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
                isoforms[k]['rid_to_corrections'][i] = [str(tint['reads'][i]['data'][j]) for j in range(M)]
                for j in range(0,M):
                    if not informative[j]:
                        isoforms[k]['rid_to_corrections'][i][j] = '-'
                    elif C[i][j] == 1 and OBJ[i][j][k].getAttr(GRB.Attr.X) > 0.9:
                        isoforms[k]['rid_to_corrections'][i][j] = 'X'
    return status,isoforms

def output_isoforms(tint, out_file):
    reads = tint['reads']
    output = list()
    output.append('#{}'.format(tint['chr']))
    output.append(str(tint['id']))
    output.append(','.join([str(s[0]) for s in tint['segs']]+[str(tint['segs'][-1][1])]))
    out_file.write('{}\n'.format('\t'.join(output)))
    for iid,isoform in enumerate(tint['isoforms']):
        output = list()
        output.append('isoform_{}'.format(iid))
        output.append(str(tint['id']))
        output.append(''.join(map(str, isoform['exons'])))
        out_file.write('{}\n'.format('\t'.join(output)))
        for i,corrections in isoform['rid_to_corrections'].items():
            output = list()
            output.append(str(reads[i]['id']))
            output.append(reads[i]['name'])
            output.append(reads[i]['chr'])
            output.append(reads[i]['strand'])
            output.append(str((reads[i]['tint'])))
            output.append(str((reads[i]['partition'])))
            output.append(str((reads[i]['poly_tail_category'])))
            output.append(str(iid))
            output.append(''.join(map(str,corrections)))
            exon_strs = [str(x) for x in corrections]
            for (j1,j2),l in reads[i]['gaps'].items():
                exon_strs[j1]+='({})'.format(l)
            output.extend(exon_strs)
            for k,v in sorted(reads[i]['poly_tail'].items()):
                output.append('{}:{}'.format(k,v))
            out_file.write('{}\n'.format('\t'.join(output)))
    for i in tint['garbage_rids']:
        output = list()
        output.append(str(reads[i]['id']))
        output.append(reads[i]['name'])
        output.append(reads[i]['chr'])
        output.append(reads[i]['strand'])
        output.append(str((reads[i]['tint'])))
        output.append(str((reads[i]['partition'])))
        output.append(str((reads[i]['poly_tail_category'])))
        output.append('*')
        output.append('*')
        exon_strs = [str(x) for x in reads[i]['data']]
        for (j1,j2),l in reads[i]['gaps'].items():
            exon_strs[j1]+='({})'.format(l)
        output.extend(exon_strs)
        for k,v in sorted(reads[i]['poly_tail'].items()):
            output.append('{}:{}'.format(k,v))
        out_file.write('{}\n'.format('\t'.join(output)))

def cluster_tint(cluster_args):
    tint, ilp_settings, min_isoform_size, logs_dir=cluster_args
    print('# Clustering tint {}'.format(tint['id']))
    if not logs_dir == None:
        os.makedirs('{}/{}'.format(logs_dir, tint['id']), exist_ok=True)
    preprocess_ilp(tint, ilp_settings)
    partition_reads(tint)
    print('# Paritions ({}) sizes: {}\n'.format(len(tint['partitions']), [len(p) for p in tint['partitions']]))
    tint['isoforms'] = list()
    tint['garbage_rids'] = list()
    for partition,(remaining_rids,incomp_rids) in enumerate(tint['partitions']):
        for rid in remaining_rids:
            tint['reads'][rid]['partition']=partition
        print('==========\ntint {}: Running {}-th partition...'.format(tint['id'], partition))
        for round in range(ilp_settings['max_rounds']):
            print('==========\ntint {}: Running {}-th round with {} reads...'.format(tint['id'], round, len(remaining_rids)))
            if len(remaining_rids) < min_isoform_size:
                break
            status,round_isoforms = run_ilp(
                tint           = tint,
                remaining_rids = remaining_rids,
                incomp_rids    = incomp_rids,
                ilp_settings   = ilp_settings,
                log_prefix     = '{}/{}/partition.{}.round.{}'.format(logs_dir, tint['id'], partition, round) if logs_dir != None else None,
            )
            if status != 'OPTIMAL' or max([0]+[len(i['rid_to_corrections']) for i in round_isoforms.values()]) < min_isoform_size:
                break
            for k,isoform in round_isoforms.items():
                print('Isoform {} size: {}'.format(k, len(isoform['rid_to_corrections'])))
                if len(isoform['rid_to_corrections']) < min_isoform_size:
                    continue
                tint['isoforms'].append(isoform)
                for rid,corrections in isoform['rid_to_corrections'].items():
                    assert rid in remaining_rids
                    remaining_rids.remove(rid)
                    tint['reads'][rid]['corrections']=corrections
                    tint['reads'][rid]['isoform']=len(tint['isoforms'])-1
            print('------->')
            print('Remaining reads: {}\n'.format(len(remaining_rids)))
            print('<-------')
        tint['garbage_rids'].extend(sorted(remaining_rids))
    for rid in tint['garbage_rids']:
        tint['reads'][rid]['isoform']     = None
        tint['reads'][rid]['corrections'] = None
    return tint

def main():
    args = parse_args()

    ilp_settings=dict(
        recycle_model    = args.recycle_model,
        K                = 2,
        epsilon          = args.epsilon,
        offset           = args.gap_offset,
        timeout          = args.timeout,
        max_rounds       = args.max_rounds,
        threads          = 1,
    )
    tints = read_segment(segment_tsv=args.segment_tsv)
    cluster_args = [
        (tint, ilp_settings, args.min_isoform_size, args.logs_dir) for tint in tints.values() if tint['id']==213
    ]
    print(ilp_settings)
    out_file = open(args.output, 'w')
    with Pool(args.threads) as p:
        if args.threads == 1:
            p.close()
        for idx,tint in enumerate(p.imap_unordered(cluster_tint, cluster_args, chunksize=5)) if args.threads>1 else enumerate(map(cluster_tint, cluster_args)):
            output_isoforms(tint, out_file)
            print('Done with {}-th transcriptional multi-intervals ({}/{})'.format(tint['id'], idx+1,len(tints)))
    out_file.close()

if __name__ == "__main__":
    main()
