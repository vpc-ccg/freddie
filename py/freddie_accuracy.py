#!/usr/bin/env python3
import os
import argparse
import glob

import numpy as np
import networkx as nx
from networkx.algorithms import clique
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
# import plotly


BIP_TRUTH = 0
BIP_TOOL = 1


tool_ls = {
    'isoform' : 'solid',
    'genome' : 'dashed',
    'segment' : 'dotted',
    'stringtie' : 'solid',
    'freddie'   : 'solid',
    'flair.0.01' : 'solid',
    'flair.0.25' : 'solid',
    'flair.0.50' : 'solid',
    'flair.0.75' : 'solid',
    'flair.1.00' : 'solid',
}


tool_colors = {
    'isoform' : 'black',#'#cbc9e2',
    'genome' : 'black',#'#9e9ac8',
    'segment' : 'black',#'#6a51a3',
    'stringtie'  : '#d95f02',
    'freddie'    : '#7570b3',
    'flair.0.01' : '#ccece6',
    'flair.0.25' : '#99d8c9',
    'flair.0.50' : '#66c2a4',
    'flair.0.75' : '#2ca25f',
    'flair.1.00' : '#006d2c',
}
names = {
    'isoform' : 'Perfect clustering (isoform alignment)',
    'genome' : 'Perfect clustering (genome alignment)',
    'segment' : 'Perfect clustering (Freddie segmentation)',
    'stringtie' : 'StringTie2',
    'freddie'   : 'Freddie',
    'flair.0.01' : 'FLAIR (1%)',
    'flair.0.25' : 'FLAIR (25%)',
    'flair.0.50' : 'FLAIR (50%)',
    'flair.0.75' : 'FLAIR (75%)',
    'flair.1.00' : 'FLAIR (100%)',
}


def parse_args():
    default_thresholds = list(np.arange(0.90, .99, 0.01))
    default_thresholds = [round(t, 2) for t in default_thresholds]
    parser = argparse.ArgumentParser(
        description="Outputs accuracy plots for different tools against the simulated grount truth")
    parser.add_argument("-id",
                        "--input-dir",
                        type=str,
                        required=True,
                        help="Input directory for edges TSV and BED files")
    parser.add_argument("-od",
                        "--out-dir",
                        type=str,
                        default='freddie_accuracy/',
                        help="Path to output plots. Default: freddie_accuracy/")
    parser.add_argument("-t",
                        "--thresholds",
                        nargs="+",
                        type=float,
                        default=default_thresholds,
                        help="Space separated list of alignment similarity thresholds. Default: {}".format(
                            ' '.join(map(str, default_thresholds)))
                        )
    args = parser.parse_args()
    return args


def read_nodes_and_edges(in_dir):
    tool_nodes = dict()
    tool_edges = dict()
    for bed in glob.iglob('{}/*.bed'.format(in_dir)):
        tool = bed.split('/')[-1][:-len('.bed')]
        if tool == 'troth':
            tool = 'truth'
        tool_nodes[tool] = set()
        for l in open(bed):
            n = l.rstrip().split('\t')[3]
            t = n.split('_')[0]
            assert t == tool or (t == 'troth' and tool == 'truth')
            tool_nodes[tool].add(n)
    assert 'truth' in tool_nodes
    for edge_tsv in glob.iglob('{}/*.edges.tsv'.format(in_dir)):
        tool = edge_tsv.split('/')[-1][:-len('.edges.tsv')]
        if tool == 'troth':
            tool = 'truth'
        tool_edges[tool] = dict()
        for l in open(edge_tsv):
            n1, n2, w = l.rstrip().split('\t')
            t1 = n1.split('_')[0]
            if t1 == 'troth':
                t1 = 'truth'
            t2 = n2.split('_')[0]
            if t2 == 'troth':
                t2 = 'truth'

            assert t1 in [tool, 'truth'], (t1,tool)
            assert t2 == tool, (t2, tool)
            assert n1 in tool_nodes[tool] or tool_nodes['truth']
            assert n2 in tool_nodes[tool]
            k = (n1, n2)
            assert not k in tool_edges[tool], (k, tool, w, tool_edges[tool])
            tool_edges[tool][k] = float(w)
        # print(tool[:7], len(tool_edges[tool]), len(tool_nodes[tool]), sep='\t' )
    for tool in tool_edges:
        if tool == 'truth':
            continue
        for k, v in tool_edges['truth'].items():
            tool_edges[tool][k] = v
    tool_edges.pop('truth')
    return tool_nodes, tool_edges


def build_graph(X_nodes, Y_nodes, edges, threshold):
    G = nx.Graph()
    G.add_nodes_from(X_nodes, bipartite=BIP_TRUTH)
    G.add_nodes_from(Y_nodes, bipartite=BIP_TOOL)
    G.add_edges_from([(u, v) for (u, v), w in edges.items() if w >= threshold])
    return G


def get_stats(G):
    comp_stats = dict(
        count=0,
        TP=0,
        FN=0,
        FP=0,
        LTP=0,
        LFX=0,
    )
    isof_stats = dict(
        truth_count=0,
        tool_count=0,
        TP_truth=0,
        TP_tool=0,
        FP=0,
        FN=0,
        LTP_truth=0,
        LTP_tool=0,
        LFP=0,
        LFN=0,
    )
    for component in nx.connected_components(G):
        Gc = G.subgraph(component)
        is_clique = (Gc.size() == len(component)*(len(component)-1)/2)
        N = sum(1 for _, bip in Gc.nodes(data='bipartite')
                if bip == BIP_TRUTH)  # number of truth isoforms
        M = sum(1 for _, bip in Gc.nodes(data='bipartite')
                if bip == BIP_TOOL)   # number of tool isoforms
        assert N + M == len(component)

        kind = None
        if is_clique and N > 0 and M > 0:
            kind = 'SM'  # clique with both tool and truth isoforms
        elif N == 0 or M == 0:
            kind = 'SX'  # only truth or only tool isoforms
        else:
            matched_req = {
                BIP_TRUTH: M,  # Each tool isoform node must connect all N truth nodes
                BIP_TOOL: N,  # Each truth isoform node must connect all M tool nodes
            }
            match_flag = all(  # check if all nodes have full oppositional matching
                sum(  # count the number of outgoing edges to nodes of the opposite type of v
                    1 for (_, u) in Gc.edges(v) if Gc.nodes[u]['bipartite'] != v_bip
                ) == matched_req[v_bip] for v, v_bip in Gc.nodes(data='bipartite')
            )
            if match_flag:
                kind = 'WM'
            else:
                kind = 'WX'
        comp_stats['count'] += 1
        isof_stats['truth_count'] += N
        isof_stats['tool_count'] += M
        if kind == 'SM':
            comp_stats['TP'] += 1
            isof_stats['TP_truth'] += N
            isof_stats['TP_tool'] += M
        elif kind == 'WM':
            comp_stats['LTP'] += 1
            isof_stats['LTP_truth'] += N
            isof_stats['LTP_tool'] += M
        elif kind == 'WX':
            comp_stats['LFX'] += 1
            isof_stats['LFN'] += N
            isof_stats['LFP'] += M
        elif kind == 'SX':
            comp_stats['FN'] += N > 0
            comp_stats['FP'] += M > 0
            isof_stats['FN'] += N
            isof_stats['FP'] += M

    return comp_stats, isof_stats


def plot_components(stats, tools, tool_labels, threshold, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(25, 10), sharex=True)
    fig.suptitle('Connected components at alignment threshold={}'.format(0.97))
    group_labels = [
        ('count', 'All CCs'),
        ('TP', 'True positive\n(cliques with GTIs and PIs)'),
        ('FN', 'False positive\n(CCs with only PIs)'),
        ('FP', 'False negative\n(CCs with only GTIs)'),
        ('LTP', 'Likely true positive\n(CCs with PIs matching all GTIs)'),
        ('LFX', 'Likely false positive/negative\n(CCs with PIs not matching all GTIs)'),
    ]
    x = np.arange(len(group_labels))  # the label locations
    width = 1/(2+len(tools))  # the width of the bars
    axes[0].set_xticks(x+width*(len(tools)-1)/2)
    axes[0].set_xticklabels([l for _, l in group_labels])
    for ax, kind in zip(axes, ['abso', 'norm']):
        if kind == 'abso':
            ax.set_ylabel('# of connected components (CCs)')
            ax.set_ylim(0, max(stats[t][threshold]['count']
                               for t in tools)*1.2)
        else:
            ax.set_ylim(0.0, 1.2)
            ax.yaxis.set_major_formatter(mticker.PercentFormatter())
            ax.set_ylabel(
                '%% of each tool\'s graph connected components (CCs)')
        for idx, tool in enumerate(tools):
            if kind == 'abso':
                data = [stats[tool][threshold][k] for k, _ in group_labels]
            else:
                data = [stats[tool][threshold][k]/stats[tool]
                        [threshold]['count'] for k, _ in group_labels]
            rects = ax.bar(x + idx*width, data, width, label=tool_labels[tool])
            for rect in rects:
                height = rect.get_height()
                ax.annotate('{}'.format(height) if kind == 'abso' else '{:.1%}'.format(height),
                            rotation=60,
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            color=tool_colors[tool],
                            ha='center', va='bottom')
    fig.legend(axes[-1].get_legend_handles_labels()[1], loc='right')
    fig.savefig('{}/comp_{}.pdf'.format(out_dir, threshold))


def sort_tool_names(tools):
    sorted_tools = list()
    tool_to_name = dict()
    if 'isoform' in tools:
        sorted_tools.append('isoform')
        tool_to_name['isoform'] = 'Perfect clustering (isoform alignment)'
    if 'genome' in tools:
        sorted_tools.append('genome')
        tool_to_name['genome'] = 'Perfect clustering (genome alignment)'
    if 'segment' in tools:
        sorted_tools.append('segment')
        tool_to_name['segment'] = 'Perfect clustering (freddie segmentation)'
    if 'freddie' in tools:
        sorted_tools.append('freddie')
        tool_to_name['freddie'] = 'Freddie'
    if 'stringtie' in tools:
        sorted_tools.append('stringtie')
        tool_to_name['stringtie'] = 'StringTie2'
    for t in sorted(tools, reverse=True):
        if t.startswith('flair'):
            sorted_tools.append(t)
            tool_to_name[t] = 'FLAIR ({:.0%})'.format(float(t[6:]))
    for t in sorted(tools, reverse=True):
        if not t in sorted_tools:
            sorted_tools.append(t)
    return sorted_tools, tool_to_name


def plot_isoforms(stats, thresholds, tools, tool_labels, prefix):
    fig, axes = plt.subplots(2, 1, figsize=(25, 10), sharex=True)
    fig.suptitle('Connected components at alignment threshold={}'.format(0.97))
    plot_titles = [
        ('TP_truth', 'True positive GTIs'),
        ('TP_tool', 'True positive PIs'),
        ('TP_tool%', 'True positive PIs (%)'),
        ('FN', 'False negative GTIs'),
        ('FP', 'False positive PIs'),
        ('FP%', 'False positive PIs (%)'),
        ('LTP_truth', 'Likely true positive GTIs\n(matching all PIs)'),
        ('LTP_tool', 'Likely true positive PIs\n(matching all GTIs)'),
        ('LTP_tool%', 'Likely true positive PIs\n(matching all GTIs, %)'),
        ('LFN', 'Likely false negative GTIs\n(not matching all PIs)'),
        ('LFP', 'Likely false positive PIs\n(not matching all GTIs)'),
        ('LFP%', 'Likely false positive PIs\n(not matching all GTIs, %)'),
    ]
    fig, axes = plt.subplots(
        nrows=4,
        ncols=3,
        figsize=(28, 19),
        sharex=True,
        gridspec_kw={'hspace': 0.3, 'wspace': 0.2}
    )

    for ax, (f, title) in zip(axes.flat, plot_titles):
        print(title)
        ax.set_xticks(thresholds)
        ax.set_xticklabels(['{:2.2f}'.format(t) for t in thresholds])
        if f.endswith('_truth') or f.endswith('FN') or f.endswith('%'):
            ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter())
        if f.endswith('_%'):
            ax.yaxis.tick_right()
        ax.tick_params(axis='x', labelrotation=45)
        ax.set_title(title)
        ax.grid(True)
        for tool in tools:
            if f.endswith('%'):
                data = [
                    100*stats[tool][t][f[:-1]]/stats[tool][t]['tool_count']
                    for t in thresholds
                ]
            elif f.endswith('_truth') or f.endswith('FN'):
                data = [
                    100*stats[tool][t][f]/stats[tool][t]['truth_count']
                    for t in thresholds
                ]
            else:
                data = [
                    stats[tool][t][f]
                    for t in thresholds
                ]
            print(tool, f, data[-2])
            ax.plot(thresholds,
                    data,
                    label=tool_labels[tool],
                    ls=tool_ls[tool],
                    lw=3,
                    color=tool_colors[tool],
                    )
        if f.endswith('_truth'):
            axsec = ax.secondary_yaxis(location='right', functions=(
                lambda x: x/100*stats[tools[0]][thresholds[0]]['truth_count'], lambda x: x*100/stats[tools[0]][thresholds[0]]['truth_count']))
            axsec.xaxis.set_major_formatter(mticker.PercentFormatter())

    plt.subplots_adjust(bottom=0.15)
    axes_handles, axes_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        axes_handles,
        axes_labels,
        ncol=len(tools),
        loc='lower center',
         bbox_to_anchor=(0.5, 0.075),
    )

    for ax in axes[-1, :]:
        ax.set_xlabel('Alignment score threshold', fontsize=20)
    fig.savefig('{}isoforms.pdf'.format(prefix))



def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=False)

    nodes, edges = read_nodes_and_edges(args.input_dir)
    tools, tool_labels = sort_tool_names(edges.keys())
    comp_stats = dict()
    isof_stats = dict()
    graphs = list()
    import pickle
    try:
        print('Trying to read from pickle')
        (comp_stats,isof_stats) = pickle.load(open(f'{args.out_dir}.pickle','rb'))
    except:
        for tool in tools:
            print('Building graphs for {}'.format(tool))
            comp_stats[tool] = dict()
            isof_stats[tool] = dict()
            for threshold in args.thresholds:
                G = build_graph(
                    X_nodes=nodes['truth'],
                    Y_nodes=nodes[tool],
                    edges=edges[tool],
                    threshold=threshold,
                )
                cur_comp_stats, cur_isof_stats = get_stats(G)
                comp_stats[tool][threshold] = cur_comp_stats
                isof_stats[tool][threshold] = cur_isof_stats
        pickle.dump((comp_stats,isof_stats), open(f'{args.out_dir}.pickle','wb+'))
    # print('Plotting components')
    # for threshold in args.thresholds:
    #     plot_components(
    #         stats=comp_stats,
    #         tools=tools,
    #         tool_labels=tool_labels,
    #         threshold=threshold,
    #         out_dir=args.out_dir,
    #     )
    print('Plotting isoforms')
    plot_isoforms(
        stats=isof_stats,
        thresholds=args.thresholds,
        tools=[t for t in tools if t.startswith('flair')],
        tool_labels=tool_labels,
        prefix=f'{args.out_dir}/all_flair.',
    )
    plot_isoforms(
        stats=isof_stats,
        thresholds=args.thresholds,
        tools=[t for t in tools if not t in {'flair.0.75','flair.0.25','flair.0.01', 'genome', 'segment'}],
        tool_labels=tool_labels,
        prefix=f'{args.out_dir}/all_tools.',
    )


if __name__ == "__main__":
    main()