#include "dag_align.h"

#include <fstream>      // std::ofstream
#include <sstream>
#include <stack>

#define GENE_ID 0
#define MATCH_S 2
#define GAP_S -1
#define MISMATCH_S -1

using namespace std;

void process_gene_test() {
    sequence_t gene_test = "CATGTC";
    sequence_list_t reads_test;

    process_gene(gene_test, reads_test);
}

void print_matrix(const align_matrix_t& D, const backtrack_matrix_t& B){
    for (size_t i = 0; i < D.size(); i++) {
        for (size_t j = 0; j < D[i].size(); j++) {
            cout << D[i][j] << " (" << B[i][j].first << "," <<B[i][j].second << ")\t";
        }
        cout << endl;
    }
}

void process_gene(sequence_t gene, sequence_list_t reads) {
    in_neighbors_t parents;
    out_neighbors_t children;
    node_to_base_t node_to_base;
    node_to_read_t node_to_read;

    // cout << gene << endl;
    node_id_t last_node = add_node(parents, children, node_to_base, node_to_read, gene[0], GENE_ID);
    node_id_t current_node = last_node;
    for (size_t i = 1; i < gene.size(); i++) {
        last_node = current_node;
        current_node = add_node(parents, children, node_to_base, node_to_read, gene[i], GENE_ID);
        add_edge(parents, children, last_node, current_node);
    }

    add_edge(parents, children, 1, 5);
    current_node = add_node(parents, children, node_to_base, node_to_read, 'X', GENE_ID);
    node_to_read[current_node].push_back(1);
    add_edge(parents, children, current_node, 0);

    std::vector<node_list_t> topo_sorted = topological_sort(parents.size(), parents, children);
    // print_graph(topo_sorted, parents, children, node_to_base, node_to_read);
    generate_dot(topo_sorted, children, node_to_base, node_to_read, "test2.dot");
    align_matrix_t D;
    backtrack_matrix_t B;
    sequence_t s = "CAGT";
    // dag_local_alignment(s, topo_sorted, parents, node_to_base, D, B);
    print_matrix(D, B);
    // for (size_t i = 0; i < D.size(); i++) {
    //     for (size_t j = 0; j < D[0].size(); j++) {
    //         cout << "(" << j<< "):" << D[j] << "\t";
    //     }
    //     cout << endl;
    // }
}

void dag_local_alignment(const sequence_t seq,
                         const node_list_t& topo_sorted,
                         const in_neighbors_t& parents,
                         const node_to_base_t& node_to_base,
                         align_matrix_t& D,
                         backtrack_matrix_t& B) {
    size_t dag_size = topo_sorted.size();
    size_t seq_size = seq.size();
    vector<size_t> node_to_topo_idx(dag_size);
    for (size_t topo_idx = 0; topo_idx < dag_size; topo_idx++) {
        // cout << topo_sorted[topo_idx] << " : " << topo_idx << endl;
        node_to_topo_idx[topo_sorted[topo_idx]] = topo_idx;
    }
    // cout << endl;
    D.clear();
    B.clear();
    D = align_matrix_t(seq_size + 1, align_row_t(dag_size+1, 0));
    B = backtrack_matrix_t(seq_size + 1, backtrack_row_t(dag_size + 1, backtrack_t(-1,-1)));

    for (size_t i = 1; i < seq_size+1; i++) {
        for (size_t j = 1; j < dag_size+1; j++) {
            size_t seq_pos = i-1;
            node_id_t node = topo_sorted[j-1];
            // cout << i << " : " << j << " - " << seq_pos << " : " << node << endl;
            if (parents[node].size() != 1) {
                if (node_to_base[node]==seq[seq_pos]) {
                    D[i][j] = MATCH_S;
                }
                // cout << '\t' << parents[node].size() << " : " << match << endl;
                continue;
            }

            size_t j_parent = node_to_topo_idx[parents[node][0]];
            align_score_t match = D[i-1][j_parent];
            if (node_to_base[node]==seq[seq_pos]) {
                match += MATCH_S;
            } else {
                match += MISMATCH_S;
            }
            align_score_t dag_skip = D[i][j_parent] + GAP_S;
            align_score_t read_skip = D[i-1][j] + GAP_S;

            if (read_skip > dag_skip && read_skip > match) {
                D[i][j] = read_skip;
                B[i][j] = backtrack_t(i-1, j);
            } else if (dag_skip > match) {
                D[i][j] = dag_skip;
                B[i][j] = backtrack_t(i, j_parent);
            } else {
                D[i][j] = match;
                B[i][j] = backtrack_t(i-1, j_parent);
            }
        }
    }
}

void add_edge(in_neighbors_t &in_neighbors, out_neighbors_t &out_neighbors, node_id_t source, node_id_t target) {
    out_neighbors[source].push_back(target);
    in_neighbors[target].push_back(source);
}

node_id_t add_node(in_neighbors_t &in_neighbors, out_neighbors_t &out_neighbors, node_to_base_t &node_to_base,
                   node_to_read_t &node_to_read, nucleotide_t base, read_id_t read_id) {
    in_neighbors.push_back(node_list_t());
    out_neighbors.push_back(node_list_t());
    node_to_read.push_back(read_list_t());
    node_to_read[node_to_read.size()-1].push_back(read_id);
    node_to_base.push_back(base);
    return in_neighbors.size() - 1;
}

void print_graph(const node_list_t node_list, const in_neighbors_t &in_neighbors, const out_neighbors_t &out_neighbors,
                 const node_to_base_t &node_to_base, const node_to_read_t &node_to_read)  {
    cout << "N [" << node_list.size() << "]" << endl;
    for (node_id_t node : node_list) {
        cout << node << " (" << node_to_base[node] << ")" << endl;
        cout << "\tIN ["<< in_neighbors[node].size() <<"]: ";
        for (node_id_t neighbor : in_neighbors[node]) {
            cout << neighbor << " (" << node_to_base[neighbor] << "), ";
        }
        cout << endl;
        cout << "\tOUT ["<< out_neighbors[node].size() <<"]: ";
        for (node_id_t neighbor : out_neighbors[node]) {
            cout << neighbor << " (" << node_to_base[neighbor] << "), ";
        }
        cout << endl;
    }
}

vector<node_list_t> topological_sort(const node_id_t node_count, const in_neighbors_t& parents,
                                     const out_neighbors_t& children) {
    vector<bool> visited(node_count, false);
    vector<node_list_t> result;
    for (node_id_t node = 0; node < node_count; node++) {
        if (visited[node] == true) {
            continue;
        }
        stack<node_id_t> opened;
        opened.push(node);
        while (opened.size() > 0) {
            node_id_t current_node = opened.top();
            if (visited[current_node] == true) {
                opened.pop();
                continue;
            }
            size_t unvisited_parents = 0;
            for (node_id_t parent : parents[current_node]) {
                if (visited[parent] == false) {
                    unvisited_parents++;
                    opened.push(parent);
                }
            }
            if (unvisited_parents == 0) {
                node_list_t segment;
                segment.push_back(current_node);
                visited[current_node] = true;
                opened.pop();
                while (true) {
                    if (children[current_node].size() != 1) {
                        break;
                    }
                    node_id_t child = children[current_node][0];
                    if (parents[child].size() != 1) {
                        break;
                    }
                    visited[child] = true;
                    segment.push_back(child);
                    current_node = child;
                }
                result.push_back(segment);
            }
        }
    }
    return result;
}

void generate_dot(const vector<node_list_t>& topo_sorted, const out_neighbors_t& children,
                  const node_to_base_t& node_to_base, const node_to_read_t& node_to_read, const string output_path)  {
    stringstream output;
    output << "digraph graphname{" << endl;
    for (size_t seg_idx = 0; seg_idx < topo_sorted.size(); seg_idx++) {
        for (node_id_t node : topo_sorted[seg_idx]) {
            output << "    " << node <<" [label=\"" << seg_idx << ":" << node << ":" << node_to_base[node] << ":" <<
                node_to_read[node].size() << "\"]" << endl;
        }
    }
    output << endl;
    for (node_list_t segment : topo_sorted) {
        for (node_id_t node : segment) {
            for (node_id_t child : children[node]) {
                output << "    " << node <<"->" << child << endl;
            }
        }
    }
    output << "}";
    ofstream ofile;
    ofile.open(output_path);
    ofile << output.str();
    ofile.close();
}

//
