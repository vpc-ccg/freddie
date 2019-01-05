#include "dag_align.h"

#include <fstream>  // std::ofstream
#include <sstream>
#include <queue>
#include <algorithm> // std::reverse
#include <stdlib.h> //abort()

#define EXONIC_ID 0x7FFFFFFF
#define INTRONIC_ID 0x7FFFFFFE
#define MATCH_S 1
#define GAP_S -1
#define MISMATCH_S -1

constexpr matrix_coordinate_t INVALID_COORDINATE = {-1,-1};

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ofstream;
using std::queue;
using std::reverse;

void process_gene_test() {
    sequence_list_t reads(1, "ACGGCGTATTGCACGT");
    sequence_t gene = "GTACGCAC";
    exonic_indicator_t exonic(gene.size(), true);
    process_gene(reads, gene, exonic);
}

void process_gene(const sequence_list_t& reads,
                  const sequence_t& gene,
                  const exonic_indicator_t& exonic) {
    in_neighbors_t parents;
    out_neighbors_t children;
    node_to_reads_t node_to_read;

    for (size_t i = 0; i < gene.size(); i++) {
        append_node(parents, children, node_to_read);
    }

    generate_dot(children, gene, exonic, node_to_read, "gene.dot");
    for (size_t i = 0; i < reads.size(); i++) {
        align_matrix_t D;
        backtrack_matrix_t B;
        local_aligner_func_t local_aligner = get_local_aligner(D, B, reads[i], gene, exonic, parents);
        local_alignment(D, B, reads[i], gene, local_aligner);
        print_matrix(reads[i], gene, D, B);
        align_path_t opt_alignment;
        align_score_t opt_score;
        extract_local_alignment(opt_alignment, opt_score, D, B);
        recalc_alignment_matrix(D, B, opt_alignment, children, local_aligner);

        print_matrix(reads[i], gene, D, B);
    }

}

align_score_t match(char i, char j, bool is_exonic) {
    if (i == j) {
        return MATCH_S + MATCH_S*is_exonic;
    } else {
        return MISMATCH_S;
    }
}

local_aligner_func_t get_local_aligner(align_matrix_t& D,
                                       backtrack_matrix_t& B,
                                       const sequence_t& read,
                                       const sequence_t& gene,
                                       const exonic_indicator_t& exonic,
                                       const out_neighbors_t& parents) {
    auto set_to_max_match = [&D, &read, &gene, &exonic] (align_score_t& opt_s, matrix_coordinate_t& opt_b, const backtrack_t& row, const backtrack_t& col) {
        align_score_t cur_s = D[row][col] + match(read[row], gene[col], exonic[col]);
        if (cur_s > opt_s) {
            opt_s = cur_s;
            opt_b = {row, col};
        }
    };
    auto set_to_max_gap = [&D] (align_score_t& opt_s, matrix_coordinate_t& opt_b, const backtrack_t& row, const backtrack_t& col) {
        align_score_t cur_s = D[row][col] + GAP_S;
        if (cur_s > opt_s) {
            opt_s = cur_s;
            opt_b = {row, col};
        }
    };
    auto local_aligner = [&D, &B, &set_to_max_match, &set_to_max_gap, &parents] (const backtrack_t& i, const backtrack_t& j) {
        align_score_t opt_s = 0;
        matrix_coordinate_t opt_b = INVALID_COORDINATE;
        // Three direct parents parents
        set_to_max_match(opt_s, opt_b, i-1, j-1); // Match
        set_to_max_gap(opt_s, opt_b, i-1, j-0); // Delete
        set_to_max_gap(opt_s, opt_b, i-0, j-1); // Insert
        // Other DAG parents
        for (node_id_t parent : parents[j-1]) {
            size_t parent_j = parent + 1;
            // Three indirect parents
            set_to_max_match(opt_s, opt_b, i-1, parent_j-1); // Match
            set_to_max_gap(opt_s, opt_b, i-1, parent_j-0); // Delete
            set_to_max_gap(opt_s, opt_b, i-0, parent_j-1); // Insert
        }
        D[i][j] = opt_s;
        B[i][j] = opt_b;
    };
    return local_aligner;
}

void local_alignment(align_matrix_t& D,
                     backtrack_matrix_t& B,
                     const sequence_t& read,
                     const sequence_t& gene,
                     const local_aligner_func_t& local_aligner) {
    cout << read << endl;
    cout << gene << endl;
    size_t read_l = read.size();
    size_t gene_l = gene.size();
    D.clear();
    B.clear();
    D = align_matrix_t(read_l + 1, align_row_t(gene_l + 1, 0));
    B = backtrack_matrix_t(read_l + 1, backtrack_row_t(gene_l + 1, INVALID_COORDINATE));

    for (size_t i = 1; i <= read_l; i++) {
        for (size_t j = 1; j <= gene_l; j++) {
            local_aligner(i, j);
        }
    }
}

void recalc_alignment_matrix(align_matrix_t& D,
                             backtrack_matrix_t& B,
                             const align_path_t& opt_alignment,
                             const in_neighbors_t& children,
                             const local_aligner_func_t& local_aligner) {
    queue<matrix_coordinate_t> clearing_queue;
    auto queue_children = [&D, &B, &children, &clearing_queue] (const matrix_coordinate_t& pos) {
        auto queue_if_child = [&D, &B, &children, &clearing_queue, &pos] (const matrix_coordinate_t& descendant) {
            if (descendant.first >= D.size()) {
                return;
            }
            if (descendant.second >= D[descendant.first].size()) {
                return;
            }
            if (B[descendant.first][descendant.second] != pos) {
                return;
            }
            clearing_queue.push(descendant);
        };
        matrix_coordinate_t descendant;
        descendant = {pos.first + 1, pos.second + 1};
        queue_if_child(descendant);
        descendant = {pos.first + 0, pos.second + 1};
        queue_if_child(descendant);
        descendant = {pos.first + 1, pos.second + 0};
        queue_if_child(descendant);
        for (const node_id_t& child : children[pos.second - 1]) {
            descendant = {pos.first + 1, child + 1};
            queue_if_child(descendant);
            descendant = {pos.first + 0, child + 1};
            queue_if_child(descendant);
            descendant = {pos.first + 1, child + 0};
            queue_if_child(descendant);
        }
    };

    for (matrix_coordinate_t pos : opt_alignment) {
        D[pos.first][pos.second] = 0;
        B[pos.first][pos.second] = INVALID_COORDINATE;
        clearing_queue.push(pos);
    }
    print_matrix("ACGGCGTATTGCACGT", "GTACGCAC", D, B);

    while (clearing_queue.size() > 0) {
        matrix_coordinate_t pos = clearing_queue.front();
        queue_children(pos);
        if (B[pos.first][pos.second] != INVALID_COORDINATE) {
            cout << "(" << pos.first << "," << pos.second << ") ";
            cout << D[pos.first][pos.second] << " --> ";
            local_aligner(pos.first, pos.second);
            cout << D[pos.first][pos.second] << endl;
        }
        clearing_queue.pop();
    }
}

void extract_local_alignment(align_path_t& opt_alignment,
                             align_score_t& opt_score,
                             const align_matrix_t& D,
                             const backtrack_matrix_t& B) {
    opt_alignment.clear();
    // Find opt score
    matrix_coordinate_t opt_tail;
    opt_score = 0;
    for (size_t i = 0; i < D.size(); i++) {
        for (size_t j = 0; j < D[0].size(); j++) {
            if (D[i][j] > opt_score) {
                opt_score = D[i][j];
                opt_tail = {i,j};
            }
        }
    }
    opt_alignment.push_back(opt_tail);
    // Backtrack
    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B[cur_pos.first][cur_pos.second];
    while (D[nxt_pos.first][nxt_pos.second] > 0) {
        opt_alignment.push_back(nxt_pos);
        cur_pos = nxt_pos;
        nxt_pos = B[cur_pos.first][cur_pos.second];
    }
    reverse(opt_alignment.begin(), opt_alignment.end());
}

void add_edge(in_neighbors_t &in_neighbors,
              out_neighbors_t &out_neighbors,
              node_id_t source,
              node_id_t target) {
    if (source >= target) {
        cout << "Source can't be >= target: " << source << " -> " << target << endl;
        abort();
    }
    if (source <= 0) {
        cout << "Source can't be <= : " << source << endl;
        abort();
    }
    if (target >= in_neighbors.size()) {
        cout << "Target can't be >= in_neighbors.size(): " << target << " >= " << in_neighbors.size() << endl;
        abort();
    }
    out_neighbors[source].push_back(target);
    in_neighbors[target].push_back(source);
}

node_id_t append_node(in_neighbors_t &in_neighbors,
                      out_neighbors_t &out_neighbors,
                      node_to_reads_t &node_to_read) {
    node_id_t new_node = in_neighbors.size();
    in_neighbors.push_back(node_list_t());
    out_neighbors.push_back(node_list_t());
    node_to_read.push_back(read_list_t());
    // add_edge(in_neighbors, out_neighbors, new_node - 1, new_node);
    return new_node;
}

void generate_dot(const node_to_reads_t& node_to_read,
                  const sequence_t& gene,
                  const exonic_indicator_t& exonic,
                  const out_neighbors_t& children,
                  const string output_path) {
    stringstream output;
    output << "digraph graphname{" << endl;
    output << "    rankdir=LR;" << endl;
    for (node_id_t node = 0; node < children.size(); node++) {
        string outline = "peripheries=";
        outline +=  exonic[node] ? "2" : "1";
        output << "    " << node <<" [label=\"" << node << ":" << gene[node] << ":" << node_to_read[node].size() << "\" "<< outline <<"]" << endl;
    }
    output << endl;
    for (node_id_t node = 0; node < children.size() - 1; node++) {
        output << "    " << node <<"->" << node + 1 << endl;
    }
    for (node_id_t node = 0; node < children.size() - 1; node++) {
        for (node_id_t child : children[node]) {
            output << "    " << node <<"->" << child << endl;
        }
    }
    output << "}" << endl;
    ofstream ofile;
    ofile.open(output_path);
    ofile << output.str();
    ofile.close();
}

// void print_read_exonic_alignment(const read_exonic_alignment_t& opt_alignment,
//                                  const align_score_t& opt_score,
//                                  const sequence_t& read,
//                                  const sequence_t& gene) {
//     size_t start, end, length;
//     start = opt_alignment.first.first;
//     end = opt_alignment.first.second;
//     length = end - start + 1;
//
//     cout << "Read: " << read.substr(start, length) <<" ("<< start <<","<< end <<")" << endl;
//     cout << "Gene: ";
//     for (auto & fragment : opt_alignment.second) {
//         start = fragment.first;
//         end = fragment.second;
//         length = end - start + 1;
//         cout << gene.substr(start, length) << " ("<< start <<","<< end <<")-";
//     }
//     cout << endl;
//     cout << "Alignment score: " << opt_score << endl;
//
// }

void print_matrix(const sequence_t& read,
                  const sequence_t& gene,
                  const align_matrix_t& D,
                  const backtrack_matrix_t& B){
    cout << "\\\t";
    sequence_t gene_temp = "$" + gene;
    for (size_t j = 0; j < gene_temp.size(); j++) {
        cout << j << "(" << gene_temp[j] << ")\t";
    }
    cout << endl;
    sequence_t read_temp = "$" + read;
    for (size_t i = 0; i < D.size(); i++) {
        cout << i << "(" << read_temp[i] << ")\t";
        for (size_t j = 0; j < D[i].size(); j++) {
            cout << D[i][j] << "(" << B[i][j].first << "," << B[i][j].second << ")\t";
        }
        cout << endl;
    }
}
