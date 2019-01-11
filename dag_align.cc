#include "dag_align.h"

#include <fstream>  // std::ofstream
#include <sstream>  // std::stringstream
#include <queue>
#include <algorithm> // std::reverse
#include <functional> // std::function
#include <stdlib.h> //abort()

#define EXONIC_ID 0x7FFFFFFF
#define INTRONIC_ID 0x7FFFFFFE
#define MATCH_S 1
#define GAP_S -1
#define MISMATCH_S -1

constexpr size_t MAX_MAPPINGS = 10;
constexpr align_score_t MIN_SCORE = 3;
constexpr matrix_coordinate_t INVALID_COORDINATE = {-1,-1};

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ofstream;
using std::queue;
using std::vector;
using std::function;
using std::reverse;

void process_gene_test() {
    sequence_list_t reads(1, "AAANNNNNNNNCCC");
    sequence_t gene = "AAACCC";
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

        read_gene_mappings_t mappings;
        vector<align_score_t> scores;
        for (size_t j = 0; j < MAX_MAPPINGS; j++) {
            align_path_t opt_alignment;
            align_score_t opt_score;
            extract_local_alignment(opt_alignment, opt_score, D, B);
            recalc_alignment_matrix(D, B, opt_alignment, children, local_aligner);
            cout << "Next best local alignment path score is " << opt_score << endl;
            if (opt_score > MIN_SCORE) {
                mappings.emplace_back(get_mapping_intervals(opt_alignment));
                scores.emplace_back(opt_score);
                cout << "Added mapping:" << endl;
                print_mapping_interval(opt_score, mappings[mappings.size()-1], reads[i], gene);
            } else {
                cout << "Alignemnt too poor. Breaking..." << endl;
                break;
            }
            print_matrix(reads[i], gene, D, B);
        }
        for (size_t j = 0; j < mappings.size(); j++) {
            print_mapping_interval(scores[j], mappings[j], reads[i], gene);
        }
        read_gene_mappings_t opt_cochain = get_optimal_cochain(scores, mappings);
        cout << "Optimal cochain:" << endl;
        print_cochain(opt_cochain);
    }
}

read_gene_mappings_t get_optimal_cochain(const vector<align_score_t>& scores,
                                         const read_gene_mappings_t& mappings) {
   constexpr size_t no_parent = -1;
   vector<size_t> D(scores.size(), 0);
   vector<size_t> B(scores.size(), no_parent);
   vector<bool> done(scores.size(), false);

   auto is_parent = [&mappings] (size_t child, size_t parent) {
       if (child == parent) {
           return false;
       }
       const read_interval_t& child_read_interval = mappings[child].first;
       const read_interval_t& parent_read_interval = mappings[parent].first;
       const gene_interval_t& child_first_gene_interval = mappings[child].second[0];
       const gene_interval_t& parent_last_gene_interval = mappings[parent].second[mappings[parent].second.size()-1];
       if (child_read_interval.first <= parent_read_interval.second) {
           return false;
       }
       if (child_first_gene_interval.first <= parent_last_gene_interval.second) {
           return false;
       }
       return true;
   };
   function<void (size_t)> recurse;
   recurse = [&D, &B, &scores, &mappings, &done, &is_parent, &recurse] (size_t i) -> void {
        if (done[i]) {
            return;
        }
        size_t max_parent_value = 0;
        size_t max_parent = no_parent;
        for (size_t j = 0; j < scores.size(); j++) {
            if (!is_parent(j, i)) {
                continue;
            }
            cout << j << "<--" << i << endl;
            if (D[j] > max_parent_value) {
                max_parent_value = D[j];
                max_parent = j;
            }
        }
        D[i] = scores[i] + max_parent_value;
        B[i] = max_parent;
        done[i] = true;
   };

   size_t opt_chain_value = 0;
   size_t opt_chain_tail = -1;
   for (size_t i = 0; i < scores.size(); i++) {
       recurse(i);
       if (D[i] > opt_chain_value) {
           opt_chain_value = D[i];
           opt_chain_tail = i;
       }
   }
   read_gene_mappings_t result;
   result.push_back(mappings[opt_chain_tail]);
   while (B[opt_chain_tail] != no_parent) {
       opt_chain_tail = B[opt_chain_tail];
       result.push_back(mappings[opt_chain_tail]);
   }

   return result;
}

read_gene_mapping_t get_mapping_intervals(const align_path_t& path) {
    read_gene_mapping_t result(
        read_interval_t(path[0].first, path[path.size()-1].first),
        gene_intervals_t(1,
            gene_interval_t(path[0].second, path[0].second)
        )
    );
    gene_intervals_t& gene_intervals = result.second;
    for (const matrix_coordinate_t& pos : path) {
        gene_interval_t& cur_gene_interval = gene_intervals[gene_intervals.size()-1];
        if (pos.second - cur_gene_interval.second > 1) {
            gene_intervals.push_back(gene_interval_t(pos.second, pos.second));
        } else {
            cur_gene_interval.second = pos.second;
        }
    }
    return result;
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
    auto local_aligner = [&D, &B, &read, &gene, &exonic, &parents] (const backtrack_t& i, const backtrack_t& j) {
        auto set_to_max_match = [&D, &read, &gene, &exonic] (align_score_t& opt_s, matrix_coordinate_t& opt_b, const backtrack_t& row, const backtrack_t& col) {
            align_score_t cur_s = D[row][col];
            cur_s += match(read[row], gene[col], exonic[col]);
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
        align_score_t opt_s = 0;
        matrix_coordinate_t opt_b = INVALID_COORDINATE;
        // Three direct parents parents
        set_to_max_gap(opt_s, opt_b, i-1, j-0); // Delete
        set_to_max_gap(opt_s, opt_b, i-0, j-1); // Insert
        set_to_max_match(opt_s, opt_b, i-1, j-1); // Match
        // Other DAG parents
        for (node_id_t parent : parents[j-1]) {
            size_t parent_j = parent + 1;
            // Three indirect parents
            set_to_max_gap(opt_s, opt_b, i-1, parent_j-0); // Delete
            set_to_max_gap(opt_s, opt_b, i-0, parent_j-1); // Insert
            set_to_max_match(opt_s, opt_b, i-1, parent_j-1); // Match
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

    while (clearing_queue.size() > 0) {
        matrix_coordinate_t pos = clearing_queue.front();
        queue_children(pos);
        if (B[pos.first][pos.second] != INVALID_COORDINATE) {
            local_aligner(pos.first, pos.second);
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

void print_cochain(const read_gene_mappings_t& chain) {
    for (read_gene_mapping_t mapping : chain) {
        read_interval_t& read_interval = mapping.first;
        cout << "R" << "("<<read_interval.first<<","<<read_interval.second<<") "<<endl;
        for (gene_interval_t gene_interval : mapping.second) {
            cout << "G" << "("<<gene_interval.first<<","<<gene_interval.second<<")"<<'-';
        }
        cout << endl;
    }
}

void print_mapping_interval(const align_score_t& score,
                            const read_gene_mapping_t& read_gene_mapping,
                            const sequence_t& read,
                            const sequence_t& gene) {
    size_t start, end, length;

    start = read_gene_mapping.first.first;
    end = read_gene_mapping.first.second;
    length = end - start + 1;
    cout << score << ": ";
    cout << read.substr(start-1, length) <<"("<<start<<","<<end<<")" << " --> ";
    for (auto & fragment : read_gene_mapping.second) {
        start = fragment.first;
        end = fragment.second;
        length = end - start + 1;
        cout << gene.substr(start-1, length) <<"("<<start<<","<<end<<")-";
    }
    cout << endl;
}

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
