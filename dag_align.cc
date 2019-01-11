#include "dag_align.h"

#include <fstream>  // std::ofstream
#include <sstream>  // std::stringstream
#include <queue>
#include <algorithm> // std::reverse
#include <functional> // std::function
#include <cstdlib> // abort()

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

// Testing function
void process_gene_test() {
    sequence_list_t reads(1, "AAANNNNNNNNCCC");
    sequence_t gene = "AAACCC";
    exonic_indicator_t exonic(gene.size(), true);
    process_gene(reads, gene, exonic);
}

// An entry point function
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
        read_gene_mappings_t opt_chain = align_read_to_dag(reads[i], gene, exonic, parents, children);
        update_dag(parents, children, node_to_read, i, opt_chain);
    }
}

void update_dag(in_neighbors_t& parents,
                out_neighbors_t& children,
                node_to_reads_t& node_to_read,
                const size_t& read_id,
                const read_gene_mappings_t& opt_chain) {
    ;
}

// Align a read to the DAG
// The alignment is local alignment
// The top MAX_MAPPINGS local alignments that are at least MIN_SCORE worth are considered for chaining
// The optimal chain is be returned
read_gene_mappings_t align_read_to_dag(const sequence_t& read,
                                       const sequence_t& gene,
                                       const exonic_indicator_t& exonic,
                                       const in_neighbors_t& parents,
                                       const out_neighbors_t& children) {
    align_matrix_t D;
    backtrack_matrix_t B;
    // Get an aligner function that has matrices, gene, read, and DAG captured
    local_aligner_func_t local_aligner = get_local_aligner(D, B, read, gene, exonic, parents);
    // Perfrom first round of local alignment for the whole read on the whole DAG
    local_alignment(D, B, read, gene, local_aligner);
    print_matrix(read, gene, D, B);

    // Now, we want to extract the best local alignments from D&B alignment matrices
    read_gene_mappings_t mappings;
    vector<align_score_t> scores;
    for (size_t i = 0; i < MAX_MAPPINGS; i++) {
        align_path_t opt_alignment;
        align_score_t opt_score;
        // This will return the sorted path of the best alignment in D
        extract_local_alignment(opt_alignment, opt_score, D, B);
        // After getting the best alignment path, we need to reset that alignment path to zeroes,
        //   and realign any entries that depended on this path
        recalc_alignment_matrix(D, B, opt_alignment, children, local_aligner);
        cout << "Next best local alignment path score is " << opt_score << endl;
        if (opt_score > MIN_SCORE) {
            mappings.emplace_back(get_mapping_intervals(opt_alignment));
            scores.emplace_back(opt_score);
            cout << "Added mapping:" << endl;
            print_mapping_interval(opt_score, mappings[mappings.size()-1], read, gene);
        } else {
            cout << "Alignemnt too poor. Breaking..." << endl;
            break;
        }
        print_matrix(read, gene, D, B);
    }
    for (size_t i = 0; i < mappings.size(); i++) {
        print_mapping_interval(scores[i], mappings[i], read, gene);
    }
    // Using the extracted alignments, find the optimal co-linear chain of them
    read_gene_mappings_t opt_cochain = get_optimal_cochain(scores, mappings);
    cout << "Optimal cochain:" << endl;
    print_cochain(opt_cochain);
    return opt_cochain;
}

// Returns the optimal co-linear chain of local alignments
read_gene_mappings_t get_optimal_cochain(const vector<align_score_t>& scores,
                                         const read_gene_mappings_t& mappings) {
    constexpr size_t no_parent = -1;
    vector<size_t> D(scores.size(), 0);
    vector<size_t> B(scores.size(), no_parent);
    vector<bool> done(scores.size(), false);

    // A lamda function to check if two mappings can be in a child-parent relationship in a chain
    auto is_parent = [&mappings] (size_t child, size_t parent) {
        // No mapping can parent itself
        if (child == parent) {
            return false;
        }
        const read_interval_t& child_read_interval = mappings[child].first;
        const read_interval_t& parent_read_interval = mappings[parent].first;
        const gene_interval_t& child_first_gene_interval = mappings[child].second[0];
        const gene_interval_t& parent_last_gene_interval = mappings[parent].second[mappings[parent].second.size()-1];
        // The start of the child read interval CANNOT come before the end of the parent read inverval
        if (child_read_interval.first <= parent_read_interval.second) {
            return false;
        }
        // The start of the child first gene interval CANNOT come before the end of the parent last gene inverval
        if (child_first_gene_interval.first <= parent_last_gene_interval.second) {
            return false;
        }
        return true;
    };
    // Recursive lambda function. Computes optimal co linear chain ending at a given mapping.
    //   Recursively computes any possible parents of the mapping before computing score for the given mapping.
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
    // Find the optimal score ending with each mapping interval
    for (size_t i = 0; i < scores.size(); i++) {
        recurse(i);
        // Record the best of them
        if (D[i] > opt_chain_value) {
            opt_chain_value = D[i];
            opt_chain_tail = i;
        }
    }
    // Backtrack from the best tail
    read_gene_mappings_t result;
    result.push_back(mappings[opt_chain_tail]);
    while (B[opt_chain_tail] != no_parent) {
        opt_chain_tail = B[opt_chain_tail];
        result.push_back(mappings[opt_chain_tail]);
    }
    reverse(result.begin(), result.end());
    return result;
}

// Compresses an alignment path to a read interval and one or more gene intervals
read_gene_mapping_t get_mapping_intervals(const align_path_t& path) {
    // We know the read interval and the start of the first gene interval
    read_gene_mapping_t result(
        read_interval_t(path[0].first, path[path.size()-1].first),
        gene_intervals_t(1,
            gene_interval_t(path[0].second, path[0].second)
        )
    );
    // Compresses the gene intervals and appends them as we find a jump in the alignment path
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

// Returns a local aligner function.
//   The lambda function takes only row and column coordinates of a cell to be computed
//   All other needed information (matrices, read, gene, and DAG) are captured at instatiation time of the aligner function
local_aligner_func_t get_local_aligner(align_matrix_t& D,
                                       backtrack_matrix_t& B,
                                       const sequence_t& read,
                                       const sequence_t& gene,
                                       const exonic_indicator_t& exonic,
                                       const out_neighbors_t& parents) {
    // Lambda function. The local aligner  that will be retured
    auto local_aligner = [&D, &B, &read, &gene, &exonic, &parents] (const backtrack_t& i, const backtrack_t& j) {
        // Sub-lambda function to update the opt value when matching two characters
        auto set_to_max_match = [&D, &read, &gene, &exonic] (align_score_t& opt_s, matrix_coordinate_t& opt_b, const backtrack_t& row, const backtrack_t& col) {
            align_score_t cur_s = D[row][col];
            if (read[row] == gene[col]) {
                cur_s += MATCH_S + MATCH_S*exonic[col];
            } else {
                cur_s += MISMATCH_S;
            }
            if (cur_s > opt_s) {
                opt_s = cur_s;
                opt_b = {row, col};
            }
        };
        // Sub-lambda function to update the opt value when adding a gap
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

// Creates and computes alignment matrix using the local aligner function
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

// Given an alignemnt path, resets the path to zeros and it B entries to INVALID_COORDINATE;
// Then, for each D entry that its path traces back to this path, recomputes the alignment using the aligner function
void recalc_alignment_matrix(align_matrix_t& D,
                             backtrack_matrix_t& B,
                             const align_path_t& opt_alignment,
                             const in_neighbors_t& children,
                             const local_aligner_func_t& local_aligner) {
    queue<matrix_coordinate_t> clearing_queue;
    // A lambda function. Given a matrix coordinate, checks each of its possible children
    // If any of possible children are actual children, then adds to the queue to be processed
    auto queue_children = [&D, &B, &children, &clearing_queue] (const matrix_coordinate_t& pos) {
        // A sub-lambda function. Queues a possible child if it is an actual child
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
        // First, check the immediate possible children (right, under, and right-under corner)
        //   Note that the the third child must always be the corner since it possibly depend on the other two children
        matrix_coordinate_t descendant;
        descendant = {pos.first + 0, pos.second + 1};
        queue_if_child(descendant);
        descendant = {pos.first + 1, pos.second + 0};
        queue_if_child(descendant);
        descendant = {pos.first + 1, pos.second + 1};
        queue_if_child(descendant);
        // Then, check possible children that come through DAG edges
        for (const node_id_t& child : children[pos.second - 1]) {
            // Note that the the third child must always be the corner since it possibly depend on the other two children
            descendant = {pos.first + 0, child + 1};
            queue_if_child(descendant);
            descendant = {pos.first + 1, child + 0};
            queue_if_child(descendant);
            descendant = {pos.first + 1, child + 1};
            queue_if_child(descendant);
        }
    };

    // First, resets alignment path so it can't be used by new alignments. Queues the path nodes.
    for (matrix_coordinate_t pos : opt_alignment) {
        D[pos.first][pos.second] = 0;
        B[pos.first][pos.second] = INVALID_COORDINATE;
        clearing_queue.push(pos);
    }

    // Process progressively each entry in the queue and add its children to the back of the queue to be processed in turn.
    while (clearing_queue.size() > 0) {
        matrix_coordinate_t pos = clearing_queue.front();
        queue_children(pos);
        if (B[pos.first][pos.second] != INVALID_COORDINATE) {
            local_aligner(pos.first, pos.second);
        }
        clearing_queue.pop();
    }
}

// Extracts the best local alignment path.
void extract_local_alignment(align_path_t& opt_alignment,
                             align_score_t& opt_score,
                             const align_matrix_t& D,
                             const backtrack_matrix_t& B) {
    opt_alignment.clear();
    // First, find optimal score in D
    matrix_coordinate_t opt_tail = {0,0};
    opt_score = D[0][0];
    for (size_t i = 0; i < D.size(); i++) {
        for (size_t j = 0; j < D[0].size(); j++) {
            if (D[i][j] > opt_score) {
                opt_score = D[i][j];
                opt_tail = {i,j};
            }
        }
    }
    opt_alignment.push_back(opt_tail);
    // Then, backtrack from the opt score back to the first positive score in the path
    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B[cur_pos.first][cur_pos.second];
    while (D[nxt_pos.first][nxt_pos.second] > 0) {
        opt_alignment.push_back(nxt_pos);
        cur_pos = nxt_pos;
        nxt_pos = B[cur_pos.first][cur_pos.second];
    }
    // Make sure to have the path in the correct orientation (top to bottom, left to right)
    reverse(opt_alignment.begin(), opt_alignment.end());
}

void add_edge(in_neighbors_t &in_neighbors,
              out_neighbors_t &out_neighbors,
              const node_id_t&  source,
              const node_id_t&  target) {
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

// Adds a node to the ends of the graph vectors
node_id_t append_node(in_neighbors_t &in_neighbors,
                      out_neighbors_t &out_neighbors,
                      node_to_reads_t &node_to_read) {
    if (out_neighbors.size() != in_neighbors.size()) {
        cout << "ERR:Graph corrupted. out_neighbors.size()!=in_neighbors.size()" << endl;
        abort();
    }
    if (node_to_read.size() != in_neighbors.size()) {
        cout << "ERR:Graph corrupted. node_to_read.size()!=in_neighbors.size()" << endl;
        abort();
    }
    node_id_t new_node = in_neighbors.size();
    in_neighbors.push_back(node_list_t());
    out_neighbors.push_back(node_list_t());
    node_to_read.push_back(read_list_t());
    return new_node;
}

// A function to generate a .DOT file of the DAG
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

// Printing here is formatted to (kinda) work with "column -t" Linux program
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
