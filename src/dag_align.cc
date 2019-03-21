#include "dag.h"
#include "commandline_flags.h"
#include "utils.h"
#include "fmt/format.h"

#include <iostream> // std::cerr, std::cout
#include <sstream>  // std::stringstream
#include <queue>
#include <unordered_map>
#include <algorithm> // std::reverse, std::sort
#include <functional> // std::function
#include <cstdlib> // abort()
#include <utility> // move()
#include <iterator> // next(), std::inserter

using namespace dag_types;

constexpr align_score_t EXONIC_S = 0;
constexpr align_score_t MATCH_S = 1;
constexpr align_score_t GAP_S = -6;
constexpr align_score_t MISMATCH_S = -6;
constexpr size_t MAX_MAPPINGS = 10;
constexpr align_score_t MIN_SCORE = 30;
constexpr matrix_coordinate_t INVALID_COORDINATE = {-1,-1};
constexpr index_t COCHAINING_PERMISSIBILITY = 10; //in inches

using fmt::format;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::queue;
using std::vector;
using std::unordered_map;
using std::function;
using std::sort;
using std::reverse;
using std::move;
using std::next;
using std::set;
using std::make_pair;

void dag_aligner::local_aligner(const index_t& i, const index_t& j) {
    // Lambda function to update the opt value
    auto set_to_max = [&, this, i, j](align_score_t& opt_s, matrix_coordinate_t& opt_b, matrix_coordinate_t& source, align_score_t cur_s) {
        cur_s += this->D[source.first][source.second];
        if (cur_s > opt_s) {
            opt_s = cur_s;
            opt_b = {source.first, source.second};
        }
    };
    align_score_t opt_s = 0;
    matrix_coordinate_t opt_b = INVALID_COORDINATE;
    align_score_t matching_score = 0;
    if (this->read[i] == this->gene[j]) {
        matching_score += MATCH_S + EXONIC_S*this->exonic_indicator[j];
    } else {
        matching_score += MISMATCH_S;
    }
    matrix_coordinate_t source;
    source = {i-1, j-0};
    set_to_max(opt_s, opt_b, source, GAP_S); // Consume read
    // Direct parent column
    source = {i-0, j-1};
    set_to_max(opt_s, opt_b, source, GAP_S); // Consume gene
    source = {i-1, j-1};
    set_to_max(opt_s, opt_b, source, matching_score); // Consume both
    // Other DAG parents
    for (index_t parent : parents[j]) {
        // Parent column from DAG
        source = {i-0, parent};
        set_to_max(opt_s, opt_b, source, GAP_S); // Consume gene
        source = {i-1, parent};
        set_to_max(opt_s, opt_b, source, matching_score); // Consume both
    }
    D[i][j] = opt_s;
    B[i][j] = opt_b;
}

void dag_aligner::extract_local_alignment() {
    align_path_t opt_alignment;
    // First, find optimal score in D
    matrix_coordinate_t opt_tail = {1,1};
    align_score_t opt_score = D[1][1];
    for (index_t i = 1; i < read.size(); i++) {
        for (index_t j = 1; j < gene.size(); j++) {
            if (D[i][j] > opt_score) {
                opt_score = D[i][j];
                opt_tail = {i,j};
            }
        }
    }
    opt_alignment.push_back(opt_tail);
    if (opt_score == 0) {
        align_scores.push_back(opt_score);
        align_paths.emplace_back(opt_alignment);
        return;
    }
    // Then, backtrack from the opt score back to the first positive score in the path
    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B[cur_pos.first][cur_pos.second];
    while (nxt_pos != INVALID_COORDINATE && D[nxt_pos.first][nxt_pos.second] > 0) {
        opt_alignment.push_back(nxt_pos);
        cur_pos = nxt_pos;
        nxt_pos = B[cur_pos.first][cur_pos.second];
    }
    // Make sure to have the path in the correct orientation (top to bottom, left to right)
    reverse(opt_alignment.begin(), opt_alignment.end());
    align_scores.push_back(opt_score);
    align_paths.emplace_back(opt_alignment);
}

void dag_aligner::recalc_alignment_matrix() {
    queue<matrix_coordinate_t> clearing_queue;
    // First, resets alignment path so it can't be used by new alignments. Queues the path nodes.
    for (matrix_coordinate_t pos : align_paths[align_paths.size()-1]) {
        D[pos.first][pos.second] = 0;
        B[pos.first][pos.second] = INVALID_COORDINATE;
        clearing_queue.push(pos);
    }
    // Process progressively each entry in the queue and add its children to the back of the queue to be processed in turn.
    while (clearing_queue.size() > 0) {
        matrix_coordinate_t pos = clearing_queue.front();
        // queue_children(pos);
        {
            // A sub-lambda function. Queues a possible child if it is an actual child
            auto branches_from_pos = [&, this, pos] (const matrix_coordinate_t& descendant) {
                if (descendant.first >= this->read.size()) {
                    return false;
                }
                if (descendant.second >= this->gene.size()) {
                    return false;
                }
                if (B[descendant.first][descendant.second] != pos) {
                    return false;
                }
                return true;
            };
            // First, check the immediate possible children (right, under, and right-under corner)
            //   Note that the the third child must always be the corner since it possibly depend on the other two children
            matrix_coordinate_t descendant;
            descendant = {pos.first + 0, pos.second + 1};
            if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
            descendant = {pos.first + 1, pos.second + 0};
            if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
            descendant = {pos.first + 1, pos.second + 1};
            if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
            // Then, check possible children that come through DAG edges
            for (const index_t& child : this->children[pos.second]) {
                // Note that the the third child must always be the corner since it possibly depend on the other two children
                descendant = {pos.first + 0, child};
                if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
                descendant = {pos.first + 1, child};
                if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
            }
        }
        if (B[pos.first][pos.second] != INVALID_COORDINATE) {
            local_aligner(pos.first, pos.second);
        }
        clearing_queue.pop();
    }
}

void dag_aligner::compress_align_paths() {
    local_mappings.resize(align_paths.size());
    for (size_t i =0; i < align_paths.size(); i++) {
        const align_path_t& path = align_paths[i];
        // We know the read interval and the start of the first gene interval
        mapping_t result(
            interval_t(path[0].first, path[path.size()-1].first),
            vector<interval_t>(1,
                interval_t(path[0].second, path[0].second)
            )
        );
        // Compresses the gene intervals and appends them as we find a jump in the alignment path
        vector<interval_t>& gene_intervals = result.second;
        for (const matrix_coordinate_t& pos : path) {
            interval_t& cur_gene_interval = gene_intervals[gene_intervals.size()-1];
            if (pos.second - cur_gene_interval.second > 1) {
                gene_intervals.push_back(interval_t(pos.second, pos.second));
            } else {
                cur_gene_interval.second = pos.second;
            }
        }
        local_mappings[i] = result;
    }
}

void dag_aligner::cochain_mappings() {
    if (align_scores.size() == 0) {
        return;
    }
    constexpr size_t no_parent = -1;
    vector<size_t> D(align_scores.size(), 0);
    vector<size_t> B(align_scores.size(), no_parent);
    vector<bool> done(align_scores.size(), false);
    opt_chain_indicator = vector<bool>(align_scores.size(), false);

    // A lamda function to check if two mappings can be in a child-parent relationship in a chain
    auto is_parent = [&mappings = this->local_mappings](size_t child, size_t parent) {
        // No mapping can parent itself
        if (child == parent) {
            return false;
        }
        const interval_t& child_read_interval = mappings[child].first;
        const interval_t& parent_read_interval = mappings[parent].first;
        const interval_t& child_first_gene_interval = mappings[child].second[0];
        const interval_t& parent_last_gene_interval = mappings[parent].second[mappings[parent].second.size()-1];
        // The start of the child read interval CANNOT come before the end of the parent read inverval
        if (child_read_interval.first + COCHAINING_PERMISSIBILITY <= parent_read_interval.second) {
            return false;
        }
        // The start of the child first gene interval CANNOT come before the end of the parent last gene inverval
        if (child_first_gene_interval.first + COCHAINING_PERMISSIBILITY <= parent_last_gene_interval.second) {
            return false;
        }
        return true;
    };
    // Recursive lambda function. Computes optimal co linear chain ending at a given mapping.
    //   Recursively computes any possible parents of the mapping before computing score for the given mapping.
    function<void (size_t)> compute_fragment_opt_score;
    compute_fragment_opt_score = [&compute_fragment_opt_score, &D, &B, &done, &is_parent, &scores = this->align_scores, &mappings = this->local_mappings] (size_t fragment_id) -> void {
        if (done[fragment_id]) {
            return;
        }
        size_t max_parent_value = 0;
        size_t max_parent_id = no_parent;
        for (size_t parent_id = 0; parent_id < scores.size(); parent_id++) {
            if (!is_parent(fragment_id, parent_id)) {
                continue;
            }
            compute_fragment_opt_score(parent_id);
            if (D[parent_id] > max_parent_value) {
                max_parent_value = D[parent_id];
                max_parent_id = parent_id;
            }
        }
        D[fragment_id] = scores[fragment_id] + max_parent_value;
        B[fragment_id] = max_parent_id;
        done[fragment_id] = true;
    };

    size_t opt_chain_value = 0;
    size_t opt_chain_tail = -1;
    // Find the optimal score ending with each mapping interval
    for (size_t i = 0; i < align_scores.size(); i++) {
        compute_fragment_opt_score(i);
        // Record the best of them
        if (D[i] > opt_chain_value) {
            opt_chain_value = D[i];
            opt_chain_tail = i;
        }
    }
    // Backtrack from the best tail
    opt_chain.push_back(opt_chain_tail);
    opt_chain_indicator[opt_chain_tail] = true;
    while (B[opt_chain_tail] != no_parent) {
        opt_chain_tail = B[opt_chain_tail];
        opt_chain.push_back(opt_chain_tail);
        opt_chain_indicator[opt_chain_tail] = true;
    }
    reverse(opt_chain.begin(), opt_chain.end());
}

//// Public functions
void dag_aligner::align_read(const string& read_name_in, const string& read_in) {
    clear_read_structures();
    read = read_in;
    if (read_name_to_id.find(read_name_in) == read_name_to_id.end()) {
        read_id = read_names.size();
        read_names.push_back(read_name_in);
        read_name_to_id[read_name_in] = read_id;
    } else {
        return;
    }
    index_t read_l = read.size();
    index_t gene_l = gene.size();
    D = align_matrix_t(read_l, align_row_t(gene_l, 0));
    B = backtrack_matrix_t(read_l, backtrack_row_t(gene_l, INVALID_COORDINATE));
    for (index_t i = 1; i < read_l; i++) {
        for (index_t j = 1; j < gene_l; j++) {
            local_aligner(i, j);
        }
    }
    for (size_t i = 0; i < MAX_MAPPINGS; i++) {
        extract_local_alignment();
        if (align_scores[align_scores.size() -1] < MIN_SCORE) {
            align_paths.pop_back();
            align_scores.pop_back();
            break;
        }
        recalc_alignment_matrix();
    }
    compress_align_paths();
    cochain_mappings();
    update_dag();
}

void dag_aligner::print_last_read_alignments() {
    stringstream ss;
    for (size_t i = 0; i < align_paths.size(); i++) {
        ss << format("{}\t", read_names[read_id]); //  Query sequence name
        ss << format("{}\t", read.size()-1); //  Query sequence length
        ss << format("{}\t", local_mappings[i].first.first-1); //  Query start (0-based)
        ss << format("{}\t", local_mappings[i].first.second-1); //  Query end (0-based)
        ss << format("{}\t", "+"); //  Relative strand: "+" or "-"
        ss << format("{}\t", gene_name); //  Target sequence name
        ss << format("{}\t", gene.size()-1); //  Target sequence length
        stringstream gene_starts;
        stringstream gene_ends;
        stringstream aln_blk_len;
        for (const interval_t& exon : local_mappings[i].second) {
            if (exon == local_mappings[i].second[0]) {
                gene_starts << format("{}", exon.first-1);
                gene_ends << format("{}", exon.second-1);
                aln_blk_len << format("{}", exon.second-exon.first);
            } else {
                gene_starts << format(",{}", exon.first-1);
                gene_ends << format(",{}", exon.second-1);
                aln_blk_len << format(",{}", exon.second-exon.first);
            }
        }
        ss << format("{}\t", gene_starts.str());  //  Target start on original strand (0-based)
        ss << format("{}\t", gene_ends.str());  //  Target end on original strand (0-based)
        ss << format("{}\t", 0); //  Number of residue matches
        ss << format("{}\t", aln_blk_len.str()); //  Alignment block length
        ss << format("{}\t", 255); //  Mapping quality (0-255; 255 for missing)
        ss << format("s1:i:{:d}\t", align_scores[i]); //  Local alignment score
        ss << format("oc:c:{:d}\n", opt_chain_indicator[i]);
    }
    cout << ss.str();
}
