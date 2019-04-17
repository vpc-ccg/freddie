#include "dag.h"
#include "utils.h"       // set_to_max()
#include "fmt/format.h"  // format()

#include <queue>
#include <iostream>   // std::cerr, std::cout
#include <sstream>    // std::stringstream
#include <algorithm>  // std::reverse, std::sort, std::max, std::min
#include <functional> // std::function
#include <cstdlib>    // abort()
#include <utility>    // move(), std::pair
#include <iterator>   // next(), std::inserter
#include <limits>     // numeric_limits
#include <math.h>     // ceil

using namespace dag_types;
using utils::set_to_max;
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
using std::function;
using std::sort;
using std::reverse;
using std::move;
using std::next;
using std::set;
using std::make_pair;
using std::pair;
using std::numeric_limits;
using std::ceil;

constexpr align_score_t EXONIC_S = 0;
constexpr align_score_t MATCH_S = 1;
constexpr align_score_t GAP_S = -6;
constexpr align_score_t MISMATCH_S = -6;
constexpr size_t MAX_MAPPINGS = 10;
constexpr align_score_t MIN_SCORE = 30;
constexpr matrix_coordinate_t INVALID_COORDINATE = {-1,-1};
constexpr index_t COCHAINING_PERMISSIBILITY = 10;
constexpr double MAX_UNALN_GENE_RATIO = 1.5;

void dag_aligner::local_aligner(align_matrix_t& D, backtrack_matrix_t& B, const string& read, const index_t& i, const index_t& j) {
    align_score_t opt_s = 0;
    matrix_coordinate_t opt_b = INVALID_COORDINATE;
    align_score_t matching_score = 0;
    if (read[i] == gene[j]) {
        matching_score = MATCH_S + EXONIC_S*exonic_indicator[j];
    } else {
        matching_score = MISMATCH_S;
    }
    matrix_coordinate_t source;
    source = {i-1, j-0};
    set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D[source.first][source.second] + GAP_S); // Consume read
    // Direct parent column
    source = {i-0, j-1};
    set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D[source.first][source.second] + GAP_S); // Consume gene
    source = {i-1, j-1};
    set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D[source.first][source.second] + matching_score); // Consume both
    // Other DAG parents
    for (index_t parent : nodes[j].parents) {
        // Parent column from DAG
        source = {i-0, parent};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D[source.first][source.second] + GAP_S); // Consume gene
        source = {i-1, parent};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D[source.first][source.second] + matching_score); // Consume both
    }
    D[i][j] = opt_s;
    B[i][j] = opt_b;
}

void dag_aligner::extract_local_alignment(local_alignment_s& loc_aln, const align_matrix_t& D, const backtrack_matrix_t& B, const string& read) {
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
    loc_aln.path.push_back(opt_tail);
    loc_aln.score = opt_score;
    if (opt_score == 0) {
        return;
    }
    // Then, backtrack from the opt score back to the first positive score in the path
    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B[cur_pos.first][cur_pos.second];
    while (nxt_pos != INVALID_COORDINATE && D[nxt_pos.first][nxt_pos.second] > 0) {
        loc_aln.path.push_back(nxt_pos);
        // Consume both read and gene
        if (cur_pos.first != nxt_pos.first && cur_pos.second != nxt_pos.second) {
            if (read[cur_pos.first] == gene[cur_pos.second]) {
                loc_aln.cigar.push_back(C_MAT);
            } else {
                loc_aln.cigar.push_back(C_MIS);
            }
        }
        // Consume read only
        if (cur_pos.first != nxt_pos.first && cur_pos.second == nxt_pos.second) {
                loc_aln.cigar.push_back(C_INS);
        }
        // Consume gene only
        if (cur_pos.first == nxt_pos.first && cur_pos.second != nxt_pos.second) {
                loc_aln.cigar.push_back(C_DEL);
        }
        cur_pos = nxt_pos;
        nxt_pos = B[cur_pos.first][cur_pos.second];
    }
    // Make sure to have the path in the correct orientation (top to bottom, left to right)
    reverse(loc_aln.path.begin(), loc_aln.path.end());
    reverse(loc_aln.cigar.begin(), loc_aln.cigar.end());
}

void dag_aligner::extract_affix_alignment(local_alignment_s& loc_aln, const align_matrix_dynamic_t& D, const backtrack_matrix_dynamic_t& B, const matrix_coordinate_t& start, const matrix_coordinate_t& end, const string& read) {
    // First, find optimal score in D
    matrix_coordinate_t opt_tail = INVALID_COORDINATE;
    align_score_t opt_score = numeric_limits<align_score_t>::min();
    for (index_t i = start.first; i <= end.first; i++) {
        for (index_t j = start.second; j <= end.second; j++) {
            matrix_coordinate_t coor(i,j);
            align_score_t score = D.at(coor);
            set_to_max<matrix_coordinate_t, align_score_t>(opt_tail, opt_score, coor, score);
        }
    }
    loc_aln.path.push_back(opt_tail);
    loc_aln.score = opt_score;
    if (opt_score == 0) {
        return;
    }
    // Then, backtrack from the opt score back to the first positive score in the path
    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B.at(cur_pos);
    while (nxt_pos != INVALID_COORDINATE && D.at(nxt_pos) > 0) {
        loc_aln.path.push_back(nxt_pos);
        // Consume both read and gene
        if (cur_pos.first != nxt_pos.first && cur_pos.second != nxt_pos.second) {
            if (read[cur_pos.first] == gene[cur_pos.second]) {
                loc_aln.cigar.push_back(C_MAT);
            } else {
                loc_aln.cigar.push_back(C_MIS);
            }
        }
        // Consume read only
        if (cur_pos.first != nxt_pos.first && cur_pos.second == nxt_pos.second) {
                loc_aln.cigar.push_back(C_INS);
        }
        // Consume gene only
        if (cur_pos.first == nxt_pos.first && cur_pos.second != nxt_pos.second) {
                loc_aln.cigar.push_back(C_DEL);
        }
        cur_pos = nxt_pos;
        nxt_pos = B.at(cur_pos);
    }
    // Make sure to have the path in the correct orientation (top to bottom, left to right)
    reverse(loc_aln.path.begin(), loc_aln.path.end());
    reverse(loc_aln.cigar.begin(), loc_aln.cigar.end());
}

void dag_aligner::recalc_alignment_matrix(align_matrix_t& D, backtrack_matrix_t& B, const local_alignment_s& loc_aln, const string& read) {
    queue<matrix_coordinate_t> clearing_queue;
    // First, resets alignment path so it can't be used by new alignments. Queues the path nodes.
    for (matrix_coordinate_t pos : loc_aln.path) {
        D[pos.first][pos.second] = 0;
        B[pos.first][pos.second] = INVALID_COORDINATE;
        clearing_queue.push(pos);
    }
    // Process progressively each entry in the queue and add its children to the back of the queue to be processed in turn.
    while (clearing_queue.size() > 0) {
        matrix_coordinate_t pos = clearing_queue.front();
        // A sub-lambda function. Queues a possible child if it is an actual child
        auto branches_from_pos = [&, read, gene_l = gene.size(), pos] (const matrix_coordinate_t& descendant) {
            if (descendant.first >= read.size()) {
                return false;
            }
            if (descendant.second >= gene_l) {
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
        for (const index_t& child : nodes[pos.second].children) {
            // Note that the the third child must always be the corner since it possibly depend on the other two children
            descendant = {pos.first + 0, child};
            if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
            descendant = {pos.first + 1, child};
            if (branches_from_pos(descendant)) {clearing_queue.push(descendant);}
        }

        if (B[pos.first][pos.second] != INVALID_COORDINATE) {
            local_aligner(D, B, read, pos.first, pos.second);
        }
        clearing_queue.pop();
    }
}

void dag_aligner::compress_align_path(local_alignment_s& loc_aln) {
    // We know the read interval and the start of the first gene interval
    interval_t read_interval (loc_aln.path.front().first, loc_aln.path.back().first);
    vector<interval_t> gene_intervals (1, interval_t(loc_aln.path.front().second, loc_aln.path.front().second));
    // Compresses the gene intervals and appends them as we find a jump in the alignment path
    for (const matrix_coordinate_t& pos : loc_aln.path) {
        interval_t& cur_gene_interval = gene_intervals.back();
        if (pos.second - cur_gene_interval.second > 1) {
            gene_intervals.push_back(interval_t(pos.second, pos.second));
        } else {
            cur_gene_interval.second = pos.second;
        }
    }
    string cigar_str = "";
    size_t last_cnt = 0;
    CIGAR_OP last_op = loc_aln.cigar.front();
    for (const CIGAR_OP& cigar_op : loc_aln.cigar) {
        if (last_op != cigar_op) {
            cigar_str += format("{:d}{:c}", last_cnt, (char)last_op);
            last_cnt = 0;
            last_op = cigar_op;
        }
        last_cnt++;
    }
    cigar_str += format("{:d}{:c}", last_cnt, (char)last_op);
    loc_aln.cigar_str = move(cigar_str);
    loc_aln.read_interval = move(read_interval);
    loc_aln.gene_intervals = move(gene_intervals);
}

void dag_aligner::cochain_mappings(vector<size_t> opt_chain, vector<local_alignment_s>& loc_alns) {
    if (loc_alns.size() == 0) {
        return;
    }
    constexpr size_t no_parent = -1;
    vector<size_t> D(loc_alns.size(), 0);
    vector<size_t> B(loc_alns.size(), no_parent);
    vector<bool> done(loc_alns.size(), false);

    // A lamda function to check if two mappings can be in a child-parent relationship in a chain
    auto is_parent = [&loc_alns](size_t child, size_t parent) {
        // No mapping can parent itself
        if (child == parent) {
            return false;
        }
        const interval_t& child_read_interval = loc_alns[child].read_interval;
        const interval_t& parent_read_interval = loc_alns[parent].read_interval;
        const interval_t& child_first_gene_interval = loc_alns[child].gene_intervals.front();
        const interval_t& parent_last_gene_interval = loc_alns[parent].gene_intervals.back();
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
    compute_fragment_opt_score = [&compute_fragment_opt_score, &D, &B, &done, &is_parent, &loc_alns] (size_t fragment_id) -> void {
        if (done[fragment_id]) {
            return;
        }
        size_t max_parent_value = 0;
        size_t max_parent_id = no_parent;
        for (size_t parent_id = 0; parent_id < loc_alns.size(); parent_id++) {
            if (!is_parent(fragment_id, parent_id)) {
                continue;
            }
            compute_fragment_opt_score(parent_id);
            if (D[parent_id] > max_parent_value) {
                max_parent_value = D[parent_id];
                max_parent_id = parent_id;
            }
        }
        D[fragment_id] = loc_alns[fragment_id].score + max_parent_value;
        B[fragment_id] = max_parent_id;
        done[fragment_id] = true;
    };

    size_t opt_chain_value = 0;
    size_t opt_chain_tail = -1;
    // Find the optimal score ending with each mapping interval
    for (size_t i = 0; i < loc_alns.size(); i++) {
        compute_fragment_opt_score(i);
        // Record the best of them
        if (D[i] > opt_chain_value) {
            opt_chain_value = D[i];
            opt_chain_tail = i;
        }
    }
    // Backtrack from the best tail
    opt_chain.push_back(opt_chain_tail);
    loc_alns[opt_chain_tail].in_opt_chain = true;
    while (B[opt_chain_tail] != no_parent) {
        opt_chain_tail = B[opt_chain_tail];
        opt_chain.push_back(opt_chain_tail);
        loc_alns[opt_chain_tail].in_opt_chain = true;
    }
    reverse(opt_chain.begin(), opt_chain.end());
}

void dag_aligner::affix_aligner(align_matrix_dynamic_t& D_affix, backtrack_matrix_dynamic_t&B_affix, align_matrix_dynamic_t& M_affix, const bool& is_prefix, const index_t& i, const index_t& j, const interval_t& gene_interval, const string& read) {
    matrix_coordinate_t coor(i,j);
    const index_t& g_start = gene_interval.first;
    const index_t& g_end   = gene_interval.second;
    align_score_t opt_s = numeric_limits<align_score_t>::min();
    align_score_t max_s = numeric_limits<align_score_t>::min();

    matrix_coordinate_t opt_b = INVALID_COORDINATE;
    align_score_t matching_score = 0;
    if (read[i] == gene[j]) {
        matching_score = MATCH_S + EXONIC_S*exonic_indicator[j];
    } else {
        matching_score = MISMATCH_S;
    }

    matrix_coordinate_t source;
    if (is_prefix) {
        source = {i-1, j-0};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_affix[source] + GAP_S); // Consume read
        max_s = std::max(max_s, M_affix[source]);
        // Direct parent column
        source = {i-0, j-1};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_affix[source] + GAP_S); // Consume gene
        max_s = std::max(max_s, M_affix[source]);
        source = {i-1, j-1};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_affix[source] + matching_score); // Consume both
        max_s = std::max(max_s, M_affix[source]);
        // Other DAG parents
        for (const index_t& parent : nodes[j].parents) {
            if (parent < g_start) {
                continue;
            }
            // Parent column from DAG
            source = {i-0, parent};
            set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_affix[source] + GAP_S); // Consume gene
            max_s = std::max(max_s, M_affix[source]);
            source = {i-1, parent};
            set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_affix[source] + matching_score); // Consume both
            max_s = std::max(max_s, M_affix[source]);
        }
    } else {
        source = {i+1, j+0};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, GAP_S); // Consume read
        max_s = std::max(max_s, M_affix[source]);
        // Direct parent column
        source = {i+0, j+1};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, GAP_S); // Consume gene
        max_s = std::max(max_s, M_affix[source]);
        source = {i+1, j+1};
        set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, matching_score); // Consume both
        max_s = std::max(max_s, M_affix[source]);
        // Other DAG children
        for (index_t child : nodes[j].children) {
            if (child > g_end) {
                continue;
            }
            // Child column from DAG
            source = {i+0, child};
            set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, GAP_S); // Consume gene
            max_s = std::max(max_s, M_affix[source]);
            source = {i+1, child};
            set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, matching_score); // Consume both
            max_s = std::max(max_s, M_affix[source]);
        }
    }
    max_s = std::max(max_s, opt_s);
    D_affix[coor] = opt_s;
    B_affix[coor] = opt_b;
    M_affix[coor] = max_s;
}

void dag_aligner::extend_opt_chain(vector<local_alignment_s>& loc_alns, vector<size_t>& opt_chain, const string& read) {
    size_t prev_mapping_id = opt_chain.front();
    for (const size_t& mapping_id : opt_chain) {
        if (mapping_id == opt_chain.front()) {
            continue;
        }
        const index_t& r_start = loc_alns[prev_mapping_id].read_interval.second;
        const index_t& r_end = loc_alns[mapping_id].read_interval.first;
        const index_t& g_start = loc_alns[prev_mapping_id].gene_intervals.front().second;
        const index_t& g_end = loc_alns[mapping_id].gene_intervals.back().first;
        const index_t affix_r_len = r_end - r_start;
        const index_t affix_g_len = std::max(g_end - g_start, (index_t) ceil(affix_r_len*MAX_UNALN_GENE_RATIO));
        align_matrix_dynamic_t D_prefix, M_prefix;
        backtrack_matrix_dynamic_t B_prefix;
        align_matrix_dynamic_t D_suffix, M_suffix;
        backtrack_matrix_dynamic_t B_suffix;
        matrix_coordinate_t coor, source, opt_b, max_sum_max_b;
        align_score_t opt_s, max_s, max_sum_max_s;
        // preprocess boundries of B_prefix and D_prefix
        coor = matrix_coordinate_t(r_start,g_start);
        source = INVALID_COORDINATE;
        D_prefix[coor] = 0;
        B_prefix[coor] = source;
        M_prefix[coor] = 0;
        for (index_t i = r_start + 1; i <= r_start + affix_r_len; i++) {
            coor   = matrix_coordinate_t(i, g_start);
            source = matrix_coordinate_t(i-1, g_start);
            opt_s = D_prefix[source] + GAP_S;
            opt_b = source;
            max_s = std::max(opt_s, M_prefix[source]);
            D_prefix[coor] = opt_s;
            B_prefix[coor] = opt_b;
            M_prefix[coor] = max_s;
        }
        for (index_t j = g_start + 1; j <= g_start + affix_g_len; j++) {
            coor   = matrix_coordinate_t(r_start, j);
            source = matrix_coordinate_t(r_start, j-1);
            opt_s = D_prefix[source] + GAP_S;
            opt_b = source;
            max_s = std::max(opt_s, M_prefix[source]);
            for (const index_t& parent : nodes[j].parents) {
                if (parent < g_start) {
                    continue;
                }
                source = matrix_coordinate_t(r_start, parent);
                set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_prefix[source] + GAP_S);
                max_s = std::max(opt_s, max_s);
            }
            D_prefix[coor] = opt_s;
            B_prefix[coor] = opt_b;
            M_prefix[coor] = max_s;
        }
        for (index_t i = r_start + 1; i <= r_start + affix_r_len; i++) {
            for (index_t j = g_start + 1; j <= g_start + affix_g_len; j++) {
                affix_aligner(D_prefix, B_prefix, M_suffix, true, i, j, interval_t(g_start, g_end), read);
            }
        }
        // preprocess boundries of B_suffix and D_suffix
        coor = matrix_coordinate_t(r_end,g_end);
        source = INVALID_COORDINATE;
        D_suffix[coor] = 0;
        B_suffix[coor] = source;
        for (index_t i = r_end - 1; i >= r_end - affix_r_len; i--) {
            coor   = matrix_coordinate_t(i, g_end);
            source = matrix_coordinate_t(i+1, g_end);
            opt_s = D_suffix[source] + GAP_S;
            opt_b = source;
            max_s = std::max(opt_s, M_suffix[source]);
            D_suffix[coor] = opt_s;
            B_suffix[coor] = opt_b;
            M_suffix[coor] = max_s;
        }
        for (index_t j = g_end - 1; j >= g_end - affix_g_len; j--) {
            coor   = matrix_coordinate_t(r_end, j);
            source = matrix_coordinate_t(r_end, j+1);
            opt_s = D_suffix[source] + GAP_S;
            opt_b = source;
            max_s = std::max(opt_s, M_suffix[source]);
            for (const index_t& child : nodes[j].children) {
                if (child > g_end) {
                    continue;
                }
                source = matrix_coordinate_t(r_end, child);
                set_to_max<matrix_coordinate_t, align_score_t>(opt_b, opt_s, source, D_suffix[source] + GAP_S);
                max_s = std::max(opt_s, max_s);
            }
            D_suffix[coor] = opt_s;
            B_suffix[coor] = opt_b;
            M_suffix[coor] = max_s;
        }
        for (index_t i = r_end - 1; i >= r_end - affix_r_len; i--) {
            for (index_t j = g_end - 1; j >= g_end - affix_g_len; j--) {
                affix_aligner(D_suffix, B_suffix, M_suffix, false, i, j, interval_t(g_start, g_end), read);
            }
        }
        // Sum of M_prefix and M_suffix
        max_sum_max_b = INVALID_COORDINATE;
        max_sum_max_s = numeric_limits<align_score_t>::min();
        for (const auto& kv : M_prefix) {
            const matrix_coordinate_t& M_prefix_coor = kv.first;
            const matrix_coordinate_t  M_suffix_coor = matrix_coordinate_t(
                coor.first,
                std::max(coor.second, g_end - affix_g_len)
            );
            set_to_max<matrix_coordinate_t, align_score_t>(max_sum_max_b, max_sum_max_s, M_prefix_coor, M_prefix[M_prefix_coor]+M_suffix[M_suffix_coor]);
        }
        for (const auto& kv : M_suffix) {
            const matrix_coordinate_t  M_prefix_coor = matrix_coordinate_t(
                coor.first,
                std::min(coor.second, g_start + affix_g_len)
            );
            const matrix_coordinate_t& M_suffix_coor = kv.first;
            set_to_max<matrix_coordinate_t, align_score_t>(max_sum_max_b, max_sum_max_s, M_prefix_coor, M_prefix[M_prefix_coor]+M_suffix[M_suffix_coor]);
        }
        local_alignment_s prefix_loc_aln;
        local_alignment_s suffix_loc_aln;
        extract_affix_alignment(prefix_loc_aln, D_prefix, B_prefix, matrix_coordinate_t(r_start, g_start), max_sum_max_b, read);
        extract_affix_alignment(suffix_loc_aln, D_suffix, B_suffix, max_sum_max_b, matrix_coordinate_t(r_end, g_end), read);
        compress_align_path(prefix_loc_aln);
        compress_align_path(suffix_loc_aln);
        loc_alns.emplace_back(prefix_loc_aln);
        loc_alns.emplace_back(suffix_loc_aln);
        prev_mapping_id = mapping_id;
    }
}

void dag_aligner::update_dag(const string& read_name, const vector<local_alignment_s>& loc_alns, const vector<size_t>& opt_chain) {
    seq_id_t read_id = read_name_to_id[read_name];
    aln_read_s& aln_read = aln_reads[read_id];
    if (opt_chain.size() == 0) {
        return;
    }
    for (const size_t& mapping_id : opt_chain) {
        mapping_s mapping;
        mapping.read_interval = loc_alns[mapping_id].read_interval;
        mapping.gene_intervals = loc_alns[mapping_id].gene_intervals;
        mapping.score = loc_alns[mapping_id].score;
        mapping.cigar_str = loc_alns[mapping_id].cigar_str;
        aln_read.mappings.emplace_back(mapping);
    }
    bool first_exon = true;
    interval_t prev_exon;
    for (const mapping_s& mapping : aln_read.mappings) {
        for (const interval_t& exon : mapping.gene_intervals) {
            for (index_t node = exon.first; node <= exon.second; node++) {
                nodes[node].read_ids.push_back(read_id);
            }
            if (first_exon) {
                first_exon = false;
                prev_exon = exon;
                continue;
            }
            edge_t e(prev_exon.second, exon.first);
            add_edge(e.first, e.second);
            edge_to_reads[e].push_back(read_id);
            prev_exon = exon;
        }
    }
}

//// Public functions
void dag_aligner::align_read(const string& read_name, const string& read) {
    if (read_name_to_id.find(read_name) == read_name_to_id.end()) {
        aln_read_s aln_read;
        aln_read.name = read_name;
        aln_read.length = read.size();
        read_name_to_id[read_name] = aln_reads.size();
        aln_reads.emplace_back(aln_read);
    } else {
        cerr << format("Read {} is already aligned. Skipping...", read_name) << endl;
        return;
    }
    align_matrix_t D = align_matrix_t(read.size(), align_row_t(gene.size(), 0));
    backtrack_matrix_t B = backtrack_matrix_t(read.size(), backtrack_row_t(gene.size(), INVALID_COORDINATE));
    for (index_t i = 1; i < read.size(); i++) {
        for (index_t j = 1; j < gene.size(); j++) {
            local_aligner(D, B, read, i, j);
        }
    }
    vector<local_alignment_s> loc_alns;
    vector<size_t> opt_chain;
    for (size_t i = 0; i < MAX_MAPPINGS; i++) {
        local_alignment_s loc_aln;
        extract_local_alignment(loc_aln, D, B, read);
        if (loc_alns.back().score < MIN_SCORE) {
            cerr << format("Read {} {}-th align score ({}) is below the MIN_SCORE ({})", read_name, loc_alns.size() -1, loc_alns.back().score, MIN_SCORE) << endl;
            break;
        }
        loc_alns.emplace_back(loc_aln);
        recalc_alignment_matrix(D, B, loc_alns.back(), read);
    }
    D.clear();
    B.clear();
    for (local_alignment_s& loc_aln : loc_alns) {
        compress_align_path(loc_aln);
    }
    cochain_mappings(opt_chain, loc_alns);
    extend_opt_chain(loc_alns, opt_chain, read);
    update_dag(read_name, loc_alns, opt_chain);
    print_last_read_alignments(loc_alns);
    loc_alns.clear();
    opt_chain.clear();
}

void dag_aligner::print_last_read_alignments(const vector<local_alignment_s> &loc_alns) {
    stringstream ss;
    seq_id_t read_id = aln_reads.size() - 1;

    for (size_t i = 0; i < loc_alns.size(); i++) {
        ss << format("{}\t", aln_reads[read_id].name); //  Query sequence name
        ss << format("{:d}\t", aln_reads[read_id].length); //  Query sequence length
        ss << format("{:d}\t", loc_alns[i].read_interval.first-1); //  Query start (0-based)
        ss << format("{:d}\t", loc_alns[i].read_interval.second-1); //  Query end (0-based)
        ss << format("{}\t", "+"); //  Relative strand: "+" or "-"
        ss << format("{}\t", gene_name); //  Target sequence name
        ss << format("{:d}\t", gene.size()-1); //  Target sequence length
        stringstream gene_starts;
        stringstream gene_ends;
        stringstream aln_blk_len;
        for (const interval_t& exon : loc_alns[i].gene_intervals) {
            if (exon == loc_alns[i].gene_intervals.front()) {
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
        index_t match_cnt = 0;
        for (const CIGAR_OP& cigar_op : loc_alns[i].cigar) {
            match_cnt += (cigar_op == C_MAT);
        }
        ss << format("{:d}\t", match_cnt); //  Number of residue matches
        ss << format("{}\t", aln_blk_len.str()); //  Alignment block length
        ss << format("{}\t", 255); //  Mapping quality (0-255; 255 for missing)
        ss << format("s1:i:{:d}\t", loc_alns[i].score); //  Local alignment score
        ss << format("oc:c:{:d}\t", loc_alns[i].in_opt_chain);
        ss << format("cg:Z:{}\n", loc_alns[i].cigar_str);
    }
    cout << ss.str();
}
