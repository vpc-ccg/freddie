#include "dag_align.h"
#include "commandline_flags.h"
#include "fmt/format.h"

#include <sstream>  // std::stringstream
#include <queue>
#include <algorithm> // std::reverse
#include <functional> // std::function
#include <cstdlib> // abort()

#define TEST_GENE "AANNAANNCCNNGG"
#define TEST_READ2 "AAAGG"
#define TEST_READ1 "AAAACC"

using namespace dag_types;

constexpr align_score_t MATCH_S = 1;
constexpr align_score_t GAP_S = -1;
constexpr align_score_t MISMATCH_S = -1;
constexpr size_t MAX_MAPPINGS = 10;
constexpr align_score_t MIN_SCORE = 1;
constexpr matrix_coordinate_t INVALID_COORDINATE = {-1,-1};

using fmt::format;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ofstream;
using std::queue;
using std::vector;
using std::function;
using std::reverse;

void dag_aligner::clear_read_structures(){
    D.clear();
    B.clear();
    align_scores.clear();
    align_paths.clear();
    local_mappings.clear();
    opt_chain_indicator.clear();
    opt_chain.clear();
}

index_t dag_aligner::append_node() {
    if (parents.size() != children.size()) {
        cout << "ERR:Graph corrupted. parents.size()!=children.size()" << endl;
        abort();
    }
    if (node_to_read.size() != children.size()) {
        cout << "ERR:Graph corrupted. node_to_read.size()!=children.size()" << endl;
        abort();
    }
    index_t new_node = children.size();
    children.push_back(node_set_t());
    parents.push_back(node_set_t());
    node_to_read.push_back(read_id_list_t());
    return new_node;
}

void dag_aligner::add_edge(const index_t& source, const index_t& target) {
    if (source >= target) {
        cout << format("Source can't be >= target: {} -> {}", source, target) << endl;
        abort();
    }
    if (target >= parents.size()) {
        cout << format("Target can't be >= parents.size(): {} -> {}", target, parents.size()) << endl;
        abort();
    }
    if (target - 1 == source) {
        return;
    }
    children[source].insert(target);
    parents[target].insert(source);
}

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
        matching_score += MATCH_S + MATCH_S*this->exonic_indicator[j];
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

void dag_aligner::update_dag() {
    vector<interval_t> exons;
    for (const size_t& mapping_id : opt_chain) {
        for (const interval_t& gene_interval : local_mappings[mapping_id].second) {
            exons.push_back(gene_interval);
        }
    }
    for (const interval_t& exon : exons) {
        for (index_t node = exon.first; node <= exon.second; node++) {
            node_to_read[node].push_back(read_id);
        }
    }
    for (size_t i = 1; i < exons.size(); i++) {
        index_t source = exons[i-1].second;
        index_t target = exons[i-0].first;
        add_edge(source, target);
    }
}

void dag_aligner::init_dag(const string& gene) {
    init_dag(gene, vector<bool>(gene.size(), false));
}

void dag_aligner::init_dag(const string& gene_in, const std::vector<bool>& exonic_indicator_in) {
    read_id = 0;
    children.clear();
    parents.clear();
    node_to_read.clear();
    gene = gene_in;
    exonic_indicator = exonic_indicator_in;
    for (index_t i = 0; i < gene.size(); i++) {
        append_node();
    }
}

void dag_aligner::align_read(const string& read_in) {
    clear_read_structures();
    read = read_in;
    read_id++;
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
        print_matrix(format("{}{}_{}.mat", globals::filenames::output_prefix, read_id, i));
        recalc_alignment_matrix();
    }
    compress_align_paths();
    cochain_mappings();
    update_dag();
}

// A function to generate a .DOT file of the DAG
void dag_aligner::generate_dot(const string& output_path) {
    stringstream output;
    output << "digraph graphname{" << endl;
    output << "    rankdir=LR;" << endl;
    for (index_t node = 0; node < children.size(); node++) {
        string outline = "peripheries=";
        outline +=  exonic_indicator[node] ? "2" : "1";
        output << "    " << node <<" [label=\"" << node << ":" << gene[node] << ":" << node_to_read[node].size() << "\" "<< outline <<"]" << endl;
    }
    output << endl;
    for (index_t node = 0; node < children.size() - 1; node++) {
        output << "    " << node <<"->" << node + 1 << endl;
    }
    for (index_t node = 0; node < children.size(); node++) {
        for (index_t child : children[node]) {
            output << "    " << node <<"->" << child << endl;
        }
    }
    output << "}" << endl;
    ofstream ofile;
    ofile.open(output_path);
    ofile << output.str();
    ofile.close();
}

void dag_aligner::print_last_read_to_paf(ofstream& out_file) {
    stringstream ss;
    for (size_t i = 0; i < align_paths.size(); i++) {
        ss << format("{}\t", read_id); //  Query sequence name
        ss << format("{}\t", read.size()); //  Query sequence length
        ss << format("{}\t", local_mappings[i].first.first); //  Query start (0-based)
        ss << format("{}\t", local_mappings[i].first.second); //  Query end (0-based)
        ss << format("{}\t", "+"); //  Relative strand: "+" or "-"
        ss << format("{}\t", "gene"); //  Target sequence name
        ss << format("{}\t", gene.size()); //  Target sequence length
        stringstream gene_starts;
        stringstream gene_ends;
        for (const interval_t& exon : local_mappings[i].second) {
            if (exon == local_mappings[i].second[0]) {
                gene_starts << format("{}", exon.first);
                gene_ends << format("{}", exon.second);
            } else {
                gene_starts << format(",{}", exon.first);
                gene_ends << format(",{}", exon.second);
            }
        }
        ss << format("{}\t", gene_starts.str());  //  Target start on original strand (0-based)
        ss << format("{}\t", gene_ends.str());  //  Target end on original strand (0-based)
        ss << format("{}\t", align_scores[i]); //  Local alignment score
        // ss << format("{}\t", 0); //  Number of residue matches
        // ss << format("{}\t", 0); //  Alignment block length
        // ss << format("{}\t", 0); //  Mapping quality (0-255; 255 for missing)
        ss << format("tg:c:{:d}\n", opt_chain_indicator[i]);
    }
    out_file << ss.str();
    out_file.flush();
}

// Printing here is formatted to (kinda) work with "column -t" Linux program
void dag_aligner::print_matrix(const string& output_path){
    stringstream ss;
    ss << "\\\t";
    ss << format("0(^)");
    for (index_t j = 1; j < D[0].size(); j++) {
        ss << format("\t{}({})", j, gene[j-1]);
    }
    ss << "\n";
    ss << "0(^)";
    for (index_t j = 0; j < D[0].size(); j++) {
        if (B[0][j] == INVALID_COORDINATE) {
            ss << format("\t{}({},{})", D[0][j], "x", "x");
        } else {
            ss << format("\t{}({},{})", D[0][j], B[0][j].first, B[0][j].second);
        }
    }
    ss << "\n";
    for (index_t i = 1; i < D.size(); i++) {
        ss << format("{}({})", i, read[i-1]);
        for (index_t j = 0; j < D[i].size(); j++) {
            if (B[i][j] == INVALID_COORDINATE) {
                ss << format("\t{}({},{})", D[i][j], "x", "x");
            }
            else {
                ss << format("\t{}({},{})", D[i][j], B[i][j].first, B[i][j].second);
            }
        }
        ss << "\n";
    }
    ofstream ofile;
    ofile.open(output_path);
    ofile << ss.str();
    ofile.close();
}

//
// void print_cochain(const read_gene_mappings_t& chain) {
//     for (read_gene_mapping_t mapping : chain) {
//         read_interval_t& read_interval = mapping.first;
//         cout << "R" << "("<<read_interval.first<<","<<read_interval.second<<") "<<endl;
//         for (gene_interval_t gene_interval : mapping.second) {
//             cout << "G" << "("<<gene_interval.first<<","<<gene_interval.second<<")"<<'-';
//         }
//         cout << endl;
//     }
// }
//
// void print_mapping_interval(const size_t& interval_id) {
//     size_t start, end, length;
//
//     start = read_gene_mapping.first.first;
//     end = read_gene_mapping.first.second;
//     length = end - start + 1;
//     cout << score << ": ";
//     cout << read.substr(start-1, length) <<"("<<start<<","<<end<<")" << " --> ";
//     for (auto & fragment : read_gene_mapping.second) {
//         start = fragment.first;
//         end = fragment.second;
//         length = end - start + 1;
//         cout << gene.substr(start-1, length) <<"("<<start<<","<<end<<")-";
//     }
//     cout << endl;
// }
//
