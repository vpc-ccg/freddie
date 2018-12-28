#include "dag_align.h"

#include <fstream>  // std::ofstream
#include <sstream>
#include <stack>
#include <stdlib.h> //abort()

#define EXONIC_ID 0x7FFFFFFF
#define INTRONIC_ID 0x7FFFFFFE
#define MATCH_S 1
#define GAP_S -1
#define MISMATCH_S -1

using namespace std;

constexpr matrix_coordinate_t INVALID = {-1,-1};
void process_gene_test() {
    sequence_list_t reads(1, "ACGTACGGCGTATTGCACGT");
    sequence_t gene = "GTACGCAC";
    vector<bool> exonic(gene.size(), true);
    process_gene(reads, gene, exonic);
}

void process_gene(const sequence_list_t& reads,
                  const sequence_t& gene,
                  const vector<bool>& exonic) {
    in_neighbors_t parents;
    out_neighbors_t children;
    node_to_reads_t node_to_read;

    for (size_t i = 0; i < gene.size(); i++) {
        append_node(parents, children, node_to_read);
    }

    generate_dot(children, gene, exonic, node_to_read, "gene.dot");
    cout << gene << endl;
    for (size_t i = 0; i < reads.size(); i++) {
        cout << reads[i] << endl;
        align_matrix_t D;
        backtrack_matrix_t B;
        local_alignment(reads[i], gene, exonic, parents, D, B);
        print_matrix(reads[i], gene, D, B);
        read_exonic_alignment_t opt_alignment;
        align_score_t opt_score;
        extract_local_alignment(opt_alignment, opt_score, D, B);
        print_read_exonic_alignment(opt_alignment, opt_score, reads[i], gene);
    }

}

align_score_t match(char i, char j, bool is_exonic) {
    if (i == j) {
        return MATCH_S + MATCH_S*is_exonic;
    } else {
        return MISMATCH_S;
    }
}

void local_alignment(const sequence_t& read,
                     const sequence_t& gene,
                     const vector<bool>& exonic,
                     const out_neighbors_t& parents,
                     align_matrix_t& D,
                     backtrack_matrix_t& B) {
    size_t read_l = read.size();
    size_t gene_l = gene.size();
    D.clear();
    B.clear();
    D = align_matrix_t(read_l + 1, align_row_t(gene_l + 1, 0));
    B = backtrack_matrix_t(read_l + 1, backtrack_row_t(gene_l + 1, matrix_coordinate_t(-1,-1)));

    for (size_t i = 1; i < read_l+1; i++) {
        for (size_t j = 1; j < gene_l+1; j++) {
            align_score_t opt_s = 0;
            matrix_coordinate_t opt_b(-1,-1);

            align_score_t cur_s;
            matrix_coordinate_t cur_b;
            // direct parent matching
            cur_b = make_pair(i-1, j-1);
            cur_s  = D[i-1][j-1] + match(read[i-1], gene[j-1], exonic[j-1]);
            if (cur_s > opt_s) {
                opt_s = cur_s;
                opt_b = cur_b;
            }
            // direct parent deletion
            cur_b = make_pair(i-1, j);
            cur_s  = D[i-1][j] + GAP_S;
            if (cur_s > opt_s) {
                opt_s = cur_s;
                opt_b = cur_b;
            }
            // direct parent insertion
            cur_b = make_pair(i, j-1);
            cur_s  = D[i][j-1] + GAP_S;
            if (cur_s > opt_s) {
                opt_s = cur_s;
                opt_b = cur_b;
            }
            for (node_id_t parent : parents[j-1]) {
                size_t j = parent + 1;
                // direct parent matching
                cur_b = make_pair(i-1, j-1);
                cur_s  = D[i-1][j-1] + match(read[i-1], gene[j-1], exonic[j-1]);
                if (cur_s > opt_s) {
                    opt_s = cur_s;
                    opt_b = cur_b;
                }
                // direct parent deletion
                cur_b = make_pair(i-1, j);
                cur_s  = D[i-1][j] + GAP_S;
                if (cur_s > opt_s) {
                    opt_s = cur_s;
                    opt_b = cur_b;
                }
                // direct parent insertion
                cur_b = make_pair(i, j-1);
                cur_s  = D[i][j-1] + GAP_S;
                if (cur_s > opt_s) {
                    opt_s = cur_s;
                    opt_b = cur_b;
                }
            }
            D[i][j] = opt_s;
            B[i][j] = opt_b;
        }
    }
}

void clear_decendants(matrix_coordinate_t source,
                      align_matrix_t D,
                      backtrack_matrix_t B) {
    return;
}

void extract_local_alignment(read_exonic_alignment_t& opt_alignment,
                             align_score_t& opt_score,
                             align_matrix_t& D,
                             backtrack_matrix_t& B) {
    //Create aliases to read fragment (std::pair) and opt_exon_fragments (std::vector of std::pairs)
    read_fragment_t & opt_read_fragment = opt_alignment.first;
    exonic_fragments_t & opt_exon_fragments = opt_alignment.second;
    opt_exon_fragments.clear();

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

    // We know only the tail of the read & tail fragment
    opt_read_fragment = {0,  opt_tail.first - 1};
    exonic_fragment_t cur_opt_exon_frag = {0, opt_tail.second - 1};

    matrix_coordinate_t cur_pos = opt_tail;
    matrix_coordinate_t nxt_pos = B[cur_pos.first][cur_pos.second];
    while (D[nxt_pos.first][nxt_pos.second] > 0) {
        // Opt alignment uses a non trivial edge in the gene DAG
        if (cur_pos.second - nxt_pos.second > 1) {
            cur_opt_exon_frag.first = cur_pos.second - 1;
            opt_exon_fragments.emplace_back(cur_opt_exon_frag);
            cur_opt_exon_frag = exonic_fragment_t(0, nxt_pos.second);
        }
        cur_pos = nxt_pos;
        nxt_pos = B[cur_pos.first][cur_pos.second];
    }
    opt_read_fragment.first = cur_pos.first - 1;
    cur_opt_exon_frag.first = cur_pos.second - 1;
    opt_exon_fragments.emplace_back(cur_opt_exon_frag);

    clear_decendants(cur_pos, D, B);
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
                  const vector<bool>& exonic,
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

void print_read_exonic_alignment(const read_exonic_alignment_t& opt_alignment,
                                 const align_score_t& opt_score,
                                 const sequence_t& read,
                                 const sequence_t& gene) {
    size_t start, end, length;
    start = opt_alignment.first.first;
    end = opt_alignment.first.second;
    length = end - start + 1;

    cout << "Read: " << read.substr(start, length) <<" ("<< start <<","<< end <<")" << endl;
    cout << "Gene: ";
    for (auto & fragment : opt_alignment.second) {
        start = fragment.first;
        end = fragment.second;
        length = end - start + 1;
        cout << gene.substr(start, length) << " ("<< start <<","<< end <<")-";
    }
    cout << endl;
    cout << "Alignment score: " << opt_score << endl;

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
