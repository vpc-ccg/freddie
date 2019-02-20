#ifndef FREDDIE_DAGALIGN_H
#define FREDDIE_DAGALIGN_H

#include <vector>
#include <unordered_map>
#include <set>
#include <functional> //std::function
#include <utility> //std::pair stuff
#include <string>
#include <iostream>
#include <fstream>  // std::ofstream

namespace dag_types {
    typedef uint32_t read_id_t;
    typedef uint32_t index_t;
    typedef uint16_t tid_t;
    typedef std::set<index_t> node_set_t;
    typedef std::vector<node_set_t> neighbors_t;
    typedef std::vector<read_id_t> read_id_list_t;
    typedef std::pair<index_t, index_t> edge_t;
    struct edge_hash {
        std::size_t operator () (const std::pair<index_t,index_t> &p) const {
            auto h1 = std::hash<index_t>{}(p.first);
            auto h2 = std::hash<index_t>{}(p.second);
            return h1 ^ h2;
        }
    };
    // Alignment matrix
    typedef int32_t align_score_t;
    typedef std::pair<index_t, index_t> matrix_coordinate_t;
    typedef std::vector<align_score_t> align_row_t;
    typedef std::vector<align_row_t> align_matrix_t;
    typedef std::vector<matrix_coordinate_t> backtrack_row_t;
    typedef std::vector<backtrack_row_t> backtrack_matrix_t;
    typedef std::vector<matrix_coordinate_t> align_path_t;
    // Optimal alignment extraction
    typedef uint32_t seq_idx_t;
    typedef std::pair<seq_idx_t, seq_idx_t> interval_t;
    typedef std::pair<interval_t, std::vector<interval_t>> mapping_t;
}

class dag_aligner {
private:
    //// Private members
    // DAG
    std::string gene_name;
    std::string gene;
    std::vector<bool> exonic_indicator;
    std::string read;
    dag_types::read_id_t read_id;
    dag_types::neighbors_t parents;
    dag_types::neighbors_t children;
    std::vector<dag_types::read_id_list_t> node_to_reads;
    std::unordered_map<dag_types::edge_t, dag_types::read_id_list_t, dag_types::edge_hash> edge_to_reads;
    std::vector<std::string> read_names;
    std::unordered_map<std::string, dag_types::read_id_t> read_name_to_id;
    // transcript GTF annotaions
    dag_types::node_set_t transcript_junctions;
    std::vector<std::string> transcripts;
    std::vector<std::vector<dag_types::interval_t>> transcript_intervals;
    // Alignment matrix
    dag_types::align_matrix_t D;
    dag_types::backtrack_matrix_t B;
    std::vector<dag_types::align_score_t> align_scores;
    std::vector<dag_types::align_path_t> align_paths;
    std::vector<dag_types::mapping_t> local_mappings;
    // Co-linear chaining
    std::vector<bool> opt_chain_indicator;
    std::vector<size_t> opt_chain;
    //// Helper functions
    void clear_read_structures();
    void add_edge(const dag_types::index_t& source, const dag_types::index_t& target);
    dag_types::index_t append_node();
    void local_aligner(const dag_types::index_t& i, const dag_types::index_t& j);
    void extract_local_alignment();
    void recalc_alignment_matrix();
    void compress_align_paths();
    void cochain_mappings();
    void update_dag();
public:
    void init_dag(const std::string& gene_name, const std::string& gene);
    void init_dag(const std::string& gene_name, const std::string& gene, const std::vector<bool>& exonic_indicator);
    void align_read(const std::string& read_name, const std::string& read);
    void print_last_read_alignments();
    void generate_dot();
    void load_state(const std::string& file_path);
};
#endif //FREDDIE_DAGALIGN_H
