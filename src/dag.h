#ifndef FREDDIE_DAG_H
#define FREDDIE_DAG_H

#include <vector>
#include <unordered_map>
#include <set>
#include <utility> //std::pair stuff
#include <string>
#include <iostream>
#include <fstream>  // std::ofstream, std::ifstream
#include <limits>   // numeric_limits

namespace dag_types {
    template<typename T>
    class min_init {
        T val;
    public:
        min_init() : val(static_cast<T>(std::numeric_limits<T>::min())) { }
        min_init(T val) : val(val) { }
        operator T&() { return val; }
        operator T() const { return val; }

    };
    template<typename T>
    class invalid_init {
        T val;
    public:
        invalid_init() : val(-1) { }
        invalid_init(T val) : val(val) { }
        operator T&() { return val; }
        operator T() const { return val; }
    };

    enum CIGAR_OP {
        C_MAT = '=',
        C_MIS = 'X',
        C_INS = 'I',
        C_DEL = 'D'
    };
    typedef invalid_init<uint32_t> seq_id_t;
    typedef min_init<uint32_t> index_t;
    // Node
    typedef std::set<index_t> node_set_t;
    typedef std::vector<seq_id_t> read_id_list_t;
    struct node_s {
        node_set_t parents;
        node_set_t children;
        read_id_list_t read_ids;
    };
    // Edge
    typedef std::pair<index_t, index_t> edge_t;
    struct edge_hash {
        std::size_t operator () (const std::pair<index_t,index_t> &p) const {
            auto h1 = std::hash<uint32_t>{}((uint32_t)p.first);
            auto h2 = std::hash<uint32_t>{}((uint32_t)p.second);
            return h1 ^ h2;
        }
    };
    // Aligned reads
    typedef min_init<int32_t> align_score_t;
    typedef std::pair<index_t, index_t> interval_t;
    struct mapping_s {
        interval_t read_interval;
        std::vector<interval_t> gene_intervals;
        align_score_t score = 0;
        std::string cigar_str = "";
    };
    struct aln_read_s {
        std::string name = "";
        index_t length = 0;
        std::vector<mapping_s> mappings;
    };
    // Annotations
    struct annot_s {
        std::string name;
        std::string chrom;
        std::vector<dag_types::interval_t> intervals;
    };

}

class dag_aligner {
private:
    //// Private members
    // DAG
    std::string gene_name;
    std::string gene;
    std::vector<bool> exonic_indicator;
    std::vector<dag_types::node_s> nodes;
    std::vector<dag_types::aln_read_s> aln_reads;
    dag_types::node_set_t aln_junctions;
    std::unordered_map<dag_types::edge_t, dag_types::read_id_list_t, dag_types::edge_hash> edge_to_reads;
    std::unordered_map<std::string, dag_types::seq_id_t> read_name_to_id;
    // Annotaions
    dag_types::node_set_t t_annot_junctions;
    std::vector<dag_types::annot_s> t_annots;
    std::unordered_map<std::string, dag_types::seq_id_t> t_annot_name_to_id;
    dag_types::node_set_t r_annot_junctions;
    std::vector<dag_types::annot_s> r_annots;
    std::unordered_map<std::string, std::string> r_annot_name_to_t_annot_name;
    //// Helper functions
    void add_edge(const dag_types::index_t& source, const dag_types::index_t& target);
    void local_aligner(const dag_types::index_t& i, const dag_types::index_t& j);
    void extract_local_alignment();
    void recalc_alignment_matrix();
    void compress_align_paths();
    void cochain_mappings();
    void extend_opt_chain();
    void update_dag(const std::string& read_name_in);
public:
    void init_dag(const std::string& gene_name, const std::string& gene);
    void init_dag(const std::string& gene_name, const std::string& gene, const std::vector<bool>& exonic_indicator);
    void align_read(const std::string& read_name, const std::string& read);
    void print_last_read_alignments();
    void generate_dot();
    void load_state(const std::string& paf_path);
    void print_dag();
};
#endif //FREDDIE_DAG_H
