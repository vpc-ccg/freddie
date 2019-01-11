#ifndef FREDDIE_DAGALIGN_H
#define FREDDIE_DAGALIGN_H

#include "global.h"
#include <vector>
#include <functional> //std::function
#include <utility> //std::pair stuff


typedef uint32_t read_id_t;
typedef uint32_t node_id_t;
typedef uint32_t read_pos_t;

typedef int32_t align_score_t;
typedef uint16_t backtrack_t;
typedef std::pair<backtrack_t, backtrack_t> matrix_coordinate_t;

typedef std::vector<matrix_coordinate_t> align_path_t;

typedef std::function<void (backtrack_t, backtrack_t)> local_aligner_func_t;

typedef std::vector<node_id_t> node_list_t;
typedef std::vector<node_list_t> out_neighbors_t;
typedef std::vector<node_list_t> in_neighbors_t;
typedef std::vector<read_id_t> read_list_t;
typedef std::vector<read_list_t> node_to_reads_t;

typedef std::string sequence_t;
typedef std::vector<sequence_t> sequence_list_t;

typedef std::vector<bool> exonic_indicator_t;

typedef std::vector<align_score_t> align_row_t;
typedef std::vector<align_row_t> align_matrix_t;
typedef std::vector<matrix_coordinate_t> backtrack_row_t;
typedef std::vector<backtrack_row_t> backtrack_matrix_t;

typedef std::pair<read_pos_t, read_pos_t> read_interval_t;
typedef std::pair<node_id_t, node_id_t> gene_interval_t;
typedef std::vector<gene_interval_t> gene_intervals_t;
typedef std::pair<read_interval_t, gene_intervals_t> read_gene_mapping_t;
typedef std::vector<read_gene_mapping_t> read_gene_mappings_t;

void process_gene_test();

void process_gene(const sequence_list_t& reads,
                  const sequence_t& gene,
                  const exonic_indicator_t& exonic);

void add_edge(in_neighbors_t &in_neighbors,
              out_neighbors_t &out_neighbors,
              node_id_t source,
              node_id_t target);

node_id_t append_node(in_neighbors_t &in_neighbors,
                      out_neighbors_t &out_neighbors,
                      node_to_reads_t &node_to_read);

align_score_t match(char i, char j, bool is_exonic);

read_gene_mapping_t get_mapping_intervals(const align_path_t& path);

read_gene_mappings_t get_optimal_cochain(const std::vector<align_score_t>& scores,
                                         const read_gene_mappings_t& mappings);

void print_cochain(const read_gene_mappings_t& chain);

void print_mapping_interval(const align_score_t& score,
                            const read_gene_mapping_t& read_gene_mapping,
                            const sequence_t& read,
                            const sequence_t& gene);

void local_alignment(align_matrix_t& D,
                     backtrack_matrix_t& B,
                     const sequence_t& read,
                     const sequence_t& gene,
                     const local_aligner_func_t& local_aligner);

void recalc_alignment_matrix(align_matrix_t& D,
                             backtrack_matrix_t& B,
                             const align_path_t& opt_alignment,
                             const in_neighbors_t& children,
                             const local_aligner_func_t& local_aligner);

void extract_local_alignment(align_path_t& opt_alignment,
                             align_score_t& opt_score,
                             const align_matrix_t& D,
                             const backtrack_matrix_t& B);

local_aligner_func_t get_local_aligner(align_matrix_t& D,
                                       backtrack_matrix_t& B,
                                       const sequence_t& read,
                                       const sequence_t& gene,
                                       const exonic_indicator_t& exonic,
                                       const out_neighbors_t& parents);

void print_matrix(const sequence_t& read,
                  const sequence_t& gene,
                  const align_matrix_t& D,
                  const backtrack_matrix_t& B);

// void print_read_exonic_alignment(const read_exonic_alignment_t& opt_alignment,
//                                  const align_score_t& opt_score,
//                                  const sequence_t& read,
//                                  const sequence_t& gene);

void generate_dot(const node_to_reads_t& node_to_read,
                  const sequence_t& gene,
                  const exonic_indicator_t& exonic,
                  const out_neighbors_t& children,
                  const std::string output_path);

#endif //FREDDIE_DAGALIGN_H
