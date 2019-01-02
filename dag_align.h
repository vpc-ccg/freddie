#ifndef FREDDIE_DAGALIGN_H
#define FREDDIE_DAGALIGN_H

#include "global.h"
#include <vector>
#include <utility> //std::pair stuff

typedef uint32_t read_id_t;
typedef uint32_t node_id_t;
typedef uint32_t read_pos_t;

typedef int32_t align_score_t;
typedef uint16_t backtrack_t;
typedef std::pair<backtrack_t, backtrack_t> matrix_coordinate_t;

typedef std::pair<read_pos_t, read_pos_t> read_fragment_t;
typedef std::pair<node_id_t, node_id_t> exonic_fragment_t;
typedef std::vector<exonic_fragment_t> exonic_fragments_t;
typedef std::pair<read_fragment_t, exonic_fragments_t> read_exonic_alignment_t;
typedef std::vector<read_exonic_alignment_t> read_exonic_alignments_t;

typedef std::vector<node_id_t> node_list_t;
typedef std::vector<node_list_t> out_neighbors_t;
typedef std::vector<node_list_t> in_neighbors_t;
typedef std::vector<read_id_t> read_list_t;
typedef std::vector<read_list_t> node_to_reads_t;


typedef std::string sequence_t;
typedef std::vector<sequence_t> sequence_list_t;

typedef std::vector<align_score_t> align_row_t;
typedef std::vector<align_row_t> align_matrix_t;
typedef std::vector<matrix_coordinate_t> backtrack_row_t;
typedef std::vector<backtrack_row_t> backtrack_matrix_t;

void process_gene_test();

void process_gene(const sequence_list_t& reads,
                  const sequence_t& gene,
                  const std::vector<bool>& exonic);

void add_edge(in_neighbors_t &in_neighbors,
              out_neighbors_t &out_neighbors,
              node_id_t source,
              node_id_t target);

node_id_t append_node(in_neighbors_t &in_neighbors,
                      out_neighbors_t &out_neighbors,
                      node_to_reads_t &node_to_read);

align_score_t match(char i, char j, bool is_exonic);

void local_alignment(const sequence_t& read,
                     const sequence_t& gene,
                     const std::vector<bool>& exonic,
                     const out_neighbors_t& parents,
                     align_matrix_t& D,
                     backtrack_matrix_t& B);

void clear_descendants(matrix_coordinate_t source,
                       const in_neighbors_t& children,
                       align_matrix_t& D,
                       backtrack_matrix_t& B);

void extract_local_alignment(read_exonic_alignment_t& opt_alignment,
                             const in_neighbors_t& children,
                             align_score_t& opt_score,
                             align_matrix_t& D,
                             backtrack_matrix_t& B);

void print_matrix(const sequence_t& read,
                  const sequence_t& gene,
                  const align_matrix_t& D,
                  const backtrack_matrix_t& B);

void print_read_exonic_alignment(const read_exonic_alignment_t& opt_alignment,
                                 const align_score_t& opt_score,
                                 const sequence_t& read,
                                 const sequence_t& gene);

void generate_dot(const node_to_reads_t& node_to_read,
                  const sequence_t& gene,
                  const std::vector<bool>& exonic,
                  const out_neighbors_t& children,
                  const std::string output_path);

#endif //FREDDIE_DAGALIGN_H
