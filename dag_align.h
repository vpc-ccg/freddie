#ifndef FREDDIE_DAGALIGN_H
#define FREDDIE_DAGALIGN_H
#include "global.h"
#include <vector>
#include <utility> //std::pair stuff

typedef uint32_t read_id_t;
typedef uint32_t node_id_t;
typedef uint8_t nucleotide_t;
typedef int32_t align_score_t;
typedef std::pair<int16_t, int16_t> backtrack_t;

typedef std::vector<node_id_t> node_list_t;
typedef std::vector<node_list_t> out_neighbors_t;
typedef std::vector<node_list_t> in_neighbors_t;
typedef std::vector<read_id_t> read_list_t;
typedef std::vector<read_list_t> node_to_read_t;
typedef std::vector<nucleotide_t> node_to_base_t;

typedef std::string sequence_t;
typedef std::vector<sequence_t> sequence_list_t;

typedef std::vector<align_score_t> align_row_t;
typedef std::vector<align_row_t> align_matrix_t;
typedef std::vector<backtrack_t> backtrack_row_t;
typedef std::vector<backtrack_row_t> backtrack_matrix_t;

void process_gene_test();
void process_gene(sequence_t gene,
                  sequence_list_t reads);

void print_graph(const node_list_t node_list,
                 const in_neighbors_t &in_neighbors,
                 const out_neighbors_t &out_neighbors,
                 const node_to_base_t &node_to_base,
                 const node_to_read_t &node_to_read);

node_list_t topological_sort(const node_id_t node_count,
                             const in_neighbors_t& in_neighbors);
void add_edge(in_neighbors_t &in_neighbors,
              out_neighbors_t &out_neighbors,
              node_id_t source,
              node_id_t target);

node_id_t add_node(in_neighbors_t &in_neighbors,
                   out_neighbors_t &out_neighbors,
                   node_to_base_t &node_to_base,
                   node_to_read_t &node_to_read,
                   nucleotide_t base,
                   read_id_t read_id);

void dag_local_alignment(const in_neighbors_t& parents,
                         const node_to_base_t& node_to_base,
                         align_matrix_t& D,
                         backtrack_matrix_t& B);

void generate_dot(const node_list_t topo_sorted,
                  const out_neighbors_t& children,
                  const node_to_base_t& node_to_base,
                  const node_to_read_t& node_to_read,
                  const std::string output_path);

#endif //FREDDIE_DAGALIGN_H
