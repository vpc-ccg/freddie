#include "dag_align.h"

#include <stack>

#define GENE_ID 0

using namespace std;

void process_gene_test() {
    sequence_t gene_test = "ACTGCCGC";
    sequence_list_t reads_test;

    process_gene(gene_test, reads_test);
}

void process_gene(sequence_t gene, sequence_list_t reads) {
    in_neighbors_t parents;
    out_neighbors_t children;
    node_to_base_t node_to_base;
    node_to_read_t node_to_read;


    cout << gene << endl;
    node_id_t last_node = add_node(parents, children, node_to_base, node_to_read, gene[0], GENE_ID);
    node_id_t current_node = last_node;
    for (size_t i = 1; i < gene.size(); i++) {
        last_node = current_node;
        current_node = add_node(parents, children, node_to_base, node_to_read, gene[i], GENE_ID);
        add_edge(parents, children, last_node, current_node);
    }

    current_node = add_node(parents, children, node_to_base, node_to_read, 'X', GENE_ID);
    add_edge(parents, children, current_node, 2);

    node_list_t topo_sorted = topological_sort(parents.size(), parents);
    print_graph(topo_sorted, parents, children, node_to_base, node_to_read);
}

void add_edge(in_neighbors_t &in_neighbors,
        out_neighbors_t &out_neighbors,
        node_id_t source,
        node_id_t target) {
    out_neighbors[source].push_back(target);
    in_neighbors[target].push_back(source);
}


node_id_t add_node(in_neighbors_t &in_neighbors,
        out_neighbors_t &out_neighbors,
        node_to_base_t &node_to_base,
        node_to_read_t &node_to_read,
        nucleotide_t base,
        read_id_t read_id) {
    in_neighbors.push_back(node_list_t());
    out_neighbors.push_back(node_list_t());
    node_to_base.push_back(base);
    node_to_read.push_back(read_id);
    return in_neighbors.size() - 1;
}


void print_graph(const node_list_t node_list,
        const in_neighbors_t &in_neighbors,
        const out_neighbors_t &out_neighbors,
        const node_to_base_t &node_to_base,
        const node_to_read_t &node_to_read)  {
    cout << "N [" << node_list.size() << "]" << endl;
    for (node_id_t node : node_list) {
        cout << node << " (" << node_to_base[node] << ")" << endl;
        cout << "\tIN ["<< in_neighbors[node].size() <<"]: ";
        for (node_id_t neighbor : in_neighbors[node]) {
            cout << neighbor << " (" << node_to_base[neighbor] << "), ";
        }
        cout << endl;
        cout << "\tOUT ["<< out_neighbors[node].size() <<"]: ";
        for (node_id_t neighbor : out_neighbors[node]) {
            cout << neighbor << " (" << node_to_base[neighbor] << "), ";
        }
        cout << endl;
    }
}

node_list_t topological_sort(const node_id_t node_count, const in_neighbors_t& in_neighbors) {
    vector<bool> visited(node_count, false);
    node_list_t result;
    for (node_id_t node = 0; node < node_count; node++) {
        if (visited[node] == true) {
            continue;
        }
        stack<node_id_t> opened;
        opened.push(node);
        while (opened.size() > 0) {
            node_id_t current_node = opened.top();
            size_t unvisited_parents = 0;
            for (node_id_t parent : in_neighbors[current_node]) {
                if (visited[parent] == false) {
                    unvisited_parents++;
                    opened.push(parent);
                }
            }
            if (unvisited_parents == 0) {
                result.push_back(current_node);
                visited[current_node] = true;
                opened.pop();
            }
        }
    }
    return result;
}










//
