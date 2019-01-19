#include "commandline_flags.h"
#include "dag_align.h"
#include "fmt/format.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::getline;
using fmt::format;

int main(int argc, char *argv[]) {
    parse_flags(argc, argv);
    print_flags();
    ifstream reads_file (globals::filenames::reads_fasta);
    string line;
    vector<string> reads;
    while (getline (reads_file, line)) {
        getline (reads_file, line);
        reads.push_back(line);
    }
    cout << fmt::format("The number of reads is {}", reads.size()) << endl;
    ifstream gene_file (globals::filenames::gene_fasta);
    string gene;
    getline(gene_file, gene);
    getline(gene_file, gene);
    cout << fmt::format("The length of the gene is {}", gene.size()) << endl;

    dag_aligner my_dag = dag_aligner();
    my_dag.init_dag(gene);
    for (size_t i = 0; i < reads.size(); i++) {
        my_dag.align_read(reads[i]);
        my_dag.generate_dot(format("dag_{}.dot", i));
    }

    return 0;
}
