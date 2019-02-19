#include "commandline_flags.h"
#include "dag_align.h"
#include "fmt/format.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>  // log10

using std::string;
using std::vector;
using std::ofstream;
using fmt::format;
using namespace globals;

void align(dag_aligner& my_dag) {
    ifstream fasta (globals::filenames::gene_fasta);
    string name;
    string seq = "^";
    if (fasta_get_record(name, seq, fasta) == false) {
        cerr << "Error: empty gene FASTA file." << endl;
    }

    my_dag.init_dag(gene_seq, gene_name);
    if (filenames::reads_paf != "") {
        my_dag.load_state(format(filenames::reads_paf));
    }

    fasta = ifstream(globals::filenames::reads_fasta);
    while (fasta_get_record(name, seq, reads_file)) {
        my_dag.align_read(name, seq);
        my_dag.print_last_read_alignments();
    }
}

void plot(dag_aligner& my_dag) {
    my_dag.generate_dot();
}

int main(int argc, char *argv[]) {
    cerr << format("🕺 Freddie 🕺: {}", string(GITINFO)) << endl;
    parse_flags(argc, argv);
    print_flags();

    dag_aligner my_dag;
    if (program::align == 1) {
        align(my_dag);
    }
    if (program::plot == 1) {
        plot(my_dag);
    }

    return 0;
}
