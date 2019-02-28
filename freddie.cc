#include "commandline_flags.h"
#include "dag_align.h"
#include "utils.h"
#include "fmt/format.h"

#include <iostream>
#include <fstream>
#include <string>

#ifndef GITINFO
#define GITINFO "@ unknown git commit"
#endif

using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using fmt::format;

using namespace globals;

void align(dag_aligner& my_dag) {
    string name, seq;
    cerr << format("Reading {} fasta file...", globals::filenames::gene_fasta) << endl;
    ifstream gene_fasta(globals::filenames::gene_fasta);
    name = "";
    seq = "^";
    if (fasta_get_record(name, seq, gene_fasta) == false) {
        cerr << "Error: empty gene FASTA file." << endl;
        abort();
    }
    cerr << format("Gene {} of length {} is read.", name, seq.size()-1) << endl;
    my_dag.init_dag(name, seq);
    if (fasta_get_record(name, seq, gene_fasta) == true) {
        cerr << "Error: Gene FASTA has more than one record." << endl;
        abort();
    }
    if (filenames::reads_paf != "") {
        my_dag.load_state(format(filenames::reads_paf));
    }

    ifstream reads_fasta(globals::filenames::reads_fasta);
    while (true) {
        name = "";
        seq = "^";
        if (fasta_get_record(name, seq, reads_fasta) == false) {
            break;
        }
        my_dag.align_read(name, seq);
        my_dag.print_last_read_alignments();
    }
}

void plot(dag_aligner& my_dag) {
    my_dag.load_state(format(filenames::reads_paf));
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
