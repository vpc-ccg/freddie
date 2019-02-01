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
using std::ofstream;
using std::getline;
using fmt::format;

int main(int argc, char *argv[]) {
    cout << format("🕺 Freddie 🕺: {}", string(GITINFO)) << endl;
    parse_flags(argc, argv);
    print_flags();

    dag_aligner my_dag = dag_aligner();
    if (globals::program::load) {
        my_dag.load_state(format(globals::filenames::data));
    } else {
        ifstream gene_file (globals::filenames::gene_fasta);
        string gene_name;
        getline(gene_file, gene_name);
        gene_name += " ";
        gene_name = gene_name.substr(1, gene_name.find(" ") - 1);
        string gene_seq;
        getline(gene_file, gene_seq);
        cout << fmt::format("The length of the gene is {}", gene_seq.size()) << endl;
        my_dag.init_dag("^" + gene_seq, gene_name);
    }
    if (globals::program::align) {
        ifstream reads_file (globals::filenames::reads_fasta);
        string line;
        vector<string> reads;
        while (getline (reads_file, line)) {
            getline (reads_file, line);
            reads.push_back("^" + line);
        }
        cout << fmt::format("The number of reads is {}", reads.size()) << endl;
        ofstream paf_file;
        if (globals::program::generate_paf) {
            paf_file.open(format("{}dag.paf", globals::filenames::output_prefix), std::ios_base::app);
        }
        for (size_t i = 0; i < reads.size(); i++) {
            my_dag.align_read(reads[i]);
            if (globals::program::generate_paf) {
                my_dag.print_last_read_to_paf(paf_file);
            }
            if (i % 50 == 0) {
                if (globals::program::save) {
                    my_dag.save_state(format("{}dag.data", globals::filenames::output_prefix));
                }
                if (globals::program::generate_dot) {
                    my_dag.generate_compressed_dot(format("{}dag_comp_{}.dot", globals::filenames::output_prefix, i));
                }
                cout << fmt::format("Done with read {}", i) << endl;
            }
        }
        if (globals::program::generate_paf) {
            paf_file.close();
        }
    }


    return 0;
}
