#include "commandline_flags.h"

#include "fmt/format.h"
#include <iostream> // std::cout, std::cerr

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using fmt::format;

using namespace globals;

string filenames::gene_fasta = "";
string filenames::reads_fasta = "";
string filenames::reads_paf = "";
string filenames::transcript_tsv = "";
string filenames::sim_read_tsv = "";
int program::align = -1;
int program::plot = -1;

void validate_flags() {
    if (program::align == 1 && program::plot == -1) {
        if (filenames::gene_fasta == "") {
            cerr << format("Error: --gene-fasta/-g is required for align command") << endl;
            print_help();
            exit(-1);
        }
        if (filenames::reads_fasta == "") {
            cerr << format("Error: --reads-fasta/-r is required for align command") << endl;
            print_help();
            exit(-1);
        }
    } else if (program::align == -1 && program::plot == 1) {
        if (filenames::reads_paf == "") {
            cerr << format("Error: --reads-paf/-p is required for plot command") << endl;
            print_help();
            exit(-1);
        }
    } else {
        cerr << format("Error: Invalid command") << endl;
        print_help();
        exit(-1);
    }
}

void parse_flags(int argc, char *argv[]){
    if (argc < 2) {
        print_help();
        exit(-1);
    }
    string command = string(argv[1]);
    if (command == "align") {
        program::align = 1;
    }
    else if (command == "plot") {
        program::plot = 1;
    } else {
        cerr << format("Error: Invalid command: {}", command) << endl;
        print_help();
        exit(-1);
    }
    for (int i = 2; i < argc; i++) {
        string current_param(argv[i]);
        if (program::align == 1) {
            if ((filenames::gene_fasta == "") && (current_param == "-g" || current_param == "--gene-fasta")) {
                filenames::gene_fasta = string(argv[i+1]);
                i++;
                continue;
            }
            if ((filenames::reads_fasta == "") && (current_param == "-r" || current_param == "--reads-fasta")) {
                filenames::reads_fasta = string(argv[i+1]);
                i++;
                continue;
            }
            if ((filenames::reads_paf == "") && (current_param == "-p" || current_param == "--reads-paf")) {
                filenames::reads_paf = string(argv[i+1]);
                i++;
                continue;
            }
        }
        if (program::plot == 1) {
            if ((filenames::reads_paf == "") && (current_param == "-p" || current_param == "--reads-paf")) {
                filenames::reads_paf = string(argv[i+1]);
                i++;
                continue;
            }
            if ((filenames::transcript_tsv == "") && (current_param == "-a" || current_param == "--transcript-tsv")) {
                filenames::transcript_tsv = string(argv[i+1]);
                i++;
                continue;
            }
            if ((filenames::sim_read_tsv == "") && (current_param == "-s" || current_param == "--sim-read-tsv")) {
                filenames::sim_read_tsv = string(argv[i+1]);
                i++;
                continue;
            }
        }
        cerr << "Error: An unrecognized or repeated parameter: " << current_param << endl;
        print_help();
        exit(-1);
    }
    validate_flags();
}

void print_flags(){
    cerr << "Globals:" << endl;
    cerr << format("\t{}:\t{}", "align", program::align) << endl;
    cerr << format("\t{}:\t{}", "plot", program::plot) << endl;
    cerr << format("\t{}:\t{}", "gene_fasta", filenames::gene_fasta) << endl;
    cerr << format("\t{}:\t{}", "reads_fasta", filenames::reads_fasta) << endl;
    cerr << format("\t{}:\t{}", "reads_paf", filenames::reads_paf) << endl;
    cerr << format("\t{}:\t{}", "transcript_tsv", filenames::transcript_tsv) << endl;
    cerr << format("\t{}:\t{}", "read_tsv", filenames::sim_read_tsv) << endl;
}

void print_help(){
    cerr << "Freddie: Co-chaining of local alignments of long reads on gene DAG" << endl;
    cerr << "Usage: freddie COMMAND [--PARAMETER VALUE]" << endl;
    cerr << "Example: freddie align -r reads.fasta -g gene.fasta > reads.paf" << endl;
    cerr << "Example: freddie plot -p reads.paf -a transcripts.tsv > reads.dot" << endl;
    cerr << "Freddie's commands:" << endl;
    cerr << "    align                        (Align reads on a gene; stdout a PAF-like file)" << endl;
    cerr << "    plot                         (Takes a PAF-like file; stdout a DOT file)" << endl;
    cerr << "Freddie's align:" << endl;
    cerr << "    -g    --gene-fasta           (type: string;   REQUIRED paramter)" << endl;
    cerr << "    -r    --reads-fasta          (type: string;   REQUIRED paramter)" << endl;
    cerr << "    -p    --reads-paf            (type: string;   OPTIONAL paramter. Loads previous read alignments. Reads from -r file found here will be skipped)" << endl;
    cerr << "Freddie's plot:" << endl;
    cerr << "    -p    --reads-paf            (type: string;   REQUIRED paramter)" << endl;
    cerr << "    -a    --transcript-tsv       (type: string;   OPTIONAL parameter. Annotates the DAG visualization of transcripts)" << endl;
    cerr << "    -s    --sim-read-tsv         (type: string;   OPTIONAL parameter. Annotates the DAG visualization of simulated reads)" << endl;
}
