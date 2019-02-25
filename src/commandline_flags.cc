#include "commandline_flags.h"

#include "fmt/format.h"
#include <iostream> // std::cout, std::cerr

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using namespace globals;

string filenames::gene_fasta = "";
string filenames::reads_fasta = "";
string filenames::reads_paf = "";
string filenames::transcript_tsv = "";
int program::align = -1;
int program::plot = -1;

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
        }
        cerr << "Error: An unrecognized or repeated parameter: " << current_param << "\n";
        print_help();
        exit(-1);
    }
}

void print_flags(){
    cerr << "Globals:" << endl;
    cerr << fmt::format("\t{}:\t{}", "align", program::align) << endl;
    cerr << fmt::format("\t{}:\t{}", "plot", program::plot) << endl;
    cerr << fmt::format("\t{}:\t{}", "gene_fasta", filenames::gene_fasta) << endl;
    cerr << fmt::format("\t{}:\t{}", "reads_fasta", filenames::reads_fasta) << endl;
    cerr << fmt::format("\t{}:\t{}", "reads_paf", filenames::reads_paf) << endl;
    cerr << fmt::format("\t{}:\t{}", "transcript_tsv", filenames::transcript_tsv) << endl;
}

void print_help(){
    cerr << "Freddie: Co-chaining of local alignments of long reads on gene DAG" << "\n";
    cerr << "Usage: freddie COMMAND [--PARAMETER VALUE]" << "\n";
    cerr << "Example: freddie align -r reads.fasta -g gene.fasta > reads.paf" << "\n";
    cerr << "Example: freddie plot -p reads.paf -a transcripts.tsv > reads.dot" << "\n";
    cerr << "Freddie's commands:" << "\n";
    cerr << "    align                        (Align reads on a gene; stdout a PAF-like file)" << "\n";
    cerr << "    plot                         (Takes a PAF-like file; stdout a DOT file)" << "\n";
    cerr << "Freddie's align:" << "\n";
    cerr << "    -g    --gene-fasta           (type: string;   REQUIRED paramter)" << "\n";
    cerr << "    -r    --reads-fasta          (type: string;   REQUIRED paramter)" << "\n";
    cerr << "    -p    --reads-paf            (type: string;   OPTIONAL paramter. Loads previous read alignments. Reads from -r file found here will be skipped)" << "\n";
    cerr << "Freddie's plot:" << "\n";
    cerr << "    -p    --reads-paf            (type: string;   REQUIRED paramter)" << "\n";
    cerr << "    -a    --transcript-tsv       (type: string;   OPTIONAL parameter. Annotates the DAG visualization)" << "\n";
}
