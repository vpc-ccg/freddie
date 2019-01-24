#include "commandline_flags.h"

#include <iostream> // std::cout
#include "fmt/format.h"

using std::cout;
using std::endl;
using std::string;

string globals::filenames::reads_fasta = "";
string globals::filenames::gene_fasta = "";
string globals::filenames::output_prefix = "";
string globals::filenames::data = "";

void parse_flags(int argc, char *argv[]){
    for (int i = 1; i < argc; i++) {
        string current_param(argv[i]);
        if (current_param == "-h" || current_param == "--help") {
            print_help();
            exit(0);
        }
        if ((globals::filenames::reads_fasta == "") && (current_param == "-r" || current_param == "--reads-fasta")) {
            globals::filenames::reads_fasta = string(argv[i+1]);
            i++;
            continue;
        }
        if ((globals::filenames::gene_fasta == "") && (current_param == "-g" || current_param == "--gene-fasta")) {
            globals::filenames::gene_fasta = string(argv[i+1]);
            i++;
            continue;
        }
        if ((globals::filenames::output_prefix == "") && (current_param == "-o" || current_param == "--output-prefix")) {
            globals::filenames::output_prefix = string(argv[i+1]);
            i++;
            continue;
        }
        if ((globals::filenames::data == "") && (current_param == "-l" || current_param == "--load")) {
            globals::filenames::data = string(argv[i+1]);
            i++;
            continue;
        }
        cout << "Unrecognized parameter or repeated parameter: " << current_param << "\n";
        print_help();
        exit(-1);
    }
    if (!(globals::filenames::data == "") ^ (globals::filenames::gene_fasta == "")) {
        cout << "Provide either just gene fasta or just data file to load!\n";
        print_help();
        exit(-1);
    }
    if (globals::filenames::reads_fasta == "") {
        cout << "Missing reads input parameter!\n";
        print_help();
        exit(-1);
    }
}

void print_flags(){
    cout << "Parameters:" << endl;
    std::cout << fmt::format("\t{}:\t{}", "reads_fasta", globals::filenames::reads_fasta) << endl;
    std::cout << fmt::format("\t{}:\t{}", "gene_fasta", globals::filenames::gene_fasta) << endl;
    std::cout << fmt::format("\t{}:\t{}", "output_prefix", globals::filenames::output_prefix) << endl;
    std::cout << fmt::format("\t{}:\t{}", "data", globals::filenames::data) << endl;
}

void print_help(){
    cout << "Freddi: Co-chaining of local alignments of long reads on gene DAGs" << "\n";
    cout << "Usage: freddie [--PARAMETER VALUE]" << "\n";
    cout << "Example: freddie -r reads.fasta -g gene.fasta -o my_out." << "\n";
    cout << "Example: freddie -r reads.fasta -l gene.data -o my_out." << "\n";
    cout << "Calib's paramters arguments:" << "\n";
    cout << "\t-r\t--reads-fasta               \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-g\t--gene-fasta                \t(type: string;   REQUIRED paramter if -l is not provided)\n";
    cout << "\t-l\t--load                      \t(type: string;   REQUIRED paramter if -g is not provided)\n";
    cout << "\t-o\t--output-prefix             \t(type: string;   Default \"\")\n";
    cout << "\t-h\t--help\n";
}
