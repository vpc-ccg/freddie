#include "commandline_flags.h"

#include "fmt/format.h"
#include <iostream> // std::cout
#include <unordered_set>

using std::cout;
using std::endl;
using std::string;
using std::unordered_set;

string globals::filenames::reads_fasta = "";
string globals::filenames::gene_fasta = "";
string globals::filenames::output_prefix = "";
string globals::filenames::data = "";
int globals::program::save = -1;
int globals::program::generate_dot = -1;
int globals::program::generate_paf = -1;
int globals::program::align = -1;
int globals::program::load = -1;

const unordered_set<string> false_vals({"false", "f", "0"});
const unordered_set<string> true_vals({"true", "t", "1"});

void process_flags() {
    globals::program::load = globals::filenames::data != "" ? 1 : 0;
    globals::program::align = globals::filenames::reads_fasta != "" ? 1 : 0;

    if (!((globals::program::load == 1) ^ (globals::filenames::gene_fasta != ""))) {
        cout << "Provide either just gene fasta or just data file to load!\n";
        print_help();
        exit(-1);
    }
    if (globals::program::generate_paf == 1 && globals::program::align != 1) {
        cout << "Cannot generate PAF file without input reads FASTA file!\n";
        print_help();
        exit(-1);
    }
    if (globals::program::generate_dot == 1 && globals::program::load != 1 && globals::program::align != 1) {
        cout << "Cannot generate DOT file without loading saved DATA file or input reads FASTA file!\n";
        print_help();
        exit(-1);
    }

}

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
        if ((globals::program::save == -1) && (current_param == "-s" || current_param == "--save")) {
            globals::program::save = 1;
            if (i+1 >= argc) {
                continue;
            }
            string val = string(argv[i+1]);
            if (val[0] == '-') {
                continue;
            }
            if (false_vals.find(val) != false_vals.end()) {
                globals::program::save = 0;
                i++;
                continue;
            }
            if (true_vals.find(val) != true_vals.end()) {
                i++;
                continue;
            }
        }
        if ((globals::program::generate_dot == -1) && (current_param == "-d" || current_param == "--generate-dot")) {
            globals::program::generate_dot = 1;
            if (i+1 >= argc) {
                continue;
            }
            string val = string(argv[i+1]);
            if (val[0] == '-') {
                continue;
            }
            if (false_vals.find(val) != false_vals.end()) {
                globals::program::generate_dot = 0;
                i++;
                continue;
            }
            if (true_vals.find(val) != true_vals.end()) {
                i++;
                continue;
            }
        }
        if ((globals::program::generate_paf == -1) && (current_param == "-p" || current_param == "--generate-paf")) {
            globals::program::generate_paf = 1;
            if (i+1 >= argc) {
                continue;
            }
            string val = string(argv[i+1]);
            if (val[0] == '-') {
                continue;
            }
            if (false_vals.find(val) != false_vals.end()) {
                globals::program::generate_paf = 0;
                i++;
                continue;
            }
            if (true_vals.find(val) != true_vals.end()) {
                i++;
                continue;
            }
        }
        cout << "Unrecognized parameter or repeated parameter: " << current_param << "\n";
        print_help();
        abort();
    }
    process_flags();
}

void print_flags(){
    cout << "Parameters:" << endl;
    std::cout << fmt::format("\t{}:\t{}", "reads_fasta", globals::filenames::reads_fasta) << endl;
    std::cout << fmt::format("\t{}:\t{}", "gene_fasta", globals::filenames::gene_fasta) << endl;
    std::cout << fmt::format("\t{}:\t{}", "output_prefix", globals::filenames::output_prefix) << endl;
    std::cout << fmt::format("\t{}:\t{}", "data", globals::filenames::data) << endl;
    std::cout << fmt::format("\t{}:\t{}", "save", globals::program::save) << endl;
    std::cout << fmt::format("\t{}:\t{}", "generate_dot", globals::program::generate_dot) << endl;
    std::cout << fmt::format("\t{}:\t{}", "generate_paf", globals::program::generate_paf) << endl;
    std::cout << fmt::format("\t{}:\t{}", "align", globals::program::align) << endl;
    std::cout << fmt::format("\t{}:\t{}", "load", globals::program::load) << endl;
}

void print_help(){
    cout << "Freddie: Co-chaining of local alignments of long reads on gene DAG" << "\n";
    cout << "Usage: freddie [--PARAMETER VALUE]" << "\n";
    cout << "Example: freddie -r reads.fasta -g gene.fasta -o my_out." << "\n";
    cout << "Example: freddie -r reads.fasta -l gene.data -o my_out." << "\n";
    cout << "Freddie's paramters arguments:" << "\n";
    cout << "    -r    --reads-fasta               (type: string;   OPTIONAL paramter needed to perform alignment)\n";
    cout << "    -g    --gene-fasta                (type: string;   REQUIRED paramter if -l is not provided)\n";
    cout << "    -l    --load                      (type: string;   REQUIRED paramter if -g is not provided)\n";
    cout << "    -o    --output-prefix             (type: string;   Default: \"\")\n";
    cout << "    -s    --save                      (type: bool;     Default: false)\n";
    cout << "    -d    --generate_dot              (type: bool;     Default: false)\n";
    cout << "    -p    --generate_paf              (type: bool;     Default: false)\n";
    cout << "    -h    --help" << endl;
}
