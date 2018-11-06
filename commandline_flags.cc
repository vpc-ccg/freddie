#include "commandline_flags.h"

using namespace std;
// Parameter definitions
string input_sam = "";
string input_genes_fasta = "";
string output_prefix = "";

void parse_flags(int argc, char *argv[]){
    for (int i = 1; i < argc; i++) {
        string current_param(argv[i]);
        if (current_param == "-h" || current_param == "--help") {
            print_help();
            exit(0);
        }
        if ((input_sam == "") && (current_param == "-r" || current_param == "--input-reads-sam")) {
            input_sam = string(argv[i+1]);
            i++;
            continue;
        }
        if ((input_genes_fasta == "") && (current_param == "-g" || current_param == "--input-gene-fasta")) {
            input_genes_fasta = string(argv[i+1]);
            i++;
            continue;
        }
        if ((output_prefix == "") && (current_param == "-o" || current_param == "--output-prefix")) {
            output_prefix = string(argv[i+1]);
            i++;
            continue;
        }

        cout << "Unrecognized parameter or repeated parameter: " << current_param << "\n";
        print_help();
        exit(-1);
    }
    if (input_sam == "" || input_genes_fasta == "") {
        cout << "Missing input or output files parameters!\n";
        print_help();
        exit(-1);
    }
}

void print_flags(){
    cout << "Parameters:\n";
    cout << "\tinput_reads_sam:\t" << input_sam << "\n";
    cout << "\tinput_genes_fasta:\t" << input_genes_fasta << "\n";
    cout << "\toutput_prefix:\t" << output_prefix << "\n";
    cout << "\n";

}

void print_help(){
    cout << "Freddi: Co-chaining of local alignments of long reads on gene DAGs" << "\n";
	cout << "Usage: freddie [--PARAMETER VALUE]" << "\n";
	cout << "Example: freddie -r reads.sam -r R2.fastq -o my_out. -e 1 -l 8 -m 5 -t 2 -k 4 --silent" << "\n";
	cout << "Calib's paramters arguments:" << "\n";
    cout << "\t-r\t--input-reads-sam               \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-g\t--input-genes-fasta              \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-o\t--output-prefix                 \t(type: string;   Default \"\")\n";
    cout << "\t-h\t--help\n";
}
