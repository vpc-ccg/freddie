#ifndef FREDDIE_COMMANDLINE_H
#define FREDDIE_COMMANDLINE_H

#include <string>
#include <unordered_map>

namespace globals {
    struct filenames {
        static std::string reads_fasta;
        static std::string gene_fasta;
        static std::string transcript_tsv;
        static std::string output_prefix;
        static std::string data;
    };
    struct program {
        static int align;
        static int load;
        static int save;
        static int generate_dot;
        static int generate_paf;
    };
}

void parse_flags(int argc, char *argv[]);

void print_flags();

void print_help();

#endif //FREDDIE_COMMANDLINE_H
