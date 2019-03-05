#ifndef FREDDIE_COMMANDLINE_H
#define FREDDIE_COMMANDLINE_H

#include <string>
#include <unordered_map>

namespace globals {
    struct filenames {
        static std::string reads_fasta;
        static std::string gene_fasta;
        static std::string reads_paf;
        static std::string transcript_tsv;
        static std::string sim_read_tsv;
    };
    struct program {
        static int align;
        static int plot;
    };
}

void parse_flags(int argc, char *argv[]);

void print_flags();

void print_help();

#endif //FREDDIE_COMMANDLINE_H
