#ifndef FREDDIE_COMMANDLINE_H
#define FREDDIE_COMMANDLINE_H

#include <string>
#include <unordered_map>

namespace globals {
    struct filenames {
        static std::string reads_fasta;
        static std::string gene_fasta;
        static std::string output_prefix;
        static std::string data;
    };
}

void parse_flags(int argc, char *argv[]);

void print_flags();

void print_help();

#endif //FREDDIE_COMMANDLINE_H
