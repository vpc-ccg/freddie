#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>


#ifndef FREDDIE_COMMANDLINE_H
#define FREDDIE_COMMANDLINE_H
#include <stdint.h>

typedef uint32_t node_id_t;

extern std::string input_sam;
extern std::string input_genes_fasta;
extern std::string output_prefix;


void parse_flags(int argc, char *argv[]);
void print_flags();
void print_help();

#endif //FREDDIE_COMMANDLINE_H
