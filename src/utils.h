#ifndef FREDDIE_UTILS_H
#define FREDDIE_UTILS_H

#include <string>
#include <vector>

std::vector<std::string> split(const std::string &s, char delim);

bool fasta_get_record(std::string& name, std::string& seq, std::ifstream& fasta_file);

#endif //FREDDIE_UTILS_H
