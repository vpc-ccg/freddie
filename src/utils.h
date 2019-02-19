#ifndef FREDDIE_UTILS_H
#define FREDDIE_UTILS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

std::vector<string> split(const std::string &s, char delim) {
    std::vector<string> result;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}


bool fasta_get_record(string& name, std::string& seq, std::ifstream& fasta_file) {
    char c;
    if (fasta_file.get(c) == false) {
        return false;
    }
    if (c != '>') {
        std::cerr << "Error: not fasta format." << std::endl;
        abort();
    }
    if (std::getline(fasta_file, name)) {
        std::cerr << "Error: not fasta format." << std::endl;
        abort();
    }
    while(fasta_file.get(c)) {
        if (c == '>') {
            break;
        }
        if (c == '\n') {
            continue;
        }
        if (c < 'A' || c > 'z') {
            std::cerr << "Error: not fasta format." << std::endl;
            abort();
        }
        seq += c;
    }
    return true;
}

#endif //FREDDIE_UTILS_H
