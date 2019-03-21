#include "utils.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using std::string;
using std::vector;
using std::stringstream;
using std::getline;
using std::cerr;
using std::endl;
using std::ifstream;

vector<string> split(const string &s, char delim) {
    vector<string> result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

bool fasta_get_record(string& name, string& seq, ifstream& fasta_file) {
    char c;
    fasta_file.get(c);
    if (fasta_file.eof()) {
        return false;
    }
    if (!getline(fasta_file, name)) {
        cerr << "Error: not fasta format; orphan '>' char" << endl;
        abort();
    }
    if (c!='>') {
        name = c + name;
    }
    while(true) {
        fasta_file.get(c);
        if (fasta_file.eof() || c == '>') {
            break;
        }
        if (c == '\n') {
            continue;
        }
        if (c < 'A' || c > 'z') {
            cerr << "Error: not fasta format; found non alphabetical char" << endl;
            abort();
        }
        seq += c;
    }
    name = name.substr(0, name.find(" ", 0));
    return true;
}
