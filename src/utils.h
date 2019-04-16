#ifndef FREDDIE_UTILS_H
#define FREDDIE_UTILS_H

#include <string>
#include <vector>

namespace utils {
    std::vector<std::string> split(const std::string &s, char delim);

    bool fasta_get_record(std::string& name, std::string& seq, std::ifstream& fasta_file);

    template <class var_t, class val_t>
    void set_to_max (var_t& a, val_t& a_val, const var_t& b, const val_t& b_val) {
        if (a_val < b_val) {
            a = b;
            a_val = b_val;
        }
    }
}
#endif //FREDDIE_UTILS_H
