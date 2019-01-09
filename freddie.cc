#include "commandline_flags.h"
#include "dag_align.h"

using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
    // parse_flags(argc, argv);
    if (argc > 1) {
        cerr << argv[0] << endl;
    }
    // print_flags();
    process_gene_test();
    return 0;
}
