#include "commandline_flags.h"

using namespace std;

int main(int argc, char *argv[]) {
    parse_flags(argc, argv);

    print_flags(); 

    return 0;
}
