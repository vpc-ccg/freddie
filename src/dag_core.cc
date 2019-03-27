#include "dag.h"
#include "commandline_flags.h"
#include "utils.h"
#include "fmt/format.h"

#include <sstream>  // std::stringstream
#include <algorithm> // std::sort
#include <cstdlib> // abort()

using namespace dag_types;


using fmt::format;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::vector;
using std::sort;


void parse_tsv_line(const string& line, const long& gene_size, annot_s& annot) {
    vector<string> columns = split(line, '\t');
    if (columns.size() != 4) {
        cerr << "Error: Incorrect number of columns at:" << endl << line << endl;
        abort();
    }
    annot.name = columns[0];
    // chr = columns[1];
    // strand = columns[2];
    for (const string& s : split(columns[3], ',')) {
        if (s.size() == 0) {
            continue;
        }
        vector<string> interval = split(s, '-');
        if (interval.size() != 2) {
            cerr << "Error: Malformated interval at:" << endl << line << endl;
            abort();
        }
        int start = stoi(interval[0]) + 1;
        int end = stoi(interval[1]);
        if (start < 1 || end < 1 || start > gene_size || end > gene_size) {
            cerr << "Error: Malformated interval at:" << endl << line << endl;
            abort();
        }
        // transcript_junctions.insert((index_t) start);
        // transcript_junctions.insert((index_t) end);
        annot.intervals.push_back(interval_t(start, end));
    }
}


void dag_aligner::clear_read_structures(){
    D.clear();
    B.clear();
    align_scores.clear();
    align_paths.clear();
    local_mappings.clear();
    opt_chain_indicator.clear();
    opt_chain.clear();
}

index_t dag_aligner::append_node() {
    if (parents.size() != children.size()) {
        cerr << "Error: Graph corrupted. parents.size()!=children.size()" << endl;
        abort();
    }
    if (node_to_reads.size() != children.size()) {
        cerr << "Error:Graph corrupted. node_to_reads.size()!=children.size()" << endl;
        abort();
    }
    index_t new_node = children.size();
    children.push_back(node_set_t());
    parents.push_back(node_set_t());
    node_to_reads.push_back(read_id_list_t());
    return new_node;
}

void dag_aligner::add_edge(const index_t& source, const index_t& target) {
    if (source >= target) {
        cerr << format("Error: Source can't be >= target: {} -> {}", source, target) << endl;
        return;
    }
    if (target >= parents.size()) {
        cerr << format("Error: Target can't be >= parents.size(): {} -> {}", target, parents.size()) << endl;
        abort();
    }
    if (target - 1 == source) {
        return;
    }
    children[source].insert(target);
    parents[target].insert(source);
}

void dag_aligner::update_dag() {
    if (opt_chain.size() == 0) {
        return;
    }
    vector<interval_t> exons;
    for (const size_t& mapping_id : opt_chain) {
        for (const interval_t& gene_interval : local_mappings[mapping_id].second) {
            exons.push_back(gene_interval);
        }
    }
    for (const interval_t& exon : exons) {
        for (index_t node = exon.first; node <= exon.second; node++) {
            node_to_reads[node].push_back(read_id);
        }
    }
    for (size_t i = 1; i < exons.size(); i++) {
        edge_t e(exons[i-1].second, exons[i-0].first);
        add_edge(e.first, e.second);
        edge_to_reads[e].push_back(read_id);
    }
}

//// Public functions
void dag_aligner::init_dag(const string& gene_name_in, const string& gene_in) {
    init_dag(gene_name_in, gene_in, vector<bool>(gene_in.size(), false));
}



void dag_aligner::init_dag(const string& gene_name_in, const string& gene_in, const vector<bool>& exonic_indicator_in) {
    //// Clearing structures
    // Alignment stuff
    read_id = -1;
    children.clear();
    parents.clear();
    node_to_reads.clear();
    // Annotaion stuff
    t_annot_junctions.clear();
    t_annots.clear();
    r_annot_junctions.clear();
    r_annots.clear();

    // Initializing required structures
    gene = gene_in;
    gene_name = gene_name_in;
    exonic_indicator = exonic_indicator_in;
    for (index_t i = 0; i < gene.size(); i++) {
        append_node();
    }

    // Initializing transcripts annotations for plotting
    if (globals::filenames::transcript_tsv != "") {
        ifstream tsv(globals::filenames::transcript_tsv);
        string line;
        while (getline(tsv, line)) {
            annot_s annot;
            parse_tsv_line(line, (long) gene.size(), annot);
            for (const interval_t& interval : annot.intervals) {
                t_annot_junctions.insert(interval.first);
                t_annot_junctions.insert(interval.second);
            }
            t_annots.emplace_back(annot);
        }
    }
    if (globals::filenames::sim_read_tsv != "") {
        ifstream tsv(globals::filenames::sim_read_tsv);
        string line;
        while (getline(tsv, line)) {
            annot_s annot;
            parse_tsv_line(line, (long) gene.size(), annot);
            for (const interval_t& interval : annot.intervals) {
                r_annot_junctions.insert(interval.first);
                r_annot_junctions.insert(interval.second);
            }
            r_annots.emplace_back(annot);
        }
    }
}

void dag_aligner::load_state(const string& paf_path) {
    ifstream ifs;
    ifs.open(paf_path);
    size_t line_num = -1;
    string buffer, rname, gname;
    rname = "";
    gname = "";
    int rlen = -1;
    int glen = -1;
    vector<interval_t> exons;
    while (true) {
        line_num++;
        bool eof = !getline(ifs, buffer);
        // Empty PAF file
        if (eof && line_num == 0) {
            cerr << format("Warning: empty PAF file: {}", paf_path) << endl;
            break;
        }
        vector<string> fields;
        if (!eof) {
            fields = split(buffer, '\t');
            // Check # of fields
            if (fields.size() < 12) {
                cerr << format("Error: PAF file has invalid number of fields ({}) at: {}:{}", fields.size(), paf_path, line_num) << endl;
                abort();
            }
            // Initialize DAG on first line
            if (line_num == 0) {
                gname = fields[5];
                glen = stoi(fields[6]) + 1;
                rname = fields[0];
                rlen = stoi(fields[1]) + 1;
                if (glen < 1) {
                    cerr << format("Error: Negative gene length at {}:{}", paf_path, line_num) << endl;
                    abort();
                }
                if (gene_name.size() > 0 && gene_name != gname) {
                    cerr << format("Error: Gene name in {} does not match gene name in {}:{}", globals::filenames::gene_fasta, paf_path, line_num) << endl;
                    abort();
                }
                if (gene.size() > 0 && gene.size() != (size_t)glen) {
                    cerr << format("Error: Gene length in {} does not match gene length in {}:{}", globals::filenames::gene_fasta, paf_path, line_num) << endl;
                    abort();
                }
                if (gene_name.size() == 0) {
                    init_dag(gname, string(glen, '^'));
                }
            }
        }
        // Resolve read at end of file or at start of new read
        if (eof || fields[0] != rname) {
            read_names.push_back(rname);
            read_id = read_names.size()-1;
            read_name_to_id[rname] = read_id;
            cerr << format("    {}:{}", rname, read_name_to_id[rname]) << endl;

            sort(exons.begin(), exons.end());

            for (size_t i = 1; i < exons.size(); i++) {
                edge_t e(exons[i-1].second, exons[i-0].first);
                add_edge(e.first, e.second);
                edge_to_reads[e].push_back(read_id);
            }
            for (const interval_t& exon : exons) {
                for (index_t node = exon.first; node <= exon.second; node++) {
                    node_to_reads[node].push_back(read_id);
                }
            }
            aln_read_intervals.emplace_back(move(exons));
            exons.clear();
        }
        if (eof) {
            break;
        }
        if (fields[0] != rname) {
            rname = fields[0];
            if (read_name_to_id.find(rname) != read_name_to_id.end()) {
                cerr << format("Error: PAF file not sorted by read name. Check at {}:{}", paf_path, line_num) << endl;
                abort();
            }
            rlen = stoi(fields[1]) + 1;
            if (rlen < 1) {
                cerr << format("Error: Negative read length at {}:{}", paf_path, line_num) << endl;
                abort();
            }
        }
        if (rlen != stoi(fields[1]) + 1) {
            cerr << format("Error: Read has inconsistent length at {}:{}", paf_path, line_num) << endl;
            abort();
        }
        vector<index_t> starts, ends;
        for (const string& start : split(fields[7], ',')) {
            int val = stoi(start)+1;
            if (val < 1) {
                cerr << format("Error: Negative start ({}) at {}:{}", start, paf_path, line_num) << endl;
                abort();
            }
            starts.push_back(val);
        }
        for (const string& end : split(fields[8], ',')) {
            int val = stoi(end)+1;
            if (val < 1) {
                cerr << format("Error: Negative start ({}) at {}:{}", end, paf_path, line_num) << endl;
                abort();
            }
            ends.push_back(val);
        }
        if (starts.size() != ends.size()) {
            cerr << format("Error: unbalanced number of ends ({}) and starts ({}) at {}:{}", starts.size(), ends.size(), paf_path, line_num) << endl;
        }
        int oc_tag = -1;
        for (size_t i = 12; i < fields.size(); i++) {
            vector<string> vs = split(fields[i], ':');
            if (vs.size() != 3) {
                cerr << format("Error: Invalid tag format at {}:{}", paf_path, line_num) << endl;
                abort();
            }
            if (fields[i].substr(0,5) == "oc:c:") {
                if (oc_tag != -1) {
                    cerr << format("Error: There is more than one oc tag at {}:{}", paf_path, line_num) << endl;
                    abort();
                }
                oc_tag = stoi(vs[2]);
                if (oc_tag != 0 && oc_tag != 1) {
                    cerr << format("Error: Invalid oc tag value ({}) at {}:{}", oc_tag, paf_path, line_num) << endl;
                    abort();
                }
            }
        }
        if (oc_tag == 1) {
            for (size_t i = 0; i < starts.size(); i++) {
                exons.push_back(interval_t(starts[i], ends[i]));
            }
        }
    }
    ifs.close();
}
