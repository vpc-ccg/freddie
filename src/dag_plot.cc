#include "dag.h"
#include "commandline_flags.h"
#include "utils.h"
#include "fmt/format.h"

#include <sstream>  // std::stringstream
#include <unordered_map>
#include <map>
#include <algorithm> // std::reverse, std::sort
#include <functional> // std::function
#include <cstdlib> // abort()
#include <utility> // move()
#include <iterator> // next(), std::inserter

using namespace dag_types;

constexpr double DOT_CANVAS_WIDTH = 100.0; //in inches
constexpr size_t MAX_PADDING = 7500;
constexpr interval_t INVALID_INTERVAL = interval_t(-1,-1);
using fmt::format;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::unordered_map;
using std::map;

const vector<string> DOT_COLORS = { // http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
};

// A function and its helpers to generate a .DOT file of the DAG

void dag_aligner::generate_dot() {
    stringstream ss;
    ss << "digraph graphname{" << endl;
    ss << "    rankdir=LR;" << endl;
    vector<index_t> junctions;
    map<index_t,size_t> junction_idx;
    // Getting a vector of all DOT vertices
    for (index_t node = 1; node < nodes.size(); node++) {
        bool flag = false;
        string comment = "";
        if (t_annot_junctions.find(node) != t_annot_junctions.end()) {
            flag = true;
            comment += "T";
        }
        if (r_annot_junctions.find(node) != r_annot_junctions.end()) {
            flag = true;
            comment += "R";
        }
        if (nodes[node].read_ids.size() > 1 && nodes[node-1].read_ids.size() == 0) {
            flag = true;
            comment += "1";
        }
        if (node + 1 < nodes.size() && nodes[node].read_ids.size() > 0 && nodes[node+1].read_ids.size() == 0) {
            flag = true;
            comment += "1";
        }
        if (nodes[node].children.size() > 0) {
            flag = true;
            comment += "C";
        }
        if (nodes[node].parents.size() > 0) {
            flag = true;
            comment += "P";
        }
        if (aln_junctions.find(node) != aln_junctions.end()) {
            if (!flag) {
                cerr << "Warning" << endl;
            }
            flag = true;
            comment += "J";
        }
        if (node == nodes.size() -1) {
            flag = true;
            comment += "E";
        }
        if (flag) {
            ss << format("    {:d} [label=\"{:d}:{:d}\"] //{:d}-{}", node, node-1, nodes[node].read_ids.size(), junctions.size(), comment) << endl;
            junction_idx[node] = junctions.size();
            junctions.push_back(node);
        }
    }
    for (const auto& kv : junction_idx) {
        cerr << format("N: {}; idx: {};", kv.first, kv.second) << endl;
    }
    // Adding edges of normal coverage between vertices
    ss << format("    edge[weight=1000, arrowhead=none];") << endl;
    for (size_t i = 1; i < junctions.size(); i++) {
        index_t start  = junctions[i-1];
        index_t end  = junctions[i];
        double len = end-start-1;
        if (len == 0) {
            ss << format("    {}->{}", start, end) << endl;
        } else {
            double coverage = 0.0;
            for (index_t node = start + 1; node < end; node++) {
                coverage += nodes[node].read_ids.size();
            }
            size_t padding = (end-start)/2;
            padding = padding > MAX_PADDING ? MAX_PADDING : padding;
            string space_padding(padding, ' ');
            ss << format("    {}->{} [label=\"{}{:.2f}{}\"]", start, end, space_padding, coverage/len, space_padding) << endl;
        }
    }
    // Adding edges of jumping coverage between vertices
    ss << format("    edge[weight=1, arrowhead=normal];") << endl;
    for (size_t node = 1; node < nodes.size() - 1; node++) {
        for (index_t child : nodes[node].children) {
            edge_t e(node, child);
            stringstream comment;
            for (const index_t& rid : edge_to_reads[e]) {
                comment << format("{},", rid);
            }
            ss << format("    {}->{} [label={:d}] // {}", node, child, edge_to_reads[e].size(), comment.str()) << endl;
        }
    }
    // Adding transcript_tsv edges
    ss << format("    edge[weight=1, arrowhead=none];") << endl;
    string transcript_dot_format = "    {}->{} [style={} color=\"{}\" label=\"t{:d}{}{:d}\"]";
    for (seq_id_t tid = 0; tid < t_annots.size(); tid++) {
        for (size_t i = 0; i < t_annots[tid].intervals.size(); i++) {
            index_t start = t_annots[tid].intervals[i].first;
            index_t end   = t_annots[tid].intervals[i].second;
            for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                index_t source = junctions[j];
                index_t target = junctions[j+1];
                ss << format(transcript_dot_format, source, target, "bold", DOT_COLORS[tid%DOT_COLORS.size()], tid, "e", i) << endl;
            }
            if (i+1 == t_annots[tid].intervals.size()) {
                continue;
            }
            index_t next_start = t_annots[tid].intervals[i+1].first;
            ss << format(transcript_dot_format, end, next_start, "dotted", DOT_COLORS[tid%DOT_COLORS.size()], tid, "i", i) << endl;
        }
    }
    // Adding sim_read_tsv edges
    ss << format("    edge[weight=1, arrowhead=none];") << endl;
    string sim_read_dot_format = "    {}->{} [style={} color=\"{}\" label=\"R{:d}{}{:d}\"]";
    for (seq_id_t rid = 0; rid < r_annots.size(); rid++) {
        for (size_t i = 0; i < r_annots[rid].intervals.size(); i++) {
            index_t start = r_annots[rid].intervals[i].first;
            index_t end   = r_annots[rid].intervals[i].second;
            for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                index_t source = junctions[j];
                index_t target = junctions[j+1];
                ss << format(sim_read_dot_format, source, target, "bold", DOT_COLORS[rid%DOT_COLORS.size()], rid, "E", i) << endl;
            }
            if (i+1 == r_annots[rid].intervals.size()) {
                continue;
            }
            index_t next_start = r_annots[rid].intervals[i+1].first;
            ss << format(sim_read_dot_format, end, next_start, "dotted", DOT_COLORS[rid%DOT_COLORS.size()], rid, "I", i) << endl;
        }
    }
    // Adding sim_read_tsv aligned correspondence edges
    ss << format("    edge[weight=1, arrowhead=none];") << endl;
    string aln_read_dot_format = "    {}->{} [style={} color=\"{}\" label=\"r{:d}{}{:d}\"]";
    for (const annot_s& r_annot : r_annots) {
        if (read_name_to_id.find(r_annot.name)==read_name_to_id.end()){
            continue;
        }
        const seq_id_t& rid = read_name_to_id[r_annot.name];
        cerr << format("rid: {}; name: {}", rid, r_annot.name) << endl;
        interval_t prev_exon = INVALID_INTERVAL;
        size_t eid = 0;
        for (const mapping_s& mapping : aln_reads[rid].mappings) {
            for (const interval_t& exon : mapping.gene_intervals) {
                cerr << format("--> eid: {}; s:{}; e:{}; ps:{}; pe:{}", eid, exon.first, exon.second, prev_exon.first, prev_exon.second) << endl;
                if (prev_exon != INVALID_INTERVAL) {
                    index_t source = prev_exon.second;
                    index_t target = exon.first;
                    ss << format(aln_read_dot_format, source, target, "dashed", DOT_COLORS[rid%DOT_COLORS.size()], rid, "i", eid) << endl;
                }
                index_t start = exon.first;
                index_t end   = exon.second;
                cerr << format("--> --> js: {}; je:{}", junction_idx[start], junction_idx[end]) << endl;
                for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                    index_t source = junctions[j];
                    index_t target = junctions[j+1];
                    cerr << format("--> --> --> j: {}; s:{}; t:{}", j, source, target) << endl;
                    ss << format(aln_read_dot_format, source, target, "solid", DOT_COLORS[rid%DOT_COLORS.size()], rid, "e", eid) << endl;
                }
                prev_exon = exon;
                eid++;
            }
        }
    }
    ss << "}" << endl;
    cout << ss.str();
}

void dag_aligner::print_dag() {
    for (index_t node = 0; node < nodes.size(); node++) {
        if (nodes[node].children.size() == 0) {
            continue;
        }
        cerr << node << " -> ";
        for (const index_t& child : nodes[node].children) {
            edge_t e(node, child);
            cerr << format("{}({}) ", child, edge_to_reads[e].size());
        }
        cerr << endl;
    }
}
