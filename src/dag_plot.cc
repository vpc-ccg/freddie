#include "dag.h"
#include "commandline_flags.h"
#include "utils.h"
#include "fmt/format.h"

#include <sstream>  // std::stringstream
#include <unordered_map>
#include <algorithm> // std::reverse, std::sort
#include <functional> // std::function
#include <cstdlib> // abort()
#include <utility> // move()
#include <iterator> // next(), std::inserter

using namespace dag_types;

constexpr double DOT_CANVAS_WIDTH = 100.0; //in inches

using fmt::format;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::unordered_map;

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
    unordered_map<index_t,size_t> junction_idx;
    // Getting a vector of all DOT vertices
    for (index_t node = 1; node < children.size(); node++) {
        bool flag = false;
        string comment = "";
        if (transcript_junctions.find(node) != transcript_junctions.end()) {
            flag = true;
            comment += "TS";
        }
        if (sim_read_junctions.find(node) != sim_read_junctions.end()) {
            flag = true;
            comment += "RS";
        }
        if (node_to_reads[node].size() > 1 && node_to_reads[node-1].size() == 0) {
            flag = true;
            comment += "1";
        }
        if (node + 1 < children.size() && node_to_reads[node].size() > 0 && node_to_reads[node+1].size() == 0) {
            flag = true;
            comment += "1";
        }
        if (children[node].size() > 0) {
            flag = true;
            comment += "C";
        }
        if (parents[node].size() > 0) {
            flag = true;
            comment += "P";
        }
        if (node == children.size() -1) {
            flag = true;
            comment += "E";
        }
        if (flag) {
            junction_idx[node] = junctions.size();
            junctions.push_back(node);
            ss << format("    {:d} [label=\"{:d}:{:d}\"] //{}", node, node-1, node_to_reads[node].size(), comment) << endl;
        }
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
                coverage += node_to_reads[node].size();
            }
            string space_padding((end-start)/2, ' ');
            ss << format("    {}->{} [label=\"{}{:.2f}{}\"]", start, end, space_padding, coverage/len, space_padding) << endl;
        }
    }
    // Adding edges of jumping coverage between vertices
    ss << format("    edge[weight=1, arrowhead=normal];") << endl;
    for (size_t node = 1; node < children.size() - 1; node++) {
        for (index_t child : children[node]) {
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
    for (tid_t tid = 0; tid < transcript_intervals.size(); tid++) {
        for (size_t i = 0; i < transcript_intervals[tid].size(); i++) {
            index_t start = transcript_intervals[tid][i].first;
            index_t end = transcript_intervals[tid][i].second;
            for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                index_t source = junctions[j];
                index_t target = junctions[j+1];
                ss << format(transcript_dot_format, source, target, "bold", DOT_COLORS[tid%DOT_COLORS.size()], tid, "e", i) << endl;
            }
            if (i+1 == transcript_intervals[tid].size()) {
                continue;
            }
            index_t next_start = transcript_intervals[tid][i+1].first;
            ss << format(transcript_dot_format, end, next_start, "dotted", DOT_COLORS[tid%DOT_COLORS.size()], tid, "i", i) << endl;
        }
    }
    // Adding sim_read_tsv edges
    ss << format("    edge[weight=1, arrowhead=none];") << endl;
    string sim_read_dot_format = "    {}->{} [style={} color=\"{}\" label=\"sim{:d}{}{:d}\"]";
    string aln_read_dot_format = "    {}->{} [style={} color=\"{}\" label=\"aln{:d}{}{:d}\"]";
    for (index_t sim_rid = 0; sim_rid < sim_read_intervals.size(); sim_rid++) {
        for (size_t i = 0; i < sim_read_intervals[sim_rid].size(); i++) {
            index_t start = sim_read_intervals[sim_rid][i].first;
            index_t end = sim_read_intervals[sim_rid][i].second;
            for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                index_t source = junctions[j];
                index_t target = junctions[j+1];
                ss << format(sim_read_dot_format, source, target, "bold", DOT_COLORS[sim_rid%DOT_COLORS.size()], sim_rid, "e", i) << endl;
            }
            if (i+1 == sim_read_intervals[sim_rid].size()) {
                continue;
            }
            index_t next_start = sim_read_intervals[sim_rid][i+1].first;
            ss << format(sim_read_dot_format, end, next_start, "dotted", DOT_COLORS[sim_rid%DOT_COLORS.size()], sim_rid, "i", i) << endl;
        }
        if (read_name_to_id.find(sim_reads[sim_rid]) == read_name_to_id.end()) {
            cerr << format("Error: Could not find sim_reads[sim_rid] (sim_rid={}) sim_reads in read_name_to_id (size={}) unordered_map", sim_rid, read_name_to_id.size()) << endl;
            for (const auto& kv : read_name_to_id) {
                cerr << format("    {} : {}", kv.first, kv.second) << endl;
            }
            abort();
            continue;
        }
        index_t rid = read_name_to_id[sim_reads[sim_rid]];
        for (size_t i = 0; i < aln_read_intervals[rid].size(); i++) {
            index_t start = sim_read_intervals[rid][i].first;
            index_t end = sim_read_intervals[rid][i].second;
            for (size_t j = junction_idx[start]; j < junction_idx[end]; j++) {
                index_t source = junctions[j];
                index_t target = junctions[j+1];
                ss << format(aln_read_dot_format, source, target, "solid", DOT_COLORS[sim_rid%DOT_COLORS.size()], sim_rid, "e", i) << endl;
            }
            if (i+1 == sim_read_intervals[sim_rid].size()) {
                continue;
            }
            index_t next_start = sim_read_intervals[rid][i+1].first;
            ss << format(aln_read_dot_format, end, next_start, "dashed", DOT_COLORS[sim_rid%DOT_COLORS.size()], sim_rid, "i", i) << endl;
        }
    }
    ss << "}" << endl;
    cout << ss.str();
}

void dag_aligner::print_dag() {
    for (index_t node = 0; node < children.size(); node++) {
        if (children[node].size() == 0) {
            continue;
        }
        cerr << node << " -> ";
        for (const index_t& child : children[node]) {
            edge_t e(node, child);
            cerr << format("{}({}) ", child, edge_to_reads[e].size());
        }
        cerr << endl;
    }
}
