//
//  path_encoder_baseline.cpp
//  lossless_dbg
//
//  Created by Jan Studený on 20/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#include "path_encoder.hpp"
#include "utils.hpp"
#include <iostream>
#include <set>
#include <map>

using namespace std;


// todo find tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this
using node_index = DeBruijnGraph::node_index;
using path_id = pair<node_index,int>;

class PathDatabaseBaseline : public PathDatabase<path_id> {
public:
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaseline(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabase(graph), graph(*graph_) {
        
    }
    
    node_index starting_node(const string& sequence) {
        return graph.kmer_to_node(sequence.substr(0,graph.get_k()));
    }
    
    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering
        vector<path_id> encoded;
        for(auto& sequence : sequences) {
            // add additional bifurcation
            additional_joins.insert(starting_node(sequence));
            additional_splits.insert(graph.kmer_to_node(sequence.substr(sequence.length()-graph.get_k())));
        }
        
        for(auto& sequence : sequences) {
            int relative_order = route_sequence(sequence);
            encoded.push_back({starting_node(sequence),relative_order});
        }
        
        return encoded;
    }
    
    int route_sequence(const string& sequence) {
        auto kmer = sequence.substr(0,graph.get_k());
        auto node = graph.kmer_to_node(kmer);
        // always putting new read above all other reads
        int relative_position = offset_for_symbol(joins[node],'~');
        int relative_starting_position = joins[node]['~'];
        joins[node]['~']++;

        int kmer_position = 0;
        for(auto& base : sequence.substr(graph.get_k())) {
            if (node_is_split(node)) {
                auto& routing_table = splits[node];
                auto rt_index = routing_table.begin();
                advance(rt_index,relative_position);
                routing_table.insert(rt_index,base);
                relative_position = rank(routing_table,base,relative_position)-1;
            }
            node = graph.traverse(node,base);
            kmer_position++;
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += offset_for_symbol(joins[node],join_symbol);
                joins[node][join_symbol]++;
            }
        }
        
        auto& routing_table = splits[node];
        auto rt_index = routing_table.begin();
        advance(rt_index,relative_position);
        routing_table.insert(rt_index,'~');
        
        return relative_starting_position;
    }
    
    // returns the number of reads that were already stored and should have lower index
    int offset_for_symbol(const map<char,int>& join,char symbol) const {
        int result = 0;
        for(auto&[base,count] : join) {
            result += count;
            if (base == symbol) break;
        }
        return result;
    }
    
    bool node_is_join(node_index node) const {
        return !graph.is_single_incoming(node) or additional_joins.count(node);
    }
    bool node_is_split(node_index node) const {
        return !graph.is_single_outgoing(node) or additional_splits.count(node);
    }
    
    int rank(const routing_table_t& routing_table, char symbol, int position) const {
        assert(position < routing_table.size());
        int result = 0;
        int i = 0;
        for(auto it=begin(routing_table);i<=position;it++,i++) {
            result += *it == symbol;
        }
        return result;
    }
    
    void compress() {
        
    }
    
    size_t num_paths() const override;
    
    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;
        
        int relative_starting_position = path.second;
        int relative_position = offset_for_symbol(joins.at(node),'~')
                                - joins.at(node).at('~')
                                + relative_starting_position;
        
        int kmer_position = 0;
        char base;
        while(true) {
            if (node_is_split(node)) {
                auto& routing_table = splits.at(node);
                auto rt_index = routing_table.begin();
                advance(rt_index,relative_position);
                base = *rt_index;
                relative_position = rank(routing_table,base,relative_position)-1;
            }
            else {
                assert(graph.is_single_outgoing(node));
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
            }
            if (base == '~') break;
            node = graph.traverse(node,base);
            kmer_position++;
            sequence.append(1,base); // 1 times base
            
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += offset_for_symbol(joins.at(node),join_symbol);
            }
        }
        
        return sequence;
    }
    
    std::vector<path_id> get_paths_going_through(const std::string &str) const override;
    
    std::vector<path_id> get_paths_going_through(PathDatabase::node_index node) const override;
    
    PathDatabase::node_index get_next_node(PathDatabase::node_index node, path_id path) const override;
    
    PathDatabase::node_index get_next_consistent_node(const std::string &history) const override;
    
private:
    // denote how many reads are joining from every branch (ATCGN~) (~ denotes start of a new read)
    std::map<node_index,map<char,int>> joins;
    // denote where the reads should go (ATCGN~) (~ denodes the end of particular read)
    
    std::map<node_index,routing_table_t> splits;
    
    std::set<node_index> additional_joins;
    std::set<node_index> additional_splits;
    
    const DeBruijnGraph & graph;
    
};


int main() {

}
