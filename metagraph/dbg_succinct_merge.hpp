#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>
#include <string>


class DBG_succ;
class Config;

namespace merge {

    DBG_succ* build_chunk(const std::vector<const DBG_succ*> &graphs,
                          size_t chunk_idx,
                          size_t num_chunks,
                          size_t num_threads,
                          size_t num_bins_per_thread);

    DBG_succ* merge_chunks(const std::string &filenamebase, size_t num_chunks);

    /*
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    DBG_succ* merge(const std::vector<const DBG_succ*> &Gv);

    /**
     * Heavily borrowing from the graph sequence traversal, this function gets a graph pointer |mergeable| and merges its
     * nodes into the target graph object |target|. The edges of |mergeable| are fully traversed and nodes are added to
     * G_t if not existing yet. This function is well suited to merge small graphs into large ones.
     */
    void merge(DBG_succ *target, const DBG_succ &mergeable);

} // namespace merge

#endif // __MERGE_HPP__
