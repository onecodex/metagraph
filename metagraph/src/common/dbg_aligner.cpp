#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "bounded_priority_queue.hpp"
#include "ssw_cpp.h"

DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_top_paths, bool verbose,
                       float sw_threshold, float re_seeding_threshold,
                       float insertion_penalty, float deletion_penalty,
                       size_t num_threads) :
                            AnnotatedDBG(dbg, annotation, num_threads),
                            num_top_paths_(num_top_paths),
                            verbose_(verbose),
                            sw_threshold_(sw_threshold),
                            re_seeding_threshold_(re_seeding_threshold),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty) {
    match_score_ = 1;
    // Substitution score for each pair of nucleotides.
    // Transition and transversion mutations have different score values.
    sub_score_ = {
        {'a', {{'a', match_score_}, {'t', -2}, {'c', -2}, {'g', -1}}},
        {'g', {{'g', match_score_}, {'t', -2}, {'c', -2}, {'a', -1}}},
        {'c', {{'c', match_score_}, {'a', -2}, {'g', -2}, {'t', -1}}},
        {'t', {{'t', match_score_}, {'a', -2}, {'g', -2}, {'c', -1}}},
        {'$', {{'$', -2}}}};
}

namespace { // Helper function.
template <typename T>
float std_dev (const std::vector<T>& list, T mean) {
    float std_dev = 0;
    for (auto it = list.begin(); it != list.end(); ++it)
        std_dev += std::pow((*it) - mean, 2.0);
    return std::pow(std_dev / list.size(), 0.5);
}
}  // namespace

DBGAligner::AlignedPath DBGAligner::align(const std::string& sequence) const {
    if (sequence.size() < graph_->get_k())
        return AlignedPath(sequence.end());

    std::vector<AlignedPath> partial_paths;
    uint32_t print_frequency = 50;
    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath> queue(num_top_paths_);
    queue.push(std::move(AlignedPath(std::begin(sequence))));
    while (!queue.empty()) {
        auto path = std::move(queue.top());
        queue.pop();

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence)) {
            partial_paths.push_back(path);
            break;
        }

        if (path.get_total_score() < 0) {
            partial_paths.push_back(path);
            queue = BoundedPriorityQueue<AlignedPath>(num_top_paths_);
            assert(queue.empty());
            path = AlignedPath(path.get_sequence_it());
        }

        node_index seed = path.size() > 0 ? path.back() : graph_->npos;
        graph_->extend_from_seed(path.get_sequence_it(), std::end(sequence),
                                       [&](node_index node) {
                                        float score = (path.size() > 0 ? match_score_ :
                                                       match_score_ * graph_->get_k());
                                        if (node != graph_->npos)
                                                path.push_back(node, annotator_->get(node), score);
                                        },
                                       [&]() { return (path.size() > 0 &&
                                                (graph_->outdegree(path.back()) > 1)); },
                                       seed);

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence)) {
            partial_paths.push_back(path);
            break;
        }

        inexact_map(path, std::end(sequence),
                    [&](node_index node, std::string::const_iterator last_mapped_position) {
                        AlignedPath alternative_path(path);
                        alternative_path.set_sequence_it(last_mapped_position);
                        if (node == graph_->npos) {
                            // Seeding was not successful. To do inexact seeding,
                            // a path with incremented query_it is pushed to the queue.
                            queue.push(std::move(alternative_path));
                            return;
                        }
                        if (alternative_path.get_total_score() > sw_threshold_ * sequence.size()) {
                            // TODO: Construct sequence last char more efficiently.
                            alternative_path.push_back(node, annotator_->get(node),
                                single_char_score(*(alternative_path.get_sequence_it() +
                                    graph_->get_k() - 1), graph_->get_node_sequence(node).back()));
                        } else {
                            alternative_path.push_back(node, annotator_->get(node));
                            alternative_path.update_total_score(
                                whole_path_score(alternative_path, std::begin(sequence)));
                        }
                        DPAlignmentKey alternative_key{.node = alternative_path.back(),
                                           .query_it = alternative_path.get_sequence_it()};
                        DPAlignmentValue alternative_value = {.score = alternative_path.get_total_score()};
                        auto dp_alignment_it = dp_alignment.find(alternative_key);
                        if (dp_alignment_it == dp_alignment.end()) {
                            dp_alignment[alternative_key] = alternative_value;
                            queue.push(std::move(alternative_path));
                        }
                        else if (dp_alignment_it->second.score < alternative_path.get_total_score()) {
                            dp_alignment_it->second = alternative_value;
                            queue.push(std::move(alternative_path));
                        }
                        });

        if (verbose_ && path.size() % print_frequency == 0) {
            auto shadow_queue(queue);
            std::cout << "Loss mean and std of paths (top path size: "
                      << path.size() << "): ";
            std::vector<float> score_values;
            score_values.reserve(queue.size());
            while (!shadow_queue.empty()) {
                score_values.push_back(shadow_queue.top().get_total_score());
                shadow_queue.pop();
            }
            float mean = std::accumulate(score_values.begin(), score_values.end(), 0);
            mean /= score_values.size();
            std::cout << mean << ", " << std_dev<float>(score_values, mean) << std::endl;
        }
    }
    if (partial_paths.empty())
        return AlignedPath(std::end(sequence));

    // TODO: How to report the score for the case with a large gap???
    AlignedPath complete_path(partial_paths.front());
    for (auto partial_path = std::begin(partial_paths) + 1;
         partial_path != std::end(partial_paths); ++ partial_path) {
        // TODO: Find shortest path between the end node of a partial path
        // and the first node of the next partial path.
        complete_path.append_path(*partial_path);
    }
    complete_path.update_total_score(whole_path_score(complete_path, std::begin(sequence)));
    // TODO: Stitch partial paths together and return the corresponding path.
    // Insertions is the parts in query that aren't matched and for deletion,
    // we have to find approximation or exact distance using labeling schemes.
    return complete_path;
}

float DBGAligner::single_char_score(char char_in_query, char char_in_graph) const {
    try {
        return sub_score_.at(std::tolower(char_in_query)).at(std::tolower(char_in_graph));
    }
    catch (const std::out_of_range&) {
        return sub_score_.at('$').at('$');
    }
}

namespace {
// helper function to print an alignment based on SSW library.
static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
    std::cout << "===== SSW result =====" << std::endl;
    std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
       << "Reference start:\t" << alignment.ref_begin << std::endl
       << "Reference end:\t" << alignment.ref_end << std::endl
       << "Query start:\t" << alignment.query_begin << std::endl
       << "Query end:\t" << alignment.query_end << std::endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
       << "Number of mismatches:\t" << alignment.mismatches << std::endl
       << "Cigar: " << alignment.cigar_string << std::endl;
    std::cout << "======================" << std::endl;
}
}  // namespace

float DBGAligner::ssw_score(const AlignedPath& path, std::string::const_iterator begin) const {
    auto ref = get_path_sequence(path.get_nodes());
    std::string query(begin, path.get_sequence_it() + graph_->get_k() - 1);

    int32_t maskLen = strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    StripedSmithWaterman::Aligner aligner;
    if (!aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen)) {
        std::cout << "Failure in SSW calculation. Computing SSW inefficiently.\n";
        return whole_path_score(path, begin);
    }
    if (verbose_) {
        PrintAlignment(alignment);
    }
    return alignment.sw_score;
}

float DBGAligner::whole_path_score(const AlignedPath& path, std::string::const_iterator begin) const {
    auto path_sequence = get_path_sequence(path.get_nodes());
    std::string query_sequence(begin, path.get_sequence_it() + graph_->get_k() - 1);
    float alignment_score[query_sequence.size() + 1][path_sequence.size() + 1] = {};
    for (size_t i = 1; i <= query_sequence.size(); ++i) {
        for (size_t j = 1; j <= path_sequence.size(); ++j) {
            alignment_score[i][j] = std::max({alignment_score[i-1][j-1] +
                    single_char_score(query_sequence[i-1], path_sequence[j-1]),
                alignment_score[i][j-1] + insertion_penalty_, alignment_score[i-1][j] + deletion_penalty_});
        }
    }
    return alignment_score[query_sequence.size()][path_sequence.size()];
}

void DBGAligner::randomly_pick_strategy(std::vector<node_index> out_neighbors,
                                        const std::function<void(node_index)> &callback) const {
    // random_device creates a different seed everytime this function is called.
    std::random_device r;
    std::default_random_engine engine(r());
    std::uniform_int_distribution<int> uniform_dist(0, out_neighbors.size() - 1);
    callback(out_neighbors[uniform_dist(engine)]);
}

void DBGAligner::pick_all_strategy(std::vector<node_index> out_neighbors,
                                  const std::function<void(node_index)> &callback) const {
    for (auto &out_neighbor : out_neighbors) {
        callback(out_neighbor);
    }
}

void DBGAligner::inexact_map(const AlignedPath &path,
                             std::string::const_iterator end,
                             const std::function<void(node_index,
                                std::string::const_iterator)> &callback) const {
    if (path.get_sequence_it() + graph_->get_k() > end)
        return;
    if (path.size() == 0) {
        callback(graph_->npos, path.get_sequence_it() + 1);
        return;
    }
    std::vector<node_index> out_neighbors;
   // TODO: char is not used currently.
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char)
                                             { out_neighbors.push_back(node); });
    pick_all_strategy(out_neighbors, [&](node_index node) {
                           callback(node, path.get_sequence_it()); });
}

std::string DBGAligner::get_path_sequence(const std::vector<node_index>& path) const {
    if (path.size() < 1)
        return "";
    // TODO: Construct sequence more efficiently.
    std::string sequence = graph_->get_node_sequence(path.front());
    auto prev_node = path.front();
    for (auto path_it = std::next(std::begin(path)); path_it != std::end(path); ++path_it) {
        std::string next_seq = graph_->get_node_sequence(*path_it);
        if (graph_->traverse(prev_node, next_seq.back()) != graph_->npos)
            sequence.push_back(next_seq.back());
        else
            sequence += next_seq;
        prev_node = *path_it;
    }
    return sequence;
}
