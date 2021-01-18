#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/vector_map.hpp"
#include "common/hashers/hash.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class Seeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class LabeledDBGAligner : public SeedAndExtendAligner<Seeder, Extender> {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    LabeledDBGAligner(const AnnotatedDBG &anno_graph,
                      const DBGAlignerConfig &config,
                      size_t num_top_labels = 1)
          : anno_graph_(anno_graph), config_(config), num_top_labels_(num_top_labels) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual ~LabeledDBGAligner() {}

    const DeBruijnGraph& get_graph() const override { return anno_graph_.get_graph(); }
    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&)> SeedGenerator;
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    Seeder build_seeder() const override { return Seeder(get_graph(), config_); }
    Extender build_extender() const override {
        return Extender(anno_graph_, config_, target_column_);
    }

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths. Each path has the property that there exists
    // at least one label which is shared by all nodes
    void align_aggregate(DBGQueryAlignment &paths,
                         const AlignmentGenerator &alignment_generator) const override;

    SeedGenerator build_seed_generator(const std::string_view query,
                                       bool orientation) const override;

  private:
    const AnnotatedDBG& anno_graph_;
    const DBGAlignerConfig config_;
    size_t num_top_labels_;
    mutable uint64_t target_column_;
};

template <typename NodeType>
class LabeledColumnExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef typename Extender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename Extender<NodeType>::score_t score_t;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config),
            anno_graph_(anno_graph) {}

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          uint64_t target_column)
          : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config),
            anno_graph_(anno_graph),
            target_columns_({ target_column }) {}

    virtual ~LabeledColumnExtender() {}

  protected:
    typedef std::deque<std::pair<NodeType, char>> Edges;

    virtual Edges fork_extension(NodeType node,
                                 std::function<void(DBGAlignment&&, NodeType)> callback,
                                 score_t min_path_score) override;

  private:
    const AnnotatedDBG &anno_graph_;
    std::vector<uint64_t> target_columns_;
    tsl::hopscotch_map<NodeType, Edges> cached_edge_sets_;
};


template <class Seeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<Seeder, Extender, AlignmentCompare>
::align_aggregate(DBGQueryAlignment &paths,
                  const AlignmentGenerator &alignment_generator) const {
    AlignmentAggregator<node_index, AlignmentCompare> path_queue(
        paths.get_query(), paths.get_query_reverse_complement(), config_
    );

    alignment_generator(
        [&](DBGAlignment&& alignment) { path_queue.add_alignment(std::move(alignment)); },
        [&](const DBGAlignment &seed) { return path_queue.get_min_path_score(seed); }
    );

    path_queue.call_alignments([&](auto&& alignment) {
        assert(alignment.is_valid(get_graph(), &config_));
        paths.emplace_back(std::move(alignment));
    });
}

template <class Seeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<Seeder, Extender, AlignmentCompare>
::build_seed_generator(const std::string_view query,
                       bool orientation) const -> SeedGenerator {
    static_assert(std::is_base_of_v<ExactMapSeeder<node_index>, Seeder>
        || std::is_base_of_v<SuffixSeeder<node_index>, Seeder>);

    return [this,query,orientation](const auto &callback) {
        auto seeder = build_seeder();
        seeder.initialize(query, orientation);
        auto &query_nodes = const_cast<std::vector<node_index>&>(seeder.get_query_nodes());

        if (std::all_of(query_nodes.begin(), query_nodes.end(), [](auto i) { return !i; }))
            return;

        auto &offsets = const_cast<std::vector<uint8_t>&>(seeder.get_offsets());
        std::vector<uint8_t> *flags = nullptr;
        const std::vector<std::vector<node_index>> *alt_query_nodes = nullptr;

        if constexpr(std::is_base_of_v<MEMSeeder<node_index>, Seeder>)
            flags = &const_cast<std::vector<uint8_t>&>(seeder.get_query_node_flags());

        if constexpr(std::is_base_of_v<SuffixSeeder<node_index>, Seeder>) {
            alt_query_nodes = &seeder.get_alt_query_nodes();

            const auto &base_seeder = seeder.get_base_seeder();
            if (const auto *mem = dynamic_cast<const MEMSeeder<node_index>*>(&base_seeder))
                flags = &const_cast<std::vector<uint8_t>&>(mem->get_query_node_flags());
        }

        auto signatures = anno_graph_.get_top_label_signatures(
            query, num_top_labels_, config_.exact_kmer_match_fraction,
            alt_query_nodes ? &query_nodes : nullptr, alt_query_nodes
        );

        if (signatures.empty())
            return;

        const auto &[label, signature] = signatures[0];

        for (size_t i = 0; i < signature.size(); ++i) {
            if (!signature[i]) {
                query_nodes[i] = 0;
                offsets[i] = 0;
                if (flags)
                    (*flags)[i] = 0;
            }
        }

        target_column_ = anno_graph_.get_annotation().get_label_encoder().encode(label);

        seeder.call_seeds(callback);
    };
}

template <typename NodeType>
inline auto LabeledColumnExtender<NodeType>
::fork_extension(NodeType node,
                 std::function<void(DBGAlignment&&, NodeType)> callback,
                 score_t min_path_score) -> Edges {
    if (cached_edge_sets_.count(node))
        return cached_edge_sets_[node];

    const auto &mat = anno_graph_.get_annotation().get_matrix();

    // get set of outgoing nodes from the parent class
    auto base_edges = DefaultColumnExtender<NodeType>::fork_extension(
        node, callback, min_path_score
    );

    Edges edges;

    // first check the simple case to avoid decoding entire rows
    if (target_columns_.size() == 1) {
        for (auto&& edge : base_edges) {
            if (mat.get(anno_graph_.graph_to_anno_index(edge.first),
                        target_columns_[0])) {
                edges.emplace_back(std::move(edge));
            }
        }

        cached_edge_sets_[node] = edges;

        return edges;
    }

    // decode all rows
    std::vector<AnnotatedDBG::row_index> base_rows;
    base_rows.reserve(base_edges.size());
    for (const auto &[out_node, c] : base_edges) {
        base_rows.push_back(anno_graph_.graph_to_anno_index(out_node));
    }

    AnnotatedDBG::row_index node_row = anno_graph_.graph_to_anno_index(node);

    std::vector<Vector<uint64_t>> rows;
    if (const auto *rd = dynamic_cast<const annot::binmat::IRowDiff*>(&mat)) {
        rows.resize(base_rows.size());
        for (size_t i = 0; i < base_rows.size(); ++i) {
            if (rd->is_anchor(base_rows[i])) {
                rows[i] = rd->get_diff(base_rows[i]);
                std::sort(rows[i].begin(), rows[i].end());
            } else if (rd->is_anchor(node_row)) {
                rows[i] = mat.get_row(base_rows[i]);
            } else {
                const auto &dbg_succ = *rd->graph();
                const auto &boss = dbg_succ.get_boss();
                if (boss.get_last(dbg_succ.kmer_to_boss_index(base_edges[i].first))) {
                    // apply diff to target_columns_
                    auto diff_row = rd->get_diff(node_row);
                    std::sort(diff_row.begin(), diff_row.end());
                    std::set_difference(target_columns_.begin(), target_columns_.end(),
                                        diff_row.begin(), diff_row.end(),
                                        std::back_inserter(rows[i]));
                } else {
                    rows[i] = mat.get_row(base_rows[i]);
                }
            }
        }
    } else {
        rows = mat.get_rows(base_rows);
        for (auto &row : rows) {
            std::sort(row.begin(), row.end());
        }
    }

    // aggregate outgoing nodes by row
    tsl::hopscotch_map<std::vector<uint64_t>, Edges, utils::VectorHash> out_labels;

    for (size_t i = 0; i < rows.size(); ++i) {
        assert(std::is_sorted(rows[i].begin(), rows[i].end()));
        out_labels[{ rows[i].begin(), rows[i].end() }].emplace_back(
            std::move(base_edges[i])
        );
    }

    // copy the current extender with the new row and continue extension
    auto fork_extender = [&](std::vector<uint64_t>&& new_target_labels,
                             std::deque<std::pair<NodeType, char>>&& cur_edges) {
        if (!new_target_labels.empty()) {
            auto fork = *this;
            fork.cached_edge_sets_.clear();
            fork.cached_edge_sets_[node] = cur_edges;
            fork.target_columns_ = std::move(new_target_labels);
            fork.update_columns(node, std::move(cur_edges), min_path_score);
            fork.extend_main([&](DBGAlignment&& extension, NodeType start_node) {
                if (start_node)
                    callback(std::move(extension), start_node);
            }, min_path_score);
        }
    };

    for (auto it = out_labels.begin(); it != out_labels.end(); ++it) {
        const auto &row = it->first;
        auto &cur_edges = it.value();

        // if the current target row is empty, fork when labels are found
        if (target_columns_.empty()) {
            if (row.empty()) {
                swap(edges, cur_edges);
            } else {
                fork_extender(std::vector<uint64_t>(row), std::move(cur_edges));
            }
            continue;
        }

        std::vector<AnnotatedDBG::row_index> intersection;
        intersection.reserve(std::min(row.size(), target_columns_.size()));
        std::set_intersection(target_columns_.begin(), target_columns_.end(),
                              row.begin(), row.end(),
                              std::back_inserter(intersection));

        if (intersection.size() == target_columns_.size()) {
            // the row matches, so pass these outgoing nodes to the current extender
            swap(edges, cur_edges);
        } else {
            // create a new extender to work on the subset
            fork_extender(std::move(intersection), std::move(cur_edges));
        }
    }

    cached_edge_sets_[node] = edges;

    return edges;
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
