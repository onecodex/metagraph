#include "aligner_labeled.hpp"

#include "graph/annotated_dbg.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/vector_map.hpp"
#include "common/utils/template_utils.hpp"

namespace mtg {
namespace graph {
namespace align {

inline void reverse_bit_vector(sdsl::bit_vector &v) {
    size_t begin = 0;
    for ( ; begin + begin + 128 <= v.size(); begin += 64) {
        uint64_t a = sdsl::bits::rev(v.get_int(begin));
        uint64_t b = sdsl::bits::rev(v.get_int(v.size() - begin - 64));
        v.set_int(begin, b);
        v.set_int(v.size() - begin - 64, a);
    }

    size_t size = (v.size() % 128) / 2;
    uint64_t a = sdsl::bits::rev(v.get_int(begin, size)) >> (64 - size);
    uint64_t b = sdsl::bits::rev(v.get_int(v.size() - begin - size, size)) >> (64 - size);
    v.set_int(begin, b, size);
    v.set_int(v.size() - begin - size, a, size);
}

void
process_seq_path(const DeBruijnGraph &graph,
                 std::string_view query,
                 const std::vector<DeBruijnGraph::node_index> &query_nodes,
                 const std::function<void(AnnotatedDBG::row_index, size_t)> &callback) {
    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    if (canonical) {
        if (query_nodes.size()) {
            auto first = std::find_if(query_nodes.begin(), query_nodes.end(),
                                      [](auto i) -> bool { return i; });
            if (first == query_nodes.end())
                return;

            size_t start = first - query_nodes.begin();

            if (canonical->get_base_node(*first) == *first) {
                for (size_t i = start; i < query_nodes.size(); ++i) {
                    if (query_nodes[i] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i])
                        ), i);
                    }
                }
            } else {
                for (size_t i = query_nodes.size(); i > start; --i) {
                    if (query_nodes[i - 1] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i - 1])
                        ), i - 1);
                    }
                }
            }
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos)
                callback(AnnotatedDBG::graph_to_anno_index(query_nodes[i]), i);
        }
    } else {
        size_t i = 0;
        if (query.front() == '#') {
            std::string map_query = graph.get_node_sequence(query_nodes[0]).substr(0, graph.get_k() - 1);
            map_query += query.substr(graph.get_k() - 1);
            graph.map_to_nodes(map_query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        } else {
            graph.map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        }
        assert(i == query_nodes.size());
    }
}

ILabeledDBGAligner::ILabeledDBGAligner(const AnnotatedDBG &anno_graph,
                                       const DBGAlignerConfig &config,
                                       size_t num_top_labels)
      : anno_graph_(anno_graph),
        graph_(anno_graph_.get_graph()),
        config_(config), num_top_labels_(num_top_labels) {}

auto ILabeledDBGAligner
::map_and_label_query_batch(const QueryGenerator &generate_query) const
        -> std::pair<BatchMapping, BatchLabels> {
    // exact k-mer matchings of each query sequence
    BatchMapping query_nodes;

    // map from Annotation row indices to (i,j), indicating position j in query i
    typedef std::vector<std::pair<size_t, size_t>> RowMapping;
    VectorMap<AnnotatedDBG::row_index, std::pair<RowMapping, RowMapping>> row_to_query_idx;

    // populate maps
    generate_query([&](std::string_view, std::string_view query, bool) {
        size_t i = query_nodes.size();

        // map query sequence to the graph
        query_nodes.emplace_back(map_sequence_to_nodes(graph_, query),
                                 std::vector<node_index>{});

        // update row_to_query_idx
        process_seq_path(graph_, query, query_nodes.back().first,
                         [&](AnnotatedDBG::row_index row, size_t j) {
            row_to_query_idx[row].first.emplace_back(i, j);
        });

        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            query_nodes.back().second = query_nodes.back().first;
            std::string query_rc(query);
            reverse_complement_seq_path(graph_, query_rc, query_nodes.back().second);

            if (graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                process_seq_path(graph_, query_rc, query_nodes.back().second,
                                 [&](AnnotatedDBG::row_index row, size_t j) {
                    row_to_query_idx[row].second.emplace_back(i, j);
                });
            }
        }
    });

    // extract rows from the row index map
    std::vector<AnnotatedDBG::row_index> rows;
    rows.reserve(row_to_query_idx.size());
    for (const auto &[row, mapping] : row_to_query_idx) {
        rows.push_back(row);
    }

    // get annotations for each row
    auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

    // we now have two maps
    // 1) row indices -> column IDs
    // 2) row indices -> (query i, query_i j) i.e., row -> the k-mer at batch[i][j]
    // and wish to construct the following map
    // query -> ( (col_1, col_1_signature), (col_2, col_2_signature), ... )

    // count labels for each query
    std::vector<VectorMap<uint64_t, uint64_t>> column_counter(query_nodes.size());
    for (size_t i = 0; i < annotation.size(); ++i) {
        for (const auto &[query_id, idx] : row_to_query_idx[rows[i]].first) {
            for (uint64_t column : annotation[i]) {
                ++column_counter[query_id][column];
            }
        }
    }

    std::vector<VectorMap<uint64_t, uint64_t>> column_counter_rc;
    if (config_.forward_and_reverse_complement
            && graph_.get_mode() != DeBruijnGraph::CANONICAL) {
        column_counter_rc.resize(query_nodes.size());
        for (size_t i = 0; i < annotation.size(); ++i) {
            for (const auto &[query_id, idx] : row_to_query_idx[rows[i]].second) {
                for (uint64_t column : annotation[i]) {
                    ++column_counter_rc[query_id][column];
                }
            }
        }
    }

    // compute target columns and initialize signatures for each query
    BatchLabels target_columns(query_nodes.size());
    for (size_t i = 0; i < column_counter.size(); ++i) {
        auto &counter = const_cast<std::vector<std::pair<uint64_t, uint64_t>>&>(
            column_counter[i].values_container()
        );

        // pick the top columns for each query sequence
        size_t num_targets = std::min(counter.size(), num_top_labels_);

        std::sort(counter.begin(), counter.end(), utils::GreaterSecond());
        if (num_targets < counter.size())
            counter.resize(num_targets);

        for (const auto &[column, count] : counter) {
            if (!target_columns[i].count(column)) {
                target_columns[i].emplace(column, Signature {
                    sdsl::bit_vector(query_nodes[i].first.size(), false),
                    config_.forward_and_reverse_complement
                            && graph_.get_mode() != DeBruijnGraph::CANONICAL
                        ? sdsl::bit_vector(query_nodes[i].second.size(), false)
                        : sdsl::bit_vector()
                });
            }
        }
    }

    for (size_t i = 0; i < column_counter_rc.size(); ++i) {
        assert(config_.forward_and_reverse_complement
            && graph_.get_mode() != DeBruijnGraph::CANONICAL);
        auto &counter = const_cast<std::vector<std::pair<uint64_t, uint64_t>>&>(
            column_counter_rc[i].values_container()
        );

        // pick the top columns for each query sequence
        size_t num_targets = std::min(counter.size(), num_top_labels_);

        std::sort(counter.begin(), counter.end(), utils::GreaterSecond());
        if (num_targets < counter.size())
            counter.resize(num_targets);

        for (const auto &[column, count] : counter) {
            if (!target_columns[i].count(column)) {
                target_columns[i].emplace(column, Signature {
                    sdsl::bit_vector(query_nodes[i].first.size(), false),
                    sdsl::bit_vector(query_nodes[i].second.size(), false)
                });
            }
        }
    }

    // fill signatures for each query
    for (size_t i = 0; i < annotation.size(); ++i) {
        for (const auto &[query_id, idx] : row_to_query_idx[rows[i]].first) {
            for (uint64_t column : annotation[i]) {
                if (target_columns[query_id].count(column))
                    target_columns[query_id][column].first[idx] = true;
            }
        }

        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            if (graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                for (const auto &[query_id, idx] : row_to_query_idx[rows[i]].second) {
                    for (uint64_t column : annotation[i]) {
                        if (target_columns[query_id].count(column))
                            target_columns[query_id][column].second[idx] = true;
                    }
                }
            } else {
                for (size_t j = 0; j < target_columns.size(); ++j) {
                    for (auto it = target_columns[j].begin();
                            it != target_columns[j].end(); ++it) {
                        it.value().second = it.value().first;
                        reverse_bit_vector(it.value().second);
                    }
                }
            }
        }
    }

    // if any of the queries has no associated labels, mark them as such
    for (size_t i = 0; i < target_columns.size(); ++i) {
        if (target_columns[i].empty()) {
            assert(!query_nodes[i].first[0]);
            assert(std::equal(query_nodes[i].first.begin() + 1,
                              query_nodes[i].first.end(),
                              query_nodes[i].first.begin()));
            target_columns[i].emplace(kNTarget, Signature {
                sdsl::bit_vector(query_nodes[i].first.size(), true),
                graph_.get_mode() == DeBruijnGraph::CANONICAL
                        || config_.forward_and_reverse_complement
                    ? sdsl::bit_vector(query_nodes[i].second.size(), true)
                    : sdsl::bit_vector()
            });
        }
    }

    return { std::move(query_nodes), std::move(target_columns) };
}

template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph),
        target_columns_({ ILabeledDBGAligner::kNTarget }) {}

template <typename NodeType>
void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    DefaultColumnExtender<NodeType>::initialize(path);
    align_node_to_target_.clear();

    assert(path.target_column != std::numeric_limits<uint64_t>::max());
    size_t target_column_idx = std::find(target_columns_.begin(), target_columns_.end(),
                                         path.target_column) - target_columns_.begin();

    align_node_to_target_[{ this->graph_.max_index() + 1, '\0', 0, 0 }] = target_column_idx;

    if (target_column_idx == target_columns_.size())
        target_columns_.push_back(path.target_column);
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>::get_outgoing(const AlignNode &node) const -> Edges {
    assert(align_node_to_target_.count(node));

    if (std::get<0>(node) == this->graph_.max_index() + 1)
        return DefaultColumnExtender<NodeType>::get_outgoing(node);

    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_);

    uint64_t target_column_idx = align_node_to_target_[node];
    uint64_t target_column = target_columns_.at(target_column_idx);

    if (target_column == ILabeledDBGAligner::kNTarget) {
        assert(this->seed_->get_offset());
        assert(this->seed_->get_offset() + 1 >= std::get<3>(node));
        size_t next_offset = this->seed_->get_offset() + 1 - std::get<3>(node);

        Edges edges = DefaultColumnExtender<NodeType>::get_outgoing(node);
        if (next_offset)
            return edges;

        std::vector<AnnotatedDBG::row_index> rows;
        rows.reserve(edges.size());

        for (const auto &[next_node, c] : edges) {
            if (canonical) {
                rows.emplace_back(AnnotatedDBG::graph_to_anno_index(
                    canonical->get_base_node(next_node)
                ));
            } else if (this->graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                rows.emplace_back(AnnotatedDBG::graph_to_anno_index(next_node));
            } else {
                this->graph_.map_to_nodes(this->graph_.get_node_sequence(next_node),
                                          [&](NodeType next_canonical) {
                    rows.emplace_back(AnnotatedDBG::graph_to_anno_index(next_canonical));
                });
            }
        }

        assert(rows.size() == edges.size());

        auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

        Edges out_edges;
        for (size_t i = 0; i < rows.size(); ++i) {
            AlignNode next { edges[i].first, edges[i].second, 0, std::get<3>(node) + 1 };
            auto find = this->table_.find(edges[i].first);
            if (find != this->table_.end())
                std::get<2>(next) = find->second.first.size();

            for (uint64_t target : annotation[i]) {
                target_column_idx = std::find(target_columns_.begin(),
                                              target_columns_.end(),
                                              target) - target_columns_.begin();

                if (target_column_idx == target_columns_.size())
                    target_columns_.emplace_back(target);

                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = target_column_idx;

                out_edges.push_back(edges[i]);
                ++std::get<2>(next);
            }
        }

        return out_edges;
    }

    if (cached_edge_sets_.count(std::get<0>(node))) {
        const auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
        if (target_column_idx < edge_sets.size())
            return edge_sets[target_column_idx];
    }

    typedef std::tuple<NodeType /* parent */,
                       NodeType /* child */,
                       char /* edge label */,
                       size_t /* index in query seq */> EdgeDescriptor;
    VectorMap<AnnotatedDBG::row_index, std::vector<EdgeDescriptor>> anno_rows_to_id;

    size_t k = this->graph_.get_k();
    auto push_path = [&](std::string_view seq, const std::vector<NodeType> &path) {
        process_seq_path(this->graph_, seq, path,
                         [&](AnnotatedDBG::row_index row, size_t i) {
            if (seq[i + k - 1] != boss::BOSS::kSentinel) {
                anno_rows_to_id[row].emplace_back(
                    i ? path[i - 1] : std::get<0>(node),
                    path[i], seq[i + k - 1], i
                );
            }
        });
    };

    // prefetch the next unitig
    const auto &column = this->table_.find(std::get<0>(node))->second;
    auto [min_i, max_i] = this->get_band(node, column, this->xdrop_cutoff_);

    size_t bandwidth = max_i - min_i;
    const score_t *S = &std::get<0>(column.first[std::get<2>(node)])[min_i];
    std::string_view q_min = this->query_.substr(this->start_);

    assert(this->query_.size() >= min_i + this->start_);
    size_t max_depth = this->query_.size() - min_i - this->start_;

    call_hull_sequences(this->graph_, std::get<0>(node),
        [&](std::string_view seq, const std::vector<NodeType> &path) {
            assert(path.size());
            assert(this->graph_.traverse(
                std::get<0>(node), seq[this->graph_.get_k() - 1]) == path[0]);
            push_path(seq, path);
        },
        [&](std::string_view seq, const auto &path, size_t depth, size_t fork_count) {
            bool result = fork_count || depth > max_depth
                || (path.size() && cached_edge_sets_.count(path.back()))
                || (path.size() >= 2 && canonical
                    && (canonical->get_base_node(path.back()) == path.back())
                        != (canonical->get_base_node(path[path.size() - 2])
                            == path[path.size() - 2]));

            if (!result) {
                // compute xdrop cutoffs
                bool has_extension = false;
                std::string_view ref = seq.substr(this->graph_.get_k() - 1);
                std::string_view qu = q_min.substr(depth + 1 - ref.size());
                for (size_t i = 0; i < bandwidth; ++i) {
                    score_t ext_score = S[i] + this->config_.score_sequences(
                        ref, { qu.data() + i, ref.size() }
                    );

                    if (ext_score >= this->xdrop_cutoff_) {
                        has_extension = true;
                        break;
                    }
                }

                result = !has_extension;
            }

            if (result && depth == 1)
                push_path(seq, path);

            return result;
        }
    );

    std::vector<AnnotatedDBG::row_index> anno_rows;
    anno_rows.reserve(anno_rows_to_id.size());
    for (const auto &pair : anno_rows_to_id) {
        anno_rows.push_back(pair.first);
    }

    sdsl::bit_vector row_mask = anno_graph_.get_annotation().get_matrix().has_column(
        anno_rows, target_column
    );

    Edges edges;
    for (size_t j = 0; j < anno_rows.size(); ++j) {
        if (row_mask[j]) {
            AnnotatedDBG::row_index row = anno_rows[j];
            for (const auto &[parent_node, child_node, c, i] : anno_rows_to_id[row]) {
                assert(c != boss::BOSS::kSentinel);
                assert(this->graph_.traverse(parent_node, c) == child_node);

                if (!i) {
                    assert(parent_node == std::get<0>(node));
                    edges.emplace_back(child_node, c);
                } else {
                    auto &parent_edge_sets = cached_edge_sets_[parent_node];
                    if (target_column_idx >= parent_edge_sets.size())
                        parent_edge_sets.resize(target_column_idx + 1);

                    cached_edge_sets_[parent_node][target_column_idx].emplace_back(
                        child_node, c
                    );
                }
            }
        }
    }

    auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
    if (target_column_idx >= edge_sets.size()) {
        edge_sets.resize(target_column_idx + 1);
        edge_sets[target_column_idx] = edges;
#ifndef NDEBUG
    } else {
        assert(edge_sets[target_column_idx] == edges);
#endif
    }

    return edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg