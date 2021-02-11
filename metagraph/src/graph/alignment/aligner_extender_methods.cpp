#include "aligner_extender_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/logger.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;


template <typename NodeType>
DefaultColumnExtender<NodeType>::DefaultColumnExtender(const DeBruijnGraph &graph,
                                                       const DBGAlignerConfig &config,
                                                       std::string_view query)
      : graph_(graph), config_(config), query(query) {
    assert(config_.check_config_scores());
    partial_sums_.resize(query.size());
    std::transform(query.begin(), query.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query) == partial_sums_.front());
    assert(config_.get_row(query.back())[query.back()] == partial_sums_.back());

    for (char c : graph_.alphabet()) {
        auto &profile_score_row = profile_score.emplace(c, query.size() + 8).first.value();
        auto &profile_op_row = profile_op.emplace(c, query.size() + 8).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = Cigar::get_op_row(c);

        std::transform(query.begin(), query.end(), profile_score_row.begin(),
                       [&row](char q) { return row[q]; });
        std::transform(query.begin(), query.end(), profile_op_row.begin(),
                       [&op_row](char q) { return op_row[q]; });
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &seed) {
    // this extender only works if at least one character has been matched
    // assert(seed.get_query().size() > 1);
    // assert(query.data() <= seed.get_query().data());

    // const char *query_end = query.data() + query.size();
    // const char *seed_end = seed.get_query().data() + seed.get_query().size();
    // assert(query_end >= seed_end);

    // // size includes the last character of the seed since the upper-left corner
    // // of the score matrix is the seed score
    // size = query_end - seed_end + 2;

    // // extend_window_ doesn't include this last character
    // extend_window_ = { seed_end - 1, size + 1 };

    // match_score_begin = partial_sums_.data() + (seed_end - 1 - query.data());
    // assert(config_.get_row(query.back())[query.back()] == match_score_begin[size - 1]);
    // assert(config_.match_score(std::string_view(extend_window_.data() - 1, size))
    //     == *match_score_begin);

    // start_node = seed.back();
    path_ = &seed;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::operator()(ExtensionCallback callback,
                                                 score_t min_path_score) {
    std::ignore = min_path_score;
    typedef std::pair<NodeType, size_t> AlignNode;
    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<AlignNode> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::vector<std::tuple<ScoreVec, ScoreVec, ScoreVec, PrevVec, OpVec>> Column;
    tsl::hopscotch_map<NodeType, Column> table;

    const char *align_start = path_->get_query().data() + path_->get_query().size() - 1;
    size_t size = query.data() + query.size() - align_start + 1;
    extend_window_ = { align_start, size - 1 };
    assert(extend_window_[0] == path_->get_query().back());

    constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

    Cigar::Operator last_op = path_->get_cigar().back().first;
    assert(last_op == Cigar::MATCH || last_op == Cigar::MISMATCH);

    auto &[S, E, F, P, O] = table.emplace(
        DeBruijnGraph::npos,
        Column{ std::make_tuple(ScoreVec(size, ninf), ScoreVec(size, ninf),
                                ScoreVec(size, ninf), PrevVec(size, AlignNode{}),
                                OpVec(size, Cigar::CLIPPED)) }
    ).first.value()[0];

    S[0] = path_->get_score() - config_.get_row(extend_window_[0])[
        path_->get_sequence().back()
    ];

    score_t max_score = S[0];
    size_t max_pos = 0;
    AlignNode best_node{ DeBruijnGraph::npos, 0};
    std::vector<AlignNode> stack { best_node };
    while (stack.size()) {
        AlignNode prev = std::move(stack.back());
        stack.pop_back();

        if (prev.second > extend_window_.size() * 2)
            continue;

        std::vector<std::pair<NodeType, char>> outgoing;
        if (!prev.first) {
            outgoing.emplace_back(path_->back(), path_->get_sequence().back());
        } else {
            graph_.call_outgoing_kmers(prev.first, [&](auto&&... edge) {
                outgoing.emplace_back(std::move(edge)...);
            });
        }

        for (const auto &[next, c] : outgoing) {
            if (c == boss::BOSS::kSentinel)
                continue;

            Column &column = table[next];
            size_t depth = column.size();
            column.emplace_back(ScoreVec(size, 0), ScoreVec(size, ninf),
                                ScoreVec(size, ninf), PrevVec(size, AlignNode{}),
                                OpVec(size, Cigar::CLIPPED));

            AlignNode cur{ next, depth };

            auto &[S, E, F, P, O] = column.back();
            auto &[S_prev, E_prev, F_prev, P_prev, O_prev] = table[prev.first][prev.second];

            score_t del_score = F_prev[0] + config_.gap_extension_penalty;
            S[0] = std::max(del_score, 0);
            if (S[0] == del_score) {
                P[0] = prev;
                O[0] = Cigar::DELETION;
            }
            F[0] = std::max(del_score, S[0] + config_.gap_opening_penalty);

            for (size_t i = 1; i < S.size(); ++i) {
                score_t ins_val = E[i - 1] + config_.gap_extension_penalty;
                score_t match_val = S_prev[i - 1] + config_.get_row(c)[extend_window_[i - 1]];
                score_t del_val = F_prev[i] + config_.gap_extension_penalty;

                S[i] = std::max({ ins_val, del_val, match_val, 0 });

                if (S[i] == ins_val) {
                    P[i] = cur;
                    O[i] = Cigar::INSERTION;
                } else if (S[i] == del_val) {
                    P[i] = prev;
                    O[i] = Cigar::DELETION;
                } else if (S[i] == match_val) {
                    P[i] = prev;
                    O[i] = c == extend_window_[i - 1] ? Cigar::MATCH : Cigar::MISMATCH;
                }

                if (S[i] > max_score) {
                    max_score = S[i];
                    max_pos = i;
                    best_node = cur;
                }

                E[i] = std::max(ins_val, S[i] + config_.gap_opening_penalty);
                F[i] = std::max(del_val, S[i] + config_.gap_opening_penalty);
            }

            auto max_it = std::max_element(S.begin(), S.end());
            if (*max_it >= max_score - config_.xdrop)
                stack.emplace_back(cur);
        };
    }

    if (max_pos < 2 && best_node.first == path_->back() && !best_node.second)
        return;

    // traceback
    assert(std::get<0>(table[best_node.first][best_node.second])[max_pos] == max_score);

    Cigar cigar;
    if (max_pos + 1 < F.size())
        cigar.append(Cigar::CLIPPED, F.size() - max_pos - 1);

    size_t pos = max_pos;
    std::vector<NodeType> path;

#ifndef NDEBUG
    score_t score_track = max_score;
#endif

    while (best_node.first || best_node.second) {
        auto &[S, E, F, P, O] = table[best_node.first][best_node.second];

#ifndef NDEBUG
        if (cigar.size() && cigar.back().first != O[pos] && (cigar.back().first == Cigar::INSERTION || cigar.back().first == Cigar::DELETION)) {
            std::cout << "foo\n";
            score_track -= config_.gap_opening_penalty - config_.gap_extension_penalty;
        }

        std::cout << score_track << "\t"
                  << pos << "\t"
                  << Cigar::opt_to_char(O[pos]) << "\t"
                  << best_node.first << ","
                  << best_node.second << "\t"
                  << P[pos].first << ","
                  << P[pos].second << "\t"
                  << S[pos] << ","
                  << E[pos] << ","
                  << F[pos] << "\t"
                  << path.size() << "\t"
                  << cigar.to_string() << "\n";

        assert(score_track == S[pos]);
#endif

        switch (O[pos]) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                cigar.append(O[pos]);
                path.push_back(best_node.first);
#ifndef NDEBUG
                char last_char = graph_.get_node_sequence(best_node.first).back();
                assert((last_char == extend_window_[pos - 1]) == (O[pos] == Cigar::MATCH));
                score_track -= config_.get_row(extend_window_[pos - 1])[last_char];
#endif
                best_node = P[pos];
                --pos;
            } break;
            case Cigar::INSERTION: {
                assert(P[pos] == best_node);
                cigar.append(O[pos]);
#ifndef NDEBUG
                score_track -= config_.gap_extension_penalty;
#endif
                --pos;
            } break;
            case Cigar::DELETION: {
                path.push_back(best_node.first);
                cigar.append(O[pos]);
                best_node = P[pos];
#ifndef NDEBUG
                score_track -= config_.gap_extension_penalty;
#endif
            } break;
            case Cigar::CLIPPED: {
                best_node = AlignNode{};
                assert(!pos || S[pos] == 0);
            }
        }
    }

    assert(path.size());

    NodeType start_node;
    score_t score;
    if (pos > 1) {
        cigar.append(Cigar::CLIPPED, pos - 1);
        start_node = path.front();
        score = max_score;
    } else {
        score = max_score - path_->get_score();
        start_node = path_->back();
        while (pos < 1) {
            ++pos;
            path.pop_back();
            --cigar.back().second;
            if (!cigar.back().second)
                cigar.pop_back();
        }
    }

    std::reverse(cigar.begin(), cigar.end());

    std::reverse(path.begin(), path.end());
    std::string seq;
    spell_path(graph_, path, seq, graph_.get_k() - 1);

    Alignment<NodeType> extension({ extend_window_.data() + pos, max_pos - pos },
                                  std::move(path),
                                  std::move(seq),
                                  score,
                                  std::move(cigar),
                                  0,
                                  path_->get_orientation(),
                                  graph_.get_k() - 1);

    assert(extension.is_valid(graph_, &config_));
    callback(std::move(extension), start_node);
}


template <typename NodeType>
auto DefaultColumnExtender<NodeType>::emplace_node(NodeType, NodeType, char, size_t, size_t, size_t, size_t, size_t)
-> std::pair<typename DPTable<NodeType>::iterator, bool> { return {}; }

template <typename NodeType>
bool DefaultColumnExtender<NodeType>::add_seed(size_t) { return true; }

template <typename NodeType>
void DefaultColumnExtender<NodeType>::check_and_push(ColumnRef&&) {}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::extend_main(ExtensionCallback, score_t) {}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::update_columns(NodeType,
                    const std::deque<std::pair<NodeType, char>> &,
                    score_t) {}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>::fork_extension(NodeType, ExtensionCallback, score_t)
-> std::deque<std::pair<NodeType, char>> { return {}; }


template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
