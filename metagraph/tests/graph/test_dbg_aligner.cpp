#include <gtest/gtest.h>

#include "test_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "dbg_aligner.hpp"
#include "aligner_methods.hpp"

#include "annotate_column_compressed.hpp"
#include "reverse_complement.hpp"
#include "alphabets.hpp"


typedef DBGAligner<>::score_t score_t;

int8_t single_char_score(const DBGAlignerConfig &config, char a, int8_t b) {
    return config.get_row(a)[b];
}

void check_score_matrix(const DBGAlignerConfig &config,
                        const char* alphabet,
                        size_t alph_size) {
    ASSERT_EQ(strlen(alphabet), alph_size);

    for (size_t i = 0; i < alph_size; ++i) {
        if (i + 1 != alph_size) {
            EXPECT_LT(int8_t(0), single_char_score(config, alphabet[i], alphabet[i]));
        }

        for (size_t j = 0; j < alph_size; ++j) {
            //check if the match score is the greatest
            if (i + 1 != alph_size) {
                EXPECT_GE(single_char_score(config, alphabet[i], alphabet[i]),
                          single_char_score(config, alphabet[i], alphabet[j]));
            }

            // checking symmetry
            EXPECT_EQ(single_char_score(config, alphabet[i], alphabet[j]),
                      single_char_score(config, alphabet[j], alphabet[i]));
        }
    }
}

TEST(DBGAlignerTest, check_score_matrix_dna) {
    check_score_matrix(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2),
                       alphabets::kAlphabetDNA5,
                       alphabets::kSigmaDNA5);
}

TEST(DBGAlignerTest, check_score_matrix_protein) {
    check_score_matrix(DBGAlignerConfig::score_matrix_blosum62,
                       alphabets::kAlphabetProtein,
                       alphabets::kSigmaProtein);
}

TEST(DBGAlignerTest, check_score_matrix_dna_unit) {
    check_score_matrix(
        DBGAlignerConfig::unit_scoring_matrix(
            1, alphabets::kAlphabetDNA, alphabets::kCharToDNA
        ),
        alphabets::kAlphabetDNA5,
        alphabets::kSigmaDNA5
    );

    check_score_matrix(
        DBGAlignerConfig::unit_scoring_matrix(
            1, alphabets::kAlphabetDNA5, alphabets::kCharToDNA5
        ),
        alphabets::kAlphabetDNA5,
        alphabets::kSigmaDNA5
    );
}

TEST(DBGAlignerTest, check_score_matrix_protein_unit) {
    check_score_matrix(
        DBGAlignerConfig::unit_scoring_matrix(
            1, alphabets::kAlphabetProtein, alphabets::kCharToProtein
        ),
        alphabets::kAlphabetProtein,
        alphabets::kSigmaProtein
    );
}

const DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));


// TODO: REPLACE THIS
#if _PROTEIN_GRAPH
    const auto *alphabet = alphabets::kAlphabetProtein;
    const auto *alphabet_encoding = alphabets::kCharToProtein;
#elif _DNA_CASE_SENSITIVE_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNACaseSent;
    const auto *alphabet_encoding = alphabets::kCharToDNACaseSent;
#elif _DNA5_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNA5;
    const auto *alphabet_encoding = alphabets::kCharToDNA5;
#elif _DNA_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNA;
    const auto *alphabet_encoding = alphabets::kCharToDNA;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

bool check_extend(std::shared_ptr<const DeBruijnGraph> graph,
                  const DBGAlignerConfig &config,
                  const DBGAligner<>::DBGQueryAlignment &paths,
                  const std::string &query) {
    assert(graph.get());
    EXPECT_EQ(query, paths.get_query());

    return paths == DBGAligner<UniMEMSeeder<>>(*graph, config).align(query);
}


template <typename Graph>
class DBGAlignerTest : public DeBruijnGraphTest<Graph> {
    void SetUp() { Cigar::initialize_opt_table(alphabet, alphabet_encoding); }
};

TYPED_TEST_CASE(DBGAlignerTest, GraphTypes);

TYPED_TEST(DBGAlignerTest, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);

    EXPECT_EQ(0ull, alt_paths.size());
}

TYPED_TEST(DBGAlignerTest, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("3=", path.get_cigar().to_string());
    EXPECT_EQ(3u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TYPED_TEST(DBGAlignerTest, align_straight_forward_and_reverse_complement) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;
    reverse_complement(query.begin(), query.end());

    auto graph = build_graph_batch<TypeParam>(k, { reference });

    auto config_fwd_and_rev = config;
    config_fwd_and_rev.forward_and_reverse_complement = true;

    DBGAligner<> aligner(*graph, config_fwd_and_rev);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    auto ext_paths = DBGAligner<UniMEMSeeder<>>(*graph, config_fwd_and_rev).align(query);

    EXPECT_TRUE(std::equal(paths.begin(), paths.end(),
                           ext_paths.begin(), ext_paths.end()));

    // test copy
    auto paths_copy = const_cast<const DBGAligner<>::DBGQueryAlignment&>(paths);
    for (const auto &path : paths_copy) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }

    // test move
    auto paths_move = std::move(paths);
    for (const auto &path : paths_move) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }
}


TYPED_TEST(DBGAlignerTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("9=", path.get_cigar().to_string());
    EXPECT_EQ(9u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, align_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGA" "AT" "ATTTGTT";
    std::string reference_2 = "AGCTTCGA" "CG" "ATTTGTT";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("17=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =     "AGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()),
              path.get_score());
    EXPECT_EQ("6=", path.get_cigar().to_string());
    EXPECT_EQ(6u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, variation) {
    size_t k = 4;
    std::string reference = "AGCAA" "C" "TCGAAA";
    std::string query =     "AGCAA" "T" "TCGAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                     reference.begin()),
              path.get_score());
    EXPECT_EQ("5=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(11u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAA" "CTC" "GAAA";
    std::string reference_2 = "TTAAGCAA" "GTC" "GAAA";
    std::string query =       "TTAAGCAA" "TGG" "GAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    DBGAligner<> aligner(*graph, config);

    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_TRUE(path.get_sequence().compare(reference_1) == 0 ||
                path.get_sequence().compare(reference_2) == 0)
        << "Path: " << path.get_sequence() << std::endl
        << "Ref1: " << reference_1 << std::endl
        << "Ref2: " << reference_2 << std::endl;
    // TODO: what about other cases?
    EXPECT_EQ("8=3X4=", path.get_cigar().to_string());
    EXPECT_EQ(12u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAA" "C" "TCTCTG" "A" "A" "C" "TTGT";
    std::string query =     "ACGCAA" "T" "TCTCTG" "T" "A" "T" "TTGT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                     reference.begin()),
              path.get_score());
    EXPECT_EQ("6=1X6=1X1=1X4=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAA" "CTTTTTT";
    std::string reference_2 = "AAAA" "TTGGGGG";
    std::string query =       "AAAA" "TTTTTTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    config.num_alternative_paths = 2;
    DBGAligner<> aligner(*graph, config);

    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(2u, alt_paths.size());
    EXPECT_NE(alt_paths.front(), alt_paths.back());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                     reference_1.begin()),
              path.get_score());
    EXPECT_EQ("4=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAA" "TTTT" "TTTT",
                                           "ACAA" "TTTT" "TGTT",
                                           "ACAA" "GTTT" "TTTT",
                                           "ACAA" "GTTT" "TGTT"};
    std::string query =                    "ACAA" "CTTT" "TCTT";

    auto graph = build_graph_batch<TypeParam>(k, references);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    config.num_alternative_paths = 2;
    config.queue_size = 100;
    DBGAligner<> aligner(*graph, config);

    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(config.num_alternative_paths, alt_paths.size());
    for (const auto &path : alt_paths) {
        EXPECT_EQ("4=1X4=1X2=", path.get_cigar().to_string())
            << query << "\n" << path.get_sequence();
        EXPECT_EQ(10u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             alt_paths.get_query(),
                             alt_paths.get_query_reverse_complement());
    }

    // TODO check with extend_mapping
}

TYPED_TEST(DBGAlignerTest, align_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                     reference.begin()),
              path.get_score());
    EXPECT_EQ("4=1X9=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(19u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TYPED_TEST(DBGAlignerTest, align_multiple_misalignment_bandwidth) {
    size_t k = 4;
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });

    for (uint64_t bandwidth : std::vector<uint64_t>{ 2, 5, 10, std::numeric_limits<uint64_t>::max()}) {
        auto config_bandwidth = config;
        config_bandwidth.bandwidth = bandwidth;

        DBGAligner<> aligner(*graph, config_bandwidth);
        auto paths = aligner.align(query);
        ASSERT_FALSE(paths.empty());

        EXPECT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(query.size() - k + 1, path.size());
        EXPECT_EQ(reference, path.get_sequence());
        EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                         reference.begin()),
                  path.get_score());
        EXPECT_EQ("4=1X9=1X6=", path.get_cigar().to_string());
        EXPECT_EQ(19u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
    }
}

TYPED_TEST(DBGAlignerTest, align_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCC"     "TTGTT";
    std::string query =     "TTTCC" "A" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference.begin(), reference.end())
                + config.gap_opening_penalty,
              path.get_score());
    EXPECT_EQ("5=1I5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TYPED_TEST(DBGAlignerTest, align_delete) {
    size_t k = 4;
    std::string reference = "TTCGA" "T" "TGGCCT";
    std::string query =     "TTCGA"     "TGGCCT";
    // alt query            "TTCGA" "T"  "GGCCT"

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    ASSERT_EQ(1ull, paths.size());
    auto &path = const_cast<DBGAligner<>::DBGAlignment&>(*paths.begin());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end())
                + config.gap_opening_penalty,
              path.get_score());

    std::unordered_set<std::string> possible_cigars { "6=1D5=", "5=1D6=" };

    EXPECT_NE(possible_cigars.end(), possible_cigars.find(path.get_cigar().to_string()));
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    Cigar cigar1, cigar2;
    cigar1.append(Cigar::Operator::MATCH, 6);
    cigar1.append(Cigar::Operator::DELETION, 1);
    cigar1.append(Cigar::Operator::MATCH, 5);

    cigar2.append(Cigar::Operator::MATCH, 5);
    cigar2.append(Cigar::Operator::DELETION, 1);
    cigar2.append(Cigar::Operator::MATCH, 6);

    path.set_cigar(std::move(cigar1));
    bool check_main = check_extend(graph, aligner.get_config(), paths, query);

    path.set_cigar(std::move(cigar2));
    bool check_alt = check_extend(graph, aligner.get_config(), paths, query);

    EXPECT_TRUE(check_main || check_alt);
}

TYPED_TEST(DBGAlignerTest, align_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATA" "CCTT" "GGCGCTCTC";
    std::string query =     "TTTCTGTATA"        "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end())
                + config.gap_opening_penalty
                + score_t(3) * config.gap_extension_penalty,
              path.get_score());
    EXPECT_EQ("10=4D9=", path.get_cigar().to_string());
    EXPECT_EQ(19u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TYPED_TEST(DBGAlignerTest, align_clipping1) {
    size_t k = 4;
    std::string reference = "GGCC" "TGTTTG";
    std::string query =     "ACCC" "TGTTTG";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    DBGAligner<> aligner(*graph, config);
    auto alt_paths = aligner.align(query);
    ASSERT_FALSE(alt_paths.empty());

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(5ull, path.size());
    EXPECT_EQ(reference.substr(2), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()),
              path.get_score());
    EXPECT_EQ("2S8=", path.get_cigar().to_string())
        << reference.substr(2) << " " << path.get_sequence();
    EXPECT_EQ(8u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         alt_paths.get_query(),
                         alt_paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), alt_paths, query));
}

TYPED_TEST(DBGAlignerTest, align_clipping2) {
    size_t k = 4;
    std::string reference = "AAA" "AGCTTCGAGGCCAA";
    std::string query =      "TT" "AGCTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(11u, path.size());
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()),
              path.get_score());
    EXPECT_EQ("2S14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TYPED_TEST(DBGAlignerTest, align_clipping_min_cell_score) {
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    DBGAlignerConfig config = ::config;
    config.min_cell_score = std::numeric_limits<score_t>::min();
    config.min_path_score = std::numeric_limits<score_t>::min();
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(7u, path.size());
    EXPECT_EQ(reference.substr(5), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()),
              path.get_score());
    EXPECT_EQ("2S13=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
}

TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    {
        DBGAlignerConfig config = ::config;
        config.min_seed_length = 2;
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min();
        config.min_path_score = std::numeric_limits<score_t>::min();
        DBGAligner<SuffixSeeder<>> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(7u, path.size());
        EXPECT_EQ(reference.substr(5), path.get_sequence());
        EXPECT_EQ(config.match_score(query.begin() + 2, query.end()),
                  path.get_score());
        EXPECT_EQ("2S13=", path.get_cigar().to_string());
        EXPECT_EQ(13u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(2u, path.get_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        EXPECT_TRUE(check_extend(graph, aligner.get_config(), paths, query));
    }
    {
        DBGAlignerConfig config = ::config;
        config.min_seed_length = 1;
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min();
        config.min_path_score = std::numeric_limits<score_t>::min();
        DBGAligner<SuffixSeeder<>> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(15u, path.size()); // includes dummy k-mers
        EXPECT_EQ(reference.substr(3), path.get_sequence());
        EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                         reference.begin() + 3),
                  path.get_score());
        EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
        EXPECT_EQ(14u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(6u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        // the unimem alignment mode skips partial k-mer matches in the beginning
        EXPECT_FALSE(check_extend(graph, aligner.get_config(), paths, query));
    }
}

TEST(DBGAlignerTest, align_suffix_seed_snp) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    DBGAlignerConfig config = ::config;
    config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
    config.min_cell_score = std::numeric_limits<score_t>::min();
    config.min_path_score = std::numeric_limits<score_t>::min();
    DBGAligner<SuffixSeeder<>> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(15u, path.size()); // includes dummy k-mers
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(),
                                     reference.begin() + 3),
              path.get_score());
    EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(6u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    // the unimem alignment mode skips partial k-mer matches in the beginning
    EXPECT_FALSE(check_extend(graph, aligner.get_config(), paths, query));
}
