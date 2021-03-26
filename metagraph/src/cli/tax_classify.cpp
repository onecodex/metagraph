#include "taxo_classify.hpp"

#include "annotation/taxonomy/taxo_classifier.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
<<<<<<< HEAD:metagraph/src/cli/tax_classify.cpp
                        const mtg::annot::TaxClassifier &tax_classifier,
=======
                        const mtg::annot::TaxoClassifier &taxo_classifier,
>>>>>>> fix the issues related to multithreading:metagraph/src/cli/taxo_classify.cpp
                        std::vector<std::pair<std::string, uint64_t> > &results,
                        const Config &config) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file);
    std::mutex result_mutex;

    for (const seq_io::kseq_t &kseq : fasta_parser) {
<<<<<<< HEAD:metagraph/src/cli/tax_classify.cpp
        std::string curr_seq = std::string(kseq.seq.s);
        std::string curr_seq_name = std::string(kseq.name.s);
        thread_pool.enqueue([&](std::string curr_seq, std::string curr_seq_name){
            uint64_t taxid = tax_classifier.assign_class(graph,
                                                          curr_seq,
                                                          config.lca_coverage_threshold);
            std::unique_lock<std::mutex> lock(result_mutex);
            results.push_back({curr_seq_name, taxid});
        }, std::move(curr_seq), std::move(curr_seq_name));
=======

        if (string(kseq.name.s).find("NC_031246.1|-137") != std::string::npos ) {
            cerr << "before thread_pool.enqueue " << kseq.name.s << "\n";
        }

        std::string curr_seq = std::string(kseq.seq.s);
        std::string curr_seq_name = std::string(kseq.name.s);
        thread_pool.enqueue([&](const double lca_coverage_threshold, std::string curr_seq, std::string curr_seq_name){
            bool debug = false;
            if (string(kseq.name.s).find("NC_031246.1|-137") != std::string::npos ) {
                cerr << "found 137 in " << kseq.name.s << "\n";
                debug = true;
            }

            uint64_t taxid = taxo_classifier.assign_class(graph,
                                                          curr_seq,
                                                          lca_coverage_threshold,
                                                          debug);
            std::unique_lock<std::mutex> lock(result_mutex);
            results.push_back({curr_seq_name, taxid});
        }, config.lca_coverage_threshold, std::move(curr_seq), std::move(curr_seq_name));
>>>>>>> fix the issues related to multithreading:metagraph/src/cli/taxo_classify.cpp
    }
}

int taxonomic_classification(Config *config) {
    assert(config);

    const std::vector<std::string> &files = config->fnames;

    Timer timer;
    logger->trace("Loading TaxonomyDB ...");
    mtg::annot::TaxClassifier tax_classifier(config->taxonomic_tree);
    logger->trace("Finished loading TaxonomyDB in {}", timer.elapsed());

    timer.reset();
    logger->trace("Graph loading ...");
    auto graph = load_critical_dbg(config->infbase);
    logger->trace("Finished graph loading in {}", timer.elapsed());

    timer.reset();
    logger->trace("Processing the classification");
    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);
<<<<<<< HEAD:metagraph/src/cli/tax_classify.cpp
=======

>>>>>>> fix the issues related to multithreading:metagraph/src/cli/taxo_classify.cpp
    std::vector<std::pair<std::string, uint64_t> > results;
    for (const std::string &file: files) {
        logger->trace("Start processing file '{}'.", file);
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
<<<<<<< HEAD:metagraph/src/cli/tax_classify.cpp
                           tax_classifier,
=======
                           taxo_classifier,
>>>>>>> fix the issues related to multithreading:metagraph/src/cli/taxo_classify.cpp
                           results,
                           *config);
    }
    thread_pool.join();
<<<<<<< HEAD:metagraph/src/cli/tax_classify.cpp
=======

    std::cerr << "results size = " << results.size() << "\n";

>>>>>>> fix the issues related to multithreading:metagraph/src/cli/taxo_classify.cpp
    for (const std::pair<std::string, uint64_t> &result: results) {
        std::string output = fmt::format("Sequence '{}' was classified with Tax ID '{}'\n",
                                         result.first, result.second);
        std::cout << output << std::endl;
    }
    logger->trace("Finished processing all the classification in {}", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
