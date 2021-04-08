#include "tax_classify.hpp"

#include "annotation/taxonomy/tax_classifier.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

uint64_t SIZE_EXECUTE_BATCH = 10000;

void execute_fasta_seq(std::string curr_seq,
                         std::string curr_seq_name,
                         const mtg::graph::DeBruijnGraph &graph,
                         const mtg::annot::TaxClassifier &tax_classifier,
                         std::vector<std::pair<std::string, uint64_t> > &results,
                         std::mutex &result_mutex,
                         const Config &config) {
    uint64_t taxid = tax_classifier.assign_class(graph,
                                                 curr_seq,
                                                 config.lca_coverage_threshold,
                                                 config.discovery_fraction);
    std::unique_lock<std::mutex> lock(result_mutex);
    results.push_back({curr_seq_name, taxid});
}

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
                        const mtg::annot::TaxClassifier &tax_classifier,
                        std::vector<std::pair<std::string, uint64_t> > &results,
                        const Config &config) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file);
    std::mutex result_mutex;

    std::vector<std::pair<std::string, std::string>> seq_batch;

    uint64_t cnt_queries_executed = 0;
    for (const seq_io::kseq_t &kseq : fasta_parser) {
        seq_batch.push_back({std::string(kseq.seq.s), std::string(kseq.name.s)});

        if (seq_batch.size() != SIZE_EXECUTE_BATCH) {
            continue;
        }
        thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string>> sequences){
            std::cerr << "start batch of size " << sequences.size() << std::endl;
            for (std::pair<std::string, std::string> &kseq: sequences) {
                execute_fasta_seq(kseq.first, kseq.second, graph, tax_classifier, results, result_mutex, config);
            }
        }, std::move(seq_batch));

        cnt_queries_executed += 1;
        if (cnt_queries_executed % 100000 == 0) {
            logger->trace("Started to executed {} queries.", cnt_queries_executed);
        }
        seq_batch.clear();
    }

    thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string>> sequences){
      for (std::pair<std::string, std::string> &kseq: sequences) {
          execute_fasta_seq(kseq.first, kseq.second, graph, tax_classifier, results, result_mutex, config);
      }
    }, std::move(seq_batch));
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
    std::vector<std::pair<std::string, uint64_t> > results;
    for (const std::string &file: files) {
        logger->trace("Start processing file '{}'.", file);
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
                           tax_classifier,
                           results,
                           *config);
    }
    thread_pool.join();
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
