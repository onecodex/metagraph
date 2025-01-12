#include <json/json.h>
#include <server_http.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "query.hpp"
#include "align.hpp"
#include "server_utils.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::graph::AnnotatedDBG;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

const std::string SEQ_DESCRIPTION_JSON_FIELD = "seq_description";
const std::string SCORE_JSON_FIELD = "score";
const std::string SEQUENCE_JSON_FIELD = "sequence";
const std::string ALIGNMENT_JSON_FIELD = "alignments";
const std::string CIGAR_JSON_FIELD = "cigar";

// convert values into proper types, i.e. 'nan' -> null, strings representing numbers -> numbers
Json::Value adjust_for_types(const std::string &v) {
    if (v == "nan")
        return Json::nullValue;

    try {
        return Json::Value(std::stof(v));
    } catch(...) {}

    return Json::Value(v);
}

std::string convert_query_response_to_json(const std::string &ret_str) {
    // TODO: we are parsing back the string generated by the 'query' code, which is ugly.
    // we should have an intermediate representation which can be converted to a string (when
    // query is invoked from the command line) or to a json (string) when invoked by the server.
    std::vector<std::string> queries = utils::split_string(ret_str, "\n");

    std::vector<std::pair<long long, Json::Value>> query_results;
    query_results.reserve(query_results.size());

    for (auto qit = queries.begin(); qit != queries.end(); ++qit) {
        std::vector<std::string> parts = utils::split_string(*qit, "\t", false);

        if (parts.size() <= 2)
            continue; // no sequences found

        Json::Value res_obj;
        std::vector<std::string> query_desc_parts = utils::split_string(parts[1], ":");

        res_obj[SEQ_DESCRIPTION_JSON_FIELD] = "";
        if (!query_desc_parts.empty())
            res_obj[SEQ_DESCRIPTION_JSON_FIELD] = query_desc_parts[0];

        if (query_desc_parts.size() > 1) {
            // we aligned first, so extracting aligned sequence and score:

            res_obj[SEQUENCE_JSON_FIELD] = query_desc_parts[1];
            res_obj[SCORE_JSON_FIELD] = (int)atoi(query_desc_parts[2].c_str());
            res_obj[CIGAR_JSON_FIELD] = query_desc_parts[3];
        }

        res_obj["results"] = Json::Value(Json::arrayValue);

        for (size_t i = 2; i < parts.size(); ++i) {
            Json::Value sampleEntry;

            std::vector<std::string> entries = utils::split_string(parts[i], ":");

            std::vector<std::string> labels
                = utils::split_string(entries[0].substr(1, entries[0].size() - 2), ";");

            sampleEntry["sample"] = labels[0];

            Json::Value properties = Json::objectValue;

            for (auto lit = ++labels.begin(); lit != labels.end(); ++lit) {
                std::vector<std::string> key_value = utils::split_string(*lit, "=");
                properties[key_value[0]] = adjust_for_types(key_value[1]);
            }

            if (properties.size() > 0) {
                sampleEntry["properties"] = properties;
            }
            sampleEntry["kmer_count"] = (int)atoi(entries[1].c_str());

            res_obj["results"].append(sampleEntry);
        }

        query_results.emplace_back(atoll(parts[0].c_str()), std::move(res_obj));
    }

    std::sort(query_results.begin(), query_results.end(), utils::LessFirst());

    Json::Value root = Json::Value(Json::arrayValue);
    // output by query id
    for (const auto &[id, res_obj] : query_results) {
        root.append(res_obj);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}


std::string process_search_request(const std::string &received_message,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const Config &config_orig) {
    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];
    if (fasta.isNull())
        throw std::domain_error("No input sequences received from client");

    Config config(config_orig);
    // discovery_fraction a proxy of 1 - %similarity
    config.discovery_fraction
            = json.get("discovery_fraction", config.discovery_fraction).asDouble();

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    if (config.discovery_fraction < 0.0 || config.discovery_fraction > 1.0) {
        throw std::domain_error(
                "Discovery fraction should be within [0, 1.0]. Instead got "
                + std::to_string(config.discovery_fraction));
    }

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::domain_error(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.count_labels = true;
    config.num_top_labels = json.get("num_labels", config.num_top_labels).asInt();
    config.fast = json.get("fast", config.fast).asBool();

    std::unique_ptr<graph::align::DBGAlignerConfig> aligner_config;
    if (json.get("align", false).asBool()) {
        aligner_config.reset(new graph::align::DBGAlignerConfig(
            initialize_aligner_config(anno_graph.get_graph().get_k(), config)
        ));
    }

    std::ostringstream oss;
    std::mutex oss_mutex;

    // writing to temporary file in order to reuse query code. This is not optimal and
    // may turn out to be an issue in production. However, adapting FastaParser to
    // work on strings seems non-trivial. An alternative would be to use
    // read_fasta_from_string for non fast queries.
    utils::TempFile tf(config.tmp_dir);
    tf.ofstream() << fasta.asString();
    tf.ofstream().close();

    // dummy pool doing everything in the caller thread
    ThreadPool dummy_pool(0);
    QueryExecutor engine(config, anno_graph, std::move(aligner_config), dummy_pool);

    engine.query_fasta(tf.name(),
        [&](const std::string &res) {
            std::lock_guard<std::mutex> lock(oss_mutex);
            oss << res;
        }
    );

    return convert_query_response_to_json(oss.str());
}

std::string process_align_request(const std::string &received_message,
                                  const graph::DeBruijnGraph &graph,
                                  const Config &config_orig) {
    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];

    Json::Value root = Json::Value(Json::arrayValue);

    Config config(config_orig);

    config.alignment_num_alternative_paths = json.get(
        "max_alternative_alignments",
        (uint64_t)config.alignment_num_alternative_paths).asInt();

    if (!config.alignment_num_alternative_paths) {
        // TODO: better throw an exception and send an error response to the client
        logger->warn("[Server] Got invalid value of alignment_num_alternative_paths = {}."
                     " The default value of 1 will be used instead...", config.alignment_num_alternative_paths);
        config.alignment_num_alternative_paths = 1;
    }

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::domain_error(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    std::unique_ptr<graph::align::IDBGAligner> aligner = build_aligner(graph, config);

    // TODO: make parallel?
    seq_io::read_fasta_from_string(fasta.asString(),
                                   [&](seq_io::kseq_t *read_stream) {
        const auto paths = aligner->align(read_stream->seq.s);

        Json::Value align_entry;
        align_entry[SEQ_DESCRIPTION_JSON_FIELD] = read_stream->name.s;

        // not supporting reverse complement yet
        Json::Value alignments = Json::Value(Json::arrayValue);

        for (const auto &path : paths) {
            Json::Value a;
            a[SCORE_JSON_FIELD] = path.get_score();
            a[SEQUENCE_JSON_FIELD] = path.get_sequence();
            a[CIGAR_JSON_FIELD] = path.get_cigar().to_string();

            alignments.append(a);
        };

        align_entry[ALIGNMENT_JSON_FIELD] = alignments;

        root.append(align_entry);
    });

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::string process_column_label_request(const graph::AnnotatedDBG &anno_graph) {
    auto labels = anno_graph.get_annotation().get_all_labels();

    Json::Value root = Json::Value(Json::arrayValue);

    for (const std::string &label : labels) {
        Json::Value entry = label;
        root.append(entry);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::string process_stats_request(const graph::AnnotatedDBG &anno_graph,
                                  const std::string &graph_filename,
                                  const std::string &annotation_filename) {
    Json::Value root;

    Json::Value graph_stats;
    graph_stats["filename"] = std::filesystem::path(graph_filename).filename().string();
    graph_stats["k"] = static_cast<uint64_t>(anno_graph.get_graph().get_k());
    graph_stats["nodes"] = anno_graph.get_graph().num_nodes();
    graph_stats["is_canonical_mode"] = (anno_graph.get_graph().get_mode()
                                            == graph::DeBruijnGraph::CANONICAL);
    root["graph"] = graph_stats;

    Json::Value annotation_stats;
    const auto &annotation = anno_graph.get_annotation();
    annotation_stats["filename"] = std::filesystem::path(annotation_filename).filename().string();
    annotation_stats["labels"] = static_cast<uint64_t>(annotation.num_labels());
    annotation_stats["objects"] = static_cast<uint64_t>(annotation.num_objects());
    annotation_stats["relations"] = static_cast<uint64_t>(annotation.num_relations());

    root["annotation"] = annotation_stats;

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::thread start_server(HttpServer &server_startup, Config &config) {
    server_startup.config.thread_pool_size = std::max(1u, get_num_threads());

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;

    logger->info("[Server] Will listen on {} port {}", config.host_address,
                 server_startup.config.port);
    return std::thread([&server_startup]() { server_startup.start(); });
}

template<typename T>
bool check_data_ready(std::shared_future<T> &data, shared_ptr<HttpServer::Response> response) {
    if (data.wait_for(0s) != std::future_status::ready) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        response->write(SimpleWeb::StatusCode::server_error_service_unavailable,
                        "Server is currently initializing, please come back later.");
        return false;
    }

    return true;
}

int run_server(Config *config) {
    assert(config);

    assert(config->infbase_annotators.size() == 1);

    ThreadPool graph_loader(1, 1);

    logger->info("[Server] Loading graph...");

    auto anno_graph = graph_loader.enqueue([&]() {
        auto graph = load_critical_dbg(config->infbase);
        logger->info("[Server] Graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);

        auto anno_graph = initialize_annotated_dbg(graph, *config);
        logger->info("[Server] Annotated graph loaded too. Current mem usage: {} MiB", get_curr_RSS() >> 20);
        return anno_graph;
    });

    // defaults for the server
    config->num_top_labels = 10000;
    config->fast = true;

    // the actual server
    HttpServer server;
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &content) {
                return process_search_request(content, *anno_graph.get(), *config);
            });
        }
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &content) {
                return process_align_request(content, anno_graph.get()->get_graph(), *config);
            });
        }
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &) {
                return process_column_label_request(*anno_graph.get());
            });
        }
    };

    server.resource["^/stats"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                            shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &) {
                return process_stats_request(*anno_graph.get(), config->infbase,
                                             config->infbase_annotators.front());
            });
        }
    };

    server.default_resource["GET"] = [](shared_ptr<HttpServer::Response> response,
                                        shared_ptr<HttpServer::Request> request) {
        logger->info("Not found " + request->path);
        response->write(SimpleWeb::StatusCode::client_error_not_found,
                        "Could not find path " + request->path);
    };
    server.default_resource["POST"] = server.default_resource["GET"];

    server.on_error = [](shared_ptr<HttpServer::Request> /*request*/,
                         const SimpleWeb::error_code &ec) {
        // Handle errors here, ignoring a few trivial ones.
        if (ec.value() != asio::stream_errc::eof
                && ec.value() != asio::error::operation_aborted) {
            logger->info("[Server] Got error {} {} {}",
                         ec.message(), ec.category().name(), ec.value());
        }
    };

    std::thread server_thread = start_server(server, *config);
    server_thread.join();

    return 0;
}

} // namespace cli
} // namespace mtg
