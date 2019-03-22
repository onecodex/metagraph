//
// Created by Jan Studený on 2019-03-11.
//
#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <ProgressBar.hpp>
#include <tclap/CmdLine.h>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;

using json = nlohmann::json;
#define _DNA_GRAPH 1

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"

#include "dbg_succinct.hpp"
#include "sequence_graph.hpp"
#include "dbg_succinct_construct.hpp"
#include "dbg_hash.hpp"


#include "path_database_list_of_bifurcation_choices.hpp"
#include "samplers.hpp"
#include "utilities.hpp"

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;



int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> inputArg("i",
                                         "input",
                                         "FASTA/Q file that should be compressed",
                                         true,
                                         "",
                                         "string");
    TCLAP::ValueArg<std::string> statsArg("s",
                                          "statistics",
                                          "Filename of json file that will output statistics about compressed file.",
                                          true,
                                          "statistics.json",
                                          "string");
    TCLAP::ValueArg<std::string> decompressedArg("d",
                                          "decompressed_file",
                                          "[for debugging purposes] Output the reads back in FASTA format.",
                                          false,
                                          "",
                                          "string");
    cmd.add(inputArg);
    cmd.add(statsArg);
    cmd.add(decompressedArg);
    cmd.parse(argc, argv);
    auto input_filename = inputArg.getValue();
    auto statistics_filename = statsArg.getValue();
    auto reads = read_reads_from_fasta(input_filename);
    auto db = PathDatabaseListBC(reads);
    db.encode(reads);
    auto statistics = db.get_statistics();
    save_string(statistics.dump(4),statistics_filename);
    if (decompressedArg.isSet()) {
        auto decompressed_filename = decompressedArg.getValue();
        write_reads_to_fasta(db.get_all_reads(),decompressed_filename);
    }

    return 0;
}
