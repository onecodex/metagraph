#include "annotation_matrix.hpp"

#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/serialization.hpp"
#include "static_annotators_def.hpp"


namespace mtg {
namespace annot {

using utils::remove_suffix;


template <class BinaryMatrixType, typename Label>
std::string
StaticBinRelAnnotator<BinaryMatrixType, Label>
::file_extension() const { return kExtension; }

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::has_label(Index i, const Label &label) const {
    assert(i < num_objects());
    try {
        return matrix_->get(i, label_encoder_.encode(label));
    } catch (...) {
        return false;
    }
}

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::has_labels(Index i, const VLabels &labels) const {
    assert(i < num_objects());

    std::set<size_t> querying_codes;
    try {
        for (const auto &label : labels) {
            querying_codes.insert(label_encoder_.encode(label));
        }
    } catch (...) {
        return false;
    }
    std::set<size_t> encoded_labels;
    for (uint64_t col_id : get_matrix().get_row(i)) {
        encoded_labels.insert(col_id);
    }
    return std::includes(encoded_labels.begin(), encoded_labels.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <class BinaryMatrixType, typename Label>
void
StaticBinRelAnnotator<BinaryMatrixType, Label>
::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension,
                            std::ios::binary);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }
    label_encoder_.serialize(outstream);
    matrix_->serialize(outstream);
}

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::merge_load(const std::vector<std::string> &filenames) {
    if (filenames.size() > 1)
        std::cerr << "Warning: Can't merge static annotators."
                     " Only the first will be loaded." << std::endl;

    std::ifstream instream(remove_suffix(filenames.at(0), kExtension) + kExtension,
                           std::ios::binary);
    if (!instream.good())
        return false;

    try {
        assert(matrix_.get());
        return label_encoder_.load(instream) && matrix_->load(instream);
    } catch (...) {
        return false;
    }
}

template <class BinaryMatrixType, typename Label>
uint64_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_objects() const {
    return matrix_->num_rows();
}

template <class BinaryMatrixType, typename Label>
uint64_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_relations() const {
    return matrix_->num_relations();
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>::except_dyn() {
    throw std::runtime_error("Dynamic actions are not supported"
                             " in static representation");
}

template <class BinaryMatrixType, typename Label>
bool StaticBinRelAnnotator<BinaryMatrixType, Label>
::dump_columns(const std::string &prefix, uint64_t num_threads) const {
    size_t m = this->num_labels();
    bool success = true;

    #pragma omp parallel for num_threads(num_threads)
    for (uint64_t i = 0; i < m; ++i) {
        std::ofstream outstream(
            remove_suffix(prefix, kExtension)
                + "." + std::to_string(i)
                + ".text.annodbg"
        );

        if (!outstream.good()) {
            std::cerr << "ERROR: dumping column " << i << " failed" << std::endl;
            success = false;
            continue;
        }

        outstream << num_objects() << " ";

        auto column = matrix_->get_column(i);

        outstream << column.size() << "\n";
        for (auto pos : column) {
            outstream << pos << "\n";
        }
    }

    return success;
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                        const LabelEncoder<Label> &label_encoder)
      : matrix_(std::move(matrix)) {
    assert(matrix_.get());
    label_encoder_ = label_encoder;
}

bool merge_load_row_diff(const std::vector<std::string> &filenames,
                         const ColumnCallback &callback,
                         size_t num_threads) {
    std::atomic<bool> error_occurred = false;

    std::vector<uint64_t> offsets(filenames.size() + 1, 0);

    // load labels
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < filenames.size(); ++i) {
        if (error_occurred)
            continue;

        std::ifstream instream(filenames[i], std::ios::binary);
        if (!instream.good()) {
            common::logger->error("Can't read from {}", filenames[i]);
            error_occurred = true;
        }

        LabelEncoder<std::string> label_encoder;
        if (!label_encoder.load(instream)) {
            common::logger->error("Can't load label encoder from {}", filenames[i]);
            error_occurred = true;
        }

        offsets[i + 1] = label_encoder.size();
    }

    // compute global offsets (partial sums)
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());

    std::vector<std::string> labels(offsets.back());

    // load annotations
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < filenames.size(); ++i) {
        if (error_occurred)
            continue;

        try {
            common::logger->trace("Loading annotations from file {}", filenames[i]);
            LabelEncoder<std::string> label_encoder;
            std::ifstream instream(filenames[i], std::ios::binary);
            binmat::RowDiff<binmat::ColumnMajor> matrix;
            if (!label_encoder.load(instream) || !matrix.load(instream)) {
                common::logger->error("Can't load {}", filenames[i]);
                std::exit(1);
            }
            if (!label_encoder.size())
                common::logger->warn("No labels in {}", filenames[i]);

            std::vector<std::unique_ptr<bit_vector>> cols = matrix.diffs().release_columns();

            for (uint32_t j = 0; j < cols.size(); ++j) {
                callback(offsets[i] + j, label_encoder.decode(j), std::move(cols[j]));
            }

        } catch (...) {
            error_occurred = true;
        }
    }

    return !error_occurred;
}

// template specialization of merge_load for RowDiffColumnAnnotator
template <>
bool StaticBinRelAnnotator<binmat::RowDiff<binmat::ColumnMajor>>
::merge_load(const std::vector<std::string> &filenames) {
    std::vector<std::unique_ptr<bit_vector>> columns;

    bool no_errors = true;
    bool merge_successful = merge_load_row_diff(
        filenames,
        [&](uint64_t, const std::string &label, std::unique_ptr<bit_vector> &&column) {
            uint64_t num_set_bits = column->num_set_bits();
            common::logger->trace("RowDiff column: {}, Density: {}, Set bits: {}", label,
                                  static_cast<double>(num_set_bits) / column->size(),
                                  num_set_bits);

            #pragma omp critical
            {
                if (columns.empty() || columns.back()->size() == column->size()) {
                    size_t col = label_encoder_.insert_and_encode(label);
                    if (col != columns.size()) {
                        common::logger->error("Duplicate columns {}", label);
                        no_errors = false;
                    }
                    columns.push_back(std::move(column));
                } else {
                    common::logger->error(
                            "Column {} has {} rows, previous column has {} rows",
                            label, columns.size(), columns.back()->size());
                    no_errors = false;
                }
            }
        },
        filenames.size() > 1u ? get_num_threads() : 0
    );

    matrix_ = std::make_unique<binmat::RowDiff<binmat::ColumnMajor>>(
            nullptr, binmat::ColumnMajor(std::move(columns)));

    return no_errors && merge_successful;
}

template class StaticBinRelAnnotator<binmat::RowConcatenated<>, std::string>;

template class StaticBinRelAnnotator<binmat::Rainbowfish, std::string>;

template class StaticBinRelAnnotator<binmat::BRWT, std::string>;

template class StaticBinRelAnnotator<binmat::BinRelWT_sdsl, std::string>;

template class StaticBinRelAnnotator<binmat::BinRelWT, std::string>;

template class StaticBinRelAnnotator<binmat::UniqueRowBinmat, std::string>;

template class StaticBinRelAnnotator<binmat::Rainbow<binmat::BRWT>, std::string>;

template class StaticBinRelAnnotator<binmat::RowDiff<binmat::ColumnMajor>, std::string>;

template class StaticBinRelAnnotator<binmat::RowDiff<binmat::BRWT>, std::string>;

template class StaticBinRelAnnotator<binmat::RowSparse, std::string>;

template class StaticBinRelAnnotator<binmat::RowDiff<binmat::RowSparse>, std::string>;

} // namespace annot
} // namespace mtg
