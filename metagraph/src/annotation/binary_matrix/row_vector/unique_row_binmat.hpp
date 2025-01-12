#ifndef __UNIQUE_ROW_BINMAT_HPP__
#define __UNIQUE_ROW_BINMAT_HPP__

#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace binmat {

class UniqueRowBinmat : public BinaryMatrix {
  public:
    explicit UniqueRowBinmat(uint64_t num_rows = 0);

    UniqueRowBinmat(std::vector<SetBitPositions>&& unique_rows,
                    std::vector<uint32_t>&& row_rank,
                    uint32_t num_columns);

    UniqueRowBinmat(const std::function<void(const RowCallback &)> &call_rows,
                    uint32_t num_columns);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return row_rank_.size(); }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const { return num_relations_; }
    // matrix density
    double density() const;

  private:
    uint32_t num_columns_ = 0;
    uint32_t num_relations_ = 0;
    std::vector<SetBitPositions> unique_rows_;
    std::vector<uint32_t> row_rank_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __UNIQUE_ROW_BINMAT_HPP__
