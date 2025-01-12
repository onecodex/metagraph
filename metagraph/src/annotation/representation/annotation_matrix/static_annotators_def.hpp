#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include "annotation_matrix.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"


namespace mtg {
namespace annot {

typedef StaticBinRelAnnotator<binmat::RowConcatenated<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<binmat::RowSparse, std::string> RowSparseAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbowfish, std::string> RainbowfishAnnotator;

typedef StaticBinRelAnnotator<binmat::BRWT, std::string> MultiBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT_sdsl, std::string> BinRelWT_sdslAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT, std::string> BinRelWTAnnotator;

typedef StaticBinRelAnnotator<binmat::UniqueRowBinmat, std::string> UniqueRowAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbow<binmat::BRWT>, std::string> RbBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::ColumnMajor>, std::string> RowDiffColumnAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::BRWT>, std::string> RowDiffBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::RowSparse>, std::string> RowDiffRowSparseAnnotator;


template <>
inline const std::string RowFlatAnnotator::kExtension = ".flat.annodbg";
template <>
inline const std::string RowSparseAnnotator::kExtension = ".row_sparse.annodbg";
template <>
inline const std::string RainbowfishAnnotator::kExtension = ".rbfish.annodbg";
template <>
inline const std::string MultiBRWTAnnotator::kExtension = ".brwt.annodbg";
template <>
inline const std::string BinRelWT_sdslAnnotator::kExtension = ".bin_rel_wt_sdsl.annodbg";
template <>
inline const std::string BinRelWTAnnotator::kExtension = ".bin_rel_wt.annodbg";
template <>
inline const std::string UniqueRowAnnotator::kExtension = ".unique_row.annodbg";
template <>
inline const std::string RbBRWTAnnotator::kExtension = ".rb_brwt.annodbg";
template <>
inline const std::string RowDiffColumnAnnotator::kExtension = ".row_diff.annodbg";
template <>
inline const std::string RowDiffBRWTAnnotator::kExtension = ".row_diff_brwt.annodbg";
template <>
inline const std::string RowDiffRowSparseAnnotator::kExtension = ".row_diff_sparse.annodbg";
} // namespace annot
} // namespace mtg

#endif // __STATIC_ANNOTATOR_DEFS_HPP__
