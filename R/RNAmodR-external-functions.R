#' @include RNAmodR.R
NULL

# for validation
.valid.SimpleSplitDataFrameList <- IRanges:::.valid.SimpleSplitDataFrameList
.valid.DataFrame <- S4Vectors:::.valid.DataFrame

# Gviz
.sequenceTrackInfo <- Gviz:::.sequenceTrackInfo

# S4Vectors
make_rownames_for_RectangularData_display  <- S4Vectors:::make_rownames_for_RectangularData_display
Vector_window <- S4Vectors:::Vector_window
labeledLine <- S4Vectors:::labeledLine
make_zero_col_DataFrame <- S4Vectors:::make_zero_col_DataFrame
set_unlisted_names <- S4Vectors:::set_unlisted_names
