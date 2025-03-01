% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience_visiumHD.R
\name{importVisiumHD}
\alias{importVisiumHD}
\title{Import a Visium HD assay}
\usage{
importVisiumHD(
  visiumHD_dir = NULL,
  expression_source = "raw",
  gene_column_index = 2,
  barcodes = NULL,
  array_subset_row = NULL,
  array_subset_col = NULL,
  pxl_subset_row = NULL,
  pxl_subset_col = NULL,
  shape = "hexagon",
  shape_size = 400
)
}
\arguments{
\item{visiumHD_dir}{Visium HD output directory (e.g. square_016um)}

\item{expression_source}{character. Raw or filter expression data. Defaults 
to "raw"}

\item{gene_column_index}{numeric. Expression column to use for gene names
1 = Ensembl and 2 = gene symbols}

\item{barcodes}{character vector. (optional) Use if you only want to load
a subset of the pixel barcodes}

\item{array_subset_row}{numeric vector. (optional) Vector with min and max 
values to subset based on array rows}

\item{array_subset_col}{numeric vector. (optional) Vector with min and max 
values to subset based on array columns}

\item{pxl_subset_row}{numeric vector. (optional) Vector with min and max 
values to subset based on row pixels}

\item{pxl_subset_col}{numeric vector. (optional) Vector with min and max 
values to subset based on column pixels}
}
\value{
VisiumHDReader object
}
\description{
Giotto import functionalities for Visium HD datasets. This function generates
a `VisiumHDReader` instance that has convenient reader functions for 
converting individual pieces of Visium HD data into Giotto-compatible 
representations when the param `visiumHD_dir` is provided.
A function that creates the full `giotto` object is also available.
These functions should have all param values provided as defaults, but
can be flexibly modified to do things such as look in alternative
directories or paths.
}
\details{
Loading functions are generated after the `visiumHD_dir` is added.
}
\examples{
# Create a `VisiumHDReader` object
reader <- importVisiumHD()

\dontrun{
# Set the visiumHD_dir
reader$visiumHD_dir <- "path to visium HD dir"
readerHD$visiumHD_dir <- visiumHD_dir

# Load tissue positions or create cell metadata
tissue_pos <- readerHD$load_tissue_position()
metadata <- readerHD$load_metadata()

Load matrix or create expression object
matrix <- readerHD$load_matrix()
expression_obj <- readerHD$load_expression()

Load transcript data (cell metadata, expression object, and transcripts per 
pixel)
my_transcripts <- readerHD$load_transcripts(array_subset_row = c(500, 1000),
                                            array_subset_col = c(500, 1000))

# Create a `giotto` object and add the loaded data
# TODO
}
}
