% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wnn.R
\name{runWNN}
\alias{runWNN}
\title{Multi omics integration with WNN}
\usage{
runWNN(
  gobject,
  spat_unit = "cell",
  feat_types = c("rna", "protein"),
  reduction_methods = c("pca", "pca"),
  reduction_names = c("rna.pca", "protein.pca"),
  k = 20,
  integrated_feat_type = NULL,
  matrix_result_name = NULL,
  w_names = c(NULL, NULL),
  verbose = FALSE
)
}
\arguments{
\item{gobject}{A Giotto object with individual PCA feat_types pre-calculated}

\item{spat_unit}{spatial unit}

\item{feat_types}{feature types to integrate. Default = c("rna", "protein")}

\item{reduction_methods}{reduction methods for each feature type. Default = c("pca", "pca")}

\item{reduction_names}{names of the reduction methods to use. Default = c("rna.pca", "protein.pca")}

\item{k}{number of neighbors to calculate cell distances, default = 20}

\item{integrated_feat_type}{integrated feature type (e.g. 'rna_protein')}

\item{matrix_result_name}{Default = 'theta_weighted_matrix'}

\item{w_names}{optional. Names for the weighted matrices. If NULL, automatic names composed by w_feat_type will be created.}

\item{verbose}{be verbose}
}
\value{
A Giotto object with a new multiomics slot containing the theta_weighted_matrix and individual weight matrices.
}
\description{
Multi omics integration with WNN
}
