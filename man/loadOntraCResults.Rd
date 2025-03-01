% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ONTraC_wrapper.R
\name{loadOntraCResults}
\alias{loadOntraCResults}
\title{loadOntraCResults}
\usage{
loadOntraCResults(
  gobject,
  ontrac_results_dir = getwd(),
  preprocessing_dir = file.path(ontrac_results_dir, "preprocessing_dir"),
  GNN_dir = file.path(ontrac_results_dir, "GNN_dir"),
  NTScore_dir = file.path(ontrac_results_dir, "NTScore_dir"),
  NTScore_reverse = FALSE
)
}
\arguments{
\item{gobject}{giotto object}

\item{ontrac_results_dir}{the directory where the ONTraC results are saved}

\item{preprocessing_dir}{the directory to save the preprocessing results.
Default is file.path(ontrac_results_dir, "preprocessing_dir")}

\item{GNN_dir}{the directory to save the GNN results. Default is
file.path(ontrac_results_dir, "GNN_dir")}

\item{NTScore_dir}{the directory to save the NTScore results. Default is
file.path(ontrac_results_dir, "NTScore_dir")}

\item{NTScore_reverse}{whether to reverse the NTScore. Default is FALSE}
}
\value{
gobject with ONTraC results
}
\description{
load ONTraC results
}
\details{
This function loads the ONTraC results into the giotto object.
}
