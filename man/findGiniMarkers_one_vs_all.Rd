% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression.R
\name{findGiniMarkers_one_vs_all}
\alias{findGiniMarkers_one_vs_all}
\title{findGiniMarkers_one_vs_all}
\usage{
findGiniMarkers_one_vs_all(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  min_expr_gini_score = 0.5,
  min_det_gini_score = 0.5,
  detection_threshold = 0,
  rank_score = 1,
  min_feats = 4,
  min_genes = NULL,
  verbose = TRUE
)
}
\arguments{
\item{gobject}{giotto object}

\item{feat_type}{feature type}

\item{spat_unit}{spatial unit}

\item{expression_values}{feat expression values to use}

\item{cluster_column}{clusters to use}

\item{subset_clusters}{selection of clusters to compare}

\item{min_expr_gini_score}{filter on minimum gini coefficient on expression}

\item{min_det_gini_score}{filter on minimum gini coefficient on detection}

\item{detection_threshold}{detection threshold for feat expression}

\item{rank_score}{rank scores for both detection and expression to include}

\item{min_feats}{minimum number of top feats to return}

\item{min_genes}{deprecated, use min_feats}

\item{verbose}{be verbose}
}
\value{
data.table with marker feats
}
\description{
Identify marker feats for all clusters in a one vs all manner
based on gini detection and expression scores.
}
\examples{
g <- GiottoData::loadGiottoMini("visium")

findGiniMarkers_one_vs_all(g, cluster_column = "leiden_clus")
}
\seealso{
\code{\link{findGiniMarkers}}
}
