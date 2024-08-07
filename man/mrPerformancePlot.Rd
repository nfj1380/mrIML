% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrPerformancePlot.R
\name{mrPerformancePlot}
\alias{mrPerformancePlot}
\title{Plot Model Performance Comparison}
\usage{
mrPerformancePlot(
  ModelPerf1 = NULL,
  ModelPerf2 = NULL,
  mod_names = c("combined", "Xonly_model"),
  mode = "classification"
)
}
\arguments{
\item{ModelPerf1}{Dataframe of model performance metrics for the first model to compare.}

\item{ModelPerf2}{Dataframe of model performance metrics for the second model to compare.}

\item{mod_names}{Character vector of model names. Default is \code{c('combined', 'Xonly_model')}.}

\item{mode}{Character string indicating whether the mode is 'classification' or 'regression'. Default is 'classification'.}
}
\value{
A list containing:
\item{p1}{A ggplot object for the boxplot of model performance metrics.}
\item{p2}{A ggplot object for the barplot of differences in performance metrics.}
\item{wide_df}{A dataframe with the wide format of model performance metrics and their differences.}
}
\description{
Create visualizations to compare the performance of two models based on their performance metrics.
}
\examples{
plots <- mrPerformancePlot(ModelPerf1 =ModelPerf_lm, ModelPerf2 = ModelPerf_rf, mod_names=c('linear_reg','rand_forest'), mode='regression' )

}
