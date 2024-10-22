% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipa_fisher.test.R
\name{lipa_fisher.test}
\alias{lipa_fisher.test}
\title{lipa_fisher.test}
\usage{
lipa_fisher.test(output, counts, alt = "two.sided")
}
\arguments{
\item{output}{Output of lipa_expressionmarkers}

\item{counts}{Output of originalcount}

\item{alt}{"two.sided", "greater", "less"}
}
\value{
A dataframe of the fisher's test output for Lipa in each
subpopulations of GFP+/GFP-.
}
\description{
(LipaTg) Fisher’s test on each of the clusters
}
\examples{
lipa_fisher.test(lipa.40pc_0.2, counts_40pc_0.2)



}
