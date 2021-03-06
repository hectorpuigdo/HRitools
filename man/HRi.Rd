% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRi.R
\name{HRi}
\alias{HRi}
\title{HRi estimators}
\usage{
HRi(x)
}
\arguments{
\item{x}{Dataframe consisting of recombination, adaptation and 0-fold sites count data.}
}
\value{
A summary of the analysis of the dataframe consisting of: 
the curvilinear model, 
LHRi of the model, 
the cutoff value where LHRi is below 0.05 from the initial LHRi value,
the Mean Ideal Adaptation or MIA,
the Maximum Detrimental Value or MDV,
the Curvature of the function or Recovery Strength,
the formula of the model as a Formula object,
the remaining LHRi values when recombination increases,
a graph with a representation of these values,
a general graph with the model
and a NLS object with the results of appliying the NLS computation to the data.
}
\description{
It allows you to calculate HRi curvilinear functions over an adapatation-recombination dataset ponderating with 0-fold sites.
}
\details{
This is a first release. It might not work with really complicated datasets, as initial parameters still need to be improved.
}
\keyword{Hill-Robertson}
