% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModulonTA_functions.R
\name{target.analysis}
\alias{target.analysis}
\title{Target Analysis}
\usage{
target.analysis(net, modulons, cc)
}
\arguments{
\item{net}{A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.}

\item{modulons}{A list with as many elements as modulons/clusters containing the constituent elements.}

\item{cc}{Connected components generated with find.connected.components() or regulatory cores as the output of core()}
}
\value{
List of the modulon constituent elements sharing targets with a given modulon connected component.
}
\description{
Find common targets between a given modulon connected component and modulon constituent elements out of the connected component.
}
\details{
This function is searches for the modulon constituent elements sharing targets with a given modulon connected component, for each connected component and modulon in a given list of connected components split by their source modulon.
}
\examples{
\dontrun{
if(interactive()){
target.analysis(net = network.TILs,
mod = modulons.TILs,
cc = cc.TILs
)
}
}
}
