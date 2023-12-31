% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spellLengths.R
\name{spellLengths}
\alias{spellLengths}
\title{Spell length calculation}
\usage{
spellLengths(s)
}
\arguments{
\item{s}{A binary vector of 0 dry and 1 wet only.}
}
\value{
Returns a list object containing a vector of dry spell lengths and a vector of wet spell lengths.
}
\description{
Function to calculate the length (duration in years) of wet or dry periods.
}
\examples{

 #use 0 for dry and 1 for wet years
 spells = c(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0)

 spellLengths(spells)

}
