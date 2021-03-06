% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/renameSignature.R
\name{renameSignature}
\alias{renameSignature}
\title{Rename selected signatures from an existing set of signatures.}
\usage{
renameSignature(sigmatrix,selectSig,renameSig)
}
\arguments{
\item{sigmatrix}{An object of class \emph{mutsig} describing the existing set of signatures.}

\item{selectSig}{A vector containing name(s) of the signature(s) to be removed.}

\item{renameSig}{A vector containing updated name(s) of the signatures.}
}
\value{
An object of \emph{mutsig} class.
}
\description{
Rename selected signatures from an existing signature-set.
}
\details{
Renaming of selected mutation signatures to an existing signature set, upon which the updated signature set is returned.
}
\examples{
renamedSigMatrix=renameSignature(sigmatrix=signatures.cosmic, selectSig=c(2,5), 
renameSig=c('two', 'five')) # rename signature 2 and 5.
}
\seealso{
\code{\link{addSignature}}, \code{\link{deleteSignature}} and \code{\link{mergeSignature}}.
}
\keyword{mutation}
\keyword{rename}
\keyword{signatures.}
