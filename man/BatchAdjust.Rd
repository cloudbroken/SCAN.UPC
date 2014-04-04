\name{BatchAdjust}
\alias{BatchAdjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper function for applying the ComBat batch-adjusting methodology to an ExpressionSet object
}
\description{
  TBD
}
\usage{
BatchAdjust(expressionSet, batchVariableName, covariateVariableNames=c())
}
\arguments{
  \item{expressionSet}{
    TBD
}
  \item{batchVariableName}{
    TBD
}
  \item{covariateVariableNames}{
    TBD
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An ExpressionSet object indicating expression levels after batch adjustment.
}
\references{
  TBD
}
\author{
  Stephen R. Piccolo
}
\note{
  TBD
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%\dontrun{
%pkgName = InstallBrainArrayPackage(celFilePath, "17.0.1", "hs", "entrezg")
%}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
