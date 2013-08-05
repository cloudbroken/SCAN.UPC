\name{UPC_RNASeq}
\alias{UPC_RNASeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Universal Probability of expression Codes (UPC) for RNA-Seq data
}
\description{
This function is used to derive UPC values for RNA-Seq data. It requires an input file that specifies a read count for each genomic region (e.g., gene). This file should list a unique identifier for each region in the first column and corresponding read counts (not RPKM/FPKM values) in the second column.

This function also can correct for the GC content and length of each genomic region. Users who wish to enable this correction must provide a separate annotation file. This tab-separated file should contain a row for each genomic region. The first column should contain a unique identifier that corresponds to identifiers from the read-count input file. The second column should indicate the length of the genomic region. And the third column should specify the number of G or C bases in the region.
}
\usage{
UPC_RNASeq(inFilePattern, annotationFilePath, outFilePath = NA, modelType = "nn", convThreshold = 0.01, ignoreZeroes = FALSE, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). This is the only required parameter.
}
  \item{annotationFilePath}{
     Absolute or relative path where the annotation file is located. This parameter is optional.
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. This is optional.
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln'') and negative-binomial (``nb'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge, it may be useful to adjust this value. (This parameter is optional.)
}
  \item{ignoreZeroes}{
     Whether to ignore read counts equal to zero when performing UPC calculations. Default is FALSE.
}
  \item{verbose}{
     Whether to output more detailed status information as files are normalized. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An \code{ExpressionSet} object that contains a row for each probeset/gene/transcript and a column for each input file.
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
}
\author{
  Stephen R. Piccolo
}
\note{
  RNA-Seq data by nature have a lot of zero read counts. Samples with an excessive number of zeroes may lead to error messages because genes cannot be allocated properly to bins. The user can specify \code{ignoreZeroes=TRUE} to avoid this error. In practice, we have seen that the resulting UPC values are similar with either approach.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
upc3 = UPC_RNASeq("ReadCounts.txt", "Annotation.txt")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
