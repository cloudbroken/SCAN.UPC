\name{UPC_RNASeq}
\alias{UPC_RNASeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Universal Probability of expression Codes (UPC) for RNA-Seq data
}
\description{
This function is used to normalize RNA-Sequencing data via the UPC method. A necessary preliminary step is to generate read counts for each gene (or transcript or exon). A popular approach for accomplishing this is to use the Tophat short-read aligner (http://tophat.cbcb.umd.edu/), followed by application of the htseq-count tool (http://www-huber.embl.de/users/anders/HTSeq/).

Gene (or transcript or exon) values can be converted to UPCs using this function. The tab-separated data file should contain a row for each gene. The first column should contain the gene ID. The second column should contain the read counts (not RPKM/FPKM values). For example:

AAB    31
AAC    255

Most users will want to correct for length and GC content. In this case, a tab-separated annotation file that specifies length and number of GC bases should be included. The first column should contain the gene ID. The second column should contain the length of the gene. The third column should specify the number of number of G or C bases in the gene.

\begin{tabular}{lll}
AAB & 1767 & 640 \\
AAC & 654 & 333 \\
AAD & 4644 & 2039 \\
AAE & 2629 & 1011 \\
\end{tabular}
}
\usage{
UPC_RNASeq(inFilePattern, annotationFilePath, outFilePath = NA, modelType = "nn", convThreshold = 0.01, verbose = TRUE)
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

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
