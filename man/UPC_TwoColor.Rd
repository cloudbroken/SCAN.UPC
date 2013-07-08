\name{UPC_TwoColor}
\alias{UPC_TwoColor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Universal Probability of expression Codes (UPC) for two-channel microarrays
}
\description{
  This function is used to normalize two-channel expression microarrays (from Agilent) using the Universal Probability of expression Codes (UPC) approach. In raw form, such microarray data come in the form of tab-separate data files.
}
\usage{
UPC_TwoColor(inFilePattern, outFilePath = NA, modelType="nn", convThreshold=0.01, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). Alternatively, a Gene Expression Omnibus identifier (e.g., GSE39655 or GSM1072833) can be specified.(This is the only required parameter.)
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. (This parameter is optional.)
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln'') and negative-binomial (``nb'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge, it may be useful to adjust this value. (This parameter is optional.)
}
  \item{verbose}{
     Whether to output more detailed status information as files are processed. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list is returned, containing two elements: a matrix containing UPC values and a vector of probe names that correspond to each row of the matrix. The matrix will contain two columns---one corresponding to each channel---for each sample. When the array design uses duplicate probe names (this is common for control probes), the vector of probe names will also contain duplicates.
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
}
\author{
  Stephen R. Piccolo
}
\note{
  If a Gene Expression Omnibus (GEO) identifier is specified for the \code{inFilePattern} parameter, an attempt will be made to download the sample(s) directly from GEO. If a study identifier (e.g., GSE39655) is specified, all CEL files from that study will be downloaded. If a sample identifier (e.g., GSM1072833) is specified, only that sample will be downloaded.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# Download an example file and save it as a temporary local file
tmpDir = tempdir()
getGEOSuppFiles("GSM1072833", makeDirectory=FALSE, baseDir=tmpDir)
inFilePath = file.path(tmpDir, "GSM1072833.CEL.gz")

# Normalize the file and save output to a file
result = UPC_TwoColor(inFilePath, "output_file.txt")

# Download a file directly from GEO and normalize it in one step
result = UPC_TwoColor("GSM1072833")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
