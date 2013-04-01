\name{SCAN_TwoColor}
\alias{SCAN_TwoColor}
\title{
  Single-Channel Array Normalization (SCAN) for two-channel Agilent microarrays
}
\description{
  This function is used to normalize two-channel expression microarrays (from Agilent) via the SCAN method. In raw form, such microarray data are tab-separate data files.
}
\usage{
SCAN_TwoColor(inFilePattern, outFilePath = NA, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). This is the only required parameter.
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. This is optional.
}
  \item{verbose}{
     Whether to output more detailed status information as files are normalized. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list is returned, containing two elements: a matrix containing normalized data values and a vector of probe names that correspond to each row of the matrix. The matrix will contain two columns---one corresponding to each channel---for each sample. When the array design contains duplicate probe names (this is common for control probes), the vector of probe names will also contain duplicates.
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
}
\author{
  Stephen R. Piccolo
}
%\note{
%  In most cases, individual probes are designed to target a specific genomic region that belongs within a particular gene.
%}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# Download an example file and save it as a temporary local file
inFilePath = file.path(tempdir(), "Vignette_Example.txt.gz")
download.file("http://www.ncbi.nlm.nih.gov/geosuppl/?acc=GSM1072833&file=GSM1072833%5F252747810638%5F020510%5FS01%5FCGH%5F105%5FDec08%5F1%5F2%2Etxt%2Egz", inFilePath)

# Normalize the file and save output to a file
result = SCAN_TwoColor(inFilePath, "output_file.txt")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
