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
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). (This is the only required parameter.)
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
By default, \code{UPC} uses the default mappings between probes and genes that have been provided by the manufacturer. However, these mappings may be outdated or may include problematic probes (for example, those that cross hybridize). The default mappings also may produce multiple summary values per gene. Alternative mappings, such as those provided by the BrainArray resource (see \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp}), allow \code{UPC} to produce a single value per gene and to use updated gene definitions. Users can specify alternative mappings using the \code{probeSummaryPackage} parameter. If specified, this package must conform to the standards of the \code{AnnotationDbi} package. The BrainArray packages can be downloaded from \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp}. When using BrainArray, be sure to download the R source package for probe-level mappings (see vignette for more information).
}

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
result = UPC_TwoColor(inFilePath, "output_file.txt")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
