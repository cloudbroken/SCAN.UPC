\name{SCAN}
\alias{SCAN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Single-Channel Array Normalization (SCAN) for two-channel microarrays
}
\description{
  This function is used to normalize two-channel expression microarrays (such as those from Agilent) via the SCAN method. In raw form, such microarray data come in the form of tab-separate data files.
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
  A list object that contains two values: a row for each probeset/gene/transcript and a column for each input file.
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
}
\author{
  Stephen R. Piccolo
}
\note{
By default, \code{SCAN} uses the default mappings between probes and genes that have been provided by the manufacturer. However, these mappings may be outdated or may include problematic probes (for example, those that cross hybridize). The default mappings also may produce multiple summary values per gene. Alternative mappings, such as those provided by the BrainArray resource (see \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp}), allow \code{SCAN} to produce a single value per gene and to use updated gene definitions. Users can specify alternative mappings using the \code{probeSummaryPackage} parameter. If specified, this package must conform to the standards of the \code{AnnotationDbi} package. The BrainArray packages can be downloaded from \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp}. When using BrainArray, be sure to download the R source package for probe-level mappings (see vignette for more information).
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
# Download an example CEL file and save it as a temporary local file
celFilePath = file.path(tempdir(), "Vignette_Example.CEL.gz")
download.file("http://www.ncbi.nlm.nih.gov/geosuppl/?acc=GSM555237&file=GSM555237.CEL.gz", celFilePath)

# Normalize a CEL file
normalized = SCAN(celFilePath)

# Normalize a CEL file and save output to a file
normalized = SCAN(celFilePath, "output_file.txt")

\dontrun{
# Normalize a CEL file and summarize at the gene level using BrainArray
# mappings for Entrez Gene
probeFilePath = file.path(tempdir(), "hgu95ahsentrezgprobe_15.0.0.tar.gz")
download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/15.0.0/entrezg.download/hgu95ahsentrezgprobe_15.0.0.tar.gz", probeFilePath)
install.packages(probeFilePath, repos=NULL, type="source")
library(hgu95ahsentrezgprobe)
normalized = SCAN(celFilePath, probeSummaryPackage=hgu95ahsentrezgprobe)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
