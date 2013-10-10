##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

shouldDownloadFromGEO = function(inFilePattern)
{
  substr(inFilePattern, 1, 3) %in% c("GSE", "GSM") & length(grep("\\.", inFilePattern)) == 0 & length(grep("\\*", inFilePattern)) == 0
}

downloadFromGEO = function(inFilePattern)
{
  tmpDir=tempdir()
  message(paste("Downloading ", inFilePattern, " directly from GEO to ", tmpDir, ".", sep=""))

  getGEOSuppFiles(inFilePattern, makeDirectory=FALSE, baseDir=tmpDir)

  if (substr(inFilePattern, 1, 3) == "GSE")
  {
    tarFilePath = file.path(tmpDir, paste(inFilePattern, "_RAW.tar", sep=""))

    if (!file.exists(tarFilePath))
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))

    individualDir = file.path(tmpDir, inFilePattern, "Files", sep="")
    dir.create(individualDir, recursive=TRUE)
    untar(tarFilePath, exdir=individualDir)
    inFilePattern = file.path(individualDir, "*", sep="")
  }

  if (substr(inFilePattern, 1, 3) == "GSM")
  {
    downloadedFiles = list.files(path=tmpDir, full.names=TRUE)

    if (length(downloadedFiles) == 0)
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))

    inFilePattern = file.path(tmpDir, basename(downloadedFiles))
  }

  inFilePattern
}

InstallBrainArrayPackage = function(celFilePath, version, organism, annotationSource)
{
  platform = cleancdfname(read.celfile.header(celFilePath, info="full")$cdfName)
  platform = sub("cdf", "", platform)

  packageName=paste(platform, organism, annotationSource, "probe", sep="")

  #if (!(packageName %in% rownames(installed.packages())))
  tmpDir = tempdir()

  packageFileName=paste(packageName, "_", version, ".tar.gz", sep="")
  tempPackageFilePath = paste(tmpDir, packageFileName, sep="")
  packageUrl = paste("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/", version, "/", annotationSource, ".download/", packageFileName, sep="")

  download.file(packageUrl, tempPackageFilePath)
  install.packages(tempPackageFilePath, repos=NULL, type="source")

  return(packageFileName)
}

ParseMetaFromGtfFile = function(gtfFilePath, fastaFilePattern, outFilePath, featureTypes=c("protein_coding"), attributeType="gene_id")
{
  fastaFilePaths = list.files(path=dirname(fastaFilePattern), pattern=glob2rx(basename(fastaFilePattern)), full.names=TRUE)

  if (length(fastaFilePaths) == 0)
    stop(paste("No FASTA file could be found at", fastaFilePattern))

  fastaStrings = readDNAStringSet(fastaFilePaths)
  fastaChromosomes = names(fastaStrings)

  message("Saving GTF data to temporary files")
  tmpDir = tempdir()
  chromosomes = NULL

  f = file(gtfFilePath, "r")
  stop = FALSE
  lineCount = 0

  while(!stop)
  {
    lineCount = lineCount + 1
    if (lineCount %% 10000 == 0)
      message(paste("Done parsing ", as.character(lineCount), " lines from ", gtfFilePath, sep=""))

    line = readLines(f, n = 1)

    if(length(line) == 0)
    {
      stop = TRUE
      next
    }

    lineItems = strsplit(line, "\t")[[1]]
    chromosome = lineItems[1]
    featureType = lineItems[2]
    subFeatureType = lineItems[3]

    if (chromosome %in% fastaChromosomes && subFeatureType == "exon" && length(intersect(featureType, featureTypes)) > 0)
    {
      chromosomes = unique(c(chromosomes, chromosome))

      startPos = as.integer(lineItems[4])
      stopPos = as.integer(lineItems[5])
      attributes = strsplit(lineItems[9], ";")[[1]]
      attributes = sub("^ ", "", attributes)
      featureID = attributes[which(grepl(paste("^", attributeType, sep=""), attributes, ignore.case=TRUE))]
      featureID = sub(paste(attributeType, " ", sep=""), "", featureID)
      featureID = gsub("\"", "", featureID)

      outRow = c(featureID, startPos, stopPos)
      write.table(matrix(outRow, nrow=1), paste(tmpDir, "/", chromosome, sep=""), sep="\t", row.names=F, col.names=F, quote=F, append=T)
    }
  }

  close(f)

  chromosomes = sort(chromosomes)
  for (chromosome in chromosomes)
  {
    message(paste("Summarizing feature values for chromosome ", chromosome, sep=""))
    ProcessGtfSubset(paste(tmpDir, "/", chromosome, sep=""), fastaStrings, chromosome, outFilePath, chromosome!=chromosomes[1])
  }
}

ProcessGtfSubset = function(tmpFilePath, fastaStrings, chromosome, outFilePath, appendToOut)
{
  data = read.table(tmpFilePath, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F)

  featureCoordsList = list()

  f = file(tmpFilePath, "r")
  stop = FALSE

  while (!stop)
  {
    line = readLines(f, n = 1)

    if(length(line) == 0)
    {
      stop = TRUE
      next
    }

    lineItems = strsplit(line, "\t")[[1]]
    featureID = lineItems[1]
    startPos = as.integer(lineItems[2])
    stopPos = as.integer(lineItems[3])
    thisRange = IRanges(startPos, stopPos)

    if (is.null(featureCoordsList[[featureID]]))
    {
      featureCoordsList[[featureID]] = thisRange
    } else {
      featureCoordsList[[featureID]] = union(featureCoordsList[[featureID]], thisRange)
    }
  }

  close(f)

  outData = NULL
  for (featureID in names(featureCoordsList))
  {
    featureData = as.matrix(as.data.frame(featureCoordsList[[featureID]]))
    featureData = cbind(featureData, NA)

    for (i in 1:nrow(featureData))
    {
      start = featureData[i,1]
      stop = featureData[i,2]

      subSequence = subseq(fastaStrings[chromosome], start, stop)
      gcCount = sum(letterFrequency(subSequence, c("C", "G")))

      featureData[i,4] = gcCount
    }

    if (nrow(featureData) == 1)
    {
      outRow = featureData[1,3:4]
    } else {
      outRow = apply(featureData[,3:4], 2, sum)
    }

    outData = rbind(outData, c(featureID, outRow))
  }

  write.table(outData, outFilePath, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=appendToOut)
}
