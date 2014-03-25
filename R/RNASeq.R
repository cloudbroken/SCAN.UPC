##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

UPC_RNASeq = function(inFilePattern, annotationFilePath=NA, outFilePath=NA, modelType="nn", convThreshold=0.01, ignoreZeroes=FALSE, verbose=TRUE)
{
  if (is.na(annotationFilePath))
  {
    annotationData = NA
  } else {
    if (!(file.exists(annotationFilePath)))
      stop(paste("No annotation file exists at ", annotationFilePath, ".", sep=""))

    print("Importing annotations...")
    annotationData = read.table(annotationFilePath, sep="\t", header=FALSE, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE)
  }

  inFilePaths = list.files(path=dirname(inFilePattern), pattern=glob2rx(basename(inFilePattern)), full.names=TRUE)

  if (length(inFilePaths) == 0)
    stop("No data files that match the pattern ", inFilePattern, " could be located.")

  outData = NULL
  rowNamesFromFirstFile = NULL

  for (inFilePath in inFilePaths)
  {
    message(paste("Processing", inFilePath))
    data = read.table(inFilePath, sep="\t", header=FALSE, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE)

    if (is.null(rowNamesFromFirstFile))
      rowNamesFromFirstFile = rownames(data)

    if (any(is.na(annotationData)))
    {
      # This is necessary to make the data consistent when no annotation data are present
      data = cbind(rownames(data), data[,1])
    } else {
      numOverlappingAnnotations = length(intersect(rownames(data), rownames(annotationData)))
      if (numOverlappingAnnotations < (nrow(data) * 0.5))
        stop("Less than half of the annotations overlap with the data in ", inFilePath, ".", sep="")

      data = merge(data, annotationData, by=0, sort=FALSE)
    }

    if (ignoreZeroes)
    {
      zeroIndices = which(data[,2] == 0)
      zeroData = NULL

      if (length(zeroIndices) > 0)
      {
        zeroData = as.matrix(data[zeroIndices,1:2])
        data = data[-zeroIndices,]
      }
    }

    naData = NULL
    naIndicator = is.na(as.numeric(data[,2]))
    if (any(naIndicator))
    {
#      if (sum(naIndicator) == 1)
#      {
        naData = data[which(naIndicator),1:2]
#      } else {
#        naData = data[which(naIndicator),1:2]
#      }

      data = data[-which(naIndicator),]
    }

    counts = as.numeric(data[,2])

    if (ncol(data) == 2)
    {
      lengths = NULL
      gc = NULL
    } else {
      if (any(is.na(data[,3])))
      {
        warning("Could not adjust for length because at least one NA value was present.")
      } else {
        if (length(unique(data[,3]))==1)
        {
          warning("Could not adjust for length because all values were the same.")
        } else {
          lengths = as.numeric(data[,3])

          if (any(is.na(lengths)))
            stop("At least one of the length values was not numeric (or NA).")
        }
      }

      if (ncol(data) == 3)
      {
        warning("Could not adjust for GC content because no GC data were present.")
      } else {
        if (any(is.na(data[,4])))
        {
          warning("Could not adjust for GC content because at least one NA value was present.")
        } else {
          if (length(unique(data[,4]))==1)
          {
            warning("Could not adjust for GC content because all values were the same.")
          }
          else {
            gc = as.numeric(data[,4])

            if (any(is.na(gc)))
              stop("At least one of the gcContent values was not numeric (or NA).")
          
            gc = gc / lengths
          }
        }
      }
    }

    upc = round(UPC_Generic(counts, lengths=lengths, gcContent=gc, modelType=modelType, convThreshold=convThreshold), 6)

    outSampleData = cbind(data[,1], upc)

    if (ignoreZeroes)
      outSampleData = rbind(outSampleData, zeroData)

    outSampleData = rbind(outSampleData, naData)

    if (is.null(outData))
    {
      outData = outSampleData
    } else {
      if (verbose)
        message(paste("Merging results for ", inFilePath, "."), sep="")

      if (length(intersect(outData[,1], outSampleData[,1])) != nrow(outData))
      {
        print(setdiff(outData[,1], outSampleData[,1]))
        message(paste(inFilePath, " has at least one feature that differs from the other file(s). Features that do not overlap across the files will default to NA."))
      }

      outData = suppressWarnings(merge(outData, outSampleData, by=1, sort=FALSE, all=TRUE))

      if (nrow(outData) == 0)
        stop("After merging the data files, no features were common across all files.")
    }
  }

  colnames(outData) = basename(c("Feature", inFilePaths))
  outData = outData[match(rowNamesFromFirstFile, outData[,1]),]

  if (!is.na(outFilePath))
  {
    if (verbose)
      message(paste("Saving results to", outFilePath))

    write.table(outData, outFilePath, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }

  featureNames = outData[,1]

  if (ncol(outData) == 2)
  {
    outData = matrix(as.numeric(outData[,2]), ncol=1)
  } else {
    outData = apply(outData[,2:ncol(outData)], 2, as.numeric)
  }

  return(list(featureNames=featureNames, data=outData))
}
