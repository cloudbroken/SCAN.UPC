##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

SCAN_TwoColor = function(inFilePattern, outFilePath=NA, verbose=TRUE)
{
  ProcessTwoColor(inFilePattern, outFilePath=outFilePath, verbose=verbose)
}

UPC_TwoColor = function(inFilePattern, outFilePath=NA, modelType="nn", convThreshold=0.01, verbose=TRUE)
{
  ProcessTwoColor(inFilePattern, outFilePath=outFilePath, upcModelType=modelType, upcConv=convThreshold, verbose=verbose)
}

ProcessTwoColor = function(inFilePattern, outFilePath=NA, upcModelType=NA, upcConv=0.01, verbose=TRUE)
{
  inFilePaths = list.files(path=dirname(inFilePattern), pattern=glob2rx(basename(inFilePattern)), full.names=TRUE)
  outData = NULL

  for (inFilePath in inFilePaths)
  {
    message(paste("Normalizing", inFilePath))
    normData = normalizeArray(filename = inFilePath)

    if (!is.na(upcModelType))
    {
      if (verbose)
        message(paste("Performing UPC transformation (", upcModelType, ") for first channel of ", inFilePath, sep=""))

      normData$norm[,1] = round(UPC_Transform(2^normData$norm[,1], modelType=upcModelType, conv=upcConv), 6)

      if (verbose)
        message(paste("Performing UPC transformation (", upcModelType, ") for second channel of ", inFilePath, sep=""))

      normData$norm[,2] = round(UPC_Transform(2^normData$norm[,2], modelType=upcModelType, conv=upcConv), 6)
    }

    outSampleData = cbind(normData$probeID[,1], normData$norm)
    colnames(outSampleData) = c("ProbeID", paste(basename(inFilePath), c("Channel1", "Channel2"), sep="_"))

    if (is.null(outData))
    {
      outData = outSampleData
    } else {
      if (verbose)
        message(paste("Merging results for ", inFilePath, "."), sep="")

      if (all(outData[,1] == outSampleData[,1]))
      {
        outData = cbind(outData, outSampleData[,2:3])
      }
      else {
        warning(paste(inFilePath, " has different features than the previous file(s), so an attempt will be made to match by feature name. However, because there are duplicate feature names on these arrays, proceed with caution."), immediate. = TRUE)
        outData = merge(outData, outSampleData, by=1, sort=FALSE)
      }
    }
  }

  if (!is.na(outFilePath))
  {
    if (verbose)
      message(paste("Saving results to", outFilePath))

    write.table(outData, outFilePath, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }

  return(list(featureNames=outData[,1], data=apply(outData[,2:ncol(outData)], 2, as.numeric)))
}
