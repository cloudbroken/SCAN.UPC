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
