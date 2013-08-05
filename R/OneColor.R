##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

SCAN = function(celFilePattern, outFilePath=NA, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, verbose=TRUE)
{
  return(processCelFiles(celFilePattern=celFilePattern, outFilePath=outFilePath, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=FALSE, exonArrayTarget=exonArrayTarget, verbose=verbose))
}

UPC = function(celFilePattern, outFilePath=NA, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, verbose=TRUE)
{
  return(processCelFiles(celFilePattern=celFilePattern, outFilePath=outFilePath, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=TRUE, exonArrayTarget=exonArrayTarget, verbose=verbose))
}

processCelFiles = function(celFilePattern, outFilePath=NA, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, UPC=FALSE, exonArrayTarget=NA, verbose=TRUE)
{
  if (convThreshold >= 1)
    stop("The convThreshold value must be lower than 1.0.")

  if (shouldDownloadFromGEO(celFilePattern))
    celFilePattern = downloadFromGEO(celFilePattern)

  celFilePaths = list.files(path=dirname(celFilePattern), pattern=glob2rx(basename(celFilePattern)), full.names=TRUE)

  if (length(celFilePaths) == 0)
    stop("No CEL files that match the pattern ", celFilePattern, " could be located.")

  if (!is.na(outFilePath))
    createOutputDir(dirPath=dirname(outFilePath), verbose=verbose)

  if (!is.na(probeLevelOutDirPath))
    createOutputDir(dirPath=probeLevelOutDirPath, verbose=verbose)

  summarized = NULL

  for (celFilePath in celFilePaths)
  {
    probeLevelOutFilePath = NA
    if (!is.na(probeLevelOutDirPath))
      probeLevelOutFilePath = paste(probeLevelOutDirPath, "/", basename(celFilePath), ".txt", sep="")

    celSummarized = processCelFile(celFilePath=celFilePath, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, UPC=UPC, convThreshold=convThreshold, probeLevelOutFilePath=probeLevelOutFilePath, exonArrayTarget=exonArrayTarget, verbose=verbose)

    if (is.null(celSummarized))
      next

    if (is.null(summarized))
    {
      summarized = celSummarized
      colnames(summarized) = basename(celFilePath)
    } else {
      previousColNames = colnames(summarized)
      summarized = merge(summarized, celSummarized, sort=FALSE, by=0)
      rownames(summarized) = summarized[,1]
      summarized = summarized[,-1]
      colnames(summarized) = c(previousColNames, basename(celFilePath))
    }
  }

  if (is.null(summarized))
    return(ExpressionSet())

  if (!is.na(outFilePath))
  {
    write.table(round(summarized, 8), outFilePath, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
    message("Saved output to ", outFilePath, sep="")
  }

  return(ExpressionSet(as.matrix(summarized)))
}

createOutputDir = function(dirPath, verbose=TRUE)
{
  if (!file.exists(dirPath))
  {
    if (!dir.create(dirPath))
      stop("An output file cannot be saved because a directory does not exist at ", dirPath, " and could not be created.")
  }
}

processCelFile = function(celFilePath, annotationPackageName, probeSummaryPackage, UPC, probeLevelOutFilePath, exonArrayTarget, nbins=25, binsize=5000, convThreshold=0.01, verbose=TRUE)
{
  if (is.na(annotationPackageName))
  {
    affyExpressionFS <- read.celfiles(celFilePath)
  } else {
    affyExpressionFS <- read.celfiles(celFilePath, pkgname=annotationPackageName)
  }

  if (is.na(exonArrayTarget))
  {
    pmSeq = pmSequence(affyExpressionFS)
  } else {
    pmSeq = pmSequence(affyExpressionFS, target=exonArrayTarget)
  }

  shouldUseProbes = width(pmSeq)==25

  if (!is.na(probeLevelOutFilePath) && file.exists(probeLevelOutFilePath))
  {
    data = read.table(probeLevelOutFilePath, sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=1)
    y_norm = data[,1]
    gam = data[,2]
    xyCoord = rownames(data)
  } else
  {
    numSequences = length(pmSeq)

    chunkSize = 100000
    chunkStartIndex = 1
    chunkEndIndex = chunkStartIndex + chunkSize - 1
    mxStartIndex = 1

    mx = matrix(nrow=sum(shouldUseProbes), ncol=80)

    while (chunkStartIndex <= numSequences)
    {
      if (chunkEndIndex > numSequences)
        chunkEndIndex = numSequences

      mxChunk = buildDesignMatrix(pmSeqs=pmSeq[chunkStartIndex:chunkEndIndex], verbose=verbose)
      mx[mxStartIndex:(mxStartIndex + nrow(mxChunk) - 1),] = mxChunk

      chunkStartIndex = chunkStartIndex + chunkSize
      chunkEndIndex = chunkEndIndex + chunkSize
      mxStartIndex = mxStartIndex + nrow(mxChunk)
    }

    if (is.na(exonArrayTarget))
    {
      pint = pm(affyExpressionFS)[which(shouldUseProbes),]
    } else {
      pint = pm(affyExpressionFS, target=exonArrayTarget)[which(shouldUseProbes),]
    }

    if ((sum(pint==0) / length(pint)) > 0.01)
    {
      message(paste(celFilePath, " has a disproportionate number of zero values, so it cannot be processed.", sep=""))
      return(NULL)
    }

    my = log2(pint)
    nGroups = length(my) / binsize
    samplingProbeIndices = getSampleIndices(total=length(pint), verbose=verbose)

    mixResult = EM_vMix(y=my[samplingProbeIndices], X=mx[samplingProbeIndices,], nbins=nbins, convThreshold=convThreshold, verbose=verbose, demo=length(grep("Vignette_Example", basename(celFilePath))) > 0)

    m1 = mx %*% mixResult$b1
    m2 = mx %*% mixResult$b2

    index = order(m1)
    y_norm = rep(0, length(my))
    for (i in 1:nGroups)
    {
      tmp = index[(binsize * i):min(binsize * i + binsize, length(my))]
      tmpSd = sig(y=my[tmp], m=m1[tmp], verbose=verbose)
      y_norm[tmp] = ((my[tmp] - m1[tmp]) / tmpSd)
    }

    bin = assign_bin(y=m1, nbins=nbins, verbose=verbose)
    gam = vresp(y=my, X=mx, bin=bin, p=mixResult$p, b1=mixResult$b1, s1=mixResult$s1, b2=mixResult$b2, s2=mixResult$s2, verbose=verbose)[,2]

    y_norm = round(y_norm, 8)
    gam = round(gam, 8)

    if (!is.na(probeLevelOutFilePath))
    {
      message("Outputting probe-level values to ", probeLevelOutFilePath)
      normOutput = cbind(xyCoord, y_norm, gam, round(my, 8), round(m1, 8), round(m2, 8))
      write.table(normOutput, probeLevelOutFilePath, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    }
  }

  if (!any(is.na(probeSummaryPackage)))
  {


    probeSetNames = get("Probe.Set.Name", probeSummaryPackage)
    probeSummaryXCoord = get("x", probeSummaryPackage)
    probeSummaryYCoord = get("y", probeSummaryPackage)
    probeSummaryXyCoord = paste(probeSummaryXCoord, probeSummaryYCoord, sep="_")

    probeSummaryData = matrix(probeSetNames, ncol=1)
    rownames(probeSummaryData) = probeSummaryXyCoord

    xyCoord = paste(getX(affyExpressionFS, type="pm")[which(shouldUseProbes)], getY(affyExpressionFS, type="pm")[which(shouldUseProbes)], sep="_")
    dataProbeIndices = which(xyCoord %in% probeSummaryXyCoord)
    keepXyCoord = xyCoord[dataProbeIndices]

    probeNames = as.character(probeSummaryData[keepXyCoord,1])
  } else
  {
    dataProbeIndices = 1:length(y_norm)

    if (is.na(exonArrayTarget))
    {
      probeNames = probeNames(affyExpressionFS)[which(shouldUseProbes)]
    } else {
      probeNames = probeNames(affyExpressionFS, target=exonArrayTarget)[which(shouldUseProbes)]
    }
  }

  if (UPC)
    return(as.matrix(tapply(gam[dataProbeIndices], probeNames, FUN=mean, trim=0.1)))

  return(as.matrix(tapply(y_norm[dataProbeIndices], probeNames, FUN=mean, trim=0.1)))
}

getSampleIndices = function(total, verbose=TRUE)
{
  interval = floor(total / 50000)
  if (interval <= 1)
    interval = 1

  seq(1, total, interval)
}

assign_bin = function(y, nbins, verbose=TRUE)
{
  quans = sort(y)[floor(length(y) * 1:nbins / nbins)]
  bins = sapply(y, function(x) { sum(x>quans) }) + 1

  if (length(table(bins)) != nbins)
  {
    if (verbose)
      message("The values were not separated into enough bins, so a tiny amount of noise will be added to make this possible.")

    set.seed(1)
    noise = rnorm(length(y)) / 10000000
    bins = assign_bin(y + noise, nbins, verbose)
  }

  bins
}

buildDesignMatrix = function(pmSeqs, verbose=TRUE)
{
  keepIndices = which(width(pmSeqs)==25)
  seqs = as.character(pmSeqs[keepIndices])

  mx = sequenceDesignMatrix(seqs)
  numA = apply(mx[,1:25], 1, sum)
  numC = apply(mx[,26:50], 1, sum)
  numG = apply(mx[,51:75], 1, sum)
  numT = 25 - (numA + numC + numG)

  mx = cbind(numT, mx, as.integer(numA^2), as.integer(numC^2), as.integer(numG^2), as.integer(numT^2))
  mx = apply(mx, 2, as.integer)
}

EM_vMix = function(y, X, nbins, convThreshold=.01, verbose=TRUE, maxIt=100, demo=FALSE)
{
  if (verbose)
    message("Starting EM")

  quan = sort(y)[floor(0.5 * length(y)) - 1]
  gam <<- cbind(as.integer(y <= quan), as.integer(y > quan))

  p = apply(gam, 2, mean)

  b1 = mybeta(y=y, X=X, gam=gam[,1], verbose=verbose)
  b2 = mybeta(y=y, X=X, gam=gam[,2], verbose=verbose)
  bin = assign_bin(y=y, nbins=nbins, verbose=verbose)
  s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
  s2 = vsig(y=y, X=X, b=b2, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)

  theta_old=c(p, b1, s1, b2, s2)

  it = 0
  conv = 1

  while (conv > convThreshold & it < maxIt)
  {
    # Expectation Step:
    gam = vresp(y=y, X=X, bin=bin, p=p, b1=b1, s1=s1, b2=b2, s2=s2, verbose=verbose)

    #M-Step
    p = apply(gam, 2, mean)
    b1 = vbeta(y=y, X=X, bin=bin, gam=gam[,1], s2=s1, prof=TRUE, verbose=verbose)
    bin = assign_bin(y=(X %*% b1), nbins=nbins, verbose=verbose)
    b2 = vbeta(y=y, X=X, bin=bin, gam=gam[,2], s2=s2, prof=FALSE, verbose=verbose)
    s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
    s2 = vsig(y=y, X=X, b=b2, gam=gam[,2], bin=bin, nbins=nbins, verbose=verbose)

    theta = c(p, b1, s1, b2, s2)
    conv = max(abs(theta - theta_old) / theta_old)
    theta_old = theta
    it = it + 1

    if (verbose)
      message("Attempting to converge...iteration ", it, ", c = ", round(conv, 6))

    if (demo & (it == 3))
      break
  }

  if (verbose)
  {
    if (demo)
    {
      message("Convergence process halted early while in demo mode.")
    } else
    {
      if (it == maxIt)
      {
        message("Reached convergence limit...", it, " iterations. Proportion of background probes: ", round(p[1], 6))
      } else {
        message("Converged in ", it, " iterations. Proportion of background probes: ", round(p[1], 6))
      }
    }
  }

  list(p=p, b1=b1, b2=b2, s1=s1, s2=s2, bin=bin)
}

mybeta = function(y, X, gam, verbose=TRUE)
{
  sqgam = sqrt(gam)
  Xw = sqgam * X
  yw = sqgam * y

  z = t(Xw) %*% Xw
  a = solve(z)

  b = a %*% t(Xw)
  as.numeric(b %*% yw)
}

vbeta = function(y, X, bin, gam, s2, prof, verbose=TRUE)
{
  vars = sqrt(s2[bin])
  sqgam = sqrt(gam)
  vars_sqgam = vars * sqgam

  Xw = 1 / vars * sqgam * X
  yw = 1 / vars * sqgam * y

  tXw = t(Xw)
  tXwXw = tXw %*% Xw
  stXwXw = solve(tXwXw)
  stXwXwtXw = stXwXw %*% tXw
  result = stXwXwtXw %*% yw

  result
}

vresp = function(y, X, bin, p, b1, s1, b2, s2, verbose=TRUE)
{
  vars0 = s1[bin]
  L0 = dn(y=y, m=(X %*% b1), s2=vars0, verbose=verbose)
  vars1 = s2[bin]
  L1 = dn(y=y, m=(X %*% b2), s2=vars1, verbose=verbose)

  gam1 = p[1] * L0 / (p[1] * L0 + p[2] * L1)
  gam2 = 1 - gam1
  cbind(gam1, gam2)
}

dn = function(y, m, s2, verbose=TRUE)
{
  1 / (sqrt(2 * pi * s2)) * exp(-1 / (2 * s2) * (y - m)^2)
}

vsig = function(y, X, b, gam, bin, nbins, verbose=TRUE)
{
  s2 = NULL

  for (i in 1:nbins)
  {
    ystar = y[bin==i]
    Xstar = X[bin==i,]
    gamstar = gam[bin==i] + .01
    resid = as.numeric(ystar - Xstar %*% b)

    s2 = c(s2, ((resid * gamstar) %*% resid) / sum(gamstar))
  }

  s2
}

sig = function(y, m, verbose=TRUE)
{
  resid = y - m
  sqrt((resid %*% resid) / length(y))
}
