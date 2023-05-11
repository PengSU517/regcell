newcellmap = function (D, R, indcells = NULL, indrows = NULL, standOD = NULL,
                       showVals = NULL, rowlabels = "", columnlabels = "", mTitle = "",
                       rowtitle = "", columntitle = "", showrows = NULL, showcolumns = NULL,
                       nrowsinblock = 1, ncolumnsinblock = 1, autolabel = TRUE,
                       columnangle = 90, sizetitles = 1.1, adjustrowlabels = 1,
                       adjustcolumnlabels = 1, colContrast = 1, outlyingGrad = TRUE,
                       darkestColor = sqrt(qchisq(0.999, 1)), drawCircles = TRUE)
{
  funcSqueeze = function(Xin, n, d, ncolumnsinblock, nrowsinblock,
                         colContrast) {
    Xblock = matrix(0, nrow = n, ncol = d)
    Xblockgrad = matrix(0, nrow = n, ncol = d)
    for (i in seq_len(n)) {
      for (j in seq_len(d)) {
        Xsel = Xin[(1 + ((i - 1) * nrowsinblock)):(i *
                                                     nrowsinblock), (1 + ((j - 1) * ncolumnsinblock)):(j *
                                                                                                         ncolumnsinblock)]
        seltable = tabulate(Xsel, nbins = 4)
        if (sum(seltable) > 0) {
          indmax = which(seltable == max(seltable))[1]
          cntmax = seltable[indmax]
          gradmax = (cntmax/(ncolumnsinblock * nrowsinblock))^(1/colContrast)
        }
        else {
          indmax = 0
          gradmax = 1
        }
        Xblock[i, j] = indmax
        Xblockgrad[i, j] = gradmax
      }
    }
    return(list(X = Xblock, Xgrad = Xblockgrad))
  }
  variable <- rownr <- rescaleoffset <- x <- y <- NULL
  type = "cell"
  n = nrow(R)
  d = ncol(R)
  blockMap = FALSE
  if (ncolumnsinblock > 1 | nrowsinblock > 1) {
    blockMap = TRUE
    if (ncolumnsinblock > d)
      stop("Input argument ncolumnsinblock cannot be larger than d")
    if (nrowsinblock > n)
      stop("Input argument nrowsinblock cannot be larger than n")
    if (!is.null(showVals))
      warning("The option showVals=D or showVals=R cannot\n
              be combined with ncolumnsinblock or nrowsinblock greater than 1,\n
              so showVals is set to NULL here.")
    showVals = NULL
  }
  if (!blockMap) {
    if (!all(dim(R) == dim(D)))
      stop("Dimensions of D and R must match")
  }
  if (!(blockMap & autolabel == FALSE)) {
    if (length(columnlabels) > 0 & length(columnlabels) !=
        d) {
      stop(paste("Number of columnlabels does not match d = ",
                 d, sep = ""))
    }
    if (length(rowlabels) > 0 & length(rowlabels) != n) {
      stop(paste("Number of rowlabels does not match n = ",
                 n, sep = ""))
    }
  }
  if (!is.null(showVals)) {
    if (!showVals %in% c("D", "R")) {
      stop(paste("Invalid \"showVals\" argument. Should be one of: NULL, \"D\", \"R\""))
    }
  }
  if (is.null(indcells))
    indcells = which(abs(R) > sqrt(qchisq(0.99, 1)))
  if (!(is.null(showcolumns) & is.null(showrows))) {
    if (is.null(showcolumns)) {
      showcolumns = seq_len(d)
    }
    else {
      if (!(all(showcolumns %in% seq_len(d))))
        stop(" showcolumns goes out of bounds")
    }
    if (is.null(showrows)) {
      showrows = seq_len(n)
    }
    else {
      if (!(all(showrows %in% seq_len(n))))
        stop(" showrows goes out of bounds")
    }
    tempMat = matrix(0, n, d)
    tempMat[indcells] = 1
    tempMat = tempMat[showrows, showcolumns]
    indcells = which(tempMat == 1)
    tempVec = rep(0, n)
    tempVec[indrows] = 1
    tempVec = tempVec[showrows]
    indrows = which(tempVec == 1)
    rm(tempMat, tempVec)
    if (!blockMap)
      D = D[showrows, showcolumns]
    R = R[showrows, showcolumns]
    if (!(blockMap & autolabel == FALSE))
      columnlabels = columnlabels[showcolumns]
    if (!(blockMap & autolabel == FALSE))
      rowlabels = rowlabels[showrows]
    n = nrow(R)
    d = ncol(R)
    if (!is.null(standOD))
      standOD = standOD[showrows]
  }
  if (type == "residual")
    outlyingGrad = 1
  X = matrix(0, n, d)
  Xrow = matrix(0, n, 1)
  Xrow[indrows, 1] = 3
  if (type == "cell" | blockMap) {
    pcells = indcells[indcells %in% which(R >= 0)]
    ncells = indcells[indcells %in% which(R < 0)]
  }
  else {
    pcells = which(R >= 0)
    ncells = which(R < 0)
  }
  X[ncells] = 1
  X[pcells] = 2
  X[is.na(R)] = 4
  if (blockMap) {
    n = floor(n/nrowsinblock)
    d = floor(d/ncolumnsinblock)
    result = funcSqueeze(X, n, d, ncolumnsinblock, nrowsinblock,
                         colContrast)
    X = result$X
    Xgrad = result$Xgrad
    result = funcSqueeze(Xrow, n, 1, 1, nrowsinblock, colContrast)
    Xrowgrad = result$Xgrad
    Xrowgrad[result$X == 0] = 0
    if (autolabel == TRUE) {
      if (ncolumnsinblock > 1 & length(columnlabels) >
          0) {
        labx = columnlabels
        columnlabels = rep(0, d)
        for (ind in seq_len(d)) {
          columnlabels[ind] = paste(labx[(1 + ((ind -1) * ncolumnsinblock))], "-", labx[(ind *ncolumnsinblock)], sep = "")
        }
      }
      if (nrowsinblock > 1 & length(rowlabels) > 0) {
        laby = rowlabels
        rowlabels = rep(0, n)
        for (ind in seq_len(n)) {
          rowlabels[ind] = paste(laby[(1 + ((ind - 1) *
                                              nrowsinblock))], "-", laby[(ind * nrowsinblock)])
        }
      }
    }
    else {
      if (length(columnlabels) > 0 & length(columnlabels) !=
          d) {
        stop(paste(" autolabel=FALSE and number of columnlabels is ",
                   length(columnlabels), " but should be ", d,
                   sep = ""))
      }
      if (length(rowlabels) > 0 & length(rowlabels) !=
          n) {
        stop(paste(" autolabel=FALSE and number of rowlabels is ",
                   length(rowlabels), " but should be ", n, sep = ""))
      }
    }
    Xdf = data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) = c("rownr", seq(1, d, 1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX = melt(Xdf, id.var = "rownr", value.name = "CatNr")
    Xgraddf = data.frame(cbind(seq(1, n, 1), Xgrad))
    colnames(Xgraddf) = c("rownr", seq(1, d, 1))
    rownames(Xgraddf) = NULL
    Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n, 1,
                                                     -1)))
    mXgrad = melt(Xgraddf, id.var = "rownr", value.name = "grad")
    mX$grad = mXgrad$grad
    mX$rescaleoffset = mXgrad$grad + 10 * mX$CatNr
    mXrow = data.frame(rownr = seq_len(n), rescaleoffset = Xrowgrad +
                         10 * 3)
    scalerange = c(0, 1)
    gradientends = scalerange + rep(c(0, 10, 20, 30, 40),
                                    each = 2)
    if (type == "cell")
      colorends = c("white", "white", "white", "blue",
                    "white", "red", "white", "black", "white",
                    "white")
    if (type == "residual")
      colorends = c("white", "white", "white", "blue",
                    "white", "red", "white", "black", "white", "white")
  }
  else {
    Ddf = data.frame(cbind(seq(1, n, 1), D))
    colnames(Ddf) = c("rownr", seq(1, d, 1))
    rownames(Ddf) = NULL
    Ddf$rownr = with(Ddf, reorder(rownr, seq(n, 1, -1)))
    mD = melt(Ddf, id.var = "rownr")
    Rdf = data.frame(cbind(seq(1, n, 1), R))
    colnames(Rdf) = c("rownr", seq(1, d, 1))
    rownames(Rdf) = NULL
    Rdf$rownr = with(Rdf, reorder(rownr, seq(n, 1, -1)))
    mR = melt(Rdf, id.var = "rownr")
    Xdf = data.frame(cbind(seq(1, n, 1), X))
    colnames(Xdf) = c("rownr", seq(1, d, 1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n, 1, -1)))
    mX = melt(Xdf, id.var = "rownr", value.name = "CatNr")
    if (!is.null(showVals)) {
      if (showVals == "D")
        mX$data = mD$value
      if (showVals == "R")
        mX$data = mR$value
    }
    if (!outlyingGrad) {
      mX$rescaleoffset = 10 * mX$CatNr
      scalerange = c(0, 1)
      gradientends = scalerange + rep(c(0, 10, 20, 30,
                                        40), each = 2)
      gradientends
      colorends = c("white", "white", "blue", "blue",
                    "red", "red", "white", "black", "white", "white")
    }
    else {
      Xgrad = matrix(NA, n, d)
      if (type == "cell") {
        Xgrad[indcells] = abs(R[indcells])
        limL = sqrt(qchisq(0.9, 1))
      }
      else {
        Xgrad = abs(R)
        limL = 0
      }
      limH = darkestColor
      Xgrad[Xgrad > limH] = limH
      Xgrad = ((Xgrad - limL)/(limH - limL))^colContrast
      Xgrad[is.na(Xgrad)] = 0
      Xgraddf = data.frame(cbind(seq(1, n, 1), Xgrad))
      colnames(Xgraddf) = c("rownr", seq(1, d, 1))
      rownames(Xgraddf) = NULL
      Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n,
                                                       1, -1)))
      mXgrad = melt(Xgraddf, id.var = "rownr", value.name = "grad")
      mX$grad = mXgrad$grad
      mX$rescaleoffset = mXgrad$grad + 10 * mX$CatNr
      scalerange = c(0, 1)
      gradientends = scalerange + rep(c(0, 10, 20, 30,
                                        40), each = 2)
      if (type == "cell")
        colorends = c("white", "white", "white", "blue",
                      "white", "red", "white", "black", "white",
                      "white")
      if (type == "residual")
        colorends = c("white", "white", "white", "blue",
                      "white", "red", "white", "black", "white",
                      "white")
    }
    tempVec = rep(0, n)
    tempVec[indrows] = 1
    mXrow = data.frame(rownr = seq_len(n), rescaleoffset = 40 -
                         (10 * tempVec))
    rm(tempVec)
    if (is.null(standOD)) {
      mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] +
        1
    }
    else {
      limL = 1
      limH = 3
      standOD[standOD > limH] = limH
      standOD = ((standOD - limL)/(limH - limL))^colContrast
      mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] +
        standOD[indrows]
    }
  }
  rowlabels = rev(rowlabels)
  base_size = 10
  columnlabels = c(columnlabels, "", "")
  circleFun = function(centerx, centery, r, npoints) {
    tt = seq(0, 2 * pi, length.out = npoints)
    xx = centerx + r * cos(tt)
    yy = centery + r * sin(tt)
    return(c(xx, yy))
  }
  if (drawCircles) {
    centerx = d + 1
    centery = n:1
    radius = 0.4
    npoints = 100
    circlePoints = mapply(circleFun, centerx, centery, radius,
                          npoints)
    positions = data.frame(rownr = rep(seq_len(n), each = npoints),
                           x = c(circlePoints[seq_len(npoints), ]), y = c(circlePoints[(npoints +
                                                                                          1):(2 * npoints), ]))
    datapoly = merge(mXrow, positions, by = c("rownr"))
  }
  ggp = ggplot(data = mX, aes(variable, rownr)) + {
    if (blockMap)
      geom_tile(aes(fill = scales::rescale(rescaleoffset, from = range(gradientends))),
                color = "gray")
  } + {
    if (!blockMap & outlyingGrad)
      geom_tile(aes(fill = scales::rescale(rescaleoffset, from = range(gradientends))),
                color = "gray")
  } + {
    if (!blockMap & !outlyingGrad)
      geom_tile(aes(fill = scales::rescale(rescaleoffset, from = range(gradientends))),
                colour = "gray")
  } + {
    if (drawCircles)
      geom_polygon(data = datapoly,
                   aes(x = x, y = y, fill = scales::rescale(rescaleoffset, from = range(gradientends)), group = rownr),
                   colour = "white")
  } + scale_fill_gradientn(colours = colorends, values = scales::rescale(gradientends),
                           rescaler = function(x, ...) x, oob = scales::squish) +
    ggtitle(mTitle) +
    coord_fixed() +
    theme_classic(base_size = base_size *1) +
    labs(x = columntitle, y = rowtitle) +
    scale_x_discrete(expand = c(0, 0), limits = as.factor(seq(1, d + 2, 1)), labels = columnlabels) +
    scale_y_discrete(expand = c(0, 0), labels = rowlabels) +
    theme(legend.position = "none", axis.ticks = element_blank(),
          plot.title = element_text(size = base_size * 2, hjust = 0.5,  vjust = 1, face = "bold"),
          axis.text.x = element_text(size = base_size * 1.8,
                                     angle = columnangle,
                                     hjust = adjustcolumnlabels,
                                     vjust = 0.5, colour = "black"),
          axis.text.y = element_text(size = base_size *
                                                                                                                                                       1.8, angle = 0, hjust = adjustrowlabels, colour = "black"),
          axis.title.x = element_text(colour = "black", size = base_size *sizetitles, vjust = 1),
          axis.title.y = element_text(colour = "black",
                                      size = base_size * sizetitles,
                                      vjust = 0),
          axis.line.x = element_blank(),
          panel.border = element_blank()) +
    annotate(geom = "segment", x = 0.5, xend = d + 0.5, y = 0.5, yend = 0.5) +
    annotate(geom = "segment",
             x = 0.5, xend = d + 0.5, y = n + 0.5, yend = n + 0.5) +
    annotate(geom = "segment", x = d + 0.5, xend = d + 0.5,
             y = 0.5, yend = n + 0.5)
  if (!is.null(showVals)) {
    txtcol = mX$CatNr
    txtcol[txtcol == 0] = "black"
    txtcol[txtcol == 1] = "white"
    txtcol[txtcol == 2] = "white"
    if (type == "residual") {
      txtcol[] = "black"
      txtcol[mXgrad$grad > 0.5] = "white"
    }
    txtcol[txtcol == 4] = "black"
    ggp = ggp + geom_text(aes(label = ifelse(is.na(data),
                                             sprintf("%1.0f", data), round(data, 1))), size = base_size *
                            0.5, colour = txtcol, na.rm = TRUE)
  }
  return(ggp)
}

