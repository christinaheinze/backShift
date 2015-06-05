
plotGraph <- function(A, main="",labels=NULL,layoutfunction=layout.circle,...){
    if(is(A, "dgCMatrix")) A <- as.matrix(A)
    if(!is.matrix(A)) stop("A needs to be a matrix")
    if(nrow(A)!=ncol(A)) stop("A needs to have as many rows as columns")
    if(is.null(labels)) labels <- if( !is.null(cc <- colnames(A)))
        cc else as.character(1:ncol(A))
    G <- graph.adjacency(A,mode="directed",weighted="a")
    if(is.null(layoutfunction)) layoutfunction <-  layout.circle
    layout <- layoutfunction(G)

    optionals <- list(...)
    if(is.null(optionals$vertex.label.cex)) optionals$vertex.label.cex <- 1.5
    if(is.null(optionals$vertex.label.color)) optionals$vertex.label.color <-rgb(0.8,0.1,0.1,0.7)
    if(is.null(optionals$vertex.color)) optionals$vertex.color <-"white"
    if(is.null(optionals$vertex.frame.color)) optionals$vertex.frame.color <-rgb(0.8,0.1,0.1,0.5)
    if(is.null(optionals$edge.color)) optionals$edge.color <-rgb(0.1,0.1,0.1,0.5)
    if(is.null(optionals$edge.arrow.size)) optionals$edge.arrow.size <-0.7
    if(is.null(optionals$edge.arrow.width)) optionals$edge.arrow.width <-2
    if(is.null(optionals$vertex.size)) optionals$vertex.size <-30
    if(is.null(optionals$vertex.label.dist)) optionals$vertex.label.dist <-0
    if(is.null(optionals$vertex.label.degree)) optionals$vertex.label.degree <-  -pi/2

    plot(G, layout=layout, main=main, vertex.label=labels,vertex.shape="circle",vertex.label.cex=optionals$vertex.label.cex, vertex.label.color=optionals$vertex.label.color,vertex.color=optionals$vertex.color ,vertex.frame.color=optionals$vertex.frame.color , edge.color=optionals$edge.color, edge.arrow.size=optionals$edge.arrow.size ,edge.arrow.width=optionals$edge.arrow.width ,vertex.size=optionals$vertex.size,vertex.label.dist=optionals$vertex.label.dist,vertex.label.degree=optionals$vertex.label.degree)
    return(layout)
}

#' Plots graph from adjacency matrix.
#'
#' @param pointEst
#' @param labels
#' @param thres.point
#' @param adjacencyEst
#' @param plotStabSelec
#' @param thres.stab
#' @param main
#' @param ...
plotGraphEdgeAttr <- function(pointEst, labels, thres.point, adjacencyEst, 
                              plotStabSelec, thres.stab, main, ...){
  
  # stability selection output
  df.As <- melt(as.matrix(adjacencyEst))
  colnames(df.As) <- c("from", "to", "edge")
  
  # point estimate
  colnames(pointEst) <- rownames(pointEst) <- colnames(adjacencyEst) <- rownames(adjacencyEst) <- labels
  df.Ap <- melt(as.matrix(pointEst))
  colnames(df.Ap) <- c("from", "to", "edge")
  
  is.empty <- sum(adjacencyEst) == 0 || sum(pointEst) == 0
  if(plotStabSelec & !is.empty){
      combined <- df.As
      combined$estimate <- df.Ap$edge
      combined.wo.zeros <- combined[-which(combined$edge == 0),]
      
      if(nrow(combined.wo.zeros) == 0){
        is.empty <- TRUE
      }else{
        # draw graph
        bio.network <-graph.data.frame(combined.wo.zeros, directed=TRUE, vertices = colnames(adjacencyEst))
        layoutfunction <-  layout.circle
        layout <- layoutfunction(bio.network)
        
        # using the absolute value of the point estimate for the transparency of the chosen color
        color.attribute <- abs(E(bio.network)$estimate)
        alphas.to.use <- convert.given.min(color.attribute, 0, max(color.attribute), range.min = 0.2, range.max = 0.8) 
        mycolors <- vary.alpha(alphas.to.use, "blue")
        E(bio.network)$color <- mycolors
        
        # using the number of times selected in the stability selection as the width
        width.attribute <- abs(E(bio.network)$edge)
        E(bio.network)$width <- convert.given.min(width.attribute, thres.stab*100, 100, range.min = 3, range.max = 8) 
      }
      
    }else if(!is.empty){
      combined <- df.Ap
      combined.wo.zeros <- combined[-which(abs(combined$edge) <= thres.point),]
      
      if(nrow(combined.wo.zeros) == 0){
        is.empty <- TRUE
        pointEst[abs(pointEst) <= thres.point] <- 0
      }else{
        # draw graph
        bio.network <-graph.data.frame(combined.wo.zeros, directed=TRUE, vertices = colnames(pointEst))
        layoutfunction <-  layout.circle
        layout <- layoutfunction(bio.network)
        
        # using the absolute value of the point estimate for the transparency of the chosen color
        color.attribute <- abs(E(bio.network)$edge)
        rmin <- 0.2
        rmax <- 0.5
        alphas.to.use <- convert.given.min(color.attribute, 0, max(color.attribute), range.min = rmin, range.max = rmax) 
        if(any(is.na(alphas.to.use))) alphas.to.use <- rep(rmax, length(color.attribute))
        mycolors <- vary.alpha(alphas.to.use, "blue")
        E(bio.network)$color <- mycolors
        
        # using the number of times selected in the stability selection as the width
        width.attribute <- abs(E(bio.network)$edge)
        rmin.w <- 5
        rmax.w <- 5
        widths.to.use <- convert(width.attribute, range.min = 5, range.max = 5) 
        if(any(is.na(widths.to.use))) widths.to.use <- rep(rmax.w, length(width.attribute))
        E(bio.network)$width <- widths.to.use
      }
    }
  
  optionals <- list(...)
  if(is.null(optionals$vertex.label.cex)) optionals$vertex.label.cex <-  2.5 
  if(is.null(optionals$vertex.label.color)) optionals$vertex.label.color <-"black" 
  if(is.null(optionals$vertex.color)) optionals$vertex.color <-"white"
  if(is.null(optionals$vertex.frame.color)) optionals$vertex.frame.color <-"black" 
  if(is.null(optionals$edge.arrow.size)) optionals$edge.arrow.size <-1
  if(is.null(optionals$edge.arrow.width)) optionals$edge.arrow.width <- 1.5
  if(is.null(optionals$vertex.size)) optionals$vertex.size <-30
  if(is.null(optionals$vertex.label.dist)) optionals$vertex.label.dist <-0
  if(is.null(optionals$vertex.label.degree)) optionals$vertex.label.degree <-  -pi/2
  
  if(!is.empty){
    plot(bio.network, layout=layout, main = main,
         vertex.shape="circle",
         vertex.label.cex=optionals$vertex.label.cex, 
         vertex.label.color=optionals$vertex.label.color,
         vertex.color=optionals$vertex.color,
         vertex.frame.color=optionals$vertex.frame.color,
         vertex.size=optionals$vertex.size,
         vertex.label.dist=optionals$vertex.label.dist,
         vertex.label.degree=optionals$vertex.label.degree,
         edge.arrow.size = optionals$edge.arrow.size, 
         edge.arrow.width = optionals$edge.arrow.width,
         edge.curved=autocurve.edges2(bio.network, start = 0.25))
  }else if(plotStabSelec){
    plotGraph(adjacencyEst,main=main,
              vertex.shape="circle",
              vertex.label.cex=optionals$vertex.label.cex, 
              vertex.label.color=optionals$vertex.label.color,
              vertex.color=optionals$vertex.color,
              vertex.frame.color=optionals$vertex.frame.color,
              vertex.size=optionals$vertex.size,
              vertex.label.dist=optionals$vertex.label.dist,
              vertex.label.degree=optionals$vertex.label.degree,
              edge.arrow.size = optionals$edge.arrow.size, 
              edge.arrow.width = optionals$edge.arrow.width)
  }else{
    plotGraph(pointEst,main=main,
              vertex.shape="circle",
              vertex.label.cex=optionals$vertex.label.cex, 
              vertex.label.color=optionals$vertex.label.color,
              vertex.color=optionals$vertex.color,
              vertex.frame.color=optionals$vertex.frame.color,
              vertex.size=optionals$vertex.size,
              vertex.label.dist=optionals$vertex.label.dist,
              vertex.label.degree=optionals$vertex.label.degree,
              edge.arrow.size = optionals$edge.arrow.size, 
              edge.arrow.width = optionals$edge.arrow.width)
  }
}

