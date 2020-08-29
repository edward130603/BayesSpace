#' Plot cluster labels in spatial context
#'
#' After running \code{spatialCluster()} or \code{spatialEnhance()}, we can
#' visualize the cluster assignments of each spot using \code{clusterPlot()} or
#' \code{enhancePlot()}, respectively.
#' 
#' @param sce A SingleCellExperiment object containing the spatial data and
#'   cluster assignments in \code{sce$spatial.cluster}.
#' @param sce.enhanced,sce.ref The enhanced subspots and original spots,
#'   respectively.
#' 
#' @return Both functions return a ggplot object. 
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce$spatial.cluster <- floor(runif(ncol(sce), 1, 4))
#' clusterPlot(sce)
#'
#' @name spatialPlot
NULL

palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
             "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9")

#' @importFrom ggplot2 ggplot aes_ geom_point scale_color_manual coord_equal labs theme_void
#'
#' @export
#' @rdname spatialPlot
clusterPlot <- function(sce, platform=c("Visium", "ST")) {
    # TODO: add platform/lattice to sce metadata instead of passing
    platform <- match.arg(platform)
    
    cdata <- data.frame(colData(sce))
    message("Making spots")
    if (platform == "Visium") {
        vertices <- .make_hex_spots(cdata)
    } else {
        vertices <- .make_square_spots(cdata)
    }
    
    # TODO: redo with geom_polygon/square for ST
    # TODO: redo with geom_polygon/hex so it scales and matches enhancePlot
    message("plotting spots")
    splot <- ggplot(data=vertices, 
                    aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) + 
        geom_polygon() +
        # scale_fill_manual() +
        labs(fill="Cluster") +
        coord_equal() +
        theme_void()
    
    splot
}

.make_hex_spots <- function(cdata) {
   NULL
}

.make_square_spots <- function(cdata, fill.col="spatial.cluster") {
    # Square spots are referenced by top left vertex
    spot_positions <- cdata[, c("col", "row", fill.col)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
    spot_positions$spot <- rownames(spot_positions)
    
    vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                                  y.offset=c(0, 0, 1, 1))
    
    spot_vertices <- merge(spot_positions, vertex_offsets)
    
    spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
    spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
    
    as.data.frame(spot_vertices)
}

#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_fixed labs theme_void
#'
#' @export
#' @rdname spatialPlot
enhancePlot <- function(sce.enhanced, sce.ref) {
## TODO add arbitrary column (e.g. marker expression) instead of cluster
## TODO maybe better to re-compute original positions from enhanced instead of passing sce.ref

    ## TODO store these as attributes when running spatialEnhance
    inputs <- .prepare_inputs(sce.ref)
    
    positions_plot <- .polygon_positions(inputs$positions, inputs$xdist*1.01, 
        inputs$ydist*1.01, data.frame(sce.enhanced$spatial.cluster))
    colnames(positions_plot) <- c("x", "y", "group", "spatial.cluster")

    eplot <- ggplot(positions_plot, aes_(x=~x, y=~(-y), group=~group, fill=~factor(spatial.cluster))) +
        geom_polygon() +
        # scale_fill_manual(values=palette) +
        coord_fixed(ratio=sqrt(1)) +
        labs(fill="Cluster") +
        theme_void()
        
    eplot
}

## (Written by Ed) Make polygon borders for each subspot
## 
## @param positions Dataframe of spot positions
## @param xdist,ydist Interspot distance
## @param cols Subspot level values (e.g. cluster label)
#' @importFrom stringr str_split
.polygon_positions <- function(positions, xdist, ydist, cols) {
    colnames(positions) <- c("x", "y")
    shift <- data.frame(Var1=c(0, 1,   0,   -1,  -1,   0,    1),
                       Var2=c(0, 1/3, 2/3, 1/3, -1/3, -2/3, -1/3))
    shift2 <- shift[c(1,2,3, #each row is one polygon
                     1,3,4,
                     1,6,7,
                     1,5,6,
                     1,2,7,
                     1,4,5),]
    group <- c(rep(seq(1, 6), each=3)) #six triangles
    
    n0 <- nrow(positions)
    positions_long <- positions[rep(seq_len(n0), 18), ]
    shift2 <- t(t(shift2)*c(xdist, ydist))
    shift_long <- shift2[rep(seq_len(18), each=n0), ]
    
    positions_long[,"x"] <- positions_long[,"x"] + shift_long[,"Var1"]
    positions_long[,"y"] <- positions_long[,"y"] + shift_long[,"Var2"]
    group_long <- paste(rep(seq_len(n0), 18),rep(rep(seq(1, 6), each=3), each=n0), sep="_")
    group_long2 <- data.frame(apply(str_split(group_long, "_", simplify=TRUE), 2, as.numeric))
    
    #add columns
    colnames(group_long2) <- c("j0", "poly")
    cols <- as.matrix(cols)
    cols_long <- matrix(nrow=nrow(positions_long), ncol=ncol(cols))
    colnames(cols_long) <- colnames(cols)
    for(i in seq_len(nrow(positions_long))){
        # print(i)
        # print(group_long2$j0[i]+(group_long2$poly[i]-1)*n0)
        cols_long[i,] <- cols[group_long2$j0[i], ]#+(group_long2$poly[i]-1)*n0,]
    }
    data.frame(positions_long, group=group_long, cols_long)
}
