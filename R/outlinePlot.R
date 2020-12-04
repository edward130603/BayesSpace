#' Add outlines to a cluster/feature plot
#' TODO add `outline` and `outline_groups` (with default all) to c/fPlot
#' 
#' @param sce SingleCellExperiment.
#' @param base_plot Base cluster/feature plot to add outlines to.
#' @param outline_group Group to outline; i.e. a categorical column in
#'   \code{colData(sce)}.
#' @param outline_values Subset of values in \code{outline_group} to outline; by
#'   default all groups/clusters will be outlined.
#' @param outline_title Title of outline groups in plot legend.
#' @param outline_size Width of outline.
#' @param outline_nudge Scale the distance each vertex in the outline is nudged
#'   towards the interior of the cluster.
#' 
#' @return ggplot object
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
#' @importFrom purrr map
#' @importFrom ggplot2 geom_polygon aes_ labs
#' @export
outlinePlot <- function(sce, base_plot, 
                        outline_group="spatial.cluster",
                        outline_values=NULL,
                        outline_title="Cluster",
                        outline_size=1.2,
                        outline_nudge=0.5) {
    
    ## Get vertices like in clusterPlot(), except fill by component
    cluster_components <- .find_connected_components(sce, outline_group)
    vertices <- .make_vertices(sce,
                               cluster_components,
                               .bsData(sce, "platform"),
                               .bsData(sce, "is.enhanced"))
    
    ## Clarify memberships
    vertices$component_group <- vertices$fill
    vertices$outline_group <- colData(sce)[vertices$spot, outline_group]
   
    ## Subset components to outline if necessary
    if (is.null(outline_values)) {
        components_to_outline <- unique(cluster_components)
    } else {
        components_to_outline <- vertices %>%
            filter(.data$outline_group %in% outline_values) %>%
            pull(.data$component_group) %>%
            unique()
    }
   
    ## Get the spot vertices associated with each component,
    ## and obtain their outlines
    .make_each_outline <- function(component_label) {
        component <- vertices %>% filter(.data$component_group == component_label)
        boundary <- .make_outline(component)
        boundary <- .nudge_outline(boundary, outline_nudge)
        
        ## Add grouping labels for plotting with geom_polygon()
        boundary$component_group <- component_label
        boundary$outline_group <- as.character(component$outline_group)[[1]]
        
        boundary
    }
    
    ## Get boundary around each component
    outlines <- map(components_to_outline, .make_each_outline)
    outlines <- do.call(rbind, outlines)
    
    plot <- base_plot + 
        geom_polygon(aes_(x=~x_nudge, y=~y_nudge, group=~component_group, color=~outline_group),
                     data=outlines,
                     size=outline_size,
                     fill=NA) +
        labs(color=outline_title)
    
    plot
}

#' Nudge outline into interior of component, to avoid overlaps between clusters
#' 
#' Moves each boundary vertex \code{nudge} of the distance towards the midpoint
#' of its neighbors. Distance is negated if the midpoint is outside the boundary
#' 
#' @param boundary Table of (x, y) coordinates
#' @param nudge Fraction of distance to midpoint of neighbors
#' @return Table of nudged coordinates at \code{(x_nudge, y_nudge)}
#' 
#' @importFrom sp point.in.polygon
.nudge_outline <- function(boundary, nudge=0.25) {

    ## Get coordinates of neighboring spots and compute distance to move
    ## towards midpoint between neighbors
    .get_nudge_dist <- function(coord) {
        c_prev <- .roll(coord, 1)
        c_next <- .roll(coord, -1)
        mid <- (c_prev + c_next) / 2
        dist <- nudge * (mid - coord)
    }
    boundary$x_dist <- .get_nudge_dist(boundary$x_rd)
    boundary$y_dist <- .get_nudge_dist(boundary$y_rd)

    ## Check if nudged vertices are inside boundary,
    ## and flip direction of nudge if they aren't
    boundary$is_interior <- point.in.polygon(boundary$x_rd + boundary$x_dist,
                                             boundary$y_rd + boundary$y_dist,
                                             boundary$x_rd,
                                             boundary$y_rd)

    .flip_sign <- function(coord, dist, loc) {
        ifelse(loc == 0, coord - dist, coord + dist)
    }
    boundary$x_nudge <- .flip_sign(boundary$x_rd, boundary$x_dist, boundary$is_interior)
    boundary$y_nudge <- .flip_sign(boundary$y_rd, boundary$y_dist, boundary$is_interior)

    boundary
}

#' Find connected components of each cluster in spatial context
#' 
#' Consider spots connected if they are neighbors and belong to same cluster
#' 
#' @param sce SingleCellExperiment
#' @param group Categorical column in \code{colData(sce)} to group spots by
#' 
#' @return Vector of component memberships
#' 
#' @keywords internal
#' 
#' @importFrom purrr imap
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#' @importFrom igraph graph_from_edgelist components
.find_connected_components <- function(sce, group="spatial.cluster") {
    if (.bsData(sce, "is.enhanced")) {
        ## TODO
        ## neighbors <- find_neighbors()
        stop("not yet implemented")
    } else {
        neighbors <- .find_neighbors(sce, .bsData(sce, "platform"), verbose=FALSE)
    }
    
    ## Convert list of lists of neighbors to dataframe of (spot, neighbor) pairs
    ## Need to return empty dataframes for neighborless-spots to preserve
    ## spot-indexing
    to_frame <- function(nbrs, idx) {
        if (length(nbrs) > 0) {
            data.frame(spot=idx, neighbor=nbrs)
        } else {
            data.frame(spot=c(), neighbor=c())
        }
    } 
    edges <- do.call(rbind, imap(neighbors, to_frame))
    
    ## Neighbor indices are zero-indexed for C++
    edges$neighbor <- edges$neighbor + 1
    
    ## Add self edges to ensure any spot without neighbors is included
    self_edges <- data.frame(spot=seq_len(ncol(sce)), 
                             neighbor=seq_len(ncol(sce)))
    edges <- rbind(edges, self_edges)
    
    ## Filter to neighbors belonging to the same cluster
    edges$spot_group <- colData(sce)[edges$spot, group]
    edges$neighbor_group <- colData(sce)[edges$neighbor, group]
    edges <- edges %>% 
        filter(.data$spot_group == .data$neighbor_group) %>%
        select(.data$spot, .data$neighbor)

    ## Find connected components 
    graph <- graph_from_edgelist(as.matrix(edges), directed=FALSE)
    comps <- components(graph)
    
    comps$membership
}

#' Obtain vertices outlining a connected component
#' 
#' @param component Table with \code{(x.vertex, y.vertex)} coordinates
#' @return Vertices from boundary of \code{component}, ordered for plotting with
#'   \code{geom_polygon()}
#'
#' @keywords internal
.make_outline <- function(component) {
    boundary <- .get_boundary_vertices(component)
    path <- .trace_path(boundary)

    ## Reorder boundary vertices for geom_path()
    boundary <- boundary[path, ]
}

#' Select boundary vertices from a connected component
#' 
#' A vertex appears once in each spot outline. Since adjacent spots coincide in
#' our plotting functions, we can determine whether the vertex is at the
#' boundary of the component by counting the number of spots it appears in. A
#' vertex that appears in the maximum possible number of spots at an
#' intersection (i.e. 3 for Visium and 4 for ST) is in the interior of the
#' component; any vertex that appears in fewer than the maximum is on the
#' boundary.
#' 
#' @param component Subset of vertices from \code{.make_vertices()}
#' @param max_neighbors Maximum number of neighboring spots for a spot to be
#'   considered on the boundary.
#'   
#' @return A data.frame of x and y coordinates for each boundary vertex
#' 
#' @keywords internal
#'
#' @importFrom magrittr %>% 
#' @importFrom dplyr mutate group_by summarise n filter select ungroup
.get_boundary_vertices <- function(vertices, max_neighbors=2) {
    vertices %>% 
        mutate(x_rd=round(.data$x.vertex, 3),
               y_rd=round(.data$y.vertex, 3)) %>%
        group_by(.data$x_rd, .data$y_rd) %>%
        summarise(n=n()) %>%
        filter(.data$n <= max_neighbors) %>% 
        select(.data$x_rd, .data$y_rd) %>%
        ungroup()  # necessary, otherwise roll/lag misbehave
}


#' Shift a vector forwards or backwards
#' 
#' Like \code{dplyr::lag()} or \code{dplyr::lead()}, but values from the
#' head/tail of the vector are rolled over to the end/start instead of being
#' replaced with a default value.
#' 
#' @param x Vector to roll.
#' @param n Number of elements to shift vector by. If positive, elements are
#'   moved forwards; if negative, backwards.
#'   
#' @return
#' A shifted vector \code{y}, where \code{y[i] = x[(i - n) \% n]}.
#' 
#' @keywords internal
#' 
#' @importFrom vctrs vec_size vec_slice vec_c
.roll <- function(x, n=1L) {
    xlen <- vec_size(x)
    n <- n %% xlen
    
    if (n == 0)
        return(x)
    
    vec_c(
        vec_slice(x, seq(xlen - n + 1, xlen)),
        vec_slice(x, seq_len(xlen - n))
    )
}

#' Trace a path clockwise through a set of vertices
#' 
#' Used to arrange boundary vertices for \code{geom_path()}
#' 
#' @param X matrix of (x, y) coordinates for vertices
#' @param max_dist Maximum distance for two vertices to be considered neighbors.
#'   Defaults to 0.6 for Visium; TODO write for ST
#' 
#' @return A vector of indices into the rows of \code{X}, ordered such that the
#'   coordinates of \code{X} may be traced into a polygon.
#'   
#' @keywords internal
#' 
#' @importFrom stats dist
#' @importFrom igraph graph_from_adjacency_matrix neighbors
.trace_path <- function(X, max_dist=0.6) {
    X <- as.matrix(X)
    
    ## TODO: add row-wise filter to choose neighbors closer then second nearest dist
    adj <- as.matrix((dist(X)))
    adj <- (adj < max_dist)
    graph <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
    
    ## Initialize start with only two neighbors 
    for (i in seq_len(nrow(X))) {
        if (length(neighbors(graph, i)) == 2) {
            src <- i
            break
        }
    }
    
    ## Initialize start
    ## To proceed CW, choose neighbor to right of start and preceding neighbor
    n1 <- neighbors(graph, src)[1]
    n2 <- neighbors(graph, src)[2]
    prev <- src
    curr <- ifelse(.is_left(X[n1, ], X[prev, ], X[n2, ]), n1, n2)
    path <- c(prev, curr)
    
    iters <- 0
    while (curr != src) {
        nbrs <- neighbors(graph, curr)
        if (length(nbrs) == 2) {
            if (nbrs[1] == prev) {
                prev <- curr
                curr <- nbrs[2]
            } else {
                prev <- curr
                curr <- nbrs[1]
            }
            path <- c(path, unname(curr))
        } else {
            for (neighbor in nbrs) {
                if (neighbor == prev) {
                    next
                }
                
                if (.is_left(X[prev, ], X[curr, ], X[neighbor, ])) {
                    prev <- curr
                    curr <- neighbor
                    path <- c(path, curr)
                    break
                }
            }
        }
        
        iters <- iters + 1
        if (iters > nrow(X)) {
            warning("Too many iterations")
            return(path)
        }
    }
    
    path
}


#' Check if \code{pt} is to left of line determined by \code{x1} and \code{x2}
#' 
#' Used to determine if a vertex in the boundary is to the left or right of the
#' preceding path.
#' 
#' @param x1,x2,pt (x,y) coordinates
#' 
#' @return \code{TRUE} if \code{pt} is to the left of \code{(x1, x2)} from the
#'   perspective of \code{x1}. Otherwise \code{FALSE}.
#' 
#' @details
#' \code{pt} is to the left of \code{(x1, x2)} if the sign of the determinant of
#' \code{(x2 - x1)} and \code{(pt - x1)} (the lines drawn from \code{x1} to each
#' of its neighbors) is positive.
#' 
#' @keywords internal
.is_left <- function(x1, x2, pt) {
    mat <- cbind(x2 - x1, pt - x1)
    det(mat) > 0
}
