#' Add outlines to a cluster/feature plot
#' TODO add `outline` and `outline_groups` (with default all) to c/fPlot
#' 
#' @param sce SingleCellExperiment
#' @param base_plot Base cluster/feature plot to add outlines to
#' @param outline Group to outline; i.e. a categorical column in
#'   \code{colData(sce)}
#' @param outline_values Optional subset of groups to outline; by default all
#'   groups will be outlined
#' 
#' @return ggplot object
#' 
#'
outlinePlot <- function(sce, base_plot, 
                        outline="spatial.cluster",
                        outline_values=NULL) {
    
    cluster_components <- .find_connected_components(sce, outline)
    vertices <- .make_vertices(sce, cluster_components,
                               .bsData(sce, "platform"), 
                               .bsData(sce, "is.enhanced"))
    
    ## TODO do for all clusters/components
    ## TODO add nudge
    boundary <- .make_outline(vertices %>% filter(fill == 1))
    
    base_plot + 
        geom_point(boundary, )
    
    
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
#' @importFrom purrr imap keep
#' @importFrom igraph graph_from_edgelist components
.find_connected_components <- function(sce, group="spatial.cluster") {
    if (.bsData(sce, "is.enhanced")) {
        ## TODO
        ## neighbors <- find_neighbors()
        stop("not yet implemented")
    } else {
        neighbors <- .find_neighbors(sce, .bsData(sce, "platform"))
    }
    
    ## Convert list of lists of neighbors to dataframe of (spot, neighbor) pairs
    edges <- keep(neighbors, ~ length(.x) > 0)
    edges <- imap(edges, ~ data.frame(spot=.y, neighbor=.x))
    edges <- do.call(rbind, edges)
    
    ## Neighbor indices are zero-indexed for C++
    edges$neighbor <- edges$neighbor + 1
    
    ## Add self edges to ensure any spot without neighbors is included
    self_edges <- data.frame(spot=seq_len(ncol(sce)), 
                             neighbor=seq_len(ncol(sce)))
    edges <- rbind(edges, self_edges)
    
    ## match groups
    edges$spot_group <- colData(sce)[edges$spot, group]
    edges$neighbor_group <- colData(sce)[edges$neighbor, group]
    edges <- edges %>% 
        filter(spot_group == neighbor_group) %>%
        select(spot, neighbor)

    ## Find connected components 
    graph <- graph_from_edgelist(as.matrix(edges), directed=FALSE)
    comps <- components(graph)
    
    comps$membership
}

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
#' @param component Subset of vertices from \link{\code{.make_vertices()}}
#' @param max_neighbors Maximum number of neighboring spots for a spot to be
#'   considered on the boundary.
#'   
#' @return A data.frame of x and y coordinates for each boundary vertex
#' 
#' @keywords internal
#' 
#' @importFrom dplyr mutate group_by summarise filter select ungroup
.get_boundary_vertices <- function(component, max_neighbors=2) {
    vertices %>% 
        mutate(x_rd=round(x.vertex, 3),
               y_rd=round(y.vertex, 3)) %>%
        group_by(x_rd, y_rd) %>%
        summarise(n=n()) %>%
        filter(n <= max_neighbors) %>% 
        select(x_rd, y_rd) %>%
        ungroup()  ## necessary otherwise roll/lag misbehave
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
#' @examples
#' roll(1:5)
#' roll(1:5, -1)
#' 
#' @keywords internal
#' 
#' @importFrom vctrs vec_size vec_slice vec_c
roll <- function(x, n=1L) {
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
#' @keywords internal
#' 
#' @importFrom igraph graph_from_adjacency_matrix neighbors
.trace_path <- function(X) {
    X <- as.matrix(X)
    
    ## TODO: add row-wise filter to choose neighbors closer then second nearest dist
    adj <- as.matrix((dist(X)))
    adj <- (adj < 0.6)
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
            message("Too many iterations")
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
