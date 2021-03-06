% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialPlot.R
\name{.make_spot_vertices}
\alias{.make_spot_vertices}
\title{Compute vertex coordinates for each spot in frame of plot}
\usage{
.make_spot_vertices(spot_positions, vertex_offsets)
}
\arguments{
\item{spot_positions}{Center for hex, top left for square}

\item{vertex_offsets}{Data frame of (x, y) offsets wrt spot position for each
vertex of spot}
}
\value{
Cartesian product of positions and offsets, with coordinates
  computed as (pos + offset)
}
\description{
Compute vertex coordinates for each spot in frame of plot
}
\keyword{internal}
