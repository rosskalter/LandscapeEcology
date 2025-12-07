convexHull <- function(pts){
  hull1.idx <- chull(pts)
  hull <- Polygon(pts[hull1.idx,])
  hull.list <- Polygons(list(hull), ID="hull")
  hull.SP <- SpatialPolygons(list(hull.list))
  return(hull.SP)
}
