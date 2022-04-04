library(sf)


cv <- function(df_values, na.rm = T){
  (sd(df_values, na.rm = na.rm)/mean(df_values, na.rm=na.rm))*100
}
'%notin%' <- Negate('%in%')


#Function that automatically gets initial values for the variogram fitting
get_initial_vals = function(v_org, data_variance = NULL, method = 'default'){

  if(method == 'default'){

    #Calculate inital nugget

    nugget_init = v_org$gamma[1] - ((v_org$dist[1]/(v_org$dist[2]-v_org$dist[1])*(v_org$gamma[2] - v_org$gamma[1])))

    nugget = max(c(0, nugget_init))

    number_of_lags = nrow(v_org)
    sill = (v_org$gamma[number_of_lags]+v_org$gamma[number_of_lags-1]+v_org$gamma[number_of_lags-2])/3

    sv_dst = v_org$dist[number_of_lags]/2
  }

  if(method == 'spline'){
    #v_org = vars[[2]]
    lst = lapply(1:(nrow(v_org)), function(x) mean(v_org$gamma[(x + -2:2)[(x + -2:2)>0]], na.rm =T))

    #Create a smoothed variogram
    smoothing_median = do.call(rbind, lst)
    exp_v_org_smoothed = v_org
    exp_v_org_smoothed$gamma = smoothing_median
    exp_v_org_smoothed = exp_v_org_smoothed[exp_v_org_smoothed$np > 100,]



    #From the smoothed variogram generate a spline which will be used to find a good inflection point
    t= smooth.spline(exp_v_org_smoothed$dist, exp_v_org_smoothed$gamma, spar=0.4)

    #plot(t)


    #Autofit has a tendency to fit badly some data with high nugget with the
    # providing starting values


    #Obtain the inflection point

    second_dev = diff(diff(predict(t)$y))


    infl <- c(FALSE,diff(sign(second_dev))!=0)

    #It is a bit rare, but if any inflection is found the minimon absolute value (closest to 0)
    #is used for range and psill
    if(all(infl == F)){
      pos_in_filtered = which.min(abs(diff(predict(t)$y)))
    }else{
      pos_in_filtered = min(which(infl == TRUE))
    }

    #Obtain range and sill
    sv_dst = v_org$dist[as.numeric(rownames(exp_v_org_smoothed[pos_in_filtered,]))]

    sill=v_org$gamma[as.numeric(rownames(exp_v_org_smoothed[pos_in_filtered,]))]

    #Nugget is obtained based on the same calculation used on QGIS Smart-Map plugin
    nugget = (exp_v_org_smoothed$gamma[2]*exp_v_org_smoothed$dist[1]-exp_v_org_smoothed$gamma[1]*exp_v_org_smoothed$dist[2])/(exp_v_org_smoothed$dist[1]-exp_v_org_smoothed$dist[2])
    #nugget = (v_org$gamma[2]*v_org$dist[1]-v_org$gamma[1]*v_org$dist[2])/(v_org$dist[1]-v_org$dist[2])

    #In case of negative nugget or nugget higher than psill,
    #it tries to find another nugget
    if(nugget <0|nugget > psill){
      nugget = exp_v_org_smoothed$gamma[1]
      if(nugget > psill){
        nugget = exp_v_org_smoothed$gamma[2]
      }
      if(nugget >= psill){
        nugget = as.numeric(0.2*varinit)
      }
      if(nugget >= psill){
        nugget = min(exp_v_org_smoothed$gamma)
      }
    }
  }

  if(method == 'variance'){

    nugget = v_org$gamma[1]
    sill = data_variance
    range_gamma = which.min(abs(v_org$gamma - sill))
    sv_dst = v_org$dist[range_gamma]
  }


  result_l = c(nugget, sill, sv_dst)
  names(result_l) = c('Nugget', 'Sill', 'Range')
  return(result_l)
}

#' Transform any projection to UTM for a sf object:
#'
#' This take in any sf geometries and returns
#'  a new sf geometry with the right UTM zone:
#' @param sf_obj  of class 'sf'
#' @keywords UTM, Projection, Simple Features, sf
#' @export
#' @examples
#' \dontrun{
#' to_utm(fields)
#' }
to_utm <- function(sf_obj) {
  # Function to get UTM Zone from mean longitude:
  long2UTM <- function(long) {
    (floor((long + 180) / 6) %% 60) + 1
  }

  # Check if the object class is 'sf':
  obj_c <- class(sf_obj)[1]
  if (obj_c == "sf") {
    # In case the object has no projectin assigned,
    #  assume it to geographic WGS84 :
    if (is.na(sf::st_crs(sf_obj))) {
      sf::st_crs(sf_obj) <- sf::st_crs(4326)
    }

    # Get the center longitude in degrees:
    bb <- sf::st_as_sfc(sf::st_bbox(sf_obj))
    bb <- sf::st_transform(bb, sf::st_crs(4326))

    # Get UTM Zone from mean longitude:
    utmzone <- long2UTM(mean(sf::st_bbox(bb)[c(1, 3)]))

    # Get the hemisphere based on the latitude:
    NS <- 100 * (6 + (mean(sf::st_bbox(bb)[c(2, 4)]) < 0))

    # Add all toghether to get the EPSG code:
    projutm <- sf::st_crs(32000 + NS + utmzone)

    # Reproject data:
    sf_obj <- sf::st_transform(sf_obj, projutm)
    return(sf_obj)
  } else {
    options(error = NULL)
    stop("Object class is not 'sf', please insert a sf object!")
  }
}


#' Calculates Root Mean Squared Error (RMSE):
#'
#' This take in any two list of observations and calculate RMSE
#'
#' @param observed list of ground truth values
#' @param predited list of predicted values
#' @keywords math, rmse
#' @export
#' @examples
#' \dontrun{
#'
#' }

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}


#' Function to get distances between consecutive points:
#'
#' This goes through every point in sf_obj and sums the distances between each consecutive point. If parameter max is FALSE sf_obj should
#' be a point feature. If max = TRUE, sf_obj should be a polygon and it will return the diagonal distance of the field.
#' @param sf_obj object of class 'sf'
#' @param max conditional variable - default = FALSE
#' @keywords Distance, Simple Features, sf
#' @export
#' @examples
#'
calc_distance <- function(sf_obj, max = FALSE) {
  if(max == TRUE){
    pts = sf::st_cast(sf::st_as_sfc(sf::st_bbox(sf_obj)), "POINT")
    return(as.numeric(sf::st_distance(pts[1], pts[3])))
  }else{
    ab_dist <- function(a, b) {
      return(sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2))
    }
    coords <- sf::st_coordinates(sf_obj)
    result <- as.numeric(sapply(c(2:nrow(coords)), function(x) {
      ab_dist(coords[x - 1, ], coords[x, ])
    }))
    return(c(result[1], result))
  }

}
