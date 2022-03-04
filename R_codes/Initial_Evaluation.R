#Import libraries for spatial analysis
library(sf) #for opening and managing spatial data
library(gstat) #for gestatiscal analysis
library(fasterize) #for raster file types
library(tmap) #library to plot spatial data
library(tmaptools) #complementary library for tmap
library(raster) #deals with raster
library(optimParallel) #does the parallel processing
library(parallel)
library(ggpubr)
library(bestNormalize)
library(dplyr)

get_initial_vals = function(v_org){
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

  psill=v_org$gamma[as.numeric(rownames(exp_v_org_smoothed[pos_in_filtered,]))]

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
  result_l = c(nugget, psill, sv_dst)
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


#set our working environment
main_dir = 'H:/HyperLayer_Data/2020/Field_15-16/Soil_Samples/20200218/Processed'

#Search an open field boundary
field_boundary = read_sf(list.files(file.path('H:/PhD_Thesis/Chapt_6/Paper_JEEG/Raw_Data/Boundary'), pattern = '.shp', full.names = T))
field_boundary = to_utm(field_boundary)

#List files that contain the spatial information
sp_files = list.files(main_dir, pattern = '.gpkg', full.names = T, recursive = T)

#Open the data - work with 0-6in
sp_file_6in = grep('SoilSamples_6in', sp_files, value = T) #List the soil samples

soilsamples_6 = read_sf(sp_file_6in) #Reads the soil sample data


#Plot a map with both information
t = tmap::tm_shape(soilsamples_6, bbox = sf::st_bbox(sf::st_buffer(field_boundary, 100))) +
  tmap::tm_dots('K', n = 5, style = "quantile", palette = tmaptools::get_brewer_pal("RdYlGn", n = 5, plot = F), border.col = 'transparent', size = 0.5) +
  tmap::tm_layout('PH Soil Samples and Validation', asp =1.2, frame = FALSE, legend.show = T, title.position = c("left", "top"))+
  tmap::tm_shape(field_boundary) + tmap::tm_borders(lwd = 2)+
  tmap::tm_compass(position = c("right", "top"), size = 3) +
  tmap::tm_scale_bar(width = 0.2, position = c("right", "bottom"), text.size  = 0.7)
print(t)
dev.off()

hist(soilsamples_6$K)
shapiro.test(soilsamples_6$K)

#Box-cox transformation
soilsamples_6$K = bestNormalize::boxcox(soilsamples_6$K, standardize = F)$x.t

hist(soilsamples_6$K)
shapiro.test(soilsamples_6$K)


#Create our grid
buff_field = st_buffer(field_boundary, 30)

grid_res = 15

bbox = sf::st_bbox(buff_field)

#Based on the grid, create a raster and use fasterize (faster way to deal with
#raster fields) - All the other layers will be sampled to this one
orig <- bbox[1:2]

frst <- fasterize::raster(buff_field, orig, res = grid_res)
frst$id = 1:length(frst)


#Calculate the diagonal to the field
max_dst = sf::st_distance(sf::st_point(c(bbox[1], bbox[2])), sf::st_point(c(bbox[3], bbox[4])))[1]

#Number of samples taken
n_samples = nrow(soilsamples_6)


#Data variance... can be used for calculating the initial values for the variogram
varinit = var(soilsamples_6$K)

area = as.numeric(st_area(field_boundary))


###############################################################
##################### Initial values ##########################
###############################################################

#Minimum number of points - NEVER SMALLER THAN 4
min_numb_neigh = 4

range = max_dst

anot_range = sqrt(area)
lagdst_2 = sqrt(area/n_samples)


#Generated the first variogram, will be used for generating good starting values for nugget, range and psill
gOK = gstat::gstat(NULL,"K", K ~ 1, soilsamples_6, maxdist = max_dst)

models = c('Exp', 'Sph', 'Gau')
model = 'Sph'
models_l = list()
n = 1
for (model in models){
    #v = gstat::variogram(gOK, cutoff = var_option[1], width = var_option[2], alpha = c(anis_direction, anis_direction+90))
    v = gstat::variogram(gOK, cutoff = anot_range, width = lagdst_2)
    plot(v)

    #v_ = v[v$dir.hor == anis_direction+90,]
    #v_anis =v[v$dir.hor == anis_direction,]
    #rownames(v_)=NULL
    #rownames(v_anis)=NULL

    #init_ = get_initial_vals(v_)
    init_v = get_initial_vals(v)

    #anis_coeff = init_['Range']/init_vanis['Range']

    # if(anis_coeff<0){anis_coeff = 0}
    # if(anis_coeff>1){anis_coeff = 1}

    m = gstat::vgm(psill = init_v['Sill']-init_v['Nugget'], model = model, range = init_v['Range'],
                   nugget = init_v['Nugget'])

    m = fit.variogram(v, m, fit.method = 7,fit.sills = TRUE)
    cv_kr = krige.cv(K~1, soilsamples_6, m, nfold=10)
    cv_kr$SErr = cv_kr$residual^2
    plot(v, m)
    Fold_RMSE = cv_kr %>%
      group_by(fold) %>%
      summarise(RMSE = sqrt(mean(SErr, na.rm =T)))

    models_l[[paste0(model, n)]][['RMSE']] = mean(Fold_RMSE$RMSE)
    models_l[[paste0(model, n)]][['model']] = m
    models_l[[paste0(model, n)]][['variogram']] = v

    n = n+1
}
errors_var = sapply(models_l, function(x) return(attr(x$model, 'SSErr')))
interpol_RMSE = sapply(models_l, function(x) return(x$RMSE))


#Variogram error export
error_l = list(errors_var, interpol_RMSE)
names(error_l) = c('Variogram_SSErr', 'CV_RMSE')

error_df = as.data.frame(do.call(rbind, error_l))

var_min = which.min(errors_var)
RMSE_min = which.min(interpol_RMSE)


selected_model = models_l[[var_min]]


m = selected_model$model
v = selected_model$variogram


# pts_sp = as(pts, 'Spatial')
# auto <- HLDP::autofitVarHLDP(formula = V1 ~ 1,input_data = pts_sp,GLS.model = c('Sph', 'Exp', 'Gau'),
#                                start_vals = c(nugget, sv_dst, psill), max_dst = max_dst)
#
# v = auto$exp_var
# m = auto$var_model
#
# plot(auto$exp_var)


preds = gstat::variogramLine(m, maxdist = max(v$dist))
v$dir.hor = round(v$dir.hor, 0)
p = ggplot2::ggplot() +
  ggplot2::geom_point(data = v, ggplot2::aes(x = dist, y = gamma, colour = np), size = 3) +
  ggplot2::geom_line(data = preds, ggplot2::aes(x = dist, y = gamma))+
  viridis::scale_color_viridis(direction = -1)+
  ggplot2::labs(colour = 'Number of Points', x = "Distance", y = "Semivariance", title = "K")+
  ggplot2::annotate("text", -Inf, Inf , label=paste0(
    'Model = ', m$model[2], '\n',
    'Nugget = ', round(m$psill[1], 4), '\n',
    'Sill = ',  round(m$psill[2]+m$psill[1], 4), '\n',
    'Range = ',  round(m$range[2], 2), '\n'
  ),color="dark gray",hjust = -0.7, vjust = 3)+
  expand_limits(x = 0, y = 0)
p

if(any(duplicated(sf::st_geometry(soilsamples_6)))){
  soilsamples_6 <- soilsamples_6[!duplicated(sf::st_geometry(soilsamples_6)),]
}

grd = stars::st_as_stars(frst)

# Calculate the number of cores
no_cores <- parallel::detectCores() - 1
size = round(nrow(grd)/no_cores, 0)
#no_cores = round(nrow(grd)/size, 0)

i=1

splits_l = list()
for (core in 1:no_cores){

  if(core == no_cores){
    x = i:nrow(grd)
  }else{
    x = i:(i+size-1)
  }
  i = i+size
  splits_l[[core]] = x
}


grd_splits = lapply(1:length(splits_l), function(x){
  grd_slice = dplyr::slice(grd, "x", splits_l[[x]])
  sf::st_crs(grd_slice) = sf::st_crs(soilsamples_6)
  return(grd_slice)
})

#tictoc::tic()
OK_l = future.apply::future_lapply(grd_splits, FUN = gstat::krige, locations = soilsamples_6, formula = K ~ 1, model = m
                                   , future.seed=TRUE)
#tictoc::toc()
#print('Done Interpolation')

OK = do.call(c, args = c(OK_l, along = c('x')))

#print('Done converting')
#plot(OK)

sf::st_crs(OK) = sf::st_crs(soilsamples_6)

fname = paste0(tempfile(), ".tif")
stars::write_stars(OK['var1.pred'], fname)
from = fname

OK_f = raster::stack(from)
OK_f = raster::readAll(OK_f)

raster::extent(OK_f) = raster::extent(frst)


plot(OK_f)

#Inverse Distance Interpolation

#Guarantee that the raster file and the vector files have the same
# crs.
raster::crs(frst) = sf::st_crs(soilsamples_6)$wkt

n_samples = nrow(soilsamples_6)

variable ='PH'

f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat::gstat(formula=as.formula(paste0(variable,'~1')), locations=train, nmax=nmx, set=list(idp=idp))
  p <- raster::predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(as.list(sf::st_drop_geometry(test)[,variable])[[variable]], p)
}

i <- sample(nrow(soilsamples_6), 0.2 * nrow(soilsamples_6))
tst <- soilsamples_6[i,]
trn <- soilsamples_6[-i,]

cl <- parallel::makeCluster(spec=parallel::detectCores(), outfile="")
parallel::clusterExport(cl=cl, list("variable", "RMSE"),
                        envir=environment())
parallel::clusterEvalQ(cl,c(library("sf"), library(gstat)))
parallel::setDefaultCluster(cl=cl)
opt <- optimParallel::optimParallel(c(4, 1), f1, test=tst, train=trn, parallel=list(forward=TRUE))
setDefaultCluster(cl=NULL); stopCluster(cl)


# Setting nmax to 5 and idp to 1 is an inverse weighted interpolation:
gs = gstat::gstat(formula=as.formula(paste0(variable,'~1')),
                  locations=soilsamples_6, nmax=opt$par[1], set=list(idp = opt$par[2]))

nv = raster::interpolate(frst, gs)
plot(nv)


val_analysis = raster::extract(nv, as(validation_6, 'Spatial'), exact = T)

validation_6$PH_pred = val_analysis

plot(validation_6$PH, validation_6$PH_pred)
cor.test(validation_6$PH, validation_6$PH_pred)


ggscatter(validation_6, x = "PH", y = "PH_pred",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Observed pH", ylab = "Predicted pH")
