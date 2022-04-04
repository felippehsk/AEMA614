#########################################################################
#########################################################################
##                                                                     ##
##           AEMA 614: Temporal and Spatial Statistics                 ##
##    Case Study: Evaluation of Kriging Methods for Soil Samples       ##
##                                                                     ##
#########################################################################
#########################################################################


#########################################################################
######################## Importing Libraries ############################
#########################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('./functions.R')

library(sf) #for opening and managing spatial data
library(gstat) #for geostatiscal analysis
library(fasterize) #for raster file types
library(tmap) #library to plot spatial data
library(tmaptools) #complementary library for tmap
library(raster) #deals with raster
library(optimParallel) #does the parallel processing
library(ggpubr)
library(bestNormalize)#BoxCox
library(dplyr)

#########################################################################
############################ Set Directories ############################
#########################################################################

#Saving Directory
s_dir = '../results'
dir.create(s_dir, showWarnings = F)

#########################################################################
####################### Setting Initial Parameters ######################
#########################################################################

#Grid Size
grid_res = 15

#Variables to interpolate and their respective names to the maps
variables = c("K", "p_1", 'AL_M3', 'OM')
variables_press = c('K_ppm', 'P_ppm', 'Al_ppm', 'OM_Percent')

#Get value for local neighbors
local_scenarios = read.csv('../R_codes/Combinations_Local.csv')

#List of selected models
models = c('Exp', 'Sph', 'Gau')


#########################################################################
############################# Open Data #################################
#########################################################################

#Search an open field boundary
field_boundary = read_sf(list.files(file.path('../data/boundary'), pattern = '.shp', full.names = T))
field_boundary = to_utm(field_boundary)

#List files that contain the spatial information
sp_files = list.files('../data/soil_samples', pattern = '.gpkg$', full.names = T, recursive = T)

#Open the data - work with 0-6in
sp_file_6in = grep('SoilSamples_6in', sp_files, value = T) #List the soil samples

soilsamples_6 = read_sf(sp_file_6in) #Reads the soil sample data

#Checks and remove duplicates
if(any(duplicated(sf::st_geometry(soilsamples_6)))){
  soilsamples_6 <- soilsamples_6[!duplicated(sf::st_geometry(soilsamples_6)),]
}

#########################################################################
####################### Set Interpolation Grid ##########################
#########################################################################

#Adds a buffer to the boundary before creating grid
buff_field = st_buffer(field_boundary, 30)

#Gets bounding box
bbox = sf::st_bbox(buff_field)

#Based on the grid, create a raster and use fasterize (faster way to deal with
#raster) - All the kriging predictions will be performed over this raster
# Size is set at the parameters
orig <- bbox[1:2]

frst <- fasterize::raster(buff_field, orig, res = grid_res)
frst$id = 1:length(frst)

#Transfor the grid to stars obejct - important for predicting with gstat
grd = stars::st_as_stars(frst)


#########################################################################
############ Calculate Max. Dst. and Lag for Variogram  #################
#########################################################################

#Maximum distance for variogram
max_dst = sqrt(as.numeric(st_area(field_boundary)))/1.2

#Calculate the distance between points and set it as lag
lag_dst = median(calc_distance(soilsamples_6))


#########################################################################
################# Main Loop - Different Variables #######################
#########################################################################

var_pos = 3
#Create the lists that will receive the results
results_l = list()
lambda_bc_l = list()
models_chosen_l = list()

for(var_pos in 1:length(variables)){

  #Select the variable from list
  var = variables[var_pos]
  #Here we use a set of variables to better look on graphs
  var_press = variables_press[var_pos]

  #Select only the data related to the variable
  var_data = soilsamples_6[,var]

  #Remove spatial from data frame - for calculations
  df_var = as.data.frame(st_drop_geometry(soilsamples_6[,var]))

  #Some summary results for map
  mean_var = mean(df_var[,var], na.rm = T)
  coef_var = cv(df_var[,var], na.rm = T)
  std_dev = sd(df_var[,var], na.rm = T)

  #Map title
  map_title = paste0('Sample Results: ', gsub('_', ' ', var_press), '\nMean: ', round(mean_var,2), '\nCV(%): ', round(coef_var,2), '\nSD: ',
                     round(std_dev,2))

  #Create map directory
  point_map_dir = file.path(s_dir, 'Point_Data')
  dir.create(point_map_dir, showWarnings = F)

  #Plot a map with both information
  t = tmap::tm_shape(soilsamples_6, bbox = sf::st_bbox(sf::st_buffer(field_boundary, 100))) +
    tmap::tm_dots(var, n = 5, style = "quantile", palette = tmaptools::get_brewer_pal("RdYlGn", n = 5, plot = F), border.col = 'transparent', size = 0.5, title = gsub('_', ' ', var_press)) +
    tmap::tm_layout(map_title, asp =1.2, frame = FALSE, legend.show = T, title.position = c("left", "top"))+
    tmap::tm_shape(field_boundary) + tmap::tm_borders(lwd = 2)+
    tmap::tm_compass(position = c("right", "top"), size = 3) +
    tmap::tm_scale_bar(width = 0.2, position = c("right", "bottom"), text.size  = 0.7)

  tmap::tmap_save(t, file.path(point_map_dir, paste0(var_press, '_Point_Map.png')))


  #Box-cox transformation

  #Run statistical test to know if normal or not
  norm_t_before = shapiro.test(df_var[,var])
  is_norm = norm_t_before$p.value > 0.05


  NORM = F
  STAND = F

  #If non normal (no significant shapiro test) run BoxCox otherwise just
  #standardize
  if(!is_norm){
    NORM = T
    print(paste0('BoxCoxing for ', var))

    #Calculate best box cox lambda and apply transformation with stand.
    bc_trans = bestNormalize::boxcox(df_var[,var], standardize = F)
    #Substitute the new values to the data set
    var_data[,var] = bc_trans$x.t
    #Save lambda value for analysis
    lambda_bc_l[[var]] = bc_trans$lambda
  }
  # else{
  #   print(paste0('STDing for ', var))
  #   STAND = T
  #   #If not boxcox, just stand.
  #   std_data = scale(df_var[,var])
  #   #Substitute the new values to the data set
  #   var_data[,var] = as.data.frame(std_data)
  # }


  #########################################################################
  ################ Sub Process - Universal and Ordinary ###################
  #########################################################################

  #Create a list of types of kriging and their formulas
  UOK_l =  list()
  UOK_l[['UK']] = as.formula(paste0(var, '~ cos((2*pi*x)/350)+y'))
  UOK_l[['OK']] = as.formula(paste0(var, '~ 1'))

  formula_pos = 1
  #For each one of the kriging methods run all the available scenarios
  for(formula_pos in 1:length(UOK_l)){

    #Get kriging type
    krig_type = names(UOK_l)[formula_pos]
    #Get kriging type formula
    formula_ = UOK_l[[formula_pos]]

    #If UK add the x and y values to the dataframe
    if(krig_type == 'UK'){
      var_data$x = st_coordinates(var_data)[,1]
      var_data$y = st_coordinates(var_data)[,2]
    }

    #trend<-krige(formula_, var_data, grd, model=NULL)
    #plot(trend)

    # if(var == "AL_M3" & krig_type == 'OK'){
    #   max_dst = sqrt(as.numeric(st_area(field_boundary)))/2
    # }else{
    #   max_dst = sqrt(as.numeric(st_area(field_boundary)))
    # }

    #Set the gstat object with the data, formula and max distance for variogram
    gOK = gstat::gstat(NULL,var, formula_, var_data, maxdist = max_dst)

    #Calculate variogram
    v = gstat::variogram(gOK, cutoff = max_dst, width = lag_dst)
    # plot(v)
    # dev.off()

    #Obtain the initial values - same approach as SAS - Jian, Olea and Yu. (1996)
    init_v = get_initial_vals(v)

    #########################################################################
    ##################  Loop - Find Best Fitting Model ######################
    #########################################################################
    model = 'Per'
    #Create list to save models results
    models_l = list()
    for (model in models){
      #Create a model object using the inital values
      m = gstat::vgm(psill = init_v['Sill']-init_v['Nugget'], model = model, range = init_v['Range'],
                     nugget = init_v['Nugget'])

      #Use the initial values in the fitting function  - method used is WLS
      #when no convergence - save the warning and set error as NA
      m = tryCatch(fit.variogram(v, m, fit.method = 7,fit.sills = TRUE), warning=function(w) w )

      if(length(grep('No convergence|singular model', unlist(m))) == 0){
        print(paste0('Fitted ', model, ' for variable ', var))
        preds = gstat::variogramLine(m, maxdist = max(v$dist))
        preds$Model = model

        #Save the model
        models_l[[paste0(model)]][['model']] = m
        models_l[[paste0(model)]][['line']] = preds
      }

    }#End of best model fitting loop

    #Re-arrange the presentation of the SSErr for the modes in a list
    errors_var = sapply(models_l, function(x) return(attr(x$model, 'SSErr')))
    models_lines_plot = do.call(rbind,lapply(models_l, function(x) return(x$line)))
    models_lines_plot$Model = as.factor(models_lines_plot$Model)

    #Select the minimum
    var_min = which.min(errors_var)

    if(length(var_min) == 0){
      #Create variogram directory
      vars_dir = file.path(s_dir, krig_type, 'Variograms')
      dir.create(vars_dir, showWarnings = F, recursive = T)

      #Plot variogram and the models
      p = ggplot2::ggplot() +
        ggplot2::geom_point(data = v, ggplot2::aes(x = dist, y = gamma, fill = np), shape = 21, size = 3) +
        viridis::scale_fill_viridis(direction = -1)+
        expand_limits(x = 0, y = 0)+
        ggplot2::labs(x = "Distance", y = "Semivariance", title = paste0(krig_type," Variogram: ", gsub('_', ' ', var_press)))+
        ggplot2::annotate("text",label=paste("NO FIT"), x = max(v$dist)/2, y = max(v$gamma)/3)


      #Save plots
      ggplot2::ggsave(file.path(vars_dir, paste0(var_press, '_variogram.png')), p)
    }else{
      selected_model = models_l[[var_min]]
      m = selected_model$model

      #Save the best model type to a list
      models_chosen_l[[var]][[krig_type]] = as.character(m$model[2])

      #Create variogram directory
      vars_dir = file.path(s_dir, krig_type, 'Variograms')
      dir.create(vars_dir, showWarnings = F, recursive = T)

      #Plot variogram and the models
      p = ggplot2::ggplot() +
        ggplot2::geom_point(data = v, ggplot2::aes(x = dist, y = gamma, fill = np), shape = 21, size = 3) +
        viridis::scale_fill_viridis(direction = -1)+
        expand_limits(x = 0, y = 0)+
        ggplot2::geom_line(data = models_lines_plot, ggplot2::aes(x = dist, y = gamma, colour = Model))+
        ggplot2::labs(x = "Distance", y = "Semivariance", title = paste0(krig_type," Variogram: ", gsub('_', ' ', var_press)))+
        ggplot2::annotate("text",label=paste('Chosen Model',
                                             '\nModel = ', m$model[2], '\nNugget = ', round(m$psill[1], 4), '\nSill = ',  round(m$psill[2]+m$psill[1], 4), '\nRange = ',  round(m$range[2], 2)),color="dark gray",x = 200, y=(m$psill[2]+m$psill[1])/3)


      #Save plots
      ggplot2::ggsave(file.path(vars_dir, paste0(var_press, '_variogram.png')), p)


      #########################################################################
      ##################  Loop - Local and Global Scenarios ###################
      #########################################################################

      i = 1
      for(i in 1:nrow(local_scenarios)){

        local_ = local_scenarios[i,]

        #LOOCV for neighbor value
        cv_kr = krige.cv(formula_, var_data, m, nfold=nrow(var_data), nmax = local_)

        #Perform back transformation on predictions to help with understanding and
        #save to new df
        if(NORM){
          results_df = data.frame(Predicted = predict(bc_trans, cv_kr$var1.pred,  inverse = T))
          results_df$Observed = predict(bc_trans, cv_kr$observed,  inverse = T)
          results_df$fold = cv_kr$fold
        }else{
          results_df = data.frame(Predicted = cv_kr$var1.pred)
          results_df$Observed = cv_kr$observed
          results_df$fold = cv_kr$fold
        }


        #Linear Regression for observed vs. predicted
        t = lm(Predicted~Observed, results_df)
        stats_ = summary(t)

        Slope = stats_$coefficients[2,1]
        Intercept = stats_$coefficients[1,1]

        #Test result for H0 => a=0 - if p value <0.05 we reject the H0
        IsIntercept0 = stats_$coefficients[7]>0.05

        # Compute t-student H0: slope=1. The estimation of coefficients and their s.d. are in sfit$coefficients
        tstats <- (1-stats_$coefficients[2,1])/stats_$coefficients[2,2]
        # Calculates two tailed probability
        pval<- 2 * pt(abs(tstats), df = df.residual(t), lower.tail = FALSE)
        IsSlope1 = pval>0.05

        #Save R squared
        Rsquared = stats_$r.squared
        #Save R squared significance
        IsRsqrSignificant = stats_$coefficients[8]<0.05

        #Calculate RMSE
        results_df$SErr = (results_df$Predicted -results_df$Observed)^2

        Fold_RMSE = results_df %>%
          group_by(fold) %>%
          summarise(RMSE = sqrt(mean(SErr, na.rm =T)))
        RMSE_r = mean(Fold_RMSE$RMSE)

        #Calculate MAPE
        results_df$APE = abs((results_df$Observed - results_df$Predicted)/results_df$Observed)*100

        MAPE = mean(results_df$APE)

        #Save all results to an analysis df
        analysis = data.frame(Rsqr = Rsquared,Slope = Slope, Intercept = Intercept, IsIntercept0 = IsIntercept0, IsSlope1 = IsSlope1, IsRsqrSignificant = IsRsqrSignificant,
                              RMSE = RMSE_r, MAPE = MAPE, Scenario = paste0(krig_type, '_Max_', local_))

        #Add scenario to CV Results DF
        results_df$Scenario = paste0(krig_type, '_Max_', local_)

        #Save analysis and df
        results_l[[var]][['Analysis']][[paste0(krig_type, '_Max_', local_)]] = analysis
        results_l[[var]][['CV_df']][[paste0(krig_type, '_Max_', local_)]] = results_df

        #Krigging
        K = krige(locations = var_data, formula = formula_, model = m, newdata = grd, nmax = local_)

        #Set projection
        sf::st_crs(K) = sf::st_crs(var_data)

        #Convert Stars to Raster
        fname = paste0(tempfile(), ".tif")
        stars::write_stars(K['var1.pred'], fname)
        from = fname

        K_f = raster::stack(from)
        K_f = raster::readAll(K_f)

        raster::extent(K_f) = raster::extent(frst)

        #Perform the back transformation
        if(NORM){
          raster::values(K_f) = predict(bc_trans, raster::values(K_f),  inverse = T)
        }
        if(STAND){
          raster::values(K_f) = raster::values(K_f)* attr(std_data, 'scaled:scale') + attr(std_data, 'scaled:center')
        }

        #Rename variable
        names(K_f) = var

        #Raster to Polygons
        K_pol = st_as_sf(raster::rasterToPolygons(K_f))

        #Clip to field boundary
        K_pol = st_intersection(K_pol, st_geometry(field_boundary))

        #Set map name
        map_title = paste0('Sample Results ', krig_type, ': ', gsub('_', ' ', var_press), '\nMax: ', local_)

        #Get legend names
        qq = quantile(st_drop_geometry(K_pol)[,var], c(0.20, 0.4, 0.6, 0.8))

        #Plot map
        t = tmap::tm_shape(K_pol, bbox = sf::st_bbox(sf::st_buffer(field_boundary, 100))) +
          tmap::tm_polygons(var, n = 5, style = "fixed", breaks = c(-Inf, qq[1], qq[2], qq[3], qq[4], Inf), palette = tmaptools::get_brewer_pal("RdYlGn", n = 5, plot = F), border.col = 'transparent', size = 0.06,  title = gsub('_', ' ', var_press)) +
          tmap::tm_layout(map_title, asp =1.2, frame = FALSE, legend.show = T, title.position = c("left", "top"))+
          tmap::tm_shape(field_boundary) + tmap::tm_borders(lwd = 2)+
          tmap::tm_shape(soilsamples_6) + tmap::tm_dots(var, n = 5, style = "fixed", breaks = c(-Inf, qq[1], qq[2], qq[3], qq[4], Inf), palette = tmaptools::get_brewer_pal("RdYlGn", n = 5, plot = F), border.col = 'black', size = 0.5,title = gsub('_', ' ', var_press)) +
          tmap::tm_compass(position = c("right", "top"), size = 3) +
          tmap::tm_scale_bar(width = 0.2, position = c("right", "bottom"), text.size  = 0.7)

        #Create map directory
        ok_map_dir = file.path(s_dir, krig_type, 'Maps')
        dir.create(ok_map_dir, showWarnings = F, recursive = T)
        tmap::tmap_save(t, file.path(ok_map_dir, paste0(var_press, '_', krig_type, '_neigh_', local_, '_Map.png')))
      }#Finish Loop for Scenarios
    }
  }#End of Sub Process - Universal and Ordinary
}#End of main loop

#Unnest analysis results
analysis = lapply(1:length(results_l),function(x){
  nutrient = names(results_l)[x]
  df_s =do.call(rbind, results_l[[x]]$Analysis)

  df = as.data.frame(df_s)
  df$Nutrient = nutrient
  return(df)} )
analysis = do.call(rbind, analysis)

#Unnest CV Results DF
df = lapply(1:length(results_l),function(x){
  nutrient = names(results_l)[x]
  df_s =do.call(rbind, results_l[[x]]$CV_df)

  df = as.data.frame(df_s)
  df$Nutrient = nutrient
  return(df)} )

df = do.call(rbind, df)

#Save Unnested df
write.csv(analysis, file.path(s_dir, 'All_Kriging_Scenarios_Analysis.csv'))
write.csv(df, file.path(s_dir, 'All_Kriging_Scenarios_DF.csv'))


#Unnest final models
models_chosen = lapply(1:length(models_chosen_l),function(x){
  nutrient = names(models_chosen_l)[x]
  df_s = t(do.call(cbind, models_chosen_l[[x]]))

  df = as.data.frame(df_s)
  names(df) = nutrient
  return(df)} )

models_chosen = do.call(cbind, models_chosen)

write.csv(models_chosen, file.path(s_dir, 'Best_Models.csv'))

lambdas = do.call(rbind, lambda_bc_l)

