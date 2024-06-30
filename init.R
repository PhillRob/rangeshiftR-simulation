###################### Libraries#####################  
## Install RangeShiftR from GitHub:
# devtools::install_github("RangeShifter/RangeShiftR-package", ref="main")
library(devtools)
library(RangeShiftR)
# devtools::install_github("PhillRob/PopulationGrowthR", ref = "main")
# library(PopulationGrowthR)
library(sf)
library(PROJ)
library(terra)
library(raster)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(here)
# install.packages("googledrive")
library(googledrive)
# install.packages("googlesheets4")
library(googlesheets4)
library(dplyr)
library(ggplot2)
library(config)
library(readr)
# system("brew install openssl")
# install.packages("curl")
library(arrow)
library(disk.frame)

library(future.apply)
library(purrr)
library(tools)

config = config::get()
# source(here("sampling.R"))

## relative path from working directory:
dirpath = here()
dir.create("Inputs", showWarnings = TRUE)
dir.create("Outputs", showWarnings = TRUE)
dir.create("Output_Maps", showWarnings = TRUE)

# authenticate google drive
auth_g_drive = function(email) {
  drive_auth(email = email)
  gs4_auth(email = email)
}
# drive_auth(
#   email = gargle::gargle_oauth_email(),
#   path = NULL,
#   scopes = "https://www.googleapis.com/auth/drive",
#   cache = gargle::gargle_oauth_cache(),
#   use_oob = gargle::gargle_oob_default(),
#   token = NULL
# )
# ____________________________

###################### convert to Raster#####################  
csv_to_raster = function(map_data, output_name, nrows, ncols, map_extent, climate_map = TRUE) {
  # convert to raster
  target_crs <-"+proj=lcc +lat_0=0 +lon_0=134 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
  colnames(map_data)[1:2] = c("lon", "lat") # set coordinates to "lon", "lat"
  coordinates(map_data) = ~ lon + lat
  proj4string(map_data)=CRS(target_crs)

  e <- map_extent # crop boundary to Australia
  r <- raster(e, ncols = ncols, nrows = nrows)

  # rasterize and average suitability probability when there are multiple points per cell
  map_data_a <- rasterize(map_data, r, fun = mean)
  map_data2a = raster::resample(map_data_a[[2]],
                                raster(ext= e, resolution = 900, crs=projection(map_data)))

  # save as Raster in "inputs/" folder
  raster::writeRaster(round(map_data2a, 4), format="ascii", filename = here("Inputs/", output_name),
                      NAflag = -9, overwrite=FALSE, bylayer = T, datatype = ifelse(climate_map, "FLT4S", "INT4S"))
}
# _____________________________________

# OUTPUT_cont.tif from drive
drive_files = function(folder) {
  files = drive_ls(folder)
  return(files)
}
maxent_suitability = drive_files(folder = "MaxENT")
if(!"OUTPUT_cont.asc" %in% list.files(here("Inputs")) | !"OUTPUT_uniform.asc" %in% list.files(here("Inputs"))) {
  
  if(!"OUTPUT_cont.tif" %in% maxent_suitability$name) {
    paste("'OUTPUT_cont.tif' file not present in 'drive folder = MaxENT' ") |> print()
  } else {
    paste("Downloading suitability 'OUTPUT_cont.tif' from 'drive folder = MaxENT' ") |> print()
    drive_download(
      file = maxent_suitability[maxent_suitability$name == "OUTPUT_cont.tif", ],
      path = here("Inputs", "OUTPUT_cont.tif"))
  }
  
  map_data = raster(here("Inputs", "OUTPUT_cont.tif")) 
  
  # map_resolution = res(map_data)
  # ncells = list(nrows = nrow(map_data), ncols = ncol(map_data))
  target_crs <-"+proj=lcc +lat_0=0 +lon_0=134 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
  
  ncells <- list(nrows = nrow(map_data), ncols = ncol(map_data))
  map_extent <- extent(map_data)
  map_resolution <- res(map_data)
  
  res <- 900  # desired pixel size in meters
  
  projected_raster <- projectRaster(map_data, crs = target_crs, res=c(res, res))
  res(projected_raster)
  
  if(!"OUTPUT_cont.asc" %in% list.files(here("Inputs"))) {
    projected_raster |>
      round(2)|>
      # calc(fun=function(x){x * 100}) |>
      writeRaster(here("Inputs","OUTPUT_cont.asc"), format="ascii",overwrite=T,datatype="FLT4S",NAflag = -9999)
  }
  
  
  if(!"OUTPUT_uniform.asc" %in% list.files(here("Inputs"))) {
    projected_raster |>
      calc(fun=function(x){x /x*0.95})|>
      writeRaster(here("Inputs","OUTPUT_uniform.asc"), format="ascii",overwrite=T,datatype="FLT4S",NAflag = -9999)
  }
  
  # file.remove(here("Inputs", "OUTPUT_cont.tif"))
}

# 
# 
# writeRaster(x,here("Inputs","OUTPUT_uniform.asc"), format="ascii",overwrite=T,datatype="INT2S",NAflag = -9999)

# _____________________


# upload output folder to drive
output_folder_upload = function(index, base_climate_map) {
  auth_g_drive(email = config$GMAIL_STORE2)
  
  if(!paste0("simulated-populations-", base_climate_map, "-climate") %in% drive_ls("Simulation-Results")$name) {
    drive_mkdir(paste0("simulated-populations-", base_climate_map, "-climate"))
    drive_mv(file = paste0("simulated-populations-", base_climate_map, "-climate"), path = "Simulation-Results/")
  }
  
  if(!paste0("output_", index) %in% drive_ls(paste0("simulated-populations-", base_climate_map, "-climate"))$name) {
    drive_mkdir(paste0("output_", index), path = paste0("simulated-populations-", base_climate_map, "-climate/"))
  }
  
  for (media in list.files(here("Outputs"))) {
    # convert output file to rds
    if(endsWith(media, ".csv")) {
      read_csv(here("Outputs", media)) |>
        saveRDS(file = here("Outputs", paste0(media, ".rds")))
    } else if(endsWith(media, ".txt")){
      read_tsv(here("Outputs", media)) |>
        saveRDS(file = here("Outputs", paste0(media, ".rds")))
    }
    
    # delete local csv file
    file.remove(here("Outputs", media))
  }
  
  for(media in list.files(here("Outputs"), pattern = ".rds")) {
    if(!media %in% drive_ls(path = paste0("simulated-populations-", base_climate_map, "-climate/output_", index, "/"))$name) {
      drive_upload(media = here("Outputs", media), path = paste0("simulated-populations-", base_climate_map, "-climate/output_", index, "/"))
      # delete local csv file
      file.remove(here("Outputs", media))
    }
  }
}

# ____________________

# update lag analysis output
update_lags = function(lag_table, base_climate_map) {
  auth_g_drive(email = config$GMAIL_STORE2)
  
  if(!paste0("simulated-populations-", base_climate_map, "-climate") %in% drive_ls("Simulation-Results")$name) {
    drive_mkdir(paste0("simulated-populations-", base_climate_map, "-climate"))
    drive_mv(file = paste0("simulated-populations-", base_climate_map, "-climate"), path = "Simulation-Results/")
  }
  
  files = drive_ls(paste0("simulated-populations-", base_climate_map, "-climate"))
  lag_table_present = paste0("lag_table_", base_climate_map, "_climate") %in% files$name
  
  if(lag_table_present) {
    file_id = files[which(files$name == paste0("lag_table_", base_climate_map, "_climate")), "id"]
  }
  
  if (!lag_table_present) {
    file = gs4_create(name = paste0("lag_table_", base_climate_map, "_climate"), sheets = lag_table)
    drive_mv(file = file, path = paste0("simulated-populations-", base_climate_map, "-climate/"))
  } else {
    sheet_append(ss = file_id$id, data = lag_table, sheet = "lag_table")
  }
}


# _________________


# upload lag plots to drive
lag_plots_folder_upload = function(base_climate_map) {
  if(!paste0("simulated-populations-", base_climate_map, "-climate") %in% drive_ls("Simulation-Results")$name) {
    drive_mkdir(paste0("simulated-populations-", base_climate_map, "-climate"))
    drive_mv(file = paste0("simulated-populations-", base_climate_map, "-climate"), path = "Simulation-Results/")
  }
  
  if(!paste0("lag_plots") %in% drive_ls(paste0("simulated-populations-", base_climate_map, "-climate"))$name) {
    drive_mkdir(paste0("lag_plots"), path = paste0("simulated-populations-", base_climate_map, "-climate/"))
  }
  
  uploaded_plots  = drive_ls(path = paste0("simulated-populations-", base_climate_map, "-climate/lag_plots/"))
  
  for(media in list.files(here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots")), pattern = ".png")) {
    if(!media %in% uploaded_plots$name) {
      drive_upload(media = here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots"), media), path = paste0("simulated-populations-", base_climate_map, "-climate/lag_plots/"))
      # delete local csv file
      file.remove(here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots"), media))
    }
  }
}


# ________________________


################ plot populations ################
plot_pops = function(id, range_df, pop_df, base_climate_map, lagresults,  n_plots = 10) {
  
  ## uniform Population
  all_reps = (range_df |> dplyr::select(Rep) |> unique() |> collect())$Rep
  lag_lengths = lagresults$Laglength |> na.omit()
  lag_reps = lagresults[which(!is.na(lagresults$Laglength)), "Rep"]
  non_lag_reps = all_reps[which(! all_reps %in% lag_reps)]
  
  if(length(lag_reps) > n_plots) {
    lag_reps = lag_reps[1:n_plots]
  }
  
  if(length(non_lag_reps) > n_plots) {
    non_lag_reps = non_lag_reps[1:n_plots]
  }
  
  # non_lag_reps = ifelse((non_lag_reps |> length()) > 50, non_lag_reps[1:50], non_lag_reps) # if more tha  50 non lag replicates, plot only first 50
  
  # number of occupied cells
  occupancy_data = open_tsv_dataset(here("Outputs", "Batch1_Sim1_Land1_Occupancy_Stats.txt")) |>
    collect() |>
    ggplot(aes(x = Year, y = Mean_OccupSuit)) + geom_line() + 
    labs(title = paste("MEAN RATIO BETWEEN OCCUPIED AND SUITABLE CELLS\n IN SIMULATION:", id), 
         subtitle = paste0("(", base_climate_map |> tools::toTitleCase(), " CLIMATE) INTRODUCTION CELL SUITABILITY: ", lagresults$intro_suitability |> unique())) +
    theme_linedraw() 
  
  if(!paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots") %in% list.files(here())) {
    dir.create(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots"))
  }
  
  occupancy_data |> ggsave(filename = paste0("simulation_", id, "_Ratio of occupied to suitable cells.png"), 
                           path = here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots")))
  
  
  # plots
  paste("lag maps") |> print()
  
  if(sum(!is.na(lag_reps)) > 0) {
    lag_plots <-  1:length(lag_reps) |>
      map(~ pop_df |>  
            filter(Rep == lag_reps[.x]) |>
            group_by(Year) |>
            summarise(NInd = sum(NInd)) |>
            ungroup() |> collect() |>
            ggplot(aes(x = Year, y = NInd)) +
            geom_line() +
            geom_vline(xintercept = lag_lengths[.x], col = "red") +
            theme_light() +
            theme(legend.position = "none") +
            labs(title = paste("SIMULATION:", id, "REPLICATE:", lag_reps[.x], "(", base_climate_map, "CLIMATE)",
                               "\n INTRODUCTION CELL SUITABILITY:", lagresults$intro_suitability),
                 subtitle = paste("SIMULATED RESIDENCY:", lagresults$Year, "ACTUAL RESIDENCY:",
                                  max((pop_df |> dplyr::filter(Rep == lag_reps[.x]) |> dplyr::select(Year) |> collect())$Year),
                                  "LAG:", TRUE))
      )
    
    # save plots
    1:length(lag_plots) |> purrr::map(~ggsave(filename = paste0("simulation_", id, "_replicate_", lag_reps[.x], ".png"),
                                              plot = lag_plots[[.x]],
                                              path = here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots"))))
    rm(list = "lag_plots")
  }
  
  
  
  paste("non - lag maps") |> print()
  
  if(sum(!is.na(non_lag_reps)) > 0) {
    non_lag_plots = 1:length(non_lag_reps) |>
      purrr::map(~ pop_df |> filter(Rep %in% non_lag_reps[.x]) |> group_by(Year) |> collect() |> 
                   summarise(NInd = sum(NInd)) |>
                   ungroup() |> collect() |> ggplot(aes(x = Year, y = NInd)) + geom_line() +
                   theme_light() + theme(legend.position = "none") + 
                   labs(title = paste("SIMULATION:", id, "REPLICATE:", non_lag_reps[.x], "(", base_climate_map, "CLIMATE)",
                                      "\n INTRODUCTION CELL SUITABILITY:", lagresults$intro_suitability), 
                        subtitle = paste("SIMULATED RESIDENCY:", lagresults$Year, "ACTUAL RESIDENCY:",
                                         max((pop_df |> filter(Rep == non_lag_reps[.x]) |> dplyr::select(Year) |> collect())$Year),
                                         "LAG:", FALSE)))
    
    # save plots
    1:length(non_lag_plots) |> purrr::map(~ggsave(filename = paste0("simulation_", id, "_replicate_", non_lag_reps[.x], ".png"),
                                                  plot = non_lag_plots[[.x]],
                                                  path = here(paste0(base_climate_map |> toTitleCase(), "-Simulated-Lag-Plots"))))
    
    rm(list = "non_lag_plots")
  }
  
  lag_plots_folder_upload(base_climate_map = base_climate_map)
}

# _____________________


######################   lag assessment#####################  
lag_assessment = function(parameters_table, id, base_climate_map) {
  
  # range
  # range_df = readRange(simulate, dirpath)
  range_df = open_tsv_dataset(here("Outputs", "Batch1_Sim1_Land1_Range.txt"))
  
  actual_residency = range_df |>
    group_by(Rep) |> summarise(actual_residency = max(Year)) |>
    rename(Species = Rep) |> collect()
    

  
  range_1000_years = range_df |> dplyr::select(Rep) |> collect() |> table() |> 
    data.frame() |> 
    rename(Replicate = Rep, Max_Years = Freq) |> 
    mutate(Max_Years = Max_Years -1)
  range_1000_years2 = range_1000_years[order(range_1000_years$Max_Years, decreasing = T), ]
  
  # lag detection
  # devtools::install_github("PhillRob/PopulationGrowthR", force = F)
  library(PopulationGrowthR)
  
  # data prep
  range_df2 <- range_df[, c(1, 2, 4)] |>
    rename(Species = Rep, Year = Year, Frequency = NInds)
  freYear <-
    aggregate(Frequency ~ Year, range_df2, function(x)
      cumsum(x))
  
  fdata1 <- cbind(range_df2 |> collect(), as.vector(unlist(freYear[2])))
  colnames(fdata1) <- c("Species", "Year", "Frequency", "Specimens")
  yeardata1 <- aggregate(Frequency ~ Year, fdata1, function(x)
    sum(x))
  
  colnames(yeardata1)<-  c("Year", "Specimens")
  
  # detection
  fitall = lagfit(data = fdata1, yeardata = yeardata1)
  lagresults<-fitall[["fitdata"]] |> as.data.frame()
  
  variables = c("index", "Rep", "intro_location", "intro_suitability",
                "carrying_capacity", "max_offspring", 
                "competition_coefficient", "short_range_distance",
                "long_range_distance", "short_range_probability",
                "emigration_probability", "no_introduced_per_cell",
                "Year", "Species", "Laglength", "actual_residency", "Lagspecies", 
                "% Lag Species", "% Lag Length of Years Assessed")
  
  if(!is.null(lagresults)) {
    # both scenes starting with "0+", "00" and "0-" slopes are identified as lag events. We are interested in lags before exponential growth i.e "0+"
    lagresults <- lagresults %>% 
      mutate(lags = Lag,
             Species = as.numeric(Species)) |>
      left_join(actual_residency) |>
      dplyr::select(Species, Laglength, actual_residency)  
    
    lag_species = lagresults |>
      dplyr::filter(!is.na(Laglength)) |> nrow()
    
    if(nrow(lagresults) > 0) {
      lagresults =(lagresults |> rename(Rep = Species) |> mutate(parameters_table)) |>
        mutate(Lagspecies = lag_species, 
               `% Lag Species` = (100 * Lagspecies / Species), 
               `% Lag Length of Years Assessed` = (100 * Laglength / Year))
      
      
      lagresults = lagresults |> dplyr::select(all_of(variables))
      # save lag output data to DB
      update_lags(lag_table = lagresults, base_climate_map = base_climate_map)
      
    } else {
      lagresults = parameters_table |> mutate(Laglist = NA, Lagspecies = NA, Laglength = NA, `% Lag Species` = NA, `% Lag Length of Years Assessed` = NA)
      lagresults = lagresults |> dplyr::select(all_of(variables))
      
      # save lag output data to DB
      update_lags(lag_table = lagresults, base_climate_map = base_climate_map)
    }
  }
  
  # population data
  pop_df = open_tsv_dataset(here("Outputs", "Batch1_Sim1_Land1_Pop.txt"))
  
  # plot populations:
  plot_pops(id = parameters_table$index, range_df = range_df, pop_df = pop_df, lagresults = lagresults, base_climate_map = base_climate_map, n_plots = 10)
  gc()
  
  # merge range data
  update_pop_data(pop_df = pop_df, id = parameters_table$index, base_climate_map = base_climate_map)
  gc()
  
  # upload raw files to drive 
  output_folder_upload(index = parameters_table$index, base_climate_map = base_climate_map)
}


#################### new introduction points and their suitability scores ####################
introduction = rbind(
  tibble(lon = 1565728, lat = -3949947, locations = "Sydney"),
  tibble(lon = 1884177, lat = -3285164, locations = "Brisbane"),
  tibble(lon = -1694855, lat = -3751189, locations = "Perth"),
  # tibble(lon = -2742416, lat = 1695296, locations = "Yeppon"),
  tibble(lon = 988155, lat = -4324335, locations = "Melbourne")
) |> mutate(Individuals = 1) |>
  dplyr::select(lon, lat, Individuals, locations)

introduction$suitability = terra::extract(raster(here("Inputs", "OUTPUT_cont.asc")),
                                   introduction[, c("lon", "lat")])
# introduction$suitability
# __________________________
# create raster from introduction points
for(location in introduction$locations) {
  projected_raster = raster(here("Inputs", "OUTPUT_cont.asc")) 
  if(!paste0(location, ".asc") %in% list.files(here("Inputs"))) {
    data = introduction |> 
      dplyr::filter(locations == location) |>
      dplyr::select(lon, lat, Individuals) |>
      csv_to_raster(
        output_name = location,
        nrows = nrow(projected_raster),
        ncols = ncol(projected_raster),
        map_extent = extent(projected_raster),
        climate_map = FALSE)
    
  }
  
}

# ______________________

#####################  update populations data for sampling #####################  
update_pop_data = function(pop_df, id, base_climate_map) {
  paste("Joining pop tables:", id) |> print()
  
  pop_df = pop_df |> collect() |>
    dplyr::mutate(Simulation = id) |>
    mutate(Simulation = as.integer(Simulation)) |>
    arrow::as_arrow_table()
  
  lst = drive_ls(paste0("simulated-populations-", base_climate_map, "-climate"))
  
  if(!paste0("joined_", base_climate_map, "_pop_data.parquet") %in% list.files(here("WWW"))) {
    if(paste0("joined_", base_climate_map, "_pop_data.parquet") %in% lst$name) {
      
      lst = lst[lst$name == paste0("joined_", base_climate_map, "_pop_data.parquet"),]
      drive_download(lst, path = here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
      pop_df = read_parquet(here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")), as_data_frame = FALSE) |>
        full_join(pop_df) 
      
      file.remove(here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
      write_parquet(pop_df |> as_arrow_table(), sink = here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
      
    } else {
      write_parquet(pop_df |> as_arrow_table(), sink = here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
    }
  } else {
    pop_df = read_parquet(here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")), as_data_frame = FALSE) |>
      full_join(pop_df) 
    
    file.remove(here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
    write_parquet(pop_df |> as_arrow_table(), sink = here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")))
  }
  
  drive_upload(here("WWW", paste0("joined_", base_climate_map, "_pop_data.parquet")), 
               path = paste0("simulated-populations-", base_climate_map, "-climate/"),
               overwrite = TRUE)
  
  
  Trash = drive_find(trashed = TRUE)
  drive_rm(Trash)
}

# _________________________


# Previous parameters and lag values
if(! "old_lag_results.csv" %in% list.files(here("Inputs"))) {
  drive_download(file = "MaxENT/old_lag_results.csv", path = here("Inputs", "old_lag_results.csv"))
  prior_values = read.csv(here("Inputs", "old_lag_results.csv"))
}

# seq(50, 500, by = 50)
#################### initialization parameters #########################
# set.seed(1300)
# intro_point = rep(introduction$locations[2],10)#sample(introduction$locations, 10, replace = T) #sample(introduction$locations, 100, replace = TRUE))
# competition = rep(0.05, 10) #runif(n = 100, min = 0.049, max = 0.051) # nrow(prior_values)
# short_range_d = rep(1100, 10)#runif(n = 100, min = 910, max = 1100) |> round(0)
# long_range_d = rep(9000, 10)#runif(n = 100, min = 1000, max = 10000) |> round(0)
# introductions = rep(50,10)#runif(n = 10, min = 50, max = 500) |> round(0)#seq(50, 500, by = 50)#rep(5,4)#runif(n = 100, min = 20, max = 100) |> round(0)
# emigration = rep(0.1, 10)#runif(n = 100, min = 0.2, max = 0.4)
# carrying_capacity = rep(10000, 10)#runif(n = 100, min = 5000, max = 10000) |> round(0)
# max_offspring = rep(50, 10)# runif(n = 10, min = 50, max = 500) |> round(0)# seq(50, 500, by = 50)# rep(50,4)#runif(n = 100, min = 5, max = 50) |> round(0)
# short_range_probability = seq(0.05, 0.95, by = .1) #rep(0.9, 10)# runif(n = 100, min = 0.7, max = 0.9)
# Year = 50 # runif(n = 1, min = 100, max = 500) |> round(0)
# replicates = 5 # runif(n = 1, min = 100, max = 500) |> round(0)


set.seed(1300)
intro_point = rep(introduction$locations[4],100)#sample(introduction$locations, 10, replace = T) #sample(introduction$locations, 100, replace = TRUE))
competition = rep(0.05, 100) #runif(n = 100, min = 0.049, max = 0.051) # nrow(prior_values)
short_range_d = runif(n = 100, min = 910, max = 1100) |> round(0)
long_range_d = runif(n = 100, min = 9100, max = 11000) |> round(0)
introductions = runif(n = 100, min = 5, max = 100) |> round(0)#seq(50, 500, by = 50)#rep(5,4)#runif(n = 100, min = 20, max = 100) |> round(0)
emigration = runif(n = 100, min = 0.05, max = 0.5)
carrying_capacity = rep(10000, 100)#runif(n = 100, min = 5000, max = 10000) |> round(0)
max_offspring = runif(n = 100, min = 50, max = 500) |> round(0)# seq(50, 500, by = 50)# rep(50,4)#runif(n = 100, min = 5, max = 50) |> round(0)
short_range_probability = runif(n = 100, min = 0.5, max = 0.9)
Year = 50 # runif(n = 1, min = 100, max = 500) |> round(0)
replicates = 10 # runif(n = 1, min = 100, max = 500) |> round(0)


parameters <- function (index) {
  set.seed(NULL)
  tibble::tibble(index = index,
                 intro_location = intro_point[index],
                 intro_suitability = introduction[introduction$locations == intro_point[index], "suitability"]$suitability,
                 carrying_capacity = carrying_capacity[index], 
                 max_offspring = max_offspring[index], # prior_values$max_offspring[index]
                 competition_coefficient = competition[index],
                 short_range_distance = short_range_d[index], # prior_values$short_range_distance[index]
                 long_range_distance = long_range_d[index], # prior_values$short_range_distance[index]
                 short_range_probability = short_range_probability[index],
                 emigration_probability = emigration[index], # sample(seq(0.005, 0.05, by = 0.001), 1)
                 no_introduced_per_cell = introductions[index], 
                 Year = Year,
                 Species = replicates)
}


