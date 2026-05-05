# define some functions

# Compute the "rank-normalized R-hat" by Vehtari et al. (2021) from jagsUI
# output.
# Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667–718.
# https://doi.org/10.1214/20-BA1221
# 
# The first input is MCMC samples. If the jagsUI output is jm, this is jm$samples.
# The second input is a string of regular expression. This is a bit
# complicated. For example, to select "BF.Fixed" and all K parameters, which 
# are indexed, Use "^BF\\.Fixed|^K\\[" A '^' specifies that the following letter
# is the beginning of a string. '\\.' specifies a literal period, which needs
# to be "escaped" by two backslashes (\\). A square bracket needs to be escaped
# with two backslashes as well. The pipe (|) indicates 'or'.  
# 

rank.normalized.R.hat <- function(samples, params, MCMC.params){
  library(posterior)
  library(coda)

  col.names <- grep(params, varnames(samples), value = TRUE, perl = TRUE)
  subset.mcmc.samples <- samples[, col.names]

  subset.mcmc.array <- as_draws_array(subset.mcmc.samples, .nchains = MCMC.params$n.chains)

  rhat.values <- apply(subset.mcmc.array, 
                       MARGIN = 3, 
                       FUN = posterior::rhat)
  
  return(rhat.values)
}

# This is for the rank-normalized R-hat (Vehtari et al. 2021)
high.Rhat <- function(x, start.year){
  return(data.frame(idx = which(x > 1.01),
                    start.year = start.year[which(x > 1.01)],
                    Rhat = x[which(x > 1.01)]))
}

# Create an observer list with new data:
# The input is only for "new" data that are not part of Laake's dataset.
create.observer.list <- function(sightings){
  years <- 1 + unique(sightings$Start.year)
  Observer <- read.csv(file = "Data/Observer_Laake.csv")
  
  # Find Observers in Laake's observer list who are also in the new observer list
  # There are multiple initials per person in some cases. ID numbers should be the
  # same for these initials
  Observer %>%
    filter(is.na(Observer) & !is.na(Initials)) %>% #-> tmp
    droplevels() %>%
    mutate(Identifier = as.factor(Initials)) %>%
    select(ID, Identifier, Name) -> Observer.NA
  
  Observer %>%
    filter(is.na(Initials) & !is.na(Observer)) %>% #-> tmp #
    droplevels() %>%
    mutate(Identifier = as.factor(Observer)) %>%
    select(ID, Identifier, Name) -> Observer.Initial.NA
  
  Observer %>%
    filter(!is.na(Observer) & !is.na(Initials)) %>%
    droplevels() %>%
    mutate(Identifier = as.factor(Initials)) %>%
    select(ID, Identifier, Name) -> Observer.not.NA
  
  Observer %>%
    filter(!is.na(Observer) & !is.na(Initials)) %>%
    droplevels() %>%
    mutate(Identifier = as.factor(Observer)) %>%
    select(ID, Identifier, Name) -> Observer.not.NA.2
  
  Observer.1 <- rbind(Observer.not.NA, 
                      Observer.NA, 
                      Observer.Initial.NA,
                      Observer.not.NA.2) 
  
  uniq.Observer.1 <- data.frame(Name = unique(Observer.1$Name),
                                ID.new = seq(1:length(unique(Observer.1$Name))))
  
  uniq.Observer.1  %>%
    left_join(Observer.1, by = "Name") %>% #-> tmp # 
    rename(ID.old = ID) %>%
    select(-Name) %>%
    rename(ID = ID.new) %>%
    mutate(data = "pre2009") -> Observer.2
  
  # Find observers from more recent data (2009/2010 to present)
  Observers.new <- unique(sightings$Observer)
  new.observers <- data.frame(ID.old = seq(1, length(Observers.new)),
                              Identifier = Observers.new,
                              ID = seq((max(Observer.2$ID)+1), 
                                       (max(Observer.2$ID) + length(Observers.new))),
                              data = "new.data")
  
  # Combine "old" and "new" observer lists
  # Find observers who showed up in both datasets.
  # Keep the original observer ID. 
  Observer.2 %>%
    left_join(new.observers, by = "Identifier") %>% #-> tmp #
    filter(!is.na(ID.y)) %>% #-> tmp#
    transmute(ID = ID.x,
              ID.new = ID.y,
              Identifier = Identifier,
              data = data.x) -> tmp
  
  # ID.new needs to be replaced with the "new" Laake-based ID (ID)
  for (k in 1:nrow(tmp)){
    new.observers[new.observers$ID == tmp$ID.new[k], "ID"] <- tmp$ID[k]
    
  }

  tmp.4 <- rbind(Observer.2, new.observers) %>%
    mutate(Observer = Identifier) %>%
    select(-Identifier) %>%
    na.omit() 
  
  # Fix one observer because "ARV/AVS" is not usable 
  ARV.ID <- tmp.4[tmp.4$Observer == "ARV/AVS", "ID"]
  ARV.ID.old <- tmp.4[tmp.4$Observer == "ARV/AVS", "ID.old"]
  tmp.4 %>%
    filter(Observer != "ARV/AVS") %>%  
    droplevels()  -> tmp.5
  
  # add those initials back in with the same ID
  all.observers <- rbind(tmp.5, 
                         data.frame(ID.old = c(ARV.ID.old, ARV.ID.old),
                                    ID = c(ARV.ID, ARV.ID), 
                                    Observer = c("ARV", "AVS"),
                                    data = "pre2009")) %>%
    mutate(obs = Observer) 
    
  all.observers %>%
    distinct(obs, .keep_all = T) -> uniq.observers 
  
  observers.list <- list(all = all.observers,
                         unique = uniq.observers)
  
  return(observers.list)
}

# Compute new Rhat statistic of Vehtari et al. (2021)
new.R.hat <- function(jm, params){
  library(posterior)
  samples <- jm$sims.list[params]
  draws <- as_draws_array(samples)
  new.r.hat <- rhat(draws)
}

# retrieve BUGS results
get.results.BUGS <- function(BUGS.file.name){
  out <- readRDS(paste0("RData/", BUGS.file.name))
  out$BUGS.out$summary %>%
    as.data.frame() %>%
    rownames_to_column(var = "parameter") %>%
    filter(grepl("Corrected", parameter)) -> summary.Nhat
  
  file.name.parts <- strsplit(BUGS.file.name, "_") %>% unlist()
  watch.dur <- strsplit(file.name.parts[4], "min") %>% unlist() %>% as.numeric()
  
  WinBUGS.Nhats.df <- data.frame(start.year = out$BUGS.input$all.years - 1,
                                 Mean = summary.Nhat$mean,
                                 SE = summary.Nhat$sd,
                                 LCL = summary.Nhat$`2.5%`,
                                 UCL = summary.Nhat$`97.5%`,
                                 model = "BUGS",
                                 min.watch = watch.dur[2],
                                 Method = "Durban",
                                 data.set = paste0("2007to", end.year))
  return(WinBUGS.Nhats.df) 
  
}


# retrieve jags results.
get.results.jags <- function(file.name){
  
  out <- readRDS(paste0("RData/", file.name))
  max.Rhats <- lapply(out$jm$Rhat, max, na.rm = T)
  
  # To compute LOOIC, need to turn zeros into NAs when there were no second station:
  jags.data <- out$jags.input$jags.data
  data.array <- jags.data$n
  data.array[,2,which(jags.data$n.station == 1)] <- NA
  
  LOOIC.n <- compute.LOOIC(loglik.array = out$jm$sims.list$log.lkhd,
                           data.array = data.array,
                           MCMC.params = out$MCMC.params)
  
  bad.Pareto <- length(which(LOOIC.n$loo.out$diagnostics$pareto_k>0.7))/length(LOOIC.n$loo.out$diagnostics$pareto_k)
  max.Pareto <- max(LOOIC.n$loo.out$diagnostics$pareto_k)
  
  # model version
  filename.parts <- unlist(strsplit(out$jags.model, "_"))
  model.v <- filename.parts[length(filename.parts)] %>% strsplit(".txt") %>% unlist()
  
  # watch duration - I forgot to extract minimum watch duration for Laake data
  # analysis... It has been fixed but rather than re-running the code, I extract
  # this information from the file name
  watch.dur <- out$jags.input$min.dur
  if (is.null(watch.dur)){
    file.name.parts <- unlist(strsplit(file.name, "min"))
    watch.dur <- unlist(strsplit(file.name.parts[2], "_"))[1] %>% as.numeric()
  }
  
  # Find out estimates - these are different among which datasets were used
  if (stringr::str_detect(file.name, paste0("min", watch.dur, "_NoBUGS_"))){
    all.start.years <- c(out$jags.input$jags.input.Laake$all.start.year,
                         out$jags.input$jags.input.new$start.years)
    data.set <- "_NoBUGS_"
    
  } else if (stringr::str_detect(file.name, "AllYears")){
    all.start.years <- out$jags.input$start.years  
    data.set <- "AllYears"
  } else if (stringr::str_detect(file.name, "Since2006")){
    all.start.years <- out$jags.input$start.years
    data.set <- "Since2006"
  } else if (stringr::str_detect(file.name, "LaakeData")){
    all.start.years <- out$jags.input$all.start.year
    data.set <- "LaakeData"
  } else if (stringr::str_detect(file.name, "_Since2010_NoBUGS")){
    all.start.years <- out$jags.input$start.years
    data.set <- "Since2010_NoBUGS"
  }
  
  model.fit <- data.frame(model = model.v,
                          data.set = data.set,
                          watch.dur = watch.dur,
                          bad.Pareto = bad.Pareto,
                          max.Pareto = max.Pareto,
                          max.Rhats = max(unlist(max.Rhats), na.rm = T))
  
  Nhats.df <- data.frame(start.year = all.start.years,
                         Mean = out$jm$mean$Corrected.Est,
                         SE = out$jm$sd$Corrected.Est,
                         LCL = out$jm$q2.5$Corrected.Est,
                         UCL = out$jm$q97.5$Corrected.Est,
                         model = model.v,
                         min.watch = watch.dur,
                         data.set = data.set)
  
  return(out.list <- list(model.fit = model.fit,
                          Nhats = Nhats.df))
  
}


# Runs a Richards' function with Pois-Binom model on datasets without using
# WinBUGS input, i.e., creating input data from output of Extract_Data_All_v2.Rmd.
# It also uses Laake's data.
# 
# Requires to provide the minimum watch duration, a version of the model (see 
# Jags models for details and differences - they are
# different in which parameters of Richards' function are time specific),
# years for which the model is fit, a directory name where data are stored,
# parameters to be monitored, and MCMC parameters as a list.
# 
# max.day = the maximum number of days to be included. It is counted from 1 December.
#           Default is 100 because there were some observations beyond day 90 in old 
#           (Laake's) dataset
# obs.n.min = the number of observation records per observer to be included. Observers with
#             less than obs.n.min records are pooled as "others" for calculating the 
#             random effects of observers. Default is 10.
# Results are saved in a subdirectory under the current working directory named "RData."
# 
# Example:
# jags.params <- c("OBS.RF", "BF.Fixed",
#                   "VS.Fixed",
#                   "mean.prob", "mean.N", "Max",
#                   "Corrected.Est", "Raw.Est", "N",
#                   "K", "S1", "S2", "P",
#                   "Max.alpha", "Max.beta",
#                   "S1.alpha", "S2.alpha",
#                   "S1.beta", "S2.beta",
#                   "P.alpha", "P.beta",
#                   "K.alpha", "K.beta",
#                   "N.alpha",
#                   "log.lkhd")
# 
# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)
#                     
# NoBUGS_Richards_fcn(min.dur = 30, 
#                     ver = "v3", 
#                     years = c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024), 
#                     data.dir = "RData/V2.1_Nov2024", 
#                     jags.params = jags.params, 
#                     MCMC.params = MCMC.params,
#                     max.day = 100,
#                     obs.n.min = 10,
#                     N.obs = 10)     
#                     
NoBUGS_Richards_fcn <- function(min.dur, years, data.dir, jags.params, MCMC.params, max.day = 100, obs.n.min = 10, N.obs = 10,Run.date = Sys.Date(), model.name, ext = ".jags"){
  
  # N.obs is the number of "top" observers who sighted the most whales among
  # all observers. 
  #Run.date <- Sys.Date()
  #model.name <- paste0(model.name.root, ver) 
  print(paste0("Starting NoBUGS_Richards_fcn at ", Sys.time(), " for Model: ", model.name))
  
  jags.model <- paste0("models/", model.name)
  
  model.name.part.1 <- strsplit(model.name, split = ext)[[1]]
  model.name.part <- strsplit(model.name.part.1, split = "model_")[[1]][2]
  
  out.file.name <- paste0("RData/JAGS_", model.name.part, 
                          "_1968to", max(years), 
                          "_min", min.dur,
                          "_NoBUGS.rds")
  
  if (!file.exists(out.file.name)){
    jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir, max.day)                        
    #jags.input.list$jags.data["N"] <- NULL
    # Modify jags data to rearrange days and provide zeros for t = 1 and t = max.day
    jags.data <- jags.input.list$jags.data
    
    # Code observers so that only observers with the minimum sample size are kept
    obs.vec <- as.vector(jags.data$obs) 
    data.frame(obs = obs.vec) %>%
      mutate(obs.f = as.factor(obs)) %>%
      group_by(obs.f) %>%
      summarize(n = n(),
                obs = first(obs)) -> obs.summary
    
    obs.too.few <- obs.summary %>% filter(n < obs.n.min)  
    
    obs.to.keep <- obs.summary %>% filter(n >= obs.n.min) 
    obs.to.keep$new.ID <- seq(1, dim(obs.to.keep)[1])
    obs.others <- max(obs.to.keep$new.ID)
    
    obs <- jags.data$obs
    new.no.obs <- obs.others + 1
    old.no.obs <- max(obs.to.keep$obs)
    for (k in 1:nrow(obs.too.few)){
      obs[obs == obs.too.few$obs[k]] <- NA
    }
    
    for (k in 1:(nrow(obs.to.keep)-1)){
      obs[obs == obs.to.keep$obs[k]]  <- obs.to.keep$new.ID[k] 
    }
    
    obs[is.na(obs)] <- obs.others
    obs[obs == old.no.obs] <- new.no.obs
    
    jags.data$obs <- obs
    jags.data$n.obs <- max(obs) - 1
    
    # Create a fixed observer effect input
    obs_counts <- table(obs.vec)
    top_obs_names.df <- data.frame(old.ID = names(sort(obs_counts, 
                                                       decreasing = TRUE))[1:N.obs],
                                   new.ID = c(1:N.obs))
                                   
    # Initialize the array with NAs or 0
    # Dimensions: Max Days, Max Stations, Total Years
    n_days_max <- dim(obs)[1] # or however you defined 'd'
    n_stations_max <- dim(obs)[2]
    n_years <- dim(obs)[3]
    
    obs_array <- array((N.obs+1), 
                       dim = c(n_days_max, 
                               n_stations_max, 
                               n_years))
    
    # Fill the array
    # Assuming your dataframe 'data' has columns: year_index, station_index, day_index _ this doesn't work
    #START HERE 2026-02-08!!
    for (i in 1:nrow(top_obs_names.df)){
      obs_array[obs == top_obs_names.df[i,"old.ID"]] <- top_obs_names.df[i, "new.ID"]
    }
      
    jags.data$obs.fixed <- obs_array
    jags.data$n.obs.fixed <- (max(top_obs_names.df$new.ID)) + 1
    
    #jags.data$scaled.day <- jags.data$day-(max.day/2)
    #jags.data["N"] <- NULL
    ###  ###  ###
    
    # center and scale VS and BF
    vs.std <- jags.data$vs
    for (k in 1:dim(vs.std)[3]){
      vs.std[,1,k] <- (vs.std[,1,k] - mean(vs.std[,1,k], na.rm = T))/sqrt(var(vs.std[,1,k], na.rm = T))
      vs.std[,2,k] <- (vs.std[,2,k] - mean(vs.std[,2,k], na.rm = T))/sqrt(var(vs.std[,2,k], na.rm = T))
    }
    
    bf.std <- jags.data$bf
    for (k in 1:dim(bf.std)[3]){
      bf.std[,1,k] <- (bf.std[,1,k] - mean(bf.std[,1,k], na.rm = T))/sqrt(var(bf.std[,1,k], na.rm = T))
      bf.std[,2,k] <- (bf.std[,2,k] - mean(bf.std[,2,k], na.rm = T))/sqrt(var(bf.std[,2,k], na.rm = T))
    }
    
    jags.data$bf.1 <- jags.data$bf
    jags.data$bf <- bf.std
    jags.data$vs.1 <- jags.data$vs
    jags.data$vs <- vs.std
    
    jags.data$start.years <- c(jags.input.list$jags.input.Laake$all.start.year,
                               jags.input.list$jags.input.new$start.years)
    
    jags.data$year.index <- jags.data$start.years - mean(jags.data$start.years)
    
    jags.input <- list(jags.data = jags.data,
                       min.dur = min.dur, 
                       jags.input.Laake = jags.input.list$jags.input.Laake,
                       jags.input.new = jags.input.list$jags.input.new,
                       jags.original.data = jags.input.list$jags.data,
                       data.dir = data.dir,
                       obs.summary = obs.summary)
    
    Start_Time<-Sys.time()
    
    jm <- jagsUI::jags(jags.data,
                       inits = NULL,
                       parameters.to.save= jags.params,
                       model.file = jags.model,
                       n.chains = MCMC.params$n.chains,
                       n.burnin = MCMC.params$n.burnin,
                       n.thin = MCMC.params$n.thin,
                       n.iter = MCMC.params$n.samples,
                       DIC = T,
                       parallel=T)
    
    Run_Time <- Sys.time() - Start_Time
    
    post <- posterior::as_draws(jm$samples)
    summary.posterior <- posterior::summarise_draws(post)
    
    jm.out <- list(jm = jm,
                   jags.input = jags.input,
                   #start.year = all.start.year,
                   jags.params = jags.params,
                   jags.model = jags.model,
                   MCMC.params = MCMC.params,
                   posterior.summary = summary.posterior,
                   Run_Time = Run_Time,
                   Run_Date = Run.date,
                   out.file.name = out.file.name,
                   Sys.env = Sys.getenv(),
                   new.run = TRUE)
    
    saveRDS(jm.out,
            file = out.file.name)
    
  } else {
    print("Previously saved results were read.")
    jm.out <- readRDS(file = out.file.name)
    jm.out$new.run <- FALSE
  }
  
  print(paste0("Finishing NoBUGS_Richards_fcn at ", Sys.time(), 
               " for Model: ", model.name))
  return(jm.out)
  
}


# Converts Granite Canyon count data to WinBUGS inputs. All raw (i.e., edited) data files 
# should be treated by Extract_Data_All_v2.Rmd. All output files should be in
# one directory (data.dir), e.g., V2.1_Nov2024. 
# years refer to years with raw data. I don't have raw data for 2006/2007 and 
# 2007/2008 and use WinBUGS input. 
# FOUND ERROR IN OUTPUT DATA WHERE SECONDARY OBSERVATIONS ARE ALL ZEROS 2024-12-12 -> This has been fixed.
# WinBUGS code does not allow differing watch lengths between primary and 
# secondary observations. So, time has to be assigned to secondary
data2WinBUGS_input <- function(data.dir, years, min.dur){
  
  library(abind)
  library(tidyverse)
  
  # this file contains all necessary inputs for 2006 - 2019 from Josh Stewart
  # We need 206/2007 and 2007/2008 data as we don't have the original data files
  data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")
  
  # e.g., 2007 refers to 2006/2007
  all.years <- c(2007, 2008, years)
  seasons <- sapply(all.years, 
                    FUN = function(x) paste0(x-1, "/", x))
  
  # These are extracted data files (Extract_Data_All_v2.Rmd)
  out.v2 <- lapply(years, 
                   FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
                                                    "_min", min.dur, 
                                                    "_Tomo_v2.rds")))
  
  all.Final_Data <- lapply(out.v2, FUN = function(x) x$Final_Data) 
  
  begin.primary <- lapply(out.v2, FUN = function(x){
    begin <- x$Final_Data %>%
      filter(station == "P") %>%
      select(begin) %>%
      pull()
    return(begin)
  })
  
  end.primary <- lapply(out.v2, FUN = function(x) {
    end <- x$Final_Data %>%
      filter(station == "P") %>%
      select(end) %>%
      pull()
    return(end)
  })
  
  begin.secondary <- lapply(out.v2, FUN = function(x){
    begin <- x$Final_Data %>%
      filter(station == "S") %>%
      select(begin) %>%
      pull()
    return(begin)
  })
  
  end.secondary <- lapply(out.v2, FUN = function(x) {
    end <- x$Final_Data %>%
      filter(station == "S") %>%
      select(end) %>%
      pull()
    return(end)
  })
  
  # Number of watch periods in each year's survey - before the 2019/2020 season
  # plus the new ones
  # This has been changed for 2024. I now have edited data for 2010 - 2024 seasons.
  # So, I can just use the first two (2006/2007 and 2007/2008). The two numbers for 
  # 2009/2010 and 2010/2011 don't match. In Durban's analysis, they were 164 and 178,
  # respectively.  136, 135
  periods <-c(data.0$periods[1:2],
              lapply(begin.primary, FUN = function(x) length(x)) %>% unlist)
  
  x <- length(periods)
  
  # out.file.name <- paste0("RData/WinBUGS_", all.years[1], "to", 
  #                         all.years[length(all.years)], "_v2_min", 
  #                         min.duration, "_", run.date, ".rds")
  
  Watch.Length. <- list()
  
  for (k in 1:length(begin.primary)){
    Watch.Length.[[k]] <- end.primary[[k]] - begin.primary[[k]]
  }
  
  # I don't have edited data for 2006/2007 and 2007/2008. So, they need to be
  # used as they were given in the WinBUGS input file from Durban. That's why
  # I take the first two columns. 
  # I may be able to use Laake's data for 2006/2007.  
  # 
  # whale count
  n <- data.0$n[,,1:2]  # take the first two years (2006/2007 and 2007/2008)
  n <- abind(n, array(0, replace(dim(n), 1, max(periods) - nrow(n))), along = 1)  %>%
    labelled::remove_attributes("dimnames")
  
  # the u data is whether there were observers on watch. 
  # 0 counts are often associated with years/shifts with 
  # no second observer. So if u=0, it will fix observation probability at 0
  # the second column for each year is for the second station - not the second
  # observer.
  u <- data.0$u[,,1:2]
  u <- abind(u, 
             array(0, replace(dim(u), 1, max(periods) - nrow(u))), along = 1) %>%
    labelled::remove_attributes("dimnames")
  
  # #visibility
  vs. <- lapply(out.v2, FUN = function(x) x$Final_Data$vs)
  
  vs <- rbind(data.0$vs[,1:2],
              array(NA, dim = c(max(periods) - nrow(data.0$vs), 2))) %>%
    labelled::remove_attributes("dimnames")
  
  # Beaufort
  bf. <- lapply(out.v2, FUN = function(x) x$Final_Data$bf)
  bf <- rbind(data.0$bf[,1:2],
              array(NA, dim = c(max(periods) - nrow(data.0$bf), 2)))%>%
    labelled::remove_attributes("dimnames")
  
  # A new observer list is created as new data are added. The new observer list
  # is saved in the Data directory. The list from the previous year is updated
  obs.list <- read.csv(file = paste0("Data/ObserverList", years[length(years)-1], ".csv"))
  
  obs.new <- unique(out.v2[[length(years)]]$Complete_Data$obs)
  new.obs <- obs.new[!c(obs.new %in% obs.list$obs)]
  obs.list <- rbind(obs.list,
                    data.frame(obs = new.obs,
                               ID = seq(max(obs.list$ID) + 1,
                                        max(obs.list$ID) + length(new.obs))))

  write.csv(obs.list, file = paste0("Data/ObserverList", max(years), ".csv"),
            row.names = FALSE)
  
  obs <- data.0$obs[,,1:2]
  
  # Obs==36 is no observer
  obs <- abind(obs, 
               array(36, replace(dim(obs), 1, max(periods) - nrow(obs))), along = 1)
  
  Watch.Length <- rbind(data.0$Watch.Length[,1:2],
                        array(NA, dim = c(max(periods) -
                                            nrow(data.0$Watch.Length), 2))) %>%
    labelled::remove_attributes("dimnames")
  
  day <- rbind(data.0$day[,1:2],
               array(NA, dim = c(max(periods) - nrow(data.0$day), 2))) %>%
    labelled::remove_attributes("dimnames")
  
  n.stations <- vector(mode = "numeric", length = length(begin.primary))
  k <- 5
  for (k in 1:length(begin.primary)){
    # need to pull out primary and secondary sightings if there were two stations
    Final_Data_k <- out.v2[[k]]$Final_Data %>% 
      mutate(f_station = as.factor(station))
    
    n.stations[k] <- length(levels(Final_Data_k$f_station))
    
    obs.year <- data.frame(obs = Final_Data_k$obs) %>% 
      left_join(obs.list, by = "obs")
    
    n.rows <- c(dim(Final_Data_k %>% filter(f_station == levels(Final_Data_k$f_station)[1]))[1],
                 dim(Final_Data_k %>% filter(f_station == levels(Final_Data_k$f_station)[2]))[1])
    
    n.k <- u.k <- obs.k <- array(data = 0, dim = c(max(n.rows), 2, 1))
    obs.k <- array(data = 36, dim = c(max(n.rows), 2, 1))

    Final_Data_k %>% 
      filter(f_station == "P") %>% 
      select(begin, end, dur, bf, vs, n, obs) %>%
      left_join(obs.list, by = "obs") -> data_P
    
    if (n.stations[k] == 1){
      n.k[,1,1] <- data_P$n
      u.k[,1,1] <- 1
      obs.k[,1,1] <- data_P$ID
      
    } else if (n.stations[k] == 2){
      
      Final_Data_k %>% 
        filter(f_station == "S") %>% 
        select(begin, end, dur, bf, vs, n, obs)  %>%
        left_join(obs.list, by = "obs")-> data_S
      
      # Find where the secondary observations need to be placed
      # Index for the closest time between primary and secondary
      idx.S <- lapply(data_S$begin, FUN = function(x){
        d.begin <- abs(data_P$begin - x) 
        which(min(d.begin) == d.begin)
        
      }) %>% unlist()
      
      # Closest differences between the primary and secondary
      dif.S <- lapply(data_S$begin, FUN = function(x){
        min(abs(data_P$begin - x) )
  
      }) %>% unlist()
      
      # if the difference is more than 30 minutes (0.02083 days), remove the 
      # secondary observation because there is no matching primary observation
      data_S <- data_S[dif.S < 0.02083,]
      idx.S <- idx.S[dif.S < 0.02083]
      
      #obs.P <- obs.year[Final_Data_k$f_station == "P", "ID"]
      
      n.k[1:n.rows[1], 1, 1] <- data_P$n
      n.k[idx.S, 2, 1] <- data_S$n
      
      u.k[1:n.rows[1], 1, 1] <- 1
      u.k[idx.S, 2, 1] <- 1
      
      obs.k[1:n.rows[1], 1, 1] <- data_P$ID
      obs.k[idx.S, 2, 1] <- data_S$ID
      
      
    }

    # if the new season has less rows than previous maximum
    if (nrow(n.k) < nrow(n)) {
      n.k.1 <- abind(n.k, 
                     array(0, dim = c(dim(n)[1] - dim(n.k)[1],
                                      2, 1)),
                     along = 1)
      n <- abind(n, n.k.1, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
      u.k.1 <- abind(u.k, 
                     array(0, dim = c(dim(n)[1] - dim(n.k)[1],
                                      2, 1)),
                     along = 1)
      u <- abind(u, u.k.1, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
      obs.k.1 <- abind(obs.k, 
                     array(36, dim = c(dim(n)[1] - dim(n.k)[1],
                                      2, 1)),
                     along = 1)
      obs <- abind(obs, obs.k.1, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
    } else if (nrow(n) < nrow(n.k)){
      n.1 <- abind(n, 
                   array(0, dim = c(dim(n.k)[1] - dim(n)[1],
                                    2, 1)),
                   along = 1)
      
      n <- abind(n.1, n.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
      u.1 <- abind(u, 
                   array(0, dim = c(dim(n.k)[1] - dim(n)[1],
                                    2, 1)),
                   along = 1)
      u <- abind(u.1, u.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
      obs.1 <- abind(obs, 
                   array(36, dim = c(dim(n.k)[1] - dim(n)[1],
                                    2, 1)),
                   along = 1)
      obs <- abind(obs.1, obs.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
    } else {
      n <- abind(n, n.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      u <- abind(u, u.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      obs <- abind(obs, obs.k, along = 3)  %>%
        labelled::remove_attributes("dimnames")
      
    }
    
    vs <- cbind(vs, c(data_P$vs, 
                      rep(NA, times = max(periods - length(data_P$vs))))) %>%
      labelled::remove_attributes("dimnames")
    
    
    bf <- cbind(bf, c(data_P$bf,
                      rep(NA, times = max(periods - length(data_P$bf))))) %>%
      labelled::remove_attributes("dimnames")
    
    Watch.Length <- cbind(Watch.Length,
                          c(Watch.Length.[[k]], 
                            rep(NA, 
                                times = max(periods) - length(Watch.Length.[[k]])))) %>%
      labelled::remove_attributes("dimnames")
    
    day <- cbind(day, 
                 c(floor(begin.primary[[k]]), 
                   rep(NA, times = max(periods) - length(begin.primary[[k]])))) %>%
      labelled::remove_attributes("dimnames")
    
  }
  
  #we're going to make N a partially observed data object with anchor points at day 1 and 90
  # TE: I don't know how these numbers were created... they are generally 2x n (not all)
  # N_inits <- as.matrix(read.table("Data/Initial Values/N_inits.txt",
  #                                 header=T))
  
  N_inits1 <- n[, 1,] * 2 + 2
  N_inits2 <- n[, 2,] * 2 + 2 
  
  N_inits <- N_inits1
  N_inits[N_inits1 < N_inits2] <- N_inits2[N_inits1 < N_inits2]
  
  N_inits <- rbind(N_inits,
                   matrix(data = NA, nrow = 2, ncol = length(periods)))
  
  for (k in 1:length(periods)){
    N_inits[(periods[k]+1):nrow(N_inits), k] <- NA  
  }
  
  #The 'data' has to be the inverse of the inits, 
  # with NAs for all of the estimated Ns, and 0s for the days 1 and 90
  N <- matrix(NA, nrow=max(periods)+2, ncol=length(periods)) 
  
  for(i in 1:length(periods)){
    N[(periods[i]+1):(periods[i]+2),i] <- 0 #True number of whales passing fixed at 0 for day 1 and 90
  }
  
  end <- Watch.Length + day
  
  #t <- round((begin+end)/2)
  t <- round((day+end)/2)
  
  # #Add a couple of extra rows of NAs to the end of the day index reference to 
  # match up with the fixed 0s in N (above), assigning them to days 1 and 90
  day <- rbind(as.matrix(day),
               matrix(NA, nrow=2, ncol=length(periods)))
  
  
  for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
    day[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
  }
  
  t <- rbind(as.matrix(t),
             matrix(NA,nrow=2,ncol=length(periods)))
  for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
    t[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
  }
  
  #Place 36s for 'no observer' for the two periods following the end of true 
  #watches (this is for the day 1 and day 90 zero-whale anchor points)
  #this will force it to the mean observation probability with no observer effect
  
  Watch.Length <- rbind(as.matrix(Watch.Length),
                        matrix(NA, nrow=2, ncol=length(periods)))
  
  for(i in 1:length(periods)){
    Watch.Length[(periods[i]+1):(periods[i]+2),i] <- 1
  }
  
  BUGS.data <- list(n = n[1:max(periods[1:x]),,1:x],
                    n.com = n[1:max(periods[1:x]),,1:x],
                    n.sp = n[1:max(periods[1:x]),,1:x],
                    n.station = dim(n[1:max(periods[1:x]),,1:x])[2],
                    n.year = dim(n[1:max(periods[1:x]),,1:x])[3],
                    n.obs = max(obs[1:max(periods[1:x]),,1:x], na.rm = T),
                    periods = periods[1:x],
                    obs = obs[1:max(periods[1:x]),,1:x],
                    #Watch.Length = 0.0625,
                    u = u[1:max(periods[1:x]),,1:x],
                    vs = vs[1:max(periods[1:x]),1:x],
                    bf = bf[1:max(periods[1:x]),1:x],
                    #day=day,
                    day = t[1:(max(periods[1:x])+2),1:x],
                    N = N[,1:x],
                    N.com = N[,1:x],
                    N.sp = N[,1:x],
                    knot = c(-1.46,-1.26,-1.02,-0.78,
                             -0.58,-0.34,-0.10,0.10,
                             0.34,0.57,0.78,1.02,1.26,1.46),
                    n.knots=14,
                    #begin=begin,
                    #end=end,
                    Watch.Length=Watch.Length[1:(max(periods[1:x])+2), 1:x])
  
  BUGS.inits <- function() list(mean.prob = 0.5,
                                BF.Fixed = 0,
                                VS.Fixed = 0,
                                mean.prob.sp = 0.5,
                                BF.Fixed.sp = 0,
                                VS.Fixed.sp = 0,
                                mean.prob.com = 0.5,
                                BF.Fixed.com = 0,
                                VS.Fixed.com = 0,
                                mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                                beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
                                BF.Switch = 1,
                                VS.Switch = 1,
                                OBS.Switch = 1,
                                sigma.Obs = 1,
                                BF.Switch.sp = 1,
                                VS.Switch.sp = 1,
                                OBS.Switch.sp = 1,
                                sigma.Obs.sp = 1,
                                BF.Switch.com = 1,
                                VS.Switch.com = 1,
                                OBS.Switch.com = 1,
                                sigma.Obs.com = 1,
                                N = N_inits,
                                N.com = N_inits,
                                N.sp = N_inits,
                                #z = matrix(1,nrow=90,ncol=6),
                                beta.sp = array(data=0, dim=c(2,x)),
                                sd.b.sp = rep(1, times = x), #c(1,1,1,1,1,1),
                                z = matrix(1, nrow=90, ncol= x))
  
  return(out.list <- list(data = BUGS.data,
                          inits = BUGS.inits,
                          all.years = all.years,
                          seasons = seasons,
                          min.dur = min.dur,
                          All_Data = out.v2,
                          Final_Data = all.Final_Data))  
}

# Create Jags input for Laake's data.
# min.dur is the minimum effort duration in minutes to be included in the analysis
LaakeData2JagsInput <- function(min.dur, max.day = 100){
  
  if (!file.exists("RData/Jags_Richards_LaakeData.rds")){
    
    library(tidyverse)
    
    # From example code in the ERAnalysis library
    # 
    # The recent survey data 1987 and after are stored in ERSurveyData and those data
    # are processed by the ERAbund program to produce files of sightings and effort.
    # The sightings files are split into Primary observer and Secondary observer sightings.
    # Primary observer sightings are whales that are not travelling North and are defined by
    # those when EXPERIMENT==1 (single observer) or a designated LOCATION when EXPERIMENT==2.
    #  For surveys 2000/2001 and 2001/2002, the primary observer was at LOCATION=="N"
    # and for all other years, LOCATION=="S".
    #
    # Based on the projected timing of the passage of the whale (t241) perpendicular to the 
    # watch station, the sighting was either contained in the watch (on effort) or not (off effort).
    # The dataframe Primary contains all of the on effort sightings and PrimaryOff contains all
    # of the off-effort sightings.  
    #
    #data(PrimaryOff)   # off-effort sightings
    load("Data/Primary.rda")      # on-effort sightings from 1987 to 2007
    load("Data/ERSurveyData.rda")
    load("Data/Observer.rda")
    
    # The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
    # are less than or equal to 4. Below the counts are shown for the 2 dataframes for
    # recent surveys since 1987/88.
    #table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
    #load("Data/PrimarySightings.rda")
    load("Data/PrimaryEffort.rda")
    load("Data/SecondaryEffort.rda")
    
    # Likewise, the secondary sightings are those with EXPERIMENT==2 but the LOCATION that
    # is not designated as primary.  there is no effort data for the secondary sightings... 
    # so, can't use it for BUGS/jags - ignore it for now.
    # data(SecondarySightings)
    
    # Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
    # exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
    # Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
    # beaufort exceeded 4 during the watch.
    # 
    # Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
    # Note that for most years <1987 there are no records with Use==FALSE because the filtered records
    # were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
    # one watch (5 hours) was missing a beaufort value so it was excluded.
    # tapply(PrimaryEffort$effort,
    #        list(PrimaryEffort$Use,
    #             PrimaryEffort$Start.year),
    #        sum)*24
    #        
    
    # Filter effort and sightings and store in dataframes Effort and Sightings
    Laake_PrimaryEffort <- PrimaryEffort[PrimaryEffort$Use,]  
    
    Laake_SecondaryEffort <- SecondaryEffort[SecondaryEffort$Use,]
    
    # Define shifts and combine effort within each shift:
    # Define shifts, group by date and shift, calculate cumulative effort, then
    # find the maximum cumulative effort for each shfit.
    #Primary.new.date <- fractional_Day2YMDhms(Laake_PrimaryEffort$end,
    #                                  Laake_PrimaryEffort$Start.year + 1)
    
    Laake_PrimaryEffort %>%
      mutate(Shift = shift.definition(fractional_Day2YMDhms(begin, 
                                                            Start.year + 1)$YMD,
                                      fractional_Day2YMDhms(begin, 
                                                            Start.year + 1)$hms)) %>%
      group_by(Date, Shift) %>% 
      mutate(cumu.effort = cumsum(effort),
             effort.per.shift = max(cumu.effort)) -> Laake_PrimaryEffort.all
    
    Laake_SecondaryEffort %>%
      mutate(Shift = shift.definition(fractional_Day2YMDhms(begin, 
                                                            Start.year + 1)$YMD,
                                      fractional_Day2YMDhms(begin, 
                                                            Start.year + 1)$hms)) %>%
      group_by(Date, Shift) %>% 
      mutate(cumu.effort = cumsum(effort),
             effort.per.shift = max(cumu.effort)) -> Laake_SecondaryEffort.all
    
    # Filter to minimum duration. min.dur is given in the unit of minutes, which
    # is converted to the unit of days:
    
    Laake_PrimaryEffort.all %>%
      filter(effort.per.shift >= min.dur / (60 * 24) ) -> Laake_PrimaryEffort
    
    Laake_SecondaryEffort.all %>%
      filter(effort.per.shift >= min.dur / (60 * 24) ) -> Laake_SecondaryEffort
    
    # Sightings = PrimarySightings
    # Sightings$seq = 1:nrow(Sightings)
    # Sightings = merge(Sightings, subset(Laake_PrimaryEffort, select=c("key")))
    # Sightings = Sightings[order(Sightings$seq),]
    
    # filter off-effort sightings and high Beaufort/vis lines from secondary sightings
    # but... there is no effort data for the secondary sightings... so, can't use it
    # for BUGS/jags - ignore it for now.
    # SecondarySightings %>% 
    #   filter(vis < 5, beaufort < 5, is.na(off)) -> secondary.sightings
    
    # For jags and WinBugs code, what I need are
    # 1. observed number of whales per day n[d, s, y], where d = # days since 12/1,
    # s = station (1 = primary, 2 = secondary), y = year. For Laake's data, maximum 
    # number of days per year was 94. 
    # 2. Beaufort sea state bf[d,y]
    # 3. Visibility code vs[d,y]
    # 4. observer code obs[d,s,y]
    # 5. the proportion of watch duration per day (fractional day, NOT hours per 
    # maximum observation period) watch.prop[d,y]
    # 6. index of survey day, i.e., the number of days since 12/1 day[d,y]
    
    # Need to count the number of days since 12-01 for each year. But 12-01 is 1.
    # Then... 
    # Count the number of whales per day and daily effort
    # In early years, surveys were conducted 10 hrs. So, the watch proportion
    # can be > 1.0, because we have used 9 hrs as maximum. 
    
    # Summarizing by day worked fine but the model requires counts per observation
    # period. Needs to be redone. 2023-09-15 DONE.
    
    Laake_PrimaryEffort %>% 
      mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
             dt = as.numeric(as.Date(Date) - Day1) + 1,
             obs = Observer) %>%
      select(Start.year, nwhales, effort, vis, 
             beaufort, obs, dt, Date, Shift, effort.per.shift) %>%
      group_by(Start.year) %>%
      mutate(shift.hr = effort.per.shift * 24,  # effort in hours
             shift.prop = shift.hr/9) -> Effort.by.period
    
    # observers - rather than using entries of Observer.rda, I pull out all observers
    # from the effort objects (primary and secondary) and re-number them. 
    obs <- c(unique(Laake_PrimaryEffort$Observer),
             unique(Laake_SecondaryEffort$Observer)) %>% unique()
    
    Laake.obs <- data.frame(obs = c(as.vector(obs), "No obs"),
                            ID.new = seq(1:(length(obs) + 1)))
    
    no.obs.ID <- max(Laake.obs$ID.new)
    
    # effort is in the unit of days. It is "effort" and not "effort.per.shift"
    # because a "period" is defined with Beaufort/Vis code and the "effort" goes
    # with the pair. "effort.per.shift" was used to extract long enough "shift"
    # for including all observations. 
    # For each season, day, and shift, the sightings need to be pooled
    # for the same condition, in order to make Laake's dataset comparable
    # to mine. Currently, each sighting is listed as one "row." This 
    # inflates the number of sightings under the same condition and 
    # observer. This has been fixed in the following. 2026-02-11 
    
    new.primary <- list()
    Laake_PrimaryEffort %>%
      count(Start.year, Date, Shift, vis, beaufort) -> Laake.unique.combos 
    k<- 3285
    for (k in 1:nrow(Laake.unique.combos)){
      Laake_PrimaryEffort %>%
        filter(Start.year == Laake.unique.combos$Start.year[k],
               Date == Laake.unique.combos$Date[k],
               Shift == Laake.unique.combos$Shift[k],
               vis == Laake.unique.combos$vis[k],
               beaufort == Laake.unique.combos$beaufort[k]) -> tmp
      
      if (Laake.unique.combos$n[k] == 1){
        tmp %>%
          select(Date, Start.year, nwhales, Observer, vis,
                 beaufort, effort.per.shift, Shift,
                 effort, begin, end, npods) %>%
          #group_by(Start.year) %>%
          mutate(Day = as.Date(Date)-as.Date(paste0(Start.year, "-12-01"))+1,
                 Season = paste0(Start.year, "/", Start.year+1),
                 shift.prop = (effort.per.shift*24)/9) %>%
          rename(Year = Start.year,
                 n = nwhales,
                 obs = Observer,
                 vs = vis,
                 bf = beaufort,
                 shift.length = effort.per.shift,
                 watch.length = effort) %>%
          select(Date, Year, n, Day, Season, obs, vs, bf,
                 shift.length, watch.length, shift.prop, begin,
                 end, npods, Shift) -> new.primary[[k]] 
      } else {
        tmp %>%
          select(Date, Start.year, nwhales, Observer, vis,
                 beaufort, effort.per.shift, Shift,
                 effort, begin, end, npods) %>%
          #group_by(Start.year)%>%
          summarize(Date = first(Date),
                    Year = first(Start.year),
                    n = sum(nwhales),
                    Day = as.Date(Date)-as.Date(paste0(Year, "-12-01"))+1,
                    Season = paste0(Year, "/", Year+1),
                    obs = first(Observer),
                    vs = first(vis),
                    bf = first(beaufort),
                    shift.length = first(effort.per.shift),
                    watch.length = sum(effort),
                    shift.prop = (first(effort.per.shift) * 24)/9,
                    begin = first(begin),
                    end = last(end),
                    npods = sum(npods),
                    Shift = first(Shift),
                    .groups = "drop") -> new.primary[[k]]
      }
    }  
    
    Laake.primary.counts <- do.call("rbind", new.primary)
    
    # Laake_PrimaryEffort %>% 
    #   group_by(Start.year) %>%
    #   reframe(Date = Date,
    #           Year = first(Start.year),
    #           n = nwhales,
    #           Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")) + 1,
    #           Season = paste0(Year, "/", Year + 1),
    #           obs = Observer,
    #           vs = vis,
    #           bf = beaufort,
    #           shift.length = effort.per.shift, # length of each shift
    #           watch.length = effort, # watch length of each row
    #           shift.prop = (effort.per.shift * 24)/9) -> Laake.primary.counts
    
    new.secondary <- list()
    Laake_SecondaryEffort %>%
      count(Start.year, Date, Shift, vis, beaufort) -> Laake.unique.combos 
    k<- 1
    for (k in 1:nrow(Laake.unique.combos)){
      Laake_SecondaryEffort %>%
        filter(Start.year == Laake.unique.combos$Start.year[k],
               Date == Laake.unique.combos$Date[k],
               Shift == Laake.unique.combos$Shift[k],
               vis == Laake.unique.combos$vis[k],
               beaufort == Laake.unique.combos$beaufort[k]) -> tmp
      
      if (Laake.unique.combos$n[k] == 1){
        tmp %>%
          select(Date, Start.year, nwhales, Observer, vis,
                 beaufort, effort.per.shift, Shift,
                 effort, begin, end, npods) %>%
          #group_by(Start.year) %>%
          mutate(Day = as.Date(Date)-as.Date(paste0(Start.year, "-12-01"))+1,
                 Season = paste0(Start.year, "/", Start.year+1),
                 shift.prop = (effort.per.shift*24)/9) %>%
          rename(Year = Start.year,
                 n = nwhales,
                 obs = Observer,
                 vs = vis,
                 bf = beaufort,
                 shift.length = effort.per.shift,
                 watch.length = effort) %>%
          select(Date, Year, n, Day, Season, obs, vs, bf,
                 shift.length, watch.length, shift.prop, begin,
                 end, npods, Shift) -> new.secondary[[k]] 
      } else {
        tmp %>%
          select(Date, Start.year, nwhales, Observer, vis,
                 beaufort, effort.per.shift, Shift,
                 effort, begin, end, npods) %>%
          #group_by(Start.year) %>%
          summarize(Date = first(Date),
                    Year = first(Start.year),
                    n = sum(nwhales),
                    Day = as.Date(Date)-as.Date(paste0(Year, "-12-01"))+1,
                    Season = paste0(Year, "/", Year+1),
                    obs = first(Observer),
                    vs = first(vis),
                    bf = first(beaufort),
                    shift.length = first(effort.per.shift),
                    watch.length = sum(effort),
                    shift.prop = (first(effort.per.shift) * 24)/9,
                    begin = first(begin),
                    end = last(end),
                    npods = sum(npods),
                    Shift = first(Shift),
                    .groups = "drop") -> new.secondary[[k]]
      }
    }  
    
    Laake.secondary.counts <- do.call("rbind", new.secondary)
    
    # Laake_SecondaryEffort %>% 
    #   group_by(Start.year) %>%
    #   reframe(Date = Date,
    #           Year = first(Start.year),
    #           n = nwhales,
    #           Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")) + 1,
    #           Season = paste0(Year, "/", Year + 1),
    #           obs = Observer,
    #           vs = vis,
    #           bf = beaufort,
    #           shift.length = effort.per.shift, # length of each shift
    #           watch.length = effort, # watch length of each row
    #           shift.prop = (effort.per.shift * 24)/9) -> Laake.secondary.counts
    
    # Laake.primary.counts %>% 
    #   group_by(Date) %>%
    #   reframe(Year = first(Start.year),
    #           Date = Date,
    #           Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")) + 1,
    #           Season = paste0(Year, "/", Year + 1),
    #           Daily.n = sum(n)) -> Laake.primary.daily.counts
    
    # Find double observer years
    double.obs.year <- unique(Laake_SecondaryEffort$Start.year) 
    all.year <- unique(Laake_PrimaryEffort$Start.year)
    
    n.station <- rep(1, length(all.year))
    n.station[all.year %in% double.obs.year] <- 2
    
    Laake.primary.counts %>%
      group_by(Year) %>%
      summarize(Season = first(Season),
                periods = first(n())) -> Laake.primary.periods
    
    Laake.secondary.counts %>%
      group_by(Year) %>%
      summarize(Season = first(Season),
                periods = first(n())) -> Laake.secondary.periods
    
    Laake.primary.periods %>% 
      left_join(Laake.secondary.periods, by = "Season") %>%
      select(Year.x, Season, periods.x, periods.y) %>%
      transmute(Start.year = Year.x,
                Season = Season,
                periods.1 = periods.x,
                periods.2 = periods.y) -> Laake.periods
    
    bf <- bf.centered <- vs <- vs.centered <- watch.prop <-  watch.length <- shift.prop <- array(dim = c(max(Laake.primary.periods$periods),
                                                                                                         2, length(all.year)))
    
    obs.input <- n.Laake <- day <- array(dim = c(max(Laake.primary.periods$periods) + 2,
                                                 2, length(all.year)))
    y <- 1
    y2 <- 1
    for (y in 1:length(all.year)){
      Laake.primary.counts %>%
        filter(Year == all.year[y]) %>%
        arrange(Day) %>%
        left_join(Laake.obs, by = "obs") -> temp.data
      
      n.Laake[1:(Laake.primary.periods$periods[y]+2), 1, y] <- c(0, temp.data$n, 0)
      obs.input[1:(Laake.primary.periods$periods[y]+2), 1, y] <- c(no.obs.ID, temp.data$ID.new, no.obs.ID)
      shift.prop[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$shift.prop
      watch.length[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$watch.length
      #shift.length[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$shift.length
      bf[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$bf 
      vs[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$vs 
      bf.centered[1:Laake.primary.periods$periods[y], 1, y] <- scale(temp.data$bf, scale = F)
      vs.centered[1:Laake.primary.periods$periods[y], 1, y] <- scale(temp.data$vs, scale = F)
      
      day[1:(Laake.primary.periods$periods[y]+2), 1, y] <- c(1, as.numeric(temp.data$Day), max.day)
      
      # fill in the secondary observations
      if (isTRUE(all.year[y] %in% double.obs.year)){
        temp.data <- Laake.secondary.counts %>%
          filter(Year == all.year[y])%>%
          left_join(Laake.obs, by = "obs")
        
        n.Laake[1:(Laake.secondary.periods$periods[y2]+2), 2, y] <- c(0,temp.data$n,0)
        obs.input[1:(Laake.secondary.periods$periods[y2]+2), 2, y] <- c(no.obs.ID, temp.data$ID.new, no.obs.ID)
        shift.prop[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$shift.prop
        watch.length[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$watch.length
        #shift.length[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$shift.length
        bf[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$bf #scale(temp.data$bf)
        vs[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$vs #scale(temp.data$vs)
        bf.centered[1:Laake.secondary.periods$periods[y2], 2, y] <- scale(temp.data$bf, scale = F)
        vs.centered[1:Laake.secondary.periods$periods[y2], 2, y] <- scale(temp.data$vs, scale = F)
        day[1:(Laake.secondary.periods$periods[y2]+2), 2, y] <- c(1, as.numeric(temp.data$Day), max.day)
        
        y2 <- y2 + 1
        
      }
    }
    
    N <- matrix(data = NA, nrow = max.day, ncol = dim(n.Laake)[3])
    N[1,] <- 0
    N[max.day,] <- 0
    
    jags.data <- list(  n = n.Laake,
                        n.station = n.station,
                        n.year = length(all.year),
                        n.obs = nrow(Laake.obs),
                        periods = Laake.periods %>% 
                          select(periods.1, periods.2) %>% simplify2array(),
                        obs = obs.input,
                        vs = vs,
                        bf = bf,
                        vs.centered = vs.centered,
                        bf.centered = bf.centered,
                        watch.prop = watch.prop,
                        watch.length = watch.length,
                        day = day,
                        n.days = max(day, na.rm = T),
                        N = N)
    
    out.list <- list(jags.data = jags.data,
                     all.start.year = all.year,
                     double.obs.year = double.obs.year,
                     obs = Laake.obs,
                     no.obs.ID = no.obs.ID,
                     min.dur = min.dur,
                     primary.counts = Laake.primary.counts,
                     secondary.counts = Laake.secondary.counts,
                     seasons = paste0(all.year, "/", all.year + 1),
                     primary.all = Laake_PrimaryEffort.all,
                     secondary.all = Laake_SecondaryEffort.all,
                     primary.filtered = Laake_PrimaryEffort,
                     secondary.filtered = Laake_SecondaryEffort)
    
    saveRDS(out.list,
            file = "RData/Jags_Richards_LaakeData.rds")
  } else {
    
    out.list <- readRDS("RData/Jags_Richards_LaakeData.rds")
  }
  return(out.list)  
}


# Create data input for Jags without running WinBUGS. Only works if
# there are raw data files. 
data2Jags_input_NoBUGS <- function(min.dur,
                                   years,
                                   data.dir,
                                   max.day = 100){
  
  seasons <- sapply(years, FUN = function(x) paste0(x-1, "/", x))
  start.years <- years - 1
  
  # These are extracted data files (Extract_Data_All_v2.Rmd)
  out.v2 <- lapply(years, 
                   FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
                                                    "_min", min.dur, "_Tomo_v2.rds")))
  
  # Get all sightings and make a data frame
  primary.sightings.list <- lapply(out.v2,
                                   FUN = function(x) {
                                     tmp <- do.call(rbind, x$sightings_primary)
                                     tmp %>%
                                       mutate(Observer = toupper(Observer),
                                              year = year(as.Date(Date, format = "%m/%d/%Y")),
                                              Start.year = min(year),
                                              station = "P") -> tmp.1
                                     return(tmp.1)
                                   })
  
  primary.sightings.df <- do.call(rbind, primary.sightings.list) 
  
  n.stations <- lapply(out.v2,
                       FUN = function(x){
                         dat <- x$Final_Data %>%
                           mutate(f.station = as.factor(station))
                         return(levels(dat$f.station) %>%
                                  length() %>%
                                  unlist())
                       }) %>%
    unlist()
  
  if (max(n.stations) > 1){
    secondary.sightings.list <- lapply(out.v2,
                                       FUN = function(x) {
                                         tmp <- do.call(rbind, x$sightings_secondary)
                                         return(tmp)})
    # Secondary doesn't exist for all years
    secondary.idx.null <- lapply(secondary.sightings.list, is.null) %>%
      unlist()
    
    secondary.sightings.df <- do.call(rbind, secondary.sightings.list[!secondary.idx.null]) %>%
      mutate(Observer = toupper(Observer),
             year = year(as.Date(Date, format = "%m/%d/%Y")),
             Start.year = min(year),
             station = "S")
    # observers - include all observers from Laake's data
    # Some observers were only found in the secondary stations. So, they 
    # need to be added. 
    obs.list <- create.observer.list(rbind(primary.sightings.df,
                                           secondary.sightings.df) %>%
                                       select(Start.year, Observer))
    
    all.observers <- obs.list$unique
    sightings.df <- rbind(primary.sightings.df, secondary.sightings.df)
    
  } else {
    obs.list <- create.observer.list(primary.sightings.df %>%
                                       select(Start.year, Observer))
    all.observers <- obs.list$unique
    
    sightings.df <- primary.sightings.df
  }

  begin. <- lapply(out.v2, FUN = function(x) x$Final_Data$begin)
  end. <- lapply(out.v2, FUN = function(x) x$Final_Data$end)
  
  periods <- lapply(out.v2, 
                    FUN = function(x){
                      x$Final_Data %>%
                        filter(station == "P") %>%
                        nrow() -> n.row
                      return(n.row)
                    }) %>% unlist()
  
  x <- length(periods)
  
  # out.file.name <- paste0("RData/JAGS_", years[1], "to", 
  #                         years[length(years)], "_v2_min", 
  #                         min.dur, "_", run.date, ".rds")

  # watch length in fractional days  
  Watch.Length <- list()

  # watch lengths in the unit of "days"
  for (k in 1:length(begin.)){
    Watch.Length[[k]] <- end.[[k]] - begin.[[k]]
  }
  
  day <- array(NA, 
               dim = c(max(periods)+2, 
                       2, length(years)))
  
  watch.length <- vs <- bf <- array(data = NA, 
                                    dim = c(max(periods), 
                                            2, length(years)))
  
  periods.mat <- array(data = NA,
                       dim = c(length(periods), 2))
  
  no.obs.ID <- max(all.observers$ID) + 1
  
  n <- array(NA, dim = c(max(periods) + 2, 
                         2, length(years)))
  
  
  obs <- array(NA, dim = dim(n))
  
  k <- 1
  for (k in 1:length(begin.)){
    Final_data <- out.v2[[k]]$Final_Data %>%
      mutate(watch.length = end - begin,
             f.station = as.factor(station))
    
    if (n.stations[k] == 2){
      Final_data.P <- Final_data %>% 
        filter(f.station == "P")
      Final_data.S <- Final_data %>%
        filter(f.station == "S")
      
      # fill in the primary 
      n.row.P <- (length(Final_data.P$n) + 2)
      n[1:n.row.P, 1, k] <- c(0, Final_data.P$n, 0)
      vs[1:(n.row.P-2), 1, k] <- c(Final_data.P$vs)
      bf[1:(n.row.P-2), 1, k] <- c(Final_data.P$bf)
      
      obs.year.P <- data.frame(obs = Final_data.P$obs) %>% 
        left_join(all.observers, by = "obs")
      
      obs[1:n.row.P,1,k] <- c(no.obs.ID, obs.year.P$ID, no.obs.ID)
      # 0.375 is 9 hrs and watch.prop = 1.0
      watch.length[1:(n.row.P-2),1,k] <- c(Final_data.P$watch.length)
      day[1:n.row.P,1,k] <- c(1, floor(Final_data.P$begin), max.day)
      periods.mat[k,1] <- nrow(Final_data.P) 
      
      # fill in the secondary 
      n.row.S <- length(Final_data.S$n) + 2
      n[1:n.row.S,2,k] <- c(0, Final_data.S$n, 0)
      vs[1:(n.row.S-2),2,k] <- c(Final_data.S$vs)
      bf[1:(n.row.S-2),2,k] <- c(Final_data.S$bf)
      
      obs.year.S <- data.frame(obs = Final_data.S$obs) %>% 
        left_join(all.observers, by = "obs")
      
      obs[1:n.row.S,2,k] <- c(no.obs.ID, obs.year.S$ID, no.obs.ID)
      watch.length[1:(n.row.S-2),2,k] <- c(Final_data.S$watch.length)
      day[1:n.row.S,2,k] <- c(1, floor(Final_data.S$begin), max.day)
      
      periods.mat[k,2] <- nrow(Final_data.S)
      
    } else {
      n.row <- length(Final_data$n) + 2
      n[1:n.row,1,k] <- c(0, Final_data$n, 0)
      vs[1:(n.row-2),1,k] <- c(Final_data$vs)
      bf[1:(n.row-2),1,k] <- c(Final_data$bf)
      
      obs.year <- data.frame(obs = Final_data$obs) %>% 
        left_join(all.observers, by = "obs")
      
      obs[1:n.row,1,k] <- c(no.obs.ID, obs.year$ID, no.obs.ID)
      watch.length[1:(n.row-2),1,k] <- c(Final_data$watch.length)
      day[1:n.row,1,k] <- c(1, floor(Final_data$begin), max.day)
      periods.mat[k,1] <- nrow(Final_data)
      
    }
    
  }
  
  # Renumber observer IDs because some observers in the list never
  # showed up in the data. The maximum ID number can be greater than
  # the total number of observers, which returns error when running
  # JAGS
  uniq.obs.1 <- apply(obs[,1,], FUN = unique, MARGIN = 2) %>% 
    unlist() %>% 
    unique() %>%
    sort()
  
  uniq.obs.2 <- apply(obs[,2,], FUN = unique, MARGIN = 2) %>% 
    unlist() %>% 
    unique() %>%
    sort()

  uniq.obs <- unique(c(uniq.obs.1, uniq.obs.2)) %>% sort
  
  uniq.obs.df <- data.frame(ID = uniq.obs, 
                            new.ID = seq(1, length(uniq.obs)))
  
  all.observers %>%
    filter(ID %in% uniq.obs) %>%
    arrange(ID) %>% #-> tmp #
    left_join(uniq.obs.df, by = "ID") %>% #-> tmp
    select(-c(data)) %>%
    transmute(ID = ID, obs = obs, ID.new = new.ID) -> new.obs.df
  
  no.obs.ID.2 <- max(new.obs.df$ID.new) + 1
  new.obs.df <- rbind(new.obs.df, 
                      data.frame(ID = no.obs.ID,
                                 obs = "No obs", 
                                 ID.new = no.obs.ID.2))
  
  obs[is.na(obs)] <- no.obs.ID
  # use "replace" to swap the original observer IDs with new IDs
  # using the look up table (new.obs.df)
  # c(new, x)[match(x, c(old, x))], where x is data
  # https://stackoverflow.com/questions/16228160/multiple-replacement-in-r
  obs.new <- array(data = NA, #no.obs.ID.2,
                   dim = dim(obs))
  
  obs.new[,1,] <- apply(obs[,1,], 
                        FUN = function(x) c(as.vector(new.obs.df$ID.new), 
                                            x)[match(x, as.vector(new.obs.df$ID), x)], 
                        MARGIN = 2)
  
  obs.new[,2,] <- apply(obs[,2,], 
                        FUN = function(x) c(as.vector(new.obs.df$ID.new), 
                                            x)[match(x, as.vector(new.obs.df$ID), x)], 
                        MARGIN = 2)
  
  # N is a partially observed (by assumption)
  # N <- array(dim = dim(day))
  # N[day == 1] <- 0
  # N[day == max.day] <- 0
  # 
  # Or a 2D
  N <- matrix(data = NA, nrow = max.day, ncol = dim(n)[3])
  N[1,] <- 0
  N[max.day,] <- 0

  # Find the maximum observed number of whales per day
  max.n.day <- array(data = 0, dim = c(100, length(years)))
  for (k in 1:length(years)){
    n.days <- unique(c(day[,1,k], day[,2,k])) %>% na.omit() %>% sort()
    n.days <- n.days[n.days > 1 & n.days < 100]
    for (k1 in 1:length(n.days)){
      n.1 <- n[day[,2,k] == n.days[k1], 1, k]
      max.1 <- ifelse(sum(!is.na(n.1)) > 0,
                      max(n.1, na.rm = T),
                      NA)
      
      n.2 <- n[day[,2,k] == n.days[k1], 2, k]
      max.2 <- ifelse(sum(!is.na(n.2)) > 0,
                      max(n.2, na.rm = T),
                      NA)
      
      max.n.day[n.days[k1], k] <- max(c(max.1, max.2, na.rm = T))
      
    }
    
  }
  
  jags.data <- list(  n = n,
                      n.station = n.stations,
                      n.year = length(years),
                      n.obs = length(uniq.obs) - 1,
                      periods = periods.mat,
                      n.days = max(day, na.rm = T),
                      obs = obs.new,
                      vs = vs,
                      bf = bf,
                      vs.centered = apply(vs, MARGIN = c(2,3), 
                                 FUN = scale, scale = FALSE),
                      bf.centered = apply(bf, MARGIN = c(2,3), 
                                 FUN = scale, scale = FALSE),
                      watch.length = watch.length,
                      watch.prop = (watch.length * 24)/9,
                      day = day,
                      N = N,
                      max.n.day = max.n.day)
  
  out.list <- list(jags.data = jags.data,
                   original.obs = obs,
                   obs = new.obs.df,
                   no.obs.ID = no.obs.ID.2,
                   min.dur = min.dur, 
                   seasons = seasons, 
                   years = years,
                   start.years = start.years,
                   data.dir = data.dir,
                   sightings = sightings.df,
                   obs.df = new.obs.df)
  
  return(out.list)
  
}

# Create jags data input for all years without uring WinBUGS input
# AllData2JagsInput_NoBUGS(min.dur, years, data.dir)
# min.dur = the minimum duration of a shift to be included; 30, 60, 85, etc.
#           data need to be extracted using the value and results saved in a 
#           directory (data.dir)
# years = the years over which the analysis is conducted. Specify it in a vector.
#         For example, c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
#         Note that these years are "ending" years. For example, 2015 means the 
#         2014/2015 season. 
#         
AllData2JagsInput_NoBUGS <- function(min.dur, years, data.dir, max.day = 90){
  library(abind)
  # this converts Laake's data into jags input
  jags.input.Laake <- LaakeData2JagsInput(min.dur,
                                          max.day = max.day)
  
  # This pulls out data from 2010 onward and create jags input
  jags.input.new <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                           years = years,
                                           data.dir = data.dir,
                                           max.day = max.day)
  
  obs.list <- create.observer.list(jags.input.new$sightings %>%
                                     select(Start.year, Observer))
  
  new.no.obs.ID <- max(obs.list$all$ID) + 1
  
  # Observer IDs need to be combined and consistent between the two datasets
  # To do that, I need to create a look-up table using the obs.list$all 
  # dataframe
  jags.input.Laake$obs %>% 
    left_join(obs.list$all, by = "obs") %>% #-> tmp
    select(obs, ID.new, ID, data) %>%
    rename(ID.old = ID.new) %>%
    filter(data != "new.data") %>%
    rbind(c("No obs", jags.input.Laake$no.obs.ID, new.no.obs.ID, "pre2009")) -> Laake.observers
  
  jags.input.new$obs %>%
    left_join(obs.list$all, by = "obs") %>%
    filter(ID.old < jags.input.new$no.obs.ID) %>%
    select(obs, ID.old, ID.x, data) %>%
    rename(ID = ID.x) %>%
    filter(data == "new.data") %>%
    rbind(c("No obs", jags.input.new$no.obs.ID, new.no.obs.ID, "new.data")) -> new.data.observers
  
  obs.all <- rbind(Laake.observers, new.data.observers) %>%
    mutate(ID.num = as.numeric(ID)) %>%
    arrange(by = ID.num)
  
  # Renumber all IDs so that there are no skipped numbers
  uniq.ID.num <- unique(obs.all$ID.num)
  for (k in 1:length(uniq.ID.num)){
    obs.all$new.ID[obs.all$ID.num == uniq.ID.num[k]] <- k
  }
  
  # separte the dataset again:
  obs.all %>%
    filter(data == "pre2009") -> Laake.observers.new
  obs.all %>%
    filter(data == "new.data") -> new.data.observers.new
  
  no.obs.ID.new <- Laake.observers.new %>% 
    filter(obs == "No obs") %>%
    pull(new.ID)
  
  # Replace observer IDs in each dataset with new IDs (obs.all.new$ID.all)
  jags.obs.Laake <- jags.input.Laake$jags.data$obs
  
  # Replace NAs with no observer ID for Laake data
  Laake.no.ID <- Laake.observers.new %>%
    filter(obs == "No obs") %>%
    pull(ID.old) %>%
    as.numeric()
  
  jags.obs.Laake[is.na(jags.obs.Laake)] <- Laake.no.ID
  
  # use "match" to swap the original observer IDs with new IDs
  # using the look up table (obs.all.new)
  # c(new, x)[match(x, c(old, x))], where x is data
  # https://stackoverflow.com/questions/16228160/multiple-replacement-in-r
 
  # Laake observer list first - this list doesn't change
  # Laake data doesn't have no-observer ID 
  jags.obs.Laake.new <- array(data = no.obs.ID.new,
                              dim = dim(jags.obs.Laake))
  
  jags.obs.Laake.new[,1,] <- apply(jags.obs.Laake[,1,], 
                                   FUN = function(x) c(as.vector(Laake.observers.new$new.ID), 
                                                       x)[match(x, as.vector(Laake.observers.new$ID.old), x)], 
                                   MARGIN = 2)
  
  jags.obs.Laake.new[,2,] <- apply(jags.obs.Laake[,2,], 
                                   FUN = function(x) c(as.vector(Laake.observers.new$new.ID), 
                                                       x)[match(x, as.vector(Laake.observers.new$ID.old), x)], 
                                   MARGIN = 2)
  
  # new observers next - this list may change according to which years are used
  # for a particular run. 
  # 
  # jags input has no-observer ID (No.obs.ID.new) for day 1 and max.day
  jags.obs.new <- jags.input.new$jags.data$obs 
  
  # new.data.no.ID <- new.data.observers.new %>%
  #   filter(obs == "No obs") %>%
  #   pull(ID.old) %>%
  #   as.numeric()
  # 
  jags.obs.new.new <- array(data = no.obs.ID.new,
                            dim = dim(jags.obs.new))
  
  jags.obs.new.new[,1,] <- apply(jags.obs.new[,1,], 
                                 FUN = function(x) c(as.vector(new.data.observers.new$new.ID), 
                                                     x)[match(x, as.vector(new.data.observers.new$ID.old), x)], 
                                 MARGIN = 2)
  
  jags.obs.new.new[,2,] <- apply(jags.obs.new[,2,], 
                                 FUN = function(x) c(as.vector(new.data.observers.new$new.ID), 
                                                     x)[match(x, as.vector(new.data.observers.new$ID.old), x)], 
                                 MARGIN = 2)
  
  # all data arrays need to be combined between Laake and new datasets
  # First make the first dimension the same size between the two arrays:
  n.row.dif <- dim(jags.obs.Laake.new)[1] - dim(jags.obs.new.new)[1]
  
  if (n.row.dif > 0){
    # observers:
    jags.obs.new.new.1 <- abind(jags.obs.new.new,
                                array(data = no.obs.ID.new, 
                                      dim = c(n.row.dif, 2, 
                                              dim(jags.obs.new.new)[3])),
                                along = 1)
    
    # Combine the arrays.
    jags.obs.all <- abind(jags.obs.Laake.new, 
                          jags.obs.new.new.1,
                          along = 3)
    
    

    # day
    jags.day.new <- abind(jags.input.new$jags.data$day,
                          array(data = NA, 
                                dim = c(n.row.dif, 2, 
                                        dim(jags.obs.new.new)[3])),
                          along = 1)
    
    jags.day.all <- abind(jags.input.Laake$jags.data$day,
                          jags.day.new,
                          along = 3)
    
  } else {
    n.row.dif <- abs(n.row.dif)
    jags.obs.Laake.new.1 <- abind(jags.obs.Laake.new,
                                array(data = no.obs.ID.new, 
                                      dim = c(n.row.dif, 2, 
                                              dim(jags.obs.Laake.new)[3])),
                                along = 1)
    
    jags.obs.all <- abind(jags.obs.Laake.new.1, 
                          jags.obs.new.new,
                          along = 3)
    
    
    # day
    jags.day.Laake <- abind(jags.input.Laake$jags.data$day,
                          array(data = NA, 
                                dim = c(n.row.dif, 2, 
                                        dim(jags.obs.Laake.new)[3])),
                          along = 1)
    
    jags.day.all <- abind(jags.day.Laake,
                          jags.input.new$jags.data$day,
                          along = 3) 
  }
  
  # whale counts:
  n.row.dif <- dim(jags.input.Laake$jags.data$n)[1] - dim(jags.input.new$jags.data$n)[1]
  
  if (n.row.dif > 0){

    # whale counts:
    jags.n.new <- abind(jags.input.new$jags.data$n,
                        array(data = NA, 
                              dim = c(n.row.dif, 2, 
                                      dim(jags.input.new$jags.data$n)[3])),
                        along = 1)
    
    jags.n.all <- abind(jags.input.Laake$jags.data$n,
                        jags.n.new,
                        along = 3)
    
  } else {
    n.row.dif <- abs(n.row.dif)
    
    # whale counts:
    jags.n.Laake <- abind(jags.input.Laake$jags.data$n,
                          array(data = NA, 
                                dim = c(n.row.dif, 2, 
                                        dim(jags.input.Laake$jags.data$n)[3])),
                          along = 1)
    
    jags.n.all <- abind(jags.n.Laake,
                        jags.input.new$jags.data$n,
                        along = 3)
    
  }
  
  
  # visibility/Beaufort/Watch length and proportion
  # Visibility and Beaufort do not have extra two rows. So need to recalculate
  # the difference in the numbers of rows
  n.row.dif <- dim(jags.input.Laake$jags.data$vs)[1] - dim(jags.input.new$jags.data$vs)[1]
  if (n.row.dif > 0){
    jags.vs.new <- abind(jags.input.new$jags.data$vs,
                         array(data = NA, 
                               dim = c(n.row.dif, 2, 
                                       dim(jags.obs.new.new)[3])),
                         along = 1)
    
    jags.vs.all <- abind(jags.input.Laake$jags.data$vs,
                         jags.vs.new,
                         along = 3)
    
    # Beaufort
    jags.bf.new <- abind(jags.input.new$jags.data$bf,
                         array(data = NA, 
                               dim = c(n.row.dif, 2, 
                                       dim(jags.obs.new.new)[3])),
                         along = 1)
    
    jags.bf.all <- abind(jags.input.Laake$jags.data$bf,
                         jags.bf.new,
                         along = 3)    
    
    # Watch proportion
    jags.watch.prop.new <- abind(jags.input.new$jags.data$watch.prop,
                                 array(data = NA, 
                                       dim = c(n.row.dif, 2, 
                                               dim(jags.obs.new.new)[3])),
                                 along = 1)
    
    jags.watch.prop.all <- abind(jags.input.Laake$jags.data$watch.prop,
                                 jags.watch.prop.new,
                                 along = 3)
    
    # Watch length
    jags.watch.length.new <- abind(jags.input.new$jags.data$watch.length,
                                   array(data = NA, 
                                         dim = c(n.row.dif, 2, 
                                                 dim(jags.obs.new.new)[3])),
                                   along = 1)
    
    jags.watch.length.all <- abind(jags.input.Laake$jags.data$watch.length,
                                   jags.watch.length.new,
                                   along = 3)
    
  } else {
    n.row.dif <- abs(n.row.dif)
    jags.vs.Laake <- abind(jags.input.Laake$jags.data$vs,
                         array(data = NA, 
                               dim = c(n.row.dif, 2, 
                                       dim(jags.obs.Laake.new)[3])),
                         along = 1)
    
    jags.vs.all <- abind(jags.vs.Laake,
                         jags.input.new$jags.data$vs,
                         along = 3)
    
    # Beaufort
    jags.bf.Laake <- abind(jags.input.Laake$jags.data$bf,
                         array(data = NA, 
                               dim = c(n.row.dif, 2, 
                                       dim(jags.obs.Laake.new)[3])),
                         along = 1)
    
    jags.bf.all <- abind(jags.bf.Laake,
                         jags.input.new$jags.data$bf,
                         along = 3)
    
    # Watch proportion
    jags.watch.prop.Laake <- abind(jags.input.Laake$jags.data$watch.prop,
                                   array(data = NA, 
                                         dim = c(n.row.dif, 2, 
                                                 dim(jags.obs.Laake.new)[3])),
                                   along = 1)
    
    jags.watch.prop.all <- abind(jags.watch.prop.Laake,
                                 jags.input.new$jags.data$watch.prop,
                                 along = 3)
    
    # Watch length
    jags.watch.length.Laake <- abind(jags.input.Laake$jags.data$watch.length,
                                     array(data = NA, 
                                           dim = c(n.row.dif, 2, 
                                                   dim(jags.obs.Laake.new)[3])),
                                     along = 1)
    
    jags.watch.length.all <- abind(jags.watch.length.Laake,
                                   jags.input.new$jags.data$watch.length,
                                   along = 3)
    
  }

  # Create partially observed N.
  N <- matrix(data = NA, nrow = max.day, ncol = dim(jags.n.all)[3])
  N[1,] <- 0
  N[max.day,] <- 0
  
  log_watch <- log(jags.watch.length.all)
  mean_log_watch <- mean(log_watch, na.rm = TRUE)
  
  jags.data <- list(  n = labelled::remove_attributes(jags.n.all, "dimnames"),
                      n.station = c(jags.input.Laake$jags.data$n.station,
                                    jags.input.new$jags.data$n.station),
                      n.year = dim(jags.bf.all)[3],
                      n.obs = nrow(obs.all) - 1,
                      periods = labelled::remove_attributes(rbind(as.matrix(jags.input.Laake$jags.data$periods), 
                                      jags.input.new$jags.data$periods), "dimnames"),
                      n.days = max(jags.input.Laake$jags.data$n.days,
                                   jags.input.new$jags.data$n.days),
                      obs = labelled::remove_attributes(jags.obs.all, "dimnames"),
                      vs = labelled::remove_attributes(jags.vs.all, "dimnames"),
                      bf = labelled::remove_attributes(jags.bf.all, "dimnames"),
                      watch.prop = labelled::remove_attributes(jags.watch.prop.all, "dimnames"),
                      watch.length = labelled::remove_attributes(jags.watch.length.all, "dimnames"),
                      mean.log.watch = mean_log_watch,
                      day = labelled::remove_attributes(jags.day.all, "dimnames"),
                      N = N)
  
  return(list(jags.data = jags.data,
              jags.input.Laake = jags.input.Laake,
              jags.input.new = jags.input.new,
              all.observers = obs.all,
              min.dur = min.dur))
  
}


# Create trace plots without using bayesplot. Sometimes bayesplot complains about
# Error in out[, i, ] <- x[[i]] : 
#number of items to replace is not a multiple of replacement length
# If selecting variable names with brackets, use two backslash escape:
# e.g., to plot all variables with S[, plot.trace.dens("S\\[", jm)

plot.trace.dens <- function(jm, var.name){
  par.names <- unlist(dimnames(jm$samples[[1]])[2])
  
  col.idx <- grep(var.name, par.names, perl = T)    

  samples.list <- list()
  
  n.samples <- nrow(jm$samples[[1]])
  n.chains <- length(jm$samples)
  samples <- lapply(jm$samples, FUN = function(x) x[, col.idx])
  if (length(col.idx) > 1){
    for (k in 1:length(col.idx)){
      samples.list[[k]] <- unlist(lapply(samples, FUN = function(x) x[,k]))
    }
    
    # Order the facet variable numerically and create a factor variable
    # Suggested by Gemini coding partner. numeric index is used to order
    # parameters with square brackets correctly. For example, without this line,
    # 10 becomes before 2. For parameter names without brackets, numeric index
    # is all NAs
    samples.df <- data.frame(seq = rep(1:n.samples, 
                                       times = n.chains),
                             sample = unlist(samples.list),
                             par.name = rep(par.names[col.idx], 
                                            each = length(samples.list[[1]])),
                             chain = rep(1:n.chains, 
                                         each = n.samples)) %>%
      mutate(numeric.index = as.numeric(str_extract(par.name, "(?<=\\[)\\d+(?=\\])"))) %>%
      mutate(par.name.ordered = factor(par.name, levels = unique(par.name[order(numeric.index)])))
    p.trace <- ggplot(samples.df) +
      geom_line(aes(x = seq, y = sample, color = chain)) +
      facet_wrap(~ par.name.ordered) +
      theme(legend.position = "none")
    
    p.dens <- ggplot(samples.df) +
      geom_density(aes(x = sample)) +
      xlab("") +
      facet_wrap(~ par.name.ordered)
  } else {
    samples.vec <- unlist(samples)
    samples.df <- data.frame(seq = rep(1:n.samples),
                             sample = samples.vec,
                             par.name = factor(par.names[col.idx]),
                             chain = rep(1:n.chains, 
                                         times = n.samples))
    p.trace <- ggplot(samples.df) +
      geom_line(aes(x = seq, y = sample, color = chain)) +
      theme(legend.position = "none") +
      theme_bw()
    
    p.dens <- ggplot(samples.df) +
      geom_density(aes(x = sample)) +
      xlab("") +
      theme_bw()
    
  }

 
  
  return(list(df = samples.df,
              p.trace = p.trace,
              p.dens = p.dens))
}

shift.definition <- function(date, time){
  dur.since.midnight <- difftime(paste(date, time), paste(date, "00:00:00"), units = "hours")
  shift.id <- ifelse(dur.since.midnight <= 7.5, 0,
                     ifelse(dur.since.midnight <= 9, 1,
                            ifelse(dur.since.midnight <= 10.5, 2,
                                   ifelse(dur.since.midnight <= 12, 3,
                                          ifelse(dur.since.midnight <= 13.5, 4,
                                                 ifelse(dur.since.midnight <= 15, 5,
                                                        ifelse(dur.since.midnight <= 16.5, 6, 7)))))))
  return(shift.id)
}

shift.definition.2010 <- function(time.dec.hrs){
  
  shift.id <- ifelse(time.dec.hrs <= 7.5, 0,
                     ifelse(time.dec.hrs <= 9, 1,
                            ifelse(time.dec.hrs <= 10.5, 2,
                                   ifelse(time.dec.hrs <= 12, 3,
                                          ifelse(time.dec.hrs <= 13.5, 4,
                                                 ifelse(time.dec.hrs <= 15, 5,
                                                        ifelse(time.dec.hrs <= 16.5, 6, 7)))))))
  return(shift.id)
}

shift.hrs <- data.frame(shift = c(1:6),
                        begin.hr = c(7.5,9,10.5,12,13.5,15),
                        end.hr = c(9,10.5,12,13.5,15,16.5))

compute.LOOIC <- function(loglik.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  n.samples <- dim(loglik.array)[1]
  
  # Create an empty list to hold the log likelihood values
  #ll.list <- vector("list", length = n.samples)
  
  loglik.list <- lapply(1:n.samples, function(s) {
    # Extract the slice for the current sample 's'
    slice <- loglik.array[s, , , ]
    # Return the non-NA values from that slice
    slice[!is.na(slice)]
  })
  
  #loglik.list <- apply(loglik.array, 1, function(slice) slice[!is.na(slice)])
  
  # loop through each MCMC sample and collect log likelihood values
  # for (s in 1:n.samples){
  #   ll.cube <- loglik.array[s,,,]
  #   ll.list[[s]] <- ll.cube[!is.na(loglik.array[s,,,])]
  # }
  
  loglik.mat <- do.call(rbind, loglik.list)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = MCMC.params$n.chains)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = MCMC.params$n.chains, 
                           k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}

# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Richards_fcn <- function(d, S1, S2, K, P, min, max){
  K <- abs(K)
  if (S1 > 0) S1 <- -S1
  if (S2 < 0) S2 <- -S2
  
  M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
  M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}


Richards_fcn_1 <- function(d, S1, S2, K1, K2, P1, P2, min, max){
  K1 <- abs(K1)
  K2 <- abs(K2)
  
  if (S1 > 0) S1 <- -S1
  if (S2 < 0) S2 <- -S2
  
  if (P1 > P2){
    P1.1 <- P1
    P1 <- P2
    P2 <- P1.1
  }
  
  M1 <- (1 + (2 * exp(K1) - 1) * exp((1/S1) * (P1 - d))) ^ (-1/exp(K1))
  M2 <- (1 + (2 * exp(K2) - 1) * exp((1/S2) * (P2 - d))) ^ (-1/exp(K2))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}

# A function to get one data file from selected directory
# Inputs are data directory name, year of survey (2021/2022 is 2022), and
# which file to be extracted (sequential number from 1 to length(files)).
get.data <- function(dir, YEAR, FILES, ff){
  # 2023-03-02 Commented the following line and added FILES input because all data for 2023
  # were combined in one file (EditedDataAll_2023.dat), which was parsed out to day-specific
  # files so they will be the same as other years. The combined file also was stored in the
  # same folder.
  
  all.lines <- read_lines(file = paste0(dir, "/", YEAR, "/", FILES[ff]))
  input.file.name <- FILES[ff]   # specify file name
  
  # look at the first three letters of the first line
  first.3 <- str_sub(all.lines[1], start = 1, end = 3)
  
  if (is.na(as.numeric(first.3)))
    all.lines <- all.lines[2:length(all.lines)]
  
  # look at all event code
  event.code <- str_sub(all.lines, start = 5, end = 5)
  
  # are there any comments?
  COMMENTS <- which(event.code == "C")
  
  # if there were comments, remove them
  if(length(COMMENTS)>0){
    data <- read.table(text = all.lines[-COMMENTS],
                       fill=T,
                       na.strings = "",
                       stringsAsFactors = F,
                       col.names = paste0("V",1:16))
  }else{
    data <- read.table(text = all.lines,
                       fill=T,
                       na.strings = "",
                       stringsAsFactors = F,
                       col.names = paste0("V",1:16))
    
  }
  
  #data <- data[, colSums(is.na(data)) < nrow(data)]
  
  # Files with extensive comments in sightings are problematic because they get split into multiple lines. 
  # We need to pull out lines that contain only numeric V1 (when they are converted into numeric)
  # data %>% 
  #   mutate(line.num = as.numeric(V1)) -> data
  
  ## Inserted on 2026-05-05
  # In 2008 and 2010, time was recorded as fractional hours... so, they need to be converted
  # into HH:MM:SS. 
  if (is.numeric(data$V4)){
    new.time <- fractional_Hr2HMS(data$V4)
    data$V4 <- new.time
  }
  ## End 
  
  data <- data[!is.na(as.numeric(data$V1)),]
  
  Starts <- which(data$V2=="B") #Find all start times
  Ends <- which(data$V2=="E") #Find all end times
  
  # if there is no "E" at the end, Add "E" with time equal to 1 second after
  # the last entry.
  if (length(Ends) == 0 | max(Ends) != nrow(data)){
    row.num <- as.numeric(data[nrow(data), 1])
    row.num.char <- ifelse(row.num < 100, 
                           paste0("0", as.character(row.num+1)),
                           as.character(row.num + 1))
    
    #if (YEAR != 2010){
    tmp <- hour(hms(data[nrow(data),4])) + 
        minute(hms(data[nrow(data),4]))/60 + 
        (second(hms(data[nrow(data),4])) + 1)/3600 
      
    #} else {
    #  tmp <- data[nrow(data), 4]
    #}
    
    HMS <- fractional_Hr2HMS(as.numeric(tmp))
    # h <- trunc(tmp)
    # m <- trunc((tmp - h) * 60)
    # s <- (((tmp - h) * 60) - m) * 60
    data <- rbind(data, c(row.num.char, "E", data[nrow(data), 3],
                          HMS,
                          rep(NA, times = 12)))
  }
  
  if(length(Starts)>0 & length(Ends)>0){
    #Make an array to hold the time differences of starts and ends
    Diffs <- matrix(NA, ncol=length(Ends), nrow=length(Starts)) 
    
    for(t in 1:length(Starts)){
      #Subtract all end times from each start time (if an end time is less than 90 seconds 
      # before a start time, that's probably an error)
      # TE: I added [t] to Ends in the following line. I think it's needed. NO... 
      # Ends does not need the subscript. 
      if (YEAR != 2010){
        Diffs[t,] <- seconds(hms(data[Starts[t],4])) - seconds(hms(data[Ends,4]))         
      } else {
        Diffs[t,] <- (as.numeric(data[Starts[t], 4]) - as.numeric(data[Ends, 4])) * 24 * 60 * 60
      }

    }
    
    #Select the differences that are likely errors (End, <90 seconds, Start. Oops!)
    Oops <- which(Diffs >=0 & Diffs<91,arr.ind = T) 
    
    if(length(Oops)>0){ #If there are any suspect errors, remove them from the data file
      data <- data[-c(Starts[Oops[1]], (Starts[Oops[1]]+1), Ends[Oops[2]]),]
    }
  }
  
  BeginDay <- mdy(data$V3) - mdy(paste0("11/30/", (YEAR - 1)))
  # Decimal hour of shift start time
  BeginHr <- hour(hms(data[, 4])) + 
                    minute(hms(data[, 4]))/60 + 
                    second(hms(data[, 4]))/3600 
					
  #if (YEAR != 2010){
  #  BeginHr <- hour(hms(data$V4)) + minute(hms(data$V4))/60 + second(hms(data$V4))/3600    
  #} else {
  #  BeginHr <- as.numeric(data$V4)
  #}

  data %>% 
    mutate(begin = as.numeric(BeginDay) + BeginHr/24,
           shift = cumsum(V2=="P")) -> data
  
  # People like to enter comments on "E" lines... remove all extra comments:
  data[data$V2 == "E", c("V5", "V6", "V7", "V8", "V9", "V10", 
                         "V11", "V12", "V13", "V14", "V15", "V16")] <- NA
  
  data$ff <- input.file.name
  
  return(data)
}

# 
# A function to extract one shift from a data file. Use get.data first and
# use the output of get.data in this function. i indicates a shift number within
# the day, which can be different from the defined shifts which are identified
# below. The first shift of a survey day may start at any time of the day due 
# to the environmental conditions. 
# 
# Defined shifts are;
# 1: 7:30:01 - 9:00:00
# 2: 9:00:01 - 10:30:00
# 3: 10:30:01 - 12:00:00
# 4: 12:00:01 - 13:30:00
# 5: 13:30:01 - 15:00:00
# 6: 15:00:01 - 16:30:00
# 
# These shifts are indicated by Shift with the capital S. 

get.shift <- function(YEAR, data, i){
  ff <- data$ff[1]   # file name
  
  # Fracitional time in 2010 and 2008 were converted to m/d/Y in the get.data
  # function (Sept 2025). So, the if-else statement is unnecessary.
  #if (YEAR == 2010){
  #  data$Shift <- shift.definition.2010(data$V4)
  #} else {
  data$Shift <- shift.definition(as.Date(data$V3, format = "%m/%d/%Y"), data$V4)
  #}

  # Each shift always begins with "P" - change in observers
  # Note that this can happen in the middle of an official shift... 
  shifts.begin <- which(data$V2 %in% "P")
  shifts.begin.df <- data.frame(event = "P",
                                shift = 1:length(shifts.begin),
                                row = shifts.begin)
  # But each shift does not always have an explicit end. 
  # The end of data file does not always contain "E" either.
  shifts.end <- which(data$V2 %in% "E")
  shifts.end.df <- data.frame(event = "E",
                              shift = NA,
                              row = shifts.end)
  
  shifts.df <- arrange(rbind(shifts.begin.df, shifts.end.df), row)
  
  max.shifts <- length(shifts.begin)
  #Only use the first observer to model random effect
  Observer <- data[shifts.begin[i], 5] %>% toupper()
  # Days since Nov 30th of the previous year because 12/1 is 1. 
  BeginDay <- mdy(data[shifts.begin[i], 3]) - mdy(paste0("11/30/", (YEAR - 1)))
  # Decimal hour of shift start time - need to add seconds because sometimes
  # the last sighting and next shift starts within one minute. This happened
  # in a 2020 data file (file 41, 2020-02-04)
  
  # Beginning hr of the shift. For 2010, time was recorded in decimal hours
  # They are now converted in the m/d/Y format in the get.data function.
  # In some rare occasions, a shift started a few seconds before the expected time
  BeginHr <- hour(hms(data[shifts.begin[i], 4])) + 
    minute(hms(data[shifts.begin[i], 4]))/60 + 
    second(hms(data[shifts.begin[i], 4]))/3600
        
  # Decimal hour of next shift start time
  if (i < max.shifts){
    event.idx <- which(shifts.df$shift %in% i)
    next.event <- shifts.df[event.idx + 1,]
    if (next.event$event == "P"){
      NextBeginHr <- hour(hms(data[next.event$row, 4])) + 
                            minute(hms(data[next.event$row, 4]))/60 + 
                            second(hms(data[next.event$row, 4]))/3600
                     
      
      EndHr <- NextBeginHr - 0.0002
	  
	  #NextBeginHr <- ifelse(YEAR != 2010,
      #                      (hour(hms(data[next.event$row, 4])) + 
      #                         (minute(hms(data[next.event$row, 4]))/60)
      #                       + (second(hms(data[next.event$row, 4]))/3600)),
      #                      as.numeric(data[next.event$row, 4]))
      
      #EndHr <- NextBeginHr - 0.00001
    } else {  # if the event is "E"
      next.P <- shifts.df[event.idx+2,]
      NextBeginHr <- hour(hms(data[next.P$row, 4])) + 
        minute(hms(data[next.P$row, 4]))/60 + 
        second(hms(data[next.P$row, 4]))/3600    
      
      EndHr <- hour(hms(data[next.event$row, 4])) + 
        minute(hms(data[next.event$row, 4]))/60 + 
        second(hms(data[next.event$row, 4]))/3600
		
	  #NextBeginHr <- ifelse(YEAR != 2010, 
      #                      (hour(hms(data[next.P$row, 4])) + 
      #                         (minute(hms(data[next.P$row, 4]))/60) 
      #                       + (second(hms(data[next.P$row, 4]))/3600)),
      #                      as.numeric(data[next.P$row, 4]))
      
      #EndHr <- ifelse(YEAR != 2010,
      #                (hour(hms(data[next.event$row, 4])) + 
      #                   (minute(hms(data[next.event$row, 4]))/60) + 
      #                   (second(hms(data[next.event$row, 4]))/3600)),
      #                as.numeric(data[next.event$row, 4]))
    }

    # Find next end hr to find the next shift to figure out spillovers
    event.idx2 <- which(shifts.df$shift %in% (i+1))
    next.event2 <- shifts.df[event.idx2+1,]
    
    if (next.event2$event == "P"){   # if the next event is also "P"
       NextEndHr <- hour(hms(data[next.event2$row, 4])) + 
         minute(hms(data[next.event2$row, 4]))/60 + 
         second(hms(data[next.event2$row, 4]))/3600 - 0.0002
		 
	   #NextEndHr <- ifelse(YEAR != 2010,
       #                    (hour(hms(data[next.event2$row, 4])) + 
       #                       (minute(hms(data[next.event2$row, 4]))/60) + 
       #                       (second(hms(data[next.event2$row, 4]))/3600) ) - 0.00001,
       #                    as.numeric(data[next.event2$row, 4]) - 0.00001)
       # 
    } else {  # if the event is "E"
       NextEndHr <- hour(hms(data[next.event2$row, 4])) + 
         minute(hms(data[next.event2$row, 4]))/60 + 
         second(hms(data[next.event2$row, 4]))/3600
		 
	   #NextEndHr <- ifelse(YEAR != 2010,
       #                    (hour(hms(data[next.event2$row, 4])) + 
       #                       (minute(hms(data[next.event2$row, 4]))/60) 
       #                     + (second(hms(data[next.event2$row, 4]))/3600)),
       #                    as.numeric(data[next.event2$row, 4]))
    }
    
  } else {    # for the last shift
    event.idx <- which(shifts.df$shift %in% max.shifts)
    next.event <- shifts.df[event.idx + 1,]  # This has to be E
    if (length(next.event) == 0){
      end.row <- nrow(data)
    } else {
      end.row <- next.event$row
    }
    EndHr <- hour(hms(data[end.row, 4])) + 
      minute(hms(data[end.row, 4]))/60 +
      second(hms(data[end.row, 4]))/3600
	  
	#EndHr <-  ifelse(YEAR != 2010,
    #                 (hour(hms(data[end.row, 4])) + 
    #                    (minute(hms(data[end.row, 4]))/60) + 
    #                    (second(hms(data[end.row, 4]))/3600)),
    #                 as.numeric(data[end.row, 4]))
  }
  
  # End time is just before next start time (replicating J Durban's calculations)
  # TE: This is incorrect. If there was an "E", we should use it. 
  #EndHr <- NextBeginHr - 0.00001 
  # Beginning time as a decimal day - NextBeginHr needs to include seconds for the
  # rare occasions when a sighting happens within a minute of the start of a
  # shift. 
  Begin <- as.numeric(BeginDay) + (BeginHr/24) 
  #End time as a decimal day
  End <- as.numeric(BeginDay) + (EndHr/24)
  
  data.shift <- data %>% filter(begin >= Begin & 
                                  begin <= End)
  
  # Remove those that were moving north:
  # Changed V14 != "North" to tolower(V14) != "north" because of entries using
  # uppercase "NORTH" in 2015 (file 11), which was not filtered out correctly in V2. 
  # This was further changed to "grep(north)" because "northbound" or "NORTHBOUND" were
  # used in some years. 
  # 
  # find those that were moving north
  north.idx <- grep("north", tolower(data.shift$V14))
  if (length(north.idx) > 0) data.shift <- data.shift[-north.idx,]
  
  if (i < max.shifts){
    # when there are multiple Es in one file: Take the first of positive values
    if (length(NextBeginHr) > 1){
      dif.BeginHr <- NextBeginHr - BeginHr
      NextBeginHr <- NextBeginHr[dif.BeginHr>0] %>% first()
    }
    
    # following shift for finding out spillovers.
    data.shift2 <- data %>% 
      filter(begin >= as.numeric(BeginDay) + NextBeginHr/24 & 
                                     begin <= as.numeric(BeginDay) + NextEndHr/24)
    
    north.idx <- grep("north", tolower(data.shift2$V14))
    if (length(north.idx) > 0) data.shift2 <- data.shift2[-north.idx,]
    
  } else {
    data.shift2 <- NA
  }

  # This really removes information from "V" entries because there is no
  # V12 when "V" is entered (i.e., NA). I think this is better because
  # there was a case (1/12/2015 9:00 - 10:30) where visibility changed from 4
  # to 5 at the end of a period (30 sec from the end). I think those sightings
  # should be included, rather than excluded. The same is true for VS below. 
  #
  # However... there were shifts (e.g., 2/5/2015 9:00-10:30) where BF changed to
  # 5 in the middle of the shift (1 hr in), so that shift should be removed. But
  # by using just "S" entries, it was not removed. So, I need to fix that. 2022-03-30
  # 
  # This is dealt by having a minimum duration of a shift to be included, which
  # is currently 85 minutes. 
  # 
  # These issues become moot if we use the method of Laake et al., in which data
  # are pooled by a constant continuing viewing condition, rather than arbitrary
  # 90 minute chunks.

  # extract all Beaufort and visibility information from "V" and "S"
  # BF based on "S" entries
  BFs.S <- data.shift %>% 
    filter(V2 == "S") %>% 
    dplyr::select(V12, begin) %>%
    transmute(BF = as.numeric(V12),
              time = begin)
  
  # BFs based on "V" entries
  BFs.V <- data.shift %>% 
    filter(V2 == "V") %>% 
    dplyr::select(V5, begin) %>%
    transmute(BF = as.numeric(V5),
              time = begin)
  
  # Combine them together and find out time from the beginning of the shift
  BFs.dt <- rbind(BFs.S, BFs.V) %>%
    arrange(time) %>%
    mutate(dt = (time - min(time)) * (24*60))  # dt in minutes
  
  # if BF changed to >4 within the last 5 minutes, keep the entire period (max(BF) < 5)
  # but if BF changed to 5 before then, make the max BF = 5.
  # This is a very clumsy way of dealing with changes in the environment. 
  BFs.dt %>% 
    filter(BF > 4) -> High.BFs
  
  if (nrow(High.BFs) == 0) {
    BFs <- BFs.dt$BF
  } else {
    if (min(High.BFs$dt) > 85){   # when the first change happened within 5 min of the shift change
      BFs <- BFs.dt$BF[BFs.dt$BF < 5]
    } else {
      BFs <- BFs.dt$BF
    }
  }
  
  # The entire shift is given the maximum BF value... 
  if (sum(!is.na(BFs)) == 0){
    BF <- NA
  } else {
    BF <- max(BFs, na.rm=T)
  }
  
  # Do the same with visibility conditions (VS)
  VSs.S <- data.shift %>% 
    filter(V2 == "S") %>% 
    dplyr::select(V13, begin) %>%
    transmute(VS = as.numeric(V13),
              time = begin)

  VSs.V <- data.shift %>% 
    filter(V2 == "V") %>% 
    dplyr::select(V6, begin) %>%
    transmute(VS = as.numeric(V6),
              time = begin)
  
  VSs.dt <- rbind(VSs.S, VSs.V) %>%
    arrange(time) %>%
    mutate(dt = (time - min(time)) * (24*60))  # dt in minutes
  
  # if VS changed to 5 within the last 5 minutes, keep the entire period (max(VS) < 5)
  # but if VS changed to 5 before then, make the max VS = 5.
  VSs.dt %>% 
    filter(VS > 4) -> High.VSs
  
  if (nrow(High.VSs) == 0) {
    VSs <- VSs.dt$VS
  } else {
    if (min(High.VSs$dt) > 85){   # when the first change happened within 5 min of the shift change
      VSs <- VSs.dt$VS[VSs.dt$VS < 5]
    } else {
      VSs <- VSs.dt$VS
    }
  }
  
  if (sum(!is.na(VSs)) == 0){
    VS <- NA
  } else {
    VS <- max(VSs, na.rm=T)
  }

  # if still NA, take the first "V" entry
  if (is.na(BF)) {BF <- data[shifts.begin[i]+1, 5]}
  if (is.na(VS)) {VS <- data[shifts.begin[i]+1, 6]}
  
  # Finding groups that were sighted over two shifts. We take the later sighting.
  Spillover <- vector(length = 0)
  # No spillover for the last shift
  if (i < max.shifts){
    # Groups = Observers. Only the first (primary) observer is considered (V5)
    # Group numbers from this watch period
    GroupsThisWatch <- data.shift %>% 
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    
    # Group numbers from next watch period
    GroupsNextWatch <- data.shift2 %>% 
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    
    # Which groups from watch i were also observed in watch i+1? 
    # They should be excluded from i and counted in i+1
    Spillover <- GroupsThisWatch[GroupsThisWatch %in% GroupsNextWatch] 
    
  }
  
  
  # v4 is only time and End is the number of days. So, it doesn't matter what year
  # I use as the starting point. I use 2022-12-01 00:00:00
  # If the last one is not E, add an E line. 
  if (data.shift[nrow(data.shift), "V2"] != "E"){
    #if (YEAR == 2010){
      # This Shift is the defined shift IDs - not based on changes in observers
    #  Shift.End <- shift.definition.2010(data.shift$V4[nrow(data.shift)])
    #  if (NextBeginHr > (BeginHr + 1.5)){
    #    V4 <- shift.hrs %>% 
    #      filter(shift == Shift.End) %>% 
    #      select(end.hr) %>%
    #      pull()        
    #  } else {
    #    V4 <- NextBeginHr - 0.00001
    #	}
    #} else {
    Shift.End <- shift.definition(as.Date(data.shift$V3[nrow(data.shift)], 
                                            format = "%m/%d/%Y"), 
                                    data.shift$V4[nrow(data.shift)])
    if (NextBeginHr > (BeginHr + 1.5)){
      V4 <- format(as.POSIXct(as.Date("2022-12-01 00:00:00") + End),
                  format = "%H:%M:%S")
    } else {
	  V4 <- fractional_Hr2HMS(NextBeginHr - 0.0002)
      #V4 <- NextBeginHr - 0.00001
    }
    #}
    
    # Add one line with "E" as the event
    data.shift <- rbind(data.shift, 
                        data.frame(V1 = max(as.numeric(data.shift$V1), na.rm = T) + 1, 
                                   V2 = "E", 
                                   V3 = data.shift[1, "V3"],
                                   V4 = V4,
                                   V5 = NA,
                                   V6 = NA,
                                   V7 = NA,
                                   V8 = NA,
                                   V9 = NA,
                                   V10 = NA,
                                   V11 = NA,
                                   V12 = NA,
                                   V13 = NA,
                                   V14 = NA,
                                   V15 = NA,
                                   V16 = NA,
                                   begin = End,
                                   shift = i, 
                                   ff = ff,
                                   Shift = Shift.End))
  }

  # Add the "key" variable, which defines a segment with constant environmental 
  # data like visibility and wind force (beaufort) and observer. It is in the format of 
  # Date_shift_ID. ID is the sequential identification number within the shift.
  # Find changes in the viewing condition
  idx.V <- which(data.shift$V2 == "V")
  bft <- c(NA, as.numeric(data.shift$V5[idx.V]))
  vis <- c(NA, as.numeric(data.shift$V6[idx.V]))
  
  # max(idx.V) should be always less than the number of rows of data.shift because
  # I added an "E" row at the end above. 
  idx.V <- c(idx.V, nrow(data.shift))
  
  bft.num <- vis.num <- key.num <- vector(mode = "numeric", length = nrow(data.shift))
  k1 <- 1
  k2 <- 0
  for (k in 1:length(idx.V)){
    key.num[k1:(idx.V[k]-1)] <- k2
    bft.num[k1:(idx.V[k]-1)] <- bft[k]
    vis.num[k1:(idx.V[k]-1)] <- vis[k]
    k1 <- idx.V[k]
    k2 <- k2 + 1
  }
  
  key.num[last(idx.V):length(key.num)] <- max(key.num)
  bft.num[last(idx.V):length(key.num)] <- last(bft)
  vis.num[last(idx.V):length(key.num)] <- last(vis)
  
  data.shift$key <- key.num
  data.shift$effort <- NA
  data.shift$start <- NA
  data.shift$end <- NA
  data.shift$time <- NA
  
  # Compute effort for each block of constant environment and observer
  k <- 1
  for (k in 1:max(key.num)){
    tmp.1 <- data.shift %>%
      filter(key == k) %>%
      summarise(time = first(begin),
                .groups = "drop") %>%
      pull(time) %>% as.numeric()
    
    if (k < max(key.num)){
      tmp.2 <- data.shift %>%
        filter(key == (k + 1)) %>%
        summarise(time = first(begin),
                  .groups = "drop") %>%
        pull(time) %>% as.numeric()
      
    } else {
      tmp.2 <- data.shift$begin[nrow(data.shift)]
    }
    
    data.shift$effort[data.shift$key == k] <- tmp.2 - tmp.1
    data.shift$start[data.shift$key == k] <- tmp.1
    data.shift$end[data.shift$key == k] <- tmp.2
    data.shift$time[data.shift$key == k] <- tmp.1 + (tmp.2-tmp.1)/2
  }
  
  # Need to remove the spillover groups in order to count npods and nwhales
  # if there was a spillover
  if(length(Spillover > 0)){ #if there are groups that spill over into following watch, 
    is.spillover <- T
    # figure out if there were any sightings that need to be considered:
    sub.data <- data.shift %>% 
      filter(V2 == "S", !(V5 %in% Spillover))  
    
    # 2023-11-29, In the following if statement, non-spillover groups are counted.
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        dplyr::select(V5, V9) %>%
        #summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        summarize(N = last(as.numeric(V9))) %>% 
		#summarize(N = last(as.numeric(V9)),
        #          .groups = "drop") %>% 
        dplyr::select(N)  %>% sum()
      npods <- length(unique(sub.data$V5))
      
    } else {
      N <- 0
      npods <- 0
    }
    
  } else {   # if there were no spillover
    is.spillover <- F
    sub.data <- data.shift %>% 
      filter(V2 == "S")
    
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        dplyr::select(V5, V9) %>%
        #summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        summarize(N = last(as.numeric(V9))) %>% 
		#summarize(N = last(as.numeric(V9)),
        #          .groups = "drop") %>% 
        dplyr::select(N)  %>% sum()
      npods <- length(unique(sub.data$V5))
      
    } else {
      N <- 0
      npods <- 0
    }
    
  }
  
  # shift and Shift are the same if the first shift of a day started before 0900 
  # and observations didn't stop until 1630. 
  # key is the sequential number within each shift that defines an equal 
  # environmental condition.
  
  # When there were at least one sighting
  #if (length(which(data.shift$V2 == "S")) != 0){
  
  # There were at least one spillover to the next shift and at least one group
  # was recorded within the shift other than the spilled over groups
  #if (is.spillover & nrow(sub.data) > 0){
  if (nrow(sub.data) > 0){
    # sightings summary
    sub.data %>%
      transmute(Date = V3, 
                Time = V4, 
                Group_ID = as.numeric(V5), 
                n = as.numeric(V9), 
                bft = as.numeric(V12), 
                vis = as.numeric(V13),
                Bearing = as.numeric(V6),
                Reticle = as.numeric(V7),
                Distance = as.numeric(V8),
                Observer = V10,
                shift = shift, 
                key = key, 
                begin = start,
                end = end,
                time = time,
                effort = effort,
                Shift = Shift) %>%
      group_by(Group_ID) %>%
      summarise(Date = first(Date),
                Time = first(Time),
                n = max(n, na.rm = T),
                bft = first(bft),
                vis = first(vis),
                Bearing = first(Bearing[n == max(n)]),
                Reticle = first(Reticle[n == max(n)]),
                Distance = first(Distance[n == max(n)]),
                Observer = first(Observer),
                shift = first(shift),
                key = first(key),
                begin = first(begin),
                end = first(end),
                time = first(time),
                effort = first(effort),
                Shift = first(Shift)) -> sub.data.shift
                #.groups = "drop") -> sub.data.shift
    
    # effort summary
    # Rare occasions when observer changes within a shift... this needs to be
    # separated because "observer" is a covariate
    
    sub.data.shift %>%
      filter(key > 0) %>%
      group_by(key, Observer) %>%
      summarise(Date = first(Date),
                npods = n(),
                nwhales = sum(n, na.rm = T),
                bft = first(bft),
                vis = first(vis),
                shift = first(shift),
                Observer = first(Observer),
                key = first(key),
                begin = first(begin),
                end = first(end),
                time = first(time),
                effort = first(effort),
                Shift = first(Shift))  -> data.shift.effort
                #.groups = "drop") -> data.shift.effort
  
  } else {
    # There were only spillover sightings, which means no sightings were
    # recorded for this shift, meaning (is.spillover & nrow(sub.data) == 0)
    # or !is.spillover & nrow(sub.data == 0). All differing sighting conditions have
    # to be separated
    data.shift %>%
      filter(key > 0, V2 != "E") %>%
      transmute(Group_ID = NA,
                Date = V3, 
                Time = V4, 
                n = 0, 
                bft = NA, 
                vis = NA,
                Bearing = NA,
                Reticle = NA,
                Distance = NA,
                Observer = Observer,
                shift = shift, 
                key = key, 
                begin = start,
                end = end,
                time = time,
                effort = effort,
                Shift = Shift) -> sub.data.shift
    
    # fix beaufort and visibility 
    for (k1 in 1:nrow(BFs.dt)){
      sub.data.shift$bft[sub.data.shift$begin == BFs.dt$time[k1]] <- BFs.dt$BF[k1]
      sub.data.shift$vis[sub.data.shift$begin == VSs.dt$time[k1]] <- VSs.dt$VS[k1]      
    }
    
    sub.data.shift %>%
      group_by(key, Observer) %>%
      summarise(Date = first(Date),
                npods = 0,
                nwhales = 0,
                bft = first(bft),
                vis = first(vis),
                shift = first(shift),
                Observer = first(Observer),
                key = first(key),
                begin = first(begin),
                end = first(end),
                time = first(time),
                effort = first(effort),
                Shift = first(Shift)) -> data.shift.effort
                #.groups = "drop") -> data.shift.effort
    
  }
  
  out.list <- list(out.df = data.frame(begin = as.numeric(Begin),
                                       end = as.numeric(End),
                                       dur = as.numeric(End) - as.numeric(Begin),
                                       bf = as.numeric(BF),
                                       vs = as.numeric(VS),
                                       n = N,
                                       npods = npods,
                                       obs = as.character(Observer),
                                       ff = ff,
                                       i = i,
                                       BeginHr = BeginHr,
                                       BeginDay = BeginDay),
                   data = sub.data.shift,
                   data.shift = data.shift,
                   shift.effort = data.shift.effort,
                   data.next.shift = data.shift2)
  return( out.list )
}

# converts a fractional date into YMD and hms
# The second input should be the second year
# of a season. For example, for the 1978/1979
# season, it should be 1979. 
fractional_Day2YMDhms <- function(x, YEAR){
  n.days <- floor(x)
  dec.hr <- (x - n.days) * 24
  hr <- floor(dec.hr)
  dec.min <- (dec.hr - hr) * 60
  m <- floor(dec.min)
  s <- floor((dec.min - m) * 60)
  

  mdy <- n.days + as.Date(paste0((YEAR-1), "-11-30"))          

  return(list(YMD = mdy,
              hms = paste(ifelse(hr < 10, paste0("0", hr), hr), 
                          ifelse(m < 10, paste0("0", m), m), 
                          ifelse(s < 10, paste0("0", s), s), sep = ":")))  
}

# Converts fractional hour to the H:M:S format.
fractional_Hr2HMS <- function(tmp){
  tmp[tmp>24] <- tmp[tmp>24] %% 24
  
  h <- trunc(tmp)
  m <- trunc((tmp - h) * 60)
  s <- round((((tmp - h) * 60) - m) * 60)
  
  HMS <- paste(ifelse(h < 10, paste0("0", h), h), 
               ifelse(m < 10, paste0("0", m), m), 
               ifelse(s < 10, paste0("0", s), s), sep = ":")
  
  return(HMS)

}

# This function compares Ver1.0 and Ver2.0 data extraction code for using
# raw data files (starting the 2020 season). The raw data files should be 
# analyzed using Formatting GC Data TE.R for Ver1.0 (saves in a .rds file)
# and Extract_Data_All_v2.Rmd for Ver2.0 (saves in a .rds file).
compare.V0.V2.raw <- function(YEAR, obs.list){
  v0.out <- readRDS(paste0("RData/out_", YEAR, "_Joshs.rds"))
  v2.out <- readRDS(paste0("RData/V2.1_Sep2023/out_", YEAR, "_min85_Tomo_v2.rds"))
  
  v2.out$Final_Data %>% 
    mutate(v = "V2") %>% # -> tmp
    dplyr::select(-dur) %>% 
    left_join(obs.list, by = "obs") -> FinalData.v2
  
  FinalData.v0 <- v0.out$FinalData %>% mutate(v = "V0") 
  
  # find if there is NA in ID - not in the look up table  
  ID.NA <- filter(FinalData.v2, is.na(ID))
  
  unique.ID.NA <- unique(ID.NA$obs)

  if (length(unique.ID.NA) > 0){
    new.obs <- data.frame(obs = NA, ID = NA)
    for (k in 1:length(unique.ID.NA)){
      FinalData.v2$ID[FinalData.v2$obs == unique.ID.NA[k]] <- max(obs.list$ID) + k
      new.obs[k,] <- c(unique.ID.NA[k], max(obs.list$ID)+k)
    }
    obs.list <- rbind(obs.list, new.obs)
    
  }
  
  
  # replace column names
  FinalData.v2 %>% dplyr::select(-obs) %>%
    mutate(obs = ID) %>%
    dplyr::select(-ID) -> FinalData.v2
  
  # rearrange the columns to match v0
  FinalData.v2 <- FinalData.v2[, names(FinalData.v0)]
  FinalData.Both <- rbind(FinalData.v2, FinalData.v0)
  
  min.begin <- min(floor(FinalData.Both$begin))
  max.begin <- max(ceiling(FinalData.Both$begin))
  
  time.steps <- min.begin:max.begin
  difs <- data.frame(begin = double(),
                     end = double(),
                     min.begin = double(), 
                     max.end = double(), 
                     n.periods = integer(), 
                     max.bf = integer(), 
                     max.vs = integer(), 
                     total.whales = integer(),
                     time.step = integer(),
                     stringsAsFactors = F)
  
  c <- k <- 1
  for (k in 1:(length(time.steps)-1)){
    tmp <- filter(FinalData.Both, begin >= time.steps[k] & begin < time.steps[k+1])
    if (nrow(tmp) > 0){
      
      tmp %>% filter(v == "V0") -> tmp.1
      tmp %>% filter(v == "V2") -> tmp.2
      
      if (nrow(tmp.1) > 0 & nrow(tmp.2) > 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      min(tmp.1$begin) - min(tmp.2$begin), 
                      max(tmp.1$end) - max(tmp.2$end),
                      nrow(tmp.1) - nrow(tmp.2),
                      max(tmp.1$bf) - max(tmp.1$bf),
                      max(tmp.1$vs) - max(tmp.1$vs),
                      sum(tmp.1$n) - sum(tmp.2$n),
                      time.steps[k])
        
        
      } else if (nrow(tmp.1) > 0 & nrow(tmp.2) == 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      NA, NA, NA, 
                      max(tmp.1$bf) - max(tmp.1$bf),
                      max(tmp.1$vs) - max(tmp.1$vs),
                      NA, time.steps[k])
      } else if (nrow(tmp.1) == 0 & nrow(tmp.2) > 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      NA, 
                      NA,
                      NA,
                      NA,
                      NA,
                      NA,
                      time.steps[k])
        
      }
      
      c <- c + 1
      
    }
    
  }  
  
  difs %>% filter(n.periods != 0 | total.whales != 0) -> difs.1
  FinalData.Both %>% mutate(time.steps = floor(FinalData.Both$begin)) -> FinalData.Both
  
  v2.out$Data_Out %>% 
    mutate(time.steps = floor(v2.out$Data_Out$begin)) -> Data_Out.v2 
  
  v2.out$Correct_Length %>%
    mutate(time.steps = floor(v2.out$Correct_Length$begin)) -> CorrectLength.v2 
  
  return(out.list <- list(difs = difs,
                          difs.1 = difs.1,
                          FinalData.Both = FinalData.Both,
                          Data_Out.v2 = Data_Out.v2,
                          CorrectLength.v2 = CorrectLength.v2,
                          v0.out = v0.out,
                          v2.out = v2.out,
                          obs.list = obs.list))
}

# This function compares outputs from data extraction codes. It uses BUGS input
# data (for data before the 2020 season because the old version (Ver1.0) does not
# work for old files) and Ver2.0. The raw data files should be 
# analyzed using Extract_Data_All_v2.Rmd for Ver2.0 (saves in a .rds file).
compare.V0.V2.BUGSinput <- function(YEAR, idx.yr, periods, obs.list){
  
  #Watch start times, as fraction of a day - stored in a different file
  begin <- as.matrix(read.table("Data/begin.txt", 
                                header=T, 
                                nrows = max(periods)))
  
  #watch end times
  end <- as.matrix(read.table("Data/end.txt", 
                              header=T,
                              nrows = max(periods)))
  
  
  # this file contains all input data for WinBUGS.
  V0.out <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")
  
  # Pull out the information for V0 dataset
  periods.V0 <- V0.out$periods[idx.yr]
  n.V0 <- V0.out$n[1:periods.V0,,idx.yr]
  n.com.V0 <- V0.out$n.com[1:periods.V0,,idx.yr]
  n.sp.V0 <- V0.out$n.sp[1:periods.V0,,idx.yr]
  obs.V0 <- V0.out$obs[1:periods.V0,,idx.yr]
  
  vs.V0 <- V0.out$vs[1:periods.V0,idx.yr]
  bf.V0 <- V0.out$bf[1:periods.V0,idx.yr]
  day.V0 <- V0.out$day[1:periods.V0,idx.yr]
  
  FinalData.V0 <- data.frame(begin = begin[1:periods[idx.yr], idx.yr],
                             end = end[1:periods[idx.yr], idx.yr],
                             bf = bf.V0,
                             vs = vs.V0,
                             n = n.V0[,1],
                             obs = obs.V0[,1],
                             BeginDay = day.V0,
                             v = "V0")
  
  # This contains the results from my version
  v2.out <- readRDS(paste0("RData/V2.1_Sep2023/out_", YEAR, "_min85_Tomo_v2.rds"))
  FinalData.v2 <- v2.out$Final_Data %>% 
    mutate(v = "V2") %>% 
    left_join(obs.list, by = "obs") %>%
    dplyr::select(-c(dur, ff, i, BeginHr)) 
  
  # find if there is NA in ID - not in the look up table  
  ID.NA <- filter(FinalData.v2, is.na(ID))
  
  unique.ID.NA <- unique(ID.NA$obs)
  
  if (length(unique.ID.NA) > 0){
    new.obs <- data.frame(obs = NA, ID = NA)
    
    for (k in 1:length(unique.ID.NA)){
      FinalData.v2[FinalData.v2$obs == unique.ID.NA[k], "ID"] <- max(obs.list$ID) + k
      new.obs[k,] <- c(unique.ID.NA[k], as.numeric(max(obs.list$ID)+k))
    }
    obs.list <- rbind(obs.list, new.obs)
    
  }
  
  # replace column names
  FinalData.v2 %>% 
    dplyr::select(-obs) %>%
    mutate(obs = ID) %>%
    dplyr::select(-ID) -> FinalData.v2
  
  # rearrange the columns to match V0
  FinalData.v2 <- FinalData.v2[, names(FinalData.V0)]
  FinalData.Both <- rbind(FinalData.v2, FinalData.V0)
  
  v2.out$Data_Out %>% 
    mutate(time.steps = floor(v2.out$Data_Out$begin)) -> Data_Out.v2 
  
  v2.out$Correct_Length %>%
    mutate(time.steps = floor(v2.out$Correct_Length$begin)) -> CorrectLength.v2 
  
  min.begin <- min(floor(FinalData.Both$begin))
  max.begin <- max(ceiling(FinalData.Both$begin))
  
  time.steps <- min.begin:max.begin
  difs <- data.frame(begin = double(),
                     end = double(),
                     min.begin = double(), 
                     max.end = double(), 
                     n.periods = integer(), 
                     max.bf = integer(), 
                     max.vs = integer(), 
                     total.whales = integer(),
                     time.step = integer(),
                     stringsAsFactors = F)
  
  c <- k <- 1
  for (k in 1:(length(time.steps)-1)){
    tmp <- filter(FinalData.Both, begin >= time.steps[k] & begin < time.steps[k+1])
    
    if (nrow(tmp) > 0){
  
      tmp %>% filter(v == "V0") -> tmp.1
      tmp %>% filter(v == "V2") -> tmp.2
      
      if (nrow(tmp.1) > 0 & nrow(tmp.2) > 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      min(tmp.1$begin) - min(tmp.2$begin), 
                      max(tmp.1$end) - max(tmp.2$end),
                      nrow(tmp.1) - nrow(tmp.2),
                      max(tmp.1$bf) - max(tmp.1$bf),
                      max(tmp.1$vs) - max(tmp.1$vs),
                      sum(tmp.1$n) - sum(tmp.2$n),
                      time.steps[k])

        
      } else if (nrow(tmp.1) > 0 & nrow(tmp.2) == 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      NA, NA, NA, 
                      max(tmp.1$bf) - max(tmp.1$bf),
                      max(tmp.1$vs) - max(tmp.1$vs),
                      NA, time.steps[k])
      } else if (nrow(tmp.1) == 0 & nrow(tmp.2) > 0){
        difs[c,] <- c(min(tmp$begin), 
                      max(tmp$end),
                      NA, 
                      NA,
                      NA,
                      NA,
                      NA,
                      NA,
                      time.steps[k])
        
      }
      c <- c + 1
      
    }
    #Sys.sleep(1.5)  
  }
  
  
  difs %>% filter(n.periods != 0 | total.whales != 0) -> difs.1
  FinalData.Both %>% mutate(time.steps = floor(FinalData.Both$begin)) -> FinalData.Both
  
  return(out.list <- list(difs = difs,
                          difs.1 = difs.1,
                          FinalData.Both = FinalData.Both,
                          FinalData.V0 = FinalData.V0,
                          Data_Out.v2 = Data_Out.v2,
                          CorrectLength.v2 = CorrectLength.v2,
                          v0.out = V0.out,
                          v2.out = v2.out,
                          obs.list = obs.list))
  
}

# This function compares whale counts, Beaufort, and visibility for
# all shifts that were recorded by Ver1.0 and Ver2.0 and creates
# a table that is sorted by the beginning time of each shift. 
# The table is returned. 2022-04-01
n.comparison <- function(FinalData.Both, difs.1, idx, YEAR){
  FinalData.Both %>% 
    filter(time.steps == difs.1[idx, "time.step"]) %>%
    mutate(begin.time = fractional_Day2YMDhms(begin, YEAR)$hms,
           end.time = fractional_Day2YMDhms(end, YEAR)$hms) %>% 
    dplyr::select(begin.time, end.time, n, bf, vs, v) -> tmp
  
  tmp %>%
    filter(v == "V0") -> tmp.0
  tmp %>%
    filter(v == "V2") -> tmp.2
  
  tmp.0 %>% 
    full_join(tmp.2, by = "begin.time") %>%
    arrange(begin.time) %>%
    transmute(begin.time = begin.time,
              end.time.Ver1.0 = end.time.x,
              end.time.Ver2.0 = end.time.y,
              n.Ver1.0 = n.x,
              n.Ver2.0 = n.y,
              vs.Ver1.0 = vs.x,
              vs.Ver2.0 = vs.y,
              Bf.Ver1.0 = bf.x,
              Bf.Ver2.0 = bf.y) -> tmp.0.2
  
  total <- c(NA, NA, "Total", 
             sum(tmp.0.2$n.Ver1.0, na.rm = T), 
             sum(tmp.0.2$n.Ver2.0, na.rm = T), 
             NA, NA, NA, NA)
  
  return(rbind(tmp.0.2, total))  
}


