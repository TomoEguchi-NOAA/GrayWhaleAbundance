
# define some functions

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
                                 data.set = "2007to2024")
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


# Runs Richards' function with Pois-Bino model on data since 2010 without
# using the WinBUGS input list, i.e., creating input data from output of 
# Extract_Data_All_v2.Rmd. Requires to provide the minimum watch duration,
# a version of the model (see Jags models for details and differences - they are
# different in which parameters of Richards' funcion are time specific),
# years for which the model is fit, a directory name where data are stored,
# parameters to be monitored, and MCMC parameters as a list. 
# 
# Results are saved in a subdirectory under the current working directory named "RData."
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
# Jags_Richards_Since2010_fcn(min.dur = 30, 
#                             ver = "v3", 
#                             years = c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024), 
#                             data.dir = "RData/V2.1_Nov2024", 
#                             jags.params = jags.params, 
#                             MCMC.params = MCMC.params)                    

Jags_Richards_Since2010_fcn <- function(min.dur, ver, years, data.dir, jags.params, MCMC.params){
  print("Starting Jags_Richards_Since2010_fcn")
  
  Run.date <- Sys.Date() #"2024-12-05" #
  
  # Minimum length of observation periods in minutes
  model.name <- paste0("Richards_pois_bino_", ver) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                          "_Since2010_NoBUGS_",
                          Run.date, ".rds")
  
  jags.input <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                       years = years,
                                       data.dir = data.dir)
  
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Run.date,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  

}

# Runs a Richards' function with Pois-Binom model on datasets without using
# WinBUGS input, i.e., creating input data from output of Extract_Data_All_v2.Rmd.
# It also uses Laake's data.
# 
# Requires to provide the minimum watch duration, a version of the model (see 
# Jags models for details and differences - they are
# different in which parameters of Richards' funcion are time specific),
# years for which the model is fit, a directory name where data are stored,
# parameters to be monitored, and MCMC parameters as a list.
# 
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
#                     MCMC.params = MCMC.params)     
#                     
NoBUGS_Richards_fcn <- function(min.dur, ver, years, data.dir, jags.params, MCMC.params){
  print("Starting NoBUGS_Richards_fcn")
  
  Run.date <- Sys.Date()
  model.name <- paste0("Richards_pois_bino_", ver) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                          "_NoBUGS_",
                          Run.date, ".rds")
  
  jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir)                        

  jags.input <- list(jags.data = jags.input.list$jags.data,
                     min.dur = min.dur, 
                     jags.input.Laake = jags.input.list$jags.input.Laake,
                     jags.input.new = jags.input.list$jags.input.new,
                     data.dir = data.dir)
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input.list$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Run.date,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
  
}


# Runs a Richards' function with Pois-Binom model on datasets for data since 2006. 
# It does not use an input list of WinBUGS from an output file from a WinBUGS run. Data
# should be extracted using Extract_Data_All_v2.Rmd. 
# Requires to provide the minimum watch duration, a version of the model (see 
# Jags models for details and differences - they are
# different in which parameters of Richards' funcion are time specific), 
# a directory name where data are stored, parameters to be monitored, and MCMC parameters as a list.
# 
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
# Jags_Richards_NoLaakeData_fcn(min.dur = 85, 
#                               ver = "v4", 
#                               years = c(2007, 2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024), 
#                               data.dir = "RData/V2.1_Nov2024", 
#                               jags.params, MCMC.params)
#                     
Jags_Richards_NoLaakeData_fcn <- function(min.dur, ver, years, data.dir, jags.params, MCMC.params){
  print("Starting Jags_Richards_NoLaakeData_fcn")
  
  Run.date <- Sys.Date()
  model.name <- paste0("Richards_pois_bino_", ver) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name,
                          "_min", min.dur,
                          "_Since2006_",
                          Run.date, ".rds")
  
  jags.input<- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                      years = years,
                                      data.dir = data.dir)
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Run.date,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
  
}


# Runs a Richards' function with Pois-Binom model on datasets for data since 2006. 
# It uses an input list of WinBUGS from an output file from a WinBUGS run. Data
# for other years should be extracted using Extract_Data_All_v2.Rmd. 
# Requires to provide the minimum watch duration, a version of the model (see 
# Jags models for details and differences - they are
# different in which parameters of Richards' funcion are time specific), the name
# and location of WinBUGS output file (a string), 
# years for which the model is fit, the number of station per year (a vector), 
# a directory name where data are stored, parameters to be monitored, and MCMC parameters as a list.
# 
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
# Jags_Richards_AllData_fcn(min.dur = 85, 
#                           ver = "v4", 
#                           WinBUGS.outfile = "RData/WinBUGS_2007to2024_v2_min85_85000.rds", 
#                           WinBUGS.years = c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024), 
#                           WinBUGS.n.stations = c(1, 1, 2, 2, rep(1, times = 6)), 
#                           data.dir = "RData/V2.1_Nov2024", 
#                           jags.params, MCMC.params)
#                           
# "AllData" is not a good name for this because it only runs years that were
# analyzed using WinBUGS.    
Jags_Richards_AllData_fcn <- function(min.dur, ver, WinBUGS.out.file, WinBUGS.years, WinBUGS.n.stations, data.dir, jags.params, MCMC.params){
  print("Starting Jags_Richards_AllData_fcn")
  
  Run.date <- Sys.Date()
  model.name <- paste0("Richards_pois_bino_", ver) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                          "_AllYears_",
                          Run.date, ".rds")
  
  # WinBUGS.years <- c(2010, 2011, 2015, 2016, 
  #                    2020, 2022, 2023, 2024)
  all.years <- c(2007, 2008, WinBUGS.years)
  
  jags.input <- AllData2JagsInput(min.dur = min.dur, 
                                  WinBUGS.years = WinBUGS.years, 
                                  WinBUGS.n.stations = WinBUGS.n.stations, 
                                  WinBUGS.out.file = WinBUGS.out.file,
                                  data.dir)
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Run.date,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
}


# Runs a Richards' function with Pois-Binom model on datasets for data in Laake's 
# analysis (ERAnalysis). 
# Requires to provide the minimum watch duration, a version of the model (see 
# Jags models for details and differences - they are
# different in which parameters of Richards' funcion are time specific), 
# 
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
# Jags_Richards_LaakeData_fcn(min.dur = 85, 
#                             ver = "v4", 
#                             jags.params, MCMC.params)
#                           
#    
Jags_Richards_LaakeData_fcn <- function(min.dur, ver, jags.params, MCMC.params){
  print("Starting Jags_Richards_LaakeData_fcn")
  
  Run.date <- Sys.Date()
  
  model.name <- paste0("Richards_pois_bino_", ver) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                          "_LaakeData_",
                          Run.date, ".rds")
  
  jags.input.Laake <- LaakeData2JagsInput(min.dur)
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input.Laake$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input.Laake,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Run.date,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
  
}


# Creating Jags input list from WinBUGS input, which does not include Laake's 
# data. WinBUGS has to be run first and results saved in a .rds file using
# WinBUGS Ver2.Rmd. 
# min.dur is the minimum duration of a watch period to be 
# included in the analysis. 
# years in this function refers to the years for which raw data were available
# and Extract_Data_All_v2.Rmd was run. The 2006/2007 and 2007/2008 data were not
# available. c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024)
# n.stations is the number of watch stations per year. There were two years (2009/2010
# and 2010/2011) when two independent stations were used. So for the first 10 years,
# this should be n.stations =  = c(1, 1, 2, 2, rep(1, times = 6)) - this input is
# now extracted from existing WinBUGS input (via data2WinBUGS_input)
# data.dir refers to the output directory of Extract_Data_All_v2.Rmd
# 
WinBUGSinputSince2006toJagsInput <- function(min.dur, 
                                             WinBUGS.out.file = "RData/WinBUGS_2007to2024_v2_min85_2024-11-23.rds",
                                             years,
                                             data.dir){

  # e.g., 2007 refers to 2006/2007
  all.years <- c(2007, 2008, years)
  seasons <- sapply(all.years, 
                    FUN = function(x) paste0(x-1, "/", x))
  
  start.years <- all.years - 1
  
  # v2 refers to v2 data extraction. 
  WinBUGS.out <- readRDS(WinBUGS.out.file)
  
  # Create WinBUGS input from raw data files.
  WinBUGS.inputs <- data2WinBUGS_input(data.dir = data.dir,
                                       years = years,
                                       min.duration = min.dur)
  
  data.WinBUGS <- WinBUGS.inputs$data
  
  # Beaufor and visibility are assumed equal between primary and secondary stations in
  # WinBUGS code. But not in Jags. So, I duplicate the secondary watch effort
  bf <- vs <- array(dim = c(max(data.WinBUGS$periods),
                            2, 
                            length(data.WinBUGS$periods)))
  
  # watch lengths are assumed equal between primary and secondary stations in
  # WinBUGS code. But not in Jags. So, I duplicate the secondary watch effort
  watch.prop <- day <- array(dim = c(dim(data.WinBUGS$Watch.Length)[1],
                                     2, 
                                     dim(data.WinBUGS$Watch.Length)[2]))
  
  bf[,1,] <- data.WinBUGS$bf
  bf[,2,] <- data.WinBUGS$bf
  
  vs[,1,] <- data.WinBUGS$vs
  vs[,2,] <- data.WinBUGS$vs
  
  day[,1,] <- data.WinBUGS$day
  day[,2,] <- data.WinBUGS$day
  
  watch.prop[,1,] <- data.WinBUGS$Watch.Length
  watch.prop[,2,] <- data.WinBUGS$Watch.Length
  
  jags.data <- list(  n = data.WinBUGS$n,
                      n.station = data.WinBUGS$n.station,
                      n.year = length(seasons),
                      n.obs = data.WinBUGS$n.obs,
                      #Daily.N = daily.N,
                      periods = cbind(data.WinBUGS$periods,
                                      data.WinBUGS$periods),
                      n.days = max(day, na.rm = T),
                      #first.day = unlist(as.vector(first.day)),
                      obs = data.WinBUGS$obs,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.prop,
                      #watch.prop = (data.WinBUGS$Watch.Length*24*60)/540,
                      day = day)
  
  out.list <- list(jags.data = jags.data,
                   min.dur = min.dur, 
                   seasons = seasons, 
                   WinBUGS.out.file = WinBUGS.out.file,
                   years = years,
                   start.years = start.years,
                   data.dir = data.dir)
  return(out.list)
}

# # Creates WinBUGS input lists from raw data plus 2006/2007 and 2007/2008 from
# # original data from the Durban era.
# data2WinBUGS_input_v2 <- function(data.dir, years, min.dur){
#   library(abind)
#   library(tidyverse)
#   
#   # this file contains all necessary inputs for 2006 - 2019:
#   data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")
#   
#   # e.g., 2007 refers to 2006/2007
#   all.years <- c(2007, 2008, years)
#   seasons <- sapply(all.years, 
#                     FUN = function(x) paste0(x-1, "/", x))
#   
#   # These are extracted data files (Extract_Data_All_v2.Rmd)
#   out.v2 <- lapply(years, 
#                    FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
#                                                     "_min", min.dur, 
#                                                     "_Tomo_v2.rds")))
# 
#   # Use the jags input and convert it to BUGS input:
#   jags.input <- data2Jags_input_NoBUGS(min.dur = min.dur,
#                                        years = years,
#                                        data.dir = data.dir)
#   
#   BUGS.data <- list(n = n[1:max(periods[1:x]),,1:x],
#                     n.com = n[1:max(periods[1:x]),,1:x],
#                     n.sp = n[1:max(periods[1:x]),,1:x],
#                     n.station = dim(n[1:max(periods[1:x]),,1:x])[2],
#                     n.year = dim(n[1:max(periods[1:x]),,1:x])[3],
#                     n.obs = max(obs[1:max(periods[1:x]),,1:x], na.rm = T),
#                     periods = periods[1:x],
#                     obs = obs[1:max(periods[1:x]),,1:x],
#                     #Watch.Length = 0.0625,
#                     u = u[1:max(periods[1:x]),,1:x],
#                     vs = vs[1:max(periods[1:x]),1:x],
#                     bf = bf[1:max(periods[1:x]),1:x],
#                     #day=day,
#                     day = t[1:(max(periods[1:x])+2),1:x],
#                     N = N[,1:x],
#                     N.com = N[,1:x],
#                     N.sp = N[,1:x],
#                     knot = c(-1.46,-1.26,-1.02,-0.78,
#                              -0.58,-0.34,-0.10,0.10,
#                              0.34,0.57,0.78,1.02,1.26,1.46),
#                     n.knots=14,
#                     #begin=begin,
#                     #end=end,
#                     Watch.Length=Watch.Length[1:(max(periods[1:x])+2), 1:x])
#   
#   BUGS.inits <- function() list(mean.prob = 0.5,
#                                 BF.Fixed = 0,
#                                 VS.Fixed = 0,
#                                 mean.prob.sp = 0.5,
#                                 BF.Fixed.sp = 0,
#                                 VS.Fixed.sp = 0,
#                                 mean.prob.com = 0.5,
#                                 BF.Fixed.com = 0,
#                                 VS.Fixed.com = 0,
#                                 mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
#                                 beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
#                                 BF.Switch = 1,
#                                 VS.Switch = 1,
#                                 OBS.Switch = 1,
#                                 sigma.Obs = 1,
#                                 BF.Switch.sp = 1,
#                                 VS.Switch.sp = 1,
#                                 OBS.Switch.sp = 1,
#                                 sigma.Obs.sp = 1,
#                                 BF.Switch.com = 1,
#                                 VS.Switch.com = 1,
#                                 OBS.Switch.com = 1,
#                                 sigma.Obs.com = 1,
#                                 N = N_inits,
#                                 N.com = N_inits,
#                                 N.sp = N_inits,
#                                 #z = matrix(1,nrow=90,ncol=6),
#                                 beta.sp = array(data=0, dim=c(2,x)),
#                                 sd.b.sp = rep(1, times = x), #c(1,1,1,1,1,1),
#                                 z = matrix(1, nrow=90, ncol= x))
#   
#   return(out.list <- list(data = BUGS.data,
#                           inits = BUGS.inits,
#                           all.years = all.years,
#                           seasons = seasons))  
#   
# }


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

    # if the new saeson has less rows than previous maximum
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
                          min.dur = min.dur))  
}

# Create Jags input for Laake's data.
# min.dur is the minimum effort duration in minutes to be included in the analysis
LaakeData2JagsInput <- function(min.dur){
  
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
  load("Data/Primary.rda")      # on-effort sightings
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
  
  # Filter to minimum duration:
  Laake_PrimaryEffort %>%
    filter(effort >= min.dur / (60 * 24) ) -> Laake_PrimaryEffort
  
  Laake_SecondaryEffort %>%
    filter(effort >= min.dur / (60 * 24) ) -> Laake_SecondaryEffort
  
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
    select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
    group_by(Start.year) %>%
    mutate(effort.min = effort * 24 * 60,
           watch.prop = effort) -> Effort.by.period
  
  # observers
  obs <- c(unique(Laake_PrimaryEffort$Observer),
           unique(Laake_SecondaryEffort$Observer)) %>% unique()
  
  Laake.obs <- data.frame(obs = obs,
                          obs.ID = seq(1:length(obs)))
  
  Laake_PrimaryEffort %>% 
    group_by(Start.year) %>%
    reframe(Date = Date,
            Year = first(Start.year),
            n = nwhales,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")) + 1,
            Season = paste0(Year, "/", Year + 1),
            obs = Observer,
            vs = vis,
            bf = beaufort,
            watch.prop = effort) -> Laake.primary.counts
  
  Laake_SecondaryEffort %>% 
    group_by(Start.year) %>%
    reframe(Date = Date,
            Year = first(Start.year),
            n = nwhales,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")) + 1,
            Season = paste0(Year, "/", Year + 1),
            obs = Observer,
            vs = vis,
            bf = beaufort,
            watch.prop = effort) -> Laake.secondary.counts
  
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
    group_by(Start.year) %>%
    summarize(Season = first(Season),
              periods = first(n())) -> Laake.primary.periods
  
  Laake.secondary.counts %>%
    group_by(Start.year) %>%
    summarize(Season = first(Season),
              periods = first(n())) -> Laake.secondary.periods
  
  Laake.primary.periods %>% 
    left_join(Laake.secondary.periods, by = "Season") %>%
    select(Start.year.x, Season, periods.x, periods.y) %>%
    transmute(Start.year = Start.year.x,
              Season = Season,
              periods.1 = periods.x,
              periods.2 = periods.y) -> Laake.periods
  
  day <- bf <- vs <- watch.prop <- obs.input <- n.Laake <- array(dim = c(max(Laake.primary.periods$periods),
                                                                         2, length(all.year)))
  
  y <- 3
  y2 <- 1
  for (y in 1:length(all.year)){
    temp.data <- Laake.primary.counts %>%
      filter(Start.year == all.year[y]) %>%
      arrange(Day) %>%
      left_join(Laake.obs, by = "obs")
    
    n.Laake[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$n
    obs.input[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$obs.ID
    watch.prop[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$watch.prop
    bf[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$bf #scale(temp.data$bf)
    vs[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$vs #scale(temp.data$vs)
    day[1:Laake.primary.periods$periods[y], 1, y] <- as.numeric(temp.data$Day)
    
    # fill in the secondary observations
    if (isTRUE(all.year[y] %in% double.obs.year)){
      temp.data <- Laake.secondary.counts %>%
        filter(Start.year == all.year[y])%>%
        left_join(Laake.obs, by = "obs")
      
      n.Laake[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$n
      obs.input[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$obs.ID
      watch.prop[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$watch.prop
      bf[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$bf #scale(temp.data$bf)
      vs[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$vs #scale(temp.data$vs)
      day[1:Laake.secondary.periods$periods[y2], 2, y] <- as.numeric(temp.data$Day)
      
      y2 <- y2 + 1
      
    }
  }
  
  jags.data <- list(  n = n.Laake,
                      n.station = n.station,
                      n.year = length(all.year),
                      n.obs = nrow(Laake.obs),
                      periods = Laake.periods %>% 
                        select(periods.1, periods.2) %>% simplify2array(),
                      obs = obs.input,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.prop,
                      day = day,
                      n.days = max(day, na.rm = T))
  
  out.list <- list(jags.data = jags.data,
                   all.start.year = all.year,
                   double.obs.year = double.obs.year,
                   obs = Laake.obs,
                   min.dur = min.dur)
  
  return(out.list)  
}


# Create Jags input data for all years using Laake data and WinBUGS data
AllData2JagsInput <- function(min.dur, 
                              WinBUGS.years, WinBUGS.n.stations, WinBUGS.out.file,
                              data.dir){
  library(tidyverse)
  
  load("Data/PrimaryEffort.rda")
  Laake.jags.input <- LaakeData2JagsInput(min.dur = min.dur)
  Laake.jags.data <- Laake.jags.input$jags.data
  Laake.start.years <- unique(PrimaryEffort$Start.year)
  
  # Convert WinBUGS input to Jags input
  Jags.input.2006<- WinBUGSdata2Jags_input(min.dur = min.dur,
                                           years = WinBUGS.years,
                                           WinBUGS.out.file,
                                           n.stations = WinBUGS.n.stations,
                                           data.dir = data.dir)
  
  .data <- Jags.input.2006$jags.data
  # .start.years <- lapply(Jags.input.2006$seasons, 
  #                        FUN = str_split, pattern = "/") %>% 
  #   lapply(FUN = function(x) {tmp <- unlist(x); tmp[1]}) %>%
  #   unlist()
  
  all.start.years <- c(Laake.start.years, Jags.input.2006$start.years)
  all.start.years <- all.start.years[!duplicated(all.start.years)]
  
  # Both datasets contain 2006. Take it out from the recent one.
  bf <- vs <- all.n <- all.obs <- array(dim = c(max(dim(Laake.jags.data$n)[1], 
                                                    dim(.data$n)[1]), 
                                                2, 
                                                (dim(Laake.jags.data$n)[3] + dim(.data$n)[3]-1)))
  
  day <- watch.prop <- array(dim = c((max(dim(Laake.jags.data$n)[1], 
                                         dim(.data$n)[1]) + 2), 
                                     2, 
                                     (dim(Laake.jags.data$n)[3] + dim(.data$n)[3]-1)))
  
  c2 <- 2
  c1 <- c <- 1
  for (k in 1:length(all.start.years)){
    if (all.start.years[k] < 2007){
      n.rows <- dim(Laake.jags.data$n)[1]
      all.n[1:n.rows, 1, c] <- Laake.jags.data$n[, 1, c1]
      all.n[1:n.rows, 2, c] <- Laake.jags.data$n[, 2, c1]
      
      all.obs[1:n.rows, 1, c] <- Laake.jags.data$obs[, 1, c1]
      all.obs[1:n.rows, 2, c] <- Laake.jags.data$obs[, 2, c1]
      
      watch.prop[1:n.rows, 1, c] <- Laake.jags.data$watch.prop[, 1, c1]
      watch.prop[1:n.rows, 2, c] <- Laake.jags.data$watch.prop[, 2, c1]
      
      day[1:n.rows, 1, c] <- Laake.jags.data$day[, 1, c1]
      day[1:n.rows, 2, c] <- Laake.jags.data$day[, 2, c1]
      
      bf[1:n.rows, 1, c] <- Laake.jags.data$bf[, 1, c1]
      bf[1:n.rows, 2, c] <- Laake.jags.data$bf[, 2, c1]
      
      vs[1:n.rows, 1, c] <- Laake.jags.data$vs[, 1, c1]
      vs[1:n.rows, 2, c] <- Laake.jags.data$vs[, 2, c1]
      
      c <- c + 1
      c1 <- c1 + 1
    } else {
      n.rows <- dim(.data$n)[1]
      all.n[1:n.rows, 1, c] <- .data$n[, 1, c2]
      all.n[1:n.rows, 2, c] <- .data$n[, 2, c2]
      
      all.obs[1:n.rows, 1, c] <- .data$obs[, 1, c2]
      all.obs[1:n.rows, 2, c] <- .data$obs[, 2, c2]
      
      watch.prop[1:(n.rows+2), 1, c] <- .data$watch.prop[1:(n.rows+2), 1, c2]
      watch.prop[1:(n.rows+2), 2, c] <- .data$watch.prop[1:(n.rows+2), 2, c2]
      
      day[1:(n.rows+2), 1, c] <- .data$day[, 1, c2]
      day[1:(n.rows+2), 2, c] <- .data$day[, 2, c2]
      
      bf[1:n.rows, 1, c] <- .data$bf[, 1, c2]
      bf[1:n.rows, 2, c] <- .data$bf[, 2, c2]
      
      vs[1:n.rows, 1, c] <- .data$vs[, 1, c2]
      vs[1:n.rows, 2, c] <- .data$vs[, 2, c2]
      c <- c + 1
      c2 <- c2 + 1
    }
    
  }
  
  n.station <- c(Laake.jags.data$n.station, .data$n.station[2:length(.data$n.station)])
  n.year <- length(all.start.years)
  n.obs <- length(unique(as.vector(all.obs))) 
  
  periods <- rbind(Laake.jags.data$periods, .data$periods[2:length(.data$n.station),])
  
  jags.data <- list(n = all.n,
                    n.station = n.station,
                    n.year = n.year,
                    n.obs = n.obs,
                    periods = periods,
                    obs = all.obs,
                    vs = vs,
                    bf = bf,
                    watch.prop = watch.prop,
                    day = day,
                    n.days = max(day, na.rm = T))
  
  out.list <- list(jags.data = jags.data,
                   min.dur = min.dur, 
                   #seasons = seasons, 
                   #WinBUGS.out.file = WinBUGS.out.file,
                   years = years,
                   start.years = all.start.years)
  return(out.list)
}

# Create data input for Jags without running WinBUGS. Only works if
# there are raw data files. 
data2Jags_input_NoBUGS <- function(min.dur,
                                   years,
                                   data.dir){
  
  seasons <- sapply(years, FUN = function(x) paste0(x-1, "/", x))
  start.years <- years - 1
  
  # These are extracted data files (Extract_Data_All_v2.Rmd)
  out.v2 <- lapply(years, 
                   FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
                                                    "_min", min.dur, "_Tomo_v2.rds")))
  
  n.stations <- lapply(out.v2,
                       FUN = function(x){
                         dat <- x$Final_Data %>%
                           mutate(f.station = as.factor(station))
                         return(levels(dat$f.station) %>%
                                  length() %>%
                                  unlist())
                       }) %>%
    unlist()
  
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
  
  Watch.Length <- list()

  for (k in 1:length(begin.)){
    Watch.Length[[k]] <- end.[[k]] - begin.[[k]]
  }
  
  vs <- bf <- array(data = NA, 
                    dim = c(max(periods), 
                            2, length(years)))
  
  periods.mat <- array(data = NA,
                       dim = c(length(periods), 2))
  # observers
  # need to convert observer initials into numbers for all years
  # This needs to be redone... 
  obs.list <- read.csv(file = "Data/ObserverList2023.csv")
  
  obs.2024 <- unique(out.v2[[length(years)]]$Complete_Data$obs)
  new.obs <- obs.2024[!c(obs.2024 %in% obs.list$obs)]
  obs.list <- rbind(obs.list,
                    data.frame(obs = new.obs,
                               ID = seq(max(obs.list$ID) + 1, 
                                        max(obs.list$ID) + length(new.obs))))
  
  watch.length <- n <- day <- array(NA, dim = c(max(periods)+2, 
                                                2, length(years)))
  
  obs <- array(obs.list %>%
                 filter(obs == "No obs") %>%
                 select(ID) %>%
                 pull(), dim = dim(n))
  
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
      n[1:(length(Final_data.P$n)+2),1,k] <- c(Final_data.P$n, 0,0)
      vs[1:length(Final_data.P$n),1,k] <- Final_data.P$vs
      bf[1:length(Final_data.P$n),1,k] <- Final_data.P$bf
      
      obs.year <- data.frame(obs = Final_data.P$obs) %>% 
        left_join(obs.list, by = "obs")
      
      obs[1:length(Final_data.P$n),1,k] <- obs.year$ID
      watch.length[1:(length(Final_data.P$n)+2),1,k] <- c(Final_data.P$watch.length, 0,0)
      day[1:(length(Final_data.P$n)+2),1,k] <- c(floor(Final_data.P$begin), 1, 90)
      periods.mat[k,1] <- nrow(Final_data.P)
      
      # fill in the secondary 
      n[1:(length(Final_data.S$n)+2),2,k] <- c(Final_data.S$n, 0,0)
      vs[1:length(Final_data.S$n),2,k] <- Final_data.S$vs
      bf[1:length(Final_data.S$n),2,k] <- Final_data.S$bf
      
      obs.year <- data.frame(obs = Final_data.S$obs) %>% 
        left_join(obs.list, by = "obs")
      
      obs[1:length(Final_data.S$n),2,k] <- obs.year$ID
      watch.length[1:(length(Final_data.S$n)+2),2,k] <- c(Final_data.S$watch.length, 0,0)
      day[1:(length(Final_data.S$n)+2),2,k] <- c(floor(Final_data.S$begin), 1, 90)
      
      periods.mat[k,2] <- nrow(Final_data.S)
      
    } else {
      n[1:(length(Final_data$n)+2),1,k] <- c(Final_data$n, 0,0)
      vs[1:length(Final_data$n),1,k] <- Final_data$vs
      bf[1:length(Final_data$n),1,k] <- Final_data$bf
      
      obs.year <- data.frame(obs = Final_data$obs) %>% 
        left_join(obs.list, by = "obs")
      
      obs[1:length(Final_data$n),1,k] <- obs.year$ID
      watch.length[1:(length(Final_data$n)+2),1,k] <- c(Final_data$watch.length, 0,0)
      day[1:(length(Final_data$n)+2),1,k] <- c(floor(Final_data$begin), 1, 90)
      periods.mat[k,1] <- nrow(Final_data)
      
    }
    
  }
  
  # Renumber observer IDs because some observers in the list never
  # showed up in the data. The maximum ID number can be greater than
  # the total number of observers, which returns error when running
  # JAGS
  uniq.obs <- apply(obs[,1,], FUN = unique, MARGIN = 2) %>% 
    unlist() %>% 
    unique() %>%
    sort()
  
  obs.list %>%
    filter(ID %in% uniq.obs) %>%
    arrange(ID) %>%
    rownames_to_column(var = "ID.tmp") %>%
    mutate(ID.new = as.numeric(ID.tmp)) %>%
    select(-ID.tmp) -> new.obs.df
  
  # use "replace" to swap the original observer IDs with new IDs
  # using the look up table (new.obs.df)
  # c(new, x)[match(x, c(old, x))], where x is data
  # https://stackoverflow.com/questions/16228160/multiple-replacement-in-r
  obs.new <- array(data = new.obs.df %>% 
                     filter(obs == "No obs") %>% 
                     select(ID.new) %>% pull(),
                   dim = dim(obs))
  
  obs.new[,1,] <- apply(obs[,1,], 
                        FUN = function(x) c(as.vector(new.obs.df$ID.new), 
                                            x)[match(x, as.vector(new.obs.df$ID), x)], 
                        MARGIN = 2)
  
  obs.new[,2,] <- apply(obs[,2,], 
                        FUN = function(x) c(as.vector(new.obs.df$ID.new), 
                                            x)[match(x, as.vector(new.obs.df$ID), x)], 
                        MARGIN = 2)
  
  jags.data <- list(  n = n,
                      n.station = n.stations,
                      n.year = length(years),
                      n.obs = length(uniq.obs),
                      periods = periods.mat,
                      n.days = max(day, na.rm = T),
                      obs = obs.new,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.length,
                      day = day)
  
  out.list <- list(jags.data = jags.data,
                   original.obs = obs,
                   obs = new.obs.df,
                   min.dur = min.dur, 
                   seasons = seasons, 
                   years = years,
                   start.years = start.years,
                   data.dir = data.dir)
  return(out.list)
  
}

# Create jags data input for all years without uring WinBUGS input
AllData2JagsInput_NoBUGS <- function(min.dur, years, data.dir){
  
  # this converts Laake's data into jags input
  jags.input.Laake <- LaakeData2JagsInput(min.dur)
  
  # This pulls out data from 2010 onward and create jags input
  jags.input.new <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                           years = years,
                                           data.dir = data.dir)
  
  # Combine all observers from both datasets and give them new ID numbers
  obs.Laake <- jags.input.Laake$obs %>%
    mutate(data = "Laake") %>%
    rename(ID = obs.ID)
  
  # In Laake's dataset, there is no ID for "No observers" and NAs are used
  # to fill in when there were no observers. 
  obs.Laake <- rbind(obs.Laake, 
                     data.frame(obs = "No obs", 
                                ID = nrow(obs.Laake)+1, 
                                data = "Laake"))
  
  # Note the no observer ID has been changed from 36 to whatever the number
  # that was assigned
  obs.new <- jags.input.new$obs %>%
    mutate(data = "new")
  
  obs.Laake %>% 
    full_join(obs.new, by = "obs") %>%
    rename(ID.Laake = ID.x) -> obs.all
  
  obs.all %>%
    mutate(ID.all = seq(1:nrow(obs.all))) %>%
    select(obs, ID.Laake, ID.new, ID.all) -> obs.all.new
  
  # Replace observer IDs in each dataset with new IDs (obs.all.new$ID.all)
  jags.obs.Laake <- jags.input.Laake$jags.data$obs
  
  # Replace NAs with no observer ID
  jags.obs.Laake[is.na(jags.obs.Laake)] <- obs.Laake[nrow(obs.Laake), "ID"]
  
  # use "match" to swap the original observer IDs with new IDs
  # using the look up table (obs.all.new)
  # c(new, x)[match(x, c(old, x))], where x is data
  # https://stackoverflow.com/questions/16228160/multiple-replacement-in-r
  
  No.obs.ID <- obs.all.new %>% 
    filter(obs == "No obs") %>% 
    select(ID.all) %>% pull()
  
  jags.obs.Laake.new <- array(data = No.obs.ID,
                              dim = dim(jags.obs.Laake))
  
  jags.obs.Laake.new[,1,] <- apply(jags.obs.Laake[,1,], 
                                   FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                       x)[match(x, as.vector(obs.all.new$ID.Laake), x)], 
                                   MARGIN = 2)
  
  jags.obs.Laake.new[,2,] <- apply(jags.obs.Laake[,2,], 
                                   FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                       x)[match(x, as.vector(obs.all.new$ID.Laake), x)], 
                                   MARGIN = 2)
  
  jags.obs.new <- jags.input.new$jags.data$obs
  
  jags.obs.new.new <- array(data = No.obs.ID,
                            dim = dim(jags.obs.new))
  
  jags.obs.new.new[,1,] <- apply(jags.obs.new[,1,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.new), x)], 
                                 MARGIN = 2)
  
  jags.obs.new.new[,2,] <- apply(jags.obs.new[,2,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.new), x)], 
                                 MARGIN = 2)
  
  # all data arrays need to be combined between Laake and new datasets
  # First make the first dimension the same size between the two arrays:
  n.row.dif <- dim(jags.obs.Laake.new)[1] - dim(jags.obs.new.new)[1]
  
  if (n.row.dif > 0){
    # observers:
    jags.obs.new.new.1 <- abind(jags.obs.new.new,
                                array(data = No.obs.ID, 
                                      dim = c(n.row.dif, 2, 
                                              dim(jags.obs.new.new)[3])),
                                along = 1)
    
    # Combine the arrays.
    jags.obs.all <- abind(jags.obs.Laake.new, 
                          jags.obs.new.new.1,
                          along = 3)
    
    # whale counts:
    jags.n.new <- abind(jags.input.new$jags.data$n,
                        array(data = NA, 
                              dim = c(n.row.dif, 2, 
                                      dim(jags.obs.new.new)[3])),
                        along = 1)
    
    jags.n.all <- abind(jags.input.Laake$jags.data$n,
                        jags.n.new,
                        along = 3)
    
    # Watch length
    jags.watch.prop.new <- abind(jags.input.new$jags.data$watch.prop,
                                 array(data = NA, 
                                       dim = c(n.row.dif, 2, 
                                               dim(jags.obs.new.new)[3])),
                                 along = 1)
    
    jags.watch.prop.all <- abind(jags.input.Laake$jags.data$watch.prop,
                                 jags.watch.prop.new,
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
                                array(data = No.obs.ID, 
                                      dim = c(n.row.dif, 2, 
                                              dim(jags.obs.Laake.new)[3])),
                                along = 1)
    
    jags.obs.all <- abind(jags.obs.Laake.new.1, 
                          jags.obs.new.new,
                          along = 3)
    
    # whale counts:
    jags.n.Laake <- abind(jags.input.Laake$jags.data$n,
                        array(data = NA, 
                              dim = c(n.row.dif, 2, 
                                      dim(jags.obs.Laake.new)[3])),
                        along = 1)
    
    jags.n.all <- abind(jags.n.Laake,
                        jags.input.new$jags.data$n,
                        along = 3)
    
    # Watch length
    jags.watch.prop.Laake <- abind(jags.input.Laake$jags.data$watch.prop,
                                 array(data = NA, 
                                       dim = c(n.row.dif, 2, 
                                               dim(jags.obs.Laake.new)[3])),
                                 along = 1)
    
    jags.watch.prop.all <- abind(jags.watch.prop.Laake,
                                 jags.input.new$jags.data$watch.prop,
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
  
  # visibility
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
  }

  jags.data <- list(  n = jags.n.all,
                      n.station = c(jags.input.Laake$jags.data$n.station,
                                    jags.input.new$jags.data$n.station),
                      n.year = dim(jags.bf.all)[3],
                      n.obs = nrow(obs.all.new),
                      periods = rbind(as.matrix(jags.input.Laake$jags.data$periods), 
                                      jags.input.new$jags.data$periods),
                      n.days = max(jags.input.Laake$jags.data$n.days,
                                   jags.input.new$jags.data$n.days),
                      obs = jags.obs.all,
                      vs = jags.vs.all,
                      bf = jags.bf.all,
                      watch.prop = jags.watch.prop.all,
                      day = jags.day.all)
  
  return(list(jags.data = jags.data,
              jags.input.Laake = jags.input.Laake,
              jags.input.new = jags.input.new))
  
}

# Create data input for Jags from WinBUGS input.
WinBUGSdata2Jags_input <- function(min.dur, 
                                   WinBUGS.out.file,
                                   years,
                                   n.stations,
                                   data.dir){
  
  all.years <- c(2007, 2008, years)
  seasons <- sapply(all.years, FUN = function(x) paste0(x-1, "/", x))
  
  start.years <- all.years - 1
  
  WinBUGS.out <- readRDS(WinBUGS.out.file)
  
  WinBUGS.inputs <- data2WinBUGS_input(data.dir = data.dir,
                                       years = years,
                                       min.dur = min.dur)
  
  data.WinBUGS <- WinBUGS.inputs$data
  # watch lengths are assumed equal between primary and secondary stations in
  # WinBUGS code. But not in Jags. So, I duplicate the secondary watch effort
  
  bf <- vs <- array(dim = c(max(data.WinBUGS$periods),
                            2, 
                            length(data.WinBUGS$periods)))
  
  watch.prop <- day <- array(dim = c(dim(data.WinBUGS$Watch.Length)[1],
                                     2, 
                                     dim(data.WinBUGS$Watch.Length)[2]))
  
  bf[,1,] <- data.WinBUGS$bf
  bf[,2,] <- data.WinBUGS$bf
  
  vs[,1,] <- data.WinBUGS$vs
  vs[,2,] <- data.WinBUGS$vs
  
  day[,1,] <- data.WinBUGS$day
  day[,2,] <- data.WinBUGS$day
  
  watch.prop[,1,] <- data.WinBUGS$Watch.Length
  watch.prop[,2,] <- data.WinBUGS$Watch.Length
  
  jags.data <- list(  n = data.WinBUGS$n,
                      n.station = n.stations,
                      n.year = length(seasons),
                      n.obs = data.WinBUGS$n.obs,
                      #Daily.N = daily.N,
                      periods = cbind(data.WinBUGS$periods,
                                      data.WinBUGS$periods),
                      n.days = max(day, na.rm = T),
                      #first.day = unlist(as.vector(first.day)),
                      obs = data.WinBUGS$obs,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.prop,
                      #watch.prop = (data.WinBUGS$Watch.Length*24*60)/540,
                      day = day)
  
  out.list <- list(jags.data = jags.data,
                   min.dur = min.dur, 
                   seasons = seasons, 
                   WinBUGS.out.file = WinBUGS.out.file,
                   years = years,
                   start.years = start.years,
                   data.dir = data.dir)
  return(out.list)
}


# Create trace and density plots for MCMC samples from jagsUI
plot.trace.dens <- function(param, jags.out){
  if (!is.character(jags.out)) stop("jags.out input has to be a character string")
  
  n.param <- ncol(eval(parse(text = paste0(jags.out, "$sims.list$", param))))
  samples <- eval(parse(text = paste0(jags.out, "$samples")))
  
  if (!is.null(n.param)){
    
    par.idx <- c(1:n.param)
    p.trace <- bayesplot::mcmc_trace(samples, paste0(param, "[", par.idx, "]"))
    p.dens <- bayesplot::mcmc_dens(samples, paste0(param, "[", par.idx, "]"))
    
  } else {
    p.trace <- bayesplot::mcmc_trace(samples, param)
    p.dens <- bayesplot::mcmc_dens(samples, param)
    
  }
  
  return(list(trace = p.trace,
              dens = p.dens))
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


compute.LOOIC <- function(loglik.array, data.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  # Convert the log-likelihood and data arrays into vectors
  loglik.vec <- c(loglik.array)
  data.vec <- c(data.array)
  
  # Convert the log-likelihood vector into a 2D array
  n.data <- length(data.vec)
  loglik.mat <- array(loglik.vec, 
                      c((n.per.chain * MCMC.params$n.chains), n.data))
  
  # remove the log-likelihood values that correspond to NA data points
  loglik.mat <- loglik.mat[, !is.na(data.vec)]
  
  # Also, some of data points (0s) were unobserved and no log likelihood values
  # exist
  colsums.loglik <- colSums(loglik.mat)
  loglik.mat <- loglik.mat[, !is.na(colsums.loglik)]
  
  # remove NAs in the data vector
  #data.vec <- data.vec[!is.na(data.vec)]
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = 4, k_threshold = 0.7)
  
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
    
    if (YEAR != 2010){
      tmp <- hour(hms(data[nrow(data),4])) + 
        minute(hms(data[nrow(data),4]))/60 + 
        (second(hms(data[nrow(data),4])) + 1)/3600 
      
    } else {
      tmp <- data[nrow(data), 4]
    }
    
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
  if (YEAR != 2010){
    BeginHr <- hour(hms(data$V4)) + minute(hms(data$V4))/60 + second(hms(data$V4))/3600    
  } else {
    BeginHr <- as.numeric(data$V4)
  }

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
  
  if (YEAR == 2010){
    data$Shift <- shift.definition.2010(data$V4)
  } else {
    data$Shift <- shift.definition(as.Date(data$V3, format = "%m/%d/%Y"), data$V4)
  }

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
  BeginHr <- ifelse(YEAR != 2010,
                    (hour(hms(data[shifts.begin[i], 4])) + 
                       (minute(hms(data[shifts.begin[i], 4]))/60) 
                     + (second(hms(data[shifts.begin[i], 4]))/3600)),
                    as.numeric(data[shifts.begin[i], 4]))
    
  # Decimal hour of next shift start time
  if (i < max.shifts){
    event.idx <- which(shifts.df$shift %in% i)
    next.event <- shifts.df[event.idx + 1,]
    if (next.event$event == "P"){
      NextBeginHr <- ifelse(YEAR != 2010,
                            (hour(hms(data[next.event$row, 4])) + 
                               (minute(hms(data[next.event$row, 4]))/60)
                             + (second(hms(data[next.event$row, 4]))/3600)),
                            as.numeric(data[next.event$row, 4]))
      
      EndHr <- NextBeginHr - 0.00001
    } else {  # if the event is "E"
      next.P <- shifts.df[event.idx+2,]
      NextBeginHr <- ifelse(YEAR != 2010, 
                            (hour(hms(data[next.P$row, 4])) + 
                               (minute(hms(data[next.P$row, 4]))/60) 
                             + (second(hms(data[next.P$row, 4]))/3600)),
                            as.numeric(data[next.P$row, 4]))
      
      EndHr <- ifelse(YEAR != 2010,
                      (hour(hms(data[next.event$row, 4])) + 
                         (minute(hms(data[next.event$row, 4]))/60) + 
                         (second(hms(data[next.event$row, 4]))/3600)),
                      as.numeric(data[next.event$row, 4]))
    }

    # Find next end hr to find the next shift to figure out spillovers
    event.idx2 <- which(shifts.df$shift %in% (i+1))
    next.event2 <- shifts.df[event.idx2+1,]
    
    if (next.event2$event == "P"){   # if the next event is also "P"
       NextEndHr <- ifelse(YEAR != 2010,
                           (hour(hms(data[next.event2$row, 4])) + 
                              (minute(hms(data[next.event2$row, 4]))/60) + 
                              (second(hms(data[next.event2$row, 4]))/3600) ) - 0.00001,
                           as.numeric(data[next.event2$row, 4]) - 0.00001)
       # 
    } else {  # if the event is "E"
       NextEndHr <- ifelse(YEAR != 2010,
                           (hour(hms(data[next.event2$row, 4])) + 
                              (minute(hms(data[next.event2$row, 4]))/60) 
                            + (second(hms(data[next.event2$row, 4]))/3600)),
                           as.numeric(data[next.event2$row, 4]))
    }
    

  } else {    # for the last shift
    event.idx <- which(shifts.df$shift %in% max.shifts)
    next.event <- shifts.df[event.idx + 1,]  # This has to be E
    if (length(next.event) == 0){
      end.row <- nrow(data)
    } else {
      end.row <- next.event$row
    }
    EndHr <-  ifelse(YEAR != 2010,
                     (hour(hms(data[end.row, 4])) + 
                        (minute(hms(data[end.row, 4]))/60) + 
                        (second(hms(data[end.row, 4]))/3600)),
                     as.numeric(data[end.row, 4]))
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
    if (YEAR == 2010){
      # This Shift is the defined shift IDs - not based on changes in observers
      Shift.End <- shift.definition.2010(data.shift$V4[nrow(data.shift)])
      if (NextBeginHr > (BeginHr + 1.5)){
        V4 <- shift.hrs %>% 
          filter(shift == Shift.End) %>% 
          select(end.hr) %>%
          pull()        
      } else {
        V4 <- NextBeginHr - 0.00001
      }

    } else {
      Shift.End <- shift.definition(as.Date(data.shift$V3[nrow(data.shift)], 
                                            format = "%m/%d/%Y"), 
                                    data.shift$V4[nrow(data.shift)])
      if (NextBeginHr > (BeginHr + 1.5)){
        V4 <- format(as.POSIXct(as.Date("2022-12-01 00:00:00") + End),
                   format = "%H:%M:%S")
      } else {
        V4 <- NextBeginHr - 0.00001
      }
    }
    
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
      summarise(time = first(begin)) %>%
      pull(time) %>% as.numeric()
    
    if (k < max(key.num)){
      tmp.2 <- data.shift %>%
        filter(key == (k + 1)) %>%
        summarise(time = first(begin)) %>%
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
                Shift = first(Shift)) -> data.shift.effort
    
    
    
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


