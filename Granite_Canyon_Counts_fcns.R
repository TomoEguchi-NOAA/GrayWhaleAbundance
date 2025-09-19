
# define some functions

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
                                 data.set = "2007to2024")
  return(WinBUGS.Nhats.df) 
  
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
                          Final_Data = all.Final_Data))  
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


