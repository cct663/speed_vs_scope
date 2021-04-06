### Script for speed vs. scope manuscript
# Written by Conor Taff; cct63@cornell.edu; cct663@gmail.com


## Load packages ----
    pacman::p_load(ggplot2, here, MASS, MBESS, tidyverse, viridis, gridExtra, sp, rptR)

## Create simulation function ----
  # The first function simulates cort response parameters for a population of animals of size n. We
    # can think of these as the 'true' phenotypes of these individuals. For purposes of this simulation
    # the cort response is described with five simulated parameters that can covary or not.
    # The parameters estimated are: 
          #base cort
          #time to reach max cort
          #max cort
          #time at max cort level
          #time of return to base cort level 
          cort_sim1 <- function(n = 20, 
                               base_mu = 5, base_sd = 2, base_min = 0.5,
                               slope_mu = 15, slope_sd = 2,
                               fastpct_mu = 0.55, fastpct_sd = 0.08, 
                               speed_mu = 25, speed_sd = 5, speed_min = 10,
                               max_mu = log(45), max_sd = 0.21, 
                               maxtime_mu = 10, maxtime_sd = 3, maxtime_min = 5,
                               return_mu = 90, return_sd = 15, 
                               extra_n = 500,
                               cor_base_speed = 0, cor_base_max = 0, cor_base_maxtime = 0, cor_base_return = 0, cor_base_slope = 0, cor_base_fastpct = 0,
                               cor_speed_max = 0, cor_speed_maxtime = 0, cor_speed_return = 0, cor_speed_slope = 0, cor_speed_fastpct = 0,
                               cor_max_maxtime = 0, cor_max_return = 0, cor_max_slope = 0, cor_max_fastpct = 0,
                               cor_maxtime_return = 0, cor_maxtime_slope = 0, cor_maxtime_fastpct = 0,
                               cor_return_slope = 0, cor_return_fastpct = 0,
                               cor_slope_fastpct = 0)
              {
                  # Build a full variance-covariance matrix based on the values for each parameter and correlation between parameters
                      cort_vcov <- matrix(nrow = 7, ncol = 7)
                      cort_vcov[1, ] <- c(base_sd^2, cor_base_speed * base_sd * speed_sd, cor_base_max * base_sd * max_sd, 
                                          cor_base_maxtime * base_sd * maxtime_sd, cor_base_return * base_sd * return_sd,
                                          cor_base_slope * base_sd * slope_sd, cor_base_fastpct * base_sd * fastpct_sd)
                      cort_vcov[2, ] <- c(cor_base_speed * base_sd * speed_sd, speed_sd^2, cor_speed_maxtime * speed_sd * max_sd, 
                                          cor_speed_maxtime * speed_sd * maxtime_sd, cor_speed_return * speed_sd * return_sd,
                                          cor_speed_slope * speed_sd * slope_sd, cor_speed_fastpct * speed_sd * fastpct_sd)
                      cort_vcov[3, ] <- c(cor_base_max * base_sd * max_sd, cor_speed_max * speed_sd * max_sd, max_sd^2, 
                                          cor_max_maxtime * max_sd * maxtime_sd, cor_max_return * max_sd * return_sd,
                                          cor_max_slope * max_sd * slope_sd, cor_max_fastpct * max_sd * fastpct_sd)
                      cort_vcov[4, ] <- c(cor_base_maxtime * base_sd * maxtime_sd, cor_speed_maxtime * speed_sd * maxtime_sd, 
                                          cor_max_maxtime * max_sd * maxtime_sd, maxtime_sd^2, cor_maxtime_return * maxtime_sd * return_sd,
                                          cor_maxtime_slope * maxtime_sd * slope_sd, cor_maxtime_fastpct * maxtime_sd * fastpct_sd)
                      cort_vcov[5, ] <- c(cor_base_return * base_sd * return_sd, cor_speed_return * speed_sd * return_sd, 
                                          cor_max_return * max_sd * return_sd, cor_maxtime_return * maxtime_sd * return_sd, return_sd^2,
                                          cor_return_slope * return_sd * slope_sd, cor_return_fastpct * return_sd * fastpct_sd)
                      cort_vcov[6, ] <- c(cor_base_slope * base_sd * slope_sd, cor_speed_slope * speed_sd * slope_sd, 
                                          cor_max_slope * max_sd * slope_sd, cor_maxtime_slope * maxtime_sd * slope_sd,
                                          cor_return_slope * return_sd * slope_sd, slope_sd^2, cor_slope_fastpct * slope_sd * fastpct_sd)
                      cort_vcov[7, ] <- c(cor_base_fastpct * base_sd * fastpct_sd, cor_speed_fastpct * speed_sd * fastpct_sd, 
                                          cor_max_fastpct * max_sd * fastpct_sd, cor_maxtime_fastpct * maxtime_sd * fastpct_sd,
                                          cor_return_fastpct * return_sd * fastpct_sd, cor_slope_fastpct * slope_sd * fastpct_sd, fastpct_sd^2)
                 
                  # Make a vector of the mean values of each parameter     
                      mu_list <- c(base_mu, speed_mu, max_mu, maxtime_mu, return_mu, slope_mu, fastpct_mu)
                  
                  # Create a dataset by sampling from a multivariate normal distribution
                        msam <- mvrnorm(n + extra_n, mu = mu_list, Sigma = cort_vcov)
                      
                      # Make a dataframe from the matrix above with columns named by paramter  
                        sim_dat <- data.frame(
                          animal = paste("id", 1:(n + extra_n), sep = "_"),
                          base = msam[, 1],
                          tmax = msam[, 2],
                          max = exp(msam[, 3]) + msam[, 1],
                          atmax = msam[, 4],
                          return = msam[, 5],
                          slope = msam[, 6],
                          fastpct = msam[, 7]
                        )
                      
                # calculate return time as return latency plus speed plus time at max
                      sim_dat$return <- sim_dat$tmax + sim_dat$atmax + sim_dat$return
                          
                # replace base, maxtime, and speed values with minimum allowed (avoids negative numbers or values that don't make sense)
                    sim_dat[which(sim_dat$base < base_min), "base"] <- base_min
                    sim_dat[which(sim_dat$tmax < speed_min), "tmax"] <- speed_min
                    sim_dat[which(sim_dat$atmax < maxtime_min), "atmax"] <- maxtime_min
                    for(j in 1:nrow(sim_dat)){
                      if(sim_dat$slope[j] > 0.8 * sim_dat$tmax[j]){
                        sim_dat$slope[j] <- 0.8 * sim_dat$tmax[j]
                      }
                    }
                  
                  # Convert simulated data to long form x/y dataframe
                        sim_dat2 <- data.frame(
                          animal = rep(sim_dat$animal[1:n], 5),
                          x = c(rep(0, n),
                                sim_dat$slope[1:n],
                                sim_dat$tmax[1:n],
                                sim_dat$atmax[1:n] + sim_dat$tmax[1:n],
                                sim_dat$return[1:n]),
                          y = c(sim_dat$base[1:n],
                                (sim_dat$max[1:n] - sim_dat$base[1:n]) * sim_dat$fastpct[1:n] + sim_dat$base[1:n],
                                sim_dat$max[1:n],
                                sim_dat$max[1:n],
                                sim_dat$base[1:n]),
                          group = c(rep("initial", n),
                                    rep("slopeturn", n),
                                    rep("reachmax", n),
                                    rep("startdecline", n),
                                    rep("return", n))
                        )
                        
                  # Make a data set of 'extra' response values that can be sampled from to add noise in second simulation
                        sim_dat3 <- data.frame(
                          animal = rep(sim_dat$animal[(n + 1):(n + extra_n)], 5),
                          x = c(rep(0, extra_n),
                                sim_dat$slope[(n + 1):(n + extra_n)],
                                sim_dat$tmax[(n + 1):(n + extra_n)],
                                sim_dat$atmax[(n + 1):(n + extra_n)] + sim_dat$tmax[(n + 1):(n + extra_n)],
                                sim_dat$return[(n + 1):(n + extra_n)]),
                          y = c(sim_dat$base[(n + 1):(n + extra_n)],
                                (sim_dat$max[(n + 1):(n + extra_n)] - sim_dat$base[(n + 1):(n + extra_n)]) * sim_dat$fastpct[(n + 1):(n + extra_n)] + sim_dat$base[(n + 1):(n + extra_n)],
                                sim_dat$max[(n + 1):(n + extra_n)],
                                sim_dat$max[(n + 1):(n + extra_n)],
                                sim_dat$base[(n + 1):(n + extra_n)]),
                          group = c(rep("initial", extra_n),
                                    rep("slopeturn", extra_n),
                                    rep("reachmax", extra_n),
                                    rep("startdecline", extra_n),
                                    rep("return", extra_n))
                        )
                        
                        dataset <- list(sim_dat2 = sim_dat2,
                                        cort_vcov = cort_vcov,
                                        extra_sims = sim_dat3)
                        return(dataset)
      }
      
    # The second function takes input from the first function (or by default, generates a dataset with the first function).
        # These parameters are used to fit and estimate full time courses of cort responses for each animal. Noise is added
        # to each of the 'true' parameters and the amount of noise can be specified. This should be considered within individual
        # variation in expression of underlying 'true' values. A second parameter sets assay error as a percentage. This is like
        # the CV of your lab assay + any other technical issues that add to measurement error. The measurement error only 
        # effects the downsampled dataset with 'bleed' data. The same animal can be sampled
        # multiple times, with noise added separately each time. The output is a list with three dataframes that represent
        # a simulated dataset with samples at specified time points, the full time course with cort at every time point, and a 
        # rank converted estimate. 
          cort_sim2 <- function(data = cort_sim1(),
                                bleed_times = c(1, 15, 30),  # Time in minutes to save 'bleed' samples at
                                base_error = 0.8,            # Variation in base from 'true' value as a percentage of true value
                                speed_error = 0.6,           # Variation in time to reach max from 'true value as percentage
                                max_error = 0.5,             # Variation in max cort from 'true' value as percentage
                                maxtime_error = 0.5,         # Variation in time at max cort from 'true' value as percentage
                                return_error = 0.5,          # Variation in time to return to base from 'true' value as percentage
                                slope_error = 0.5,           # Variation in length of fast slope initial increase from 'true' value as percentage
                                fastpct_error = 0.5,         # Variation in percent of max reached during initial fast increase from 'true' as percentage
                                sample_times = 2,            # Number of samples to draw for each animal
                                assay_error = 0.2,          # Measurement error after sampling as percent of value normal error mean = 0
                                timecourse_max = 170,        # Number of minutes for the full time course
                                performance_contributions = c(0, 0, 0, 0, 0, 1, 0, 1), 
                                        # Relative contributions of base, speed, max, maxtime, return, slope, fastpct, and random error to performance
                                sm_span = 0.5              # Smoothing parameter for loess regression. Smaller = more wiggly.
                                ) 
              {
            
            # Here is some ugly data wrangling that is just getting everything into the right format to proceed
                # Increases length of the 'true' data to account for multiple samples taken from each individual animal
                    data2 <- data$sim_dat2[rep(seq_len(nrow(data$sim_dat2)), sample_times), ]
                # Adds a column to identify sample number when animals are sampled more than once
                    data2$sample <- rep(seq(1, sample_times, 1), each = nrow(data$sim_dat2))
                # Makes a new identifier that includes both animal id and sample number
                    data2$animal_sample <- paste(data2$animal, data2$sample, sep = "_")
                # Pivots to a wide format with all paramters for each unique animal sample combination in one row
                    data3 <- as.data.frame(pivot_wider(data2, id_cols = animal_sample, names_from = group, values_from = c(x, y)))
                    
            # Calculating the performance metric
                # Sums and standardizes performance contributions to make them relative
                    per_conts <- performance_contributions / sum(performance_contributions)
                # Scales each parameter to the same z scored scale
                    data3$maxtime <- data3$x_startdecline - data3$x_reachmax
                    data3$maxtime <- scale(data3$maxtime)
                    data3$y_initial <- scale(data3$y_initial)
                    data3$x_reachmax <- scale(data3$x_reachmax)
                    data3$y_reachmax <- scale(data3$y_reachmax)
                    data3$x_return <- scale(data3$x_return)
                    data3$x_slopeturn <- scale(data3$x_slopeturn)
                    data3$fastpct <- (data3$y_slopeturn - data3$y_initial) / (data3$y_reachmax - data3$y_initial)
                    data3$fastpct <- scale(data3$fastpct)
                    for(k in 1:nrow(data3)){
                      data3$performance[k] <- sum(per_conts[1:7] * 
                            data3[k, c("y_initial", "x_reachmax", "y_reachmax", "maxtime", "x_return", "x_slopeturn", "fastpct")]) +
                            (rnorm(1, mean = 0, sd = 1) * per_conts[8])
                    }
                # Join performance measure back to the initial data frame
                    data2 <- plyr::join(data2, data3[, c("animal_sample", "performance")], "animal_sample", "left", "first")
                    
                    
            
            # Adding noise to 'true' values to represent within individual variation in response expression
                # This is best thought of as within individual variation. It is set by changing the values to a percent of the observed value
                  # that is attributable to a random draw from the population distribution (error) and a percent that is attributable to the 
                  # true underlying value of the individual (1 - error).
                    
                  # Sample an 'extra' cort response curve to add noise from. This maintains the population level covariance structure and
                    # means so that the parameters all vary together.
                      id <- sample(unique(data$extra_sims$animal), sample_times * (nrow(data$sim_dat2) / 5))
                      noise <- subset(data$extra_sims, data$extra_sims$animal %in% id)
                    
                    
                  # Noise for base cort
                    data2[which(data2$group == "initial"), "y"] <- (1 - base_error) * data2[which(data2$group == "initial"), "y"] +
                      base_error * noise[which(noise$group == "initial"), "y"]
                      
                    # Old version that didn't maintain correlation structure. Has been replaced.  
                      #rnorm(n = length(data2[which(data2$group == "initial"), "y"]),
                       #                  mean = mean(data2[which(data2$group == "initial"), "y"]),
                         #                sd = sd(data2[which(data2$group == "initial"), "y"]))
                    
                  # Noise for speed to response
                    data2[which(data2$group == "slopeturn"), "x"] <- (1 - slope_error) * data2[which(data2$group == "slopeturn"), "x"] +
                      slope_error * noise[which(noise$group == "slopeturn"), "x"]
                      
                        # Old version that didn't maintain correlation structure. Has been replaced.
                            #rnorm(n = length(data2[which(data2$group == "slopeturn"), "x"]),
                            #    mean = mean(data2[which(data2$group == "slopeturn"), "x"]), 
                            #    sd = sd(data2[which(data2$group == "slopeturn"), "x"]))
                    
                  # noise for fast percent
                    data2[which(data2$group == "slopeturn"), "y"] <- (1 - fastpct_error) * data2[which(data2$group == "slopeturn"), "y"] +
                      fastpct_error * noise[which(noise$group == "slopeturn"), "y"]
                      
                            # Old version that didn't maintain correlation structure. Has been replaced.
                                      #rnorm(n = length(data2[which(data2$group == "slopeturn"), "y"]),
                                      #      mean = mean(data2[which(data2$group == "slopeturn"), "y"]),
                                      #      sd = sd(data2[which(data2$group == "slopeturn"), "y"]))

                  # Noise for speed to reach maxcort
                    data2[which(data2$group == "reachmax"), "x"] <- (1 - speed_error) * data2[which(data2$group == "reachmax"), "x"] +
                      speed_error * noise[which(noise$group == "reachmax"), "x"]
                      
                                # Old version that didn't maintain correlation structure. Has been replaced.
                                      #rnorm( n = length(data2[which(data2$group == "reachmax"), "x"]),
                                       #    mean = mean(data2[which(data2$group == "reachmax"), "x"]),
                                        #   sd = sd(data2[which(data2$group == "reachmax"), "x"]))
                    
                  # Noise for maxcort
                    data2[which(data2$group == "reachmax"), "y"] <- (1 - max_error) * data2[which(data2$group == "reachmax"), "y"] +
                      max_error * noise[which(noise$group == "reachmax"), "y"]
                                
                                # Old version that didn't maintain correlation structure. Has been replaced.
                                    #rnorm(n = length(data2[which(data2$group == "reachmax"), "y"]),
                                    #    mean = mean(data2[which(data2$group == "reachmax"), "y"]),
                                    #    sd = sd(data2[which(data2$group == "reachmax"), "y"]))
                        
                  # change start decline to match new max cort
                    data2[which(data2$group == "startdecline"), "y"] <- data2[which(data2$group == "reachmax"), "y"]
                    
                  # Noise for max time
                    data2[which(data2$group == "startdecline"), "x"] <- (1 - maxtime_error) * data2[which(data2$group == "startdecline"), "x"] +
                      maxtime_error * noise[which(noise$group == "startdecline"), "x"]
                      
                              # Old version that didn't maintain correlation structure. Has been replaced.
                                  #rnorm(n = length(data2[which(data2$group == "startdecline"), "x"]),
                                   #         mean = mean(data2[which(data2$group == "startdecline"), "x"]),
                                    #        sd = sd(data2[which(data2$group == "startdecline"), "x"]))
                      
                  # Noise for return
                    data2[which(data2$group == "return"), "x"] <- (1 - return_error) * data2[which(data2$group == "return"), "x"] +
                      return_error * noise[which(noise$group == "return"), "x"]
                      
                                  # Old version that didn't maintain correlation structure. Has been replaced.
                                        #rnorm(n = length(data2[which(data2$group == "return"), "x"]),
                                         #  mean = mean(data2[which(data2$group == "return"), "x"]),
                                          # sd = sd(data2[which(data2$group == "return"), "x"]))
                    
                  # noise for return basecort level
                    data2[which(data2$group == "return"), "y"] <- (1 - base_error) * data2[which(data2$group == "return"), "y"] +
                      base_error * noise[which(noise$group == "return"), "y"]
                      
                                # Old version that didn't maintain correlation structure. Has been replaced.
                                      #rnorm(n = length(data2[which(data2$group == "return"), "y"]),
                                       #    mean = mean(data2[which(data2$group == "return"), "y"]),
                                        #   sd = sd(data2[which(data2$group == "return"), "y"]))  
                
            # Make a time sequence with filler points every 5 minutes to make loess fit smoother. I did this because otherwise you end
                # up with different length (of time) gaps on the x axis and the curves can fit in weird ways with certain parameters.
                    
                  # Set up an object with slots every minute to record the full time series data
                        time_seq <- data.frame(time = seq(0, timecourse_max, 1))
                        preds <- as.data.frame(matrix(nrow = length(unique(data2$animal_sample)), ncol = nrow(time_seq)))
                        colnames(preds) <- time_seq$time
                        rownames(preds) <- unique(data2$animal_sample)
                        
                  # Set up an object with slots every 5 minutes that will be used for fitting the loess curve. This could be
                      # tried with different time intervals and might give some flexibility on how wiggly the line is. Would
                      # be easy to add that as an option in the function if useful.
                          tseq <- data.frame(time = seq(0, timecourse_max, 1), y = NA)
                          preds2 <- as.data.frame(matrix(nrow = length(unique(data2$animal_sample)), ncol = nrow(tseq)))
                          colnames(preds2) <- tseq$time
                          rownames(preds2) <- unique(data2$animal_sample)
                          
                  # Loop through each animal and fill in based on straight lines between the sampled points where their cort value
                      # would be every five minutes. Then fit a loess curve of those straight lines with adjustible smoothing and
                      # save the value every one minute for the output.
                            for(i in 1:length(unique(data2$animal_sample))){
                                # make a temporary dataset with just one animal sample combination
                                    temp <- subset(data2, data2$animal_sample == unique(data2$animal_sample)[i])
                                # The first time (minute 0) = the initial baseline parameter
                                    tseq$y[1] <- temp$y[1]
                                # For times between 0 and slopeturn, samples are taken from the line drawn between those points
                                    for(k in 2:nrow(tseq)){
                                        if(tseq$time[k] < temp$x[2]){
                                          slope <- (temp$y[2] - temp$y[1]) / (temp$x[2] - temp$x[1])
                                          tseq$y[k] <- tseq$time[k] * slope + tseq$y[1]
                                        }
                                    }
                                # For times between slopeturn and reaching max a different line is used
                                    for(k in 2:nrow(tseq)){
                                        if(tseq$time[k] < temp$x[3]){
                                          if(tseq$time[k] > temp$x[2]){
                                            slope <- (mean(temp$y[3:4]) - temp$y[2]) / (temp$x[3] - temp$x[2])
                                            tseq$y[k] <- (tseq$time[k] - temp$x[2]) * slope + temp$y[2]
                                          }
                                        }
                                    }
                                # For times at the plateau, the max value is used.    
                                    for(k in 2:nrow(tseq)){
                                      if(tseq$time[k] > temp$x[3]){
                                        if(tseq$time[k] < temp$x[4]){
                                          tseq$y[k] <- mean(temp$y[3:4])
                                        }
                                      }
                                    }
                               # For times in the decline, a straight slope is used. This one needs a counter to account for the time correctly.     
                                    counter <- 0
                                    for(k in 2:nrow(tseq)){
                                      if(tseq$time[k] > temp$x[4]){
                                        if(tseq$time[k] < temp$x[5]){
                                          slope <- (temp$y[5] - mean(temp$y[3:4])) / (temp$x[5] - temp$x[4])
                                          tseq$y[k] <- mean(temp$y[3:4]) + (tseq$time[k] - tseq$time[k - 1] + counter) * slope
                                          counter <- counter + 1
                                        }
                                      }
                                    }
                                # Finally, after reaching decline the rest are filled in with the return value.
                                    for(k in 2:nrow(tseq)){
                                      if(tseq$time[k] > temp$x[5]){
                                        tseq$y[k] <- temp$y[5]
                                      }
                                    }
                                    
                                # Now a loess curve is fit based on that 5 minute time series with smoothing and the predicted values
                                  # for a one minute time course are saved.
                                    suppressWarnings(
                                      fits <- predict(loess(tseq$y ~ tseq$time, span = sm_span), newdata = time_seq$time)
                                    )
                                    preds[i, ] <- fits
                            }
                        # Convert any negative numbers to 0
                          preds[preds < 0] <- 0
            
            # Make simulated data set sampled just at the bleed times
                    preds2 <- preds[, bleed_times + 1]
                    preds2$animal_sample <- rownames(preds2)
                    preds2 <- plyr::join(preds2, data2[, c("animal_sample", "animal", "sample", "performance")], "animal_sample", "left", "first")
                    preds2_long <- as.data.frame(pivot_longer(preds2, cols = seq(2, ncol(preds2) - 3), values_to = "cort", names_to = "time"))
                    preds2_long$time <- as.numeric(preds2_long$time)
                  # Add in assay error. To avoid this set to 0 in function call.
                    for(i in 1:nrow(preds2_long)){
                      error_amt <- assay_error * preds2_long$cort[i]
                      preds2_long$cort[i] <- preds2_long$cort[i] + runif(1, - error_amt, error_amt)
                    }
    
            # Make full time course. This is just data wrangling from the object created above.
                  pred_time <- preds
                  for(j in 1:nrow(pred_time)){
                    temp <- na.omit(t(pred_time[j, ]))[length(na.omit(t(pred_time[j, ])))]
                    replace <- t(is.na(pred_time[j, ]) == TRUE)
                    pred_time[j, replace] <- temp
                  }
                  pred_time$animal_sample <- rownames(pred_time)
                  pred_time <- plyr::join(pred_time, data2[, c("animal_sample", "animal", "sample", "performance")], "animal_sample", "left", "first")
                  pred_time_long <- as.data.frame(pivot_longer(data = pred_time, cols = seq(2, ncol(pred_time) - 3), values_to = "cort",
                                                 names_to = "time"))
                  pred_time_long$time <- as.numeric(pred_time_long$time)
              
            # make rank time course
                pred_rank <- as.data.frame(preds)
                pred_rank[] <- lapply(-pred_rank, rank, ties.method = "min")
                pred_rank <- as.data.frame(t(pred_rank))
                pred_rank$time <- time_seq$time
                long_rank <- as.data.frame(pivot_longer(pred_rank, cols = starts_with("id"), names_to = "animal_sample", values_to = "rank"))
                n_s <- ncol(pred_rank) - 1
                colors <- long_rank[1:n_s, ]
                colors <- colors[order(colors$rank), ]  
                colors$virid <- viridis(n_s)
                colors <- colors[, c("animal_sample", "virid")]
                long_rank <- plyr::join(long_rank, colors, "animal_sample", "left", "all")
                long_rank <- long_rank[order(long_rank$time, long_rank$rank), ]
                long_rank$animal_sample2 <- factor(long_rank$animal_sample,
                                            levels = long_rank$animal_sample[1:n_s])
                
            # True values
                true <- as.data.frame(pivot_wider(data$sim_dat2, names_from = group, values_from = c(x, y)))
                true$baseline <- true$y_initial
                true$slope <- true$x_slopeturn
                true$fastpct <- (true$y_slopeturn - true$y_initial) / (true$y_reachmax - true$y_initial)
                true$speed <- true$x_reachmax
                true$max <- true$y_reachmax
                true$atmax <- true$x_startdecline - true$x_reachmax
                true$return <- true$x_return - true$x_startdecline
                true <- true[, c("animal", "baseline", "slope", "fastpct", "speed", "max", "atmax", "return")]
                true <- plyr::join(true, preds2_long[, c("animal", "performance")], "animal", "left", "first")
                
            # Calculate area under the curve for the true values
                for(u in 1:nrow(true)){
                    box_coords <- matrix(c(0, 0,
                                            0, true$baseline[u],
                                           true$slope[u], true$fastpct[u] * true$max[u],
                                           true$speed[u], true$max[u],
                                           true$speed[u] + true$atmax[u], true$max[u],
                                           true$return[u] + true$speed[u] + true$baseline[u], true$baseline[u],
                                           true$return[u] + true$speed[u] + true$baseline[u], 0,
                                           0, 0),
                                          nrow = 8, ncol = 2, byrow = TRUE)
                    box_coords2 <- box_coords[c(2:6, 2), ]
                    true$AUCg[u] <- Polygon(box_coords, hole = FALSE)@area / (true$return[u] + true$speed[u] + true$baseline[u])
                    true$AUCi[u] <- Polygon(box_coords2, hole = FALSE)@area / (true$return[u] + true$speed[u] + true$baseline[u])
                }
                
            # Calculate area under the curve for observed curves
                obs_AUC <- data.frame(animal_sample = unique(preds2_long$animal_sample))
                for(u in 1:nrow(obs_AUC)){
                  sub <- subset(preds2_long, preds2_long$animal_sample == obs_AUC$animal_sample[u])
                  box1 <- matrix(nrow = nrow(sub) + 2, ncol = 2)
                  box1[1, ] <- c(sub$time[1], sub$cort[1])
                  for(k in 2:nrow(sub)){
                    box1[k, ] <- c(sub$time[k], sub$cort[k])
                  }
                  box1[nrow(sub) + 1, ] <- c(sub$time[nrow(sub)], sub$cort[1])
                  box1[nrow(sub) + 2, ] <- box1[1, ]
                  
                  box2 <- box1
                  box2[1, 2] <- 0
                  box2[nrow(box2) - 1, 2] <- 0
                  box2[nrow(box2), 2] <- 0
                  
                  obs_AUC$sample_AUCi[u] <- Polygon(box1, hole = FALSE)@area / (sub$time[nrow(sub)] - sub$time[1])
                  obs_AUC$sample_AUCg[u] <- Polygon(box2, hole = FALSE)@area / (sub$time[nrow(sub)] - sub$time[1])
                }
                
                for(u in 1:nrow(obs_AUC)){
                  sub <- subset(pred_time_long, pred_time_long$animal_sample == obs_AUC$animal_sample[u])
                  box1 <- matrix(nrow = nrow(sub) + 2, ncol = 2)
                  box1[1, ] <- c(sub$time[1], sub$cort[1])
                  for(k in 2:nrow(sub)){
                    box1[k, ] <- c(sub$time[k], sub$cort[k])
                  }
                  box1[nrow(sub) + 1, ] <- c(sub$time[nrow(sub)], sub$cort[1])
                  box1[nrow(sub) + 2, ] <- box1[1, ]
                  
                  box2 <- box1
                  box2[1, 2] <- 0
                  box2[nrow(box2) - 1, 2] <- 0
                  box2[nrow(box2), 2] <- 0
                  
                  obs_AUC$full_AUCi[u] <- Polygon(box1, hole = FALSE)@area / (sub$time[nrow(sub)] - sub$time[1])
                  obs_AUC$full_AUCg[u] <- Polygon(box2, hole = FALSE)@area / (sub$time[nrow(sub)] - sub$time[1])
                  
                  sub2 <- subset(sub, sub$time < max(preds2_long$time) & sub$time > min(preds2_long$time))
                  
                  box1 <- matrix(nrow = nrow(sub2) + 2, ncol = 2)
                  box1[1, ] <- c(sub2$time[1], sub2$cort[1])
                  for(k in 2:nrow(sub2)){
                    box1[k, ] <- c(sub2$time[k], sub2$cort[k])
                  }
                  box1[nrow(sub2) + 1, ] <- c(sub2$time[nrow(sub2)], sub2$cort[1])
                  box1[nrow(sub2) + 2, ] <- box1[1, ]
                  
                  box2 <- box1
                  box2[1, 2] <- 0
                  box2[nrow(box2) - 1, 2] <- 0
                  box2[nrow(box2), 2] <- 0
                  
                  obs_AUC$subfull_AUCi[u] <- Polygon(box1, hole = FALSE)@area / (sub$time[nrow(sub2)] - sub2$time[1])
                  obs_AUC$subfull_AUCg[u] <- Polygon(box2, hole = FALSE)@area / (sub$time[nrow(sub2)] - sub2$time[1])
                }
                
                obs_AUC <- plyr::join(obs_AUC, pred_time_long[, c("animal_sample", "animal", "sample")], "animal_sample", "left", "first")
                true_tem <- true[, c("animal", "AUCg", "AUCi")]
                colnames(true_tem)[2:3] <- c("true_AUCg", "true_AUCi")
                obs_AUC <- plyr::join(obs_AUC, true_tem, "animal", "left")
              
            # output all five dataframes plus the initial 'true' values
              dataset2 <- list(simulated_dataset_long = preds2_long,
                              timecourse_long = pred_time_long,
                              rank_timecourse = long_rank,
                              AUC_measures = obs_AUC,
                              true_values = true)
              return(dataset2)
          }
          
    ## Repeatability functions
          # This function will take a simulated dataset where each individual is observed at least two times and calculate repeatability
            # of cort for each time measure, for AUCi & AUCg (from full simulated timecourse and from sub-sampled points),
            # and profile repeatability (Reed & Romero). 
        
        # This is the main function that will calculate repeatability for several different AUCs and for each timepoint as well
          # as producing three plots that are saved in a list.
            cort_repeat <- function(data = cort_sim2(sample_times = 6),       # data input from cort simulation functions above
                                    boots = 100,                              # number of bootstraps for repeatability calculations
                                    perms = 0                                 # number of permutations for repeatability calculations; 0 = NA
                                    )
              {
              
              # Calculate profile repeatability as described in Reed & Romero 2019. This is a really tedious calculation that involves
                # counting the number of times lines cross etc. The code could probably be cleaner. One PR is calculated for each individual
                # and then an average can be taken for all of them together.
              
                    # How many times are cort values saved at
                        times <- length(unique(data$simulated_dataset_long$time))
                        
                    # new data frame with one row for each animal
                        pr_animal <- data.frame(animal = unique(data$simulated_dataset_long$animal))
                        
                    # Loop through each animal and calculate PR
                            for(u in 1:nrow(pr_animal)){
                                    sub <- subset(data$simulated_dataset_long, data$simulated_dataset_long$animal == pr_animal$animal[u])
                                    pr_animal$possible_cross[u] <- sum(seq(1:(length(unique(sub$sample)) - 1))) * (times - 1)
                                    obs_cross <- 0
                                
                                  # data wrangling to count the number of crosses observed
                                    sub2 <- as.data.frame(pivot_wider(sub, names_from = time, values_from = cort))
                                    subx <- sub2
                                    colnames(subx) <- paste0(colnames(sub2), "_rep")
                                    colnames(subx)[2] <- "animal"
                                    sub3 <- plyr::join(sub2, subx, "animal", "left")
                                    sub3 <- subset(sub3, sub3$sample != sub3$sample_rep)
                                    for(k in 1:nrow(sub3)){
                                      or <- order(c(sub3$sample[k], sub3$sample_rep[k]))
                                      sam <- c(sub3$sample[k], sub3$sample_rep[k])
                                      sub3$pair[k] <- paste(sam[or[1]], sam[or[2]], sep = "_")
                                    }
                                    sub3 <- subset(sub3, !duplicated(sub3$pair))
                                    
                                    for(w in 1:nrow(sub3)){
                                      for(t in 1:(times - 1)){
                                        if(sub3[w, t + 4] > sub3[w, times + 7 + t]){
                                          if(sub3[w, t + 5] < sub3[w, times + 8 + t]){obs_cross <- obs_cross + 1}
                                        }
                                        if(sub3[w, t + 4] < sub3[w, times + 7 + t]){
                                          if(sub3[w, t + 5] > sub3[w, times + 8 + t]){obs_cross <- obs_cross + 1}
                                        }
                                      }  
                                    }  
                                    pr_animal$obs_cross[u] <- obs_cross 
                                    
                                # calculting the other numbers needed for PR
                                    sam_times <- unique(data$simulated_dataset_long$time)
                                    vars <- rep(NA, length(sam_times))
                                    for(m in 1:length(sam_times)){
                                      subvar <- subset(sub, sub$time == sam_times[m])
                                      vars[m] <- var(subvar$cort)
                                    }
                                    
                                    pr_animal$maxvar[u] <- max(vars)
                                    pr_animal$avevar[u] <- mean(vars)
                    
                            }
                      
                        # Saving the numbers and making final PR calculation for each animal
                        # If you want to save PR for each animal need to modify to return this 'pr_animal' object
                            pr_animal$cross_score <- pr_animal$maxvar * round(10 * pr_animal$obs_cross / pr_animal$possible_cross) / 5
                            pr_animal$base <- ((pr_animal$maxvar + pr_animal$avevar + pr_animal$cross_score) / 100) - 5
                            pr_animal$PR <- 1 - (1 / (1 + exp(- pr_animal$base)))
                            
                        # Avarage PR for the whole dataset
                            avg_PR <- mean(pr_animal$PR)
                            
              # Calculate repeatability of the other metrics using rptR and simple LMMs
                    # Repeatability of AUC measures
                            suppressWarnings(
                              r_sample_AUCi <- rpt(sample_AUCi ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            suppressWarnings(
                              r_sample_AUCg <- rpt(sample_AUCg ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            suppressWarnings(
                              r_full_AUCi <- rpt(full_AUCi ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            suppressWarnings(
                              r_full_AUCg <- rpt(full_AUCg ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            suppressWarnings(
                              r_subfull_AUCi <- rpt(subfull_AUCi ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            suppressWarnings(
                              r_subfull_AUCg <- rpt(subfull_AUCg ~ (1|animal), grname = "animal", data = data$AUC_measures, nboot = boots, npermut = perms))
                            
              # Calculate repeatability for measures at each timepoint
                  temp <- as.data.frame(pivot_wider(data$simulated_dataset_long, names_from = time, values_from = cort))
                  colnames(temp)[5:(ncol(temp))] <- paste0("time", colnames(temp)[5:(ncol(temp))])
                  save_rs <- list(NA)
                  for(a in 1:(ncol(temp) - 4)){
                    temp$response <- temp[, a + 4]
                    save_rs[[a]] <- rpt(response ~ (1|animal), grname = "animal", data = temp, nboot = boots, npermut = perms, datatype = "Gaussian")
                  }
                  
                  temp_rs <- data.frame(repeatability = rep(NA, length(save_rs)), lowCI = rep(NA, length(save_rs)), 
                                        highCI = rep(NA, length(save_rs)), pLRT = rep(NA, length(save_rs)), pPERM = rep(NA, length(save_rs)))
                  for(t in 1:length(save_rs)){
                    temp_rs$repeatability[t] <- save_rs[[t]]$R
                    temp_rs$lowCI[t] <- save_rs[[t]]$CI_emp[1]
                    temp_rs$highCI[t] <- save_rs[[t]]$CI_emp[2]
                    temp_rs$pLRT[t] <- save_rs[[t]]$P[1]
                    temp_rs$pPERM[t] <- save_rs[[t]]$P[2]
                  }
                  
              # Build a data frame with all the repeatability measures
                  cort_repeat <- data.frame(measure = c("PR", "sample_AUCi", "sample_AUCg", "full_AUCi", "full_AUCg", "subfull_AUCi", "subfull_AUCg",
                                                        colnames(temp)[5:(ncol(temp) - 1)]),
                                            repeatability = c(avg_PR, unlist(r_sample_AUCi$R), unlist(r_sample_AUCg$R), 
                                                              unlist(r_full_AUCi$R), unlist(r_full_AUCg$R), 
                                                              unlist(r_subfull_AUCi$R), unlist(r_subfull_AUCg$R), unlist(temp_rs$repeatability[1:nrow(temp_rs)])),
                                            lowCI = c(NA, unlist(r_sample_AUCi$CI_emp[1]), unlist(r_sample_AUCg$CI_emp[1]), 
                                                                unlist(r_full_AUCi$CI_emp[1]), unlist(r_full_AUCg$CI_emp[1]), 
                                                                unlist(r_subfull_AUCi$CI_emp[1]), unlist(r_subfull_AUCg$CI_emp[1]), unlist(temp_rs$lowCI[1:nrow(temp_rs)])),
                                            highCI = c(NA, unlist(r_sample_AUCi$CI_emp[2]), unlist(r_sample_AUCg$CI_emp[2]), 
                                                      unlist(r_full_AUCi$CI_emp[2]), unlist(r_full_AUCg$CI_emp[2]), 
                                                      unlist(r_subfull_AUCi$CI_emp[2]), unlist(r_subfull_AUCg$CI_emp[2]), unlist(temp_rs$highCI[1:nrow(temp_rs)])),
                                            pLRT = c(NA, unlist(r_sample_AUCi$P[1]), unlist(r_sample_AUCg$P[1]), 
                                                      unlist(r_full_AUCi$P[1]), unlist(r_full_AUCg$P[1]), 
                                                      unlist(r_subfull_AUCi$P[1]), unlist(r_subfull_AUCg$P[1]), unlist(temp_rs$pLRT[1:nrow(temp_rs)])),
                                            pPERM = c(NA, unlist(r_sample_AUCi$P[2]), unlist(r_sample_AUCg$P[2]), 
                                                      unlist(r_full_AUCi$P[2]), unlist(r_full_AUCg$P[2]), 
                                                      unlist(r_subfull_AUCi$P[2]), unlist(r_subfull_AUCg$P[2]), unlist(temp_rs$pPERM[1:nrow(temp_rs)]))
                                            )
                            
              # Make three different plots
                  p1 <- ggplot(data$simulated_dataset_long, mapping = aes(x = time, y = cort, by = animal_sample, color = animal)) +
                    geom_line() + theme_classic() + scale_color_viridis(discrete = TRUE) +
                    guides(color = FALSE)
                  
                  p2 <- p1 + facet_wrap(~ animal)
                  
                  p3 <- ggplot(data = data$simulated_dataset_long, mapping = aes(x = animal, y = cort, color = animal, fill = animal)) +
                    geom_boxplot(alpha = 0.2) + scale_fill_viridis(discrete = TRUE) + 
                    geom_point() + facet_wrap(~time, scales = "free") + scale_color_viridis(discrete = TRUE) + 
                    guides(color = FALSE, fill = FALSE) +
                    theme_bw() + theme(axis.text.x = element_text(angle = 90))
                  
              # Combine plot and data
                  output <- list(repeat_est = cort_repeat,
                                 plot_all = p1,
                                 plot_panels = p2,
                                 plot_dots = p3)
                  
                  return(output)
                            
            }
            
      
    # Plot functions. The datasets generated in cort_sim2 can be used to make any number of different plots or models. But this
        # function creates a three panel summary plot that illustrates the data. It's only input is a list with the three data
        # frames output by 'cort_sim2'. It will by default generate a new dataset when called if one isn't provided.
          plot_cort_sim <- function(data = cort_sim2())
              {
            
            spoints <- unique(data$simulated_dataset_long$time)
            p1 <- ggplot(data = data$simulated_dataset, mapping = aes(x = time, y = cort, by = animal_sample)) +
              geom_line(size = 0.5, color = "coral3", alpha = 0.7) + 
              guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") +
              geom_vline(xintercept = spoints, linetype = "dashed")
            p2 <- ggplot(data = data$timecourse_long, mapping = aes(x = time, y = cort, color = animal, by = animal_sample)) +
              #stat_smooth(geom = "line", method = "loess", span = 0.3, alpha = 0.7, se = FALSE, color = "coral3") + 
              geom_line(size = 0.5, color = "coral3", alpha = 0.7) + 
              guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") + coord_cartesian(xlim = c(0, 60)) +
              geom_vline(xintercept = spoints, linetype = "dashed")
            p3 <- ggplot(data = data$rank_timecourse, mapping = aes(x = time, y = rank, color = animal_sample2)) +
              geom_line(size = 1.6, alpha = 0.8) + guides(color = FALSE) + theme_classic() + 
              scale_color_viridis(discrete = TRUE) + xlim(0, 60) +
              geom_vline(xintercept = spoints, linetype = "dashed")
            grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 2, 2), c(3, 3, 3)))
          }
          plot_cort_sim_1 <- function(data = cort_sim2(), x_max = 60, y_max = 85)
              {
            spoints <- unique(data$simulated_dataset_long$time)
            p <- ggplot(data = data$timecourse_long, mapping = aes(x = time, y = cort, by = animal_sample)) +
              #stat_smooth(geom = "line", method = "loess", span = 0.7, size = 0.5, alpha = 0.6, se = FALSE, color = "coral3") + 
              geom_line(size = 0.5, alpha = 0.6, color = "coral3") +
              guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") + coord_cartesian(xlim = c(0, x_max)) +
              ylim(0, y_max) #+
              #geom_vline(xintercept = spoints, linetype = "dashed")
            p
          }


## Troubleshooting ----      
    # testing span. The function smooths the 

          
## Initial demonstration ----
      set.seed(100)   # makes it reproducible by using same random seed
      demo <- cort_sim2()
      p <- plot_cort_sim(demo)
      ggsave(here::here("3_r_scripts/demo.png"), p, device = "png", width = 10.2, height = 7.5)
      
      ggsave(here::here("3_r_scripts/demo_fullx.png"), 
        plot_cort_sim_1(demo, x_max = 160, y_max = 100),
        device = "png", width = 9, height = 5, units = "in")
      
## Variation in speed vs. scope ----
    set.seed(955)
    d1 <- cort_sim2(cort_sim1(speed_mu = 25, speed_sd = 0, max_sd = 0.2))  
    d2 <- cort_sim2(cort_sim1(max_sd = 0.2))
    d3 <- cort_sim2(cort_sim1(max_sd = 0, speed_sd = 10))
    p1 <- plot_cort_sim_1(d1)
    p1 <- p1 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "A")
    p2 <- plot_cort_sim_1(d2)
    p2 <- p2 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "B")
    p3 <- plot_cort_sim_1(d3)
    p3 <- p3 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "C")
    ggsave(here::here("3_r_scripts/relative_variation.png"),
           ggpubr::ggarrange(p1, p2, p3, nrow = 1),
           device = "png", width = 9.2, height = 3, units = "in")
                    
## Simulate data ----
    #plot_cort_sim(cort_sim2(cort_sim1(n = 60, speed_mu = 18, speed_sd = 3, base_mu = 9, base_sd = 3, cor_base_max = 0.3, cor_base_speed = -0.5, cor_speed_max = -0.7, cor_max_return = 0.4), 
    #                        bleed_times = c(1, 15, 30), sample_times = 1))
          
    # Speed vs. scope illustration
          d1 <- cort_sim2(cort_sim1(cor_speed_max = -0.9, speed_sd = 9, return_sd = 0, maxtime_sd = 0, slope_sd = 2, fastpct_sd = .1))
          d2 <- cort_sim2(cort_sim1(cor_speed_max = 0, speed_sd = 9, return_sd = 0, maxtime_sd = 0, slope_sd = 2, fastpct_sd = .1))
          d3 <- cort_sim2(cort_sim1(cor_speed_max = 0.9, speed_sd = 9, return_sd = 0, maxtime_sd = 0, slope_sd = 2, fastpct_sd = .1))
          p1 <- plot_cort_sim_1(data = d1)
          p1 <- p1 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "A")
          p2 <- plot_cort_sim_1(data = d2)
          p2 <- p2 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "B")
          p3 <- plot_cort_sim_1(data = d3)
          p3 <- p3 + annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "C")
          
          
        ggsave(here::here("3_r_scripts/sp_sc_correlation.png"),
          ggpubr::ggarrange(p1, p2, p3, nrow = 1),
          device = "png", width = 9.2, height = 3, units = "in")
          
    # Only variation in speed. None in max
          set.seed(82)
          d1 <- cort_sim2(cort_sim1(max_sd = 0, speed_mu = 30, speed_sd = 3, base_sd = 0.1), sample_times = 1, sm_span = 0.2)
          d2 <- cort_sim2(cort_sim1(max_sd = 0, speed_mu = 30, speed_sd = 6, base_sd = 0.1), sample_times = 1, sm_span = 0.2)
          d3 <- cort_sim2(cort_sim1(max_sd = 0, speed_mu = 30, speed_sd = 12, base_sd = 0.1), sample_times = 1, sm_span = 0.2)
          p1 <- plot_cort_sim_1(data = d1)
          p2 <- plot_cort_sim_1(data = d2)
          p3 <- plot_cort_sim_1(data = d3)
          
          ggsave(here::here("3_r_scripts/var_in_speed_only.png"),
                 ggpubr::ggarrange(p1, p2, p3, nrow = 1),
                 device = "png", width = 9.2, height = 3, units = "in")
          
    
## Simulate RWBL data ----
     rwde <- subset(rwd_long, rwd_long$LH.stage == "early-breeding")  
      rwde_wide <- subset(rw_dat, rw_dat$LH.stage == "early-breeding")
      rwde <- subset(rwde, rwde$band != "1372-33452")   # This one has an unusually high cort response
    set.seed(165)
     sim_rwd <- cort_sim2(cort_sim1(
       n = 23,
       base_mu = 3.5,
       base_sd = 3.5,
       slope_mu = 17,
       slope_sd = 0.7,
       fastpct_mu = 0.7,
       fastpct_sd = 0.03,
       max_mu = log(56),
       max_sd = 0.5,
       speed_mu = 65,
       speed_sd = 5,
       return_mu = 120,
       maxtime_mu = 25,
       cor_max_slope = -0.7
     ), 
     sample_times = 1,
     bleed_times = c(1, 3, 8, 15, 30, 45, 60),
     assay_error = 0.2,
     sm_span = 0.2)         
     
     
     sim_rwd$simulated_dataset_long$type <- "simulated"
     sim_rwd2 <- sim_rwd$simulated_dataset_long[, c("animal_sample", "time", "cort", "type")]
     rwde2 <- rwde[, c("band", "tsec", "cort")]
     rwde2$tsec <- rwde2$tsec / 60
     rwde2$type <- "red-winged blackbird"
     colnames(rwde2) <- c("animal_sample", "time", "cort", "type")
     
     comp_sim <- rbind(rwde2, sim_rwd2)
     
     p1 <- ggplot(comp_sim, mapping = aes(x = time, y = cort, by = animal_sample, color = type)) +
       geom_line(alpha = 0.7, size = 0.7) +
       theme_classic() +
       coord_cartesian(xlim = c(0, 60), ylim = c(0, 125)) +
       xlab("Time (minutes)") + ylab("Corticosterone (ng/ul)") +
       scale_color_manual(values = c("slateblue", "orange")) +
       guides(color = FALSE) +
       annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "A")
     
     p2 <- ggplot(comp_sim, mapping = aes(x = time, y = cort, color = type, fill = type)) +
       geom_smooth() +
       theme_classic() + 
       coord_cartesian(xlim = c(0, 60), ylim = c(0, 125)) +
       xlab("Time (minutes)") + ylab("Corticosterone (ng/ul)") +
       scale_color_manual(values = c("slateblue", "orange")) +
       scale_fill_manual(values = c("slateblue", "orange")) +
       theme(legend.position = c(0.8, 0.15)) +
       guides(color = guide_legend(title = ""), fill = guide_legend(title = "")) +
       annotate(geom = "text", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, label = "B")

      ggpubr::ggarrange(p1, p2, nrow = 1)     
     
     ggsave(here::here("3_r_scripts/rwbl_sim.png"),
            ggpubr::ggarrange(p1, p2, nrow = 1),
            device = "png", units = "in", width = 9, height = 4)
     
     
## Scenario 1
     sc1 <- cort_sim2(cort_sim1(
       n = 100, base_mu = 4, base_sd = 3, slope_mu = 18, slope_sd = 0.5, fastpct_mu = 0.78, max_mu = log(61),
       max_sd = 0.3, speed_mu = 60, speed_sd = 5, return_mu = 120, maxtime_mu = 25, cor_max_slope = 0
     ), 
     sample_times = 1, bleed_times = c(15, 30, 60), assay_error = 0.2, sm_span = 0.2,
     performance_contributions = c(0, 0, 1, 0, 0, 0, 0, 3))   
     
     # Relative contributions of base, speed, max, maxtime, return, slope, fastpct, and random error
     
     sc1d <- subset(sc1$simulated_dataset_long, sc1$simulated_dataset_long$time == 15)
     
     ggplot(sc1d, mapping = aes(x = cort, y = performance)) + 
       geom_point() + 
       geom_smooth(method = "lm") +
       theme_classic() + xlab("Corticosterone at 30 minutes") +
       ylab("Performance")
     ggplot(sc1$true_values, mapping = aes(x = max, y = 1 * performance)) + 
       geom_point() +
       theme_classic() +
       geom_smooth(method = "lm")
     
          
     
     
     