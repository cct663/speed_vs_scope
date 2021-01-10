### Script for speed vs. scope manuscript
# Written by Conor Taff; cct63@cornell.edu; cct663@gmail.com


## Load packages ----
    pacman::p_load(ggplot2, here, MASS, MBESS, tidyverse, viridis, gridExtra)

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
                               base_mu = 5, base_sd = 1.5, base_min = 0.5,
                               speed_mu = 30, speed_sd = 5, speed_min = 10,
                               max_mu = log(50), max_sd = log(1.2), 
                               maxtime_mu = 10, maxtime_sd = 3, maxtime_min = 4,
                               return_mu = 90, return_sd = 15, 
                               cor_base_speed = 0, cor_base_max = 0, cor_base_maxtime = 0, cor_base_return = 0,
                               cor_speed_max = 0, cor_speed_maxtime = 0, cor_speed_return = 0,
                               cor_max_maxtime = 0, cor_max_return = 0,
                               cor_maxtime_return = 0)
              {
                  # make variance covariance matrix
                      cort_vcov <- matrix(nrow = 5, ncol = 5)
                      cort_vcov[1, ] <- c(base_sd^2, cor_base_speed * base_sd * speed_sd, cor_base_max * base_sd * max_sd, 
                                          cor_base_maxtime * base_sd * maxtime_sd, cor_base_return * base_sd * return_sd)
                      cort_vcov[2, ] <- c(cor_base_speed * base_sd * speed_sd, speed_sd^2, cor_speed_maxtime * speed_sd * max_sd, 
                                          cor_speed_maxtime * speed_sd * maxtime_sd, cor_speed_return * speed_sd * return_sd)
                      cort_vcov[3, ] <- c(cor_base_max * base_sd * max_sd, cor_speed_max * speed_sd * max_sd, max_sd^2, 
                                          cor_max_maxtime * max_sd * maxtime_sd, cor_max_return * max_sd * return_sd)
                      cort_vcov[4, ] <- c(cor_base_maxtime * base_sd * maxtime_sd, cor_speed_maxtime * speed_sd * maxtime_sd, 
                                          cor_max_maxtime * max_sd * maxtime_sd, maxtime_sd^2, cor_maxtime_return * maxtime_sd * return_sd)
                      cort_vcov[5, ] <- c(cor_base_return * base_sd * return_sd, cor_speed_return * speed_sd * return_sd, 
                                          cor_max_return * max_sd * return_sd, cor_maxtime_return * maxtime_sd * return_sd, return_sd^2)
                      
                      mu_list <- c(base_mu, speed_mu, max_mu, maxtime_mu, return_mu)
                  
                  # Create a dataset
                        msam <- mvrnorm(n, mu = mu_list, Sigma = cort_vcov)
                      
                      sim_dat <- data.frame(
                        animal = paste("id", 1:n, sep = "_"),
                        base = msam[, 1],
                        tmax = msam[, 2],
                        max = exp(msam[, 3]),
                        atmax = msam[, 4],
                        return = msam[, 5]
                      )
                      
                # calculate return time as return latency plus speed plus time at max
                      sim_dat$return <- sim_dat$tmax + sim_dat$atmax + sim_dat$return
                          
                      # replace base, maxtime, and speed values with minimum allowed (avoids negative numbers)
                          sim_dat[which(sim_dat$base < base_min), "base"] <- base_min
                          sim_dat[which(sim_dat$tmax < speed_min), "tmax"] <- speed_min
                          sim_dat[which(sim_dat$atmax < maxtime_min), "atmax"] <- maxtime_min
                  
                  # Convert simulated data to long form x/y dataframe
                        sim_dat2 <- data.frame(
                          animal = rep(sim_dat$animal, 4),
                          x = c(rep(0, n),
                                sim_dat$tmax,
                                sim_dat$atmax + sim_dat$tmax,
                                sim_dat$return),
                          y = c(sim_dat$base,
                                sim_dat$max,
                                sim_dat$max,
                                sim_dat$base),
                          group = c(rep("initial", n),
                                    rep("reachmax", n),
                                    rep("startdecline", n),
                                    rep("return", n))
                        )
                        return(sim_dat2)
      }
      
    # The second function takes input from the first function (or by default, generates a dataset with the first function).
        # These parameters are used to fit and estimate full time courses of cort responses for each animal. Noise is added
        # to each of the 'true' parameters and the amount of noise can be specified. This can be considered measurement error
        # or an individual deviation from their 'true' parameters (or a combination of both). The same animal can be sampled
        # multiple times, with noise added separately each time. The output is a list with three dataframes that represent
        # a simulated dataset with samples at specified time points, the full time course with cort at every time point, and a 
        # rank converted estimate. 
          cort_sim2 <- function(data = cort_sim1(),
                                bleed_times = c(2, 15, 30),
                                base_error = 0.1,
                                speed_error = 0.1,
                                max_error = 0.1,
                                maxtime_error = 0.1,
                                return_error = 0.1,
                                sample_times = 1,
                                timecourse_max = 170,
                                performance_contributions = c(0, 0, 0, 0, 0, 1),
                                sm_span = .6) #relative contribution of base, speed, max, maxtime, return, and error to performance
              {
            data2 <- data[rep(seq_len(nrow(data)), sample_times), ]
            data2$sample <- rep(seq(1, sample_times, 1), each = nrow(data))
            data2$animal_sample <- paste(data2$animal, data2$sample, sep = "_")
            data3 <- as.data.frame(pivot_wider(data2, id_cols = animal_sample, names_from = group, values_from = c(x, y)))
            per_conts <- performance_contributions / sum(performance_contributions)
            data3$maxtime <- data3$x_startdecline - data3$x_reachmax
            data3$maxtime <- scale(data3$maxtime)
            data3$y_initial <- scale(data3$y_initial)
            data3$x_reachmax <- scale(data3$x_reachmax)
            data3$y_reachmax <- scale(data3$y_reachmax)
            data3$x_return <- scale(data3$x_return)
            for(k in 1:nrow(data3)){
              data3$performance[k] <- sum(per_conts[1:4] * data3[k, c("y_initial", "x_reachmax", "y_reachmax", "maxtime", "x_return")]) 
                    + rnorm(1, mean = 0, sd = 1) * per_conts[6]
            }
            data2 <- plyr::join(data2, data3[, c("animal_sample", "performance")], "animal_sample", "left", "first")
            
            # Add noise (could be individual variation or measurement error)
              # Noise for base cort
                data2[which(data2$group == "initial"), "y"] <- data2[which(data2$group == "initial"), "y"] +
                  rnorm(n = length(data2[which(data2$group == "initial"), "y"]),
                        mean = 0,
                        sd = base_error * data2[which(data2$group == "initial"), "y"])
              # Noise for speed to response
                data2[which(data2$group == "reachmax"), "x"] <- data2[which(data2$group == "reachmax"), "x"] +
                  rnorm(n = length(data2[which(data2$group == "reachmax"), "x"]),
                        mean = 0,
                        sd = speed_error * data2[which(data2$group == "reachmax"), "x"])
              # Noise for maxcort
                data2[which(data2$group == "reachmax"), "y"] <- data2[which(data2$group == "reachmax"), "y"] +
                  rnorm(n = length(data2[which(data2$group == "reachmax"), "y"]),
                        mean = 0,
                        sd = max_error * data2[which(data2$group == "reachmax"), "y"])
              # Noise for max time
                data2[which(data2$group == "startdecline"), "x"] <- data2[which(data2$group == "startdecline"), "x"] +
                  rnorm(n = length(data2[which(data2$group == "startdecline"), "x"]),
                        mean = 0,
                        sd = maxtime_error * data2[which(data2$group == "startdecline"), "x"])
              # Noise for return
                data2[which(data2$group == "return"), "x"] <- data2[which(data2$group == "return"), "x"] +
                  rnorm(n = length(data2[which(data2$group == "return"), "x"]),
                        mean = 0,
                        sd = return_error * data2[which(data2$group == "return"), "x"])
             
            
            time_seq <- data.frame(time = seq(0, timecourse_max, 1))
            preds <- as.data.frame(matrix(nrow = length(unique(data2$animal_sample)), ncol = nrow(time_seq)))
            colnames(preds) <- time_seq$time
            rownames(preds) <- unique(data2$animal_sample)
            for(i in 1:length(unique(data2$animal_sample))){
              temp <- subset(data2, data2$animal_sample == unique(data2$animal_sample)[i])
              suppressWarnings(
                fits <- predict(loess(temp$y ~ temp$x, span = sm_span), newdata = time_seq$time)
              )
              preds[i, ] <- fits
            }
            
            # make simulated dataset
              preds2 <- preds[, bleed_times + 1]
              preds2$animal_sample <- rownames(preds2)
              preds2 <- plyr::join(preds2, data2[, c("animal_sample", "animal", "sample", "performance")], "animal_sample", "left", "first")
              preds2_long <- as.data.frame(pivot_longer(preds2, cols = seq(2, ncol(preds2) - 3), values_to = "cort", names_to = "time"))
              preds2_long$time <- as.numeric(preds2_long$time)
    
            # make full time course
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
              
            # output all three dataframes  
              dataset <- list(simulated_dataset_long = preds2_long,
                              timecourse_long = pred_time_long,
                              rank_timecourse = long_rank,
                              true_values = data)
              return(dataset)
      }
      
    # Plot function. The datasets generated in cort_sim2 can be used to make any number of different plots or models. But this
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
          plot_cort_sim_1 <- function(data = cort_sim2(), x_max = 60){
            spoints <- unique(data$simulated_dataset_long$time)
            p <- ggplot(data = data$timecourse_long, mapping = aes(x = time, y = cort, by = animal_sample)) +
              #stat_smooth(geom = "line", method = "loess", span = 0.7, size = 0.5, alpha = 0.6, se = FALSE, color = "coral3") + 
              geom_line(size = 0.5, alpha = 0.6, color = "coral3") +
              guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") + coord_cartesian(xlim = c(0, x_max)) +
              geom_vline(xintercept = spoints, linetype = "dashed")
            p
          }


## Troubleshooting ----      
    # testing span. The function smooths the 
          plot_cort_sim()
          
          set.seed(95)
          tt <- cort_sim2(cort_sim1(n = 5))
          plot_cort_sim(data = tt)
          
          ggplot(data = tt$timecourse_long, mapping = aes(x = time, y = cort, color = animal_sample)) +
            #stat_smooth(geom = "line", method = "loess", span = 0.3, alpha = 0.7, se = FALSE, color = "coral3") + 
            geom_line(size = 0.5, alpha = 0.7) + 
            guides(color = FALSE) +
            theme_classic() + xlab("Time") + ylab("Corticosterone") + coord_cartesian(xlim = c(0, 120)) +
            geom_vline(xintercept = spoints, linetype = "dashed") + geom_point() +
            geom_point(data = tt$true_values, mapping = aes(x = x, y = y, shape = animal), color = "blue")
          
## Initial demonstration ----
      set.seed(100)   # makes it reproducible by using same random seed
      demo <- cort_sim2()
      p <- plot_cort_sim(demo)
      ggsave(here::here("3_r_scripts/demo.png"), p, device = "png", width = 10.2, height = 7.5)
      
## Variation in speed vs. scope ----
    set.seed(955)
    d1 <- cort_sim2(cort_sim1(speed_mu = 25, speed_sd = 0))  
    d2 <- cort_sim2()
    d3 <- cort_sim2(cort_sim1(max_sd = 0))
    p1 <- plot_cort_sim_1(d1)
    p2 <- plot_cort_sim_1(d2)
    p3 <- plot_cort_sim_1(d3)
    ggpubr::ggarrange(p1, p2, p3, nrow = 1)
                    
## Simulate data ----
    #plot_cort_sim(cort_sim2(cort_sim1(n = 60, speed_mu = 18, speed_sd = 3, base_mu = 9, base_sd = 3, cor_base_max = 0.3, cor_base_speed = -0.5, cor_speed_max = -0.7, cor_max_return = 0.4), 
    #                        bleed_times = c(1, 15, 30), sample_times = 1))
          
    # Speed vs. scope illustration
          set.seed(80)
          d1 <- cort_sim2(cort_sim1(cor_speed_max = -0.8))
          d2 <- cort_sim2(cort_sim1(cor_speed_max = 0))
          d3 <- cort_sim2(cort_sim1(cor_speed_max = 0.8))
          p1 <- plot_cort_sim_1(data = d1)
          p2 <- plot_cort_sim_1(data = d2)
          p3 <- plot_cort_sim_1(data = d3)
          
          ggpubr::ggarrange(p1, p2, p3, nrow = 1)
          
    # Only variation in speed. None in max
          set.seed(80)
          d1 <- cort_sim2(cort_sim1(max_sd = 0, speed_sd = 0))
          d2 <- cort_sim2(cort_sim1(max_sd = 0, speed_sd = 6))
          d3 <- cort_sim2(cort_sim1(max_sd = 0, speed_sd = 12))
          p1 <- plot_cort_sim_1(data = d1)
          p2 <- plot_cort_sim_1(data = d2)
          p3 <- plot_cort_sim_1(data = d3)
          
          ggpubr::ggarrange(p1, p2, p3, nrow = 1)
          
    