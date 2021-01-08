### Script for speed vs. scope manuscript
# Written by Conor Taff; cct63@cornell.edu; cct663@gmail.com


## Load packages ----
    pacman::p_load(ggplot2, here, MASS, MBESS, tidyverse, viridis, gridExtra)

## Load RWBL field data ----
        rw_dat <- read.delim(here::here("1_raw_data/rwbl_data.txt"))      
        rw_dat <- subset(rw_dat, rw_dat$Bird_no != 10)    # This bird is missing most data
    
    ## Replacing labeling on a couple rows to make them conform to level names
        rw_dat$LH.stage<-gsub("mid-breeding","early-breeding",rw_dat$LH.stage)	
        rw_dat$LH.stage[38]<-"early-breeding"
    
    ## Make a list of colors in order of rows to plot LH stages in different colors
        collist<-rep(NA,nrow(rw_dat))
        for(i in 1:nrow(rw_dat)){
          if(rw_dat$LH.stage[i]=="early-breeding"){collist[i]<-"lightblue"}
          if(rw_dat$LH.stage[i]=="late-breeding"){collist[i]<-"gray60"}
          if(rw_dat$LH.stage[i]=="molt"){collist[i]<-"red"}
        }

## Analysis of RWBL field data ----
        ## Loop through data frame to calculate slope and max cort achieved for each bird's samples
            for(i in 1:nrow(rw_dat)){
              # Set up dataframe for each bird with just the first three samples (less than 10 minutes)
                calc <- as.data.frame(matrix(nrow = 3, ncol = 3))
                colnames(calc) <- c("time", "cort", "lcort")
                calc[, 1] <- t(rw_dat[i, c("tube1_tsec", "tube2_tsec", "tube3_tsec")])
                calc[, 2] <- t(rw_dat[i, c("tube1_cort", "tube2_cort", "tube3_cort")])
                calc[, 3] <- log(calc[, 2] + 1)
                rownames(calc) <- c(1, 2, 3)
                calc <- subset(calc, calc$time < 600 & is.na(calc$cort) == FALSE)
              # Calculate the slope over the first 8 minutes from a simple linear regression. Both normal and log scale.   
                calc$time2 <- calc$time / 60
                ifelse(nrow(calc) > 1, slope <- coef(lm(calc$cort ~ calc$time2))[2], slope <- NA)
                rw_dat$slope[i] <- slope
                ifelse(nrow(calc) > 1, rw_dat$slope.l[i] <- coef(lm(calc$lcort ~ calc$time2))[2], rw_dat$slope.l[i] <- NA)
              # Calculate max cort reached on log and normal scale from loess regression
                values <- t(rw_dat[i, c("tube1_cort", "tube2_cort", "tube3_cort", "tube4_cort", "tube5_cort",
                                  "tube6_cort", "tube7_cort", "tube8_cort")])
                values.l <- log(values + 1)
                rw_dat$max[i] <- max(na.omit(values))
                rw_dat$max.l[i] <- max(na.omit(values.l))
                values.l <- gsub("NA", 0, values.l)
                rw_dat$max.t[i] <- which.max(values.l)
            } 
        
        ### Make a long version of the data for easier plotting
              rwd_long1 <- as.data.frame(pivot_longer(rw_dat, cols = tube1_tsec:tube8_cort,
                                       names_to = "label", values_to = "values"))
              for(i in 1:nrow(rwd_long1)){
                rwd_long1$tube[i] <- strsplit(rwd_long1$label[i], "_")[[1]][1]
                rwd_long1$type[i] <- strsplit(rwd_long1$label[i], "_")[[1]][2]
              }
              rwd_long <- as.data.frame(pivot_wider(rwd_long1, id_cols = c(tube, band, type, Cap_Time_s, sex, 
                                                            mass, tarsus, wing, head.bill, LH.stage,
                                                             fat, age, slope, slope.l, max, max.l, max.t, thirtypct, time95,
                                                             time.max, max.loess), 
                                      values_from = values, names_from = type))

## Descriptive RWBL plots ----                            
        # Plot general pattern of cort response for each group    
              pa <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = cort, by = band)) + 
                geom_point(alpha = 0.6) + 
                geom_smooth(method = "loess", span = 1.2, se = FALSE, color = "gray60", alpha = 0.7, size = 0.6) + 
                facet_wrap(~ LH.stage, scale = "free_y") + theme_classic() + xlab("Time (minutes)") +
                ylab("Corticosterone (ng/ul)") +
                geom_smooth(data = rwd_long, mapping = aes(x = tsec / 60, y = cort, by = LH.stage), col = "coral3", fill = "coral3")
              ggsave(here::here("3_r_scripts/basic_rwbl_profile.png"),
                     pa, device = "png", width = 9.7, height = 3.75, units = "in")
              
        # Same plot as above but on log scale
              pal <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = log(cort), by = band)) + 
                geom_point(alpha = 0.6) + 
                geom_smooth(method = "loess", span = 1.2, se = FALSE, color = "gray60", alpha = 0.7, size = 0.6) + 
                facet_wrap(~ LH.stage, scale = "free_y") + theme_classic() + xlab("Time (minutes)") +
                ylab("Corticosterone (ng/ul)") +
                geom_smooth(data = rwd_long, mapping = aes(x = tsec / 60, y = log(cort), by = LH.stage), col = "coral3", fill = "coral3")
              ggsave(here::here("3_r_scripts/basic_rwbl_profile_log.png"),
                     pal, device = "png", width = 9.7, height = 3.75, units = "in")
              
        # On one plot normal scale
              pb <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = cort, color = LH.stage, by = band)) +
                geom_point(alpha = 0.6) +
                #geom_smooth(method = "loess", span = 1.2, se = FALSE, alpha = 0.7, size = 0.6) +
                #stat_smooth(geom = "line", method = "loess", span = 1.5, alpha = 0.6, size = 0.7, se = FALSE) +
                geom_line(alpha = 0.6, size = 0.7) +
                scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) +
                theme_classic() + xlab("Time (minutes)") + ylab("Corticosterone (ng/ul)") +
                guides(fill = FALSE, color = FALSE) + coord_cartesian(y = c(0, 120))
              
        # On one plot log scale
              pb2 <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = log(cort), color = LH.stage, by = band)) +
                geom_point(alpha = 0.6) +
                geom_line(alpha = 0.6, size = 0.7) +
                #geom_smooth(method = "loess", span = 1.2, se = FALSE, alpha = 0.7, size = 0.6) +
                #stat_smooth(geom = "line", method = "loess", span = 1.5, alpha = 0.6, size = 0.7, se = FALSE) +
                scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) +
                theme_classic() + xlab("Time (minutes)") + ylab("log Corticosterone (ng/ul)") +
                guides(fill = FALSE, color = FALSE)
          
          ggsave(here::here("3_r_scripts/summary_lines.png"),        
                grid.arrange(pb, pb2, nrow = 1),
                device = "png", width = 9.7, height = 4.8)

## Shape of response by stage ----        
        # Plot shape by life history stage as percent of maximum 
              
              rwd_long$pct_max <- rwd_long$cort / rwd_long$max
              p1 <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = pct_max, color = LH.stage, fill = LH.stage)) + 
                geom_smooth(method = "loess", span = 1.2, se = TRUE) + 
                scale_fill_viridis(discrete = TRUE) +
                scale_color_viridis(discrete = TRUE) + 
                xlab("Time (minutes)") +
                ylab("Corticosterone percent of maximum") + 
                theme_classic() + ylim(0, 1.12) + 
                coord_cartesian(xlim = c(0, 60)) +
                theme(legend.position = c(0.8, 0.2)) +
                guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)
              
              p2 <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = pct_max, color = LH.stage, fill = LH.stage, by = band)) + 
                stat_smooth(geom = "line", method = "loess", span = 1.5, alpha = 0.6, size = 0.7, se = FALSE) + 
                scale_fill_viridis(discrete = TRUE) +
                scale_color_viridis(discrete = TRUE) + 
                xlab("Time (minutes)") +
                ylab("Corticosterone percent of maximum") + 
                theme_classic() + ylim(0, 1.12) +
                coord_cartesian(xlim = c(0, 60)) +
                theme(legend.position = c(0.8, 0.2)) +
                guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)
              
              ggsave(here::here("3_r_scripts/rwbl_s_profiles.png"),
                  ggpubr::ggarrange(p2, p1),
                  device = "png", width = 9, height = 4, units = "in")
              
              
          # Plot shape by life history stage as percent of from log scale 
              
              rwd_long$pct_maxl <- log(rwd_long$cort + 1) / log(exp(rwd_long$max.l) + 1)
              p1 <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = pct_maxl, color = LH.stage, fill = LH.stage)) + 
                geom_smooth(method = "loess", span = 1.2, se = TRUE) + 
                scale_fill_viridis(discrete = TRUE) +
                scale_color_viridis(discrete = TRUE) + 
                xlab("Time (minutes)") +
                ylab("Corticosterone percent of maximum") + 
                theme_classic() + ylim(0, 1.12) + 
                coord_cartesian(xlim = c(0, 60)) +
                theme(legend.position = c(0.8, 0.2)) +
                guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)
              
              p2 <- ggplot(data = rwd_long, mapping = aes(x = tsec / 60, y = pct_maxl, color = LH.stage, fill = LH.stage, by = band)) + 
                stat_smooth(geom = "line", method = "loess", span = 1.5, alpha = 0.6, size = 0.7, se = FALSE) + 
                scale_fill_viridis(discrete = TRUE) +
                scale_color_viridis(discrete = TRUE) + 
                xlab("Time (minutes)") +
                ylab("Corticosterone percent of log maximum") + 
                theme_classic() + ylim(0, 1.12) +
                coord_cartesian(xlim = c(0, 60)) +
                theme(legend.position = c(0.8, 0.2)) +
                guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)
              
              ggsave(here::here("3_r_scripts/rwbl_s_profiles_log.png"),
                     ggpubr::ggarrange(p2, p1),
                     device = "png", width = 9, height = 4, units = "in")
            
## Amount of variation ----
        vari <- rw_dat[, c("LH.stage", "band", "tube1_cort", "max", "max.l", "slope", "slope.l", "time95")]      
        colnames(vari)[3] <- "base"
        vari$base.l <- log(vari$base + 1)
        s_vari <- vari %>%
          group_by(LH.stage) %>%
          summarize(base = sd(na.omit(base)) / mean(na.omit(base)) * 100,
                    log_base = sd(na.omit(base.l)) / mean(na.omit(base.l)) * 100,
                    max = sd(na.omit(max)) / mean(na.omit(max)) * 100,
                    log_max = sd(na.omit(max.l)) / mean(na.omit(max.l)) * 100,
                    slope = sd(na.omit(slope)) / mean(na.omit(slope)) * 100,
                    log_slope = sd(na.omit(slope.l)) / mean(na.omit(slope.l)) * 100,
                    time95 = sd(na.omit(time95)) / mean(na.omit(time95)) * 100)
        saveRDS(s_vari, here::here("3_r_scripts/cv_table.rds"))
              
## Slope vs. max RWBL ----
     # Plot on normal scale
          r_all <- round(summary(lm(max ~ slope, data = rw_dat))$r.squared, 2)
          p_all <- round(summary(lm(max ~ slope, data = rw_dat))$coefficients[2, 4], 2)
          p1 <- ggplot(data = rw_dat, mapping = aes(x = slope, y = max)) +
                geom_point(mapping = aes(color = LH.stage)) + 
                geom_smooth(method = "lm", color = "gray40", alpha = 0.2) +
                scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
                theme_classic() + xlab("Slope of initial increase (ng/ul/min)") +
                ylab("Maximal corticosterone achieved (ng/ul)") +
                #coord_cartesian(ylim = c(-10, 200)) +
                theme(legend.position = c(0.2, 0.8)) +
                guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)+
            annotate(geom = "text", label = paste0("R2 = ", r_all), x = 5, y = 7) +
            annotate(geom = "text", label = paste0("P = ", p_all), x = 4.9, y = -5)
              
          cols <- viridis(n = 3)
          
          r_eb <- round(summary(lm(max ~ slope, data = subset(rw_dat, rw_dat$LH.stage == "early-breeding")))$r.squared, 2)    
          p_eb <- round(summary(lm(max ~ slope,data = subset(rw_dat, rw_dat$LH.stage == "early-breeding")))$coefficients[2, 4], 2)
          p2 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "early-breeding"), mapping = aes(x = slope, y = max)) +
            geom_point(color = cols[1]) + 
            geom_smooth(method = "lm", color = cols[1], alpha = 0.2, fill = cols[1]) +
            theme_classic() + xlab("") +
            ylab("") +
            guides(fill = FALSE, color = FALSE) +
            annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
            annotate(geom = "text", label = paste0("R2 = ", r_eb), x = 5, y = 28) +
            annotate(geom = "text", label = paste0("P = ", p_eb), x = 4.9, y = -15)
          
          r_lb <- round(summary(lm(max ~ slope, data = subset(rw_dat, rw_dat$LH.stage == "late-breeding")))$r.squared, 2)    
          p_lb <- round(summary(lm(max ~ slope,data = subset(rw_dat, rw_dat$LH.stage == "late-breeding")))$coefficients[2, 4], 2)
          p3 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "late-breeding"), mapping = aes(x = slope, y = max)) +
            geom_point(color = cols[2]) + 
            geom_smooth(method = "lm", color = cols[2], alpha = 0.2, fill = cols[2]) +
            theme_classic() + xlab("") +
            ylab("") +
            guides(fill = FALSE, color = FALSE) +
            annotate(geom = "text", label = "C", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
            annotate(geom = "text", label = paste0("R2 = ", r_lb), x = 0.8, y = 7) +
            annotate(geom = "text", label = paste0("P = ", p_lb), x = 0.8, y = 4)
          
          r_m <- round(summary(lm(max ~ slope, data = subset(rw_dat, rw_dat$LH.stage == "molt")))$r.squared, 2)    
          p_m <- round(summary(lm(max ~ slope, data = subset(rw_dat, rw_dat$LH.stage == "molt")))$coefficients[2, 4], 2)
          p4 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "molt"), mapping = aes(x = slope, y = max)) +
            geom_point(color = cols[3]) + 
            geom_smooth(method = "lm", color = cols[3], alpha = 0.2, fill = cols[3]) +
            theme_classic() + xlab("") +
            ylab("") +
            guides(fill = FALSE, color = FALSE) +
            annotate(geom = "text", label = "D", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
            annotate(geom = "text", label = paste0("R2 = ", r_m), x = 0.38, y = 6) +
            annotate(geom = "text", label = paste0("P = ", p_m), x = 0.375, y = 3)
          
          ggsave(here::here("3_r_scripts/slope_vs_max.png"), 
                 grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3), c(1, 1, 4))),
                 device = "png", width = 9, height = 6.2, units = "in")
     
              
          # Plot on log scale
                r_all <- round(summary(lm(max.l ~ slope.l, data = rw_dat))$r.squared, 2)
                p_all <- round(summary(lm(max.l ~ slope.l, data = rw_dat))$coefficients[2, 4], 2)
                p1 <- ggplot(data = rw_dat, mapping = aes(x = slope.l, y = max.l)) +
                  geom_point(mapping = aes(color = LH.stage)) + 
                  geom_smooth(method = "lm", color = "gray40", alpha = 0.2) +
                  scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
                  theme_classic() + xlab("Slope of initial increase log scale (ng/ul/min)") +
                  ylab("Maximal corticosterone achieved (log ng/ul)") +
                  #coord_cartesian(ylim = c(-10, 200)) +
                  theme(legend.position = c(0.2, 0.8)) +
                  guides(fill = guide_legend(title = "stage"), color = guide_legend(title = "stage")) +
                  annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5)+
                  annotate(geom = "text", label = paste0("R2 = ", r_all), x = 0.35, y = 1.8) +
                  annotate(geom = "text", label = paste0("P = ", p_all), x = 0.35, y = 1.2)
                
                cols <- viridis(n = 3)
                
                r_eb <- round(summary(lm(max.l ~ slope.l, data = subset(rw_dat, rw_dat$LH.stage == "early-breeding")))$r.squared, 2)    
                p_eb <- round(summary(lm(max.l ~ slope.l,data = subset(rw_dat, rw_dat$LH.stage == "early-breeding")))$coefficients[2, 4], 2)
                p2 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "early-breeding"), mapping = aes(x = slope.l, y = max.l)) +
                  geom_point(color = cols[1]) + 
                  geom_smooth(method = "lm", color = cols[1], alpha = 0.2, fill = cols[1]) +
                  theme_classic() + xlab("") +
                  ylab("") +
                  guides(fill = FALSE, color = FALSE) +
                  annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
                  annotate(geom = "text", label = paste0("R2 = ", r_eb), x = .35, y = 3.5) +
                  annotate(geom = "text", label = paste0("P = ", p_eb), x = .35, y = 3)
                
                r_lb <- round(summary(lm(max.l ~ slope.l, data = subset(rw_dat, rw_dat$LH.stage == "late-breeding")))$r.squared, 2)    
                p_lb <- round(summary(lm(max.l ~ slope.l,data = subset(rw_dat, rw_dat$LH.stage == "late-breeding")))$coefficients[2, 4], 2)
                p3 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "late-breeding"), mapping = aes(x = slope.l, y = max.l)) +
                  geom_point(color = cols[2]) + 
                  geom_smooth(method = "lm", color = cols[2], alpha = 0.2, fill = cols[2]) +
                  theme_classic() + xlab("") +
                  ylab("") +
                  guides(fill = FALSE, color = FALSE) +
                  annotate(geom = "text", label = "C", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
                  annotate(geom = "text", label = paste0("R2 = ", r_lb), x = 0.22, y = 1.7) +
                  annotate(geom = "text", label = paste0("P = ", p_lb), x = 0.22, y = 1.5)
                
                r_m <- round(summary(lm(max.l ~ slope.l, data = subset(rw_dat, rw_dat$LH.stage == "molt")))$r.squared, 2)    
                p_m <- round(summary(lm(max.l ~ slope.l, data = subset(rw_dat, rw_dat$LH.stage == "molt")))$coefficients[2, 4], 2)
                p4 <- ggplot(data = subset(rw_dat, rw_dat$LH.stage == "molt"), mapping = aes(x = slope.l, y = max.l)) +
                  geom_point(color = cols[3]) + 
                  geom_smooth(method = "lm", color = cols[3], alpha = 0.2, fill = cols[3]) +
                  theme_classic() + xlab("") +
                  ylab("") +
                  guides(fill = FALSE, color = FALSE) +
                  annotate(geom = "text", label = "D", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5) +
                  annotate(geom = "text", label = paste0("R2 = ", r_m), x = 0.18, y = 2.1) +
                  annotate(geom = "text", label = paste0("P = ", p_m), x = 0.18, y = 1.9)
                
                ggsave(here::here("3_r_scripts/slope_vs_max_log.png"), 
                       grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3), c(1, 1, 4))),
                       device = "png", width = 9, height = 6.2, units = "in")       
      
## Slope vs. time95 RWBL ----
        p1 <- ggplot(data = rw_dat, mapping = aes(x = slope, y = time95 / 60, color = LH.stage, fill = LH.stage)) +
            geom_point() + geom_smooth(method = "lm") +
            scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
            theme_classic() + xlab("Initial slope") + ylab("Time reaching 95% of maximal") +
            facet_wrap(~ LH.stage, scale = "free") +
                      guides(fill = FALSE, color = FALSE)
        ggsave(here::here("3_r_scripts/slope_vs_95.png"),
               p1, device = "png", width = 8.4, height = 3.1, units = "in")
                
        p2 <- ggplot(data = rw_dat, mapping = aes(x = slope.l, y = time95 / 60, color = LH.stage, fill = LH.stage)) +
          geom_point() + geom_smooth(method = "lm") +
          scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
          theme_classic() + xlab("Initial slope (log scale)") + ylab("Time reaching 95% of maximal") +
          facet_wrap(~ LH.stage, scale = "free") +
          guides(fill = FALSE, color = FALSE)
        ggsave(here::here("3_r_scripts/slope_vs_95_log.png"),
               p2, device = "png", width = 8.4, height = 3.1, units = "in")
                
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
          cort_sim1 <- function(n = 25, 
                               base_mu = 5, base_sd = 1.5, base_min = 0.05,
                               speed_mu = 40, speed_sd = 10,
                               max_mu = log(50), max_sd = log(1.2), 
                               maxtime_mu = 4, maxtime_sd = 1,
                               return_mu = 120, return_sd = 15, 
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
                      # replace low base values with detection limit
                          sim_dat[which(sim_dat$base < base_min), ] <- base_min
                  
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
                                bleed_times = c(2, 30),
                                base_error = 0.1,
                                speed_error = 0.1,
                                max_error = 0.1,
                                maxtime_error = 0.1,
                                return_error = 0.1,
                                sample_times = 2,
                                timecourse_max = 170,
                                performance_contributions = c(0, 40, 40, 0, 0, 20)) #relative contribution of base, speed, max, maxtime, return, and error to performance
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
                fits <- predict(loess(temp$y ~ temp$x, span = 0.6), newdata = time_seq$time)
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
                              rank_timecourse = long_rank)
              return(dataset)
      }
      
    # Plot function. The datasets generated in cort_sim2 can be used to make any number of different plots or models. But this
        # function creates a three panel summary plot that illustrates the data. It's only input is a list with the three data
        # frames output by 'cort_sim2'. It will by default generate a new dataset when called if one isn't provided.
          plot_cort_sim <- function(data = cort_sim2())
              {
            
            spoints <- unique(data$simulated_dataset_long$time)
            p1 <- ggplot(data = data$simulated_dataset, mapping = aes(x = time, y = cort, color = animal, by = animal_sample)) +
              geom_line(size = 0.5) + guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") 
            p2 <- ggplot(data = data$timecourse_long, mapping = aes(x = time, y = cort, color = animal, by = animal_sample)) +
              geom_smooth(method = "loess", span = 0.3, size = 0.5, se = FALSE) + guides(color = FALSE) +
              theme_classic() + xlab("Time") + ylab("Corticosterone") + coord_cartesian(xlim = c(0, 60)) +
              geom_vline(xintercept = spoints, linetype = "dashed")
            p3 <- ggplot(data = data$rank_timecourse, mapping = aes(x = time, y = rank, color = animal_sample2)) +
              geom_line(size = 1.6, alpha = 0.8) + guides(color = FALSE) + theme_classic() + 
              scale_color_viridis(discrete = TRUE) + xlim(0, 60) +
              geom_vline(xintercept = spoints, linetype = "dashed")
            grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 2, 2), c(3, 3, 3)))
          }
      
## Simulate data ----
    plot_cort_sim(cort_sim2(cort_sim1(n = 60, speed_mu = 18, speed_sd = 3, base_mu = 9, base_sd = 3, cor_base_max = 0.3, cor_base_speed = -0.5, cor_speed_max = -0.7, cor_max_return = 0.4), 
                            bleed_times = c(1, 15, 30), sample_times = 1))
    