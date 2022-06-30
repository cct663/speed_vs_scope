### Script for speed vs. scope manuscript
# Written by Conor Taff; cct63@cornell.edu; cct663@gmail.com


## Load packages ----
    pacman::p_load(ggplot2, here, tidyverse, viridis, gridExtra)

## Load RWBL field data ----
        rw_dat <- read.delim(here::here("1_raw_data/rwbl_data.txt"))      
        rw_dat <- subset(rw_dat, rw_dat$Bird_no != 10)    # This bird is missing most data
        rw_dat <- subset(rw_dat, rw_dat$band != "1372-33452")
    
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
        
        # Calculate max from loess regression
          for(i in 1:nrow(rw_dat)){	
            yy <- as.vector(t(rw_dat[i, 38:45]))	
            xx <- as.vector(t(rw_dat[i, 30:37]))
            dd <- as.data.frame(cbind(yy, xx))
            colnames(dd) <- c("yy", "xx")
            ss <- loess(yy ~ xx, span = 1.5) 
            max <- max(predict(ss, seq(min(na.omit(dd$xx)), max(na.omit(dd$xx)))))
            pct <- predict(ss, seq(min(na.omit(dd$xx)), max(na.omit(dd$xx)))) / max
            se <- seq(min(na.omit(dd$xx)), max(na.omit(dd$xx)))
            rw_dat$thirtypct[i] <- pct[which(1800 == se)]
            rw_dat$time95[i] <- se[min(which(pct > 0.95))]
            rw_dat$time.max[i] <- se[which.max(se)]
            pred <- predict(ss,seq(min(na.omit(dd$xx)), max(na.omit(dd$xx))))
            rw_dat$max.loess[i] <- max(pred)
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
        
## Time95 vs. max RWBL ----
        p1 <- ggplot(data = rw_dat, mapping = aes(x = time95 / 60, y = max.loess, color = LH.stage, fill = LH.stage)) +
          geom_point() + geom_smooth(method = "lm") +
          scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
          theme_classic() + xlab("Time to 95% of max") + ylab("Max cort (ng/ul)") +
          facet_wrap(~ LH.stage, scale = "free") +
          guides(fill = FALSE, color = FALSE)
        ggsave(here::here("3_r_scripts/time95_vs_max.png"),
               p1, device = "png", width = 8.4, height = 3.1, units = "in")
        
        p2 <- ggplot(data = rw_dat, mapping = aes(x = time95 / 60, y = max.l, color = LH.stage, fill = LH.stage)) +
          geom_point() + geom_smooth(method = "lm") +
          scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
          theme_classic() + xlab("Time reaching 95% of maximal") + ylab("log max cort (ng/ul)") +
          facet_wrap(~ LH.stage, scale = "free") +
          guides(fill = FALSE, color = FALSE)
        ggsave(here::here("3_r_scripts/time95_vs_maxl.png"),
               p2, device = "png", width = 8.4, height = 3.1, units = "in")
                

    