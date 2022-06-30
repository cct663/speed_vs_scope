### Script for speed vs. scope manuscript
# Written by Conor Taff; cct63@cornell.edu; cct663@gmail.com

## Load packages ----
    pacman::p_load(ggplot2, here, MASS, MBESS, tidyverse, viridis, gridExtra, rptR, simcoRt, sp)

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
     
          
     
     
     