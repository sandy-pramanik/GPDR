


# Fig 1: KDE vs GPDR vs BDR ----
rm(list = ls())
library(tidyverse)

## KDE and GP fit
simulation.output.path = "/Users/sandipanpramanik/mystaff/JHU/research/Distribution regression/code/github/3-31-24/simulation_output"    # specify path to the GPDR_sourcecode folder
load(file.path(simulation.output.path, "simulation_results_alpha_25.RData"))
xx = as.data.frame(simRes)
# xx = read.table("simulation_results", header=TRUE)

ss = xx %>% group_by(V1, V2) %>% summarise(
  exp_l2_mean = mean(V4) / 100,
  exp_l2_std = sqrt(var(V4) / n()) / 100,
  kde_l2_mean = mean(V5) / 100,
  kde_l2_std = sqrt(var(V5) / n()) / 100
)

## BDR fit
load(file.path(simulation.output.path, "BDR_simulation_results_alpha_25.RData"))
bb = as.data.frame(BDR_simRes)
# bb = read.table("BDR/bdr_simulation_results", header=TRUE)

ss2 = bb %>% group_by(V1, V2) %>% summarise(
  bdr_k10_l2_mean = mean(V4) / 100,
  bdr_k10_l2_std = sqrt(var(V4) / n()) / 100,
  bdr_k50_l2_mean = mean(V5) / 100,
  bdr_k50_l2_std = sqrt(var(V5) / n()) / 100
)


## plot
data.frame(
  n = c(ss$V1, ss$V1, ss2$V1, ss2$V1),
  m = c(ss$V2, ss$V2, ss2$V2, ss2$V2),
  method = c(rep("GPDR", nrow(ss)), rep("KDE", nrow(ss)), 
             rep("BDR\n (k=10)", nrow(ss2)), rep("BDR\n (k=50)", nrow(ss2))),
  L2 = c(ss$exp_l2_mean, ss$kde_l2_mean, ss2$bdr_k10_l2_mean, ss2$bdr_k50_l2_mean), 
  std = c(ss$exp_l2_std, ss$kde_l2_std, ss2$bdr_k10_l2_std, ss2$bdr_k50_l2_std)
) %>% ggplot(aes(x=m, y=L2, color=method)) +
  geom_line() + 
  geom_errorbar(aes(ymin = L2 - 1.96*std, ymax = L2 + 1.96*std),
                position = "identity", width = 0.05) +
  facet_wrap(vars(n), nrow=2, ncol=3, labeller=label_both) +
  scale_x_log10() + 
  scale_y_log10() +
  ylab("square L2") +
  theme_bw()

data.frame(
  n = c(ss$V1, ss$V1, ss2$V1, ss2$V1),
  m = c(ss$V2, ss$V2, ss2$V2, ss2$V2),
  method = c(rep("GPDR", nrow(ss)), rep("KDE", nrow(ss)), 
             rep("BDR\n(k=10)\n", nrow(ss2)), rep("BDR\n(k=50)\n", nrow(ss2))),
  L2 = c(ss$exp_l2_mean, ss$kde_l2_mean, ss2$bdr_k10_l2_mean, ss2$bdr_k50_l2_mean), 
  std = c(ss$exp_l2_std, ss$kde_l2_std, ss2$bdr_k10_l2_std, ss2$bdr_k50_l2_std)
) %>% dplyr::filter(m > 10) %>% ggplot(aes(x=n, y=L2, color=method)) +
  geom_line(linewidth=0.75) + geom_point(aes(shape=method)) + 
  # geom_errorbar(aes(ymin = L2 - 1.96*std, ymax = L2 + 1.96*std),
  #               position = "identity", width = 0.05) +
  facet_wrap(vars(m), nrow=2, ncol=3, labeller=label_both) +
  scale_x_log10() + 
  scale_y_log10() +
  # scale_x_continuous(limits=c(50, 100, 200, 300, 400)) +
  ylab("Posterior Risk") +
  # ggplot2::theme(
  #   axis.title.x = ggplot2::element_text(size=15,
  #                                        margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
  #   axis.title.y = ggplot2::element_text(size=15,
  #                                        margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
  #   axis.text.x = ggplot2::element_text(color = "black", size = 13,
  #                                       angle = 0, hjust = .5, vjust = 1),
  #   axis.text.y = ggplot2::element_text(color = "black", size = 13),
  #   axis.ticks.x = ggplot2::element_line(linewidth = .5),
  #   axis.ticks.length.x = ggplot2::unit(.2, "cm"),
  #   axis.ticks.y = ggplot2::element_line(linewidth = .5),
#   axis.ticks.length.y = ggplot2::unit(.2, "cm"),
#   panel.background = ggplot2::element_blank(),
#   panel.border = ggplot2::element_rect(color='black', linetype = "solid",
#                                        fill = NA, linewidth = 1),
#   panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
#                                            colour = "grey90"),
#   panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
#                                            colour = "grey90"),
#   strip.text.x = ggplot2::element_text(
#     size = 15,
#     face = "bold"
#   ),
#   strip.text.y = ggplot2::element_text(
#     size = 15,
#     face = "bold"
#   ),
#   strip.background = ggplot2::element_rect(color="black", linewidth=1),
#   # legend.title = ggplot2::element_blank(),
#   legend.key.width = ggplot2::unit(1, "cm"),
#   legend.key.height = ggplot2::unit(.75, "cm"),
#   # legend.key.size = ggplot2::unit(.5, "cm"),
#   legend.spacing.x = ggplot2::unit(.5, 'cm'),
#   legend.text=ggplot2::element_text(size=13),
#   legend.position = 'bottom'
# )
theme_bw()  # 8 X 5
# theme(legend.position = c(0.9, 0.4), 
#       legend.text = element_text(size=5))

# ggsave(filename="plots/risk_by_m.png", p2, width=8, height=6)



# Fig 2: empirical cdf of DP(.01) vs DP(25) ----
rm(list = ls())

# sourcecode.path = xxx    # specify path to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_functions.R"))

nSamp = 100
set.seed(8)
dpi.df = rbind.data.frame(data.frame('dp_samples' = rdp(n = nSamp, alpha = .1)$samples,
                                     'alpha_dp' = 0.1),
                          data.frame('dp_samples' = rdp(n = nSamp, alpha = 25)$samples,
                                     'alpha_dp' = 25))
head(dpi.df)

dpi.df$alpha_dp = factor(x = dpi.df$alpha_dp,
                         labels = c(paste(expression(alpha*' = 0.1')),
                                    paste(expression(alpha*' = 25'))),
                         levels = c(0.1, 25))

alpha_dp.label = c('0.1' = paste(expression(alpha*' = 0.1')),
                   '25' = paste(expression(alpha*' = 25')))
ggplot2::ggplot(data = dpi.df) +
  ggplot2::coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  ggplot2::facet_grid(.~alpha_dp,
                      labeller = label_parsed
                      # labeller = ggplot2::label_bquote(cols = alpha .(alpha_dp))
                      # labeller = ggplot2::labeller(alpha_dp = alpha_dp.label)
  ) +
  ggplot2::stat_ecdf(ggplot2::aes(x = dp_samples),
                     geom = "step", linewidth = 1) +
  # ggplot2::geom_line(
  #   ggplot2::aes(x = roc_fp, y = roc_sens),
  #   linewidth = .7, color = 4
  # ) +
  # ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  # ggplot2::geom_label(
  #   data = plotdf.auc,
  #   ggplot2::aes(x = 75, y = 25, label = auc_label), fontface = "bold"
  # ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size=22, face="bold"),
                 plot.subtitle = ggplot2::element_text(size=20),
                 axis.title.x = ggplot2::element_text(size=20,
                                                      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                 # axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size=20,
                                                      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                 axis.text.x = ggplot2::element_text(color = "black", size = 15,
                                                     angle = 0, hjust = .5, vjust = 1),
                 # axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(color = "black", size = 15),
                 # axis.ticks.x = ggplot2::element_line(linewidth = 1),
                 axis.ticks.length.x = ggplot2::unit(.2, "cm"),
                 # axis.ticks.x = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                      fill = NA, linewidth = 1.5),
                 panel.background = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                          colour = "grey90"),
                 panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                          colour = "grey90"),
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold"),
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold"),
                 strip.background = ggplot2::element_rect(color="black", linewidth=1),
                 legend.title = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"), 
                 legend.key.size = ggplot2::unit(30, "cm"),
                 legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=15),
                 legend.position = 'bottom') +
  # ggplot2::guides(color = ggplot2::guide_legend(
  #   nrow = 2, byrow=FALSE#,
  #   # override.aes = list(
  #   #   linetype = "solid",
  #   #   shape = c(rep(16, length(levels(plotdf$method))-1), NA))
  # )) +
  ggplot2::labs(#title = 'International Trade Data Analysis', subtitle = 'Trade Occurence Prediction',
    x = 'x',
    y = 'Empirical CDF')   # 11 X 5.5



# Fig 4: risks for alpha=0.01: KDE vs GPDR ----
rm(list = ls())
library(tidyverse)

simulation.output.path = "/Users/sandipanpramanik/mystaff/JHU/research/Distribution regression/code/github/3-31-24/simulation_output"    # specify path to the GPDR_sourcecode folder
load(file.path(simulation.output.path, "simulation_results_alpha_0.1.RData"))
xx = as.data.frame(simRes)
head(xx)
class(xx)

ss = xx %>% group_by(V1, V2) %>% summarise(
  exp_l2_mean = mean(V4) / 100,
  exp_l2_std = sqrt(var(V4) / n()) / 100,
  kde_l2_mean = mean(V5) / 100,
  kde_l2_std = sqrt(var(V5) / n()) / 100
)
head(ss)

postrisk.df = rbind.data.frame(data.frame(n=ss$V1, m=ss$V2, method="GPDR", 
                                          L2=ss$exp_l2_mean, std=ss$exp_l2_std),
                               data.frame(n=ss$V1, m=ss$V2, method="KDE", 
                                          L2=ss$kde_l2_mean, std=ss$kde_l2_std))
postrisk.df$m = as.factor(postrisk.df$m)
postrisk.df$method = factor(x=postrisk.df$method,
                            levels = c('GPDR', 'KDE'))
postrisk.df=postrisk.df[postrisk.df$m %in% c(10, 100, 500, 1000),]

method.color = scales::hue_pal()(4)[3:4]
method.pointshape = c(15, 3)
names(method.color) = names(method.pointshape) = c('GPDR', 'KDE')

ggplot(data = postrisk.df) +
  facet_wrap(vars(m), nrow=1, ncol=7, labeller=label_both) +
  geom_line(aes(x=n, y=L2, color=method),
            linewidth=1) +
  geom_point(aes(x=n, y=L2, shape=method, color=method),
             size = 3) + 
  scale_color_manual(values = method.color) +
  scale_shape_manual(values = method.pointshape) +
  scale_x_log10() + 
  scale_y_log10() +
  # ylab("square L2") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size=15,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=15,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 13,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 13),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 15,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 15,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(1, "cm"),
    legend.key.height = ggplot2::unit(.75, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.spacing.x = ggplot2::unit(.5, 'cm'),
    legend.text=ggplot2::element_text(size=13),
    legend.position = 'bottom'
  ) +
  # theme_bw() +
  labs(#title = 'alpha=0.01',
       x = 'n', y = 'Posterior Risk')   # 4 X 11



# Fig 6: risk ratios: KDE vs GPDR ----
rm(list = ls())
library(tidyverse)

simulation.output.path = "/Users/sandipanpramanik/mystaff/JHU/research/Distribution regression/code/github/3-31-24/simulation_output"    # specify path to the GPDR_sourcecode folder

sim.results = c('simulation_results_alpha_0.1',
                'simulation_results_alpha_25',
                'simulation_results_cp',
                'simulation_results_cp_dep')
sim_setting_name = c('DP (0.1)',
                     'DP (25)',
                     'Independent\nContinuous Covariates',
                     'Dependent\nContinuous Covariates')

postriskratio.df = NULL
for(l in 1:length(sim.results)){
  
  load(file.path(simulation.output.path, 
                 paste0(sim.results[l],".RData")))
  xx = as.data.frame(simRes)
  head(xx)
  class(xx)
  
  ss = xx %>% group_by(V1, V2) %>% summarise(
    exp_l2_mean = mean(V4) / 100,
    exp_l2_std = sqrt(var(V4) / n()) / 100,
    kde_l2_mean = mean(V5) / 100,
    kde_l2_std = sqrt(var(V5) / n()) / 100
  )
  head(ss)
  
  postriskratio.df = rbind.data.frame(postriskratio.df,
                                      data.frame(
                                        n = ss$V1,
                                        m = ss$V2,
                                        L2_ratio = ss$exp_l2_mean/ss$kde_l2_mean,
                                        sim_setting = sim_setting_name[l]
                                      ))
  
}
postriskratio.df$n = as.factor(postriskratio.df$n)
postriskratio.df$sim_setting = factor(x = postriskratio.df$sim_setting,
                                      levels = sim_setting_name)

ggplot(data = postriskratio.df,
       aes(x=m, y=L2_ratio, color=n)) +
  # facet_grid(.~sim_setting, scales="free_y") +
  facet_wrap(.~sim_setting, scales="free", nrow=1, ncol=4) +
  # facet_wrap(vars(sim_setting), nrow=1, ncol=4, labeller=label_both) +
  geom_line(linewidth=1.5) + 
  geom_point(shape=20, size=7) + 
  scale_x_log10() + 
  scale_y_log10() +
  # ylab("square L2") +
  # theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 25,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 25),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    # legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(2, "cm"),
    legend.key.height = ggplot2::unit(2, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.spacing.x = ggplot2::unit(.5, 'cm'),
    legend.text=ggplot2::element_text(size=30),
    legend.title = ggplot2::element_text(size=30),
    legend.position = 'bottom'
  ) +
  labs(#title = 'alpha=0.01',
    x = 'm', y = 'Ratio of Posterior Risk of\nGPDR With Respect To KDE')   # 25 X 7.5



# Fig 7: GPDR risks for varying degree of dependence ----
## indep vs rho=0 ====
rm(list = ls())
library(tidyverse)

simulation.output.path = "/Users/sandipanpramanik/mystaff/JHU/research/Distribution regression/code/4-2-24-1/simulation"    # specify path to the GPDR_sourcecode folder

sim.results = c('indep' = 'simulation_results_cp_indep',
                'rho0' = 'simulation_results_cp_dep_rho0')
sim_setting.label = c('Independent', 'AR(1) Coefficient 0')
names(sim_setting.label) = names(sim.results)

postrisk.df = NULL
for(l in 1:length(sim.results)){
  
  load(file.path(simulation.output.path, 
                 paste0(sim.results[l],".RData")))
  xx = as.data.frame(simRes)
  head(xx)
  class(xx)
  
  ss = xx %>% group_by(V1, V2) %>% summarise(
    exp_l2_mean = mean(V4) / 100,
    exp_l2_std = sqrt(var(V4) / n()) / 100,
    kde_l2_mean = mean(V5) / 100,
    kde_l2_std = sqrt(var(V5) / n()) / 100
  )
  head(ss)
  
  postrisk.df = rbind.data.frame(postrisk.df,
                                 data.frame(
                                   n = ss$V1,
                                   m = ss$V2,
                                   L2 = ss$exp_l2_mean,
                                   sim_setting = names(sim.results)[l]
                                 ))
  
}
postrisk.df$m = as.factor(postrisk.df$m)
postrisk.df$sim_setting = factor(x = postrisk.df$sim_setting,
                                 levels = names(sim.results))
# postrisk.df = postrisk.df[postrisk.df$m!=50,]

m.label = paste0('m: ', levels(postrisk.df$m))
names(m.label) = levels(postrisk.df$m)

sim_setting.color = scales::hue_pal()(length(sim.results))
sim_setting.pointshape = c(3, 4)
names(sim_setting.color) = names(sim_setting.pointshape) =
  names(sim.results)

ggplot(data = postrisk.df) +
  facet_grid(.~m, scales="fixed",
             labeller = ggplot2::labeller(m = m.label)) +
  # facet_wrap(.~m, scales="free", nrow=1, ncol=4, labeller=label_both) +
  # facet_wrap(vars(sim_setting), nrow=1, ncol=4, labeller=label_both) +
  geom_line(aes(x=n, y=L2, color=sim_setting), linewidth=1) + 
  geom_point(aes(x=n, y=L2, color=sim_setting, shape=sim_setting),
             size=5, stroke=1) + 
  scale_color_manual(values = sim_setting.color,
                     labels = sim_setting.label) +
  scale_shape_manual(values = sim_setting.pointshape,
                     labels = sim_setting.label) +
  scale_x_log10() + 
  scale_y_log10() +
  # ylab("square L2") +
  # theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 20,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 20),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(2, "cm"),
    legend.key.height = ggplot2::unit(2, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.spacing.x = ggplot2::unit(.5, 'cm'),
    legend.text=ggplot2::element_text(size=25),
    legend.position = 'bottom'
  ) +
  labs(#title = 'alpha=0.01',
    x = 'n', y = 'Posterior Risk')   # 13 X 6


## rho=0, rho=0.5, rho=0.8 ====
rm(list = ls())
library(tidyverse)

simulation.output.path = "/Users/sandipanpramanik/mystaff/JHU/research/Distribution regression/code/4-2-24-1/simulation"    # specify path to the GPDR_sourcecode folder

sim.results = c('rho0' = 'simulation_results_cp_dep_rho0',
                'rho0.5' = 'simulation_results_cp_dep_rho0.5',
                'rho0.8' = 'simulation_results_cp_dep_rho0.8')
sim_setting.label = c(0, 0.5, 0.8)
names(sim_setting.label) = names(sim.results)

postrisk.df = NULL
for(l in 1:length(sim.results)){
  
  load(file.path(simulation.output.path, 
                 paste0(sim.results[l],".RData")))
  xx = as.data.frame(simRes)
  head(xx)
  class(xx)
  
  ss = xx %>% group_by(V1, V2) %>% summarise(
    exp_l2_mean = mean(V4) / 100,
    exp_l2_std = sqrt(var(V4) / n()) / 100,
    kde_l2_mean = mean(V5) / 100,
    kde_l2_std = sqrt(var(V5) / n()) / 100
  )
  head(ss)
  
  postrisk.df = rbind.data.frame(postrisk.df,
                                      data.frame(
                                        n = ss$V1,
                                        m = ss$V2,
                                        L2 = ss$exp_l2_mean,
                                        sim_setting = names(sim.results)[l]
                                      ))
  
}
postrisk.df$m = as.factor(postrisk.df$m)
postrisk.df$sim_setting = factor(x = postrisk.df$sim_setting,
                                 levels = names(sim.results))
# postrisk.df = postrisk.df[postrisk.df$m!=50,]

m.label = paste0('m: ', levels(postrisk.df$m))
names(m.label) = levels(postrisk.df$m)

sim_setting.color = scales::hue_pal()(length(sim.results))
sim_setting.pointshape = c(4, 15, 16)
names(sim_setting.color) = names(sim_setting.pointshape) =
  names(sim.results)

ggplot(data = postrisk.df) +
  facet_grid(.~m, scales="fixed",
             labeller = ggplot2::labeller(m = m.label)) +
  # facet_wrap(.~m, scales="free", nrow=1, ncol=4, labeller=label_both) +
  # facet_wrap(vars(sim_setting), nrow=1, ncol=4, labeller=label_both) +
  geom_line(aes(x=n, y=L2, color=sim_setting), linewidth=1) + 
  geom_point(aes(x=n, y=L2, color=sim_setting, shape=sim_setting),
             size=5, stroke=2) + 
  scale_color_manual(values = sim_setting.color,
                     labels = sim_setting.label) +
  scale_shape_manual(values = sim_setting.pointshape,
                     labels = sim_setting.label) +
  scale_x_log10() + 
  scale_y_log10() +
  # ylab("square L2") +
  # theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=25,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 20,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 20),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 25,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(2, "cm"),
    legend.key.height = ggplot2::unit(2, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.spacing.x = ggplot2::unit(.5, 'cm'),
    legend.text=ggplot2::element_text(size=25),
    legend.position = 'bottom'
  ) +
  labs(#title = 'alpha=0.01',
    x = 'n', y = 'Posterior Risk')   # 13 X 6


