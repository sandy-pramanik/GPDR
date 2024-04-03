

rm(list = ls())


# sourcing library and codes ----
sourcecode.path = ...    # specifies path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"


# data generation specifics ----
n=100
m=10
alpha_DP = .1
k = 10
xseq = seq(0, 1, by = 0.01)
beta.truth = beta1(xseq)


# specifics for KDE ----
set.seed(0)
data = generate_data_dp(beta1, 100, 100)
K = MK(nu=Inf, sigma=1, l=0.25)
## x grid and final evaluate grid are always the same
## for all m and n, therefore we only need to calculate once
Kxx = K(data$kde[[1]]$x, data$kde[[1]]$x)
Kxc = K(xseq, data$kde[[1]]$x)
Kcc = K(xseq, xseq)


# simulate one data using DP(alpha) and fit GPDR and KDE ----
## simulate
set.seed(1)
data = generate_data_dp(beta = beta1, n = n, m = m, alpha = alpha_DP)

## model fitting
fit0 = GPfit(data$X, data$y, nu=Inf, l=0.25, k=k, sigma=0.1, verbose=FALSE)
fit1 = kde_fit(data, Kxx, Kxc, Kcc, sig=0.1, newx=xseq)

## estimated regression functions from GPDR and KDE
pred0 = predict(fit0, xseq, mean.only=TRUE)$E
pred1 = fit1$E

## empirical squared error risk
(exp_l2 = sum((beta.truth - pred0)^2))
(kde_l2 = sum((beta.truth - pred1)^2))


## plotting to compare regression function estimates from GPDR and KDE with truth ====
plotdf = rbind.data.frame(data.frame('x' = xseq, 'fx' = beta.truth, 'est_type' = 'Truth'),
                          data.frame('x' = xseq, 'fx' = pred0, 'est_type' = 'GPDR'),
                          data.frame('x' = xseq, 'fx' = pred1, 'est_type' = 'KDE'))
head(plotdf)

plotdf$est_type = factor(x = plotdf$est_type,
                         levels = c('Truth', 'GPDR', 'KDE'))

est_type.label = c('Truth' = 'Truth',
                   'GPDR' = paste0('GPDR (Risk: ', round(exp_l2, 2), ')'), 
                   'KDE' = paste0('KDE (Risk: ', round(kde_l2, 2), ')'))

est_type.colors = c('Truth' = 'red3', 'GPDR' = 1, 'KDE' = 4)

fig =
  ggplot2::ggplot(data = plotdf) +
  # ggplot2::coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  ggplot2::geom_line(
    ggplot2::aes(x = x, y = fx, color = est_type),
    linewidth = 1
  ) +
  ggplot2::scale_color_manual(values = est_type.colors,
                              labels = est_type.label) +
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
                                                      fill = NA, linewidth = 1),
                 panel.background = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                          colour = "grey90"),
                 panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                          colour = "grey90"),
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold"),
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold"),
                 strip.background = ggplot2::element_rect(color="black", linewidth=1),
                 legend.title = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(1, "cm"),
                 legend.key.height = ggplot2::unit(1, "cm"), 
                 # legend.key.size = ggplot2::unit(30, "cm"),
                 legend.spacing.y = ggplot2::unit(5, 'cm'), 
                 legend.text=ggplot2::element_text(size=15),
                 legend.position = 'bottom') +
  # ggplot2::guides(color = ggplot2::guide_legend(
  #   nrow = 2, byrow=FALSE#,
  #   # override.aes = list(
  #   #   linetype = "solid",
  #   #   shape = c(rep(16, length(levels(plotdf$method))-1), NA))
  # )) +
  ggplot2::labs(title = paste0('n=', n, ', m=', m), 
                # subtitle = paste0('n=', n),
                x = 'x', y = 'Regression Function')   # 8 X 6
fig

# fig.m10 = fig
# fig.m100 = fig
# 
# ggpubr::ggarrange(fig.m10, fig.m100, nrow = 1, ncol = 2,
#                   legend = 'bottom')  # 6.5 X 13
