library(tidyverse)
library(latex2exp)

all_profs <- readRDS('ppomp_profs.rds') %>% 
  as.data.frame() %>% 
  select(prof_var, prof_val, starts_with("ll"))

rho_profs <- readRDS("ppomp_rho_prof.rds") %>% 
  as.data.frame() %>%
  select(rho_prof, starts_with("ll")) %>% 
  mutate(prof_var = 'rho') %>% 
  rename(prof_val = rho_prof)

all_profs <- bind_rows(all_profs, rho_profs)
prof_params <- unique(all_profs$prof_var)


results_long1 <- all_profs %>%
  rename(logLik = ll_total, logLik_se = ll_se) %>%
  # filter(logLik > 11400) %>%
  # filter(logLik_se < 2) %>%
  # filter(variable == prof_var) %>%
  group_by(prof_var, prof_val) %>%
  slice_max(order_by = logLik, n = 1) %>% 
  ungroup()

results_long2 <- all_profs %>%
  rename(logLik_se = ll_se) %>%
  # filter(logLik > 11400) %>%
  # filter(logLik_se < 2) %>%
  # filter(variable == prof_var) %>%
  group_by(prof_var, prof_val) %>%
  mutate(logLik = max(ll_Brazil) + max(ll_India) + max(ll_Indonesia) + max(ll_Mexico)) %>%
  slice_max(order_by = logLik, n = 1) %>% 
  ungroup()

all_mcaps <- list()
all_ci <- list()
all_mcaps2 <- list()
all_ci2 <- list()
for (p in prof_params) {
  
  # else if (p == 'theta') {
  #     mcap_tmp <- mcap(
  #       logLik = results_long %>% filter(prof_var == p) %>% pull(logLik),
  #       parameter = results_long %>% filter(prof_var == p) %>% pull(value),
  #       span = 0.4
  #     )
  #   } else if (p == 'rho') {
  #     mcap_tmp <- mcap(
  #       logLik = results_long %>% filter(prof_var == p) %>% pull(logLik),
  #       parameter = results_long %>% filter(prof_var == p) %>% pull(value),
  #       span = 0.45
  #     )
  #   } else {
  mcap_tmp <- mcap(
    logLik = results_long1 %>% filter(prof_var == p) %>% pull(logLik),
    parameter = results_long1 %>% filter(prof_var == p) %>% pull(prof_val)
  )
  
  mcap_tmp2 <- mcap(
    logLik = results_long2 %>% filter(prof_var == p) %>% pull(logLik),
    parameter = results_long2 %>% filter(prof_var == p) %>% pull(prof_val)
  )
  
  all_mcaps[[p]] <- mcap_tmp
  tmp_ci <- mcap_tmp$ci
  names(tmp_ci) <- c("lower", "upper")
  
  all_ci[[p]] <- c(tmp_ci, 'mle' = mcap_tmp$mle)
  
  all_mcaps2[[p]] <- mcap_tmp2
  tmp_ci2 <- mcap_tmp2$ci
  names(tmp_ci2) <- c("lower", "upper")
  
  all_ci2[[p]] <- c(tmp_ci2, 'mle' = mcap_tmp2$mle)
}

all_mcaps <- purrr::map_df(
  prof_params,
  function(x) all_mcaps[[x]]$fit %>% mutate(prof_var = x)
)

all_mcaps2 <- purrr::map_df(
  prof_params,
  function(x) all_mcaps2[[x]]$fit %>% mutate(prof_var = x)
)

all_ci <- all_ci |> as.data.frame() |> t() |> as.data.frame()
all_ci$prof_var <- rownames(all_ci)

all_ci2 <- all_ci2 |> as.data.frame() |> t() |> as.data.frame()
all_ci2$prof_var <- rownames(all_ci2)

# all_mles <- all_results %>% 
#   slice_max(order_by = ll) %>% 
#   select(ll, all_of(prof_params)) %>% 
#   rename(logLik = ll) %>%
#   # mutate(V_0 = sqrt(V_0)) %>%
#   pivot_longer(
#     cols = -logLik,
#     names_to = 'prof_var',
#     values_to = 'value'
#   )

all_specific_best <- readRDS("ppomp_specific.rds") %>% 
  summarize(best_ll = max(ll_Brazil) + max(ll_India) + max(ll_Indonesia) + max(ll_Mexico)) %>% 
  unlist()

# prof_params <- c("mu", 'theta', 'kappa', 'V_0', 'xi', 'rho', 'lambda')
h1_plots <- list()
for (i in 1:(length(prof_params))) {
  p = prof_params[i]
  if (p == 'mu') {
    my_lab <- TeX('$\\mu$')
  } else if (p == 'theta') {
    my_lab <- TeX('$\\theta$')
  } else if (p == 'kappa') {
    my_lab <- TeX('$\\kappa$')
  } else if (p == 'V_0') {
    my_lab <- TeX('$\\sqrt{V_0}$')
  } else if (p == 'xi') {
    my_lab <- TeX('$\\xi$')
  } else if (p == 'rho') {
    my_lab <- TeX('$\\rho')
  }
  
  gg_tmp <- ggplot() +
    # geom_point(
    #   data = results_long1 %>% filter(prof_var == p), 
    #   aes(x = prof_val, y = logLik)
    # ) +
    geom_point(
      data = results_long2 %>% filter(prof_var == p), 
      aes(x = prof_val, y = logLik, col = (logLik >= all_specific_best - qchisq(0.95, 3) / 2))
      # col = 'red'
    ) +
    # geom_line(
    #   data = all_mcaps %>% filter(prof_var == p), 
    #   aes(x = parameter, y = smoothed), 
    #   col = 'blue'
    # ) +
    geom_line(
      data = all_mcaps2 %>% filter(prof_var == p), 
      aes(x = parameter, y = smoothed, col = (smoothed >= all_specific_best - qchisq(0.95, 3) / 2)),
      col = 'blue'
    ) +
    # geom_line(data = all_mcaps, aes(x = parameter, y = quadratic), col = 'red') +
    # geom_vline(
    #   data = all_ci %>% filter(prof_var == p), 
    #   aes(xintercept = lower), linetype = 'dashed'
    # ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = lower), linetype = 'dashed'
      # col = 'red'
    ) +
    # geom_vline(
    #   data = all_ci %>% filter(prof_var == p), 
    #   aes(xintercept = upper), 
    #   linetype = 'dashed'
    # ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = upper), 
      linetype = 'dashed'
      # col = 'red'
    ) +
    # geom_vline(
    #   data = all_ci %>% filter(prof_var == p), 
    #   aes(xintercept = mle), col = 'blue'
    # ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = mle),
      col = 'blue'
    ) +
    geom_hline(yintercept = all_specific_best) +
    geom_hline(yintercept = all_specific_best - qchisq(0.95, 3) / 2, linetype = 'dashed') +
    # This code chunk below will plot the MLE along with the smoothed marginal 
    # MLEs and the confidence interval. 
    # geom_point(
    #   data = all_mles %>% filter(prof_var == p),
    #   aes(x = value, y = logLik), col = 'red'
    # ) +
    labs(y = 'Log-Likelihood', title = my_lab) +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    ) +
    theme_bw() +
    scale_y_continuous(limits = c(27750, 27852)) +
    theme(
      axis.title.x = element_blank(), 
      plot.title = element_text(
        hjust = 0.5, size = 8, margin = margin(t = 0, unit = "pt")
      ),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 6.5),
      plot.margin = margin(5.5, 3, 0, -1, unit = 'pt')
    ) + 
    scale_color_manual(values = c("black", 'red')) + 
    guides(col = FALSE)
  
  if (i %in% c(1, 4)) {
    h1_plots[[i]] <- gg_tmp
  } else {
    h1_plots[[i]] <- gg_tmp + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  }
}

cowplot::plot_grid(
  plotlist = h1_plots,
  align = 'h', ncol = 3,
  rel_widths = c(1.25, 1, 1)
)


h1_plots <- list()
for (i in 1:(length(prof_params))) {
  p = prof_params[i]
  if (p == 'mu') {
    my_lab <- TeX('$\\mu$')
  } else if (p == 'theta') {
    my_lab <- TeX('$\\theta$')
  } else if (p == 'kappa') {
    my_lab <- TeX('$\\kappa$')
  } else if (p == 'V_0') {
    my_lab <- TeX('$\\sqrt{V_0}$')
  } else if (p == 'xi') {
    my_lab <- TeX('$\\xi$')
  } else if (p == 'rho') {
    my_lab <- TeX('$\\rho')
  }
  
  gg_tmp <- ggplot() +
    geom_point(
      data = results_long1 %>% filter(prof_var == p),
      aes(x = prof_val, y = logLik)
    ) +
    geom_point(
      data = results_long2 %>% filter(prof_var == p), 
      aes(x = prof_val, y = logLik, col = (logLik >= all_specific_best - qchisq(0.95, 3) / 2)),
      col = 'red'
    ) +
    geom_line(
      data = all_mcaps %>% filter(prof_var == p),
      aes(x = parameter, y = smoothed),
      col = 'blue'
    ) +
    geom_line(
      data = all_mcaps2 %>% filter(prof_var == p), 
      aes(x = parameter, y = smoothed, col = (smoothed >= all_specific_best - qchisq(0.95, 3) / 2)),
      col = 'red'
    ) +
    # geom_line(data = all_mcaps, aes(x = parameter, y = quadratic), col = 'red') +
    geom_vline(
      data = all_ci %>% filter(prof_var == p),
      aes(xintercept = lower), linetype = 'dashed'
    ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = lower), linetype = 'dashed',
      col = 'red'
    ) +
    geom_vline(
      data = all_ci %>% filter(prof_var == p),
      aes(xintercept = upper),
      linetype = 'dashed'
    ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = upper), 
      linetype = 'dashed',
      col = 'red'
    ) +
    geom_vline(
      data = all_ci %>% filter(prof_var == p),
      aes(xintercept = mle), col = 'blue'
    ) +
    geom_vline(
      data = all_ci2 %>% filter(prof_var == p), 
      aes(xintercept = mle),
      col = 'red'
    ) +
    geom_hline(yintercept = all_specific_best) +
    geom_hline(yintercept = all_specific_best - qchisq(0.95, 3) / 2, linetype = 'dashed') +
    # This code chunk below will plot the MLE along with the smoothed marginal 
    # MLEs and the confidence interval. 
    # geom_point(
    #   data = all_mles %>% filter(prof_var == p),
    #   aes(x = value, y = logLik), col = 'red'
    # ) +
    labs(y = 'Log-Likelihood', title = my_lab) +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    ) +
    theme_bw() +
    scale_y_continuous(limits = c(27750, 27850)) +
    theme(
      axis.title.x = element_blank(), 
      plot.title = element_text(
        hjust = 0.5, size = 8, margin = margin(t = 0, unit = "pt")
      ),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 6.5),
      plot.margin = margin(5.5, 3, 0, -1, unit = 'pt')
    ) + 
    scale_color_manual(values = c("black", 'red')) + 
    guides(col = FALSE)
  
  if (i %in% c(1, 4)) {
    h1_plots[[i]] <- gg_tmp
  } else {
    h1_plots[[i]] <- gg_tmp + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  }
}

cowplot::plot_grid(
  plotlist = h1_plots,
  align = 'h', ncol = 3,
  rel_widths = c(1.25, 1, 1)
)

all_specific <- readRDS("ppomp_specific2.rds")

all_specific_best <- readRDS("ppomp_specific2.rds") %>% 
  summarize(best_ll = max(ll_Brazil) + max(ll_India) + max(ll_Indonesia) + max(ll_Mexico))


rho_tops1 <- rho_res %>% 
  select(ll_total, ll_se, rho) %>% 
  rename(logLik = ll_total, logLik_se = ll_se) %>% 
  filter(rho >= -0.8) %>%
  # filter(logLik > 11000) %>% 
  group_by(rho) %>% 
  slice_max(order_by = logLik, n = 1)

mcap_rho1 <- mcap(
  logLik = rho_tops1$logLik, parameter = rho_tops1$rho,
  span = 0.5
)

rho_tops2 <- rho_res %>% 
  select(ll_total, ll_se, rho, ll_Brazil, ll_India, ll_Indonesia, ll_Mexico) %>% 
  rename(logLik = ll_total, logLik_se = ll_se) %>% 
  filter(rho >= -0.8) %>%
  # filter(logLik > 11000) %>% 
  group_by(rho) %>% 
  mutate(best_ll = max(ll_Brazil) + max(ll_Indonesia) + max(ll_India) + max(ll_Mexico)) %>% 
  slice_max(order_by = best_ll, n = 1)

mcap_rho2 <- mcap(
  logLik = rho_tops2$best_ll, parameter = rho_tops2$rho,
  span = 0.5
)