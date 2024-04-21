all_profs <- readRDS('ppomp_profs.rds') %>% 
  as.data.frame()

rho_profs <- readRDS("ppomp_rho_prof.rds") %>% 
  as.data.frame() %>%
  # select(rho_prof, starts_with("ll")) %>% 
  mutate(prof_var = 'rho') %>% 
  rename(prof_val = rho_prof)

all_profs <- bind_rows(all_profs, rho_profs)

best_mods <- all_profs %>% 
  group_by(prof_var, prof_val) %>%
  summarize(
    ll_Mexico = ll_Mexico[which.max(ll_Mexico)],
    ll_India = ll_India[which.max(ll_India)],
    ll_Indonesia = ll_Indonesia[which.max(ll_Indonesia)],
    ll_Brazil = ll_Brazil[which.max(ll_Brazil)], 
    `mu[Mexico]` = `mu[Mexico]`[which.max(ll_Mexico)],
    `mu[India]` = `mu[India]`[which.max(ll_India)],
    `mu[Indonesia]` = `mu[Indonesia]`[which.max(ll_Indonesia)],
    `mu[Brazil]` = `mu[Brazil]`[which.max(ll_Brazil)],
    `kappa[Mexico]` = `kappa[Mexico]`[which.max(ll_Mexico)],
    `kappa[India]` = `kappa[India]`[which.max(ll_India)],,
    `kappa[Indonesia]` = `kappa[Indonesia]`[which.max(ll_Indonesia)],
    `kappa[Brazil]` = `kappa[Brazil]`[which.max(ll_Brazil)],
    `rho[Mexico]` = `rho[Mexico]`[which.max(ll_Mexico)],,
    `rho[India]` = `rho[India]`[which.max(ll_India)],,
    `rho[Indonesia]` = `rho[Indonesia]`[which.max(ll_Indonesia)],
    `rho[Brazil]` = `rho[Brazil]`[which.max(ll_Brazil)],
    `xi[Mexico]` = `xi[Mexico]`[which.max(ll_Mexico)],,
    `xi[India]` = `xi[India]`[which.max(ll_India)],,
    `xi[Indonesia]` = `xi[Indonesia]`[which.max(ll_Indonesia)],
    `xi[Brazil]` = `xi[Brazil]`[which.max(ll_Brazil)],
    `v0[Mexico]` = `v0[Mexico]`[which.max(ll_Mexico)],,
    `v0[India]` = `v0[India]`[which.max(ll_India)],,
    `v0[Indonesia]` = `v0[Indonesia]`[which.max(ll_Indonesia)],
    `v0[Brazil]` = `v0[Brazil]`[which.max(ll_Brazil)],
    `theta[Mexico]` = `theta[Mexico]`[which.max(ll_Mexico)],,
    `theta[India]` = `theta[India]`[which.max(ll_India)],,
    `theta[Indonesia]` = `theta[Indonesia]`[which.max(ll_Indonesia)],
    `theta[Brazil]` = `theta[Brazil]`[which.max(ll_Brazil)]
  ) %>%
  ungroup()


for (i in 1:nrow(best_mods)) {
  best_mods[i, is.na(best_mods[i, ] |> unlist())] <- best_mods[i, 'prof_val']
}


best_mod_results <- best_mods %>% 
  mutate(best_ll = ll_Mexico + ll_Brazil + ll_India + ll_Indonesia) %>%
  group_by(prof_var) %>% 
  slice_max(order_by = best_ll, n = 1)

global_search_results <- readRDS("ppomp_specific.rds") %>% 
  summarize(
    ll_Mexico = ll_Mexico[which.max(ll_Mexico)],
    ll_India = ll_India[which.max(ll_India)],
    ll_Indonesia = ll_Indonesia[which.max(ll_Indonesia)],
    ll_Brazil = ll_Brazil[which.max(ll_Brazil)], 
    `mu[Mexico]` = `mu[Mexico]`[which.max(ll_Mexico)],
    `mu[India]` = `mu[India]`[which.max(ll_India)],
    `mu[Indonesia]` = `mu[Indonesia]`[which.max(ll_Indonesia)],
    `mu[Brazil]` = `mu[Brazil]`[which.max(ll_Brazil)],
    `kappa[Mexico]` = `kappa[Mexico]`[which.max(ll_Mexico)],
    `kappa[India]` = `kappa[India]`[which.max(ll_India)],,
    `kappa[Indonesia]` = `kappa[Indonesia]`[which.max(ll_Indonesia)],
    `kappa[Brazil]` = `kappa[Brazil]`[which.max(ll_Brazil)],
    `rho[Mexico]` = `rho[Mexico]`[which.max(ll_Mexico)],,
    `rho[India]` = `rho[India]`[which.max(ll_India)],,
    `rho[Indonesia]` = `rho[Indonesia]`[which.max(ll_Indonesia)],
    `rho[Brazil]` = `rho[Brazil]`[which.max(ll_Brazil)],
    `xi[Mexico]` = `xi[Mexico]`[which.max(ll_Mexico)],,
    `xi[India]` = `xi[India]`[which.max(ll_India)],,
    `xi[Indonesia]` = `xi[Indonesia]`[which.max(ll_Indonesia)],
    `xi[Brazil]` = `xi[Brazil]`[which.max(ll_Brazil)],
    `v0[Mexico]` = `v0[Mexico]`[which.max(ll_Mexico)],,
    `v0[India]` = `v0[India]`[which.max(ll_India)],,
    `v0[Indonesia]` = `v0[Indonesia]`[which.max(ll_Indonesia)],
    `v0[Brazil]` = `v0[Brazil]`[which.max(ll_Brazil)],
    `theta[Mexico]` = `theta[Mexico]`[which.max(ll_Mexico)],,
    `theta[India]` = `theta[India]`[which.max(ll_India)],,
    `theta[Indonesia]` = `theta[Indonesia]`[which.max(ll_Indonesia)],
    `theta[Brazil]` = `theta[Brazil]`[which.max(ll_Brazil)]
  ) %>% 
  mutate(
    best_ll = ll_India + ll_Mexico + ll_Indonesia + ll_Brazil,
    prof_var = 'none', prof_val = NA
  )

final_results <- bind_rows(best_mod_results, global_search_results) %>% 
  arrange(desc(best_ll))

saveRDS(final_results, 'best_model_params_and_likelihoods.rds')

# Get Smoothed MLE --------------------------------------------------------

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

results_long2 <- all_profs %>%
  rename(logLik_se = ll_se) %>%
  # filter(logLik > 11400) %>%
  # filter(logLik_se < 2) %>%
  # filter(variable == prof_var) %>%
  group_by(prof_var, prof_val) %>%
  mutate(logLik = max(ll_Brazil) + max(ll_India) + max(ll_Indonesia) + max(ll_Mexico)) %>%
  slice_max(order_by = logLik, n = 1) %>% 
  ungroup()

all_mcaps2 <- list()
all_ci2 <- list()
for (p in prof_params) {
  
  mcap_tmp2 <- mcap(
    logLik = results_long2 %>% filter(prof_var == p) %>% pull(logLik),
    parameter = results_long2 %>% filter(prof_var == p) %>% pull(prof_val)
  )
  
  all_mcaps2[[p]] <- mcap_tmp2
  tmp_ci2 <- mcap_tmp2$ci
  names(tmp_ci2) <- c("lower", "upper")
  
  all_ci2[[p]] <- c(tmp_ci2, 'mle' = mcap_tmp2$mle)
}


all_mcaps2 <- purrr::map_df(
  prof_params,
  function(x) all_mcaps2[[x]]$fit %>% mutate(prof_var = x)
)

all_ci2 <- all_ci2 |> as.data.frame() |> t() |> as.data.frame()
all_ci2$prof_var <- rownames(all_ci2)


