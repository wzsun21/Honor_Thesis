library(pomp)

sp500_statenames <- c("V", "S", 'S1', 'R')
sp500_rp_names <- c("mu", "kappa", "theta", "xi", "rho")
sp500_ivp_names <- c("V_0")
sp500_parameters <- c(sp500_rp_names, sp500_ivp_names)
sp500_covarnames <- "covaryt"

rproc_filter <- "
  double dWv, dZ, dWs, rt;

  rt = covaryt;

  dWs = (rt-mu+0.5*V)/(sqrt(V));
  dZ = rnorm(0, 1);

  dWv = rho * dWs + sqrt(1 - rho * rho) * dZ;

  S1 = S1 * rt; // given exactly
  S += S * (mu + sqrt(fmax(V, 0.0)) * dWs);  // given dynamically, or stochastically
  V += kappa * (theta - V) + xi * sqrt(fmax(V, 0.0)) * dWv;

  if (V<=0) {
    V=1e-32;
  }
"

# This is for simulation
rproc_sim <- "
  double dWv, dZ, dWs;

  dWs = rnorm(0, 1);
  dZ = rnorm(0, 1);

  dWv = rho * dWs + sqrt(1 - rho * rho) * dZ;

  S += S * (mu + sqrt(fmax(V, 0.0)) * dWs);  // given dynamically, or stochastically
  R = (mu - 0.5 * V) + sqrt(V) * dWs;
  V += kappa * (theta - V) + xi * sqrt(fmax(V, 0.0)) * dWv;

  if (V<=0) {
    V=1e-32;
  }
"

sp500_rinit <- "
  V = V_0;
  S = 1106.75; // TODO: change this for the first return, start fitting on second?
  S1= 1106.75;
  R = 0;
"

# For fitting, neither rmeasure matters.
#   sp500_rmeasure_filt <- "
#   y=covaryt;
# "

  sp500_rmeasure <- "
  y = R;
"

sp500_dmeasure <- "
   lik=dnorm(y, mu-0.5*V, sqrt(V), give_log);
"

my_ToTrans <- "
   T_xi = log(xi);
   T_kappa = log(kappa);
   T_theta = log(theta);
   T_V_0 = log(V_0);
   T_rho = log((rho + 1) / (1 - rho));
"

my_FromTrans <- "
  kappa = exp(T_kappa);
  theta = exp(T_theta);
  xi = exp(T_xi);
  V_0 = exp(T_V_0);
  rho = -1 + 2 / (1 + exp(-T_rho));
"

sp500_partrans <- parameter_trans(
  toEst = Csnippet(my_ToTrans),
  fromEst = Csnippet(my_FromTrans)
)

sp500 <- read.csv('SPX.csv') %>%
  # mutate(date = as.Date(Date)) %>%
  # mutate(diff_days = difftime(date, min(date), units = 'day')) %>%
  # mutate(time = as.numeric(diff_days)) %>%
  mutate(time = (-1):(n() - 2)) %>%
  mutate(y = log(Close / lag(Close))) %>%
  select(time, y) %>%
  drop_na()

fit_data <- sp500 %>% select(time, y) %>% filter(time > 0)

# Filter
sp500.filt <- pomp(
  data=fit_data,
  statenames = sp500_statenames,
  paramnames = sp500_parameters,
  covarnames = sp500_covarnames,
  times = "time",
  t0=0,
  covar=covariate_table(
    time=sp500 %>% pull(time),
    covaryt=sp500 %>% pull(y),
    times = "time",
    order = 'constant'
  ),
  # rmeasure = Csnippet(sp500_rmeasure_filt),
  # rmeasure = Csnippet(sp500_rmeasure_sim),
  dmeasure = Csnippet(sp500_dmeasure),
  rprocess = discrete_time(step.fun = Csnippet(rproc_filter), delta.t = 1),
  rinit = Csnippet(sp500_rinit),
  partrans = sp500_partrans,
  params = unlist(results %>% slice_max(order_by = ll, n = 1) %>% select(rho, mu, theta, kappa, V_0, xi)) # Need to specify parameters
)

# Confidence interval of the ending time spot volatility
pfilt_obj <- pfilter(sp500.filt, Np = 1000, save.states = TRUE)
plot(pfilt_obj)
quantile(sqrt(pfilt_obj@saved.states[[3522]] %>% t() %>% as.data.frame() %>% pull(V)) * sqrt(252), probs = c(0.025, 0.975, 0.5))

# Plot volatility 
pfilt_obj <- pfilter(sp500.filt, Np = 1000, filter.mean = TRUE)
pfilt_obj@filter.mean %>% t() %>% as.data.frame() %>% pull(V) %>% 
  plot(x = 1:3522, y = ., type = 'l')

# Date vector
sp500_date <- read.csv('SPX.csv') %>%
  mutate(date = as.Date(Date)) %>%
  # mutate(diff_days = difftime(date, min(date), units = 'day')) %>%
  # mutate(time = as.numeric(diff_days)) %>%
  mutate(time = (-1):(n() - 2)) %>%
  mutate(y = log(Close / lag(Close))) %>%
  select(time, y, date) %>%
  filter(time > 0)

pfilt_obj@filter.mean %>% t() %>% as.data.frame() %>% pull(V) %>% 
  plot(x = sp500_date$date, y = ., type = 'l')

# Plot estimated volatility--------------------------

emerging_vol <- matrix(nrow = 2190, ncol = 2) # specify appropriate nrow
for (i in 1:2190) {
  emerging_vol[i,1] <- mean(sqrt(pfilt_obj@saved.states[[i]] %>% t() %>% as.data.frame() %>% pull(V)) * sqrt(252))
  emerging_vol[i,2] <- mean(sqrt(pfilt_obj2@saved.states[[i]] %>% t() %>% as.data.frame() %>% pull(V)) * sqrt(252))
}
emerging_vol <- as.data.frame(emerging_vol)
emerging_vol$Date <- emerging$Date[-1]
colnames(emerging_vol)[1:2] <- c("mean_1","mean_2")

# benchmark for comparison: for SP500, we have VIX; for others, we may use this. Take mexico model as an example.
#r <- 0.0001615378
#ref <- mexico_logrt-r
#emerging_vol$ref <- sqrt(ref^2)

# This is the plot for Mexico model
# Not for SP500
# Just serves as an example

# ggplot(data=emerging_vol)+
#   geom_line(aes(x=Date,y=mean_1, col="Estimated Volatility from Shared-Rho Model"),alpha=0.5)+
#   geom_line(aes(x=Date,y=mean_2, col="Estimated Volatility from Shared-Kappa Model"),alpha=0.5)+
#   geom_line(aes(x=Date,y=ref, col="Reference"))+
#   geom_hline(aes(yintercept=0.1672, col = "Mean Reversion Level from Shared-Rho Model"),linetype='dashed')+
#   geom_hline(aes(yintercept=0.1567, col = "Mean Reversion Level from Shared-Kappa Model"), linetype='dashed')+
#   theme_bw()+
#   theme(legend.position="bottom")+
#   xlab("Time")+
#   ylab("Volatility Estimate")+
#   scale_color_manual(values=c("#F8766D","#00BA38","#619CFF","red","blue"))+
#   guides(color = guide_legend(nrow = 2))
