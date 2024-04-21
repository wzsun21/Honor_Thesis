library(panelPomp)
library(pomp)
library(doParallel)
library(foreach)
library(doRNG)
library(tidyverse)

# Data Import and Manipulation, and Parallel Computing Setup --------------

brazil <- read.csv("brazil.csv")
brazil <- brazil[,c(1,6)]
brazil$Date <- as.Date(brazil$Date)

india <- read.csv("india.csv")
india <- india[,c(1,5)]
india$Date <- as.Date(india$Date, format = "%d-%b-%y")

indonesia <- read.csv("indonesia.csv")
indonesia <- indonesia[,c(1,6)]
indonesia$Date <- as.Date(indonesia$Date)

mexico <- read.csv("mexico.csv")
mexico <- mexico[,c(1,6)]
mexico$Date <- as.Date(mexico$Date)

emerging <- merge(brazil, india, by = "Date", all = FALSE)
emerging <- merge(emerging, indonesia, by = "Date", all = FALSE)
emerging <- merge(emerging, mexico, by = "Date", all = FALSE)

colnames(emerging)[2:5] <- c("brazil","india","indonesia","mexico")

for (i in 2:5){
  emerging[,i] <- as.numeric(emerging[,i])
}

emerging <- na.omit(emerging)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()  
registerDoParallel(cores)
registerDoRNG(34118892)

brazil_logrt <- log(emerging$brazil[-1]/emerging$brazil[-nrow(emerging)])
india_logrt <- log(emerging$india[-1]/emerging$india[-nrow(emerging)])
indonesia_logrt <- log(emerging$indonesia[-1]/emerging$indonesia[-nrow(emerging)])
mexico_logrt <- log(emerging$mexico[-1]/emerging$mexico[-nrow(emerging)])
panel_data <- cbind(brazil_logrt,india_logrt,indonesia_logrt,mexico_logrt)
colnames(panel_data) <- c("Brazil","India","Indonesia","Mexico")

# Create POMP Object ------------------------------------------------------
create_pomp_object <- function(unit, rproc_type = 'filter') { 
  # unit by name
  
  panel_statenames <- c("V", "S")
  panel_rp_names <- c("mu", "kappa", "theta", "xi", "rho")
  panel_ivp_names <- c("v0")
  panel_parameters <- c(panel_rp_names, panel_ivp_names)
  panel_covarnames <- "covaryt"
  
  rproc1 <- "
  double dWv, dZ, dWs, rt;
  
  rt=covaryt;
  dWs = (rt-mu+0.5*V)/(sqrt(V));
  dZ = rnorm(0, 1);
  
  dWv = rho * dWs + sqrt(1 - rho * rho) * dZ;

  S += S * (mu + sqrt(fmax(V, 0.0)) * dWs);
  V += kappa * (theta - V) + xi * sqrt(fmax(V, 0.0)) * dWv;
  
  if (V<=0) {
    V=1e-32;
  } 
"
  if (unit == 'Brazil'){
    panel_rinit <- "
    V = v0;
    S = 51127;"
    }
  
  if (unit == 'India') {
    panel_rinit <- "
    V = v0;
    S = 22715.33;"
  }
  if (unit == 'Indonesia') {
    panel_rinit <- "
    V = v0;
    S = 4765.729;"
  }
  if (unit == 'Mexico') {
    panel_rinit <- "
    V = v0;
    S = 40447.96;"
  }
  
  panel_rmeasure_filt <- "
  y=exp(covaryt);
"
  
  panel_rmeasure_sim <- "
  y = (mu - 0.5 * V) + sqrt(V);
"
  
  panel_dmeasure <- "
   lik=dnorm(y, mu-0.5*V, sqrt(V), give_log);
"
  # Stock-specific data
  unit_data <- panel_data[, unit, drop=FALSE]
  
  my_ToTrans <- "
       T_xi = log(xi);
       T_kappa = log(kappa);
       T_theta = log(theta);
       T_v0 = log(v0);
       T_mu = log(mu);
       T_rho = log((rho + 1) / (1 - rho));
    "

  my_FromTrans <- "
      kappa = exp(T_kappa);
      theta = exp(T_theta);
      xi = exp(T_xi);
      v0 = exp(T_v0);
      mu = exp(T_mu);
      rho = -1 + 2 / (1 + exp(-T_rho));
    "

  panel_partrans <- parameter_trans(
    toEst = Csnippet(my_ToTrans),
    fromEst = Csnippet(my_FromTrans)
  )
  

# Create Filter -----------------------------------------------------------
  panel.filt <- pomp(data=data.frame(
    y=as.numeric(unit_data), time=1:nrow(unit_data)
  ),
  statenames = panel_statenames,
  paramnames = panel_parameters,
  covarnames = panel_covarnames,
  times = "time",
  t0=0,
  covar = covariate_table(
    time=0:nrow(unit_data),
    covaryt=c(0, as.numeric(unit_data)),
    times = "time",
    order="constant"
  ),
  rmeasure = Csnippet(panel_rmeasure_filt),
  dmeasure = Csnippet(panel_dmeasure),
  rprocess = discrete_time(step.fun = Csnippet(rproc1), delta.t = 1),
  rinit = Csnippet(panel_rinit),
  partrans = panel_partrans
  )
  
  return(panel.filt)
}

# Create PanelPOMP object -------------------------------------------------
create_panelPomp <- function(stocks, specific_names) {
  my_pomps <- list()
  for (stock in stocks) {
    my_pomps[[stock]] <- create_pomp_object(unit = stock, rproc_type = 'filter')  # Your function here
  }
  
  # Create a matrix of specific parameters, number of rows is the number of 
  # specific parameters, number of columns is the number of stocks
  specific_values <- matrix(
    nrow = length(specific_names),
    ncol = length(stocks)
  )
  
  # Name the rows and columns using the parameter and stock names, respectively
  rownames(specific_values) <- specific_names
  colnames(specific_values) <- stocks
  
  # Parameter vector containing all starting parameters (somewhat arbitrary)
  # Here We used results from previous analysis in 1-d model
  all_params <- c(
    'kappa' = 0.03626649,
    'rho' = 0.0007916, 
    'theta' = 8.411001e-05,
    'v0' = 1.508205e-05,
    'xi' = 0.001795454,
    'mu' = 3.389192e-04
  )
  
  # Get the names of the shared parameters by looking at which parameters 
  # are NOT in specific_names
  shared_names <- names(all_params)[!names(all_params) %in% specific_names]
  
  # Loop through all specific_names and name the rows of specific_values 
  # using those names. 
  for (sp in specific_names) {
    specific_values[sp, ] <- all_params[sp]
  }
  
  # Get all of the shared values 
  shared_values <- all_params[shared_names]
  
  # List containing both shared an unit-specific parameters 
  output <- list()
  
  # Add shared an unit specific parameters here. 
  output$shared <- shared_values
  output$specific <- specific_values
  
  # Create and return panelPomp object. 
  return(
    panelPomp(
      my_pomps, 
      shared = output$shared, 
      specific = output$specific
    )
  )
}

# Preparation to call MIF -------------------------------------------------

RUN_LEVEL = 3
NCORES = 36 ## SUBJECT TO CHANGE IF USE HPC 

NREP = 1000
N_STOCKS = ncol(panel_data) # minus the date column

# name of stocks; panel_data is the data of stock prices
names_vec <- colnames(panel_data)[colnames(panel_data) != 'Date']

shared_params <- c()  # Just have to change this vector to specify the shared-parameters
# In case of all-specific, leave the vector empty as here

NP <- switch(RUN_LEVEL, 25, 100, 500)
NMIF <- switch(RUN_LEVEL, 2, 20, 200)
NREPS_EVAL <- switch(RUN_LEVEL, 2, 5, 24)
NP_EVAL <- switch(RUN_LEVEL, 50, 200, 1000)

all_params <- c("mu", "kappa", "theta", "xi", "rho", "v0")
specific_params <- setdiff(all_params, shared_params)

# Create the panel pomp object
my_panel_pomp_model <- create_panelPomp(stocks = names_vec, specific_names = specific_params)

# Create unit-specific search box -----------------------------------------

lower_bounds <- c(
  'mu' = 1e-5, 'theta' = 1e-5, 'v0' = 1e-7,
  'rho' = -1e-2, 'xi'=1e-5, 'kappa'=1e-4
)

upper_bounds <- c(
  "mu" = 1e-3, 'theta' = 1e-3, 'v0' = 1e-5,
  'rho' = 0, 'xi'=1e-2, 'kappa'=1e-1
)

# specific parameter
starts_units <- runif_design(
  lower = lower_bounds[specific_params],
  upper = upper_bounds[specific_params],
  nseq = NREP * N_STOCKS
)

specific_starts_list <- list()
# For each rep, we need a NPARAMS x N_STOCKS matrix.
for (i in 1:NREP) {
  tmp <- t(starts_units)[, ((i - 1) * N_STOCKS + 1):((i) * N_STOCKS)]
  # Get the names of the panelPomp object
  colnames(tmp) <- names(my_panel_pomp_model)
  specific_starts_list[[i]] <- tmp
}

# Create shared parameter search box --------------------------------------

shared_params_start <- runif_design(
  lower = lower_bounds[shared_params],
  upper = upper_bounds[shared_params],
  nseq = NREP
)

# Define RW_SD ------------------------------------------------------------

panel_rw.sd_rp <- 0.01
panel_rw.sd_ivp <- 0.1
panel_cooling.fraction50 <- 0.5
panel_rw.sd <- rw_sd(
  mu = panel_rw.sd_rp,
  theta = panel_rw.sd_rp,
  kappa = panel_rw.sd_rp,
  xi = panel_rw.sd_rp,
  rho = panel_rw.sd_rp,
  v0 = ivp(panel_rw.sd_ivp)
)

# Do global Search --------------------------------------------------------

file_name <- paste0('RUN_LEVEL_', RUN_LEVEL, '/panelPompBPIF', 'FILE_NAME', '.rds') # specify FILE_NAME 
file_name_evals <- paste0('RUN_LEVEL_', RUN_LEVEL, 'evals_ind', 'FILE_NAME', '.rds')

# Save all bpif objects to file `file_name`
global_bpif <- bake(file = file_name, {
  foreach(i=1:NREP) %dopar% {
    tmp_shared_params <- unlist(shared_params_start[i, , drop = FALSE]) # Comment out this line if all-specific
    tmp_specific_params <- specific_starts_list[[i]]
  
    mif2(
      my_panel_pomp_model,
      Np = NP,
      Nmif = NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = panel_rw.sd,
      cooling.type = 'geometric',
      shared.start = tmp_shared_params, # Comment out this line if all-specific
      specific.start = tmp_specific_params,
      block = TRUE # use BPIF to improve efficiency
    )
  }
})

# Save the evaluations + parameter values to `file_name_evals`
  bake(
  file = file_name_evals, {

    final_pars <- as.data.frame(t(sapply(global_bpif, coef)))
    lik_evals <- foreach(i=1:length(global_bpif), .combine = rbind) %dopar% {
      reps <- t(replicate(
        NREPS_EVAL,
        unitlogLik(pfilter(global_bpif[[i]], Np = NP_EVAL))
      )) #Used t
        apply(reps, 2, logmeanexp)
      #panel_logmeanexp(reps, MARGIN = 1, se = TRUE)
    }

    # change to report each individual ll
    # Here only report the added sum of ll's: final_pars$logLike_i <- like_evals[, i]
    final_pars$logLik <- lik_evals[, 1]+lik_evals[, 2]+lik_evals[,3]+lik_evals[,4]

    final_pars
  })
  
