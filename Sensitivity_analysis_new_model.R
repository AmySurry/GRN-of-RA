# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)

# T/F statement for activating the pathway
activation = T

# Activation time period - can be a single value or a vector
t_start <- 3600   # Start time for activation
t_end <- 20000     # End time for activation

#t_start <- c(0,14400)   # Start time for activation
#t_end <- c(14400,40000)     # End time for activation


#### Creating the model for the NF-KB pathway ####
NF_KB <- function(time,state,parameters){
  
  # Defining the genes included in the model
  IKKn <- state[1] # IKK complex non-active
  IKKa <- state[2] # IKK complex active
  NFkBn <- state[3] # Nuclear NFkB 
  A20 <- state[4] # A20
  IkBa <- state[5] # IkBa protein 
  IkBat <- state[6] # IkBa mRNA 
  
  
  # Defining the parameters
  Kdeg <- parameters["Kdeg"]
  a2 <- parameters["a2"]
  k1 <- parameters["k1"]
  a3 <- parameters["a3"]
  i1a <- parameters["i1a"]
  k3 <- parameters["k3"]
  k2 <- parameters["k2"]
  delta <- parameters["delta"]
  epsilon <- parameters["epsilon"]
  Cdeg <- parameters["Cdeg"]
  C4a <- parameters["C4a"]
  C5a <- parameters["C5a"]
  C3a <- parameters["C3a"]
  
  # Activation of the pathway (on or off)
  Tr <- 0 # Receptor inactivation
  
  if(length(t_start) > 1){ # If there are multiple activation periods
    for(i in 1:length(t_start)){
      if(activation && time >= t_start[i] && time <= t_end[i]){
        Tr <- 1 # Receptor activation
        break
      }
    }
  }else{ # If there is only one activation period
    if(activation && time >= t_start && time <= t_end){
      Tr <- 1 # Receptor activation
    }else{
      Tr <- 0 # Receptor inactivation
    }
  }
  
  # The differential equations
  dIKKn <- Kdeg - (Kdeg*IKKn) - (Tr*k1*IKKn)
  dIKKa <- (Tr*k1*IKKn) - (k3+Kdeg+(Tr*k2*A20))*IKKa
  dNFkBn <- (a3*IKKa)*(1-NFkBn)*(delta/(IkBa+delta)) - (i1a*IkBa)*(NFkBn/(NFkBn+epsilon))
  dA20 <- (Cdeg*NFkBn) - (Cdeg*A20)
  dIkBa <- (C4a*IkBat) - (C5a*IkBa) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa+delta))) - (i1a*IkBa*(NFkBn/(NFkBn+epsilon)))
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat) #- (C4a*IkBat)
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat)
  
  return(list(c(dIKKn,dIKKa,dNFkBn,dA20,dIkBa,dIkBat)))
}


#### Solving the differential equations ####
# Initial condition 
# x0 <- c(A = 0, I = 0, P = 0) 
x0 <- c(IKKn = 1, # µM (Unknown)
        IKKa = 0, # µM (Unknown)
        NFkBn = 0.02, # µM (found in the literature - Hoffmann et al. (2002))
        A20 = 0.02, 
        IkBa = 0.5,
        IkBat = 0.02)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters <- c(k3 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                a2 = 0.0763, # s^-1 IκBα degradation induced by IKK complex 
                #alpha3 = 0.000372 + 0.000106 + 0.0428,  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)
                k1 = 0.00195, # s^-1 Activation of IKK complex
                a3 = 0.0946, # s^-1 IKK (IκBα|NF-κB) association
                i1a = 0.000595, # s^-1 IκBα nuclear import leading to... (repression of NF-κB)
                Kdeg = 0.000107,
                k2 = 0.0357,
                delta = 0.108,
                epsilon = 0.0428,
                Cdeg = 0.000106,
                C4a =0.00313,
                C5a = 0.0000578,
                C3a = 0.000372)


# Time sequence
times <- seq(0, 20000, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters, maxsteps=1000000) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting


# Sobol sensitivity analysis function

sobol_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val, n_iterations = 3000){
  
  # Inputs:
  # 
  # model -> the ODEs which you want to analyze
  # var_pars -> parameters whose sensitivity will be calculated
  # x0_init -> initital value for each state in the system of ODEs
  # var_min -> the smallest value the parameters can have
  # var_max -> the largest value the parameters can have (chosen randomly in the space)
  # time_val -> the time interval, which the sensitivity analysis will take place.
  # n_iterations -> the number of iterations that will be run.
  
  
  # calculating the sobol_result data.frame
  sobol_result <- ODEsobol(mod = model, # the model that is used
                           pars = var_pars, # the parameters present
                           state_init = x0_init, # initial state value
                           times = time_val, # time value
                           n = n_iterations, # number of iterations -> more is better accuracy but also takes time to run
                           rfuncs = "runif", # how to distribute the parameter values - here uniformly
                           rargs = paste0("min = ", var_min, # minimum value the parameters can have
                                          ", max = ", var_max), # maximum value the parameters can have
                           sobol_method = "Martinez", # which sobol method to be used
                           ode_method = "lsoda", # which ode method to be used (this is default in ode)
                           parallel_eval = TRUE, # parallel evaluation, not important to understand
                           parallel_eval_ncores = 2) # number of cores to be used to evaluate, also not important to understand
  
  return(sobol_result)
  
}

# Function with the differential equations for the sensitivity analysis, used as the model parameter in the sobol_sensitivity function
NF_KB_sensitivity <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    R <- 1
    
    # The differential equations
    dIKKn <- Kdeg - (Kdeg*IKKn) - (Tr*k1*IKKn)
    dIKKa <- (Tr*k1*IKKn) - (k3+Kdeg+(Tr*k2*A20))*IKKa
    dNFkBn <- (a3*IKKa)*(1-NFkBn)*(delta/(IkBa+delta)) - (i1a*IkBa)*(NFkBn/(NFkBn+epsilon))
    dA20 <- (Cdeg*NFkBn) - (Cdeg*A20)
    dIkBa <- (C4a*IkBat) - (C5a*IkBa) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa+delta))) - (i1a*IkBa*(NFkBn/(NFkBn+epsilon)))
    dIkBat <- (C3a*NFkBn) - (C3a*IkBat) #- (C4a*IkBat)
    dIkBat <- (C3a*NFkBn) - (C3a*IkBat)
    
    return(list(c(dIKKn,dIKKa,dNFkBn,dA20,dIkBa,dIkBat)))
    
  })
}

#Vector with the named parameters

bound_var <- c("Kdeg",
               "a2",
               "k1",
               "a3",
               "i1a",
               "k3",
               "k2",
               "delta",
               "epsilon",
               "Cdeg",
               "C4a",
               "C5a",
               "C3a")

#Vectors with the minimum and maximum values for the parameters
bound_min_var <- c(0.0001,
                   0.05,
                   0.001,
                   0.05,
                   0.0005,
                   0.001,
                   0.003,
                   0.1,
                   0.04,
                   0.0001,
                   0.003,
                   0.00005,
                   0.0003)
bound_max_var <- c(0.0003,
                   0.1,
                   0.003,
                   0.15,
                   0.0008,
                   0.003,
                   0.006,
                   0.3,
                   0.08,
                   0.0003,
                   0.005,
                   0.00008,
                   0.0006)


# Time sequence
times <- seq(0.01, 1000, by = 1) 

# Doing the sensitivity analysis
sensitive_sobol_final <- sobol_sensitivity(NF_KB_sensitivity,
                                           bound_var,
                                           x0,
                                           bound_min_var,
                                           bound_max_var,
                                           times)

#Plottig the sensitivity analysis for all

#For IKKn 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IKKn", main_title = "IKKn sensitivity - SOBOL", type = "l", lwd = 3)
#For IKKa 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IKKa", main_title = "IKKa sensitivity - SOBOL", type = "l", lwd = 3)
#For NFkBn
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "NFkBn", main_title = "NFkBn sensitivity - SOBOL", type = "l", lwd = 3)
#For A20 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "A20", main_title = "A20 sensitivity - SOBOL", type = "l", lwd = 3)
#For IkBa 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IkBa", main_title = "IkBa sensitivity - SOBOL", type = "l", lwd = 3)
#For IkBat
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IkBat", main_title = "IkBat sensitivity - SOBOL", type = "l", lwd = 3)

#plotting the sensitivity analysis for betas only

#For IKKn 
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "IKKn", main_title = "IKKn sensitivity - SOBOL", type = "l", lwd = 3)
#For IKKa 
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "IKKa", main_title = "IKKa sensitivity - SOBOL", type = "l", lwd = 3)
#For NFkBn
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "NFkBn", main_title = "NFkBn sensitivity - SOBOL", type = "l", lwd = 3)
#For A20 
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "A20", main_title = "A20 sensitivity - SOBOL", type = "l", lwd = 3)
#For IkBa 
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "IkBa", main_title = "IkBa sensitivity - SOBOL", type = "l", lwd = 3)
#For IkBat
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "IkBat", main_title = "IkBat sensitivity - SOBOL", type = "l", lwd = 3)
