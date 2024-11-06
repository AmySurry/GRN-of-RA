rm(list=ls())


#Getting Libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)

#Installing the package if not already installed
#install.packages("ODEsensitivity")
#install.packages("devtools")
#install.packages("ggplot2")


#NF_KB model
NF_KB <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    # Receptor activation
    if (activation == T) {
      R <- 1 
    } else {
      R <- 0
    }
    
    # The differential equations
    dA <- R * beta1 - alpha1 * A
    dI <- A * beta3 - alpha2 * I
    dP <- A * beta2 - beta4 * I - alpha3 * P
    
    return(list(c(dA, dI, dP)))
    
  })
}


# Initials
x0 <- c(A = 0, I = 0, P = 0) 

# Parameter values from literature
alpha1 <- 0.00145
alpha2 <- 0.0763
alpha3 <- 0.000372 + 0.000106 + 0.0428
beta1 <- 0.00195
beta2 <- 0.0946
beta3 <- 0.0946
beta4 <- 0.000595

# Parameter vector
parameters <- c(alpha1 = alpha1, 
                alpha2 = alpha2, 
                alpha3 = alpha3,
                beta1 = beta1, 
                beta2 = beta2, 
                beta3 = beta3, 
                beta4 = beta4,
                activation = T)

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB)


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
    dA <- R * beta1 - alpha1 * A
    dI <- A * beta3 - alpha2 * I
    dP <- A * beta2 - beta4 * I - alpha3 * P
    
    return(list(c(dA, dI, dP)))
    
  })
}

#Vector with the named parameters

bound_var <- c("alpha1",
               "alpha2",
               "alpha3",
               "beta1",
               "beta2",
               "beta3",
               "beta4")

#Vectors with the minimum and maximum values for the parameters
bound_min_var <- c(0.001,
                   0.05,
                   0.01,
                   0.0005,
                   0.05,
                   0.05,
                   0.0001)
bound_max_var <- c(0.003,
                   0.1,
                   0.06,
                   0.003,
                   0.1,
                   0.1,
                   0.1)


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

#For A (Ikk complex)
plot(sensitive_sobol_final, pars_plot = bound_var[1:7], state_plot = "A", main_title = "A sensitivity - SOBOL", type = "l", lwd = 3)
#For I (IkBa)
plot(sensitive_sobol_final, pars_plot = bound_var[1:7], state_plot = "I", main_title = "I sensitivity - SOBOL", type = "l", lwd = 3)
#For P (NF-KB)
plot(sensitive_sobol_final, pars_plot = bound_var[1:7], state_plot = "P", main_title = "P sensitivity - SOBOL", type = "l", lwd = 3)


#plotting the sensitivity analysis for betas only

#For A (Ikk complex)
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "A", main_title = "A sensitivity - SOBOL", type = "l", lwd = 3)
#For I (IkBa) 
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "I", main_title = "I sensitivity - SOBOL", type = "l", lwd = 3)
#For P (NF-KB)
plot(sensitive_sobol_final, pars_plot = bound_var[4:7], state_plot = "P", main_title = "P sensitivity - SOBOL", type = "l", lwd = 3)


