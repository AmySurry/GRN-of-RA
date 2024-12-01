# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)
library(RColorBrewer)
library(pbapply)


# Sobol sensitivity analysis function
sobol_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val, n_iterations = 3000){
  
  # Inputs:

  # model -> the ODEs which you want to analyze
  # var_pars -> parameters whose sensitivity will be calculated
  # x0_init -> initial value for each state in the system of ODEs
  # var_min -> the smallest value the parameters can have
  # var_max -> the largest value the parameters can have (chosen randomly in the space)
  # time_val -> the time interval, which the sensitivity analysis will take place.
  # n_iterations -> the number of iterations that will be run.
  
  # Print a message to indicate the start of the process
  cat("Starting Sobol sensitivity analysis...\n")
  
  # Progress bar for tracking iterations
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
                           parallel_eval_ncores = 2, # number of cores to be used to evaluate, also not important to understand
                           progress_bar = TRUE)  # Add progress bar
  
  cat("Sensitivity analysis complete!\n")
  
  return(sobol_result)
  
}

# Function with the differential equations for the sensitivity analysis, used as the model parameter in the sobol_sensitivity function
NF_KB_sensitivity <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    Tr <- 1
    
    # The differential equations
    dIKKn <- Kdeg - (Kdeg*IKKn) - (Tr*k1*IKKn)
    dIKKa <- (Tr*k1*IKKn) - (k3+Kdeg+(Tr*k2*A20))*IKKa
    dNFkBn <- (a3*IKKa)*(1-NFkBn)*(delta/(IkBa+delta)) - (i1a*IkBa)*(NFkBn/(NFkBn+epsilon))
    dA20 <- (Cdeg*NFkBn) - (Cdeg*A20)
    dIkBa <- (C4a*IkBat) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa+delta))) - (i1a*IkBa*(NFkBn/(NFkBn+epsilon)))
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

# Initial values
x0 <- c(IKKn = 1, 
        IKKa = 0, 
        NFkBn = 0,
        A20 = 0, 
        IkBa = 0,
        IkBat = 0)

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
times <- seq(0.01, 500, by = 1) 

# Doing the sensitivity analysis
sensitive_sobol_final <- sobol_sensitivity(NF_KB_sensitivity,
                                           bound_var,
                                           x0,
                                           bound_min_var,
                                           bound_max_var,
                                           times)

# Define a color palette (e.g., using 13 colors to match the number of parameters)
# color_palette <- brewer.pal(n = 11, name = "Spectral")  # For up to 11 colors
# color_palette <- c(color_palette, "darkgreen", "hotpink")  # Adding extra distinct colors if needed

color_palette <- c("red", "blue", "green", "purple", "darkorange", "black", "hotpink", "cyan", "darkgreen","darkred", "darkcyan", "grey", "coral1")

# Plotting the sensitivity analysis for all
# plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "all", main_title = "All sensitivity - SOBOL", type = "l", lwd = 3)

#For IKKn 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IKKn", main_title = "IKKn sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
#For IKKa 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IKKa", main_title = "IKKa sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
#For NFkBn
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "NFkBn", main_title = "NFkBn sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
#For A20 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "A20", main_title = "A20 sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
#For IkBa 
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IkBa", main_title = "IkBa sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
#For IkBat
plot(sensitive_sobol_final, pars_plot = bound_var[1:13], state_plot = "IkBat", main_title = "IkBat sensitivity - SOBOL", type = "l", lwd = 3, col = color_palette)
