rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)
library(gridExtra)

# T/F statement for activating the pathway
activation = T

# Activation time period - can be a single value or a vector
t_start <- 3600   # Start time for activation
t_end <-  32400     # End time for activation

# t_start <- c(0,5400, 10800)   # Start time for activation
# t_end <- c(2700,8100,13500)     # End time for activation


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
  Kprod <- parameters["Kprod"]
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
  dIkBa <- (C4a*IkBat) - (C5a*IkBa) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa+delta))) - (i1a*IkBa)*(NFkBn/(NFkBn+epsilon))
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat) #- (C4a*IkBat)
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat)
  
  return(list(c(dIKKn,dIKKa,dNFkBn,dA20,dIkBa,dIkBat)))
}


#### Solving the differential equations ####
# Initial condition 
# x0 <- c(A = 0, I = 0, P = 0) 
x0 <- c(IKKn = 1, 
        IKKa = 0, 
        NFkBn = 0,
        A20 = 0, 
        IkBa = 0,
        IkBat = 0)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters <- c(k3 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                a2 = 0.0763, # s^-1 IκBα degradation induced by IKK complex 
                #alpha3 = 0.000372 + 0.000106 + 0.0428,  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)
                k1 = 0.00195, # s^-1 Activation of IKK complex
                a3 = 0.0946, # s^-1 IKK (IκBα|NF-κB) association
                i1a = 0.000595, # s^-1 IκBα nuclear import leading to... (repression of NF-κB)
                Kdeg = 0.000125,
                Kprod = 0.000025,
                k2 = 0.0357,
                delta = 0.108,
                epsilon = 0.0428,
                Cdeg = 0.000106,
                C4a =0.00313,
                C5a = 0.0000578,
                C3a = 0.000372)


# Time sequence
times <- seq(0, 32400, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters, maxsteps=1000000) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting



#### Plotting the results ####
p_NFkBn <- ggplot(res_NF_KB, aes(x = time/3600, y = NFkBn)) +
  geom_line(linewidth = 1, col = "blue") +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme_minimal()


trajectory <- ggplot(res_NF_KB, aes(x = NFkBn, y = IkBa)) +
  geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

grid.arrange(p_NFkBn,trajectory, nrow = 1)
