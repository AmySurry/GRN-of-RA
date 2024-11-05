rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)

# T/F statement for activating the pathway
activation = T

# Differential equations
NF_KB <- function(time,state,parameters){
  A <- state[1] # IKK complex
  I <- state[2] # IKBa
  P <- state[3] # NF-KB
  
  a1 <- parameters["alpha1"]
  b1 <- parameters["beta1"]
  a2 <- parameters["alpha2"]
  b2 <- parameters["beta2"]
  a3 <- parameters["alpha3"]
  b3 <- parameters["beta3"]
  b4 <- parameters["beta4"]
  
  if(activation == T){
    R <- 1 # Receptor activation
  }else{
    R <- 0
  }
  
  # The differential equations
  dA <- R*b1-a1*A
  dI <- A*b3-a2*I
  dP <- A*b2-b4*I-a3*P
  
  return(list(c(dA,dI,dP)))
}


# Solving the differential equations
# Initial condition 
x0 <- c(A = 0, I = 0, P = 0) 

# Parameter values from Jaruszewicz-Błońska et al. (2023)
alpha1 <- 0.00145 # s^-1 Spontaneous inactivation of IKK complex
alpha2 <- 0.0763 # s^-1 IκBα degradation induced by IKK complex
alpha3 <- 0.000372+0.000106 + 0.0428  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)

beta1 <- 0.00195 # s^-1 Activation of IKK complex
beta2 <- 0.0946 # s^-1 IKK (IκBα|NF-κB) association
beta3 <- 0.0946 # s^-1 IKK (IκBα|NF-κB) association
beta4 <- 0.000595 # s^-1 IκBα nuclear import leading to... (repression of NF-κB)

# Parameter vector
parameters <- c(alpha1 = alpha1, 
                alpha2 = alpha2, 
                alpha3 = alpha3,
                beta1 = beta1, 
                beta2 = beta2, 
                beta3 = beta3, 
                beta4 = beta4)

# Time sequence
times <- seq(0, 1500, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB)


# Plotting the results
ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = A, color = "IKK complex")) +
  geom_line(aes(y = I, color = "IKBa")) +
  geom_line(aes(y = P, color = "NF-KB")) +
  labs( x = "Time", y = "Concentration") +
  scale_color_manual(name = " ", 
                     values = c("IKK complex" = "blue", 
                                            "IKBa" = "orange", 
                                            "NF-KB" = "green")) +
  theme_minimal()

# Zoomin in on the first 50 seconds
ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = A, color = "IKK complex")) +
  geom_line(aes(y = I, color = "IKBa")) +
  geom_line(aes(y = P, color = "NF-KB")) +
  labs( x = "Time", y = "Concentration") +
  scale_color_manual(name = " ", 
                     values = c("IKK complex" = "blue", 
                                "IKBa" = "orange", 
                                "NF-KB" = "green")) +
  xlim(0, 50) +
  ylim(0,0.125) +
  theme_minimal()

