rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)

# T/F statement for activating the parthway
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
    R <- 1
  }else{
    R <- 0
  }
  
  dA <- R*b1-a1*A
  dI <- A*b3-a2*I
  dP <- A*b2-b4*I-a3*P
  
  return(list(c(dA,dI,dP)))
}


# Solving the differential equations
# Initial condition 
x0 <- c(A = 0, I = 0, P = 0) 

beta1 <- log(2)
beta2 <- log(2)
beta3 <- log(2)
beta4 <- log(2)
alpha1 <- 0.5
alpha2 <- 0.5
alpha3 <- 0.5

# Adjust parameters to allow comparison
parameters <- c(alpha1 = alpha1, 
                alpha2 = alpha2, 
                alpha3 = alpha3,
                beta1 = beta1, 
                beta2 = beta2, 
                beta3 = beta3, 
                beta4 = beta4)

# Time sequence
times <- seq(0, 100, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB)


# Plotting the results
ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = A, color = "A")) +
  geom_line(aes(y = I, color = "I")) +
  geom_line(aes(y = P, color = "P")) +
  labs( x = "Time", y = "Concentration") +
  scale_color_manual(name = " ", values = c("A" = "blue", "I" = "orange", "P" = "green")) +
  theme_minimal()


