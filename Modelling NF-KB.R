rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)

# T/F statement for activating the pathway
activation = T

# Activation time period - can be a single value or a vector
# t_start <- 0   # Start time for activation
# t_end <- 5000     # End time for activation

t_start <- c(0,14400)   # Start time for activation
t_end <- c(14400,40000)     # End time for activation


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
  Tr <- 1 # Receptor inactivation
  
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
  dIkBa <- (C4a*IkBat) - (C5a*IkBa) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa-delta))) - (i1a*IkBa*(NFkBn/(NFkBn+epsilon)))
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat)
  
  return(list(c(dIKKn,dIKKa,dNFkBn,dA20,dIkBa,dIkBat)))
}


#### Solving the differential equations ####
# Initial condition 
# x0 <- c(A = 0, I = 0, P = 0) 
x0 <- c(IKKn = 0.02, # µM (Unknown)
        IKKa = 0.02, # µM (Unknown)
        NFkBn = 0.02, # µM (found in the literature - Hoffmann et al. (2002))
        A20 = 0.02, 
        IkBa = 0.02,
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
times <- seq(0, 14400, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting


#### Plotting the results ####
plot <- ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = IKKn, color = "IKKn"),lwd = 1) +
  geom_line(aes(y = IKKa, color = "IKKa"),lwd = 1) +
  geom_line(aes(y = NFkBn, color = "NFkB"),lwd = 1) +
  geom_line(aes(y = A20, color = "A20"),lwd = 1) +
  geom_line(aes(y = IkBa, color = "IkBa"),lwd = 1) +
  geom_line(aes(y = IkBat, color = "IkBat"),lwd = 1) +
  labs( x = "Time (sec)", y = "Concentration (µM)") +
  scale_color_manual(name = " ", 
                     values = c("IKKn" = "blue", 
                                "IKKa" = "orange", 
                                "NFkB" = "green",
                                "A20" = "Pink",
                                "IkBa" = "Yellow",
                                "IkBat" = "Black")) +
  theme_minimal() +
  # Change the font size and style of the axis and legend text
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

# Add markers for activation and relaxation
if(t_start[1] !=res_NF_KB$time[1] && t_end[length(t_end)] != res_NF_KB$time[length(res_NF_KB$time)]){
  if(activation){
    plot <- plot + 
      geom_vline(xintercept = t_start, linetype = "dotted", color = "darkgreen") + # Dotted line for activation
      geom_vline(xintercept = t_end, linetype = "dotted", color = "darkred") + # Dotted line for relaxation
      annotate("text", x = t_start+10, y = max(res_NF_KB$A, res_NF_KB$I, res_NF_KB$P), # Add text/label to the line for activation
               label = "Activation", color = "darkgreen", vjust = -0.5, hjust = 0) +
      annotate("text", x = t_end+10, y = max(res_NF_KB$A, res_NF_KB$I, res_NF_KB$P), # Add text/label to the line for relaxation
               label = "Relaxation", color = "darkred", vjust = -0.5, hjust = 0)
  }
}

# Show the plot
print(plot)



# # Zoomin in on the first 50 seconds
# ggplot(res_NF_KB, aes(x = time)) +
#   geom_line(aes(y = A, color = "IKK complex"), lwd = 1) +
#   geom_line(aes(y = I, color = "IKBa"), lwd = 1) +
#   geom_line(aes(y = P, color = "NF-KB"), lwd = 1) +
#   labs( x = "Time", y = "Concentration") +
#   scale_color_manual(name = " ",
#                      values = c("IKK complex" = "blue",
#                                 "IKBa" = "orange",
#                                 "NF-KB" = "green")) +
#   xlim(0, 50) +
#   ylim(0, 0.4) +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 12, face = "bold"),
#         axis.text = element_text(size = 10, face = "bold"),
#         legend.text = element_text(size = 10))
