rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)

# T/F statement for activating the pathway
activation = T

# Activation time period - can be a single value or a vector
t_start <- 0   # Start time for activation
t_end <- 1000     # End time for activation

# t_start <- c(500,2000)   # Start time for activation
# t_end <- c(1000,3000)     # End time for activation


#### Creating the model for the NF-KB pathway ####
NF_KB <- function(time,state,parameters){
  
  # Defining the genes included in the model
  A <- state[1] # IKK complex
  # KB <- state[2] # NF-KB & IkBa complex
  P <- state[2] # NF-kB
  I <- state[3] # IkBa
  IN <- state[4] # IkBa nucleus
  A20 <- state[5] # A20

  
  # Defining the parameters
  b1 <- parameters["beta1"]
  b2 <- parameters["beta2"]
  b3 <- parameters["beta3"]
  b4 <- parameters["beta4"]
  b5 <- parameters["beta5"]
  b6 <- parameters["beta6"]
  b7 <- parameters["beta7"]
  a1 <- parameters["alpha1"]
  a2 <- parameters["alpha2"]
  a3 <- parameters["alpha3"]
  a4 <- parameters["alpha4"]
  a5 <- parameters["alpha5"]

  # Activation of the pathway (on or off)
  R <- 0 # Receptor inactivation
  
  if(length(t_start) > 1){ # If there are multiple activation periods
    for(i in 1:length(t_start)){
      if(activation && time >= t_start[i] && time <= t_end[i]){
        R <- 1 # Receptor activation
        break
      }
    }
  }else{ # If there is only one activation period
    if(activation && time >= t_start && time <= t_end){
      R <- 1 # Receptor activation
    }else{
      R <- 0 # Receptor inactivation
    }
  }
  
  K <- 1
  n <- 1
  
  # The differential equations
  dA <- R*b1-b2*A20-A*a1
  dI <- b3/(1+(I/K)^n)-a2*A
  dP <- b4*A-b5*IN-a3*P
  dIN <- b6*P-a4*IN
  dA20 <- b7*P-a5*A20

  return(list(c(dA, dI, dP, dIN, dA20)))
}


#### Solving the differential equations ####
# Initial condition 
# x0 <- c(A = 0, I = 0, P = 0) 
x0 <- c(A = 0.02, # µM (Unknown)
        I = 0.02, # µM (Unknown)
        P = 0.02, # µM (found in the literature - Hoffmann et al. (2002))
        IN = 0.02, # µM (Unknown)
        A20 = 0.02) # µM (Unknown)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters <- c(beta1 = 0.00195, # s^-1 Activation of IKK complex
                beta2 = 0.03557, # Inactivation of IKKa induced by A20
                beta3 = 0.00313, # s^-1 IκBα translation
                beta4 = 0.108, # s^-1 s^-1 Michaelis-type constant for NF-κB release from IκBα
                beta5 = 0.0428, # s^-1 Michaelis-Menten type constant (NF-κB inactivation) 
                beta6 = 0.000372, # s^-1 IκBα synthesis
                beta7 = 0.000106, # s^-1 Synthesis of A20
                alpha1 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                alpha2 = 0.0763, # s^-1 IκBα degradation induced by IKK complex
                alpha3 = 0.00145, # s^-1 Spontaneous inactivation of NF-κB (unknown)
                alpha4 = 0.000372, # s^-1 nuclear IkBa degradation
                alpha5 = 0.000106)  # s^-1 A20 degradation


# Time sequence
times <- seq(0, 1000, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting


#### Plotting the results ####
plot <- ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = A, color = "IKK complex"),lwd = 1) +
  geom_line(aes(y = I, color = "IkBa"),lwd = 1) +
  geom_line(aes(y = P, color = "NF-kB"),lwd = 1) +
  geom_line(aes(y = IN, color = "Nuclear IkBa"),lwd = 1) +
  geom_line(aes(y = A20, color = "A20"),lwd = 1) +
  labs( x = "Time (sec)", y = "Concentration (µM)") +
  scale_color_manual(name = " ", 
                     values = c("IKK complex" = "blue", 
                                "IkBa" = "red",
                                "Nuclear IkBa" = "orange", 
                                "NF-kB" = "green",
                                "A20" = "black")) +
  theme_minimal() +
  # Change the font size and style of the axis and legend text
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))
  # ylim(-50, 50) # Set the y-axis limits
  # xlim(250, 1000) # Set the x-axis limits

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

