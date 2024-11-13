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

t_start <- c(500,2000)   # Start time for activation
t_end <- c(1000,3000)     # End time for activation


#### Creating the model for the NF-KB pathway ####
NF_KB <- function(time,state,parameters){
  
  # Defining the genes included in the model
  A <- state[1] # IKK complex
  I <- state[2] # IkBa
  P <- state[3] # NF-kB
  
  # Defining the parameters
  a1 <- parameters["alpha1"]
  b1 <- parameters["beta1"]
  a2 <- parameters["alpha2"]
  b2 <- parameters["beta2"]
  a3 <- parameters["alpha3"]
  b3 <- parameters["beta3"]
  
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
  
  # The differential equations
  dA <- R*b1-a1*A
  dI <- A*b2-a2*I
  dP <- A*b2-b3*I-a3*P
  
  return(list(c(dA,dI,dP)))
}


#### Solving the differential equations ####
# Initial condition 
# x0 <- c(A = 0, I = 0, P = 0) 
x0 <- c(A = 0.02, # µM (Unknown)
        I = 0.02, # µM (Unknown)
        P = 0.02) # µM (found in the literature - Hoffmann et al. (2002))

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters <- c(alpha1 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                alpha2 = 0.0763, # s^-1 IκBα degradation induced by IKK complex 
                alpha3 = 0.000372 + 0.000106 + 0.0428,  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)
                beta1 = 0.00195, # s^-1 Activation of IKK complex
                beta2 = 0.0946, # s^-1 IKK (IκBα|NF-κB) association
                beta3 = 0.000595) # s^-1 IκBα nuclear import leading to... (repression of NF-κB)

# Time sequence
times <- seq(0, 5000, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting


#### Plotting the results ####
plot <- ggplot(res_NF_KB, aes(x = time)) +
  geom_line(aes(y = A, color = "IKK complex"),lwd = 1) +
  geom_line(aes(y = I, color = "IkBa"),lwd = 1) +
  geom_line(aes(y = P, color = "NF-kB"),lwd = 1) +
  labs( x = "Time (sec)", y = "Concentration (µM)") +
  scale_color_manual(name = " ", 
                     values = c("IKK complex" = "blue", 
                                "IkBa" = "orange", 
                                "NF-kB" = "green")) +
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
