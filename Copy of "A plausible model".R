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
t_end <- 18000     # End time for activation

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

  return(list(c(dIKKn,dIKKa,dNFkBn,dA20,dIkBa,dIkBat)))
}


#### Solving the differential equations ####
# Initial condition 
x0 <- c(IKKn = 1, 
        IKKa = 0, 
        NFkBn = 0,
        A20 = 0, 
        IkBa = 0,
        IkBat = 0)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters <- c(k3 = 0.00145,
                a2 = 0.0763,
                k1 = 0.00195,
                a3 = 0.0946,
                i1a = 0.000595, 
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
times <- seq(0, 18000, by = 1) 

# Solve ODEs
res_NF_KB <- ode(x0, times, NF_KB, parameters, maxsteps=1000000) 
res_NF_KB <- as.data.frame(res_NF_KB) # Convert the results to a data frame for plotting



#### Plotting the results ####
plot <- ggplot(res_NF_KB, aes(x = time / 3600)) +  # Converting seconds to hours
  geom_line(aes(y = IKKn, color = "IKKn"), lwd = 1) +
  geom_line(aes(y = IKKa, color = "IKKa"), lwd = 1) +
  geom_line(aes(y = NFkBn, color = "NFkBn"), lwd = 1) +
  geom_line(aes(y = A20, color = "A20"), lwd = 1) +
  geom_line(aes(y = IkBa, color = "IkBa"), lwd = 1) +
  geom_line(aes(y = IkBat, color = "IkBat"), lwd = 1) +
  labs(x = "Time (hours)", y = "Concentration (µM)") +
  scale_color_manual(name = " ", 
                     values = c("IKKn" = "blue", 
                                "IKKa" = "orange", 
                                "NFkBn" = "green",
                                "A20" = "pink",
                                "IkBa" = "gold1",
                                "IkBat" = "red")) +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

# Add markers for activation and relaxation
if (t_start[1] != res_NF_KB$time[1] || t_end[length(t_end)] != res_NF_KB$time[length(res_NF_KB$time)]) {
  if (activation) {
    plot <- plot + 
      geom_vline(xintercept = t_start / 3600, linetype = "dotted", color = "darkgreen") + # Converting t_start to hours
      geom_vline(xintercept = t_end / 3600, linetype = "dotted", color = "darkred") + # Converting t_end to hours
      annotate("text", x = t_start / 3600 + 0.1, 
               y = max(res_NF_KB$IkBa, res_NF_KB$A20, res_NF_KB$NFkBn, res_NF_KB$IKKa, res_NF_KB$IKKn), 
               label = "Activation", color = "darkgreen", vjust = -0.5, hjust = 0) +
      annotate("text", x = t_end / 3600 + 0.1, 
               y = max(res_NF_KB$IkBa, res_NF_KB$A20, res_NF_KB$NFkBn, res_NF_KB$IKKa, res_NF_KB$IKKn), 
               label = "Relaxation", color = "darkred", vjust = -0.5, hjust = 0)
  }
}

# Show the plot
print(plot)



p_IKKa <- ggplot(res_NF_KB, aes(x = time/3600, y = IKKa)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa")

p_NFkBn <- ggplot(res_NF_KB, aes(x = time/3600, y = NFkBn)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn")

p_IkBa <- ggplot(res_NF_KB, aes(x = time/3600, y = IkBa)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa")

p_A20 <- ggplot(res_NF_KB, aes(x = time/3600, y = A20)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20")

p_IkBat <- ggplot(res_NF_KB, aes(x = time/3600, y = IkBat)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat")

grid.arrange(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1)
