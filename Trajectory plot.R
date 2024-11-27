rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)
library(gridExtra)
library(grid)
library(ggarrow)



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
  dIkBa <- (C4a*IkBat) - (a2*IKKa*IkBa) - (a3*IKKa*(1-NFkBn)*(IkBa/(IkBa+delta))) - (i1a*IkBa)*(NFkBn/(NFkBn+epsilon))
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat) #- (C4a*IkBat)
  dIkBat <- (C3a*NFkBn) - (C3a*IkBat)
  
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
parameters_A <- c(k3 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                  a2 = 0.0762, # s^-1 IκBα degradation induced by IKK complex 
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
                  C5a = 0.000058,
                  C3a = 0.000372)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters_B <- c(k3 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                  a2 = 0.02, # s^-1 IκBα degradation induced by IKK complex 
                  #alpha3 = 0.000372 + 0.000106 + 0.0428,  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)
                  k1 = 0.00195, # s^-1 Activation of IKK complex
                  a3 = 0.0946, # s^-1 IKK (IκBα|NF-κB) association
                  i1a = 0.0001, # s^-1 IκBα nuclear import leading to... (repression of NF-κB)
                  Kdeg = 0.000125,
                  Kprod = 0.000025,
                  k2 = 0.0357,
                  delta = 0.108,
                  epsilon = 0.0428,
                  Cdeg = 0.000106,
                  C4a =0.00313,
                  C5a = 0.00001,
                  C3a = 0.000372)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters_C <- c(k3 = 0.00145, # s^-1 Spontaneous inactivation of IKK complex
                a2 = 0.01, # s^-1 IκBα degradation induced by IKK complex 
                #alpha3 = 0.000372 + 0.000106 + 0.0428,  # s^-1 IκBα synthesis and degradation of IκBα mRNA + Synthesis of A20 + Michaelis-Menten type constant (NF-κB inactivation)
                k1 = 0.00195, # s^-1 Activation of IKK complex
                a3 = 0.0946, # s^-1 IKK (IκBα|NF-κB) association
                i1a = 0.0001, # s^-1 IκBα nuclear import leading to... (repression of NF-κB)
                Kdeg = 0.000125,
                Kprod = 0.000025,
                k2 = 0.0357,
                delta = 0.108,
                epsilon = 0.0428,
                Cdeg = 0.000106,
                C4a =0.00313,
                C5a = 0.00001,
                C3a = 0.000372)


# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters_9A <- c(k3 = 0.00145, 
                a2 = 0.0763, 
                k1 = 0.00195, 
                a3 = 0.0946, 
                i1a = 0.000595, 
                Kdeg = 0.000125,
                Kprod = 0.000025,
                k2 = 0, #Changed
                delta = 0.108,
                epsilon = 0.0428,
                Cdeg = 0.000106,
                C4a =0.00313,
                C5a = 0.0000578,
                C3a = 0.000372)

# Parameter values from Jaruszewicz-Błońska et al. (2023)
parameters_9B <- c(k3 = 0.00145, 
                   a2 = 0.0763, 
                   k1 = 0.00195, 
                   a3 = 0.0946, 
                   i1a = 0.01, 
                   Kdeg = 0.000125,
                   Kprod = 0.000025,
                   k2 = 3.57, #Changed
                   delta = 0.108,
                   epsilon = 0.0428,
                   Cdeg = 0.000106,
                   C4a =0.00313,
                   C5a = 0.0000578,
                   C3a = 0.000372)


# Time sequence
times_A <- seq(0, 32400, by = 1) 
times_BC <- seq(0, 108000, by = 1) 
times_9A <- seq(0, 28800, by = 1)
times_9B <- seq(0, 54000, by = 1)

# T/F statement for activating the pathway
activation = T

# Activation time period - can be a single value or a vector
t_start <- 3600   # Start time for activation
t_end <-  32400     # End time for activation (plot A)

# Solve ODEs
res_A <- ode(x0, times_A, NF_KB, parameters_A, maxsteps=1000000) 
res_A <- as.data.frame(res_A) # Convert the results to a data frame for plotting

# Activation time period - can be a single value or a vector
t_end <-  108000     # End time for activation (plot B and C)

res_B <- ode(x0, times_BC, NF_KB, parameters_B, maxsteps=1000000) 
res_B <- as.data.frame(res_B) # Convert the results to a data frame for plotting

res_C <- ode(x0, times_BC, NF_KB, parameters_C, maxsteps=1000000) 
res_C <- as.data.frame(res_C) # Convert the results to a data frame for plotting


# Activation time period - can be a single value or a vector
t_end <-  28800     # End time for activation (plot B and C)

res_9A <- ode(x0, times_9A, NF_KB, parameters_9A, maxsteps=1000000) 
res_9A <- as.data.frame(res_9A) 

# Activation time period - can be a single value or a vector
t_end <-  54000     # End time for activation (plot B and C)

res_9B <- ode(x0, times_9B, NF_KB, parameters_9B, maxsteps=1000000) 
res_9B <- as.data.frame(res_9B) 


#### Plotting the results ####
A_NFkBn <- ggplot(res_A, aes(x = time/3600)) +
  geom_line(aes(y = NFkBn, col = "NFkBn"),linewidth = 1) +
  geom_line(aes(y = IkBa, col = "IkBa"),linewidth = 1) +
  labs( x = "Time (hours)", y = "Concentration") +
  scale_color_manual(values = c("NFkBn" = "blue", "IkBa" = "red")) +
  theme_minimal()

trajectory_A <- ggplot(res_A, aes(x = NFkBn, y = IkBa)) +
  # geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))+
  geom_arrow(arrow_mid = "head_wings", mid_place = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), col = "red")


A <- grid.arrange(A_NFkBn,trajectory_A, nrow = 1, top = textGrob("A: a2 = 0.0762, c5a = 0.000058, i1a = 0.000595", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02))
grid.arrange(A_NFkBn,trajectory_A, nrow = 1)


ggplot(res_A, aes(x = NFkBn, y = A20)) +
  # geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "A20") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))+
  geom_arrow(arrow_mid = "head_wings", mid_place = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), col = "red")


ggplot(res_A, aes(x = NFkBn, y = IKKa)) +
  # geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IKKa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))+
  geom_arrow(arrow_mid = "head_wings", mid_place = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), col = "red")


ggplot(res_A, aes(x = NFkBn, y = IkBat)) +
  # geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBat") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10)) +
  geom_arrow(arrow_mid = "head_wings", mid_place = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), col = "red")




B_NFkBn <- ggplot(res_B, aes(x = time/3600, y = NFkBn)) +
  geom_line(linewidth = 1, col = "blue") +
  labs( x = "Time (hours)", y = "NF-kBn") +
  theme_minimal()


trajectory_B <- ggplot(res_B, aes(x = NFkBn, y = IkBa)) +
  geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

B <- grid.arrange(B_NFkBn,trajectory_B, nrow = 1, top = textGrob("B: a2 =0.02, c5a = 0.00001, i1a = 0.0001", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02))


C_NFkBn <- ggplot(res_C, aes(x = time/3600, y = NFkBn)) +
  geom_line(linewidth = 1, col = "blue") +
  labs( x = "Time (hours)", y = "NF-kBn") +
  theme_minimal()


trajectory_C <- ggplot(res_C, aes(x = NFkBn, y = IkBa)) +
  geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

C <- grid.arrange(C_NFkBn,trajectory_C, nrow = 1, top = textGrob("C: a2 =0.01, c5a = 0.00001, i1a = 0.0001", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02))

grid.arrange(A,B,C, nrow = 3)




NFkBn_9A <- ggplot(res_9A, aes(x = time/3600, y = NFkBn)) +
  geom_line(linewidth = 1, col = "blue") +
  labs( x = "Time (hours)", y = "NF-kBn") +
  theme_minimal()


trajectory_9A <- ggplot(res_9A, aes(x = NFkBn, y = IkBa)) +
  geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "IkBa") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

p_9A <- grid.arrange(NFkBn_9A,trajectory_9A, nrow = 1, top = textGrob("A: k2 = 0", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02))


NFkBn_9B <- ggplot(res_9B, aes(x = time/3600, y = NFkBn)) +
  geom_line(linewidth = 1, col = "blue") +
  labs( x = "Time (hours)", y = "NF-kBn") +
  theme_minimal()


trajectory_9B <- ggplot(res_9B, aes(x = NFkBn, y = A20)) +
  geom_point(size = 0.5, col = "red") +
  labs(x = "NFkBn", y = "A20") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))

p_9B <- grid.arrange(NFkBn_9B,trajectory_9B, nrow = 1, top = textGrob("B: k2 = 3.57, i1a = 0.01 ", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02))


grid.arrange(p_9A,p_9B, nrow = 2)
