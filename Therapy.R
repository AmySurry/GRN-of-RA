rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(ODEsensitivity)
library(gridExtra)
library(grid)

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


#### Increased IKKa degradation ####
# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["k3"] <- parameters["k3"] * (1 + therapy_factor)  # Increase IKKa degradation
  parameters_therapy["k2"] <- parameters["k2"] * (1 + therapy_factor)  # Increase IKKa degradation

  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
IKKa_p <- ggplot(combined_results, aes(x = time / 3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs(x = "Time (hours)", y = "Concentration (µM)", title = "Effect of Therapy on IKKa Degradation") +
  theme_minimal() +
  theme(
    legend.position = "bottom",       # Move legend to the bottom
    legend.direction = "horizontal"     # Arrange legend items vertically
  ) 



# Function to extract the legend from a ggplot
extract_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

# Extract the shared legend
shared_legend <- extract_legend(IKKa_p)

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on IKKa Degradation", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
  )



#### Increase IkBa degradation ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["C5a"] <- parameters["C5a"] * (1 + therapy_factor)  # Increase IKKa degradation
  parameters_therapy["a2"] <- parameters["a2"] * (1 + therapy_factor)  # Increase IKKa degradation
  
  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on IkBa Degradation", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)




#### Decrease IkBa degradation ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["C5a"] <- parameters["C5a"] * (1 - therapy_factor)  # Decrease IKKa degradation
  parameters_therapy["a2"] <- parameters["a2"] * (1 - therapy_factor)  # Decrease IKKa degradation
  
  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on IkBa Degradation", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)



#### Decreased IkBat translation ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["C4a"] <- parameters["C4a"] * (1-therapy_factor)  # Decrease IKKa degradation

  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on IkBat Translation", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)



#### Decreased NFkBn release from complex ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["a3"] <- parameters["a3"] * (1-therapy_factor)  # Increase IKKa degradation
  
  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on NFkBn release from complex", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)





#### Increase NFkBn and IkBa complex formation ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["i1a"] <- parameters["i1a"] * (1+therapy_factor)  # Increase IKKa degradation
  
  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on NFkBn and IkBa complex formation", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)






#### Inhibit NFkBn binding to DNA ####

# Therapy factors: 0% (baseline), 25%, 50%, 75%, 100%
therapy_levels <- c(0, 0.25, 0.5, 0.75, 1.0)

# Store results in a list
therapy_results <- list()

# Loop through each therapy level
for (therapy_factor in therapy_levels) {
  
  # Adjusted parameters for therapy
  parameters_therapy <- parameters
  parameters_therapy["C3a"] <- parameters["C3a"] * (1-therapy_factor)  # Increase IKKa degradation
  parameters_therapy["Cdeg"] <- parameters["Cdeg"] * (1-therapy_factor)  # Increase IKKa degradation
  
  # Solve ODEs with therapy
  res_therapy <- ode(x0, times, NF_KB, parameters_therapy, maxsteps = 1000000)
  res_therapy <- as.data.frame(res_therapy)
  res_therapy$Therapy <- paste0(therapy_factor * 100, "%")  # Add therapy label
  
  # Store in results list
  therapy_results[[paste0("Therapy_", therapy_factor * 100)]] <- res_therapy
}

# Combine all results into a single data frame for plotting
combined_results <- do.call(rbind, therapy_results)

# Specify the desired order of the legend
combined_results$Therapy <- factor(combined_results$Therapy, 
                                   levels = c("100%", "75%", "50%", "25%", "0%"))

# Plot results for IKKa
p_IKKa <- ggplot(combined_results, aes(x = time/3600, y = IKKa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IKKa") +
  theme(legend.position = "none")

p_NFkBn <- ggplot(combined_results, aes(x = time/3600, y = NFkBn, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "NF-kBn") +
  theme(legend.position = "none")

p_IkBa <- ggplot(combined_results, aes(x = time/3600, y = IkBa, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBa") +
  theme(legend.position = "none")

p_A20 <- ggplot(combined_results, aes(x = time/3600, y = A20, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "A20") +
  theme(legend.position = "none")

p_IkBat <- ggplot(combined_results, aes(x = time/3600, y = IkBat, color = Therapy)) +
  geom_line() +
  labs( x = "Time (hours)", y = "Concentration (µM)", title = "IkBat") +
  theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p_IKKa, p_NFkBn, p_IkBa, p_A20, p_IkBat, nrow = 1),  # Arrange plots in one row
  shared_legend,  # Add the shared legend below
  nrow = 2,
  heights = c(10, 1),  # Adjust the relative heights (10 for plots, 1 for legend)
  top = textGrob("Effect of Therapy on NFkBn binding to DNA", gp = gpar(fontsize = 16, fontface = "bold"), hjust = 0, x=0.02)
)
