# Planning the Management of Wild Animal Populations Using a Multi-Trophic
# Mechanistic Model

rm(list=ls())

# during this practical we wiiL implement a discrete, mechanistic model for six
# distinct populations: lichen, plants, moose forage, caribou, moose and wolves. 

# lichen, plant and moose forage growth are based on a typical logistic growth

# lichen, plant and moose forage are depleted through herbivory foiLowing a
# type II functional response

# caribou and moose growth are based on type II functional response, which
# depends on forage availability.

# caribou and moose die when there is not enough food to sustain the whole
# population. 

# wolf growth is based on type II functional response, which depends on caribou
# and moose availability. 


# STEP ONE: BASE MODEL

# parameter values
sP <- 0.5 #plant annual mas growth potential
Kp <- 240 #plant carrying capacity (kg/ha)
lP <- 2.68275 # annual plant requirement by an individual caribou (kg)
hp <- 300 # maximum rate of plant intake (kg/ha)
sL <- 0.07 #lichen annual max growth potential
Kl <- 870 #lichen carrying capacity (kg/ha)
lL <- 2.68275 #annual lichen requirement by an individual caribou (kg)
hl <- 400 #maximum rate of lichen intake (kg/ha)
sH <- 0.5 #forage annual max growth potential
Kh <- 970 #forage carrying capacity (kg/ha)
lH <- 3.1682 #annual forage requirement by an individual moose (kg)
hh <- 1000 #maximum rate of moose forage intake (kg/ha)
fC <- 1 #maximum fecundity for caribou (average number of offspring per individual)
eC <- 1.85 #asymptotic kiiL rate on caribou per wolf (unit-less)
mC <- 0 # hunting rate on caribou
fM <- 1.5 #maximum fecundity for moose (average number of offspring per individual)
eM <- 0.6 # asymptotic kiiL rate on moose (animal/ha)
mM <- 0 #hunting rate on moose
mW <- 0 #hunting rate on wolf
b <- 0.8 #upper limit of annual rate of increase
g <- 0.38 #wolf death rate
dC <- 460 #depredation efficiency on caribou (animal/ha; prey density at half the asymptotic kiiL rate)
dM <- 46 #predation efficiency on moose (animal/ha; prey density at half the asymptotic kiiL rate)

# generate data frame
n_steps <- 200
population_initial <- data.frame(time=1:n_steps,P=numeric(n_steps),L=numeric(n_steps),H=numeric(n_steps),C=numeric(n_steps),M=numeric(n_steps),W=numeric(n_steps))
# add initial population data
population_initial <- within(population_initial,{
  P[1] <- 240 #plant biomass density (kg/ha)
  L[1] <- 870 #lichen biomas density (kg/ha)
  H[1] <- 970 #moose forage biomass density (kg/ha)
  C[1] <- 7 #caribou density (animals/ha)
  M[1] <- 25 #moose density (animals/ha)
  W[1] <- 8 #wolf density (animals/ha)
})


# calculate populations
for(t in 2:n_steps){
  population_initial <- within(population_initial,{
    #Primary producers - plants
    plant_birth <- sP * P[t-1] * (1-(P[t-1]/Kp))
    plant_death <- (C[t-1] * lP * P[t-1])/(hp + P[t-1])
    P[t] <- max(0,P[t-1] + plant_birth - plant_death) #plants consumed by Caribou - the max is used to avoid negative values
    
    #Primary producers - lichen
    lichen_birth <- sL * L[t-1] * (1-(L[t-1]/Kl))
    lichen_death <- (C[t-1] * lL * L[t-1])/(hl + L[t-1])
    L[t] <- max(0,L[t-1] + lichen_birth - lichen_death) #lichen consumed by Caribou
    
    forage_birth <- sH * H[t-1] * (1-(H[t-1])/Kh)
    forage_death <- H[t-1] * lH * M[t-1]/(hh + H[t-1])
    H[t] <- max(0,H[t-1] + forage_birth - forage_death) #plants consumed by Moose
    
    # First trophic level - Caribou
    caribou_growth <- C[t-1] * fC * (L[t-1]/(L[t-1]+hl)) * (P[t-1]/(P[t-1]+hp))
    caribou_death <- C[t-1] * (1 - P[t-1]/(P[t-1]+hp)) * (1 - L[t-1]/(L[t-1]+hl))
    caribou_predation <- (W[t-1]* eC * C[t-1])/(C[t-1] + dC) 
    C[t] <- max(0,C[t-1] + caribou_growth - caribou_death - caribou_predation)
    
    # First trophic level - Moose
    moose_growth <- M[t-1] * fM * H[t-1]/(H[t-1] + hh)
    moose_death <- M[t-1] * (1 - H[t-1]/(H[t-1] + hh))
    moose_predation <-  (W[t-1] * eM * M[t-1])/(M[t-1] + dM)
    M[t] <- max(0,M[t-1] + moose_growth - moose_death - moose_predation)
    
    # Third trophic level - Wolf
    wolf_growth <- W[t-1] * b * ((C[t-1]/(C[t-1] + dC)) + (M[t-1]/(M[t-1] + dM)))
    wolf_death <- W[t-1] * g
    W[t] <- max(0,W[t-1] + wolf_growth - wolf_death)
  
    # Predator
    # W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    # W.death <- g*W[t-1]
    # W[t] <- max(0,W[t-1] + W.growth - W.death)
      
  })
}

# plot outcome of calculations
library(ggplot2)
colours <- c("Plants"="darkgreen","Lichen"="green","Shrubs"="lightgreen","Caribou"="lightblue","Moose"="blue","Wolf"="maroon")
ggplot(data = population_initial)+
  geom_line(mapping=aes(x=time,y=P,color="Plants"), linetype = "dashed") +
  geom_line(mapping=aes(x=time,y=L,color="Lichen"), linetype = "dashed") +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs"), linetype = "dashed") +
  geom_line(mapping=aes(x=time,y=C,color="Caribou"), linetype = "dashed") +
  geom_line(mapping=aes(x=time,y=M,color="Moose"), linetype = "dashed") +
  geom_line(mapping=aes(x=time,y=W,color="Wolf"), linetype = "dashed") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Original model")+
  scale_color_manual(values = colours) +
  theme_bw()


# Q1. Now set the initial populations to 0 and run the model. What happens?
population_initial <- within(population_initial,{
  P[1] <- 240 #plant biomass density (kg/ha)
  L[1] <- 870 #lichen biomas density (kg/ha)
  H[1] <- 970 #moose forage biomass density (kg/ha)
  C[1] <- 0 #caribou density (animals/ha)
  M[1] <- 0 #moose density (animals/ha)
  W[1] <- 0 #wolf density (animals/ha)
})
#     All plant species remain at a constant level, the animal species never
#     have any individuals.

#Q2. Set the caribou and moose populations to the values in the table. How are
#    the two species doing?
population_initial <- within(population_initial,{
  P[1] <- 240 #plant biomass density (kg/ha)
  L[1] <- 870 #lichen biomas density (kg/ha)
  H[1] <- 970 #moose forage biomass density (kg/ha)
  C[1] <- 7 #caribou density (animals/ha)
  M[1] <- 25 #moose density (animals/ha)
  W[1] <- 0 #wolf density (animals/ha)
})
#    Plants remain consistent throughout, at a lower density than the other
#    producers. Lichen show large oscillation. Shrubs show a little decline and
#    oscillation before settling on a constant level. Both moose and caribou
#    show an initial increase, before levelling off. Caribou a bit less stable
#    than moose. More moose are supported than Caribou. 

# Q3. Now set the wolf initial population to the one given in the table. What
#     happens to the moose and caribou populations. Which one is more impacted?
population_initial <- within(population_initial,{
  P[1] <- 240 #plant biomass density (kg/ha)
  L[1] <- 870 #lichen biomas density (kg/ha)
  H[1] <- 970 #moose forage biomass density (kg/ha)
  C[1] <- 7 #caribou density (animals/ha)
  M[1] <- 25 #moose density (animals/ha)
  W[1] <- 8 #wolf density (animals/ha)
})
#     Similar k for each of the species. But grass now also shows oscillations.
#     Very little change in moose and caribou, because wolves very quickly
#     decline to extinction. But Caribou o
