# source init file
library(here)
source(here("init.R"))
source("simulation_function.R")

n = 41:100 # uniform from 7

for (i in n) {
  options(timeout = 100000000)
  print(parameters(index=i))
  my_simulation(parameters = parameters(index = i),
                base_climate_map = "maxent")
}

# intro_suitability = map_data$values |> unique()
# Ecological attributes and hypotheses include enemy release, novel weapons, empty niches, the evolution of increased competitive ability, and others (Rai 2015). 