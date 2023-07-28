#Load libraries
library(tidyverse)
library(ggpubr)
library(reshape2)


#Load (RT-)qPCR data from TaqMan Array Cards
stds_1 = read_delim("./data/raw_data/stds/2023-07-03_AL_STDC_1.txt")
stds_2 = read_delim("./data/raw_data/stds/2023-07-03_AL_STDC_2.txt")
stds_3 = read_delim("./data/raw_data/stds/2023-07-03_AL_STDC_3.txt")

#Load dPCR data from QIACuity 
starting_quant = read_csv("./data/raw_data/stds/targets_starting_quant.csv")


#Starting Quantities
##"Quantity" refers to the quantity of target within the template in copies/uL

quantity = starting_quant %>%
  mutate("10^6" = starting_quant) %>%
  mutate("10^5" = starting_quant/10) %>%
  mutate("10^4" = starting_quant/100) %>%
  mutate("10^3" = starting_quant/1000) %>%
  mutate("10^2" = starting_quant/10000) %>%
  mutate("10" =   starting_quant/100000) %>%
  mutate("1" =    starting_quant/1000000) %>%
  drop_na() %>%
  select(-plasmid, -starting_quant) %>%
  pivot_longer(!`Target Name`, names_to = "Sample Name", values_to = "quantity")



#Standard curve from Machine 1 
stds_1 = stds_1 %>%
  select('Sample Name', 'Target Name', 'CT') %>%
  left_join(quantity, by = c("Target Name","Sample Name")) %>%
  drop_na() %>%
  mutate(CT = as.numeric(CT)) %>%
  mutate(log_quant = log10(quantity)) %>%
  mutate(machine = "one")

names(stds_1) = c("sample", "target", "ct", "quant", "log_quant", "machine")


#Standard curve from Machine 2 
stds_2 = stds_2 %>%
  select('Sample Name', 'Target Name', 'CT') %>%
  left_join(quantity, by = c("Target Name","Sample Name")) %>%
  drop_na() %>%
  mutate(CT = as.numeric(CT)) %>%
  mutate(log_quant = log10(quantity)) %>%
  mutate(machine = "two")

names(stds_2) = c("sample", "target", "ct", "quant", "log_quant", "machine")


#Standard curve from Machine 3
stds_3 = stds_3 %>%
  select('Sample Name', 'Target Name', 'CT') %>%
  left_join(quantity, by = c("Target Name","Sample Name")) %>%
  drop_na() %>%
  mutate(CT = as.numeric(CT)) %>%
  mutate(log_quant = log10(quantity)) %>%
  mutate(machine = "three")

names(stds_3) = c("sample", "target", "ct", "quant", "log_quant", "machine")


#Combine standard curves 
stds = rbind(stds_1, stds_2, stds_3)

stds = stds %>%
  drop_na() %>%
  filter(target != "ascaris_lumbricoides") #This assay is not functioning


#Calculate the standard curve(s)

##Create dataframes for target and machine variables
targets = stds %>% select(target) %>% unique() 
machines = stds %>% select(machine) %>% unique()

#If you want to calculate based off of all three curves

##Set up empty dataframe
results = data.frame(matrix(ncol = 5, nrow = 0))

##Provide column names
colnames(results) = c('machine', 'target', 'slope', 'intercept', 'r.squared')


for(i in 1:nrow(targets)){
  results[i,1] = "all"
  results[i,2] = targets[i,1]
  results[i,3] = summary(lm(ct ~ log_quant, na.action=na.exclude, stds[stds$target == unlist(targets[i,]),]))$coefficients[2,1]
  results[i,4] = summary(lm(ct ~ log_quant, na.action=na.exclude, stds[stds$target == unlist(targets[i,]),]))$coefficients[1,1]
  results[i,5] = summary(lm(ct ~ log_quant, na.action=na.exclude, stds[stds$target == unlist(targets[i,]),]))$r.squared
}

  


#If you want one curve per machine

##Set up empty dataframe
results = data.frame(matrix(ncol = 5, nrow = 0))

##Provide column names
colnames(results) = c('machine', 'target', 'slope', 'intercept', 'r.squared')


for(i in 1:nrow(machines)){
  for(j in 1:nrow(targets)){
    machine = machines[i,1]
    target = targets[j,1]
    slope = summary(lm(ct ~ log_quant, 
                       na.action=na.exclude,  
                       stds[stds$target == unlist(targets[j,]) & stds$machine == unlist(machines[i,]),]))$coefficients[2,1]
    intercept = summary(lm(ct ~ log_quant, 
                       na.action=na.exclude,  
                       stds[stds$target == unlist(targets[j,]) & stds$machine == unlist(machines[i,]),]))$coefficients[1,1]
    r.squared = summary(lm(ct ~ log_quant, 
                           na.action=na.exclude,  
                           stds[stds$target == unlist(targets[j,]) & stds$machine == unlist(machines[i,]),]))$r.squared
  output = data.frame(machine, target, slope, intercept, r.squared)
  results = rbind(results, output)
    }
}


#Calculate the PCR efficiency
results = results %>% 
  mutate(efficiency = 10^-(1/slope))

#Add in the number of points used
results %>%
  left_join(
    stds %>%
      group_by(machine, target) %>%
      tally())

#Save the dataframe

write_csv(results, "./data/processed_data/standard_curves_cleaned.csv")
