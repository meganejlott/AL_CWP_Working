#LOAD LIBRARIES
library(tidyverse)
library(magrittr)
library(purrr)



#LOAD DATA

##Standard Curves
standard_curves = read_csv("./data/processed_data/standard_curves_cleaned.csv")


##Meta data
sample_key = read_csv("./data/raw_data/sample_key_1.csv")
processing_key = read_csv("./data/raw_data/processing_key_1.csv")
target_key = read_csv("./data/raw_data/target_key_1.csv")


##Raw Data
###The code below compiles all of the txt files into one dataframe

#read in file names for TAC Data
raw_files = fs::dir_ls("./data/raw_data/TAC_Raw", glob="*.txt")

#Read in the data from the txt files and combine into one dataframe.
TAC_data = raw_files %>% 
  purrr::map_dfr(read_delim, .id = "source", col_names = TRUE) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "./TAC_Raw", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, ".txt", "")) %>% 
  separate(source, into = c("drop", "card_num"), sep = "dat/AL_") %>%
  select("card_num", "Sample Name", "Target Name", "CT") 

names(TAC_data) = c("card_num", "sample_name", "target", "ct")



#Join in the meta data
meta_data = left_join(sample_key, processing_key)


TAC_data %<>% 
  mutate(card_num = as.numeric(card_num)) %>%
  left_join(meta_data)


#Join in the standard curve data
TAC_data %<>% left_join(standard_curves)



#Calculate the copies per uL of template
TAC_data %<>%
  mutate(ct = as.numeric(ct), 
         intercept = as.numeric(intercept), 
         slope = as.numeric(slope)) %>%
  mutate(copies_uL_template = 10^((ct-intercept)/slope))




col_order = c("sample_name", "sample_type", 
               "target", "ct", "copies_uL_template",
               "slope", "intercept", "r.squared", "efficiency", 
               "card_num", "position", "machine", "tac_date",
               "tac_label", "dummy_code", "prep_num", 
               "extraction_date", "extraction_group", "extraction_kit", 
               "qubit_date", "qubit_group", "qubit_dna", "qubit_rna")
TAC_data = TAC_data[, col_order]




#SAVE PROCESSED DATA 
write_csv(TAC_data, "./data/processed_data/TAC_data_cleaned.csv")
