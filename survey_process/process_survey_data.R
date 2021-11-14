# Process data from the general survey 2018 data
pacman::p_load(tidyverse, sjlabelled)


df <- read_spss("GSS2018.sav")


outcome_variables <-  c("ABINSPAY", "ABMEDGOV1", "ABHELP1", "ABHELP2", "ABHELP3", "ABHELP4", 
                        "ABMORAL", "ABSTATE1", "ABSTATE2")

demogrpahic_variables <- c("SEX", "RACE", "CLASS", "PARTYID")

sample <- df %>% select(outcome_variables, demogrpahic_variables, ID) 

# print questions
sample %>% get_label()

# print value codes
sample %>% get_label()

# recode "no" "yes" to = 0, missing as "1"
model_data_bernoulli <- sample %>% 
  mutate_at(outcome_variables, 
            list(~recode(., '1'=0, '2'=0, '3'=0, '4'=0, '5'=0))) %>% 
  mutate_at(outcome_variables, 
            list( ~replace_na(., "1"))) %>% 
  drop_na() 

# print counts & proportion for each item
map(outcome_variables, function(x){
  model_data_bernoulli %>% count_(x) %>% 
    mutate(proportion = n/sum(n))
})

# save data with demographics
write_csv(model_data_bernoulli, "model_data.csv")

