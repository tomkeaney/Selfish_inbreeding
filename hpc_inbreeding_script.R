

library(tidyverse) # for tidy style coding and plotting
library(data.table) # for efficient handling of large dataframes
#library(Rmpi) # for parallel computing on a hpc
library(snow) # complements the Rmpi package for hpc computing
cl <- makeCluster(64, type="SOCK")

# mendelian genetics function for reproduction

make_mating_table <- function(gene_location){
  
  make_offspring <- function(X, Y, type, zygote_freq, gene_location){
    tibble(Female_genotype = X,
           Male_genotype = Y,
           type,
           zygote_freq,
           locus_type = gene_location)
  }
  
  # Specify the possible offspring genotypes for all the potential crosses; we use these for the type argument in the make_offspring function
  
  # autosomal
  
  # II x II
  
  a_genotype_1 <- c("A_IA_I.Female", "A_IA_I.Male")
  
  # II x IO
  
  a_genotype_2 <- c("A_IA_I.Female", "A_IA_I.Male", 
                    "A_IA_O.Female", "A_IA_O.Male")
  
  # II x OO
  
  a_genotype_3 <- c("A_IA_O.Female", "A_IA_O.Male")
  
  # IO x IO
  
  a_genotype_4 <- c("A_IA_I.Female", "A_IA_I.Male", 
                    "A_IA_O.Female", "A_IA_O.Male", 
                    "A_OA_O.Female", "A_OA_O.Male")
  
  # IO x OO
  
  a_genotype_5 <- c("A_IA_O.Female", "A_IA_O.Male",
                    "A_OA_O.Female", "A_OA_O.Male")
  
  # OO x OO
  
  a_genotype_6 <- c("A_OA_O.Female", "A_OA_O.Male")
  
  # XY
  
  
  # II x IY_I
  
  xy_genotype_1 <- c("X_IX_I.Female", "X_IY_I.Male")
  
  # II x IY_O
  
  xy_genotype_2 <- c("X_IX_I.Female", "X_IY_O.Male")
  
  # II x OY_I
  
  xy_genotype_3 <- c("X_IX_O.Female", "X_IY_I.Male")
  
  # II x OY_O
  
  xy_genotype_4 <- c("X_IX_O.Female", "X_IY_O.Male")
  
  # IO x IY_I
  
  xy_genotype_5 <- c("X_IX_I.Female", "X_IY_I.Male",
                     "X_IX_O.Female", "X_OY_I.Male")
  
  # IO x IY_O
  
  xy_genotype_6 <- c("X_IX_I.Female", "X_IY_O.Male", 
                     "X_IX_O.Female", "X_OY_O.Male")
  
  # IO x OY_I
  
  xy_genotype_7 <- c("X_IX_O.Female", "X_IY_I.Male",
                     "X_OX_O.Female", "X_OY_I.Male")
  
  # IO x OY_O
  
  xy_genotype_8 <- c("X_IX_O.Female", "X_IY_O.Male",
                     "X_OX_O.Female", "X_OY_O.Male")
  
  # OO x IY_I
  
  xy_genotype_9 <- c("X_IX_O.Female", "X_OY_I.Male")
  
  # OO x IY_O
  
  xy_genotype_10 <- c("X_IX_O.Female", "X_OY_O.Male")
  
  # OO x OY_I
  
  xy_genotype_11 <- c("X_OX_O.Female", "X_OY_I.Male")
  
  # OO x OY_O
  
  xy_genotype_12 <- c("X_OX_O.Female", "X_OY_O.Male")
  
  # ZW
  
  # IW_I x II
  
  zw_genotype_1 <- c("Z_IZ_I.Male", "Z_IW_I.Female")
  
  # IW_I x IO
  
  zw_genotype_2 <- c("Z_IZ_I.Male", "Z_IZ_O.Male", 
                     "Z_IW_I.Female", "Z_OW_I.Male")
  
  # IW_I x OO
  
  zw_genotype_3 <- c("Z_IZ_O.Male", "Z_OW_I.Female")
  
  # IW_O x II
  
  zw_genotype_4 <- c("Z_IZ_I.Male", "Z_IW_O.Female")
  
  # IW_O x IO
  
  zw_genotype_5 <- c("Z_IZ_I.Male", "Z_IZ_O.Male",
                     "Z_IW_O.Female", "Z_OW_O.Female")
  
  # IW_O x OO
  
  zw_genotype_6 <- c("Z_IZ_O.Male", "Z_OW_O.Female")
  
  # OW_I X II
  
  zw_genotype_7 <- c("Z_IZ_O.Male", "Z_IW_I.Female")
  
  # OW_I x IO
  
  zw_genotype_8 <- c("Z_IZ_O.Male", "Z_OZ_O.Male",
                     "Z_IW_I.Female", "Z_OW_I.Female")
  
  # OW_I x OO
  
  zw_genotype_9 <- c("Z_OZ_O.Male", "Z_OW_I.Female")
  
  # OW_O X II
  
  zw_genotype_10 <- c("Z_IZ_O.Male", "Z_IW_O.Female")
  
  # OW_O x IO
  
  zw_genotype_11 <- c("Z_IZ_O.Male", "Z_OZ_O.Male",
                      "Z_IW_O.Female", "Z_OW_O.Female")
  
  # OW_O x OO
  
  zw_genotype_12 <- c("Z_OZ_O.Male", "Z_OW_O.Female")
  
  # cytoplasmic
  
  # I x I
  # I x O
  
  c_genotype_1 <- c("C_I.Female", "C_I.Male")
  
  # O x O
  # O x I
  
  c_genotype_2 <- c("C_O.Female", "C_O.Male")
  
  
  
  # Now calculate the zygote frequencies for each cross
  
  # autosomal
  
  # even frequency of two offspring genotypes
  
  freq_2 <- rep(0.5, 2)
  
  # even frequency between four offspring types
  
  freq_4 <- rep(0.25, 4)
  
  # when there are 6 offspring genotypes
  
  freq_6 <- c(0.125, 0.125,
              0.25, 0.25,
              0.125, 0.125)
  
  bind_rows(
    list(
      make_offspring("A_IA_I", "A_IA_I", a_genotype_1, freq_2, "autosomal"),
      make_offspring("A_IA_I", "A_IA_O", a_genotype_2, freq_4, "autosomal"),
      make_offspring("A_IA_I", "A_OA_O", a_genotype_3, freq_2, "autosomal"),
      make_offspring("A_IA_O", "A_IA_I", a_genotype_2, freq_4, "autosomal"),
      make_offspring("A_IA_O", "A_IA_O", a_genotype_4, freq_6, "autosomal"),
      make_offspring("A_IA_O", "A_OA_O", a_genotype_5, freq_4, "autosomal"),
      make_offspring("A_OA_O", "A_IA_I", a_genotype_3, freq_2, "autosomal"),
      make_offspring("A_OA_O", "A_IA_O", a_genotype_5, freq_4, "autosomal"),
      make_offspring("A_OA_O", "A_OA_O", a_genotype_6, freq_2, "autosomal"),
      
      make_offspring("X_IX_I", "X_IY_I", xy_genotype_1, freq_2, "XY"),
      make_offspring("X_IX_I", "X_IY_O", xy_genotype_2, freq_2, "XY"),
      make_offspring("X_IX_I", "X_OY_I", xy_genotype_3, freq_2, "XY"),
      make_offspring("X_IX_I", "X_OY_O", xy_genotype_4, freq_2, "XY"),
      make_offspring("X_IX_O", "X_IY_I", xy_genotype_5, freq_4, "XY"),
      make_offspring("X_IX_O", "X_IY_O", xy_genotype_6, freq_4, "XY"),
      make_offspring("X_IX_O", "X_OY_I", xy_genotype_7, freq_4, "XY"),
      make_offspring("X_IX_O", "X_OY_O", xy_genotype_8, freq_4, "XY"),
      make_offspring("X_OX_O", "X_IY_I", xy_genotype_9, freq_2, "XY"),
      make_offspring("X_OX_O", "X_IY_O", xy_genotype_10, freq_2, "XY"),
      make_offspring("X_OX_O", "X_OY_I", xy_genotype_11, freq_2, "XY"),
      make_offspring("X_OX_O", "X_OY_O", xy_genotype_12, freq_2, "XY"),
      
      make_offspring("Z_IW_I", "Z_IZ_I", zw_genotype_1, freq_2, "ZW"),
      make_offspring("Z_IW_I", "Z_IZ_O", zw_genotype_2, freq_4, "ZW"),
      make_offspring("Z_IW_I", "Z_OZ_O", zw_genotype_3, freq_2, "ZW"),
      make_offspring("Z_IW_O", "Z_IZ_I", zw_genotype_4, freq_2, "ZW"),
      make_offspring("Z_IW_O", "Z_IZ_O", zw_genotype_5, freq_4, "ZW"),
      make_offspring("Z_IW_O", "Z_OZ_O", zw_genotype_6, freq_2, "ZW"),
      make_offspring("Z_OW_I", "Z_IZ_I", zw_genotype_7, freq_2, "ZW"),
      make_offspring("Z_OW_I", "Z_IZ_O", zw_genotype_8, freq_4, "ZW"),
      make_offspring("Z_OW_I", "Z_OZ_O", zw_genotype_9, freq_2, "ZW"),
      make_offspring("Z_OW_O", "Z_IZ_I", zw_genotype_10, freq_2, "ZW"),
      make_offspring("Z_OW_O", "Z_IZ_O", zw_genotype_11, freq_4, "ZW"),
      make_offspring("Z_OW_O", "Z_OZ_O", zw_genotype_12, freq_2, "ZW"),
      
      make_offspring("C_I", "C_I", c_genotype_1, freq_2, "cytoplasmic"),
      make_offspring("C_I", "C_O", c_genotype_1, freq_2, "cytoplasmic"),
      make_offspring("C_O", "C_I", c_genotype_2, freq_2, "cytoplasmic"),
      make_offspring("C_O", "C_O", c_genotype_2, freq_2, "cytoplasmic")
    )) %>% 
    filter(locus_type == gene_location)
}

# load the autosomal inheritance table for faster computation

offspring_genotypes_autosome <- 
  make_mating_table(gene_location = "autosomal") %>% 
  select(1:4) %>% 
  rename(zygote_type = type) %>% 
  separate_wider_delim(zygote_type, names = c("Genotype", "Sex"), delim = ".") %>% 
  mutate(Sex = if_else(Sex == "Female", 1, 0)) %>% 
  as.data.table()

# main simulation function

continuous_time_simulation <- function(row,
                                       parameters,
                                       inheritance_scheme){
  library(data.table) # note that we must load these packages within the function for HPC to work
  library(tidyverse)
  print(paste("Doing row", row)) # this tells you which row in the parameter space is being modelled
  
  # prop_i_table <- data.table(time = numeric(), 
  #                           proportion_I = numeric(),
  #                          population_size = numeric()) # filled in as sim progresses
  
  keep_going <- TRUE # if the inbreeding allele fixes or goes extinct, this will change to false and the while loop will quit early
  
  Starting_pop_size <- parameters$Starting_pop_size[row] 
  N <- parameters$N[row] # constant
  number_mutants <- parameters$number_mutants[row] # constant at 1 
  baseline_mean_lifespan <- parameters$baseline_mean_lifespan[row] # constant at 2
  time_end <- parameters$time_end[row] # a cut-off point for each run 
  sex_expressed <- parameters$sex_expressed[row]
  chromosome <- parameters$chromosome[row]
  heterozygous_genotype <- parameters$heterozygous_genotype[row]
  homozygous_genotype <- parameters$homozygous_genotype[row]
  hemizygous_genotype <- parameters$hemizygous_genotype[row]
  #C <- parameters$C[row]
  v <- parameters$v[row]
  refractory_period <- parameters$refractory_period[row]
  D <- parameters$D[row]
  dominance <- parameters$dominance[row]
  
  # define the starting genotypes for each sex so the population table can be built
  
  Female_starting_genotype <- inheritance_scheme[.N]$Female_genotype
  
  Male_starting_genotype <- inheritance_scheme[.N]$Male_genotype
  
  # Set the number of breeding sites
  
  breeding_sites <- 0.5*Starting_pop_size
  
  # Initialize the Individual_ID  and Family_ID counters
  
  Individual_ID_counter <- Starting_pop_size
  
  Family_ID_counter <- Starting_pop_size/N # family size equals the no. offspring produced by a single female
  
  # the simulation tracks the population via a data.table
  
  # create the starting population - note that females are sex = 1 and males are sex = 0
  
  population <-
    data.table(mortality_rate = 1/baseline_mean_lifespan,
               Sex = rbinom(n = Starting_pop_size, 1, prob = 0.5),
               birth_time = 0,
               matings = 0,
               reproduced = 0,
               mated_with = "NA",
               inbred_mating = 0,
               refractory_period_end = 0
    )[, `:=` (Genotype = ifelse(Sex > 0, Female_starting_genotype, Male_starting_genotype),
              breeding = Sex,
              Individual_ID = .I,
              Family_ID = rep(1:(.N/N), each = N, length.out = .N))]
  
  # seed population with the inbreeding allele
  
  population[sample(which(str_detect(Genotype, pattern = chromosome)), size = number_mutants, replace = F),
             Genotype := str_replace(Genotype, pattern = paste0(chromosome, "_O"), replacement = paste0(chromosome, "_I"))]
  
  # Determining the next event
  
  # check when the next death occurs
  
  next_death <- rexp(n = 1, rate = sum(population[, mortality_rate]))
  
  who_died <- population[sample(.N, 1, prob = mortality_rate)]
  
  # check when the next female-male encounter occurs
  
  # create the possible encounter list
  
  number_females <- sum(population$Sex > 0)
  
  encounter_possibilities <- 
    CJ(Female_ID = population[Sex > 0]$Individual_ID,
       Male_ID = population[Sex < 1]$Individual_ID)[, encounter_rate := v/number_females]
  
  # check the time
  
  next_mating <- rexp(n = 1, rate = sum(encounter_possibilities[, encounter_rate]))
  
  which_encounter <- encounter_possibilities[sample(.N, 1, prob = encounter_rate)]
  
  
  # Initialize the timer t to the first encounter
  
  t <- pmin(next_death, next_mating)
  
  # With the initial population ready to go and the first event found, start the timer and let the simulation run. In short, time progresses as events occur. Events can trigger state changes for the individuals in the population, leading to death, mating and offspring production.
  
  while (t <= time_end & keep_going) {
    
    # what type of encounter happens at time t
    
    if(next_death < next_mating){
      
      # Male mortality events
      
      if(who_died[,Sex] < 1){
        
        # remove male from population table
        
        population <- population[!Individual_ID %chin% who_died$Individual_ID]
        
      }
      
      # Female mortality events
      
      if(who_died[,Sex] > 0){
        
        # remove female from population table and candidate list
        
        population <- population[!Individual_ID %chin% who_died$Individual_ID]
        
        # check if female mortality event frees up breeding site
        
        current_breeders <- sum(population$breeding > 0)
        
        # If there is an available breeding site, and at least one female to fill it, recruit a new breeder
        
        if(current_breeders < breeding_sites && sum(population$Sex > 0 & population$breeding < 1) > 0){
          
          # assign the new breeders
          
          population <- population[sample(which(breeding < 1 & Sex > 0),
                                          size = 1), # note that all living females have equal prob of becoming a breeder
                                   breeding := 1]
        }
      }
    } else{
      
      # Opposite sex encounters
      
      female <- population[Individual_ID %chin% which_encounter$Female_ID]
      male <- population[Individual_ID %chin% which_encounter$Male_ID]
      
      # if the encounter occurs involving an individual expressing the A_I allele, provide an opportunity for inbreeding. The hurdle requirement for an inbreeding opportunity is that an individual must stay alive long enough to meet any member of the opposite sex. We then code the simulation such that this opposite sex individual is swapped out for a full-sibling that can be mated with. This simulates a common situation in nature, where due to population viscosity, relatives live in close geographic proximity and are thus more likely to be encountered, or to be encountered early in life. 
      
      # First determine if a homogametic individual carrying one copy of the I allele will inbreed on this occasion.  
      
      heterozygote_inbreeds <- rbinom(1, 1, prob = dominance)
      
      # find males that will inbreed when they are the sex that expresses inbreeding tolerance
      
      # male heterozygotes (FALSE if there isn't a heterozygote male genotype)
      
      if(sex_expressed < 1 & 
         str_detect(male$Genotype, heterozygous_genotype) & 
         heterozygote_inbreeds > 0 &
         male$refractory_period_end < t){
        
        # find a new mate (sister) for the males
        
        sister_mating <-
          population[Family_ID %chin% male$Family_ID
          ][Sex > 0, .SD[sample(1, 1)] # assign a sibling, then check if they're receptive
          ][matings < 1][, mated_with := male[,Genotype]] 
        
        # if mating occurred, update the population
        
        if(nrow(sister_mating)>0){
          
          mates <- rbindlist(list(male, sister_mating))
          
          population[mates,
                     `:=`(matings = matings + 1,
                          inbred_mating = 1 * (Sex > 0),
                          mated_with = ifelse(Sex > 0, i.mated_with, NA),
                          refractory_period_end = (t + refractory_period * baseline_mean_lifespan) * (Sex < 1)),
                     on = .(Individual_ID)]
        }
      }
      
      # male homozygotes (also includes cases where males only ever have one copy of the chromosome e.g. X, Y, C)
      
      if(sex_expressed < 1 & 
         str_detect(male$Genotype, homozygous_genotype) &
         male$refractory_period_end < t){
        
        # find a new mate (sister) for the males
        
        sister_mating <-
          population[Family_ID %chin% male$Family_ID
          ][Sex > 0, .SD[sample(1, 1)] # assign a sibling, then check if they're receptive
          ][matings < 1][, mated_with := male[,Genotype]] 
        
        # if mating occurred, update the population
        
        if(nrow(sister_mating)>0){
          
          mates <- rbindlist(list(male, sister_mating))
          
          population[mates,
                     `:=`(matings = matings + 1,
                          inbred_mating = 1 * (Sex > 0),
                          mated_with = ifelse(Sex > 0, i.mated_with, NA),
                          refractory_period_end = (t + refractory_period * baseline_mean_lifespan) * (Sex < 1)),
                     on = .(Individual_ID)]
        }
      }
      
      # find females that will inbreed when they are the sex that expresses inbreeding tolerance
      
      # female heterozygotes (FALSE if there isn't a heterozygote female genotype)
      
      if(sex_expressed > 0 & 
         str_detect(female$Genotype, heterozygous_genotype) & 
         heterozygote_inbreeds > 0 &
         female$matings < 1){
        
        # find a new mate (brother) for the female
        
        brother_mating <-
          population[Family_ID %chin% female$Family_ID
          ][Sex < 1, .SD[sample(1, 1)] # assign a sibling, then check if they're receptive
          ][t > refractory_period_end] 
        
        # if mating occurred, update the population
        
        if(nrow(brother_mating)>0){
          
          mates <- rbindlist(list(female[, mated_with := brother_mating[,Genotype]], brother_mating))
          
          population[mates,
                     `:=`(matings = matings + 1,
                          inbred_mating = 1 * (Sex > 0),
                          mated_with = ifelse(Sex > 0, i.mated_with, NA),
                          refractory_period_end = (t + refractory_period * baseline_mean_lifespan) * (Sex < 1)),
                     on = .(Individual_ID)]
        }
      }
      
      # female homozygotes (also includes cases where females only ever have one copy of the chromosome e.g. Z, W, C)
      
      if(sex_expressed > 0 & 
         str_detect(female$Genotype, homozygous_genotype) &
         female$matings < 1){
        
        # find a new mate (brother) for the females
        
        brother_mating <-
          population[Family_ID %chin% female$Family_ID
          ][Sex < 1, .SD[sample(1, 1)] # assign a sibling, then check if they're receptive
          ][t > refractory_period_end] 
        
        # if mating occurred, update the population
        
        if(nrow(brother_mating)>0){
          
          mates <- rbindlist(list(female[, mated_with := brother_mating[,Genotype]], brother_mating))
          
          population[mates,
                     `:=`(matings = matings + 1,
                          inbred_mating = 1 * (Sex > 0),
                          mated_with = ifelse(Sex > 0, i.mated_with, NA),
                          refractory_period_end = (t + refractory_period * baseline_mean_lifespan) * (Sex < 1)),
                     on = .(Individual_ID)]
        }
      }
      
      # standard outbred mating
      
      if((female$matings < 1 & t > male$refractory_period_end) |
         (female$matings < 1 & t > male$refractory_period_end &
          sex_expressed < 1 & str_detect(male$Genotype, heterozygous_genotype) & heterozygote_inbreeds < 1) |
         (female$matings < 1 & t > male$refractory_period_end &
          sex_expressed > 0 & str_detect(female$Genotype, heterozygous_genotype) & heterozygote_inbreeds < 1)){
        
        # update the population
        
        mates <- rbindlist(list(female[, mated_with := male[,Genotype]], male))
        
        population[mates,
                   `:=`(matings = matings + 1,
                        inbred_mating = 0,
                        mated_with = ifelse(Sex > 0, i.mated_with, NA),
                        refractory_period_end = (t + refractory_period * baseline_mean_lifespan) * (Sex < 1)),
                   on = .(Individual_ID)] 
      }
      
    }
    
    # reproduction
    
    # check if a female can now produce offspring, either because they're previously mated and have secured a breeding site or because they already held a breeding site and have now mated
    
    new_mated_breeder <- population[Sex > 0 & matings > 0 & breeding > 0 & reproduced < 1, 
                                    .(Individual_ID,
                                      inbred_mating,
                                      Female_genotype = Genotype,
                                      Male_genotype = mated_with)]
    
    if(nrow(new_mated_breeder) > 0) {
      # add offspring to the population. Each mated female that holds a breeding site produces N offspring
      offspring <- 
        new_mated_breeder[inheritance_scheme, 
                          on = .(Female_genotype = Female_genotype,
                                 Male_genotype = Male_genotype), 
                          nomatch = NULL, allow.cartesian  = TRUE
        ][, .SD[sample(.N, # I think this value is very occasionally zero, causing the sim to break
                       size = N, 
                       prob = zygote_freq, 
                       replace = T)]
        ][, Family_ID := .GRP + Family_ID_counter # assign these offspring to a new family 
        ][, .(Genotype, 
              Sex, 
              inbred_mating,
              Family_ID)
        ][, `:=`(mortality_rate = ifelse(inbred_mating > 0, 
                                         1/(baseline_mean_lifespan + D), # the cost of inbreeding: D <= 0
                                         1/baseline_mean_lifespan), # outbred offspring mortality risk
                 birth_time = t,
                 breeding = 0,
                 matings = 0,
                 reproduced = 0,
                 mated_with = "NA",
                 refractory_period_end = t,
                 Individual_ID = .I + Individual_ID_counter,
                 inbred_mating = 0)]
      
      # bind the offspring table to the existing population table and update which females have reproduced 
      
      population <- rbindlist(list(population, offspring), use.names = TRUE
      )[new_mated_breeder, reproduced := 1, on = .(Individual_ID)]
      
      # update the Individual_ID counter
      Individual_ID_counter <- max(population$Individual_ID)
      Family_ID_counter <- max(population$Family_ID) 
    }
    
    
    # Calculate the frequency of the I allele, quit early if I fixes or goes extinct
    
    # calc allele freq if autosomal locus   
    if(chromosome == "A"){
      prop_i <-
        (length(population$Genotype[str_detect(population$Genotype, heterozygous_genotype)]) + 
           2*length(population$Genotype[str_detect(population$Genotype, homozygous_genotype)]))/ (nrow(population)*2)
    }
    
    # calc allele freq if hemizygous locus: W, Y or cytoplasmic
    if(chromosome == "W" | chromosome == "Y" | chromosome == "C"){
      prop_i <-
        (length(population$Genotype[str_detect(population$Genotype, hemizygous_genotype)])/ 
           length(population$Genotype[str_detect(population$Genotype, chromosome)]))
    }
    
    # calc allele freq if diploid in one sex and haploid in the other: X and Z
    if(chromosome == "X" | chromosome == "Z"){
      prop_i <-
        if(hemizygous_genotype == "X_IY_O"){
          (length(population$Genotype[str_detect(population$Genotype, heterozygous_genotype)]) + 
             2*length(population$Genotype[str_detect(population$Genotype, homozygous_genotype)]) +
             length(population$Genotype[str_detect(population$Genotype, hemizygous_genotype)]))/ 
            (nrow(population[Sex > 0])*2 + nrow(population[Sex < 1]))}
      else{
        (length(population$Genotype[str_detect(population$Genotype, heterozygous_genotype)]) + 
           2*length(population$Genotype[str_detect(population$Genotype, homozygous_genotype)]) +
           length(population$Genotype[str_detect(population$Genotype, hemizygous_genotype)]))/ 
          (nrow(population[Sex < 1])*2 + nrow(population[Sex > 0]))}
    }
    
    # this is a diagnostic to make sure the model is running well - it can be commented out when running the big simulation
    #prop_i_table <- rbindlist(list(prop_i_table, list(t, prop_i, nrow(population))))
    
    #print(paste0("Population size = ", nrow(population),
    #            ", time = ", round(t, 3)))
    
    if(prop_i > 0.95 | prop_i < 0.0001 | nrow(population) < 2) keep_going <- FALSE
    
    # Move t to next encounter
    
    # determining the next event
    
    # check when the next death occurs
    
    next_death <- rexp(n = 1, rate = sum(population[, mortality_rate]))
    
    who_died <- population[sample(.N, 1, prob = mortality_rate)]
    
    # check when the next female-male encounter occurs
    
    # create the possible encounter list
    
    number_females <- sum(population$Sex > 0)
    
    encounter_possibilities <- 
      CJ(Female_ID = population[Sex > 0]$Individual_ID,
         Male_ID = population[Sex < 1]$Individual_ID)[, encounter_rate := v/number_females]
    
    # check the time
    
    next_mating <- rexp(n = 1, rate = sum(encounter_possibilities[, encounter_rate]))
    
    which_encounter <- encounter_possibilities[sample(.N, 1, prob = encounter_rate)]
    
    # Initialize the timer t to the next encounter
    
    t <- t + pmin(next_death, next_mating)  
    
  }
  finish_time <- t
  final_pop_size = nrow(population)
  
  #prop_i_table
  
  # Print the simulation results
  #list(
  results <- 
    parameters[row, ] %>% mutate(I_frequency = prop_i,
                                 finish_time = finish_time,
                                 final_pop_size = final_pop_size) %>% 
    select(-c(heterozygous_genotype, homozygous_genotype, hemizygous_genotype,
              number_mutants, baseline_mean_lifespan, N, time_end)) %>% 
    as.data.frame()
  
  write_csv(results, paste("/lustre/miifs01/project/m2_jgu-tee/tk_output/sim_results_", row, ".csv", sep = ""))
  # population,
  #prop_i_table)
  
}

# find parameter space

resolution <- 30
starting_pop_size_autosomes <- 200 # both sexes harbour two copies of each autosomal chromosome = 400 autosomal haplotypes

parameters <-
  expand_grid(
    #C = c(2, 10),
    chromosome = c("A", "X", "Y", "Z", "W", "C"),
    v = c(1, 5, 50),
    D = seq(0, -0.99, length = resolution), # inbreeding depression
    refractory_period = seq(0, 1, length = resolution)
  ) %>% 
  full_join(tibble(chromosome = c("A", "A", "A", "A", "A", "A",
                                  "X", "X", "X", "X",
                                  "Y",
                                  "Z", "Z", "Z", "Z",
                                  "W",
                                  "C", "C"),
                   sex_expressed = c(0, 0, 0, 1, 1, 1,
                                     0, 1, 1, 1,
                                     0,
                                     0, 0, 0, 1,
                                     1,
                                     0, 1),
                   dominance = c(0, 0.5, 1, 0, 0.5, 1,
                                 1, 0, 0.5, 1, 
                                 1, 
                                 0, 0.5, 1, 1,
                                 1,
                                 1, 1)) %>% 
              mutate(heterozygous_genotype = case_when(chromosome == "A" ~ "A_IA_O",
                                                       chromosome == "X" ~ "X_IX_O",
                                                       chromosome == "Y" ~ "NA",
                                                       chromosome == "Z" ~ "Z_IZ_O",
                                                       chromosome == "W" ~ "NA",
                                                       chromosome == "C" ~ "NA"),
                     homozygous_genotype = case_when(chromosome == "A" ~ "A_IA_I",
                                                     chromosome == "X" ~ "X_IX_I",
                                                     chromosome == "Y" ~ "NA",
                                                     chromosome == "Z" ~ "Z_IZ_I",
                                                     chromosome == "W" ~ "NA",
                                                     chromosome == "C" ~ "NA"),
                     hemizygous_genotype = case_when(chromosome == "A" ~ "NA",
                                                     chromosome == "X" ~ "X_IY_O",
                                                     chromosome == "Y" ~ "X_OY_I",
                                                     chromosome == "Z" ~ "Z_IW_O",
                                                     chromosome == "W" ~ "Z_OW_I",
                                                     chromosome == "C" ~ "C_I"),
                     Starting_pop_size = case_when(chromosome == "A" ~ starting_pop_size_autosomes,
                                                   chromosome == "X" | chromosome == "Z" ~ 
                                                     starting_pop_size_autosomes + 0.33*starting_pop_size_autosomes,
                                                   chromosome == "Y" | chromosome == "W" ~ starting_pop_size_autosomes*4,
                                                   chromosome == "C" ~ starting_pop_size_autosomes*2)),
            relationship = "many-to-many", by = "chromosome") %>% 
  mutate(baseline_mean_lifespan = 1,
         N = 5, # subject to change
         number_mutants = 4, # this equates to a starting frequency of 0.01
         time_end = 1000, # with avg lifespan = 1, this is ~ roughly 1000 gens
         parameter_space_ID = row_number())

parameters_autosome <- parameters %>% filter(chromosome == "A")

# this is the part we need to run in parallel

#results <- lapply(1:2, continuous_time_simulation, 
#           parameters = parameters_autosome, 
#          offspring_genotypes_autosome)

# here's the snow version of the function, note that the cl argument is a cluster object that I haven't specified yet

results <- clusterApply(cl, 1:nrow(parameters_autosome), continuous_time_simulation,
                        parameters = parameters_autosome, 
                        offspring_genotypes_autosome)

all_results <- do.call(rbind, results)


write_csv(all_results, "/lustre/miifs01/project/m2_jgu-tee/tk_output/sim_results_complete.csv")



