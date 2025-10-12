#!/usr/bin/env Rscript
options(echo=TRUE) # if you want see commands in output file
# Expect command line args at the end.
args = commandArgs(trailingOnly = TRUE)
# Skip args[1] to prevent getting --args

params<-read.table(args[2])
row_number= as.numeric(args[3])
print(row_number)
print(params[row_number,])

library(dplyr)

sample_vec <- function(x, ...) x[sample(length(x), ...)] 

make_mating_table <- function(gene_location){
  
  make_offspring <- function(X, Y, offspring_genotype, zygote_freq, sex){
    data.frame(Female_genotype = X,
               Male_genotype = Y,
               offspring_genotype,
               zygote_freq,
               sex)
  }
  
  # Specify the possible offspring genotypes for all the potential crosses; we use these for the offspring_genotype argument in the make_offspring function
  
  # offspring genotypes
  
  offspring_genotypes_1 <- c(2,2)
  offspring_genotypes_2 <- c(1, 1, 2, 2)
  offspring_genotypes_3 <- c(1, 1)
  offspring_genotypes_4 <- c(0, 0, 1, 1, 2, 2)
  offspring_genotypes_5 <- c(0, 0, 1, 1)
  offspring_genotypes_6 <- c(0, 0)
  
  offspring_genotypes_7 <- c(1, 2)
  offspring_genotypes_8 <- c(0, 1, 1, 2)
  offspring_genotypes_9 <- c(0, 1)
  offspring_genotypes_10 <- c(1, 0) # this is diff from above bc of the order with the sexes
  offspring_genotypes_11 <- c(2, 1)
  offspring_genotypes_12 <- c(2,1,1,0)
  
  # offspring sex
  
  offspring_sex_2 <- c(0, 1)
  offspring_sex_4 <- c(0, 1, 0, 1)
  offspring_sex_6 <- c(0, 1, 0, 1, 0, 1)
  
  # even frequency of two offspring genotypes
  
  freq_2 <- rep(0.5, 2)
  
  # even frequency between four offspring types
  
  freq_4 <- rep(0.25, 4)
  
  # when there are 6 offspring genotypes
  
  freq_6 <- c(0.125, 0.125,
              0.25, 0.25,
              0.125, 0.125)
  
  if(gene_location == "A"){
    books <- rbind(
      make_offspring(2, 2, offspring_genotypes_1, freq_2, offspring_sex_2),
      make_offspring(2, 1, offspring_genotypes_2, freq_4, offspring_sex_4),
      make_offspring(2, 0, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(1, 2, offspring_genotypes_2, freq_4, offspring_sex_4),
      make_offspring(1, 1, offspring_genotypes_4, freq_6, offspring_sex_6),
      make_offspring(1, 0, offspring_genotypes_5, freq_4, offspring_sex_4),
      make_offspring(0, 2, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(0, 1, offspring_genotypes_5, freq_4, offspring_sex_4),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "X"){
    books <- rbind(
      make_offspring(2, 1, offspring_genotypes_7, freq_2, offspring_sex_2),
      make_offspring(2, 0, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(1, 1, offspring_genotypes_8, freq_4, offspring_sex_4),
      make_offspring(1, 0, offspring_genotypes_5, freq_4, offspring_sex_4),
      make_offspring(0, 1, offspring_genotypes_9, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "Y"){
    books <- rbind(
      make_offspring(0, 1, offspring_genotypes_10, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "Z"){
    books <- rbind(
      make_offspring(1, 2, offspring_genotypes_11, freq_2, offspring_sex_2),
      make_offspring(0, 2, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(1, 1, offspring_genotypes_12, freq_4, offspring_sex_4),
      make_offspring(0, 1, offspring_genotypes_5, freq_4, offspring_sex_4),
      make_offspring(1, 0, offspring_genotypes_10, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "W"){
    books <- rbind(
      make_offspring(1, 0, offspring_genotypes_9, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "C"){
    books <- rbind(
      make_offspring(1, 0, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(1, 1, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2),
      make_offspring(0, 1, offspring_genotypes_6, freq_2, offspring_sex_2)
    )
  }
  
  if(gene_location == "P"){
    books <- rbind(
      make_offspring(1, 0, offspring_genotypes_6, freq_2, offspring_sex_2),
      make_offspring(1, 1, offspring_genotypes_3, freq_2, offspring_sex_2),
      make_offspring(0, 0, offspring_genotypes_6, freq_2, offspring_sex_2),
      make_offspring(0, 1, offspring_genotypes_3, freq_2, offspring_sex_2)
    )
  }
  return(books)  
}

offspring_genotypes_autosome <- make_mating_table("A")
offspring_genotypes_X <- make_mating_table("X")
offspring_genotypes_Y <- make_mating_table("Y")
offspring_genotypes_Z <- make_mating_table("Z")
offspring_genotypes_W <- make_mating_table("W")
offspring_genotypes_C <- make_mating_table("C")
offspring_genotypes_P <- make_mating_table("P")

sample_mating_table <- function(inheritance_scheme, 
                                f,
                                mother){
  
  # cut to possible genotypes
  possibilities <- 
    inheritance_scheme[inheritance_scheme$Female_genotype == mother[4] &
                         inheritance_scheme$Male_genotype == mother[9], c(3,5)]
  # get prob of producing each genotype
  probs <- 
    inheritance_scheme[inheritance_scheme$Female_genotype == mother[4] &
                         inheritance_scheme$Male_genotype == mother[9], 4]
  # sample
  possibilities[sample(size = f,
                       x = nrow(possibilities), 
                       prob = probs,
                       replace = TRUE), ]
}


continuous_time_simulation <- function(row,
                                       parameters,
                                       inheritance_scheme){
  
  #print(paste("Doing row", row)) # this shows which row in the parameter space is being modelled
  
  Starting_pop_size <- round(parameters$Starting_pop_size[row], 0)
  f <- parameters$f[row] # fecundity constant
  mutation_time <- parameters$mutation_time[row] # introduce an I allele after family structure is established
  baseline_mean_lifespan <- parameters$baseline_mean_lifespan[row] # constant at 1
  time_end <- parameters$time_end[row] # a cut-off point for each run 
  sex_expressed <- parameters$sex_expressed[row]
  chromosome <- parameters$chromosome[row]
  v <- parameters$v[row]
  refractory_period <- parameters$refractory_period[row]
  D <- parameters$D[row]
  dominance <- parameters$dominance[row]
  parameter_space_ID <- parameters$parameter_space_ID[row]
  mutation_events <- parameters$mutation_events[row]
  
  # Set the number of breeding sites
  
  breeding_sites <- round(0.2*Starting_pop_size, 0)
  
  # what inheritance system does this run follow
  offspring_genotypes <- inheritance_scheme
  
  # Set the maximum number of I alleles that can be found in each sex
  if(chromosome == "A"){
    female_max_I <- 2
    male_max_I <- 2
  }
  
  if(chromosome == "X"){
    female_max_I <- 2
    male_max_I <- 1
  }
  
  if(chromosome == "Z"){
    female_max_I <- 1
    male_max_I <- 2
  }
  
  if(chromosome == "Y"){
    female_max_I <- 0
    male_max_I <- 1
  }
  
  if(chromosome == "W"){
    female_max_I <- 1
    male_max_I <- 0
  }
  
  if(chromosome == "C" | chromosome == "P"){
    female_max_I <- 1
    male_max_I <- 1
  }
  
  # make matrix to hold results; updated as sim progresses
  # col1 = time, col2 = prop I, col3 = pop size, col4 = prop virgin female deaths
  results_matrix <- matrix(nrow = time_end*4+2, ncol = 4) # record each time point
  
  # make matrix to hold population; updated as sim progresses
  
  # col1 = ID 
  # col2 = Family ID
  # col3 = Sex: females = 1 and males = 0
  # col4 = Genotype: 0, 1 and 2 = copies of inbreeding allele
  # col5 = mortality rate
  # col6 = encountered relative: NA = NO, 1 = YES
  # col7 = mating state: -Inf not in pop, NA = unmated, real = out, Inf = mated female
  # col8 = inbred mating: NA = NO, 1 = YES (only matters for females)
  # col9 = mated_genotype: NA = unmated, otherwise 0,1,2 (see mating table)
  # col10 = breeding site:  NA = NO, 1 = YES 
  # col11 = no. matings (only matters for males)
  # col12 = offspring produced: NA = NO, 1 = YES
  
  pop_matrix <- matrix(nrow = Starting_pop_size*2, # pop expands with initial repro pulse
                       ncol = 12)
  # ID & Family ID
  pop_matrix[1:Starting_pop_size, 1:2] <- 1:Starting_pop_size
  # assign sex
  pop_matrix[1:Starting_pop_size, 3] <- rbinom(n = Starting_pop_size, 1, prob = 0.5)
  # female_starting_genotype
  pop_matrix[pop_matrix[,3] < 1 & !is.na(pop_matrix[,3]), 4] <- 0
  # male_starting_genotype
  pop_matrix[pop_matrix[,3] > 0 & !is.na(pop_matrix[,3]), 4] <- 0
  # assign mortality rates
  pop_matrix[1:Starting_pop_size, 5] <- 1/baseline_mean_lifespan
  # set the unused rows to state -Inf 
  pop_matrix[(Starting_pop_size + 1):nrow(pop_matrix), 7] <- -Inf
  # mate count
  pop_matrix[1:Starting_pop_size, 11] <- 0
  # offspring production status
  pop_matrix[1:Starting_pop_size, 12] <- NA
  # populate breeding sites
  # the starting no. of females generally exceeds the number of breeding sites, which is starting_pop_size/f. The code below selects the initial breeding site holders
  
  if(nrow(pop_matrix[pop_matrix[,3] > 0 & !is.na(pop_matrix[,3]),]) > breeding_sites){
    
    initial_breeders <- head(pop_matrix[pop_matrix[,3] > 0 & !is.na(pop_matrix[,3]),1], breeding_sites)
    
    pop_matrix[initial_breeders,10] <- 1 # take advantage of ID = row number for initial pop
    
  } else{pop_matrix[pop_matrix[,3] > 0 & !is.na(pop_matrix[,3]), 10] <- 1}
  
  # Initialise counter for the results table
  
  next_update <- 0 # keep track of when to update the results
  next_row <- 0 # keep track of which row to update
  
  # Initialise the Individual_ID and Family_ID counters
  
  Individual_ID_counter <- Starting_pop_size
  
  Family_ID_counter <- Starting_pop_size # each individual descends from a distinct family at onset
  
  # Initialise the timer t
  
  t <- 0
  
  # Set initial pop size and freq of I allele for results table
  
  Prop_I <- 0 
  pop_size <- Starting_pop_size
  total_female_deaths <- 0
  mated_female_deaths <- 0
  
  # Start population without the I allele to generate family structure
  # Flips to 1 at mutant intro time point 
  
  mutant_introduced <- 0
  
  keep_going <- TRUE # if the inbreeding allele fixes or goes extinct, this will change to false and the while loop will quit early
  
  # With the initial population ready to go, start the timer and let the simulation run.
  
  while(t <= time_end & keep_going){
    
    print(paste0("Population size = ", pop_size, 
                 ", breeders = ", sum(pop_matrix[,10], na.rm = T), 
                 ", time = ", round(t, 3), ", Prop I =", Prop_I, ", mutation events =", mutant_introduced))
    
    # find next event 
    
    # next death: this is the sum of the mortality rates for all individuals in the population
    
    next_death <- t + rexp(n = 1, rate = sum(pop_matrix[, 5], na.rm = T))
    
    # next receptive mating encounter
    
    # find no. of females in mating pool & separate by encounter experience
    
    receptive_females_first_encounter <- 
      pop_matrix[pop_matrix[,3] > 0 &
                   is.na(pop_matrix[,6]) &
                   is.na(pop_matrix[,7]),, drop = FALSE]
    
    receptive_females_second_encounter <- 
      pop_matrix[pop_matrix[,3] > 0 &
                   !is.na(pop_matrix[,6]) &
                   is.na(pop_matrix[,7]),, drop = FALSE]
    
    # find no. of males in mating pool
    receptive_males <- pop_matrix[pop_matrix[,3] < 1 & is.na(pop_matrix[,7]),, drop = FALSE]
    
    # Find the time the next encounter occurs: plug the sum of the rates into the exponential function. 
    # The population level encounter rate is the product of the rate at which a single male finds a single female, the number of receptive females in the population, and the number of receptive males in the population
    
    next_first_encounter <- t + 
      rexp(n = 1, rate = v*nrow(receptive_females_first_encounter)*nrow(receptive_males))
    
    next_secondary_encounter <- t + 
      rexp(n = 1, rate = v*nrow(receptive_females_second_encounter)*nrow(receptive_males))
    
    # time in - Inf, Inf and Na are possible options that the code can handle 
    next_time_in <- min(pop_matrix[is.finite(pop_matrix[,7]),7])
    
    # find which event happens next and update t
    t <- pmin(next_death,
              next_time_in, 
              next_first_encounter,
              next_secondary_encounter,
              next_update, # update the population
              na.rm = TRUE) # ... if a rate is 0, NaN produced.
    
    
    if(t == next_update & !is.na(next_update)){# record time, I prop and pop size
      results_matrix[next_row+1,1] <- t
      results_matrix[next_row+1,2] <- round(Prop_I, 4)
      results_matrix[next_row+1,3] <- pop_size # popsize
      results_matrix[next_row+1,4] <- round(mated_female_deaths / total_female_deaths, 3)
      next_update <- next_update + 0.25
      next_row <- next_row + 1
      total_female_deaths <- 0 # reset the count
      mated_female_deaths <- 0 # reset the count
    }
    
    
    if(t == next_death){# remove an individual from the pop
      who_died <- 
        sample_vec(size = 1, # choose one
                   x = pop_matrix[!is.na(pop_matrix[,1]),1], # subset to current pop
                   prob = pop_matrix[!is.na(pop_matrix[,5]),5]) # weight by mortality rate
      # add a death if it was a female
      if(nrow(pop_matrix[pop_matrix[,1] == who_died &
                         !is.na(pop_matrix[,1]) &
                         pop_matrix[,3] > 0,, drop = FALSE]) > 0){total_female_deaths <- total_female_deaths + 1}
      
      # add virgin female deaths
      if(nrow(pop_matrix[pop_matrix[,1] == who_died &
                         !is.na(pop_matrix[,1]) &
                         pop_matrix[,3] > 0 &
                         is.infinite(pop_matrix[,7]),, drop = FALSE]) > 0){mated_female_deaths <- mated_female_deaths + 1}
      
      # remove individual from pop matrix
      pop_matrix[pop_matrix[,1] == who_died, 7] <- -Inf # NA means time-in here, so special edit required
      pop_matrix[pop_matrix[,1] == who_died, c(1:6, 8:12)] <- NA 
      
      # re-order to make steps like adding offspring easier later on
      pop_matrix <- pop_matrix[order(pop_matrix[,1]),]
      
    }
    
    # check if there are free breeding sites and whether females are available to fill them 
    
    current_breeders <- sum(pop_matrix[, 10], na.rm = T)
    
    # get list of IDs for floating females
    floating_females <- pop_matrix[!is.na(pop_matrix[,1]) & # alive
                                     pop_matrix[,3] > 0 & # female
                                     is.na(pop_matrix[,10]), # non-breeding
                                   1] # return the IDs only
    
    # If so, recruit a new breeder
    # All prospective females have equal probability
    
    if(current_breeders < breeding_sites & length(floating_females) > 0){
      
      # assign the new breeder
      
      new_breeder <- 
        sample_vec(size = 1, # choose one
                   x = floating_females) # subset to floaters
      
      pop_matrix[pop_matrix[,1] == new_breeder, 10] <- 1
    }
    
    if(t == next_time_in & !is.na(next_time_in)){ # a male re-enters the mating pool
      pop_matrix[pop_matrix[,7] == next_time_in, 7] <- NA # change to receptive
    }
    
    #### mating
    
    if(t == next_first_encounter &
       !is.na(next_first_encounter)){# does first encounter lead to (inbred) mating?
      
      # Determine whether a heterozygote inbreeds on this occasion. 
      # Depends on genotype if this matters
      heterozygote_inbreeds <- rbinom(1, 1, prob = dominance)
      
      # which female
      female_ID <- sample_vec(receptive_females_first_encounter[,1], 1)
      # get meta-data
      female <- subset(pop_matrix, pop_matrix[,1] == female_ID)
      # how many inbreeding alleles does she carry?
      alleles_female <- female[,4]
      
      mates <- NULL # reset this every time as a safeguard - MAYBE REMOVE?
      
      # find brothers that are in the mating pool
      brothers <-
        pop_matrix[pop_matrix[,2] == female[, 2] & # find family members
                     pop_matrix[,3] < 1 & # that are male
                     is.na(pop_matrix[,7]) & # and in the mating pool
                     !is.na(pop_matrix[,1]), # remove NAs
                   1] 
      # find the specific brother - if there aren't any, inbreeding does not happen
      if(length(brothers) > 0){# choose brother randomly
        chosen_brother <-
          subset(pop_matrix, 
                 pop_matrix[,1] == sample_vec(size = 1, x = brothers))
        # how many inbreeding alleles does he carry?
        alleles_brother <- chosen_brother[,4]
        brother_ID <- chosen_brother[,1]
      }else{alleles_brother <- 0} # we need this for the next if statement
      
      # now determine whether inbreeding occurs:
      # which individual expresses the allele
      # does that individual have the allele
      # is it expressed (depends on genomic region, no. copies and dominance)
      
      if(# female expression determines outcome
        # dominance doesn't matter
        length(brothers) > 0 & sex_expressed > 0 & female_max_I == alleles_female |
        # dominance matters
        length(brothers) > 0 & sex_expressed > 0 & 
        0 < alleles_female & alleles_female < female_max_I & heterozygote_inbreeds > 0 |
        # male expression determines outcome
        # dominance doesn't matter
        length(brothers) > 0 & sex_expressed < 1 & male_max_I == alleles_brother |
        # dominance matters
        length(brothers) > 0 & sex_expressed < 1 & 
        0 < alleles_brother & alleles_brother < male_max_I & heterozygote_inbreeds > 0){
        
        # do inbreeding
        # update the pop matrix
        # female
        pop_matrix[pop_matrix[,1] == female_ID, 6] <- 1 # relative has been encountered
        pop_matrix[pop_matrix[,1] == female_ID, 7] <- Inf # female leaves mating pool
        pop_matrix[pop_matrix[,1] == female_ID, 8] <- 1 # inbreeding occurs
        pop_matrix[pop_matrix[,1] == female_ID, 9] <- alleles_brother # mates genotype
        
        # male
        pop_matrix[pop_matrix[,1] == brother_ID, 7] <- t + refractory_period # male leaves mating pool
        pop_matrix[pop_matrix[,1] == brother_ID, 8] <- 1 # inbreeding occurs
        pop_matrix[pop_matrix[,1] == brother_ID & !is.na(pop_matrix[,1]), 11] <-
          pop_matrix[pop_matrix[,1] == brother_ID & !is.na(pop_matrix[,1]), 11] + 1
      } else{
        # inbreeding is avoided
        # females that had no receptive brother to encounter are recorded as having had their chance for inbreeding early in life. When the male refractory period != 0, this is possible but unlikely (because all siblings are produced at the same time). Most commonly, this will occur when a female produces an all-female brood (0.03125 probability when f=5)
        
        pop_matrix[pop_matrix[,1] == female_ID, 6] <- 1 # relative has been encountered
      }
    }
    
    if(t == next_secondary_encounter &
       !is.na(next_secondary_encounter)){ 
      # If the individual has already encountered a sibling, don't swap and let encounter proceed. 
      
      # which female
      female_ID <- sample_vec(receptive_females_second_encounter[,1], 1)
      # get meta-data
      female <- subset(pop_matrix, pop_matrix[,1] == female_ID)
      # how many inbreeding alleles does she carry?
      alleles_female <- female[,4]  
      
      # which male
      male_ID <- sample_vec(receptive_males[,1], 1)
      # get meta-data
      male <- subset(pop_matrix, pop_matrix[,1] == male_ID)
      # how many inbreeding alleles does he carry?
      alleles_male <- male[,4] 
      
      # If the pair happen to be siblings, check if they inbreed  
      
      # Determine whether a heterozygote inbreeds on this occasion. 
      # Depends on genotype if this matters
      heterozygote_inbreeds <- rbinom(1, 1, prob = dominance)
      
      if(
        # female expression determines outcome
        # dominance doesn't matter
        female[,2] == male[,2] & sex_expressed > 0 & female_max_I == alleles_female |
        # dominance matters
        female[,2] == male[,2] & sex_expressed > 0 & 
        0 < alleles_female & alleles_female < female_max_I & heterozygote_inbreeds > 0 |
        # male expression determines outcome
        # dominance doesn't matter
        female[,2] == male[,2] & sex_expressed < 1 & male_max_I == alleles_male |
        # dominance matters
        female[,2] == male[,2] & sex_expressed < 1 & 
        0 < alleles_male & alleles_male < male_max_I & heterozygote_inbreeds > 0){
        
        # do inbreeding
        # update the pop matrix
        # female
        pop_matrix[pop_matrix[,1] == female_ID, 7] <- Inf # female leaves mating pool
        pop_matrix[pop_matrix[,1] == female_ID, 8] <- 1 # inbreeding occurs
        pop_matrix[pop_matrix[,1] == female_ID, 9] <- alleles_male # mates genotype
        
        # male
        pop_matrix[pop_matrix[,1] == male_ID, 7] <- t + refractory_period # male leaves mating pool
        pop_matrix[pop_matrix[,1] == male_ID & !is.na(pop_matrix[,1]), 11] <-
          pop_matrix[pop_matrix[,1] == male_ID & !is.na(pop_matrix[,1]), 11] + 1
      } else{
        # do outbreeding
        # update the pop matrix
        # female
        pop_matrix[pop_matrix[,1] == female_ID, 7] <- Inf # female leaves mating pool
        pop_matrix[pop_matrix[,1] == female_ID, 9] <- alleles_male # mates genotype
        
        # male
        pop_matrix[pop_matrix[,1] == male_ID, 7] <- t + refractory_period # male leaves mating pool
        pop_matrix[pop_matrix[,1] == male_ID & !is.na(pop_matrix[,1]), 11] <-
          pop_matrix[pop_matrix[,1] == male_ID & !is.na(pop_matrix[,1]), 11] + 1
      }
    }
    
    # Consequences of death and mating: reproduction
    
    # check if a female can now produce offspring, either because they're previously mated and have secured a breeding site or because they already hold a breeding site and have now mated
    # make sure that previous breeders are excluded
    
    new_mated_breeder <- pop_matrix[is.infinite(pop_matrix[,7]) & # mated
                                      !is.na(pop_matrix[,10]) & # holds breeding site
                                      is.na(pop_matrix[,12]),, drop = FALSE] # hasn't reproduced
    
    if(nrow(new_mated_breeder) > 0){
      # add offspring to the population
      # each mated female that holds a breeding site produces f offspring
      
      # first check whether the mutant I allele should be added
      if(mutant_introduced < mutation_events & t > mutation_time){
        which_sex <- rbinom(1, 1, prob = 0.5)
        
        if(chromosome == "A" & which_sex == 1 |
           chromosome == "X" & which_sex == 1 |
           chromosome == "Z" & which_sex == 1){
          new_mated_breeder[4] <- 1
        }
        
        if(chromosome == "A" & which_sex == 0 |
           chromosome == "X" & which_sex == 0 |
           chromosome == "Z" & which_sex == 0){
          new_mated_breeder[9] <- 1
        }
        
        if(chromosome == "W"|
           chromosome == "C"){
          new_mated_breeder[4] <- 1
        }
        
        if(chromosome == "Y" |
           chromosome == "P"){
          new_mated_breeder[9] <- 1
        }
        
        mutant_introduced <- mutant_introduced + 1
      }
      
      next_row_to_fill <- length(pop_matrix[!is.na(pop_matrix[,1]),1]) + 1
      last_row_to_fill <- next_row_to_fill + f - 1
      next_ID <- Individual_ID_counter + 1
      last_ID <- Individual_ID_counter + f
      Family_ID_counter <- Family_ID_counter + 1
      
      # assign IDs
      pop_matrix[next_row_to_fill:last_row_to_fill, 1] <- next_ID:last_ID
      # assign all offspring to a single family
      pop_matrix[next_row_to_fill:last_row_to_fill, 2] <- Family_ID_counter
      # assign sex and genotype using our mating table sampling function
      offspring_genos <- 
        sample_mating_table(inheritance_scheme,
                            f, 
                            mother = new_mated_breeder)
      pop_matrix[next_row_to_fill:last_row_to_fill, 3] <- offspring_genos[,2]
      pop_matrix[next_row_to_fill:last_row_to_fill, 4] <- offspring_genos[,1]
      # assign mortality rates
      if(is.na(new_mated_breeder[8])){
        pop_matrix[next_row_to_fill:last_row_to_fill, 5] <- 1/baseline_mean_lifespan
      } else{ # apply effect of inbreeding depression
        pop_matrix[next_row_to_fill:last_row_to_fill, 5] <- 1/(baseline_mean_lifespan + D)
      }
      # fill in the mating and breeding site details - everyone starts as a floating virgin
      pop_matrix[next_row_to_fill:last_row_to_fill, 6:10] <- NA
      # mate count
      pop_matrix[next_row_to_fill:last_row_to_fill, 11] <- 0
      
      
      # update the mothers offspring production status
      
      pop_matrix[pop_matrix[,1] == new_mated_breeder[1], 12] <- 1
      
      # update the individual ID counter (redundant but more readable to do this here)
      Individual_ID_counter <- last_ID
      
    }      
    
    # Calculate the frequency of the I allele, quit early if I fixes or goes extinct
    
    pop_size <- nrow(pop_matrix[!is.na(pop_matrix[,1]),, drop = FALSE]) # use this to update the results
    n_females <- nrow(pop_matrix[!is.na(pop_matrix[,1]) &
                                   pop_matrix[,3] > 0,, drop = FALSE])
    n_males <- pop_size - n_females
    
    # calc allele freq if autosomal locus   
    if(chromosome == "A"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/(pop_size*2) # x2 because diploid
    }
    
    # calc allele freq if W locus   
    if(chromosome == "W"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/n_females
    }
    
    # calc allele freq if Y locus   
    if(chromosome == "Y"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/n_males 
    }
    
    # calc allele freq if X locus   
    if(chromosome == "X"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/(n_females*2 + n_males)
    }
    
    # calc allele freq if Z locus   
    if(chromosome == "Z"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/(n_females + n_males*2) 
    }
    
    # calc allele freq if C locus   
    if(chromosome == "C" |
       chromosome == "P"){
      Prop_I <-
        sum(pop_matrix[,4], na.rm = T)/pop_size 
    }
    
    # quit condition
    if(mutant_introduced > 0 & Prop_I > 0.9 |
       mutant_introduced > 0 & Prop_I == 0 | 
       pop_size < 2){keep_going <- FALSE}
    
  }
  
  results_matrix[next_row+1,1] <- t
  results_matrix[next_row+1,2] <- round(Prop_I, 4)
  results_matrix[next_row+1,3] <- pop_size
  results_matrix[next_row+1,4] <- round(mated_female_deaths / total_female_deaths, 3)
  results_matrix <- results_matrix[-(next_row+2:nrow(results_matrix)),]
  # save results as a csv.  
  
  results_matrix
  
  write.csv(results_matrix,
            paste("results/rowID_", 
                  parameter_space_ID, 
                  chromosome, ".csv", 
                  sep = ""))
  #write.csv(results_matrix,
  #        paste("sim_results/rowID_", 
  #            parameter_space_ID, 
  #          chromosome, ".csv", 
  #        sep = ""))
}

continuous_time_simulation(row = row_number, 
                           parameters = params, 
                           inheritance_scheme = offspring_genotypes_autosome)