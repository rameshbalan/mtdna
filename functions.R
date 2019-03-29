# this function takes the base population and returns a new population table
# that has the number of contributions each genotype makes to the gamete pool
# with this implementation we actually currently dont do anything with the
# selection fucntion.  However, I can imagine exploring how much of a fitness
# must be present to stop the spread of Medea so lets keep this framework handy
male_selection <- function(male_pop, Ne, s){
  # create an empty vector for absolute fitness
  male.abs.fit <- vector(mode = "numeric", length = 32)
  #If wanted to add selection against MEDEA
  # here we could decrement abs.fit by a value to represent the cost of
  # being medea positive.
  #males:
  
  # ("ABABM", "ABAbM", "ABaBM", "ABabM", "AbABM", "AbAbM", "AbaBM", "AbabM",
  #   "aBABM", "aBAbM", "aBaBM", "aBabM", "abABM", "abAbM", "abaBM", "ababM", 
  #   "ABABm", "ABAbm", "ABaBm", "ABabm", "AbABm", "AbAbm", "AbaBm", "Ababm",
  #   "aBABm", "aBAbm", "aBaBm", "aBabm", "abABm", "abAbm", "abaBm", "ababm")
  
  male.abs.fit[1:32] <- 1 
  
  
  # convert absolute fitness to relative or marginal fitness
  w_male <- male.abs.fit / max(male.abs.fit[male_pop != 0]) * (male_pop / (Ne/2))
  
  # use marginal fitness to draw from population
  as.vector(rmultinom(1:32, size = (Ne/2), prob = w_male))
  
}

female_selection <- function(female_pop, Ne, s){
  # create an empty vector for absolute fitness
  female.abs.fit <- vector(mode = "numeric", length = 32)
  #If wanted to add selection against MEDEA
  # here we could decrement abs.fit by a value to represent the cost of
  # being medea positive.
  #females:
  
  # ("ABABM", "ABAbM", "ABaBM", "ABabM", "AbABM", "AbAbM", "AbaBM", "AbabM",
  #   "aBABM", "aBAbM", "aBaBM", "aBabM", "abABM", "abAbM", "abaBM", "ababM", 
  #   "ABABm", "ABAbm", "ABaBm", "ABabm", "AbABm", "AbAbm", "AbaBm", "Ababm",
  #   "aBABm", "aBAbm", "aBaBm", "aBabm", "abABm", "abAbm", "abaBm", "ababm")
  
  female.abs.fit[1:32] <- 1  
  
  
  w_female <- female.abs.fit / max(female.abs.fit[female_pop != 0]) * (female_pop / (Ne/2))
  # use marginal fitness to draw from population
  as.vector(rmultinom(1:32, size = (Ne/2), prob = w_female))
}

makeGametes_male <- function(male_pop, r, Ne, U.a, U.b){
  # gamete pool: "ABM", "AbM", "aBM", "abM", "ABm", "Abm", "aBm", "abm")
  # M and m are variants of a mitochondrial gene from two different populations0
  #
  #create an empty vector for gametes
  x <- vector(length=8, mode="numeric")
  # ABM
  x[1] <- (sum(male_pop[1]*(1/4),male_pop[2]*(1/8),male_pop[3]*(1/8),male_pop[4]*((1/8)*(1-r)),male_pop[5]*(1/8),
               male_pop[7]*(1/8)*r,male_pop[9]*(1/8),male_pop[10]*(1/8)*r,male_pop[13]*((1/8)*(1-r))))/(Ne/2)
  # AbM
  x[2] <- (sum(male_pop[2]*(1/8),male_pop[4]*(1/8)*r,male_pop[5]*(1/8),male_pop[6]*(1/4),male_pop[7]*((1/8)*(1-r)),
               male_pop[8]*(1/8),male_pop[10]*((1/8)*(1-r)),male_pop[13]*(1/8)*r,male_pop[14]*(1/8)))/(Ne/2)
  # aBM
  x[3] <- (sum(male_pop[3]*(1/8),male_pop[4]*(1/8)*r,male_pop[7]*((1/8)*(1-r)),male_pop[9]*(1/8),male_pop[10]*((1/8)*(1-r)),
               male_pop[11]*(1/4),male_pop[12]*(1/8),male_pop[13]*(1/8)*r,male_pop[15]*(1/8)))/(Ne/2)
  # abM
  x[4] <- (sum(male_pop[4]*((1/8)*(1-r)),male_pop[7]*(1/8)*r,male_pop[8]*(1/8),male_pop[10]*(1/8)*r,male_pop[12]*(1/8),
               male_pop[13]*((1/8)*(1-r)),male_pop[14]*(1/8),male_pop[15]*(1/8),male_pop[16]*(1/4)))/(Ne/2)
  # ABm
  x[5] <- (sum(male_pop[17]*(1/4),male_pop[18]*(1/8),male_pop[19]*(1/8),male_pop[20]*((1/8)*(1-r)),male_pop[21]*(1/8),
               male_pop[23]*(1/8)*r,male_pop[25]*(1/8),male_pop[26]*(1/8)*r,male_pop[29]*((1/8)*(1-r))))/(Ne/2)
  # Abm
  x[6] <- (sum(male_pop[18]*(1/8),male_pop[20]*(1/8)*r,male_pop[21]*(1/8),male_pop[22]*(1/4),male_pop[23]*((1/8)*(1-r)),
               male_pop[24]*(1/8),male_pop[26]*((1/8)*(1-r)),male_pop[29]*(1/8)*r,male_pop[30]*(1/8)))/(Ne/2)
  # aBm
  x[7] <- (sum(male_pop[19]*(1/8),male_pop[20]*(1/8)*r,male_pop[23]*((1/8)*(1-r)),male_pop[25]*(1/8),male_pop[26]*((1/8)*(1-r)),
               male_pop[27]*(1/4),male_pop[28]*(1/8),male_pop[29]*(1/8)*r,male_pop[31]*(1/8)))/(Ne/2)
  # abm
  x[8] <- (sum(male_pop[20]*((1/8)*(1-r)),male_pop[23]*(1/8)*r,male_pop[24]*(1/8),male_pop[26]*(1/8)*r,male_pop[28]*(1/8),
               male_pop[29]*((1/8)*(1-r)),male_pop[30]*(1/8),male_pop[31]*(1/8),male_pop[32]*(1/4)))/(Ne/2)
  
  
  # Incorporation of mutation per gamete
  y <- vector(length=8, mode="numeric")
  # gamete pool: "ABM", "AbM", "aBM", "abM", "ABm", "Abm", "aBm", "abm")
  #ABM
  y[1] <- x[1] - U.a*x[1] - U.b*x[1] - U.a*U.b*x[1]
  #AbM
  y[2] <- x[2] - U.a*x[2] + U.b*x[1]
  #aBM
  y[3] <- x[3] - U.b*x[3] + U.a*x[1]
  #abM
  y[4] <- x[4] + U.a*x[2] + U.b*x[3] + U.a*U.b*x[1]
  #ABm
  y[5] <- x[5] - U.a*x[5] - U.b*x[5] - U.a*U.b*x[5]
  #Abm
  y[6] <- x[6] - U.a*x[6] + U.b*x[5]  
  #aBm
  y[7] <- x[7] - U.b*x[7] + U.a*x[5]
  #abm
  y[8] <- x[8]  + U.a*x[6] + U.b*x[7] + U.a*U.b*x[5]
  ## now we return this vector
  return(y)
  
}


makeGametes_female <- function(female_pop, r, Ne, U.a, U.b){
  # gamete pool: "ABM", "AbM", "aBM", "abM", "Ab'M", "a'BM", "ab'M", "a'bM", "a'b'M","ABm", "Abm", "aBm", "abm", "Ab'm", "a'Bm", "ab'm", "a'bm", "a'b'm"
  # in these cases a prime symbol indicates "poisoned" for that
  # allele.
  # M and m are variants of a mitochondrial gene from two different populations
  #
  #create an empty vector for gametes
  x <- vector(length=18, mode="numeric")
  # ABM
  x[1] <- (sum(female_pop[1]*(1/4),female_pop[2]*(1/8),female_pop[3]*(1/8),female_pop[4]*((1/8)*(1-r)),female_pop[5]*(1/8),
               female_pop[7]*(1/8)*r,female_pop[9]*(1/8),female_pop[10]*(1/8)*r,female_pop[13]*((1/8)*(1-r))))/(Ne/2)
  # AbM
  x[2] <- (sum(female_pop[6]*(1/4),female_pop[8]*(1/8),female_pop[14]*(1/8)))/(Ne/2)
  # aBM
  x[3] <- (sum(female_pop[11]*(1/4),female_pop[12]*(1/8),female_pop[15]*(1/8)))/(Ne/2)
  # abM
  x[4] <- (female_pop[16]*(1/4))/(Ne/2)
  
  # Add the 'poisened' gametes
  #Ab'M
  x[5] <- (sum(female_pop[2]*(1/8),female_pop[4]*r*(1/8),female_pop[5]*(1/8),female_pop[7]*((1/8)*(1-r)),female_pop[10]*((1/8)*(1-r)),
               female_pop[13]*(1/8)*r))/(Ne/2)
  
  # a'BM
  x[6] <- (sum(female_pop[3]*(1/8),female_pop[4]*(1/8)*r,female_pop[7]*((1/8)*(1-r)),female_pop[9]*(1/8),female_pop[10]*((1/8)*(1-r)),
               female_pop[13]*(1/8)*r))/(Ne/2)
  
  # ab'M
  x[7] <- (sum(female_pop[12]*(1/8),female_pop[15]*(1/8)))/(Ne/2)
  
  # a'bM
  x[8] <- (sum(female_pop[8]*(1/8),female_pop[14]*(1/8)))/(Ne/2)
  
  # a'b'M
  x[9] <- (sum(female_pop[4]*((1/8)*(1-r)),female_pop[7]*(1/8)*r,female_pop[10]*(1/8)*r,
               female_pop[13]*((1/8)*(1-r))))/(Ne/2)
  
  # ABm
  x[10] <- (sum(female_pop[17]*(1/4),female_pop[18]*(1/8),female_pop[19]*(1/8),female_pop[20]*((1/8)*(1-r)),female_pop[21]*(1/8),
                female_pop[23]*(1/8)*r,female_pop[25]*(1/8),female_pop[26]*(1/8)*r,female_pop[29]*((1/8)*(1-r))))/(Ne/2)
  # Abm
  x[11] <- (sum(female_pop[22]*(1/4),female_pop[24]*(1/8),female_pop[30]*(1/8)))/(Ne/2)
  # aBm
  x[12] <- (sum(female_pop[27]*(1/4),female_pop[28]*(1/8),female_pop[31]*(1/8)))/(Ne/2)
  # abm
  x[13] <- (female_pop[32]*(1/4))/(Ne/2)
  
  #Ab'm
  x[14] <- (sum(female_pop[18]*(1/8),female_pop[20]*r*(1/8),female_pop[21]*(1/8),female_pop[23]*((1/8)*(1-r)),
                female_pop[26]*((1/8)*(1-r)),female_pop[29]*(1/8)*r))/(Ne/2)
  
  # a'Bm
  x[15] <- (sum(female_pop[19]*(1/8),female_pop[20]*(1/8)*r,female_pop[23]*((1/8)*(1-r)),female_pop[25]*(1/8),
                female_pop[26]*((1/8)*(1-r)),female_pop[29]*(1/8)*r))/(Ne/2)
  
  # ab'm
  x[16] <- (sum(female_pop[28]*(1/8),female_pop[31]*(1/8)))/(Ne/2)
  
  # a'bm
  x[17] <- (sum(female_pop[24]*(1/8),female_pop[30]*(1/8)))/(Ne/2)
  
  # a'b'm
  x[18] <- (sum(female_pop[20]*((1/8)*(1-r)),female_pop[23]*(1/8)*r,female_pop[26]*(1/8)*r,
                female_pop[29]*((1/8)*(1-r))))/(Ne/2)
  
  # Incorporation of mutation per gamete
  y <- vector(length=18, mode="numeric")
  # gamete pool: "ABM", "AbM", "aBM", "abM", "Ab'M", "a'BM", "ab'M", "a'bM", "a'b'M", 
  #              "ABm", "Abm","aBm", "abm", "Ab'm", "a'Bm", "ab'm", "a'bm", "a'b'm"
  #ABM
  y[1] <- x[1] - U.a*x[1] - U.b*x[1] - U.a*U.b*x[1]
  #AbM
  y[2] <- x[2] - U.a*x[2] + U.b*x[1]
  #aBM
  y[3] <- x[3] - U.b*x[3] + U.a*x[1]
  #abM
  y[4] <- x[4] + U.a*x[2] + U.b*x[3] + U.a*U.b*x[1]
  #Ab'M
  y[5] <- x[5] - U.a*x[5]
  #a'BM
  y[6] <- x[6] - U.b*x[6]
  #ab'M
  y[7] <- x[7] + U.a*x[5]
  #a'bM
  y[8] <- x[8] + U.b*x[6]
  #a'b'M Nothing can mutate into this allele, as mutation into a poison allele is impossible
  y[9] <- x[9]
  #ABm
  y[10] <- x[10] - U.a*x[10] - U.b*x[10] - U.a*U.b*x[10]
  #Abm
  y[11] <- x[11] - U.a*x[11] + U.b*x[10]
  #aBm
  y[12] <- x[12] - U.b*x[12] + U.a*x[10]
  #abm
  y[13] <- x[13] + U.a*x[11] + U.b*x[12] + U.a*U.b*x[10] 
  #Ab'm
  y[14] <- x[14] - U.a*x[14]
  #a'Bm
  y[15] <- x[15] - U.b*x[15]
  #ab'm
  y[16] <- x[16] + U.a*x[14]
  #a'bm
  y[17] <- x[17] + U.b*x[15]
  #a'b'm Nothing can mutate into this allele, as mutation into a poison allele is impossible
  y[18] <- x[18]
  ## now we return this vector
  return(y)
  
}



# this fnx takes a gamete pool and reconstitutes a
makePop <- function(male_gametes, female_gametes,Ne, s, t, q){
  # first we calculate probs based on frequency of gametes in pool
  # then we just draw from a multinomial distribution
  mg <- male_gametes
  fg <- female_gametes
  z.freq <- c(#1 ABxABM
    mg[1]*(1-s-t) * fg[1]*(1-s-t),
    #2 ABxAbM + ABxAb'M
    mg[1]*(1-s-t) * fg[2]*(1-s) + mg[1]*(1-s-t) * fg[5]*(1-s),
    #3 ABxaBM + ABxa'BM
    mg[1]*(1-s-t) * fg[3]*(1-t) + mg[1]*(1-s-t) * fg[6]*(1-t),
    #4 ABxabM + ABxab’M + ABxa’bM + ABxa’b’M
    mg[1]*(1-s-t) * fg[4] + mg[1]*(1-s-t) * fg[7] + mg[1]*(1-s-t) * fg[8] + mg[1]*(1-s-t) * fg[9],
    #5 AbxABM
    mg[2]*(1-s)  * fg[1]*(1-s-t),
    #6 AbxAbM
    mg[2]*(1-s)  * fg[2]*(1-s),
    #7 AbxaBM + Abxa’BM
    mg[2]*(1-s) * fg[3]*(1-t) + mg[2]*(1-s) * fg[6]*(1-t),
    #8 AbxabM + Abxa’bM
    mg[2]*(1-s) * fg[4] + mg[2]*(1-s) * fg[8],
    #9 aBxABM
    mg[3]*(1-t) * fg[1]*(1-s-t),
    #10 aBxAbM + aBxAb’M
    mg[3]*(1-t) * fg[2]*(1-s) + mg[3]*(1-t) * fg[5]*(1-s),
    #11 aBxaBM
    mg[3]*(1-t) * fg[3]*(1-t),
    #12 aBxabM + aBxab’M
    mg[3]*(1-t) * fg[4] + mg[3]*(1-t) * fg[7],
    #13 abxABM
    mg[4] * fg[1]*(1-s-t),
    #14 abxAbM
    mg[4] * fg[2]*(1-s),
    #15 abxaBM
    mg[4] * fg[3]*(1-t),
    #16 abxabM
    mg[4] * fg[4],
    
    #17 ABxABm
    mg[5]*(1-s-t) * fg[10]*q*(1-s-t),
    #18 ABxAbm + ABxAb'm
    mg[5]*(1-s-t) * fg[11]*q*(1-s) + mg[5]*(1-s-t) * fg[14]*q*(1-s),
    #19 ABxaBm + ABxa'Bm
    mg[5]*(1-s-t) * fg[12]*q*(1-t) + mg[5]*(1-s-t) * fg[15]*q*(1-t),
    #20 ABxabm + ABxab’m + ABxa’bm + ABxa’b’m
    mg[5]*(1-s-t) * fg[13]*q + mg[5]*(1-s-t) * fg[16]*q + mg[5]*(1-s-t) * fg[17]*q + mg[5]*(1-s-t) * fg[18]*q,
    #21 AbxABm
    mg[6]*(1-s) * fg[10]*q*(1-s-t),
    #22 AbxAbm
    mg[6]*(1-s)* fg[11]*q*(1-s),
    #23 AbxaBm + Abxa’Bm
    mg[6]*(1-s) * fg[12]*q*(1-t) + mg[6]*(1-s) * fg[15]*q*(1-t),
    #24 Abxabm + Abxa’bm
    mg[6]*(1-s) * fg[13]*q + mg[6]*(1-s) * fg[17]*q,
    #25 aBxABm
    mg[7]*(1-t) * fg[10]*q*(1-s-t),
    #26 aBxAbm + aBxAb’m
    mg[7]*(1-t) * fg[11]*q*(1-s) + mg[7]*(1-t) * fg[14]*q*(1-s),
    #27 aBxaBm
    mg[7]*(1-t) * fg[12]*q*(1-t),
    #28 aBxabm + aBxab’m
    mg[7]*(1-t) * fg[13]*q + mg[7]*(1-t) * fg[16]*q,
    #29 abxABm
    mg[8] * fg[10]*q*(1-s-t),
    #30 abxAbm
    mg[8] * fg[11]*q*(1-s),
    #31 abxaBm
    mg[8] * fg[12]*q*(1-t),
    #32 abxabm
    mg[8] * fg[13]*q)
  
  z <- (rmultinom(1:32, Ne, prob = z.freq))
  return(z)
}


#Tracking allele frequency of M4 (B)
mon4.fnx <- function(results, max){
  count<-1
  x<-vector()
  for(j in seq.int(from=1, to=max)){
    x[count] <- ((2*sum(results[c(1,3,9,11,17,19,25,27), j]) + sum(results[c(2,4,5,7,10,12,13,15,18,20,21,23,26,28,29,31), j]))) / (2*sum(results[,j]))
    count<-count+1
  }
  return(x)
}

#Tracking allele frequency of M1 (A)
mon1.fnx <- function(results, max){
  count<-1
  y<-vector()
  for(s in seq.int(from=1, to=max)){
    y[count] <- ((2*sum(results[c(1,2,5,6,17,18,21,22), s]) + sum(results[c(3,4,7,8,9,10,13,14,19,20,23,24,25,26,29,30), s]))) / (2*sum(results[,s]))
    count<-count+1
  }
  return(y)
}

#Tracking allele frequency of mtDNA (M)
mtDNA.fnx <- function(results, max){
  count<-1
  m<-vector()
  for(s in seq.int(from=1, to=max)){
    m[count] <- (((sum(results[c(1:16), s]))) / (sum(results[,s])))
    count<-count+1
  }
  return(m)
}
