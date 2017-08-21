

# ZERO INFLATED CLUSTER WEIGHTED POISSON MODEL
# June 13th - 2017 
# Creator: Nik Pocuca 
# The following is code written using two libraries and their respective dependencies.
# If you had not yet installed them, please do so now, run the following function.

initLibrary <- function(){

  install.packages('flexCWM')
  install.packages('pscl')
  
  
}

# initFunctions()

# IMPORT OF LIBRARIES.
library(flexCWM)
library(pscl)


load_book <- function(){
    
  CONTRACTS.f <- read.csv(file = "CONTRACTS.csv")
  #BOOK#
  # factor(CONTRACTS.f$BRAND=="F",labels=c("other","F"))
  CONTRACTS.f$powerF <- factor(1*(CONTRACTS.f$POWER%in%letters[4:6])+ 2*(CONTRACTS.f$POWER%in%letters[7:8]),labels=c("other","DEF","GH"))
  CONTRACTS.f$powerF <- factor(1*(CONTRACTS.f$POWER%in%letters[4:6])+ 2*(CONTRACTS.f$POWER%in%letters[7:8]),labels=c("other","DEF","GH"))
  CONTRACTS.f$GAS <- CONTRACTS.f$Gas
  CONTRACTS.f$EXPOSURE <- CONTRACTS.f$Exposure
  CONTRACTS.f$NB  <- CONTRACTS.f$ClaimNb
  CONTRACTS.f <<- CONTRACTS.f
}
declare_g <- function(data){
  AGECAR <<- data$AGECAR
  AGEDRIVER <<- data$AGEDRIVER
  Brand <<- data$Brand
  BRAND <<- data$BRAND
  brandF <<- data$brandF
  CarAge <<- data$CarAge
  ClaimNb <<- data$ClaimNb
  Density <<- data$Density
  DENSITY <<- data$DENSITY
  DriverAge <<- data$DriverAge
  Exposure <<- data$Exposure
  EXPOSURE <<- data$EXPOSURE
  Gas <<- data$Gas 
  GAS <<- data$GAS
  NB <<- data$NB 
  Power <<- data$Power
  POWER <<- data$POWER
  powerF <<- data$powerF
  Region <<- data$Region
}
create24 <- function() {
  c_24 <<- CONTRACTS.f[CONTRACTS.f$Region=="R24",]
}

#CWM Addons, not found in CWM package



getFitted <- function(object, ...){
  best <- getBestModel(object,...)
  obj  <- best$models[[1]]
  if (!is.null(obj$GLModel)){
    lr <- lapply(seq_len(obj$k), function(i){
      par <- obj$GLModel[[i]]
      c(list(fitted=par$model$fitted.values),par[-1])
    })
    names(lr) <- paste0("GLMComp.",seq_len(obj$k))
    lr
  } else NULL
}
getResiduals <- function(object, ...){
  best <- getBestModel(object,...)
  obj  <- best$models[[1]]
  if (!is.null(obj$GLModel)){
    lr <- lapply(seq_len(obj$k), function(i){
      par <- obj$GLModel[[i]]
      c(list(resid=par$model$residuals),par[-1])
    })
    names(lr) <- paste0("GLMComp.",seq_len(obj$k))
    lr
  } else NULL
}



#Initialize
load_book()
create24()
declare_g(c_24)
fregzi <- NB ~ DriverAge + CarAge + Density + powerF



#Beginning of zwm function 
zcwm <- function(data, formulaZP, np){
  
  
  #initializations of zero space.
  data.z <<- data
  data.z$NB <<- as.integer(data.z$NB <= 0) 
  
  
  cat('Beginning Partitioning using Cwm')
  cat('--------------------------------')
  cat(' ')
  
  #First Partition of space. (Poisson) 
  
  dclareCWM <- function(){
    cwm_poisson <<- cwm(formulaY = formulaZP,
                      data= data,
                      familyY = poisson(link="log"),
                      Xnorm = cbind(DriverAge,CarAge,Density),
                      modelXnorm = "EVV",
                      k = 1:np)
  
  
  
  #Declare global for second partition. 
  declare_g(data.z)
  #Second Partition of space. (Binomial)
  
  
  cwm_binomial <<- cwm(formulaY = formulaZP,
                       data = data.z,
                       Xnorm = cbind(DriverAge,CarAge,Density),
                       familyY = binomial(link = "logit"),
                       modelXnorm = "VVI",
                       k = 1:np)
  
  }
  
  
  
  dclareCWM()
  
  
  
  nPara <- length(getParGLM(cwm_binomial)[1])
  
  c_pois <<- getCluster(cwm_poisson)
  c_bin <<- getCluster(cwm_binomial)
  lex <<- paste(c_pois, c_bin, sep= "")
  test <<- cbind(c_pois,c_bin,lex)
  partitions <<- match(lex,unique(lex))
  cuts <<- data.frame(c_pois,c_bin,lex,partitions)
  

  

  poissonVectors <- returnVectors(pois_coef)
  binomialVectors <- returnVectors(bin_coef)
  
  
  #Begin gluing partitions through zero inflated poisson. 
  data_space <<- cbind(data,cuts)
  
  # Message: July 19th 2017
  # Create a function that will generate dataspaces, with their respective vectors based on the lexicon. 
 
  
  
  # Message:  July 21st, 2017 - Nik Pocuca 
  # I need to create a method for appending subspaces based on the number of unique partitions. 
  # Should be straight forward:
  # 1. create a subspace class,
  # 2. create  a method for appending k subspaces, to the p space classes.
  # ____________________________________________________________________________________
  # 3. create a function that will extract each subspace, run the zip_mk1.R program 
  #              on each of them and store it's results in another class called zswap
  # ____________________________________________________________________________________
  # In the end, I will pass one giant class into a function that will iteratively run 
  
  
  
  partSpace <- genData(dSpace = data_space, 
          pVectors = poissonVectors,
          bVectors = binomialVectors)
  # PARTITIONED SPACE IS COMPLETED 
  
  
  
  counter <- 1
  result <<- c()
  for(i in partSpace){
  
    
    
  print(                 )
  print(                  )
  print(paste("Now attempting maximization number",counter))
  print("---------------------------")
    
  
  
    
 tryCatch({    holdModel <-  z_mk1(formula = formulaZP, 
          zipModelE = i,
          data = i$data,
          dist = "poisson")    
  
  print(summary(holdModel))
    
  
  counter <- counter + 1
  assign(paste('zeroNikT', counter,sep=''), holdModel)
 },error = function(){})
  
  
  
  
  }
  
  

  #Begin the Zero Inflated Poisson Here. 
  class(result) <- "zcwm"
  
  return(result)
} 
#End of zcwm function




runZ <- function(){

zeroTEST <- zcwm(data = c_24,
             formulaZP = fregzi,
             np = 3)

}









# APPENDIX
# Below are adjustments or extra code that I used at one point but now dont use at all. 

#pull out fitted values for poisson
#p_f <<- c(getFitted(cwm_poisson))
#p1_f <<- as.vector(p_f$GLMComp.1$fitted)
#p2_f <<- as.vector(p_f$GLMComp.2$fitted)
#  p3_f <<- as.vector(p_f$GLMComp.3$fitted)
#fit_poisson <<- cbind(p1_f,p2_f)

#pull out fitted values for zero
#b_f <<- c(getFitted(cwm_binomial))
#b1_f <<- as.vector(b_f$GLMComp.1$fitted)
#b2_f <<- as.vector(b_f$GLMComp.2$fitted)
#b3_f <<- as.vector(b_f$GLMComp.3$fitted)
#fit_binomial <<- cbind(b1_f,b2_f)


#intialize vectors
#fit <<- c()
#fit_z <<- c()

#print("Matching fitted values....")
#match fittted values for both of them.
#for(i in 1:length(fit_binomial[,1])){
#matcher_z <- fit_binomial[i,c_bin[i]]
#matcher <- fit_poisson[i,c_pois[i]]
#  fit <<- as.vector(append(fit,matcher))
#  fit_z <<- as.vector(append(fit,matcher_z))
#}
#note this takes a long time. 


#fit_cuts <<- cbind(fit,fit_z)







