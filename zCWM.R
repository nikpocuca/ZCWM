

# ZERO INFLATED CLUSTER WEIGHTED POISSON MODEL
# June 13th - 2017 
# Creator: Nik Pocuca 
# The following is code written using two libraries and their respective dependencies.
# If you had not yet installed them, please do so now, please run the following function. 

initFunctions <- function(){

  install.packages('flexCWM')
  install.packages('pscl')
  
  
}

# initFunctions()

# IMPORTATION OF LIBRARIES. 
library(flexCWM)
library(pscl)


load_book <- function(){
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

#CWM Addons. Not found in flexCWM package
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
zcwm <- function(data, np){
  
  
  #initializations of zero space.
  data.z <<- data
  data.z$NB <<- as.integer(data.z$NB <= 0) 
  
  
  print("Beginning Partitioning using Cwm")
  #First Partition of space. (Poisson) 
  cwm_poisson <<- cwm(formulaY = fregzi,
                      data= data,
                      familyY = poisson(link="log"),
                      Xnorm = cbind(DriverAge,CarAge,Density),
                      modelXnorm = "EVV",
                      k = 1:np)
  
  
  
  #Declare global for second partition. 
  declare_g(data.z)
  #Second Partition of space. (Binomial)
  cwm_binomial <<- cwm(formulaY = fregzi,
                       data = data.z,
                       Xnorm = cbind(DriverAge,CarAge,Density),
                       familyY = binomial(link = "logit"),
                       modelXnorm = "VVI",
                       k = 1:np)
  
  
  nPara <- length(getParGLM(cwm_binomial)[1])
  
  c_pois <<- getCluster(cwm_poisson)
  c_bin <<- getCluster(cwm_binomial)
  lex <<- paste(c_pois, c_bin, sep= "")
  test <<- cbind(c_pois,c_bin,lex)
  partitions <<- match(lex,unique(lex))
  cuts <<- data.frame(c_pois,c_bin,lex,partitions)
  

  
  # RESTRUCTURING COEFFICIENTS
  #| =============================================================================================|
  #| RETURN VECTORS FUNCTION                                                                      |
  #| Nik Pocuca July 19th - 2017                                                                  |
  #| Returns glm vectors in the form of a dataframe.                                              |
  #| =============================================================================================|
  returnVectors <- function(cwmGLM){
  
    # Get names 
    vNames <-  names(cwmGLM$GLMComp.1$coefficients)
  
    
                                      # Placeholder dataframe. 
                                      dataPlaceholder <-  data.frame()  
                                      for(i in cwmGLM){
                                      dataPlaceholder <- rbind(dataPlaceholder, i$coefficients)
                                      }
                                
  
  # Set Names for dataframe
  colnames(dataPlaceholder) <- vNames
  
  return(dataPlaceholder)} # END OF RETURN VECTORS FUNCTION 
  #| =============================================================================================|

  poissonVectors <- returnVectors(pois_coef)
  binomialVectors <- returnVectors(bin_coef)
  
  
  #Begin gluing partitions through zero inflated poisson. 
  data_space <<- cbind(data,cuts)
  
  # Message: July 19th 2017
  # Create a function that will generate dataspaces, with their respective vectors based on the lexicon. 
  
  # CLASS CREATION OF LEXICON VECTORS
  #| =============================================================================================|
  #| LEXICON VECTOR CLASS                                                                         |
  #| Nik Pocuca July 21th - 2017                                                                  |
  #| Definition of a Lexicon vector class, an LVC is a class that contains two vectors            |
  #| and the associated lexicon.                                                                  |
  #| =============================================================================================|
  lexVector <- setClass(Class = "Lexicon Vector",
                        slots = c(lexicon = "character",pVector = "vector",bVector = "vector"))
  
  #Set a method for the class
  setMethod("$", "Lexicon Vector", function(x, name) {
              slot(x, name)
            })
  # END OF LEXICON VECTOR CLASS 
  #| =============================================================================================|
  
  # GENERATING VECTOR OBJECTS
  #| =============================================================================================|
  #| GENVEC FUNCTION                                                                              |
  #| Nik Pocuca July 20th - 2017                                                                  |
  #| Takes in a lexicon, generates an object   with names as numbers of the corresponding vectors |
  #| =============================================================================================|
  genVec <- function(tag, pVec, bVec){
        vSplitTag <- unlist(strsplit(tag, ""))
        
        vecObject <- lexVector(lexicon = tag, 
                              pVector = t(pVec[vSplitTag[1],]),
                              bVector = t(pVec[vSplitTag[2],]))
   
 
    
    
  return(vecObject)
  } # END OF GENVEC FUNCTION
  #| =============================================================================================|
  
  # GENERATING DATASPACES
  #| =============================================================================================|
  #| GENSPACE FUNCTION                                                                            |
  #| Nik Pocuca July 20th - 2017                                                                  |
  #| Creates a data space of each cut in an object so you can pass in the data, and vectors.      |
  #| =============================================================================================|
  genData <- function(dSpace, pVectors, bVectors){
      lexicon <- as.character(unique(data_space$lex))
      
      genObject <- c()
      for(i in lexicon){
        
                        genVecObj <- genVec(tag = i, 
                                     pVec = pVectors,
                                     bVec = bVectors)
                        
                        
                        
                        # END OF FOR LOOP
                        genObject <- c(genObject,genVecObj)
      }
    
      
      
      
      genAns <- c()
      for(i in unique(data_space$partitions)){

    
      genAns <- c(genAns ,assign(paste('subspace',i,sep=''), 
             subspace(data= data_space[data_space$partitions == i,], vectors = genObject[[i]] ) #assigning value
             ))
      
        
      }
      
    #  print(genAns)
  return(genAns)
  } #END OF GENSPACE FUNCTION 
  #| =============================================================================================|
  
  
  
  # CLASS CREATION OF DATA SPACE
  #| =============================================================================================|
  #| DATA SPACE CLASS                                                                      |
  #| Nik Pocuca July 21th - 2017                                                                  |
  #| Definition of entire data space, will contain subclass of a data space which conatins k      |
  #| partitions.                                                                                  |
  #|                                                                                              |
  #| =============================================================================================|
  PSPACE <- setClass(Class = "PSPACE",
                     slots = c(
                     k = "numeric",
                     spaces = "list"
                     ))
  
  
  
  # SETTING METHOD FOR ACCESSING INFORMATION
  setMethod("$", "PSPACE", function(x, name) { 
    slot(x, name)
  })
  
  
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
  
  
  
  
  # CLASS CREATION OF SUBSPACE
  #| =============================================================================================|
  #| SUBSPACE CLASS                                                                               |
  #| Nik Pocuca July 21th - 2017                                                                  |
  #| Definition of subspace of dataspace class. Each subspace conatins the dataset with the       |
  #| referenced partition, and a coupled lexicon vector.                                          |
  #|                                                                                              |
  #| =============================================================================================|
  subspace <- setClass(Class = "subspace",
                     slots = c(
                       data = "data.frame",
                       vectors = "Lexicon Vector"
                     ))
  
  
  
  # SETTING METHOD FOR ACCESSING INFORMATION
  setMethod("$", "subspace", function(x, name) { 
    slot(x, name)
  })
  
  
  #| =============================================================================================|
  
  
  
  partSpace <- genData(dSpace = data_space, 
          pVectors = poissonVectors,
          bVectors = binomialVectors)
  # PARTITIONED SPACE IS COMPLETED 
  
 
  

  #Generate object of models to be passed in. 
  data_space_1 <<- data_space[data_space$partitions == 1,]
  
  zipModel_1 <<- list(
  data = data_space_1,
  p_vector = pois_split[data_space_1$c_pois[1]],
  b_vector = bin_split[data_space_1$c_bin[1]]
  )
  
  
  data_space_2 <<- data_space[data_space$partitions == 2,]
  zipModel_2 <<- list(
    data = data_space_2,
    p_vector = pois_split[data_space_2$c_pois[1]],
    b_vector = bin_split[data_space_2$c_bin[1]]
  )
  
  data_space_3 <<- data_space[data_space$partitions == 3,]
  zipModel_3 <<- list(
    data = data_space_3,
    p_vector = pois_split[data_space_3$c_pois[1]],
    b_vector = bin_split[data_space_3$c_bin[1]]
  )
  
  data_space_4 <<- data_space[data_space$partitions == 4,]
  zipModel_4 <<- list(
    data = data_space_4,
    p_vector = pois_split[data_space_4$c_pois[1]],
    b_vector = bin_split[data_space_4$c_bin[1]]
  )
  
  data_space_5 <<- data_space[data_space$partitions == 5,]
  zipModel_5 <<- list(
    data = data_space_5,
    p_vector = pois_split[data_space_5$c_pois[1]],
    b_vector = bin_split[data_space_5$c_bin[1]]
  )
  
  data_space_6 <<- data_space[data_space$partitions == 6,]
  zipModel_6 <<- list(
    data = data_space_6,
    p_vector = pois_split[data_space_6$c_pois[1]],
    b_vector = bin_split[data_space_6$c_bin[1]]
  )
  
  #Begin the Zero Inflated Poisson Here. 
  
  result <- list(
    space = data,
    space.z = data.z$NB,
    cwm_bin = cwm_binomial,
    cwm_pois = cwm_poisson)
  class(result) <- "zcwm"
  
  return(result)
} 
#End of zcwm function

zero <- zcwm(c_24)




zero_1 <- zeroinfl(formula = fregzi, 
                     data = zipModel_1$data,
                     dist = "poisson")


zero_2 <- zeroinfl(formula = fregzi, 
                   data = zipModel_2$data,
                   dist = "poisson")


zero_3 <- zeroinfl(formula = fregzi, 
                   data = zipModel_3$data,
                   dist = "poisson")

zero_4 <- zeroinfl(formula = fregzi, 
                   data = zipModel_4$data,
                   dist = "poisson")

zero_5 <- zeroinfl(formula = fregzi, 
                    data = zipModel_5$data,
                    dist = "poisson")

zero_6 <- zeroinfl(formula = fregzi, 
                   data = zipModel_6$data,
                   dist = "poisson")

summary(zero_1)
summary(zero_2)
summary(zero_3)
summary(zero_4)
summary(zero_5)
summary(zero_6)

#mine

zeroNik_1 <- z_mk1(formula = fregzi, 
                   zipModelE = zipModel_1,
                   data = zipModel_1$data,
                   dist = "poisson")
zeroNik_2 <- z_mk1(formula = fregzi,
                   zipModelE = zipModel_2,
                   data = zipModel_2$data,
                   dist = "poisson")
zeroNik_3 <- z_mk1(formula = fregzi,
                   zipModelE = zipModel_3,
                   data = zipModel_3$data,
                   dist = "poisson")

zeroNik_4 <- z_mk1(formula = fregzi, 
                   zipModelE = zipModel_4,
                   data = zipModel_4$data,
                   dist = "poisson")

zeroNik_5 <- z_mk1(formula = fregzi, 
                   zipModelE = zipModel_5,
                   data = zipModel_5$data,
                   dist = "poisson")

zeroNik_6 <- z_mk1(formula = fregzi, 
                   zipModelE = zipModel_6,
                   data = zipModel_6$data,
                   dist = "poisson")


summary(zeroNik_1)
summary(zeroNik_2)
summary(zeroNik_3)
summary(zeroNik_4)
summary(zeroNik_5)
summary(zeroNik_6)









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










