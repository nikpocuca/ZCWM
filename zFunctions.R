

# ZFUNCTIONS
# July 22nd - 2017
# Creator: Nik Pocuca
# The following are function definitions for outputs and calls that will be used by the zcwm
# function.






#| RESTRUCTURING COEFFICIENTS
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









#| GENERATING VECTOR OBJECTS
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












