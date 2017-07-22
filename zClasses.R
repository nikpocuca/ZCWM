
# ZCLASSES
# July 22nd - 2017
# Creator: Nik Pocuca
# The following is are class definitions for objects that will be used by the zcwm function
# in the main call.




# CLASS CREATION OF LEXICON VECTORS
#| =============================================================================================|
#| LEXICON VECTOR CLASS                                                                         |
#| Nik Pocuca July 21th - 2017                                                                  |
#| Definition of a Lexicon vector class, an LVC is a class that contains two vectors            |
#| and the associated lexicon.                                                                  |
#| =============================================================================================|
lexVector <- setClass(Class = "Lexicon Vector",
slots = c(lexicon = "character",
          pVector = "vector",
          bVector = "vector"))

#Set a method for the class
setMethod("$", "Lexicon Vector", function(x, name) {
    slot(x, name)
})
# END OF LEXICON VECTOR CLASS
#| =============================================================================================|







# CLASS CREATION OF DATA SPACE
#| =============================================================================================|
#| DATA SPACE CLASS                                                                             |
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
#  END OF DATASPACE CLASS
#| =============================================================================================|







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

#  END OF SUBSPACE CLASS
#| =============================================================================================|












