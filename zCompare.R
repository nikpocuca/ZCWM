

zeroCompare <- function() {
  
  
  
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
  
  
  
  
}