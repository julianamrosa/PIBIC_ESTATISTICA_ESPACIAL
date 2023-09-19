#Discretizando y e padronizando x
georgia_nb_std <- georgia_data
georgia_nb_std$PctBach <- as.integer(georgia_nb_std$PctBach)

#TotPop90
georgia_nb_std$TotPop90 <- (georgia_nb_std$TotPop90-
                               mean(georgia_nb_std$TotPop90))/
  sd(georgia_nb_std$TotPop90)
#PctRural
georgia_nb_std$PctRural <- (georgia_nb_std$PctRural-
                              mean(georgia_nb_std$PctRural))/
  sd(georgia_nb_std$PctRural)
#PctEld
georgia_nb_std$PctEld <- (georgia_nb_std$PctEld-
                              mean(georgia_nb_std$PctEld))/
  sd(georgia_nb_std$PctEld)
#PctFB
georgia_nb_std$PctFB <- (georgia_nb_std$PctFB-
                            mean(georgia_nb_std$PctFB))/
  sd(georgia_nb_std$PctFB)
#PctPov
georgia_nb_std$PctPov <- (georgia_nb_std$PctPov-
                           mean(georgia_nb_std$PctPov))/
  sd(georgia_nb_std$PctPov)
#PctBlack
georgia_nb_std$PctBlack <- (georgia_nb_std$PctBlack-
                            mean(georgia_nb_std$PctBlack))/
  sd(georgia_nb_std$PctBlack)

#Teste 10
startTime <- Sys.time()
mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov",
              "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
              METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic", MODEL="poisson")
endTime <- Sys.time()
#17.36 mins

#Teste 11 ou Teste 17
startTime <- Sys.time()
mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov",
              "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
#23.7 mins

#Teste 16
startTime <- Sys.time()
mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov",
              "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="cv", MODEL="poisson")
endTime <- Sys.time()
#27.48 mins

#Teste 19
startTime <- Sys.time()
mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov",
              "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="FIXED_G", BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
#1.07 horas

### Testes do Artigo ###

#gwnbr
startTime <- Sys.time()
out <- mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld",
              "PctFB", "PctPov", "PctBlack"), LAT="Y",
       LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="negbin", MGWR="no")
endTime <- Sys.time()
#1.05 mins

#gwpr
startTime <- Sys.time()
out <- mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
              XVAR=c("TotPop90", "PctRural", "PctEld",
                     "PctFB", "PctPov", "PctBlack"), LAT="Y",
              LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
              BANDWIDTH="aic", MODEL="poisson", MGWR="no")
endTime <- Sys.time()
#1.8 mins

#mgwnbr
startTime <- Sys.time()
out <- mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
       XVAR=c("TotPop90", "PctRural", "PctEld",
              "PctFB", "PctPov", "PctBlack"), LAT="Y",
       LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
#17.07 mins

#mgwpr
startTime <- Sys.time()
out <- mgwnbr(DATA=georgia_nb_std, YVAR="PctBach",
              XVAR=c("TotPop90", "PctRural", "PctEld",
                     "PctFB", "PctPov", "PctBlack"), LAT="Y",
              LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
              BANDWIDTH="aic", MODEL="poisson")
endTime <- Sys.time()
#14.66 mins

#logistic
startTime <- Sys.time()
mgwnbr(DATA=logistic, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="LOGISTIC", H=159)
endTime <- Sys.time()
