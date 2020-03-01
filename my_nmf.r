# Implementation of Non-Negative Matrix Factorization
# Theorem 1: Euclidian Distance
# Theorem 2: Kullback-Leibler Divergence

# Dataset taken from
# https://www.kaggle.com/lewisduncan93/the-economic-freedom-index




library(ggplot2)

# Clear environment and set current working directory
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# create test matrices
A <- matrix(c(1:9),byrow=T, nrow=3)
B <- matrix(c(9:1),byrow=T, nrow=3)

A
B


### v ~~ WH
### THEOREM 1
##### UPDATE RULES
theorem1 <- function(V, r) {
  
  ## euclidian distance init
  euclDistSqd <- function(A, B) {
    result <- sum((A - B)**2)
    return(result)
  }
  
  ## initialisation of n & m
  n <- nrow(V)
  m <- ncol(V)
  
  ## initialisation of W & H
  W <- matrix(c(runif(n*r)), n, r)
  H <- matrix(c(runif(r*m)), r, m)
  
  distance <- sqrt(euclDistSqd(V, W%*%H))
  
  while (TRUE) {
  ## updating H
  updateFactorH <- (t(W) %*% V) / (t(W) %*% W %*% H)
  H <- H * updateFactorH
  
  ## updating W
  updateFactorW <- (V %*% t(H)) / (W %*% H %*% t(H))
  W <- W * updateFactorW
  
  newDistance <- sqrt(euclDistSqd(V, W%*%H))
  
  #print(distance)
  #print(newDistance)
  
  # distance comparison
  if (newDistance==distance){
    break
  }
  #if (newDivergence < 1e-20) {
  #  break
  #}
  
  else {
    distance <- newDistance
  }
  }
  
  returnObj <- list("W"=W, "H"=H)
  return(returnObj)
}

result <- theorem1(A, 2)
result

# check the result
result$W %*% result$H





### v ~~ WH
### THEOREM 2
##### UPDATE RULES
theorem2 <- function(V, r) {
  
  ## kullback-leiber
  kulLeiDiv <- function(A, B) {
    result <- sum(A * log(A / B) - A + B)
    return(result)
  }
  
  ## initialisation of n & m
  n <- nrow(V)
  m <- ncol(V)
  
  ## initialisation of W & H
  W <- matrix(c(runif(n*r)), n, r)
  H <- matrix(c(runif(r*m)), r, m)
  
  
  divergence <- kulLeiDiv(V, W%*%H)
  #for (iter in c(1:100)) {
  
  while (TRUE) {
    
    ## updating H
    updateFactorH = matrix(, nrow = r, ncol = m)
    for (a in c(1:r)){
      for (mu in c(1:m)){
        updateFactorH[a, mu] <- (sum(W[,a]* (V[,mu]/(W%*%H)[,mu]))) /
                                  sum(W[,a])
      }
    } 
    
    H <- H * updateFactorH

    ### updating W
    updateFactorW = matrix(, nrow = n, ncol = r)
    for (i in c(1:n)){
      for (a in c(1:r)){
        updateFactorW[i, a] <- (sum(H[a,]* V[i,] / (W %*% H)[i,])) / sum(H[a,])
      }
    } 
    
    
    W <- W * updateFactorW
    
    
    newDivergence <- kulLeiDiv(V, W%*%H)
    
    #cat("old divergence: ", divergence, "\n")
    #cat("new divergence: ", newDivergence, "\n\n")
    # distance comparison
    if (newDivergence==divergence){
      break
    }
    #if (newDivergence < 1e-10) {
    #  break
    #}
    else {
      divergence <- newDivergence
      #print(divergence)
  }
  }
  returnObj <- list("W"=W, "H"=H)
  return(returnObj)
}



result <- theorem2(A, 2)
result

# check the result
result$W %*% result$H



#### APPLICATION
data <- read.csv("data/economic_freedom_index2019_data.csv")

# set country name as index
row.names(data) <- data$Country.Name

# drop unnecessary columns
data[c('Country.Name', 'WEBNAME', 'Region', 'Country')] <- c(rep(NULL, 4))

# remove dollar signs and thousands separator
data$GDP..Billions..PPP. <- gsub('\\$', '', data$GDP..Billions..PPP. )
data$GDP.per.Capita..PPP. <- gsub('\\$', '', data$GDP.per.Capita..PPP.)
data$GDP.per.Capita..PPP. <- gsub(',', '', data$GDP.per.Capita..PPP.)

# convert all columns to numeric
data[,c(2:ncol(data))] <- sapply(data[,c(2:ncol(data))], as.numeric)

# drop nas
data <- data[complete.cases(data),]



selectedCountries <-  c("Bahrain","Paraguay","Colombia", "Mozambique","Mongolia",
                        "Angola", "Portugal", "Tanzania", "Vanuatu", "Cyprus","Netherlands",
                        "Montenegro", "Luxembourg","Austria", "Djibouti", "Finland", "Singapore",
                        "Greece", "Cambodia")

df <- data[selectedCountries, c(2:ncol(data))]



### Comparison
result <- theorem2(t(df), 2)
result

#result <- theorem1(t(df), 2)
#result


plottingDf <- as.data.frame(t(result$H))
colnames(plottingDf) <- c("Factor1", "Factor2")
plottingDf <- scale(plottingDf)
plottingDf <- as.data.frame(plottingDf)
plottingDf$Country <- rownames(df)


a <- ggplot(plottingDf, aes(x = Factor1, y = Factor2))
a <- a + geom_point()
a <- a + geom_text(aes(label=Country),hjust=1, vjust=1, size=3)
a <- a + labs(title = "NMF Result",
              subtitle = "Kullback Leibler Implementation"
              )
a <- a + xlim(-3, 3)
a <- a + ylim(-3, 3)
a


