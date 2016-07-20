setwd("~/R Codes/Asset Allocation")

library(gdata);library(moments);library(pracma);library(forecast);library(MASS);library(tseries);
library(fGarch);library(fArma);library(Matrix);library(matrixcalc);

######################### Function fitting best GARCH model by minimizing AIC #########################

garchSearch = function(
  ts,  
  minOrder,
  maxOrder,
  trace=FALSE, dstn)
{
  bestAic = 1e9
  len = nrow( ts ) 
  for( p in minOrder[1]:maxOrder[1] ) for( q in minOrder[2]:maxOrder[2] )
  {
    if( p == 0 && q == 0 )
    {    
      next
    }    
    formula = as.formula( paste( sep="", "formula ~ garch(", p, ",", q, ")" ) )
    fit = tryCatch( garchFit( formula, data=ts, cond.dist= dstn ),
                    error=function( err ) FALSE,
                    warning=function( warn ) FALSE )
    if( !is.logical( fit ) )
    {    
      fitAic = fit@fit$ics[1]
      if( fitAic < bestAic )
      {    
        bestAic = fitAic
        bestFit = fit
        bestModel = c( p, q )
      }    
      if( trace )
      {    
        ss = paste( sep="", "(", p, ",", q, "): AIC = ", fitAic )
        print( ss ) 
      }    
    }    
    else
    {    
      if( trace )
      {    
        ss = paste( sep="", "(", p, ",", q, "): None" )
        print( ss ) 
      }    
    }    
  }
  
  if( bestAic < 1e9 )
  {
    return( list( aic=bestAic, fit=bestFit, model=bestModel ) )
  }
  
  return( FALSE )
}

######################################################################################################


###################################### Risk Estimating Function ######################################
risk.fun <- function(data, Sigma, w, a){
  
  D <- dim(data)[2]
  m<-v<-nu<-rep(0,D)
  
  for (j in 1:D)
    {  # multivariate Student-t estimated parameters
      fit <- sstdFit(data[,j]) 
      m[j] <- fit$estimate[1] 
      v[j] <- fit$estimate[2]
      nu[j] <- fit$estimate[3]
    }
  
  n <- 90000  # number of samples
  mu <- rep(0,D)
  
  # Generate a multivariate distribution following the centered t distribution with nu degrees of freedom and covariance matrix Sigma
  X <- mvrnorm(n, mu, Sigma)
  grouped_t <- U <- matrix(rep(0,n*D),nrow=n)
  
  for (i in 1:n)
    {
    for (j in 1:D)
      {
      grouped_t[i,j] = (X[i,j]+m[j])/sqrt(rchisq(1,nu[j])/nu[j]) # multivariate t
      U[i,j] <- dt(grouped_t[i,j], df=nu[j]) # Uniforms
    }
  }
  
  fitc <- sstdFit(rowSums(grouped_t*w)); mc <- fitc$estimate[1]; vc <- fitc$estimate[2]; nuc <- fitc$estimate[3]
  
  C <- dt(qt(U,nuc), nuc)  # Student-t Copula
  
  # CVaR of simulated scenarios
  returns <- rowSums(grouped_t*w); returns = sort(returns, decreasing = TRUE)
  VaR = abs(quantile(returns,1-a))
  CVaR <- abs(sum(returns[ceiling(n*a):length(returns)])/(n-ceiling(n*a))) 
  
  # CVaR of historical returns
  hist_returns <- rowSums(data*w); hist_returns = sort(hist_returns, decreasing = TRUE)
  hist_VaR = abs(quantile(hist_returns,1-a))
  hist_CVaR <- abs(sum(hist_returns[ceiling(nrow(data)*a):length(hist_returns)])/(nrow(data)-ceiling(nrow(data)*a))) 
  
  # Max Drawdown
  max_drawdown <- abs(min(hist_returns))
  
  results <- c(VaR, CVaR, hist_VaR, hist_CVaR, max_drawdown) 
  return(results)
}

####################################################################################################






SPTR = FNERTR = SPGSCITR = GDDUEAFE = LBUSTRUU = LF98TRUU = GDLEEGF = USO = SX5E = data.frame()
datalist <- list(FNERTR, GDDUEAFE, GDLEEGF, LBUSTRUU, LF98TRUU, SPGSCITR, SPTR, SX5E, USO)

for (i in 1:9){
  datalist[[i]] <- read.xls("Data_04_15.xlsx", sheet = i, header = TRUE, sep = ',')
  datalist[[i]] <- na.omit(datalist[[i]])  
}

# Store as Compound Returns 
data <- sapply(sapply(datalist, function(x){x[1:(dim(x)[1]-1),2] = diff(log(x[,2]))}), function(x){x = x[-length(x)[[1]]]})

#merge returns into one data frame which contains returns for all days in common
merged.data <- Reduce(function(...) merge(..., by ="Date"), datalist)
#order by date
merged.data <- merged.data[order(as.Date(merged.data$Date, format="%Y-%M-%D"), decreasing = FALSE),]
#omit rows containing NA's
merged.data <- merged.data[complete.cases(merged.data),]
#rename merged data frame
names(merged.data) <- c("Date", "FNERTR", "GDDUEAFE", "GDLEEGF", "LBUSTRUU", "LF98TRUU", "SPGSCITR", "SPTR", "SX5E", "USO")
merged.data <- apply(merged.data[,2:10],2,function(x)diff(log(x)));


######################### Exploratory Data Analysis ###########################

# Kurtosis/Skewness/Normality
kurtosis <- sapply(data, kurtosis); kurtosis 
skewness <- sapply(data, skewness); skewness
sapply(data, shapiro.test)

# Autocorrelation/ Partial Autocorrelation
sapply(data, acf)
sapply(data,function(x)pacf(abs(x)))

# Hurst Exponent
hurst <- sapply(data, function(x){hurstexp(x, d=50, display=TRUE);print("________________________________")})


######################### Fitting Time Series Model ###########################

# Portfolio Weights (FNERTR, GDDUEAFE, GDLEEGF, LBUSTRUU, LF98TRUU, SPGSCITR, SPTR, SX5E, USO)
w <- c(0.02, 0.02, 0.2115654059, 0.1860097137, 0.25, 0.02, 0.06961727, 0.02 , 0.20280760760)
#w <- c(rep(1/9,9)) # equally weighted
#w <- c(0,0,0,.4,0,0,.6,0,0) # traditional 60/40
ts <- rowSums(merged.data*w)
# EDA
kurtosis(ts);skewness(ts);acf(ts);pacf(ts);hurstexp(ts, d=50, display=TRUE)  # Skewed Student-t
order.param <- garchSearch(ts, minOrder=c(0,0), maxOrder=c(5,5), trace=FALSE, "sstd")
p <- order.param$model[1]; q <- order.param$model[2]
garch.fit <- garchFit(formula=~garch(1,1),data=ts,trace=F,cond.dist="sstd")


######################### Robust Covariance Matrix ###########################
cov.pca <- prcomp(cov(merged.data[,1:ncol(merged.data)]), center = TRUE, scale. = TRUE, retx = TRUE) 
summary(cov.pca) # run diagnosis

# clean up covariance matrix
pc.use <- 4 # explains 99% of variance
# reconstruct matrix function
prcomp.recon <- function(pca, pcs){
  if(is.null(pcs)) pcs <- seq(pca$sdev)
  recon <- as.matrix(pca$x[,pcs]) %*% t(as.matrix(pca$rotation[,pcs]))
  if(pca$scale[1] != FALSE){
    recon <- scale(recon , center=FALSE, scale=1/pca$scale)
  }
  if(pca$center[1] != FALSE){
    recon <- scale(recon , center=-pca$center, scale=FALSE)
  }
  recon
}

trunc <- prcomp.recon(cov.pca, c(1:pc.use))

# Force positive definiteness
Sigma <- as.matrix(nearPD(trunc, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE,
                doSym = FALSE, doDykstra = TRUE, only.values = FALSE,
                ensureSymmetry = !isSymmetric(trunc),
                eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                maxit = 100, conv.norm.type = "I", trace = FALSE)$mat)


######################### Risk Measures ###########################

# level of confidence for CVaR
a=.95

risk <- risk.fun(merged.data, Sigma, w, a)

writeLines(paste0(a*100,"% VaR: ",round(risk[1]*100,2)," % \n",a*100,"% CVaR: ",round(risk[2]*100,2),
                  " % \n", a*100,"% Historical VaR: ",round(risk[3]*100,2)," % \n", a*100,
                  "% Historical CVaR: ",round(risk[4]*100,2)," % \n", "Max.Drawdown: ",
                  round(risk[5]*100,2)," % \n "))