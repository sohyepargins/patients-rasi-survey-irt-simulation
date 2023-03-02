setwd("/Users/sohyepark/Rfiles/rfiles/edms787/project/acculturative-survey-irt-simulation")
library(plyr)
require(rjags)
library(sfsmisc) # For PPC 
#-------------------------------------------------------------------------------
# Graded model
#-------------------------------------------------------------------------------
# rasi data
rasi <- read.csv("rasionly2.csv", header = T)
#rasi<-data.matrix(data.frame(lapply(rasi, ordered, levels= c("Never", "Rarely", "Sometimes", "A lot", "Always"))))
n <- nrow(rasi)
m <- ncol(rasi)
K<-apply(rasi, 2, max)

2.5^(-2)
0.5^(-2)# SD = .5 <=> Tau = 4

# model fitting
init <- list( 
              list('kappa.star' = matrix(rnorm(14*4), 14,4),
                  .RNG.name = "base::Mersenne-Twister", .RNG.seed = 498),
              list('kappa.star' = matrix(rnorm(14*4), 14,4),
                   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 502),
              list('kappa.star' = matrix(rnorm(14*4), 14,4),
                   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 504)
              )
kappa <- t(sapply(1:m, function(j) c(rep(NA, K[j]-1), rep(0.0, max(K) - K[j]))))

#kappa.star <- matrix(rnorm(14 * 4), 14, 4)
#kappa <- t(apply(kappa.star, 1, sort) )
#init <- list( 
#              list('kappa' = kappa,
 #                  .RNG.name = "base::Mersenne-Twister", .RNG.seed = 499),
  #                .RNG.name = "base::Mersenne-Twister", .RNG.seed = 500),
   #           list('kappa' = kappa,
    #               .RNG.name = "base::Mersenne-Twister", .RNG.seed = 501) 
     #         )
# mod.grm2 <- jags.model("grm.jags",
#                       data = list('Y' = rasi, 
#                                   'n' = n,
#                                   'm' = m,
#                                   'K' = K,
#                                   'kappa' = kappa),
#                       init=init,
#                       n.chains = 3, n.adapt = 2000)
# 
# update(mod.grm2, n.iter = 4000)#burn in
# result.grm2 <- coda.samples(mod.grm2, 
#                            variable.names = c('alpha', 'theta', 'kappa'), 
#                            n.iter = 500, thin = 1)
# summary(result.grm2)



### Differen priors for alpha and kappa: mod.grm3

mod.grm4 <- jags.model("grm4.jags",
                       data = list('Y' = rasi, 
                                   'n' = n,
                                   'm' = m,
                                   'K' = K,
                                   'kappa' = kappa),
                       init=init,
                       n.chains = 3, n.adapt = 4000)

update(mod.grm4, n.iter = 10000) #burn in
result.grm4 <- coda.samples(mod.grm4, 
                            variable.names = c('alpha', 'theta', 'kappa'), 
                            n.iter = 10000, thin = 1)


acfplot(result.grm4[,51])
gelman.plot(result.grm4)  
gelman.diag(result.grm4) #PSRF
effectiveSize(result.grm4)  


summary(result.grm4)
plot(result.grm4[,77])
plot(result.grm4[, "theta[1]"])

#Plot
plot(result.grm4, density =F)
plot(result.grm4[,"theta"])
plot(result.grm4[[1]])
plot(result.grm4[[1]][,"theta"])



# Posterior predictive model checks
post.draws1 = as.matrix(result.grm4)


kappalist <- function(v) {
  ## Creates a list of kappa values from a vector of BUGS output
  p <- max(as.numeric(gsub("(kappa\\[|,[[:digit:]]+\\])", "", grep("kappa", names(v), value=TRUE))))
  lapply(1:p, function(j) v[grep(paste("kappa\\[", j, ",[[:digit:]]+\\]", sep=""), names(v))])
}
newdata <- function(theta, alpha, kappa){
  ## Generates a new data set from a graded response model
  ## theta = vector (length n) of "ability" values
  ## alpha = vector (length p) of discrimination parameters
  ##         or a single value (vector of length 1)
  ## kappa = a list of kappa values
  n <- length(theta)
  p <- length(kappa)
  Ystar <- t(sapply(theta, function(th) rlogis(p, alpha*th, 1.0)))
  Y <- t(apply(Ystar, 1, function(y) sapply(seq(p), function(j) cut(y[j], breaks=c(-Inf, kappa[[j]], Inf), labels=FALSE))))
  dimnames(Y) <- NULL
  return(Y)
}
ppktau <- function(sims){
  require(plyr)
  parnames <- colnames(sims)
  dumfun <- function(v){
    theta <- v[grep("theta", parnames)]
    alpha <- v[grep("alpha", parnames)]
    kappa <- kappalist(v)
    p <- length(kappa)
    out <- cor(newdata(theta, alpha, kappa), method="kendall")[lower.tri(diag(p))]
    return(out)
  }
  out <- adply(sims, 1, dumfun, .progress="text")
  return(out[,-1])
}
densityShade <- function(x, from, to, col=rgb(0.5, 0.5, 0.5, 0.5), border=NA, ...){
  den <- density(x)
  from <- max(from, min(den$x))
  to <- min(to, max(den$x))
  den <- density(x, from=from, to=to)
  x <- den$x
  y <- den$y
  xx <- c(x, rev(x))
  yy <- c(y, rep(0, times=length(x)))
  polygon(xx, yy, col=col, border=border, ...)
}
makeNames <- function(name, ..., symmetric.matrix=FALSE, diag=TRUE){
  if (symmetric.matrix){
    idx <- which(upper.tri(diag(max(...)), diag=diag), arr.ind=TRUE)
  } else {
    idx <- expand.grid(rev(list(...)))
  }
  name <- paste(name, "[", sep="")
  if (ncol(idx)>1){
    for (j in ncol(idx):2){
      name <- paste(name, idx[,j], ",", sep="")
    }
  }
  name <- paste(name, idx[,1], "]", sep="")
  return(name)
}
plotpppval <- function(sims, true, xlim=NULL){
  if (!is.matrix(sims))
    sims <- as.matrix(sims)
  npar <- ncol(sims)
  ##mult.fig(npar, mar=c(3, 1.0, 1.0, 1.0), oma=rep(0, 4))
  mult.fig(npar, mar=c(2.5, 0.25, 0.25, 0.25), oma=rep(0, 4))
  den <- apply(sims, 2, density)
  if (is.null(xlim)){
    xl <- range(as.vector(sapply(den, function(d) range(d$x))))
    xl <- xl + c(-1, 1)*0.04*diff(xl)
  }
  else xl <- xlim
  ylmax <- max(sapply(den, function(d) max(d$y)))
  yl <- c(0, 1.04*ylmax)
  for(i in 1:npar){
    plot(density(sims[,i]), xlim=xl, ylim=yl, tcl=-0.2, mgp=c(1, 0.1, 0.0), xlab=colnames(sims)[i], ylab="", bty="n", yaxt="n", main="")
    if (mean(true[i]<sims[,i])>0.5){
      from <- -Inf
      to <- true[i]
      pval <- formatC(mean(true[i]>sims[,i]), digits=2, format="f")
    }
    else{
      from <- true[i]
      to <- Inf
      pval <- formatC(mean(true[i]<sims[,i]), digits=2, format="f")
    }
    densityShade(sims[,i],from=from, to=to)
    text(xl[1], 0.75*ylmax, labels=pval, pos=4)
    ##abline(v=true[i], col="red")
  }
}

## Tau correlations and PPC plot
pp.eq1<-ppktau(post.draws1)
colnames(pp.eq1)<-makeNames("ktau", 1:14, 1:14, symmetric.matrix = TRUE, diag = FALSE)
ktau1<-cor(rasi, method="kendall")[lower.tri(diag(m))]
plotpppval(pp.eq1,ktau1)

#caterplot: doesn't seem to be useful
dev.new()
library(mcmcplots)
caterplot(pp.eq, collapse=FALSE)
caterpoints(ktau,pch="x", col ="red")