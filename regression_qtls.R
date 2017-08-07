#!/usr/bin/R


#====================================================================================
#
#        AUTHOR: Ashley Tehranchi, tehranchi576@gmail.com 
#  ORGANIZATION: Stanford University
#
#====================================================================================

linearRegression <-
function(pd, gt, tmppfx, indiv) {

cat("Regression... "); flush.console()
options("scipen"=100)
n.indi = indiv
postTotal = read.table(pd, sep="\t", header=T)
print(pd)
genofile = matrix(scan(file = gt, sep = "\t", what = double(), nlines = n.indi), byrow = TRUE, nrow = n.indi)
genoTotal = 1 - (0.5 * genofile)

# restrict to MAF cutoff
MAF= 0.02
a = t(genoTotal)
b= cbind(postTotal,a)
mafsub = b[apply(b[,16:75], MARGIN = 1, function(x) mean(x) >= MAF & mean(x) <= 0.98), ]
post = mafsub[,1:15]
genotypes = as.matrix(mafsub[,16:75])
postTemp = mafsub

n.snps = nrow(post)
#genotypes=as.matrix(t(genos))
weights = rep(1/n.indi,n.indi)
real.props = genotypes %*% weights
depths=  post$Depth
postProp = post$POSTfreq
estimated.props = postProp
G = genotypes[,-n.indi] - genotypes[,n.indi]
Y = postProp - genotypes[,n.indi]
eps = 0.0001
# force postfreq to be between 0.0001 and 0.9999
props.for.weights = pmin(1-eps,pmax(estimated.props,eps))
# weight = depth / (adjusted.post * (1-adjusted.post))
regression.weights = depths / (props.for.weights * (1-props.for.weights) )
good = which( (postProp>0.1) & (postProp < 0.9))
m = lm(Y[good] ~ G[good,]-1,weights=regression.weights[good])  ## run without intercept
coefs = m$coef
s = summary(m)
cov.mat = s$cov.unscaled * s$sigma^2
big.cov.mat = matrix(NA,n.indi,n.indi)
big.cov.mat[-n.indi,-n.indi] = cov.mat
big.cov.mat[n.indi,n.indi] = sum(cov.mat)
big.cov.mat[n.indi,-n.indi] = big.cov.mat[-n.indi,n.indi] = -rowSums(cov.mat)
vars = sapply(1:n.snps, function(i) genotypes[i,] %*% big.cov.mat %*% genotypes[i,])
vars[vars < 0] <- 0 
all.coefs = c(coefs, 1-sum(coefs))
preProps = genotypes %*% all.coefs 
preProps[preProps > 1] <- 1 
preVars = vars 
postProps = estimated.props
postProps[postProps > 1] <- 1 
postVars = 1/regression.weights
postVars[postVars > 1] <- 1 
Zs = (postProps - preProps)/sqrt(preVars +postVars)
p.values = 2*(1-pnorm(abs(Zs)))

analytic.pv = function(preProp,preVar,postProp,depth){
    possible.props = 0:depth / depth
    extreme.props = (0:depth)[abs(possible.props - preProp) >= abs(postProp-preProp)]
    if(preVar < 0.0001){
      return(sum(dbinom(extreme.props,depth,preProp)))
    }
    f = function(x){
      denom = pnorm(1,mean=preProp,sd=sqrt(preVar)) - pnorm(0,mean=preProp,sd=sqrt(preVar))
      dnorm(x,mean=preProp,sd=sqrt(preVar)) * sum(dbinom(extreme.props,depth,x)) / denom
    }
    integrate(f,0,1)$value
}

an.pvs<-sapply(1:n.snps, function(i) analytic.pv(preProps[i],preVars[i],postProps[i],depths[i]))
post$prechipfreq <- preProps
post$pvalue <- an.pvs
post$SNPpostfreq <- 1 - post$POSTfreq 
post$SNPprefreq <- 1 - post$prechipfreq


combined = cbind(post, genotypes)
sig = subset(combined, pvalue <= 0.05)
nonsig = subset(combined, pvalue > 0.05)

####### LD calc for sig
r2.sig = vector("numeric",length=nrow(sig))
holdnum = nrow(sig)-1
for (i in 1:holdnum){
    if (sig[i,2] - sig[i+1,2] < 200000) {
          r2.sig[i] <- cor(t(sig[i,20:79]), t(sig[i+1,20:79]), method="pearson")
  }
  else{ r2.sig[i] <- 0}
}
sig$r2 <- (r2.sig)^2
ldfileOUTsig = sig[, c("Chr", "position", "r2")]
ldfileOUTsigFinal = ldfileOUTsig[complete.cases(ldfileOUTsig),]
totaloutsig = paste(tmppfx, ".sigLD", sep="")
write.table(ldfileOUTsigFinal, file= totaloutsig, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)

####### LD calc for nonsig
r2.nonsig = vector("numeric",length=nrow(nonsig))
holdnum = nrow(nonsig)-1
for (i in 1:holdnum){
    if (nonsig[i,2] - nonsig[i+1,2] < 200000) {
          r2.nonsig[i] <- cor(t(nonsig[i,20:79]), t(nonsig[i+1,20:79]), method="pearson")
  }
  else{ r2.nonsig[i] <- 0}
}
nonsig$r2 <- (r2.nonsig)^2
ldfileOUTnonsig = nonsig[, c("Chr", "position", "r2")]
ldfileOUTnonsigFinal = ldfileOUTnonsig[complete.cases(ldfileOUTnonsig),]
totaloutnonsig = paste(tmppfx, ".nonsigLD", sep="")
write.table(ldfileOUTnonsigFinal, file= totaloutnonsig, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)


#cluster.sig = data.frame(Chr=character(),start=numeric(), end=numeric(), stringsAsFactors=FALSE) 
#holdcount = nrow(ldfileOUTsig)-1
#for (i in 1:holdcount){
#    if (ldfileOUTsig[i,3] < 0.8) {
#      cluster.sig[i,1] = ldfileOUTsig[i,1]
#      cluster.sig[i,2] = ldfileOUTsig[i,2] 
#      cluster.sig[i,3] = ldfileOUTsig[i,2] 
#    }
#    else {
#    
#    }
#

# Summary plots
png(paste(tmppfx, ".summaryPlots.png", sep=""), width=1200, height=1200, pointsize=40)
par(mfrow=c(2,2))
plot(table(round(p.values,d=2)))
qqplot(p.values,runif(n.snps), )#,xlim=c(0,0.1),ylim=c(0,0.1))
abline(0,1)
plot(table(round(an.pvs,d=2)))
qqplot(an.pvs,runif(n.snps))#,xlim=c(0,0.1),ylim=c(0,0.1))
abline(0,1)
dev.off()

png(paste(tmppfx, ".qqplot.png", sep=""), width=1200, height=1200, pointsize=40)
qqplot(runif(n.snps), an.pvs, xlab="Expected", ylab="Observed")#,xlim=c(0,0.1),ylim=c(0,0.1))
abline(0,1)
dev.off()

## remove SNPs with pre =0 or 1
print(dim(post))
temp = subset(post, prechipfreq>0 & prechipfreq<1)
## must be sorted by pvalue
x = temp[with(temp, order(pvalue)), ]
print(dim(x))
xs = x
options("scipen"=100)
xs$start <- xs$position - 1
bed = xs[,c("Chr", "start", "position")]  ## reorder by column numbers
sortedBed = bed[with(bed,order(Chr, start, position)),]
finaloutFile = paste(tmppfx, ".total.txt", sep="")
finalbed = paste(tmppfx, ".bed", sep="")
write.table(x, file= finaloutFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
write.table(sortedBed, file=finalbed, sep="\t", quote=FALSE, row.names=FALSE, col.names=F)

}

linearRegression("postdataIN", "genosoutName", "prefix", "numIndiv")
