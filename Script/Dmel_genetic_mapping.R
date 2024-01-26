library(qtl)
library(zoo)

# Goal: To map Drosophila melanogaster gene region(s) containing resistance to infection by parasitic wasp Leptopilina boulardi. Previous controlled crosses identified a resistance factor or chromosome 2. First step was broad mapping on chromosome 2. Infected 385 F2 D. melanogaster flies from an backross with a resistant and suceptible parental genotype with L. boulardi and recorded their resistance status as resistant or susceptible (1 or 2 in phenotype column below). Selected 10 insertion-deletion markers interspersed along chromosome 2 for mapping and genotyped the 385 D. melanogaster. Then, interval mapping was done using R/QTL.

## Broad mapping ##

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# First inspect the data
comb1 <- read.csv(file="qtl_chr3rest_with5cm.csv")
head(comb1)

## read genotypes from cross into rqtl
mapthis <- read.cross("csv", file="qtl_chr3rest_with5cm.csv", estimate.map=FALSE)

### initally creating a genetic map using marker segreation ratios, to see how linked markers are in current cross
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>9)) ### only keep individuals typed for more than nine markers
mapthis <- est.rf(mapthis) ## estimates recombination fraction (rf)
lg <- formLinkageGroups(mapthis, max.rf=0.5, min.lod=6) ### assign markers into linkage groups
rf <- pull.rf(mapthis) ### get rf
lod <- pull.rf(mapthis, what="lod") ### get lod scores
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score") ## plot rf v lod scores
mapthis <- calc.errorlod(mapthis, error.prob=0.01) ### calculate error probabilities
newmap <- est.map(mapthis, error.prob=0.01) ### estimate map with new error probabilities, to see how genetic map changes with greater genotyping error
plotMap(newmap,show.marker.names = TRUE)
mapthis <- replace.map(mapthis, newmap) ### replace old map with new one

#run this code to use standard drosophila genetic map
markers=colnames(read.csv("qtl_chr3rest_with5cm.csv")) ## column names contain recombination map position
markers2=markers[2:length(markers)]
markers2=gsub('Chr2L_','',markers2)
markers2=gsub('Chr2R_','',markers2)
markers2=gsub('cM','',markers2)
markers2=as.numeric(markers2)
flymap=as.list(markers2) 
flymap=newmap
flymap[[1]]<-markers2
names(flymap[[1]])=markers[2:length(markers)]
plotMap(flymap,show.marker.names = TRUE)
mapthis <- replace.map(mapthis, flymap) ### use recombination estimates from drosophila genetic map rather than calculating here using marker ratios. used further
##############################

plotGeno(mapthis, chr=1, ind=c(1:115)) # Plot the genotypes. Genotypes AA and AB are indicated by white and black circles. X's in intervals having a crossover.
mapthis <- calc.genoprob(mapthis, step=1, error.prob=0.05) ### calculate genotype probabilities along map
out.hk <- scanone(mapthis, method="hk") ## Genome scan with single qtl model using Haley-Knott (HK) regression. Easy and fast, but performs poorly in the case of selective genotyping (Feenstra et al. 2006), which is not the case here.
summary(out.hk)
operm.hk <- scanone(mapthis, method="hk", n.perm=1000) ### here carrying out a permutation test
summary(operm.hk, alpha=0.05)
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE) # Test if marker is above log threshold. Yes it is!

#get confidence intervals
cM_in_qtl=out.hk$pos[out.hk$lod>(max(out.hk$lod)-1.5)]
low=min(cM_in_qtl)
up=max(cM_in_qtl)
low
up

#composite interval mapping, see if more than one qtl exists
out.cim <- cim(mapthis,window=10)

## Combined plot with Haley-Knott single qtl model (Black line) and composite interval mapping (Red line) which blue shaded region indicating area of significance.
#pdf(file="qtl.pdf",height=5,width=5)
plot(out.hk,xlab="Chromosome 2 Map Position (cM)",cex.lab=1.4, cex.axis=1.4)
plot(out.cim,add=T,col="red")
rect(low,-5,up,20, col= rgb(0,0,1.0,alpha=0.15), border = NA)
#dev.off()
#quartz.save(type = "pdf", file="qtl.pdf",height=4,width=4)

# Idetified a single 24cM QTL from broad mapping.


## Fine mapping ## 

## Here the approach was use the same backcross as before but only genotype resistant flies and identify regions where the marker ratio departed from a 50:50 expectation (homozygous vs heterozygous). The broad mapping results from before was used to focus on the region from 3 to 27 cM on chromosome 2. The first step was to identify flies with a recombination breakpoint within 3 to 27 cM. ~1300 F2 flies were phenotyped and genotyped. But only ~300 had a breakpint between the 3 to 27 cM.

#this is the observed risk ratio
#from non-recombinant flies (adults with capsules) (same  marker genotype at 3 and 27cM)
#thre were 738 flies RR genotype and 260 SS flies
observedRR=906/282

## Choosing informative markers involved employing a χ2 drop, which was determined by simulating 1,000 datasets using the estimated risk ratio from nonrecombinant flies and the observed recombination fraction. This χ2 drop delineated a region containing the gene in 95% of the simulations.
genosim=function(nrecombs=298, risk=observedRR){
  
  #number genetically resistant recombinants
  nR=rbinom(1, nrecombs, risk/(risk+1))
  
  #recomninants at 3 and 27 cM, so 240 0.1cM intervals
  #1,0 recombinants
  mat=matrix(1,nrow=1000,ncol=240)
  breakpoint=sample(238, 1000,replace=T)+1
  for (i in 1:1000){
    mat[i,breakpoint[i]:240]=0
  }
  
  #0,1 recombinants
  mat2=matrix(0,nrow=1000,ncol=240)
  breakpoint=sample(238, 1000,replace=T)+1
  for (i in 1:1000){
    mat2[i,breakpoint[i]:240]=1
  }
  
  #shuffle the two
  mat3=rbind(mat,mat2)[sample(1:2000),]
  
  #split into resistant and susceptible at position 70 (ie 7cm from start, 10cM) 
  #assign R as genotype 1
  matR=mat3[mat3[,70]==1,]
  matS=mat3[mat3[,70]!=1,]
  
  #sample resistant and susceptible lines following risk ratio
  genotypes=rbind(matR[1:nR,],matS[1:(nrecombs-nR),])
  genotypes
}



genosim2=function(N=251, R=observedRR,ldrop=1.5,chidrop=4.5){
  
  QTL=vector()
  CIsize=vector()
  CIgene=vector()
  lodCIsize=vector()
  lodCIgene=vector()
  chiCIsize=vector()
  chiCIgene=vector()
  
  
  for(i in 1:10000){
    genos=genosim(nrecombs=N, risk=R)
    x=which(colSums(genos)==max(colSums(genos)))
    #note this is to randomly choose marker when ties
    QTL[i]=x[sample(length(x))][1]
    
    #get chi2 scores
    chi2=vector()
    for(j in 1:ncol(genos)){
      x=prop.test(sum(genos[,j]),nrow(genos),p=0.5)$statistic
      #next line changes sign depnding on whether in correct direction
      chi2[j]=ifelse((sum(genos[,j])/nrow(genos)>0.5),
                     x,-x)
    }
    
    ci_chi=which(chi2>(max(chi2)-chidrop))
    chiCIsize[i]=(max(ci_chi)-min(ci_chi))/10
    chiCIgene[i]=ifelse(min(ci_chi)<=70&max(ci_chi)>=70,"gene in CI","gene not in CI")
    
    #convert to lod 
    #            lod=chi2/(2*log(10))
    ##            ci=which(lod>(max(lod)-ldrop))
    #            lodCIsize[i]=(max(ci)-min(ci))/10
    #            lodCIgene[i]=ifelse(min(ci)<=70&max(ci)>=70,"gene in CI","gene not in CI")
    
  }
  #use this for lod drop conf intervals
  #    c(table(CIgene)[2]/length(CIgene),
  #    mean(CIsize),
  #    table(lodCIgene)[2]/length(lodCIgene),
  #    mean(lodCIsize)
  #    )
  #use this for chi-squared ci
  c(table(chiCIgene)[2]/length(chiCIgene),
    #      mean(chiCIsize),
    #      table(chiCIgene)[2]/length(chiCIgene),
    mean(chiCIsize)
  )
}
###############run simulations

#########################
##test different chi-squared drops. called lods as old script!

lods=40:60/10
effect_L=matrix(nrow=length(lods),ncol=2)

for(k in 1:length(lods)){
  effect_L[k,]=genosim2(N=251, R=observedRR,chidrop=lods[k])
  print(k)
}
effect_L2=cbind(lods,effect_L[,1:2])
colnames(effect_L2)[2:3]=c("proportion_qtl_containing_gene","mean_qtl_size")
#write.csv(effect_L2,file="QTL simulations different chi drops.csv")

effect_L2=read.csv(file="QTL simulations different chi drops.csv")

#pdf(file="QTL simulations different chi drops.pdf",height=4,width=5)
par(mfrow=c(1,1))
plot(100*effect_L2[,3]~effect_L2[,2],
     ylim=c(0,9),
     xlab="chi2 drop",
     ylab="% times gene outside 95% CI interval",
     main='risk ratio=3.21,N=298')
abline(h=5,col="red",lwd=3)

#dev.off()

#this is the drop in chi2 to get confidence intervals on location. Derived empirically from simulation script.
chi_drop=4.6


HSSH <- read.csv("HS-SH,15-5.csv")
head(HSSH)

#create matrix of marker genotypes
genotypes1=HSSH[,grep("cM",names(HSSH))]
names(genotypes1)=gsub(".cM", "", names(genotypes1))
genotypes1 <- genotypes1[-c(3)]
names(genotypes1)=gsub("X", "", names(genotypes1))
ind=order(as.numeric(names(genotypes1)))
genotypes1=data.matrix(genotypes1[,ind])


##only samples with wasp 1 amplification
genotypes1 <- genotypes1[HSSH$Wasp.DNA.amplified=="Y",]


#impute missing genotypes when flanking values same
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

genotypes2=t(apply(genotypes1,1,na.approx))
ind=rowSums(!t(apply(genotypes2,1,is.wholenumber)))==0
genotypes=genotypes2[ind,]

table(genotypes[,7] == genotypes[,9]) ## only 2 informative recombinants between 10.3 and 11

#test if genotypes differ from 50:50
genotypes[genotypes==2]=0

chi2=vector()
for(i in 1:ncol(genotypes)){
  x=prop.test(sum(genotypes[,i]),nrow(genotypes),p=0.5)$statistic
  #next line changes sign depending on whether in correct direction
  chi2[i]=ifelse((sum(genotypes[,i])/nrow(genotypes)>0.5),
                 x,-x)
}

#get location CIs

locations=as.numeric(colnames(genotypes1))
outside=which(chi2<max(chi2-4.6))
peak=which(chi2==max(chi2))
low=locations[outside[max(which(outside<peak))]]
up=locations[outside[min(which(outside>peak))]]

low
up


#pdf(file="qtl mapping in adults.pdf",height=3.8,width=4)
par(mar=c(5.1,4.3,2.1,1.2))

# Plot of chi2 (test of wheather marker ratio departs from 50:50 ratio) along chromosome. Single 1.6cM peak detected.
plot(chi2~locations,type="l",
     xlab="Chromosome 2 Map Position (cM)",
     ylab=expression(chi^2),
     las=1,bty="l",cex.lab=1.3, cex.axis=1.3)
points(chi2~locations,type="p",pch=20)
rect(low,-15,up,40, col= rgb(0,0,1.0,alpha=0.15), border = NA)
#remove this line: use for guessing where to put markers
#abline(h=max(chi2)-chi_drop,col="red",lwd=3)
#beware this is not interval mapping so pay little attention to line between markers

#dev.off()