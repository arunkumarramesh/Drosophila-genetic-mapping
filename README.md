# Drosophila-genetic-mapping
Goal: ...

**First step is to create a genetic map using marker segreation ratios, to see how linked markers are in current cross**

Set working directory and load library

```
library(qtl)
library(zoo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```
First inspect the data
```
comb1 <- read.csv(file="qtl_chr3rest_with5cm.csv")
head(comb1)
```
|   | Phenotype | Chr2L_3cM | Chr2L_5cM | Chr2L_7cM | Chr2L_10.3cM | Chr2L_17cM | Chr2L_27cM | Chr2L_34cM | Chr2L_54cM | Chr2R_67cM | Chr2R_104cM |
|---|-----------|-----------|-----------|-----------|--------------|------------|------------|------------|------------|------------|-------------|
| 1 | NA        | 1         | 1         | 1         | 1            | 1          | 1          | 1          | 1          | 1          | 1           |
| 2 | 1         | A         | A         | H         | A            | A          | A          | A          | A          | H          | A           |
| 3 | 2         | H         | H         | H         | H            | H          | H          | H          | H          | H          | A           |
| 4 | 1         | A         | A         | A         | A            | A          | A          | A          | A          | A          | A           |
| 5 | 2         | A         | <NA>      | A         | A            | A          | A          | A          | A          | A          | A           |
| 6 | 2         | H         | H         | H         | H            | H          | H          | H          | A          | A          | H           |

Read genotypes from cross into rqtl
```
mapthis <- read.cross("csv", file="qtl_chr3rest_with5cm.csv", estimate.map=FALSE)
```
only keep individuals typed for more than nine markers
```
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>9)) 
```
estimates recombination fraction (rf)
```
mapthis <- est.rf(mapthis)
```
assign markers into linkage groups
```
lg <- formLinkageGroups(mapthis, max.rf=0.5, min.lod=6) 
```
get rf
```
rf <- pull.rf(mapthis) 
```
get lod scores
```
lod <- pull.rf(mapthis, what="lod")
```
plot rf v lod scores
```
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/9a4bb452-f06a-4ba6-ab6c-fd9864652144)
calculate error probabilities
```
mapthis <- calc.errorlod(mapthis, error.prob=0.01)
```
estimate map with new error probabilities, to see how genetic map changes with greater genotyping error
```
newmap <- est.map(mapthis, error.prob=0.01) 
plotMap(newmap,show.marker.names = TRUE)
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/9667c306-ec35-4f8c-ae71-5e8de0ba178e)

replace old map with new one
```
mapthis <- replace.map(mapthis, newmap)
```

run this code to use standard drosophila genetic map. 

column names contain recombination map position. use recombination estimates from drosophila genetic map (from Flybase) rather than those calculated here using marker ratios.
```
markers=colnames(read.csv("qtl_chr3rest_with5cm.csv")) 
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
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/090f89f6-249b-4f8c-af71-3c9ce7d6a675)

Flybase map used further
```
mapthis <- replace.map(mapthis, flymap)
```
Plot the genotypes.  Genotypes AA and AB are indicated by white and black circles.  X's in intervals having a crossover.
```
plotGeno(mapthis, chr=1, ind=c(1:115))
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/c01cae8f-f617-443c-8d68-6451203385c9)

calculate genotype probabilities along map
```
mapthis <- calc.genoprob(mapthis, step=1, error.prob=0.05)
```
Genome scan with single qtl model using Haley-Knott (HK) regression. Easy and fast, but performs poorly in the case of selective genotyping (Feenstra et al. 2006), which is not the case here.
```
out.hk <- scanone(mapthis, method="hk")
summary(out.hk)
```
|          | chr | pos | lod  |
|----------|-----|-----|------|
| c1.loc12 | 1   | 15  | 17.9 |

here carrying out a permutation test. LOD thresholds (1000 permutations)
```
operm.hk <- scanone(mapthis, method="hk", n.perm=1000) 
summary(operm.hk, alpha=0.05)
```
|    | lod  |
|----|------|
| 5% | 1.46 |

Test if marker is above log threshold. Yes it is!
```
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)
```
|          | chr | pos | lod  | pval |
|----------|-----|-----|------|------|
| c1.loc12 | 1   | 15  | 17.9 | 0    |


get confidence intervals
```
cM_in_qtl=out.hk$pos[out.hk$lod>(max(out.hk$lod)-1.5)]
low=min(cM_in_qtl)
up=max(cM_in_qtl)
low
```
11
```
up
```
25

Composite interval mapping, see if more than one qtl exists
```
out.cim <- cim(mapthis,window=10)
```
Combined plot with Haley-Knott single qtl model (Black line) and composite interval mapping (Red line) which blue shaded region indicating area of significance.
```
plot(out.hk,xlab="Chromosome 2 Map Position (cM)",cex.lab=1.4, cex.axis=1.4)
plot(out.cim,add=T,col="red")
rect(low,-5,up,20, col= rgb(0,0,1.0,alpha=0.15), border = NA)
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/9a59fb08-c13c-4824-90fe-b46715bbd68e)

Idetified a single 24cM QTL from broad mapping.

Next fine mapping.

This is the oberved risk ratio from non-recombinant flies (adults with capsules) (same  marker genotype at 3 and 27cM). There were 738 flies RR genotype and 260 SS flies
```
observedRR=906/282
```
function to simulate recombinant genotypes
```
genosim=function(nrecombs=298, risk=observedRR){
  
  #number genetically resistant recombinants
  nR=rbinom(1, nrecombs, risk/(risk+1))
  
  #recombinants at 3 and 27 cM, so 240 0.1cM intervals
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
    
  }
  #use this for chi-squared ci
  c(table(chiCIgene)[2]/length(chiCIgene),
    mean(chiCIsize)
  )
}
```
run simulations

test different chi-squared drops. called lods.
```
lods=40:60/10
effect_L=matrix(nrow=length(lods),ncol=2)

for(k in 1:length(lods)){
  effect_L[k,]=genosim2(N=251, R=observedRR,chidrop=lods[k])
  print(k)
}
effect_L2=cbind(lods,effect_L[,1:2])
colnames(effect_L2)[2:3]=c("proportion_qtl_containing_gene","mean_qtl_size")
write.csv(effect_L2,file="QTL simulations different chi drops.csv")

effect_L2=read.csv(file="QTL simulations different chi drops.csv")

#pdf(file="QTL simulations different chi drops.pdf",height=4,width=5)
par(mfrow=c(1,1))
plot(100*effect_L2[,3]~effect_L2[,2],
     ylim=c(0,9),
     xlab="chi2 drop",
     ylab="% times gene outside 95% CI interval",
     main='risk ratio=3.21,N=298')
abline(h=5,col="red",lwd=3)
```
![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/be220203-f1a7-41cb-873e-5a18f0fe598a)


This is the drop in chi2 to get confidence intervals on location. Derived empirically from simulation script.
```
chi_drop=4.6
```
Read input file
```
HSSH <- read.csv("HS-SH,15-5.csv")
head(HSSH)
```
| Sample | Plate | Well | Row | Col | Egg lay date | Egg lay period | Egg transfer (ET) | Infection start date | Infection end date | Infection period | Infection start date - Egg lay date | No. vials made on that date | Yeast usage       | Cage | Total flies collected on that date | 3 cM | 27 cM | Genotype (3 and 27 cM) | Plate Label | 7 cM | 10.3 cM | 12 cM | 17 cM | Wasp primer melting temperature  | Wasp DNA amplified | 8 cM | 11 cM | 10.7 cM | 11.3 cM | 11.6 cM | 8.5 cM | 9 cM | 10 cM |
|--------|-------|------|-----|-----|--------------|----------------|-------------------|----------------------|--------------------|------------------|-------------------------------------|-----------------------------|-------------------|------|------------------------------------|------|-------|------------------------|-------------|------|---------|-------|-------|----------------------------------|--------------------|------|-------|---------|---------|---------|--------|------|-------|
| 15     | 1     | G2   | G   | 2   | 25. Jan      | Overnight      | 25. Jan           | 28. Jan              | 29.01.19           | 1 day            | 3                                   | 10                          | Plastic, no yeast | 2    | 15                                 | H    | S     | HS                     | 1_G2        | H    | H       | S     | S     | 80.82924652                      | N                  | NA   | H     | NA      | H       | H       | NA     | NA   | NA    |
| 18     | 1     | B3   | B   | 3   | 27. Jan      | Overnight      | 27. Jan           | 29. Jan              | 31.01.19           | 2 days           | 2                                   | 70                          | Plastic, no yeast | 2    | 89                                 | H    | S     | HS                     | 1_B3        | H    | H       | S     | S     | 81.12721252                      | N                  | NA   | H     | NA      | H       | H       | NA     | NA   | NA    |
| 34     | 1     | B5   | B   | 5   | 27. Jan      | Overnight      | 27. Jan           | 29. Jan              | 31.01.19           | 2 days           | 2                                   | 70                          | Plastic, no yeast | 2    | 89                                 | H    | S     | HS                     | 1_B5        | H    | H       | S     | S     | 81.87213135                      | N                  | NA   | S     | NA      | NA      | NA      | NA     | NA   | NA    |
| 50     | 1     | B7   | B   | 7   | 27. Jan      | Overnight      | 27. Jan           | 29. Jan              | 31.01.19           | 2 days           | 2                                   | 70                          | Plastic, no yeast | 2    | 89                                 | H    | S     | HS                     | 1_B7        | H    | H       | S     | S     | 78.44551086                      | Y                  | NA   | H     | NA      | H       | H       | NA     | NA   | NA    |

Create matrix of marker genotypes
```
genotypes1=HSSH[,grep("cM",names(HSSH))]
names(genotypes1)=gsub(".cM", "", names(genotypes1))
genotypes1 <- genotypes1[-c(3)]
names(genotypes1)=gsub("X", "", names(genotypes1))
ind=order(as.numeric(names(genotypes1)))
genotypes1=data.matrix(genotypes1[,ind])
```
Only use samples with evidence of parasitic wasp infection
```
genotypes1 <- genotypes1[HSSH$Wasp.DNA.amplified=="Y",]
```
Impute missing genotypes when flanking values same
```
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

genotypes2=t(apply(genotypes1,1,na.approx))
ind=rowSums(!t(apply(genotypes2,1,is.wholenumber)))==0
genotypes=genotypes2[ind,]
table(genotypes[,7] == genotypes[,9])
```
| FALSE | TRUE |
|-------|------|
| 2     | 109  |

only 2 informative recombinants between 10.3 and 11

Test if genotypes differ from 50:50 ratio
```
genotypes[genotypes==2]=0

chi2=vector()
for(i in 1:ncol(genotypes)){
  x=prop.test(sum(genotypes[,i]),nrow(genotypes),p=0.5)$statistic
  #next line changes sign depending on whether in correct direction
  chi2[i]=ifelse((sum(genotypes[,i])/nrow(genotypes)>0.5),
                 x,-x)
}
```
get location CIs
```
locations=as.numeric(colnames(genotypes1))
outside=which(chi2<max(chi2-4.6))
peak=which(chi2==max(chi2))
low=locations[outside[max(which(outside<peak))]]
up=locations[outside[min(which(outside>peak))]]
low
```
10
```
up
```
11.6

Plot of chi2 (test of whether marker ratio departs from 50:50 ratio) along chromosome. Single 1.6cM peak detected.

Beware this is not interval mapping so pay little attention to line between markers

![image](https://github.com/arunkumarramesh/Drosophila-genetic-mapping/assets/23363383/22c54fca-e76b-4c85-a023-392ab1605d32)


