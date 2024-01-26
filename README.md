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


#get confidence intervals
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


#this is the drop in chi2 to get confidence intervals on location. Derived empirically from simulation script.
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
```
