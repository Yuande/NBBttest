# NBBttest
## Introduction

NBBttest is R package for implementing negative binomial beta t-test (or called NBBt-test). NBBt-test is a statistical method for differential analysis of multiple types of RNA-seq data. NBBt-test is based on negative binomial distribution with parameters p and r: $x \backsim NB\left(p,r\right)$ where r is number of RNA sequences failed to be sequenced in RNA population (or RNA library) and p is  proportion of the sequenced RNA species in RNA library and follows beta distibution with parameters $\alpha$ and $\beta$: $𝑝 \backsim 𝑏𝑒𝑡𝑎(𝛼,𝛽)$ . We constructed an iteration algorithm to optimally estimate weight w and  $\alpha$  and  $\beta$  and used $\hat{\alpha}$  and  $\hat{𝛽}$  to estimate 𝑝 from p initial value and variance V of 𝑝. Note that different types of RNA-seq data have the initial different values of 𝑝. RNA-seq data are often generated from small samples with, for example, 3 ~ 6 replicates (or libraries). Small samples have big divergent variances and high gap probability(see Figure 1):

![image](https://user-images.githubusercontent.com/14003650/185705154-323cf284-fb98-49ac-9ee4-84c01cea82d2.png)

Figure 1. Statistical effect of small samples
The dataset was simulated from 10,000 experiments with two conditions randomly sampled from a negative binomial distribution NB(100, 50). a: variances distribution along mean. b: sorted variances along experiments. c: Probability of gap occurring between two datasets.

To handle these two statistical effects, we introduce a variable called gene-wise or isoform-wise inflation-shrinkage variable: 

![\Large \rho_{gi}=\sqrt{\varphi_{gi}\zeta{gi}}](https://latex.codecogs.com/svg.latex?\Large&space;\rho_{gi}=\sqrt{\varphi_{gi}\zeta{gi}}) 

where 


![\Large \varphi_{gi}=max\left[\frac{min(X_{gAi})}{max(X_{gBi})},\frac{min(X_{gBi})}{max(X_{gAi})}\right]](https://latex.codecogs.com/svg.latex?\Large&space;\varphi_{gi}=max\left[\frac{min(X_{gAi})}{max(X_{gBi})},\frac{min(X_{gBi})}{max(X_{gAi})}\right]) 

and


![\Large \zeta_{gi}=\ln\left(1+\frac{\bar{X_{gi}}\sigma_{gi}^2+1}{\bar{X_{gAi}}\sigma_{gAi}^2+\bar{X_{gBi}}\sigma_{gBi}^2+1}\right)](https://latex.codecogs.com/svg.latex?\Large&space;\zeta_{gi}=\ln\left(1+\frac{\bar{X_{gi}}\sigma_{gi}^2+1}{\bar{X_{gAi}}\sigma_{gAi}^2+\bar{X_{gBi}}\sigma_{gBi}^2+1}\right)) 
 

where 

![\Large \bar{X_{gi}}=\frac{1}{2}\left(\bar{X_{gAi}}+\bar{X_{gBi}}\right)](https://latex.codecogs.com/svg.latex?\Large&space;\bar{X_{gi}}=\frac{1}{2}\left(\bar{X_{gAi}}+\bar{X_{gBi}}\right)) 


and 


![\Large \bar{X_{gki}}=\frac{1}{m_k}\sum_{j=1}^{m_k}x_{gkij}](https://latex.codecogs.com/svg.latex?\Large&space;\bar{X_{gki}}=\frac{1}{m_k}\sum_{j=1}^{m_k}x_{gkij}) 


where k=1 for A and 2 for B. $\rho_{gi} > 1$ means that thera is a gap between data A and B for isoform i within gene g, $\rho_{gi}<1$  means that there is an overlap between data A and B for isoform i within gene g. 

Test for differential screens of the ith sgRNA targeting gene g using

$$t_{gi}^\alpha=\frac{\rho_g}{\omega_\alpha}t_{gi}$$

test for differential expression of the ith RNA isoform of gene g using 

$$t_{gi}^\alpha=\frac{\rho_{gi}}{\omega_\alpha}t_{gi}$$

test for differential expression or screen of gene g using

$$t_{g}^\alpha=\frac{\rho_{g}}{\omega_\alpha}t_{g}$$

where $\omega_\alpha$ is an averaged null $\rho$ under statistical significance $\alpha$, used as threshold for $\rho $. If $\rho > \omega_\alpha$, then $𝑡^\alpha > 𝑡$, meaning t-value is inflated. This case happens when there is significant difference between conditions. If $\rho < \omega_\alpha$, then $𝑡^\alpha < 𝑡$ , meaning t-value is shrunken. This case happens when there is no significant difference between two data. If $\rho =\omega_\alpha$, then $𝑡^\alpha = 𝑡$. This case occurs when sample size > 15. In current differential analysis of RNA-seq, NBBt-test has the lowest type I error rate and the best performance in small samples.

## Install NBBttest
Three ways to install NBBttest.
1. NBBttest can be installed from GitHub using devtools in R Console (Mac) or RGui(WIndows) or Rstudio:
```
library(devtools)

install_github("yuande/NBBttest")
```
2. Directly install NBBttest from R console(Mac) or RGui(window): 
Click package on R console or RGui and choose install package, choose CRAN mirrors and click OK and find NBBttest and click it.
3. Use R install.packages function to install NBBttest:
```
install.packages("NBBttest")
```
This way also needs to choose CRAM mirrors.

```
library(NBBttest)
```

## Data Preparation
When RNA-seq data have been produced from RNA sequecing experiments, user should first perform pipeline analysis of the RNA sequence read data and map RNA sequences to a reference genome. Currently many pipeline tools such as BWA, Bowtie2, tophat2, star and galaxy can be used to map and annotate RNA sequences on a reference genome. The pipeline analysis generates count matrix. The count matrix contains two parts: Annotation information and count data. Information may contain tagid (isoform_id or exon_id), geneid, gene name, chromosome, DNA strand, etc columns, depending on a pipeline tool that user used. Information columns are in the left side of the matrix. It has at least one column for geneid or tagid (isoformid). The count data contain two conditions each having several replicate libraries and must be in the right side. Here is an example:
```
data(jkttcell)

jkttcell[1:10,]
```
or using head to display the data jkttcell:
```
head(jkttcell)
```
Data jkttcell is a matrix count dataset generated from RNA-seq data due to differential polyadenylation in Jurkat T-cell betweem resting and stimulating statuses using BWA. Data jkttcell contains 7 columns for information of poly(A) sites in the left side and 10 columns for count data.

## Check data quality
Use QC to plot data of two replicates. The following code is an example to show plot of log2 data of replicates NS_A and NS_B in Jurkat T-cell RNA-seq data(jktcell):
```
data(jkttcell)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", log = "log", col = "blue", pch = 19)
```
Here "nci" is column number for data information, data begins with column 8 and ends with column 13. S1 is sample repicate1 and S2 is sample replicate2. Replicate1 is specified in column 8 and replicate2 is specified in column 9. "method" has two options: "plot" and "heatmap" and log also has two options: "none" and 
"log". This code outputs plot:
![image](https://user-images.githubusercontent.com/14003650/185797051-d296446b-2d46-426e-b921-076433a76d85.png)

Figure 2. QCplot of count data of two replicate libraries in no stimulated cells. Replicate libraries NS.A and NS.B have high correlation: All data dots follow on neighbor area of diagonal line.

```
QC(dat=jkttcell, nci=7, method = "heatmap", log = "log")}
```
This code outputs heatmap of correlation efficients between replicates:
![image](https://user-images.githubusercontent.com/14003650/185797319-3f6fd78d-0001-4811-bdc3-87215d3b7750.png)

Figure 3. Heatmap for correlation coefficients between replicate libraries. NS.A, NS.B, and NS.C are replicate libraries in no stimulated cells and 48h.A, 48h.B and 48h.C are replicate libraries 48h stimulated cells.

## Perform differential analysis
NBBttest has mbetattest function. This function performs differential expression of genes(RNA), differential splicing of RNA-isoforms, differential adenylation of RNA-isoforms (poly(A)), differential CRISPR screening.
```
data(jkttcell) 
res<-mbetattest(X=jkttcell, nci=7, na=3, nb=3, alpha=0.05, norm="yes", side="both", level="isoform", padjust_methods="fdr")
}
```
Here X is an object of NBBttest, that is, cout data, "nci" is number of columns for data information such as gene_id, gene name, exon_id or isoform_id or tag_id, chrosomome,annotation, strand etc. "na" and "nb" are numbers of replicates or libraries in conditions A and B, respectively. "alpha" is significance level, the default is 0.05. "norm" indicates whether data are normalized or not. "side" indicates one-side t-test or two-side t-test. "side" has three options: "both", "up" and "down". If side="up", then p-value is given with t-test in the left tail. If side="down", p-value is given with t-test in right tail. If side ="both", p-value is given with t-test in two sides. "level" has 6 options: "isoform", "sgRNA", "RNA", "splicing.gene", "polyA.gene", and "CRISPR.gene". If user's data come from poly(A) RNA-seq or splicing RNA-seq in which gene has multiple RNA isoform and wants differential splicing or differential adenylation, then user can set level="isoform". If user's data are CRISPR RNA-seq data or small hairpin RNA-seq (shRNA-seq), then user can set level="sgRNA". "splicing.gene", "polyA.gene", and "CRISPR.gene" are options for differential expression at gene level for data from splicing RNA-seq, poly(A) RNA-seq and CRISPR RNA-seq,respectively. If user's data are RNA data where each gene has one RNA isoform, then user can set level="RNA". 

```
exonskip<-read.table("merge_graphs_exon_skip_C3_count.txt", header=T)
dim(exonskip)
res<-mbetattest(X=exonskip,nci=10,na=11,nb=10,alpha=0.05, norm="yes",side="both", level="splicing.gene")
```
```
exonskp<-read.table("merge_graphs_exon_skip_C3_count.txt", header=T)
dim(exonskp)
exonskp[1:10,1:10]
res<-mbetattest(X=exonskp,nci=10,na=15,nb=32,alpha=0.05, norm="yes", side="both", level="isoform")
```

```
setwd("/Volumes/WD2/Band_data/Band_data_results/spladder_results/DDX39_GRCh38_93_gff3_result/")
DDX39.mutek<-read.table("merge_graphs_mutex_exons_C3_count.txt",header=T)
DDX39.mutek.isof<-mbetattest(X=DDX39.mutek, nci=12, na=3, nb=3, alpha=0.05, norm="yes",side="both", level="isoform")

```
## Annotation
Alternative splicing is detected in any element of 3'UTR, 5'UTR, exons and introns within a gene 
using RNA-seq data where RNA reads are mapped to a reference genome. As an example for annotation, 
the RNA-seq reads derived from human samples can be mapped onto human genome reference (GRCh38) 
using different methods, for example, HTSeq, spladder, rMAT, cufflinks, etc. These methods can detect 
alternative splicing sites within genes. However, none of these methods does gene annotation for users. 
Our NBBttest offers a R function for annotating genes with exons or isoforms.
```
data(DDX39_100)
data(gtfa)
DDX39_30<-annotat(infile=DDX39_100,mfile=gtfa,type="gene")
```
## Heatmap for differential expression or differential splicing or differential CRISPR screening

NBBttest has two heatmap functions: myheatmap and myheatmap2. The myheatmap uses z-score to normalize data and uses heatmap.2 to visualize differential expression of genes or isoforms detected. The myheatmap2 uses n-score to normalize the data and uses heatmap.2 to visualize differential expression of genes or isoforms detected by NBBttest.
```
data(result)
colclass=c("1","1","1","2","2","2")
oldpar <- par(no.readonly =TRUE)
par(mar=c(7.5,5.5,3.5,1.2))
par(oma=c(3,1,1,3))
myheatmap (dat=result, IDcol=1, nci=7, r=6, r1=3,r2=3, colrs="terrain.colors", rowBarColor=NULL, 
colBarColor=colclass, labrow="no", labcol="yes", rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
```
"dat" is result outputted by mbetattest and contain data information columns(in this example, data information column =7), data columns(r, r1 and r2)  and t-test result column. The t-test result columns must have "selection" or "select" column that lists "1" for DE genes (or DE isoforms) or "0" for no differentially expressed genes or isoforms. This code outputs the following heatmap plot:
![image](https://user-images.githubusercontent.com/14003650/185808658-ccad6681-be28-4f70-ae39-7c5b44fc4693.png)

Figure 4. Heatmap made with z-score for differential expressions of pol(A) RNA isoforms between stimulated and no stimulated cells.

The myheatmap2 uses selection to choose genes or isoforms in the data and then uses n-scale to normalize the selected data. Different from z-score, n-score does not follow standard normal distribution with mean = 0 and 
variance =1 for all rows but it has the same largest count in all rows and shows multiple colors for numeric 
difference between two conditions. The myheatmap2  has multiple options to select map color, distance, cluster 
and x-lab and y-lab angles. It can be able to display multiple datasets in two ways: if multiple 
datasets have the same row names or features, the these datasets are put onto the different columns separated 
with empty column named with dataset names. If multiple datasets have the same column names of the datasets, 
then put them on different rows saparated with empty rows named with dataset names or whatever names user 
specifies.  

```
data(result)
colclass=c("1","1","1","2","2","2")
oldpar <- par(no.readonly =TRUE)
par(mar=c(7.5,5.5,3.5,1.2))
par(oma=c(3,1,1,3))
myheatmap2(dat=result, IDcol=1, nci=7, r=6, colrs="bluered", rwcex=1.8, clcex=1.8, x=10, tree="both", method="euclidean", ky=1.5, rowBarColor=NULL,  colBarColor=colclass, labrow="no", labcol="yes", adjrow=c(0.2, 0.0 ), adjcol = c(1, 1) , rwangle=0, clangle=30, maptitle="My heatmap2",keyvalue="count")
```
![image](https://user-images.githubusercontent.com/14003650/185811378-9af27083-16ef-4a47-937a-4af62d309caa.png)

Figure 5. Heatmap made with n-score for differential expressions of pol(A) RNA isoforms between stimulated and no stimulated cells.

NBBttest also provides pathwayHeatmap. The pathwayHeatmap function is used to show differential expressions 
of pathways or functions between conditions. These pathways or functions were detected by function 
annotation or gene ontology methods such as David function analysis tools or Ingenuity pathway analysis or Metascape. The p-values for DE genes are used as weights for expression value across genes in a pathway or a function. The pathwayHeatmap function needs two datasets: count dataset and pathway dataset. The count dataset contain gene column, count data columns and p-value column. The pathway dataset contain pathway column and gene columns.

```
data(upGAm)
data(pathwy.A.up)
pathwayup<-pathwy.A.up
colclass=c("1","1","1","1","2","2","2","2","2","2")
par(mar=c(7.5,5.5,3.5,3))
par(oma=c(3,1,1,10))
pathwayHeatmap(dat=upGAm,pathway=pathwayup,nci=1,r1=4,r2=6,colclass=colclass,rowclass=NULL,colrs="greenred",maptitle="pathway up-expression in Group A")
```

![image](https://user-images.githubusercontent.com/14003650/185812122-007f8a93-3d13-4400-b0ee-cd9e4eb54364.png)

Figure 6. Heatmap for pathways.

## NBBplot
NBBttest provides NBBplot to display map and expression counts of exons within a specified gene across replicates in two conditions.
```
data(exondata)

NBBplot(res=exondata, gene="H2-DMb1", nci=9, na=3, nb=3, C1="WT", C2="KO")
```
![image](https://user-images.githubusercontent.com/14003650/185812883-404f75c2-c07c-4ac0-8496-0af0e1f9ecbd.png)

Figure 7. NBBplot for differential expression of exons in gene H2-DMb1. The top part is expressions of exons across replicates in two conditions KO and WT. The bottum part is phyical map of exons in gene H2-DMb1. The red boxes are DE exons.

## Reference
Tan YD and Guda C NBBt-test: a versatile method for differential analysis of multiple types of RNA-seq data. Scientific Report, 12833 (2022). (https://www.nature.com/articles/s41598-022-15762-x)




