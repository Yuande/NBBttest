# NBBttest
## Introduction
NBBttest is R package for implementing negative binomial beta t-test (or called NBBt-test). NBBt-test is a statistical method for differential analysis of multiple types of RNA-seq data. NBBt-test is based on negative binomial distribution with parameters p and r: 

 $$ x \backsim NB(p,r)$$
 
where r is number of RNA sequences failed to be sequenced in RNA population (or RNA library) and p is proportion of the sequenced RNA species in RNA library and follows beta distibution with parameters $\alpha$ and $\beta$:

$$ 𝑝 \backsim 𝑏𝑒𝑡𝑎(𝛼,𝛽)$$

We constructed an iteration algorithm to optimally estimate weight w and $\alpha$ and $\beta$ and used $\hat{\alpha}$ and $\hat{𝛽}$ to estimate 𝑝 from p initial value and variance V of 𝑝. Note that different type of RNA-seq data, the initial value of 𝑝 is different. RNA-seq data are often generated from small samples with, for example, 3 ~ 6 replicates (or libraries). Small samples have bigger divergent variances and high gap probability(see Figure 1):

![image](https://user-images.githubusercontent.com/14003650/185705154-323cf284-fb98-49ac-9ee4-84c01cea82d2.png)
Figure 1. Statistical effect of small samples
The dataset was simulated from 10,000 experiments with two conditions randomly sampled from a negative binomial distribution NB(100, 50). a: variances distribution along mean. b: sorted variances along experiments. c: Probability of gap occurring between two datasets.

To handle these two statistical effects, we introduce a variable called gene-wise or isoform-wise inflation and shrinkage variable:

$$\rho_{gi}=\sqrt{\varphi_{gi}\zeta_{gi}}$$

where 

$$\varphi_{gi}=max\left[ \frac{min(X_{gAi})}{max(X_{gBi})},\frac{min(X_{gBi})}{max(X_{gAi})}\right]$$

and

$$\zeta_{gi}=\ln\left(1+\frac{\bar{X_{gi}}\sigma_{gi}^2 +1}{\bar{X_{gAi}}\sigma_{gAi}^2 +\bar{X_{gBi}}\sigma_{gBi}^2+ 1}\right)$$

where
$$\bar{X_{gi}}=\frac{1}{2}\left(\bar{X_{gAi}}+\bar{X_{gBi}}\right)$$
and 
$$\bar{X_{gki}}=\frac{1}{m_k}\sum_{j=1}^{m_k} x_{gkij}$$ where k=1 for A and 2 for B.
$\rho_{gi} > 1$ means that thera is a gap between data A and B for isoform i within gene g, $\rho_{gi}<1$  means that there is an overlap between data A and B for isoform i within gene g. 

Test for differential screens of the ith sgRNA targeting gene g using
$$t_{gi}^\alpha=\frac{\rho_g}{\omega_\alpha}t_{gi}$$
test for differential expression of the ith RNA isoform of gene g using 
$$t_{gi}^\alpha=\frac{\rho_{gi}}{\omega_\alpha}t_{gi}$$
test for differential expression or screen of gene g using
$$t_{g}^\alpha=\frac{\rho_{g}}{\omega_\alpha}t_{g}$$
where $\omega_\alpha$ is an averaged null $\rho$ under statistical significance $\alpha$, used as threshold for $\rho $. If $\rho > \omega_\alpha$, then $𝑡^\alpha > 𝑡$, meaning t-value is inflated. This case happens when there is significant difference between conditions. If $\rho < \omega_\alpha$, then $𝑡^\alpha < 𝑡$ , meaning t-value is shrunken. This case happens when there is no significant difference between two data. If $\rho =\omega_\alpha $, then $𝑡^\alpha = 𝑡$. This case occurs when sample size > 15.

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
jkttcell is matrix count data generated from RNA-seq data due to differential polyadenylation in Jurkat T-cell betweem resting and stimulating statuses using BWA. Data jkttcell contains 7 columns for information of poly(A) sites in the left side and 10 columns for count data.

## Check data quality
Use QC to plot data of two replicates. The following code is an example to show plot of log2 data of replicates NS_A and NS_B in Jurkat T-cell RNA-seq data(jktcell):
```
data(jkttcell)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", 
log = "log", col = "blue", pch = 19)
```
Here nci is column number for data information, data begins with column 8 and ends with column 13. S1 is sample repicate1 and S2 is sample replicate2. Replicate1 is specified in column 8 and replicate2 is specified in column9. method has two options: "plot" and "heatmap" and log also has two options: "none" and 
"log". This code outputs plot:
![image](https://user-images.githubusercontent.com/14003650/185797051-d296446b-2d46-426e-b921-076433a76d85.png)
```
QC(dat=jkttcell, nci=7, method = "heatmap", log = "log")}
```
This code outputs heatmap of correlation efficients between replicates:
![image](https://user-images.githubusercontent.com/14003650/185797319-3f6fd78d-0001-4811-bdc3-87215d3b7750.png)



