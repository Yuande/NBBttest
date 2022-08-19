# NBBttest
## Introduction
NBBttest is R package for implementing negative binomial beta t-test (or called NBBt-test). NBBt-test is a statistical method for differential analysis of multiple types of RNA-seq data. NBBt-test is based on negative binomial distribution with parameters p and r: 

 $$ x \backsim NB(p,r)$$
 
where r is number of RNA sequences failed to be sequenced in RNA population (or RNA library) and p is proportion of the sequenced RNA species in RNA library and follows beta distibution with parameters $\alpha$ and $\beta$:

$$ 𝑝 \backsim 𝑏𝑒𝑡𝑎(𝛼,𝛽)$$

We constructed an iteration algorithm to optimally estimate weight w and $\alpha$ and $\beta$ and used $\hat{\alpha}$ and $\hat{𝛽}$ to estimate 𝑝 and variance V of 𝑝. However, to estimate 𝑝 by iteration, we need an initial value of 𝑝 for each gene or each isoform



![image](https://user-images.githubusercontent.com/14003650/185698478-a8ad2f85-b673-49aa-a5d0-cea217879fa6.png)

