# NBBttest
## Introduction
NBBttest is R package for implementing negative binomial beta t-test (or called NBBt-test). NBBt-test is a statistical method for differential analysis of multiple types of RNA-seq data. NBBt-test is based on negative binomial distribution with parameters p and r: 

 $$ x \backsim NB(p,r)$$
 
where r is number of RNA sequences failed to be sequenced in RNA population (or RNA library) and p is proportion of the sequenced RNA species in RNA library and follows beta distibution with parameters $\alpha$ and $\beta$:

$$ 𝑝 \backsim 𝑏𝑒𝑡𝑎(𝛼,𝛽)$$

We constructed an iteration algorithm to optimally estimate weight w and $\alpha$ and $\beta$ and used $\hat{\alpha}$ and $\hat{𝛽}$ to estimate 𝑝 from p initial value and variance V of 𝑝. Note that different type of RNA-seq data, the initial value of 𝑝 is different. RNA-seq data are often generated from small samples with, for example, 3 ~ 6 replicates (or libraries). Small samples have bigger divergent variances and high gap probability(see Figure 1):
![image](https://user-images.githubusercontent.com/14003650/185704654-2f8a34fd-515e-4da2-8ce6-80b8c1ea876f.png)




![image](https://user-images.githubusercontent.com/14003650/185698478-a8ad2f85-b673-49aa-a5d0-cea217879fa6.png)

