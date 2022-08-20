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

$\rho_{gi}\ge 1$ means a gap between data A and B, $\rho_{gi} le 1$  means overlap between data A and B. 

Test for differential screens of the ith sgRNA targeting gene g using
$$t_{gi}^\alpha=\frac{\rho_g}{\omega_\alpha}t_{gi}$$,

test for differential expression of the ith RNA isoform of gene g using 
$$t_{gi}^\alpha=\frac{\rho_{gi}}{\omega_\alpha}t_{gi}$$,



