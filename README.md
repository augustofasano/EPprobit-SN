# EPprobit-SN

This repository is associated with the conference proceedings [Fasano, Anceschi, Rebaudo and Franzolini (2023a)](https://giovannirebaudo.github.io/Publications/2023FasanoAnceschiFranzoliniRebaudo.pdf) *Efficient expectation propagation for high-dimensional probit models* and [Fasano, Anceschi, Rebaudo and Franzolini (2023b)](https://giovannirebaudo.github.io/Publications/2023CLADAGFasanoAnceschiFranzoliniRebaudo.pdf) *Efficient computation of predictive probabilities in probit models via Expectation Propagation*. The key contributions of the two papers are outlined below.

[Fasano, Anceschi, Rebaudo and Franzolini (2023a)](https://giovannirebaudo.github.io/Publications/2023FasanoAnceschiFranzoliniRebaudo.pdf):
> [...] we focus on the expectation propagation (EP) approximation of the posterior distribution in Bayesian probit regression under a multivariate Gaussian prior distribution. Adapting more general derivations in [Anceschi et al. (2023)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2023.2169150), we show how to leverage results on the extended multivariate skew-normal distribution to derive an efficient implementation of the EP routine having a per-iteration cost that scales linearly in the number of covariates.

[Fasano, Anceschi, Rebaudo and Franzolini (2023b)](https://giovannirebaudo.github.io/Publications/2023CLADAGFasanoAnceschiFranzoliniRebaudo.pdf):
> [...] we focus on the computation of posterior predictive probabilities in Bayesian probit models via Expectation Propagation (EP).
Leveraging more general results in recent literature, we show that such predictive probabilities admit a closed-form expression.

This repository provides codes to replicate the simulation studies reported in Section 4 of [Fasano, Anceschi, Rebaudo and Franzolini (2023a)](https://giovannirebaudo.github.io/Publications/2023FasanoAnceschiFranzoliniRebaudo.pdf) and in Section 4 of [Fasano, Anceschi, Rebaudo and Franzolini (2023b)](https://giovannirebaudo.github.io/Publications/2023CLADAGFasanoAnceschiFranzoliniRebaudo.pdf).

More precisely, we provide the `R` code to implement **Algorithms 1 and 2 presented in [Fasano, Anceschi, Rebaudo and Franzolini (2023a)](https://giovannirebaudo.github.io/Publications/2023FasanoAnceschiFranzoliniRebaudo.pdf)** to obtain **efficient EP approximations of the posterior moments** of the parameters in a probit model with spherical Gaussian prior distribution.
We recall that Algorithms 1 and 2 have per-iteration costs O(p<sup>2</sup>n) and O(pn<sup>2</sup>), respectively, so the former is used when p<n, the latter otherwise. After the algorithm has converged, one can then **efficiently compute the approximated EP posterior predictive probabilities** for held-out units exploiting the closed-form expressions presented in [Fasano, Anceschi, Rebaudo and Franzolini (2023b)](https://giovannirebaudo.github.io/Publications/2023CLADAGFasanoAnceschiFranzoliniRebaudo.pdf).

In addition, we also provide code to perform posterior inference with other three different methods, used for comparison purposes:

1. **i.i.d. sampling from the exact unified skew-normal distribution** [(Durante, 2019)](https://academic.oup.com/biomet/article-abstract/106/4/765/5554418)
2. **partially factorized mean-field variational Bayes (PFM-VB) approximation** [(Fasano, Durante and Zanella, 2022)](https://academic.oup.com/biomet/article-abstract/109/4/901/6581071)
3. **EP implemented via the `R` function `EPprobit` from the package `EPGLM`** [(Chopin and Ridgway, 2017)](https://projecteuclid.org/journals/statistical-science/volume-32/issue-1/Leave-Pima-Indians-Alone--Binary-Regression-as-a-Benchmark/10.1214/16-STS581.full)

Structure of the repository:

* the functions to implement the above methods can be found in the `R` source file [`functions.R`](https://github.com/augustofasano/EPprobit-SN/blob/main/functions.R)
* a tutorial with the code to reproduce the results in the papers is available at [`Illustration.md`](https://github.com/augustofasano/EPprobit-SN/blob/main/Illustration.md)
