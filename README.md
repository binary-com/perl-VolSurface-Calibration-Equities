### VolSurface::Calibration::SABR

This repository is Binary.com's equities volatility calibration - a variant of standard SABR model. This assumes that the input volatility surface is in moneyness terms, with a matrix of moneyness and tenor defined. 

Basically a Volatility Surface Calibration algorithm tries to check a Vol-Surface to make sure it satisfies some basic requirements. If not, it will change the surface to make it a valid VolSurface - satisfying arbitrage free prices across tenors and strikes. The standard SABR (stochastic alpha, beta, rho) model is used to estimate implied volatility of an instrument in the derivatives market. It is one of the most popular models in the industry, represented by :

```

σ_impl= α x log⁡(F_0⁄K)/D(ζ) x  {1+[(2γ_2-γ_(1 )^2+1⁄(F_mid^2 ))/24 x   ((σ_0 C(F_mid ))/α)^2+ (ργ_1)/4 x  (σ_0 C((F_mid ))/α+ (2-3ρ^2)/24]  x ε }

```
The above equation basically expresses implied volatility as some sort of moneyness function (the alpha x log(F/K)D(...) part). This term is adjusted by a factor in the square brackets and then added back. In our variant of the standard SABR approach, we modify the terms in the square brackets. This calibration approach is based upon modeling the term structure of ATM volatiity and Skew using exponential functions, as it is widely observed that ATM vols term structure or skew term structure is mostly convex. We have observed that this variant results in a more consistent option prices.

For optimization, we use a form of the Downhill Simplex Method or Nelder-Mead (available as the R function optim). 

#### Documentation

Further details of the calibration model is available in MS Word and pdf formats at :

https://github.com/mm-binary/perl-VolSurface-Calibration-SABR/blob/mm/initial_movement/documentation/Binary's_equities_volatility_calibration.docx

https://github.com/mm-binary/perl-VolSurface-Calibration-SABR/blob/mm/initial_movement/documentation/Binary's_equities_volatility_calibration.pdf
