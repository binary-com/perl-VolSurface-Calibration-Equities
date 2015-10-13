### VolSurface::Calibration::SABR

This repository is the Binary.com's equities volatility calibration - a variant of standard SABR model. Basically a Volatility Surface Calibration algorithm tries to check a Vol-Surface to make sure it satisfies some basic requirements. If not, it will change the surface to make it a valid VolSurface. The standard SABR (stochastic alpha, beta, rho_ model is used to estimate implied volatility of an instrument in the derivatives market. 

```

σ_impl= α x log⁡(F_0⁄K)/D(ζ) x  {1+[(2γ_2-γ_(1 )^2+1⁄(F_mid^2 ))/24 x   ((σ_0 C(F_mid ))/α)^2+ (ργ_1)/4 x  (σ_0 C((F_mid ))/α+ (2-3ρ^2)/24]  x ε }

```

The term is basically expressing implied volatility as some sort of moneyness function (the alpha x log(F/K)D(...) part). This term is adjusted by some factor in the square brackets and then added back. In our variant of the standard SABR approach, we modify the terms in the square brackets. This calibration approach is based upon modeling the term structure of ATM volatiity and Skew using exponential functions, as it is widely observed that ATM vols term structure or skew term structure is mostly convex.

#### Documentation

Further details of the calibration model is available in MS Word and pdf formats at :

https://github.com/mm-binary/perl-VolSurface-Calibration-SABR/blob/mm/initial_movement/documentation/Binary's_equities_volatility_calibration.docx

https://github.com/mm-binary/perl-VolSurface-Calibration-SABR/blob/mm/initial_movement/documentation/Binary's_equities_volatility_calibration.pdf


To read more about how this algorithm works please refer to "SABR\_RMG\_Implementation.pdf" file in the "/doc" directory.
