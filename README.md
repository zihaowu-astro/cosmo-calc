# Cosmology Calculator

A browser-based cosmology calculator for flat $\Lambda$CDM, in the spirit of
[Ned Wright's Cosmology Calculator](https://astro.ucla.edu/~wright/CosmoCalc.html). Line-of-sight integrals are evaluated with adaptive Simpson quadrature, and radiation is included as $\Omega_r = 4.165\times10^{-5}\,/\,h^2$. All quantities are validated against Astropy. Conventions follow
[Hogg (1999)](https://arxiv.org/abs/astro-ph/9905116).

## Parameters

| Set | $H_0$ | $\Omega_m$ | $\Omega_\Lambda$ |
|-----|-----|-----|-----|
| Planck 2018 | 67.66 | 0.3111 | 0.6889 |
| WMAP9 | 69.32 | 0.2865 | 0.7135 |

$H_0$ in $\mathrm{km\,s^{-1}\,Mpc^{-1}}$; both flat ($\Omega_k = 0$).

- Planck 2018 (TT,TE,EE+lowE+lensing+BAO): Planck Collaboration VI (2020),
  *A&A* **641**, A6 ([arXiv:1807.06209](https://arxiv.org/abs/1807.06209)), Table 2.
- WMAP9 (WMAP+eCMB+BAO+H₀): Hinshaw et al. (2013), *ApJS* **208**, 19
  ([arXiv:1212.5226](https://arxiv.org/abs/1212.5226)), Table 4.

