# Cosmology Calculator

A browser-based cosmology calculator for flat $\Lambda$CDM, in the spirit of
[Ned Wright's Cosmology Calculator](https://astro.ucla.edu/~wright/CosmoCalc.html). Line-of-sight integrals are evaluated with adaptive Simpson quadrature. All quantities are validated against Astropy. Conventions follow
[Hogg (1999)](https://arxiv.org/abs/astro-ph/9905116).

## Parameters

| Set | $H_0$ | $\Omega_m$ | $\Omega_\Lambda$ |
|-----|-----|-----|-----|
| Planck 2018 | 67.66 | 0.3111 | 0.6889 |
| Planck 2015 | 67.74 | 0.3089 | 0.6911 |
| Planck 2013 | 67.77 | 0.3086 | 0.6914 |
| WMAP9 | 69.32 | 0.2865 | 0.7135 |
| WMAP7 | 70.40 | 0.2720 | 0.7280 |
| $H_0 = 70$, $\Omega_m = 0.3$ | 70.00 | 0.3000 | 0.7000 |

$H_0$ in $\mathrm{km\,s^{-1}\,Mpc^{-1}}$; all flat ($\Omega_k = 0$). The quoted
$\Omega_m$ is the published value, which counts the massive neutrino as matter;
neutrinos are then treated as massless in $\Omega_r$.

- Planck 2018 (TT,TE,EE+lowE+lensing+BAO): Planck Collaboration VI (2020),
  *A&A* **641**, A6 ([arXiv:1807.06209](https://arxiv.org/abs/1807.06209)), Table 2.
- Planck 2015 (TT,TE,EE+lowP+lensing+ext): Planck Collaboration XIII (2016),
  *A&A* **594**, A13 ([arXiv:1502.01589](https://arxiv.org/abs/1502.01589)), Table 4.
- Planck 2013 (Planck+WP+highL+BAO): Planck Collaboration XVI (2014),
  *A&A* **571**, A16 ([arXiv:1303.5076](https://arxiv.org/abs/1303.5076)), Table 5.
- WMAP9 (WMAP+eCMB+BAO+H₀): Hinshaw et al. (2013), *ApJS* **208**, 19
  ([arXiv:1212.5226](https://arxiv.org/abs/1212.5226)), Table 4.
- WMAP7 (WMAP+BAO+H₀ ML): Komatsu et al. (2011), *ApJS* **192**, 18
  ([arXiv:1001.4538](https://arxiv.org/abs/1001.4538)), Table 1.
- $H_0 = 70$, $\Omega_m = 0.3$, $\Omega_\Lambda = 0.7$.
