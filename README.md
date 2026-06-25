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

## Distance modulus

We clarify a common source of confusion in distance-modulus calculations: the distinction between bolometric flux and flux density.

The usual distance modulus,

$$m_{\rm bol}=M_{\rm bol}+5\log_{10}\left(\frac{D_L}{10\,{\rm pc}}\right),$$

is defined for bolometric quantities. Here $D_L$ is the luminosity distance, which already includes the two cosmological dimming factors associated with photon-energy redshift and time dilation.

Most observed magnitudes, however, are based on a flux density, usually $f_\nu$, rather than on a bolometric flux. For a monochromatic AB-like measurement matched to the corresponding rest-frame frequency, the relation becomes

$$m_{\rm UV}=M_{\rm UV}+5\log_{10}\left(\frac{D_L}{10\,{\rm pc}}\right)-2.5\log_{10}(1+z).$$



The extra term is a bandwidth effect. A frequency interval observed with width $d\nu_{\rm obs}$ corresponds to a wider interval in the emitted frame,

$$d\nu_{\rm em}=(1+z)\,d\nu_{\rm obs}.$$

The flux density $f_\nu$ ($\mathrm{erg\,s^{-1}\,{Hz}^{-1}}$) is defined per unit observed frequency. The transformation between emitted and observed frequency intervals therefore has a factor of $1+z$:

$$f_{\nu_{\rm obs}}=\frac{(1+z)\,L_{\nu_{\rm em}}\left[(1+z)\nu_{\rm obs}\right]}{4\pi D_L^2}.$$

When written in magnitudes, this factor gives the term

$$-2.5\log_{10}(1+z).$$

In contrast, the bolometric flux ($\mathrm{erg\,s^{-1}}$)  is,

$$F_{\rm bol}=\int f_\nu\,d\nu,$$

where the bandwidth transformation is integrated over. 



This $-2.5\log_{10}(1+z)$ factor is necessary, for example, when converting observed magnitudes of high-redshift galaxies into absolute rest-frame UV magnitudes $M_{\rm UV}$.
