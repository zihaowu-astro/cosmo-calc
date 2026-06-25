# Cosmology Calculator

A small, dependency-free cosmology calculator in the spirit of
[Ned Wright's Cosmology Calculator](https://astro.ucla.edu/~wright/CosmoCalc.html).

You **cannot** edit the cosmological parameters — instead you pick one of two
standard, fixed parameter sets:

| Set | H₀ | Ω_M | Ω_Λ | Source |
|-----|----|-----|-----|--------|
| **Planck 2018** | 67.36 | 0.3153 | 0.6847 | Planck 2018 VI, TT,TE,EE+lowE+lensing+BAO |
| **WMAP9** | 69.32 | 0.2865 | 0.7135 | Hinshaw et al. 2013 (WMAP9 + eCMB + BAO + H₀) |

Both are flat ΛCDM. Radiation density (photons + 3 neutrino species) is added
following Wright's convention, Ω_R = 4.165×10⁻⁵ / h².

Enter a redshift *z* and the calculator returns:

- Age of the universe now and at *z*
- Light travel time (lookback time)
- Comoving radial distance
- Angular diameter distance and angular size scale (kpc/″)
- Luminosity distance and distance modulus
- Comoving volume enclosed within *z*

Everything runs client-side in plain JavaScript — distances are integrated with
an adaptive Simpson's rule.

## Run locally

It's a static site, so just open `index.html`, or serve the folder:

```bash
python3 -m http.server 8000
# then visit http://localhost:8000
```

## Deploy on GitHub Pages

1. Create a new GitHub repository and push these files:

   ```bash
   git init
   git add .
   git commit -m "Cosmology calculator"
   git branch -M main
   git remote add origin https://github.com/<you>/<repo>.git
   git push -u origin main
   ```

2. In the repo, go to **Settings → Pages** and set the source to
   **GitHub Actions**. The included workflow (`.github/workflows/pages.yml`)
   publishes the site automatically on every push to `main`.

   *(Alternatively, set Pages source to "Deploy from a branch → main / root" —
   no workflow needed since this is a static site.)*

Your calculator will be live at `https://<you>.github.io/<repo>/`.

## Files

- `index.html` — page structure
- `style.css` — styling
- `cosmo.js` — cosmology math + UI logic
