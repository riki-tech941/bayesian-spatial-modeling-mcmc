The project is split into two parts:

1. **Spatial Gaussian Process (Air quality):**
   - Fits a spatial GP on European PM10 data.
   - Uses a Gibbs sampler with an adaptive Metropolis-Hastings step on log(phi) to get a stable acceptance rate around 20-25%.
   - Handles covariance inversion using Cholesky factorizations and forward substitution to avoid numerical issues.
   - Includes regular grid prediction (kriging) and an analytical marginalization step to compare sampling performance.

2. **Trajectory Hidden Markov Model (Movement tracking):**
   - 2-state and 3-state HMM models fitted to 2D tracking data.
   - Inverse-Wishart priors for covariance matrices and Dirichlet priors for transition probabilities.
   - Latent paths sampled in block using Forward-Backward (Viterbi) instead of scalar Gibbs to reduce autocorrelation.

## Structure

- `src/`: R script with samplers, diagnostics, and plots.
- `report/`: PDF report with mathematical derivations, trace plots, and results.
- `data/`: input datasets (`data\\\_dogs2.csv`, `grid.csv`, `data\\\_for\\\_homework2.RData`).

## Author

Riccardo Siervo
[LinkedIn](https://linkedin.com/in/riccardo-siervo)
