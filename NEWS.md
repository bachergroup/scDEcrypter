# scDEcrypter 0.8.7

- Corrected the penalized mean update to weight each component by its
  posterior count divided by its variance.
- Replaced the accelerated mean substep with proximal-gradient steps whose
  step size guarantees descent of the penalized M-step objective.
- Aligned the observed-data log-likelihood for partly observed labels with
  the joint likelihood used by the E-step and M-step.
- Retained the converged EM iterate and exposed the penalized observed-data
  objective after every iteration as `objective_generation`.
- Kept the initialization parameters and returned weights synchronized.
