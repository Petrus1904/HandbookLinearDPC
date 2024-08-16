# Handbook Linear Data-driven Predictive Control - Examples
This repository contains MATLAB code examples of implementing different data-driven predictive controllers (DPC) on 2 different systems. 

## Requirements:
* **YALMIP** -- you could re-express the problem in a quadprog shape, but YALMIP is used for readability
* **MATLAB System Identification toolbox** -- For data acquisition and system identification of the MPC controllers.
* **OSQP** -- Optional. It is used as the solver, but this can be replaced by 'quadprog' in the solver settings.

## What it contains
**Plane Model Control**
* `MPC_planeModel` -- MPC with Luenberger observer
* `MPC_Kalman_planeModel` -- MPC with a Kalman filter
* `SPC_planeModel` -- Subspace Predictive Control
* `DeePC_planeModel` -- DeePC with l2 regularization on g
* `DeePC_PIreg_planeModel` -- DeePC with PI-regularization on g
* `gammaDDPC_planeModel` -- γ-DDPC, an LQ-decomposition DPC strategy
* `GDPC_planeModel` -- Generalized DPC that combines SPC and DeePC

**Four Tank System Control**
* `MPC_FourTank` -- MPC with Luenberger observer
* `MPC_Kalman_FourTank` -- MPC with a Kalman filter
* `SPC_FourTank` -- Subspace Predictive Control
* `DeePC_FourTank` -- DeePC with l2 regularization on g
* `gammaDDPC_FourTank` -- γ-DDPC, an LQ-decomposition DPC strategy
* `GDPC_FourTank` -- Generalized DPC that combines SPC and DeePC

All of the above mentioned scripts are standalone runs of the same control problem, but with a different control strategy. These scripts do depend on a few functions, including:
* `GetDataFourTankModel` -- Constructs system and experiment data (including the Hankel matrices) from the Four Tank model
* `GetDataPlaneModel` -- Constructs system and experiment data (including the Hankel matrices) from the Plane model
* `GetMarkovMatrix` -- Only used in the MPC scripts to build prediction matrices from state-space
* `HankeRausRegularization` -- A function that can be used to plot the Hanke-Raus regularization score. This score can be used to select most suitable hyperparameters in DeePC, γ-DDPC and GDPC. The implementation is commented in the corresponding run scripts.

Use `help [function_name]` to get detailed information on how to use the functions.

## How to Cite
**APA style:**

P.C.N. Verheijen, V. Breschi and M. Lazar. (2023). *"Handbook of linear data-driven predictive control: Theory, implementation and design."* Annual Reviews in Control, 56, https://doi.org/10.1016/j.arcontrol.2023.100914

**Bibtex:**
```
@article{HandbookDPC:Verheijen2023,
  title = {Handbook of linear data-driven predictive control: Theory, implementation and design},
  journal = {Annual Reviews in Control},
  volume = {56},
  pages = {100914},
  year = {2023},
  issn = {1367-5788},
  doi = {https://doi.org/10.1016/j.arcontrol.2023.100914},
  author = {P.C.N. Verheijen and V. Breschi and M. Lazar},
}
```


