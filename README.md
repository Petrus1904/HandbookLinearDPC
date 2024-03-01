# Handbook Linear Data-driven Predictive Control - Examples
This repository contains MATLAB code examples of implementing different data-driven predictive controllers (DPC) on 2 different systems. 

## Requirements:
* YALMIP -- you could re-express the problem in a quadprog shape, but YALMIP is used for readability
* MATLAB SysID toolbox -- For data acquisition and system identification of the MPC controllers.
* OSQP -- Optional. It is used as the solver, but this can be replaced by 'quadprog' in the solver settings.

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
