# Double-Pendulum-Simulation-MATLAB
Numerical simulation and chaos dynamics of a nonlinear double pendulum in MATLAB, evaluating Explicit Euler, RK4, and ODE45 solvers with mechanical energy conservation analysis.
# Numerical Simulation and Chaos Dynamics of a Nonlinear Double Pendulum

This repository contains MATLAB implementations for the numerical simulation, stability analysis, and chaos dynamics of a double compound pendulum system. The project models nonlinear coupled equations of motion and systematically compares the performance, stability, and energy conservation of different Ordinary Differential Equation (ODE) integration schemes.

---

## 📌 Project Overview

The double pendulum is a classical benchmark in nonlinear mechanics exhibiting extreme sensitivity to initial conditions (deterministic chaos) at elevated energy states. 

This project explores:
- Derivation and state-space representation of the coupled nonlinear equations of motion.
- Small-angle oscillations (quasi-linear periodic motion) vs. large-angle regimes (chaotic dynamics) using MATLAB's adaptive `ode45` solver.
- Sensitivity to initial conditions by perturbing starting angles by $\Delta \theta = 0.001 \text{ rad}$ to observe trajectory divergence (Lyapunov behavior).
- Implementation of fixed-step explicit numerical integrators: **Explicit Euler** and **4th-Order Runge-Kutta (RK4)** across step sizes $\Delta t \in \{0.01, 0.005, 0.001\} \text{ s}$.
- Total mechanical energy monitoring ($E = T + V$) to evaluate numerical damping and artificial energy divergence.

---

## 📐 Mathematical Formulation

The dynamic equations of motion derived via Lagrangian mechanics are formulated as a 4-dimensional first-order state-space system $\mathbf{\dot{x}} = f(t, \mathbf{x})$:
$$\mathbf{x} = \begin{bmatrix} \theta_1 & \omega_1 & \theta_2 & \omega_2 \end{bmatrix}^T$$

Where $\theta_1, \theta_2$ denote angular displacements and $\omega_1, \omega_2$ represent angular velocities of the primary and secondary links. The equations are numerically integrated in `pendulum_equations.m`.

---

## 📂 Repository Structure

```text
├── Scripts/
│   ├── pendulum_equations.m      # State-space ODE function of the double pendulum[cite: 2]
│   ├── part2.m                   # Small vs. large angle dynamics using ode45[cite: 2]
│   ├── part3.m                   # Chaos analysis & sensitivity to initial conditions[cite: 2]
│   ├── part4.m                   # Explicit Euler vs. RK4 comparative solver[cite: 2]
│   └── part5.m                   # Mechanical energy tracking (T, V, E_total)[cite: 2]
├── Plots/
│   ├── small_vs_large_angle.png  # Angular trajectory and phase space plots[cite: 2]
│   ├── chaos_divergence.png      # Trajectory separation over time[cite: 2]
│   └── energy_conservation.png   # Euler drift vs. RK4 energy conservation[cite: 2]
├── Documentation/
│   └── Final_Project_Report.pdf  # Comprehensive mathematical derivation and results[cite: 2]
└── README.md
