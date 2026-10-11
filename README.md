# Supplementary material

[![arXiv](https://img.shields.io/badge/arXiv-2607.26248-b31b1b.svg)](https://arxiv.org/abs/2607.26248)

This repository contains the MATLAB source code to reproduce the numerical simulations presented in our paper:

> **"Strong imposition of Dirichlet boundary velocities in structure-preserving discretizations of elastodynamics"**

**Authors:** Cristobal Ponce¹, Hector Ramirez¹, Yongxin Wu², Ning Liu², and Yann Le Gorrec²  
¹ *Universidad Técnica Federico Santa María, Valparaíso, Chile.*  
² *Université Marie et Louis Pasteur, SUPMICROTECH, CNRS, Institut FEMTO-ST, Besançon, France.*

## Repository structure

The source code is located in the `Matlab Code` directory. It provides a general implementation framework for linear systems, along with specific subfolders corresponding to the numerical examples discussed in the manuscript:

- `General code for 1D and 2D linear systems/`
- `1. Linear Timoshenko beam/`
- `2. Geometrically exact spatial string/` *(Note: Includes one additional example not shown in the paper)*
- `3. Nonlinear 2D Elasticity - Neo-Hookean material/`

## System requirements

- **Software:** MATLAB R2021b or newer.
- **Toolboxes:** Partial Differential Equation (PDE) Toolbox.
- **Hardware:** A standard personal computer is sufficient to run the codes. Detailed simulation times for each example are reported in the paper.

## Citation policy

If you use this code or these models in your own research, we kindly request that you cite the appropriate sources:

```bibtex

1. For the underlying linear port-Hamiltonian models:
@article{ponce2024systematic,
  title={A systematic methodology for port-Hamiltonian modeling of multidimensional flexible linear mechanical systems},
  author={Ponce, Cristobal and Wu, Yongxin and Le Gorrec, Yann and Ramirez, Hector},
  journal={Applied Mathematical Modelling},
  volume={134},
  pages={434--451},
  year={2024},
  publisher={Elsevier}
}


2. For the underlying nonlinear/hyperelastic port-Hamiltonian models:
@article{ponce2025port,
  title={A port-Hamiltonian framework for the modeling and FEM discretization of hyperelastic systems},
  author={Ponce, Cristobal and Wu, Yongxin and Le Gorrec, Yann and Ramirez, Hector},
  journal={Applied Mathematical Modelling},
  pages={116403},
  year={2025},
  publisher={Elsevier}
}

3. For the numerical implementation (Strong imposition of Dirichlet boundary velocities):
@article{ponce2026strong,
  title={Strong imposition of Dirichlet boundary velocities in structure-preserving discretizations of elastodynamics},
  author={Ponce, Cristobal and Ramirez, Hector and Wu, Yongxin and Liu, Ning and Le Gorrec, Yann},
  journal={arXiv preprint arXiv:2607.26248},
  year={2026}
}

