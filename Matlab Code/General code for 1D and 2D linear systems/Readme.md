# General code for 1D and 2D linear systems

This folder contains the MATLAB implementation for simulating 1D and 2D linear port-Hamiltonian elastodynamics. 

Specifically, it provides the numerical schemes for the strong imposition of Dirichlet boundary velocities discussed in our preprint:

*   **Standard FEM:** Based on Hamilton's principle (**Corollary 1**).
*   **Mixed-FEM:** Based on the Generalized Hamilton's principle / Hu-Washizu formulation (**Corollary 2**).

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

2. For the numerical implementation (Strong imposition of Dirichlet boundary velocities):
@article{ponce2026strong,
  title={Strong imposition of Dirichlet boundary velocities in structure-preserving discretizations of elastodynamics},
  author={Ponce, Cristobal and Ramirez, Hector and Wu, Yongxin and Liu, Ning and Le Gorrec, Yann},
  journal={arXiv preprint arXiv:2607.26248},
  year={2026}
}
