# fvm-gmsfem

Implementation of the Generalized Multiscale Finite Element Method (GMsFEM) for solution problems in heterogeneous media with finite volume approximation on the fine grid.

Papers about GMsFEM:

* Efendiev Y, Galvis J, Hou TY. Generalized multiscale finite element methods (GMsFEM) [PDF](https://arxiv.org/abs/1301.2866)
* Chung ET, Efendiev Y, Li G, Vasilyeva M. Generalized multiscale finite element methods for problems in perforated heterogeneous domains [PDF](https://arxiv.org/abs/1501.03536)
* Akkutlu IY, Efendiev Y, Vasilyeva M. Multiscale model reduction for shale gas transport in fractured media [PDF](https://arxiv.org/abs/1507.00113)
* Vasilyeva M, Chung ET, Efendiev Y, Tyrylgin A. A three-level multi-continua upscaling method for flow problems in fractured porous media [PDF](https://arxiv.org/abs/1810.01581)

Implementation of the method contains several parts:

* **fine grid system generation** (./systemT.py) - create mass and stiffness matrices and right - hand side vector
* **local domains (coarse grid) generation** (./local-domain/) - create files with coarse cells coordinates and cell indices in local domains
* **multiscale basis function calculation** (./gmsfem-basis/) - solve local spectral problems to generate and save multiscale basis functions
* **projection matrix generation** (./ms-rgen/) - load multiscale basis functions and create projection matrix R
* **fine scale and multiscale solver** (./solver/) - solve fine grid system or/and multiscale solver

Implementation based on the [FEniCS](https://fenicsproject.org) (geometry objects, functions for saving and visualization) and [PETSc](https://www.mcs.anl.gov/petsc/).

## How to use

Fine grid simulations:
1. run fenics container
  > ```docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable```
2. create folders ./data/out/, ./data/modelF/
3. generate fine grid system for a given heterogeneous permeability field in ./data/perm/
  > ```python systemT.py```
4. run fine grid solver in ./solver/
  > ```./solver F 50 1.0e-4 80 80 ../data/modelF/ ../data/out/ 1 ./ err.txt```

Multiscale simulations:
1. create folders ./data/modelMs/omega10/, ./data/modelMs/eigen/, ./data/modelMs/dof/
2. local domains generations (coarse grid) in ./local-domain/
  > ``` ./omegas 2 ../data/omega10/ 10 0.1 0.0 10 0.1 0.0```
3. multiscale basis generation in ./gmsfem-basis/
  > ``` ./run```
4. generate R in ./ms-rgen/
  > ``` ./rgen 1 1600 121 ../data/modelMs/ 8```
5. solve multiscale in ./solver/
  > ``` ./solver C 50 1.0e-4 80 80 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt```
