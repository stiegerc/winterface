
# Winterface

Winterface or 'Wannier Interface' is a code able to generate device scale Hamiltonian matrices starting from output of a tool called Wannier90 (www.wannier.org). Any Quantum Chemistry package (e.g. Quantum Espresso, VASP, ...) that is interfaced to Wannier90 can thus be interfaced to a Quantum Transport code. Since a Density Functional Theory (DFT) calculation is typically carried out on a small system, an upscaling procedure to much larger super cells is required for device scale simulations. This, along with identification of relevant interactions, is what Winterface provides. With minimal user input and supervision it is able to generate Hamiltonian matrices that ensure efficient Quantum Transport simulations whilst ensuring the relevant physics are captured appropriately. Other more advanced ways of upscaling are supported aswell.

The code is subdivided in 3 parts: A matrix library called *libmat*, a lattice library called *liblat* and a tool exposing the functionality of the former called *ltool*

Documentation generated using doxygen is included.
A publication explaining the code presented here can be found on arXiv (https://arxiv.org/abs/2007.04268)

### Prerequisites

The code depends on BLAS and LAPACK (www.netlib.org), OpenSSL (https://www.openssl.org/) and cppunit (https://www.freedesktop.org/wiki/Software/cppunit/) for its optional test suite.

### Installing

Edit the file 'def.mk' to suit your environment, after which the command 'make' should compile everything. A recent version of g++ is required able to support the C++17 standard. Other compilers may work but are not tested or supported.

## Running the tests

For 'libmat' and 'liblat', a test suite is included found in each directory. After compilation, it can be run by the command 'make run'.

## Authors

**Christian Stieger** - *ETH Zurich* - *Integrated Systems Laboratory* - 2014-2019

## License

This project is licensed under the 2-Clause BSD License (https://opensource.org/licenses/BSD-2-Clause)
