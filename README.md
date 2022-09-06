# Oscillon-Nonlinear-simulation

This program is dedicated to efficient parallel computing simulation of spherically symmetric oscillon.

## Quick Start

To run this program, install the openmpi library and execute the following command:

```g++ -fopenmp g72.cpp```

or on mac, please try:

```g++11 -fopenmp g72.cpp```

then run:

```./a.out```

the output file include ``energy.dat``, ``pdrphi.dat``, ``pdtphi.dat``, ``pdtphiflag.dat``, ``phi.dat``, ``phiflag.dat``, ``r.dat``, ``rT.dat`` and ``T.dat``.

## Details

### 1.Initial Setup and Background Introduce

The action of the oscillon we focus on can be written as
$S_6=\f{m}{2\li}\int \d^{3}x\;{\cal L} \equiv \f{m}{2\li}\int \d^{3}x\(-\frac{1}{2}\partial_\mu\phi \partial^\mu\phi-\frac{1}{2}\phi^2+\frac{1}{2}\phi^4-\frac{1}{2}g\phi^6\).$

### 2.Input parameters

``R``: Computational domain.

``zoombeta``: Presicion parameter, the higher, the smaller the ``dr``

``sz``: Number of grids.

``dr``: Space step.

``cfl``: ``dt``/``dr``.

``dt``: Time step.

``tnum``: Number of step.

``tlength``: Total compute time length;

``Dim``: Compute domain dimension.

``g``: Potential parameter.

``A``: Intial condition setup parameter.

``a``: Intial condition setup parameter.

``sigma``: Intial condition setup parameter.

``peroid``: Output peroid, when ``timestep/peroid``=``ceil(timestep/peroid)``, run output function. (NOT THE PEROID OF FIELD!!!!!!!)

``pi``: Yes, it is $\pi$.

``numofrT``: Evenly set ``numofrT`` points in the computational domain to output the field period of the points.
