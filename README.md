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
