# Oscillon-Nonlinear-simulation

This program is dedicated to efficient parallel computing simulation of spherically symmetric oscillon.

## Quick Start

To run this program, install the openmpi library and execute the following command:

```g++ -fopenmp oscillon-simulation.cpp```

or on mac, please try:

```g++-12 -fopenmp oscillon-simulation.cpp```

then run:

```./a.out```

the output file include ``energy.dat``, ``pdrphi.dat``, ``pdtphi.dat``, ``pdtphiflag.dat``, ``phi.dat``, ``phiflag.dat``, ``r.dat``, ``rT.dat`` and ``T.dat``.

## Details

### 1.Background Introduce and Intial Setup

#### Bacground (For more detials, please see arxiv:******* )

The action of the oscillon we focus on can be written as:

$$
S_6=\int d^{3}x\(-\frac{1}{2}\partial_\mu\phi \partial^\mu\phi-\frac{1}{2}\phi^2+\frac{1}{2}\phi^4-\frac{1}{2}g\phi^6\).
$$

Then the EoM becomes:

$$
\ddot\phi- \nabla^2 \phi +\phi -2\phi^3+3g\phi^5=0.
$$

#### Initialize

This code sets grids on `r=(i+0.5)*dr`, where `i=[0,sz]` to avoid the singularity at `r=0`.

If does not have `phiflag.dat` and `pdtphiflag.dat` files, please set the code line 322-323 to be:

```
	initialize(phi,pdtphi,r);
	//continue_initialize(phi,pdtphi);
```

Intial condition is setted to be:

$$
\phi=Ae^{-\frac{(r^2-a)^2}{\sigma^2}},\ \partial_t\phi=0.
$$

This is an empirically good initial condition for oscillon formation of cascading levels. You can also change into other format in `initialize` function.

If has `phiflag.dat` and `pdtphiflag.dat` files, please set the code line 322-323 to be:

```
	//initialize(phi,pdtphi,r);
	continue_initialize(phi,pdtphi);
```

In this way, the program will be initialized with the input `**flag.dat` file and start computing. This allows you to arbitrarily set initial conditions.

### 2.Input parameters

``R``: Computational domain.

``calR``: A optional total energy calculate domain, see section5-`energy.dat`.

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

### 3.Absorb Boundary Conditions

This code use absorb boundary conditions to simulate a single oscillon. For more detials see arxiv:****.

### 4.Computing

The program uses the rk4 (Runge–Kutta 4th-order) method to simulate the $\phi$.

### 5.Output File

``energy.dat``: This file outputs the total energy, the first column outputs the time, the second column outputs the full computational domain total energy, and the third column outputs the total energy in the $r<$`calR` space.

``r.dat``: This file outputs the position $r$ of all grid points.

``pdrphi.dat``: This file outputs the field $\partial_r\phi$ value on all grids, the first column of each row outputs the time, and thereafter the field $\partial_r\phi$ value at that time point.

``pdtphi.dat``: This file outputs the field $\partial_t\phi$ value on all grids, the first column of each row outputs the time, and thereafter the field $\partial_t\phi$ value at that time point.

``phi.dat``: This file outputs the field $\phi$ value on all grids, the first column of each row outputs the time, and thereafter the field $\phi$ value at that time point.

``pdtphiflag.dat``: This file outputs the breakpoint file of $\phi$, and the breakpoint file can be used to continue the calculation in the next calculation. The output format is the same as `phi.dat` but without the first column which outputs time.

``phiflag.dat``: This file outputs the breakpoint file of $\partial_t\phi$, and the breakpoint file can be used to continue the calculation in the next calculation. The output format is the same as `pdtphi.dat` but without the first column which outputs time.

``rT.dat``: This file outputs the position $r$ of the grid points which will output peroid of those points.

``T.dat``: This file outputs the $\phi$ field period corresponding to the points in the `rT.dat` file.
