# Lab 8

Consider the relativistic motion of a particle in a strong electrostatic field E(x). In normalized units the equations of motion are

<p align="center">
<img src="stuffy_stuff/f3.png" width="100">
</p>

If we consider a linear electric field E(x)=-x the particle will oscillate about the origin. For initial conditions *x(0) = x0, p(0) = 0* we can integrate the trajectory and determine the oscillation frequency. Due to the nonlinearity the frequency will depend on *x0*. It is the goal of this exercise to obtain a diagram which shows the frequency in dependence of *x0*.


----
For the classic RK4 method

<p align="center">
<img src="stuffy_stuff/rk4.png" width="150">
</p>

a third order interpolation polynomial can be used
to find function values in between two points [t_n, t_{n+1}] via

<p align="center">
<img src="stuffy_stuff/f1.png" width="400">
</p>
where

<p align="center">
<img src="stuffy_stuff/f2.png" width="250">
</p>

Use the bisection method to find the exact position at which one period is over.
