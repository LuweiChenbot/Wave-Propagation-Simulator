## Wave-Propagation-Simulator
For theoretical background and modeling details, see:
👉 [📄 WavePropagationIntroduction.pdf](./WavePropagationIntroduction.pdf)
### Introduction
Mathematical model for wave propagation are often entangled in implicit nonlinear partial
differential equations with complex spatial and temporal derivatives. Often, the complexity 
of these models increases with that of the analyzed surfaces. However, the propagation
in a stiff vibrating steel string presents a uniquely ideal case. This string can be effectively
described using up to second-order time variables and well-defined boundary conditions
Bensa (2003). This essay explores the implementation of finite-difference schemes to this
problem, benefiting from the linear nature of the proposed model. Theoretical background
and rationale for this particular methodology will be first introduced. Then, an implementation 
is presented. Finally, a discussion on outcomes and limitations are provided.

### References
Bensa, B. S. K.-M. R. . S. J. O. r., J. (2003). The simulation of piano string vibration: from
physical models to finite difference schemes and digital waveguides. 
The Journal of
the Acoustical Society of America, 114(2), 1095–1107. doi: 10.1121/1.1587146
