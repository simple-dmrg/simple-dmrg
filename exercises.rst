Day 1 Exercises
===============

1. Consider a reduced density matrix $\rho$ corresponding to a maximally mixed state in a Hilbert space of dimension md.  Compute the truncation error associated with keeping only the largest m eigenvectors of $\rho$.  Fortunately, the reduced density matrix eigenvalues for ground states of local Hamiltonians decay much more quickly!

2. Explore computing the ground state energy of the Heisenberg model using the infinite system algorithm.  The exact Bethe ansatz result in the thermodynamic limit is E/L = 0.25 - ln(2) = -0.443147.  Note the remarkable accuracy obtained with an extremely small block basis of size $m \sim 10$.  Why does the DMRG work so well in this case?

3. Entanglement entropy:

    (a) Calculate the bipartite (von Neumann) entanglement entropy at the center of the chain during the infinite system algorithm.  How does it scale with L?

    (b) Now, using the finite system algorithm, calculate the bipartite entanglement entropy for every bipartite splitting.  How does it scale with subsystem size x?

    (c) From the above, estimate the central charge of the "Bethe phase" (1D quasi-long-range Neel phase) of the 1D Heisenberg model, and in light of that, think again about your answer to the last part of exercise 2.

4. XXZ model:

    (a) Change the model to accomodate spin-exchange anisotropy: $H = \frac{J}{2}(S^+ S^- + \mathrm{h.c.}) + J_z S^z S^z$.
    
    (b) For $J_z/J > 1$ ($J_z/J < 1$), the ground state is known to be an Ising antiferromagnet (ferromagnet), and thus fully gapped.
        Verify this by investigating scaling of the entanglement entropy as in exercise 3.


Day 2 Exercises
===============

1. Calculate the spin (triplet) gap by finding the ground state energy in the $S_z=1$ sector.  How does it scale with 1/L?

2. Calculate the total weight of each sector in the enlarged system block after constructing each block of $\rho$.

3. Starting with simple_dmrg_02_finite_system.py, implement a spin-spin correlation function measurement of the free two sites at each step in the finite system algorithm, i.e., calculate $\langle\vec{S}_{i}\cdot\vec{S}_{i+1}$.

   For L=20, m=50, you should get:

4. Implement the ring term $H_\mathrm{ring} = S^z_{i} S^z_{i+1} S^z_{i+2} S^z_{i+3}$, in 10 minutes.  [Note that this term is one of the pieces of the SU(2)-invariant ring-exchange operator for sites (i, i+1, i+2, i+3).]