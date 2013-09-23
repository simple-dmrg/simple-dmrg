Exercises
=========

Day 1
-----

1.  Consider a reduced density matrix :math:`\rho` corresponding to a maximally mixed state in a Hilbert space of dimension :math:`md`.  Compute the truncation error associated with keeping only the largest m eigenvectors of :math:`\rho`.  Fortunately, the reduced density matrix eigenvalues for ground states of local Hamiltonians decay much more quickly!

2.  Explore computing the ground state energy of the Heisenberg model using the infinite system algorithm.  The exact Bethe ansatz result in the thermodynamic limit is :math:`E/L = 0.25 - \ln 2 = -0.443147`.  Note the respectable accuracy obtained with an extremely small block basis of size :math:`m \sim 10`.  Why does the DMRG work so well in this case?

3.  Entanglement entropy:

    (a) Calculate the bipartite (von Neumann) entanglement entropy at the center of the chain during the infinite system algorithm.  How does it scale with :math:`L`?

    (b) Now, using the finite system algorithm, calculate the bipartite entanglement entropy for every bipartite splitting.  How does it scale with subsystem size :math:`x`?

        .. hint::
            To create a simple plot in python::

                >>> from matplotlib import pyplot as plt
                >>> x_values = [1, 2, 3, 4]
                >>> y_values = [4, 2, 7, 3]
                >>> plt.plot(x_values, y_values)
                >>> plt.show()

    (c) From the above, estimate the central charge :math:`c` of the "Bethe phase" (1D quasi-long-range Néel phase) of the 1D Heisenberg model, and in light of that, think again about your answer to the last part of exercise 2.

        The formula for fitting the central charge on a system with open boundary conditions is:

        .. math::

            S = \frac{c}{6} \ln \left[ \frac{L}{\pi} \sin \left( \frac{\pi x}{L} \right) \right] + A

        where :math:`S` is the von Neumann entropy.

        .. hint::
            To fit a line in python::

                >>> x_values = [1, 2, 3, 4]
                >>> y_values = [-4, -2, 0, 2]
                >>> slope, y_intercept = np.polyfit(x_values, y_values, 1)

4.  XXZ model:

    (a) Change the code (ever so slightly) to accommodate spin-exchange anisotropy: :math:`H = \sum_{\langle ij \rangle} \left[ \frac{J}{2} (S_i^+ S_j^- + \mathrm{h.c.}) + J_z S_i^z S_j^z \right]`.

    (b) For :math:`J_z/J > 1` (:math:`J_z/J < -1`), the ground state is known to be an Ising antiferromagnet (ferromagnet), and thus fully gapped.
        Verify this by investigating scaling of the entanglement entropy as in exercise 3.  What do we expect for the central charge in this case?


Day 2
-----

1.  Using ``simple_dmrg_03_conserved_quantum_numbers.py``, calculate the "spin gap" :math:`E_0(S_z=1) - E_0(S_z=0)`.  How does the gap scale with :math:`1/L`?  Think about how you would go about computing the spectral gap in the :math:`S_z=0` sector:  :math:`E_1(S_z=0) - E_0(S_z=0)`, i.e., the gap between the ground state and first excited state *within* the :math:`S_z=0` sector.

2.  Calculate the total weight of each :math:`S_z` sector in the enlarged system block after constructing each block of :math:`\rho`.  At this point, it's important to fully understand *why* :math:`\rho` is indeed block diagonal, with blocks labeled by the total quantum number :math:`S_z` for the enlarged system block.

3.  Starting with ``simple_dmrg_02_finite_system.py``, implement a spin-spin correlation function measurement of the free two sites at each step in the finite system algorithm, i.e., calculate :math:`\langle\vec{S}_{i}\cdot\vec{S}_{i+1}\rangle` for all :math:`i`.  In exercise 3 of yesterday's tutorial, you should have noticed a strong period-2 oscillatory component of the entanglement entropy.  With your measurement of :math:`\langle\vec{S}_{i}\cdot\vec{S}_{i+1}\rangle`, can you now explain this on physical grounds?

    Answer:
    ``finite_system_algorithm(L=20, m_warmup=10, m_sweep_list=[10, 20, 30, 40, 40])`` with :math:`J = J_z = 1` should give :math:`\langle \vec{S}_{10} \cdot \vec{S}_{11} \rangle = -0.363847565413` on the last step.

4.  Implement the "ring term" :math:`H_\mathrm{ring} = K \sum_i S^z_{i} S^z_{i+1} S^z_{i+2} S^z_{i+3}`.  Note that this term is one of the pieces of the SU(2)-invariant four-site ring-exchange operator for sites (:math:`i`, :math:`i+1`, :math:`i+2`, :math:`i+3`), a term which is known to drive the :math:`J_1`-:math:`J_2` Heisenberg model on the two-leg triangular strip into a quasi-1D descendant of the spinon Fermi sea ("spin Bose metal") spin liquid [see http://arxiv.org/abs/0902.4210].

    Answer:
    ``finite_system_algorithm(L=20, m_warmup=10, m_sweep_list=[10, 20, 30, 40, 40])`` with :math:`K = J = 1`, should give :math:`E/L = -0.40876250668`.
