Day 1 Exercises
===============

1.  Consider a reduced density matrix :math:`\rho` corresponding to a maximally mixed state in a Hilbert space of dimension :math:`md`.  Compute the truncation error associated with keeping only the largest m eigenvectors of :math:`\rho`.  Fortunately, the reduced density matrix eigenvalues for ground states of local Hamiltonians decay much more quickly!

2.  Explore computing the ground state energy of the Heisenberg model using the infinite system algorithm.  The exact Bethe ansatz result in the thermodynamic limit is :math:`E/L = 0.25 - \ln 2 = -0.443147`.  Note the respectable accuracy obtained with an extremely small block basis of size :math:`m \sim 10`.  Why does the DMRG work so well in this case?

3.  Entanglement entropy:

    (a) Calculate the bipartite (von Neumann) entanglement entropy at the center of the chain during the infinite system algorithm.  How does it scale with :math:`L`?

    (b) Now, using the finite system algorithm, calculate the bipartite entanglement entropy for every bipartite splitting.  How does it scale with subsystem size :math:`x`?

        .. hint::
            to create a simple plot in python::

                >>> from matplotlib import pyplot as plt
                >>> x_values = [1, 2, 3, 4]
                >>> y_values = [4, 2, 7, 3]
                >>> plt.plot(x_values, y_values)
                >>> plt.show()

    (c) From the above, estimate the central charge of the "Bethe phase" (1D quasi-long-range Néel phase) of the 1D Heisenberg model, and in light of that, think again about your answer to the last part of exercise 2.

        .. hint::
            to fit a line in python::

                >>> slope, y_intercept = np.polyfit([1, 2, 3, 4], [-4, -2, 0, 2], 1)

        The formula for fitting the central charge on a system with open boundary conditions is:

        .. math::

            S = \frac{c}{6} \ln \left[ \frac{L}{\pi} \sin \left( \frac{\pi x}{L} \right) \right] + A

        where :math:`S` is the von Neumann entropy.

4.  XXZ model:

    (a) Change the code (ever so slightly) to accommodate spin-exchange anisotropy: :math:`H = \sum_{<ij>} \left[ \frac{J}{2} (S_i^+ S_j^- + \mathrm{h.c.}) + J_z S_i^z S_j^z \right]`.

    (b) For :math:`J_z/J > 1` (:math:`J_z/J < -1`), the ground state is known to be an Ising antiferromagnet (ferromagnet), and thus fully gapped.
        Verify this by investigating scaling of the entanglement entropy as in exercise 3.  What do we expect for the central charge in this case?

