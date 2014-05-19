# Competition kernels plan

## Compute alpha directly from Lotka-Volterra equations

Using either commonly used form of the competitive Lotka-Volterra equations, compute the effectve competition kernels.  This follows directly from the equations, rearranged and solved for alpha.

I think that for a multispecies system, this will involve solving a matrix equation, but that's not too bad?  The problem is for n species there are n^2 unknown alpha coefficients.  So we can do this for single species as it invades, but some work is going to be involved to work this out properly for multispecies systems.  Basically there will be a family of different alpha coefficients that will generate the same dynamics, at a given set of species densities.

## Look at dependence on species density and nonlinearity

Does it make any sense to think about competition kernels if we have to think about $alpha(x_i, x_j, N_i N_j)N_j$ (and even worse for systems with more than two species where the other species really end up in that equation.

Should $N_j$ remain outside this equation?  Probably because that means alpha is the per-capita effect.

One way of thinking about this would be to partition alpha into a density-dependent and density independent term:

$(alpha(x_i, x_j) + alpha*(x_i, x_j, N_j)) N_j$

with the second term being zero in most models.

But what density to evaluate the first model at?

Another way might be to think of the Taylor series around the equilibrium density of $N_j$, then we can write $alpha(x_i, x_j, N_j + delta)$ as

$(alpha(x_i, x_j, N_j) + delta * da/dN + delta^2/2 d2a/dN2)$, which will include the density dependence of competition naturally, but only be reasonable in the vicinity of equilibrium.

## Is there anything useful in the population competition?

Is is nice because it's the same thing as the LV model for a linear system?  For nonlinear systems, should we just be looking at a big pile of alphas?

## dN/dt vs dp/dt

For invasion fitness, dN_i/dt is OK, but what we will often want is dp_i/dt.  That requires mean fitness (the weighted sum of all dN_i/dt values).  Once we have that we can reframe the equations in terms of things like the Price equation and hopefully partition out frequency dependent and independent competition out.
