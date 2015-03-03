# History of competition models

There are a set of papers by MacArthur in the mid-late 70s (@MacArthur-1964, @MacArthur-1967, @MacArthur-1969, @MacArthur-1970)

## @MacArthur-1964

This paper is the precursor for the @MacArthur-1967 paper, but the presentation is quite different.  It seems mostly inspired by Levins' book (Evolution in changing environments; 1968) with the fine-grained/coarse-grained species and fitness sets present, but also shares a lot with Tilman's later work on resource competition (e.g. @Tilman-1980).

The basic idea remains:

> It is well known that related species often differ in either habitat or size, and thereby avoid competitive elimination...Among species which specialize on a single uniform resources, only the most effective one will survive.

The modelling proceeds much more like the R^* model, with differential equations for species and for resources ($dN_i/dt$, $dR_j/dt$), but the equations for change in resources turn out not to matter in their analysis.

## @MacArthur-1967

This is the paper that coined the term "limiting similarity", and I think the first to introduce continuous competition kernels.

Starting with Lotka Volterra equations

\begin{equation}
\frac{d N_i}{d t} = r_i N_i \frac{Ki - N_i - \sum_j \alpha_{ij} N_i}{K_i}
\end{equation}

where $\alpha_{ij}$ is the competition *experienced* by species $i$, *exerted* by $j$ (in their words "the relative depression in $r_i$, caused by an indiviual of species $j$ compared to another of species $i$").

The focus is on conditions for invasion (which is in turn governed by $K_i > \sum_j \alpha_{ij} N_j$.

The link to resources comes via figure 1 where they argume that $U_i(R)$ is a function over resource $R$ that measures that probability that an item of of resource is consumed in unit time by an individual of species $i$, and therefore that the probability of two species simultaneously trying for the same resource is $U_i(R)U_j$(R)$, so

\begin{equation}
\alpha_{ij} = \frac{\int U_i(R)U_j(R) dR}{\int (U_i(R))^2 dR}
\end{equation}

with the form being "more complicated" if the resources are not instantly replenished (but assert that will not be required in their presentation).

Because the equation above is a convolution, if the $U$ functions are Gaussian then the $\alpha$ function is also Gaussian.

% Still have to read through MacArthur-1969, MacArthur-1970

## Extensions of LV competition models

* @Roughgarden-1974: Effects of leptokurtic and platykurtic resource kernels on ability of species to invade communities; competition kernels are derived from the utilisation kernels by convolution.
* @Holt-1985: considers density independent mortality (term $-m_i N_i$ added to $dN_i/dt$), which has no qualitative effect, reducing equilibrium density by a proportional amount, and some theoretical treatment of nonlinear competition functions (that end up subsuming all of $r_i$ and $K_i$, too).
* @Ayala-1973: Alternative model to the LV model, with various extensions that change linearities of interactions.

## Criticism of LV models

* Become cumbersome, non-identifiable @Tilman-1987.  Also notes that the effect of viewing things in a LV framework leads to detecting little competition (or little apparent competition) in communities that should have widespread competition.
* Criticism of constant/linear assumptions in @Schoener-1974 (see p 336)
* Criticism of constant per-capita competition in Wilbur 1972, Ayala et al. 1973, Neill 1974, @Schoener-1974, Smith-Gill and Gill 1978, Abrams 1980 (references in @Holt-1985).

## Other assorted bits

* @Pigolotti-2007, @Pigolotti-2010 show how species clustering happens in Lotka-Volterra models
* @Scheffer-2006: species can coexist by being sufficiently different or sufficiently similar, and "lumpy" distributions can emerge.  They use a competition function that is gaussian.
* @Bender-1984: shows what experiments would be required to infer competition coefficients from data; for $n$ species, $n$ experiments will be required.
* @Schoener-1974: discusses links between resource utilisation kernels and competition coefficients.

# Recent work

* @Mirrahimi-2014: Connections bbetween Tilman style models and LV style models, but looks fairly inpenetrable.  A lot of the motivations look similar to ours.

# References
