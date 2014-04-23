% Comparison of the Tilman model (Monod function for essential
  resource) and the Lotka-Volterra model

Georges Kunstler

# Comparison of Lotka Volterra and R* at equilibrium
This is based on the Chapter 7 of Tilman 1982.

# Lotka Volterra model
\begin{equation}
\frac{dN_1}{dt} =  N_1 \times r \times (1 - N_1/K_1 -\alpha_{21}
\times N_2/K_1)
\end{equation}

\begin{equation}
\frac{dN_2}{dt} =  N_2 \times r \times (1 - N_2/K_2 -\alpha_{12}
\times N_1/K_2)
\end{equation}

at equilibrium $\frac{dN_1}{dt} = 0$ and same for $N_2$.

Thus

\begin{equation}
\bar{N_1} =  K_1 - \alpha_{21} \times \bar{N_2}
\end{equation}

and

\begin{equation}
\bar{N_2} =  K_2 - \alpha_{12} \times \bar{N_1}
\end{equation}

The basic idea of the chapter 7 of Tilman is to take the R* model at
equilibrium and rework out the equation to match the same structure as
these two equations and see what match $K$ and $\alpha$.

## R* model based on Huisman & Weissing equations:

Here I follow the notation of Huisman & Weissing (2001) because this is what we have used so far. This is different from Tilman chapter 7.

\begin{equation}
\frac{dN_i}{dt} =  N_i \times (min(p_{1i}(R_1), p_{2i}(R_2))-m)
\end{equation}

\begin{equation}
\frac{dR_j}{dt} = D \times (S_j-R_i) - \sum_{i=1}^I {C_{ji} N_i
min(p_{1i}(R_1),p_{2i}(R_2))}
\end{equation}

\begin{equation}
p_{ji}(R_j) = \frac{r R_j}{K'_{ji}+R_j}
\end{equation}
I use $K'$ to make the difference between the $K$ in Lotka-Volterra.

j and i in 1, 2.

We assume $D=m$.

A species $i$ is defined by two vectors of parameters (of length two for two resources):
- $K_{.i}$ which is the half saturation constant for each resource ($K_{1i}$ for resource $R_1$ and $K_{2i}$ for resource $R_2$).
- $C_{.i}$ which is the consumption rate (or the content of nutriment) for each resource ($C_{1i}$ for resource $R_1$ and $C_{2i}$ for resource $R_2$).

The resource requirement of species $i$ is defined by:  
$R_{11}^* = \frac{m K_{11}}{r-m}$

I assume that species 1 as the minimum $R^*$ for $R_1$ and species 2 for $R_2$.

## Equilibrium for R*

At equilibrium $\frac{dN}{dt} = 0$ and $\frac{dR}{dt} = 0$.

$min(p_{1i}(R_1), p_{2i}(R_2)) = m$.

# If $R_1$ is the limting factor for both species

The carrying capacity of $N_1$ is defined as the equilibrium population when
single species. Which is:

$0 = D(S_1 - R_{11}^*) - C_{11} \bar{N_1} m$

with 
$R_{11}^* = \frac{m K}{r-m}$

Thus $N_1$ equilibrium for single species case is:
$K_1 = \bar{N_1} = \frac{(S_1 - R_{11}^*)}{C_{11}}$

### Equilibrium with two species

$\bar{N_1} = \frac{S_1-R_{11}^*}{C_{11}} - \frac{C_{12}}{C_{11}} \times \bar{N_2}$
$\bar{N_1} = K_1 - \frac{C_{12}}{C_{11}} \times \bar{N_2}$

Thus 

$K_1 = \frac{(S_1 - R_{11}^*)}{C_{11}}$

and 

$\alpha_{21} = \frac{C_{12}}{C_{11}}$

In the same way:

$K_2 = \frac{(S_1 - R_{12}^*)}{C_{12}}$

and 

$\alpha_{12} = \frac{C_{11}}{C_{12}}$

# If $R_2$ is the limting factor for both species

Applying the same approach as for when $R_1$ is limiting gives:

$K_1 = \frac{(S_2 - R_{21}^*)}{C_{21}}$

and 

$\alpha_{21} = \frac{C_{22}}{C_{21}}$

In the same way:

$K_2 = \frac{(S_2 - R_{22}^*)}{C_{22}}$

and 

$\alpha_{12} = \frac{C_{21}}{C_{22}}$

# If $R_1$ is the limting factor for species 1 and $R_2$ for species 2.

Applying the same approach as for when $R_1$ is limiting gives:

$K_1 = \frac{(S_1 - R_{11}^*)}{C_{11}}$

and 

$\alpha_{21} = \frac{C_{12}}{C_{11}}$

In the same way:

$K_2 = \frac{(S_2 - R_{22}^*)}{C_{22}}$

and 

$\alpha_{12} = \frac{C_{21}}{C_{22}}$

# Determining when $R_1$ or $R_2$ are limiting?

For the species 1 the boundary between case where the species is limited by $R_1$ or $R_2$ is given by:
$\frac{(S_2 - R_{21}^*)}{C_{21}} = \frac{(S_1 - R_{11}^*)}{C_{11}}$

If the first part is smaller than the second part of the previous equation then $R_2$ limited if inverse $R_1$ limited.

If we develop for $R^*$ this gives:
$\frac{(S_2 - \frac{m K'_{21}}{r -m})}{C_{21}} = \frac{(S_1 - \frac{m K'_{11}}{r-m})}{C_{11}}$

To start to explore what kind of shape of competition would be predicted by these equations I will assume that a species is defined by a trait $x$ with $K'_{1i} =C_{1i} = x$ and $K'_{2i} = C_{2i} = 1 - x$. 

Replacing $K'$ and $C$ in the previous equation give the condition for $R_2$ limiting depending on $S_1, S_2$ and $x$.

Thus a species defined by the trait $x$ is $R_2$ limited if 

$\frac{S_2}{S_1} < \frac{1-x}{x}$ and $R_1$ limited if inverse. 

**RICH THIS NEEDS TO BE DOUBLE CHECKED I WILL RE-CHECK THAT WHEN BACK FROM CANBERRA!!!**.


## R function to plot alpha and beta in function of x 

**The R code is commented need to include it via knitr!**

<!-- K.and.alpha.in.fun.x <- function(x, xp, m = 0.25, r = 1, S.1 = 1, S.2 = 1){ -->
<!--   x.R.2.lim.TF   <- (S.2/S.1) < ((1-x)  / x) -->
<!--   xp.R.2.lim.TF <- (S.2/S.1) < ((1-xp)/ xp) -->
<!--   if (x.R.2.lim.TF) { -->
<!-- 	  K.x <- (S.2 - m*(1-x)/(r-m))/(1-x) -->
<!-- 	  alpha.xp.x <- (1-xp)/(1-x) -->
<!--   }else{ -->
<!-- 	  K.x <- (S.1 - m*(x)/(r-m))/(x) -->
<!-- 	  alpha.xp.x <- (xp)/(x) -->
<!--   } -->
<!--   #species x.p -->
<!--   if (xp.R.2.lim.TF) { -->
<!-- 	  K.xp <- (S.2 - m*(1-xp)/(r-m))/(1-xp) -->
<!-- 	  alpha.x.xp <- (1-x)/(1-xp) -->
<!--   }else{ -->
<!-- 	  K.xp <- (S.1 - m*(xp)/(r-m))/(xp) -->
<!-- 	  alpha.x.xp <- (x)/(xp) -->
<!--   } -->
<!--   return(data.frame(K.x = K.x, K.xp = K.xp,  -->
<!--                     alpha.xp.x = alpha.xp.x, alpha.xp.x = alpha.xp.x))  -->
<!-- } -->

<!-- x.vec <- (1:99)/100 -->

<!-- xp.fix <- 0.4 -->
<!-- K.alpha <- do.call('rbind', lapply(x.vec, K.and.alpha.in.fun.x , xp =xp.fix)) -->

<!-- par(mfrow = c(2,2)) -->
<!-- plot(x.vec, K.alpha[['alpha.xp.x']], type ='l',  -->
<!--      xlab = 'trait x', ylab = paste('competitive effect of xp on x (xp =', xp.fix,')')) -->
<!-- plot(x.vec, K.alpha[['K.x']], type ='l',  -->
<!--      xlab = 'trait x', ylab = 'K.x') -->
<!-- plot(x.vec, K.alpha[['alpha.x.xp']], type ='l',  -->
<!--      xlab = 'trait x', ylab = paste('competitive effect of x on xp (xp =', xp.fix,')')) -->
<!-- plot(x.vec, K.alpha[['K.xp']], type ='l',  -->
<!--      xlab = 'trait x', ylab = 'K.xp') -->



# References

Huisman, J., and Weissing, F.J. (2001). Biological conditions for oscillations and chaos generated by multispecies competition. Ecology 82, 2682–2695.

Tilman, D. (1982) Resource Competition and Community Structure. Princeton University Press, Princeton, NJ. 
