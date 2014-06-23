library(grImport)
library(grid)
library(Revolve)
library(vectoR)
library(tree)

source("R/util.R")
source("R/util-grid.R")
source("R/figures-intro.R")
source("R/figures-rstar.R")
source("R/figures-tree.R")
source("R/common.R")

source("R/tree-model.R")
source("R/circlepack.R")

set_cairo_fonts_kreon()

res_rstar_1 <- prep_rstar_1()
res_rstar_2 <- prep_rstar_2(S=c(.5, .5))

## Need to bring the height at maturity in a bit so that it's not so
## stupid.  We reach an alpha value of under 3 at hmat = 2.5, so
## that's a reasonable place to start.
res_tree_lma <- prep_tree_lma()
res_tree_hmat <- prep_tree_hmat()

to_pdf(fig_gaussian_competition(1),
       "figs/gaussian_competition_1.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_gaussian_competition(2),
       "figs/gaussian_competition_2.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_gaussian_competition(3),
       "figs/gaussian_competition_3.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_gaussian_competition(4),
       "figs/gaussian_competition_4.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_gaussian_competition(5),
       "figs/gaussian_competition_5.pdf", width=4.5, height=3, cairo=TRUE)

to_pdf(fig_bridge(1),
       "figs/bridge_1.pdf", width=4.5, height=3)

to_pdf(fig_rstar_alpha(res_rstar_1, 1),
       "figs/rstar1_alpha_1.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_rstar_alpha(res_rstar_1, 2),
       "figs/rstar1_alpha_2.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_rstar_alpha(res_rstar_1, 3),
       "figs/rstar1_alpha_3.pdf", width=4.5, height=3, cairo=TRUE)

## TODO: Relabel x axis as specificicity of resource use (demand more
## of A, more of B, and equal demands).
##
## TODO: Why are we always sensitive to invasion to the right?  Tweak
## the parameters to prevent this.  I wonder in this model if invasion
## is always possible when interspecific competition is less than
## intraspecific competition?
##
## Look into this tomorrow.
to_pdf(fig_rstar_alpha(res_rstar_2, 1),
       "figs/rstar2_alpha_1.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_rstar_alpha(res_rstar_2, 2),
       "figs/rstar2_alpha_2.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_rstar_alpha(res_rstar_2, 3),
       "figs/rstar2_alpha_3.pdf", width=4.5, height=3, cairo=TRUE)
to_pdf(fig_rstar_alpha(res_rstar_2, 4),
       "figs/rstar2_alpha_4.pdf", width=4.5, height=3, cairo=TRUE)

to_pdf(fig_tree_model(),
       "figs/tree_model.pdf", width=4, height=4, cairo=TRUE)

to_pdf(fig_tree(res_tree_lma, 1, "alpha"),
       "figs/tree_alpha_lma_1.pdf", width=5, height=3.5, cairo=TRUE)
to_pdf(fig_tree(res_tree_lma, 2, "alpha"),
       "figs/tree_alpha_lma_2.pdf", width=5, height=3.5, cairo=TRUE)

to_pdf(fig_tree_hmat(),
       "figs/tree_hmat.pdf", width=3, height=4, cairo=TRUE)

to_pdf(fig_tree(res_tree_hmat, 1, "alpha"),
       "figs/tree_alpha_hmat_1.pdf", width=5, height=3.5, cairo=TRUE)
to_pdf(fig_tree(res_tree_hmat, 2, "alpha"),
       "figs/tree_alpha_hmat_2.pdf", width=5, height=3.5, cairo=TRUE)

to_pdf(fig_tree(res_tree_lma, 2, "w"),
       "figs/tree_w_lma.pdf", width=5, height=3.5, cairo=TRUE)
to_pdf(fig_tree(res_tree_hmat, 2, "w"),
       "figs/tree_w_hmat.pdf", width=5, height=3.5, cairo=TRUE)
