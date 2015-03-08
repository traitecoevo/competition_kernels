library(tree2)
source("tree_new.R")

options(mc.cores=3L)
obj <- tree_competition_prepare("lma", parallel=TRUE)

ans <- tree_competition(trait_matrix(0.2, "lma"), obj, n=200)

plot(ans$x_mutant, ans$w, type="l", log="x")
plot(ans$x_mutant, ans$alpha, type="l", log="x")
