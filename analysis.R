library("rootSolve")
library("rmarkdown")
library("tinytex")
library("nleqslv")
library("numDeriv")
library("deSolve")
library("plant")
source("R/competition.R")
source("R/dd99.R")
source("R/figures.R")
source("R/plant_competition.R")
source("R/plant_figure.R")
source("R/rstar_figure.R")
source("R/rstar_model.R")
source("R/rstar_plot.R")
source("R/rstar_support.R")
source("R/shapes.R")
source("R/utils.R")
dir.create("ms/figures", FALSE, TRUE)
dir.create("ms/figures/shape", FALSE, TRUE)

pdf("ms/figures/kernel.pdf", height = 3.5, width = 5.5, pointsize = 10L)
fig_kernel()
dev.off()

pdf("ms/figures/components.pdf", height = 4L, width = 8L)
fig_components()
dev.off()

pdf("ms/figures/rstar_components_U_1.pdf", height = 4L, width = 8L)
fig_rstar_components(1)
dev.off()

pdf("ms/figures/rstar_components_U_2.pdf", height = 4L, width = 8L)
fig_rstar_components(2)
dev.off()

pdf("ms/figures/rstar_components_UC_1.pdf", height = 4L, width = 8L)
fig_rstar_components(3)
dev.off()

pdf("ms/figures/rstar_components_UC_2.pdf", height = 4L, width = 8L)
fig_rstar_components(4)
dev.off()

pdf("ms/figures/shape/constant.pdf", height = 1.3, width = 2L)
shape_constant()
dev.off()

pdf("ms/figures/shape/gaussian.pdf", height = 1.3, width = 2L)
shape_gaussian()
dev.off()

pdf("ms/figures/shape/platykurtic.pdf", height = 1.3, width = 2L)
shape_platykurtic()
dev.off()

pdf("ms/figures/shape/logistic.pdf", height = 1.3, width = 2L)
shape_logistic()
dev.off()

pdf("ms/figures/shape/gaussian_offset.pdf", height = 1.3, width = 2L)
shape_gaussian_offset()
dev.off()

pdf("ms/figures/shape/exponential.pdf", height = 1.3, width = 2L)
shape_exponential()
dev.off()

pdf("ms/figures/shape/gaussian_with_hat.pdf", height = 1.3, width = 2L)
shape_gaussian_with_hat()
dev.off()

pdf("ms/figures/shape/laplacian.pdf", height = 1.3, width = 2L)
shape_laplacian()
dev.off()

pdf("ms/figures/shape/gaussian_with_spike.pdf", height = 1.3, 
    width = 2L)
shape_gaussian_with_spike()
dev.off()

pdf("ms/figures/shape/step_asymmetric.pdf", height = 1.3, width = 2L)
shape_step_asymmetric()
dev.off()

plant_lma_shared <- plant_competition_prepare("lma", parallel = TRUE)
plant_hmat_shared <- plant_competition_prepare("hmat", parallel = TRUE)

plant_lma_1 <- plant_competition(0.17, plant_lma_shared)
plant_lma_2 <- plant_competition(0.07, plant_lma_shared)
plant_hmat_1 <- plant_competition(5, plant_hmat_shared)
plant_hmat_2 <- plant_competition(15.5, plant_hmat_shared)

pdf("ms/figures/plant_lma_components.pdf", height = 4L, width = 8L)
fig_plant_components(plant_lma_1, plant_lma_2)
dev.off()

pdf("ms/figures/plant_hmat_components.pdf", height = 4L, width = 8L)
fig_plant_components(plant_hmat_1, plant_hmat_2)
dev.off()

pdflatex("ms/competition-kernels.tex")
pdflatex("ms/competition-kernels-sm.tex")
