library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(pbapply)
library(here)
load((here("applications/data/coords_danube.Rdata")))

#### Danube data
X <- danube$data
g <- graph_from_edgelist(danube$flow_edges)
g <- graphicalExtremes:::set_graph_parameters(g)
n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
rholist <- seq(.01, .22, length.out=200)

eglearn_fit <- eglearn2(Y, rholist=rholist, reg_method="ns", complete_Gamma=F)
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, graphicalExtremes:::set_graph_parameters)
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, graphicalExtremes:::set_graph_parameters)


# plotting
pdf(file = here("applications/figures/danube.pdf"), width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(g, layout = coords_danube, edge.arrow.size=.3)
dev.off()

pdf(file =  here("applications/figures/danube_aic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$aic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube_bic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$bic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube_mbic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$mbic, layout = coords_danube)
dev.off()

