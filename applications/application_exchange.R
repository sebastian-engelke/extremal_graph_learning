library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(pbapply)
library(here)

#### Currency data
# loading the data: X contains the absolute residuals from fitting an Arma(0, 2) - Garch(1, 1)
# process to the negative log returns,
# and g is the minimum spanning tree obtained by Engelke & Volgushev
load(here("applications/data/currency.Rdata"))
load(here("applications/data/coords_exchange.Rdata"))

g <- graphicalExtremes:::set_graph_parameters(g)
n <- nrow(X)
p <- .95 # same as in the tree learning paper
Y <- data2mpareto(X, p)
rholist <- seq(.01, .5, length.out=100)

eglearn_fit <- eglearn(Y, rholist=rholist, reg_method="ns", complete_Gamma=F)
connected_tmp <- sapply(eglearn_fit$graph, is_connected)

eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(graphicalExtremes:::set_graph_parameters(gr))
})



pdf(file = here("applications/figures/exchange.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(g, layout = coords_exchange, edge.arrow.size=.3)
dev.off()

pdf(file = here("applications/figures/exchange_aic.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$aic, layout = coords_exchange)
dev.off()

pdf(file = here("applications/figures/exchange_bic.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$bic, layout = coords_exchange)
dev.off()

pdf(file = here("applications/figures/exchange_mbic.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$mbic, layout = coords_exchange)
dev.off()

