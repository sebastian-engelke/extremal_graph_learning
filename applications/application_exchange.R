library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(pbapply)
library(here)
library(Matrix)
library(latex2exp)


source(("simulations/functions_paper.R"))

#### Currency data
# loading the data: X contains the absolute residuals from fitting an Arma(0, 2) - Garch(1, 1)
# process to the negative log returns,
# and g is the minimum spanning tree obtained by Engelke & Volgushev
load(("applications/data/currency.Rdata"))
load(("applications/data/coords_exchange.Rdata"))

g <- set_graph_parameters(g)

X_all <- X
X <- X_all[,]
X_test <- X_all[3001:3790,]


n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y_all <- data2mpareto(X_all, p)
Y <- data2mpareto(X, p)
Y_test <- data2mpareto(X_test, p)
rholist <- seq(.005, 1, length.out=200)


eglearn_fit <- eglearn(Y, rholist=rholist, reg_method="ns", complete_Gamma=TRUE)
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})


plot(eglearn_fit$graph_ic$mbic, layout = coords_exchange)


logliks_eglearn_test <- sapply(
    seq_along(rholist),
    FUN = function(j) {
        Gamma <- eglearn_fit$Gamma[[j]]
        if(is.null(Gamma)) return(c(NA, NA, NA))
        loglik_HR(
            data = Y_test,
            Gamma = Gamma,
            graph = eglearn_fit$graph[[j]]
        )
    }
)


loglik_mbic <- loglik_HR(
    data = Y_test,
    Gamma = eglearn_fit$Gamma_ic$mbic,
    graph = eglearn_fit$graph_ic$mbic
)

# Fit tree graph to the data
emst_fit <- emst(data = Y, method = "vario")

# Compute likelihood/ICs, and plot fitted graph, parameters
loglik_emst <- loglik_HR(
    data = Y_test,
    Gamma = emst_fit$Gamma,
    graph = emst_fit$graph
)



gg_test <- ggplot(mapping = aes(x = rholist, y = logliks_eglearn_test['loglik', ])) +
    geom_line() +
    geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
    #geom_hline(aes(yintercept = loglik_emst['loglik']), lty = "dotted") +
    geom_hline(aes(yintercept = loglik_mbic['loglik']), lty = "solid") +
    xlab(TeX("Tuning parameter $\\rho$")) +
    ylab("Test log-likelihood") +
    scale_x_continuous(
        breaks = rholist,
        labels = round(rholist, 3),
        sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(eglearn_fit$graph, igraph::gsize),
            name = "Number of edges"
        )
    )

gg_test



pdf(file = "applications/figures/exchange_test.pdf", width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
gg_test
dev.off()




pdf(file = ("applications/figures/exchange.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(g, layout = coords_exchange, edge.arrow.size=.3)
dev.off()

# pdf(file = here("applications/figures/exchange_aic.pdf"), width = 5, height = 5)
# par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
# plot(eglearn_fit$graph_ic$aic, layout = coords_exchange)
# dev.off()
# 
# pdf(file = here("applications/figures/exchange_bic.pdf"), width = 5, height = 5)
# par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
# plot(eglearn_fit$graph_ic$bic, layout = coords_exchange)
# dev.off()

pdf(file = here("applications/figures/exchange_mbic.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$mbic, layout = coords_exchange)
dev.off()

pdf(file = here("applications/figures/exchange_35.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[70]], layout = coords_exchange)
dev.off()

pdf(file = here("applications/figures/exchange_45.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[90]], layout = coords_exchange)
dev.off()

pdf(file = here("applications/figures/exchange_53.pdf"), width = 5, height = 5)
par(cex = 0.6, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[106]], layout = coords_exchange)
dev.off()

