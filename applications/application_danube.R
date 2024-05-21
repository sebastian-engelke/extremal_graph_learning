library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(pbapply)
library(here)
library(Matrix)

source("simulations/functions_paper.R")

load("applications/data/coords_danube.Rdata")

#### Danube data
X <- danube$data_clustered[1:214,]

X_test <- danube$data_clustered[215:428,]
g <- graph_from_edgelist(danube$flow_edges)
g <- graphicalExtremes:::set_graph_parameters(g)
n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y_all <- data2mpareto(danube$data_clustered, p)
Y <- data2mpareto(X, p)
Y_test <- data2mpareto(X_test, p)
rholist <- seq(.005, .15, length.out=25)
chi_hat <- emp_chi(data = X, p = p)
Gamma_hat <- emp_vario(data = X, p = p)
chi_hat_test <- emp_chi(data = X_test, p = p)
Gamma_hat_test <- emp_vario(data = X_test, p = p)

# Run eglearn for a suitable list of penalization parameters
#rholist <- seq(0, 0.1, length.out = 11)
eglearn_fit <- eglearn(
    data = Y,
    rholist = rholist,
    complete_Gamma = TRUE
)

# Compute the corresponding likelihoods/ICs
# logliks_eglearn <- sapply(
#     seq_along(rholist),
#     FUN = function(j) {
#         Gamma <- eglearn_fit$Gamma[[j]]
#         if(is.null(Gamma)) return(c(NA, NA, NA))
#         loglik_HR(
#             data = Y,
#             Gamma = Gamma,
#             graph = eglearn_fit$graph[[j]]
#         )
#     }
# )

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

# Get flow graph from package
flow_graph <- getDanubeFlowGraph()

# Fit parameter matrix with this graph structure
flow_Gamma <- fmpareto_graph_HR(data = Y, graph = flow_graph)

# Compute likelihood/ICs and plot parameters
# flow_loglik <- loglik_HR(
#     data = Y,
#     Gamma = flow_Gamma,
#     graph = flow_graph
# )

flow_loglik_test <- loglik_HR(
    data = Y_test,
    Gamma = flow_Gamma,
    graph = flow_graph
)



# ggplot(mapping = aes(x = rholist, y = logliks_eglearn['bic', ])) +
#     geom_line() +
#     geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
#     geom_hline(aes(yintercept = flow_loglik['bic']), lty = "dashed") +
#     # geom_hline(aes(yintercept = loglik_emst['bic']), lty = "dotted") +
#     xlab("rho") +
#     ylab("BIC") +
#     scale_x_continuous(
#         breaks = rholist,
#         labels = round(rholist, 3),
#         sec.axis = sec_axis(
#             trans = ~., breaks = rholist,
#             labels = sapply(eglearn_fit$graph, igraph::gsize),
#             name = "Number of edges"
#         )
#     )


ggplot(mapping = aes(x = rholist, y = logliks_eglearn_test['loglik', ])) +
    geom_line() +
    geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
    geom_hline(aes(yintercept = flow_loglik_test['loglik']), lty = "dashed") +
    # geom_hline(aes(yintercept = loglik_emst['bic']), lty = "dotted") +
    xlab("rho") +
    ylab("BIC") +
    scale_x_continuous(
        breaks = rholist,
        labels = round(rholist, 3),
        sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(eglearn_fit$graph, igraph::gsize),
            name = "Number of edges"
        )
    )




# plot_fitted_params <- function(G0, G1, xlab = 'True', ylab = 'Fitted'){
#     return(
#         ggplot()
#         + geom_point(aes(
#             x = G0[upper.tri(G0)],
#             y = G1[upper.tri(G1)]
#         ))
#         + geom_abline(slope = 1, intercept = 0)
#         + xlab(xlab)
#         + ylab(ylab)
#     )
# }

# # Compare fitted chi to empirical one
# best_index <- 8#which.min(logliks_eglearn['bic', ])
# best_Gamma <- eglearn_fit$Gamma[[best_index]]
# best_graph <- eglearn_fit$graph[[best_index]]
# plotDanubeIGraph(graph = best_graph)
# plot_fitted_params(chi_hat, Gamma2chi(best_Gamma), xlab='Empirical')
# plot_fitted_params(chi_hat_test, Gamma2chi(best_Gamma), xlab='Empirical')

# plot_fitted_params(Gamma_hat_test, best_Gamma, xlab='Empirical')





eglearn_fit <- eglearn(Y_all, rholist=rholist, reg_method="ns", complete_Gamma=FALSE)
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
connected_sub <- sapply(1:length(rholist), function(r)
  is_connected(induced_subgraph(graph=eglearn_fit$graph[[r]], v=c(1:22, 28:31))))
eglearn_fit$graph <- lapply(eglearn_fit$graph, graphicalExtremes:::set_graph_parameters)
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, graphicalExtremes:::set_graph_parameters)


# plotting
pdf(file = here("applications/figures/danube.pdf"), width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(g, layout = coords_danube, edge.arrow.size=.3)
dev.off()

# pdf(file =  here("applications/figures/danube_aic.pdf"), width = 5, height = 5)
# par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
# plot(eglearn_fit$graph_ic$aic, layout = coords_danube)
# dev.off()
# 
# pdf(file =  here("applications/figures/danube_bic.pdf"), width = 5, height = 5)
# par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
# plot(eglearn_fit$graph_ic$bic, layout = coords_danube)
# dev.off()

pdf(file =  here("applications/figures/danube_mbic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$mbic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube_06.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[12]], layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube_09.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[18]], layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube_12.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[25]], layout = coords_danube)
dev.off()

