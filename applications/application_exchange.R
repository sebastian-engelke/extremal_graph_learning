library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(pbapply)
library(here)
library(Matrix)
library(MASS)
library(glassoFast)
library(latex2exp)

source(here("simulations/functions_paper.R"))

#### Currency data
# loading the data: X contains the absolute residuals from fitting an Arma(0, 2) - Garch(1, 1)
# process to the negative log returns,
# and g is the minimum spanning tree obtained by Engelke & Volgushev
load(here("applications/data/currency.Rdata"))
load(here("applications/data/coords_exchange.Rdata"))

n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
rholist <- seq(.005, 1, length.out=200)

eglearn_fit <- eglearn2(Y, rholist=rholist, reg_method="ns", ic="hr")
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})



# plotting graphs
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



#### model testing
# data splitting
load(here("applications/data/currency.Rdata"))
n <- nrow(X)
n_train <- ceiling(n/2)
X_all <- X
X <- X_all[1:n_train,]
X_test <- X_all[(n_train+1):n,]

# training
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
Y_test <- data2mpareto(X_test, p)
rholist <- seq(0, .3, length.out=11)
small_rholist <- c(.015, .045, .225, .255)

eglearn_test <- eglearn2(Y, rholist=sort(c(rholist, small_rholist)), reg_method="ns", ic="hr")
connected_tmp <- sapply(eglearn_test$graph, is_connected)
eglearn_test$graph <- lapply(eglearn_test$graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})
eglearn_test$graph_ic <- lapply(eglearn_test$graph_ic, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})

# calculating likelihoods of estimated models (including mbic) and the emst
logliks_eglearn_test <- pbsapply(
  seq_along(sort(c(rholist, small_rholist))),
  FUN = function(j) {
    Gamma <- eglearn_test$Gamma[[j]]
    if(any(is.na(Gamma))) return(c(NA, NA, NA))
    loglik_HR(
      data = Y_test,
      Gamma = Gamma,
      graph = eglearn_test$graph[[j]]
    )
  }
)

loglik_aic <- loglik_HR(
  data = Y_test,
  Gamma = eglearn_test$Gamma_ic$aic,
  graph = eglearn_test$graph_ic$aic
)

loglik_bic <- loglik_HR(
  data = Y_test,
  Gamma = eglearn_test$Gamma_ic$bic,
  graph = eglearn_test$graph_ic$bic
)

loglik_mbic <- loglik_HR(
  data = Y_test,
  Gamma = eglearn_test$Gamma_ic$mbic,
  graph = eglearn_test$graph_ic$mbic
)

emst_fit <- emst(data = Y, method = "vario")
loglik_emst <- loglik_HR(
  data = Y_test,
  Gamma = emst_fit$Gamma,
  graph = emst_fit$graph
)

# plotting
gg_test <- ggplot(mapping = aes(x = sort(c(rholist, small_rholist)), y = logliks_eglearn_test['loglik', ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  #geom_hline(aes(yintercept = loglik_emst['loglik']), lty = "dotted", lwd=1.25) +
  geom_hline(aes(yintercept = loglik_aic['loglik']), lty = "solid", lwd=1.25, color="pink") +
  geom_hline(aes(yintercept = loglik_bic['loglik']), lty = "solid", lwd=1.25, color="lightblue") +
  geom_hline(aes(yintercept = loglik_mbic['loglik']), lty = "solid", lwd=1.25, color="darkgrey") +
  xlab(TeX("Tuning parameter $\\rho$")) +
  ylab("Test log-likelihood") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 2),
    sec.axis = sec_axis(
      trans = ~., breaks = rholist,
      labels = sapply(rank(c(rholist, small_rholist))[1:length(rholist)], function(i) igraph::gsize(eglearn_test$graph[[i]])),
      name = "Number of edges"
    )
  )

pdf(file = here("applications/figures/exchange_test.pdf"), width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
gg_test
dev.off()




#### comparing vanilla glasso and ns
# transforming data to (approximate) marginal normality
load(here("applications/data/currency.Rdata"))
load(here("applications/data/coords_exchange.Rdata"))
n <- nrow(X)
d <- ncol(X)
# Y <- log(X)
Y <- 0*X
for(i in 1:d){
  bc <- invisible(boxcox(lm(1+X[,i] ~ 1)))
  lambda <- bc$x[which.max(bc$y)]
  Y[,i] <- ((1+X[,i])^lambda - 1)/lambda
}
# testing normality of the new variables
qqnorm(Y[,sample(1:d, 1)])

# vanilla glasso
rholist <- seq(.25, .75, length.out=100)
S <- cov(Y)
R <- diag(diag(S)^(-1/2)) %*% S %*% diag(diag(S)^(-1/2))
glasso_fit <- sapply(rholist, function(r) glassoFast(R, rho=r)$wi, simplify="array")
glasso_graph <- apply(glasso_fit, 3, function(M) igraph::graph_from_adjacency_matrix(abs(M) > 1e-8, mode="undirected", diag=FALSE))
glasso_graph <- lapply(glasso_graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})

# vanilla neighborhood selection
rholist <- seq(.2, .99, length.out=100)
ns_fit <- glasso_mb2(data=Y, samp_size=n, lambda=rholist)
ns_graph <- apply(ns_fit$adj.est, 3, function(M) igraph::graph_from_adjacency_matrix(M, mode="undirected", diag=FALSE))
ns_graph <- lapply(ns_graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})

# edge counts of the estimated graphs
unlist(lapply(eglearn_fit$graph, function(g) length(E(g))))
unlist(lapply(glasso_graph, function(g) length(E(g))))
unlist(lapply(ns_graph, function(g) length(E(g))))

# comparing estimated graphs with 31 edges
length(E(eglearn_fit$graph[[118]]))
length(E(glasso_graph[[62]]))
length(E(ns_graph[[81]]))
pdf(file =  here("applications/figures/comparison_exchange_31.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[118]], layout = coords_exchange)
plot(glasso_graph[[62]], layout = coords_exchange)
plot(ns_graph[[81]], layout = coords_exchange)
dev.off()

# comparing estimated graphs with 49 edges
length(E(eglearn_fit$graph[[87]]))
length(E(glasso_graph[[41]]))
length(E(ns_graph[[65]]))
pdf(file =  here("applications/figures/comparison_exchange_49.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[87]], layout = coords_exchange)
plot(glasso_graph[[41]], layout = coords_exchange)
plot(ns_graph[[65]], layout = coords_exchange)
dev.off()

