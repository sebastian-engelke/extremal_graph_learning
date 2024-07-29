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

load(here("applications/data/coords_danube.Rdata"))

#### Danube data
ids <- (1:31)[-(23:27)]
X <- danube$data_clustered[,ids]
d <- ncol(X)
flow_edges_26 <- danube$flow_edges[1:25,]
flow_edges_26 <- apply(flow_edges_26, 1:2, function(x){
  if(x <= 22) return(x)
  else return(x-5)
})
g <- graph_from_edgelist(flow_edges_26)
g <- set_graph_parameters(g)
coords_danube <- coords_danube[ids,]

n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
rholist <- seq(.005, .5, length.out=100)

eglearn_fit <- eglearn2(Y, rholist=rholist, reg_method="ns", ic="hr")
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, set_graph_parameters)
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, set_graph_parameters)


# plotting
pdf(file = here("applications/figures/danube26.pdf"), width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(g, layout = coords_danube, edge.arrow.size=.3)
dev.off()

pdf(file =  here("applications/figures/danube26_aic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$aic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube26_bic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$bic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube26_mbic.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph_ic$mbic, layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube26_075.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[15]], layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube26_200.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[40]], layout = coords_danube)
dev.off()

pdf(file =  here("applications/figures/danube26_405.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[81]], layout = coords_danube)
dev.off()



#### model testing
load(here("applications/data/coords_danube.Rdata"))
ids <- (1:31)[-(23:27)]
X <- danube$data_clustered[,ids]
d <- ncol(X)
flow_edges_26 <- danube$flow_edges[1:25,]
flow_edges_26 <- apply(flow_edges_26, 1:2, function(x){
  if(x <= 22) return(x)
  else return(x-5)
})
g <- graph_from_edgelist(flow_edges_26)
g <- set_graph_parameters(g)
coords_danube <- coords_danube[ids,]

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
small_rholist <- .03*(1:3)/4

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

loglik_flow <- loglik_HR(
  data = Y_test,
  Gamma = complete_Gamma(emp_vario(Y), g),
  graph = g
)

# plotting
gg_test <- ggplot(mapping = aes(x = sort(c(rholist, small_rholist)), y = logliks_eglearn_test['loglik', ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = loglik_flow['loglik']), lty = "dashed", lwd=1.25) +
  geom_hline(aes(yintercept = loglik_emst['loglik']), lty = "dotted", lwd=1.25) +
  geom_hline(aes(yintercept = loglik_aic['loglik']), lty = "solid", lwd=1.25, color="pink") +
  geom_hline(aes(yintercept = loglik_bic['loglik']), lty = "solid", lwd=1.25,  color="lightblue") +
  geom_hline(aes(yintercept = loglik_mbic['loglik']), lty = "solid", lwd=1.25, color="darkgrey") +
  xlab(TeX("Tuning parameter $\\rho$")) +
  ylab("Test log-likelihood") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 2),
    sec.axis = sec_axis(
      trans = ~., breaks = rholist,
      labels = sapply(c(1, (length(small_rholist)+2):(length(small_rholist)+length(rholist))), function(i) igraph::gsize(eglearn_test$graph[[i]])),
      name = "Number of edges"
    )
  )

pdf(file = here("applications/figures/danube26_test.pdf"), width = 5, height = 5)
par(cex = .8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
gg_test
dev.off()



#### comparing vanilla glasso and ns
# transforming data to (approximate) marginal normality
n <- nrow(X)
d <- ncol(X)
# Y <- log(X)
Y <- 0*X
for(i in 1:d){
  bc <- invisible(boxcox(lm(X[,i] ~ 1)))
  lambda <- bc$x[which.max(bc$y)]
  Y[,i] <- (X[,i]^lambda - 1)/lambda
}
# testing normality of the new variables
qqnorm(Y[,sample(1:d, 1)])

# vanilla glasso
rholist <- seq(.9, 1, length.out=100)
S <- cov(Y)
R <- diag(diag(S)^(-1/2)) %*% S %*% diag(diag(S)^(-1/2))
glasso_fit <- sapply(rholist, function(r) glassoFast(R, rho=r)$wi, simplify="array")
glasso_graph <- apply(glasso_fit, 3, function(M) igraph::graph_from_adjacency_matrix(abs(M) > 1e-8, mode="undirected", diag=FALSE))
glasso_graph <- lapply(glasso_graph, set_graph_parameters)

# vanilla neighborhood selection
rholist <- c(seq(.2, .99, length.out=100), seq(.995, 1, length.out=100))
ns_fit <- glasso_mb2(data=Y, samp_size=n, lambda=rholist)
ns_graph <- apply(ns_fit$adj.est, 3, function(M) igraph::graph_from_adjacency_matrix(M, mode="undirected", diag=FALSE))
ns_graph <- lapply(ns_graph, set_graph_parameters)

# edge counts of the estimated graphs
unlist(lapply(eglearn_fit$graph, function(g) length(E(g))))
unlist(lapply(glasso_graph, function(g) length(E(g))))
unlist(lapply(ns_graph, function(g) length(E(g))))

# comparing and plotting estimated graphs with 27 edges
length(E(eglearn_fit$graph[[80]]))
length(E(glasso_graph[[66]]))
length(E(ns_graph[[130]]))
pdf(file =  here("applications/figures/comparison_danube26_27.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[80]], layout = coords_danube)
plot(glasso_graph[[66]], layout = coords_danube)
plot(ns_graph[[130]], layout = coords_danube)
dev.off()

# comparing and plotting estimated graphs with 50 edges
length(E(eglearn_fit$graph[[21]]))
length(E(glasso_graph[[27]]))
length(E(ns_graph[[26]]))
pdf(file =  here("applications/figures/comparison_danube26_50.pdf"), width = 5, height = 5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(eglearn_fit$graph[[21]], layout = coords_danube)
plot(glasso_graph[[27]], layout = coords_danube)
plot(ns_graph[[26]], layout = coords_danube)
dev.off()

