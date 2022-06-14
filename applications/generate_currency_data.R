library(igraph)
library(graphicalExtremes)
library(timeSeries)
library(fGarch)
library(tidyverse)
library(here)

# download data from GitHub
load(url("https://github.com/sebastian-engelke/extremal_tree_learning/raw/master/application/data/BoE_exchange_rates.RData"))
data <- BoE_exchange_rates


# create the data: fit an Arma(0, 2) - Garch(1, 1) process to the negative log returns
# and extract the absolute values of the residuals
x <-  apply(data, 2, FUN = function(x) diff(log(x)))
d <-  dim(x)[2]
n <-  dim(x)[1]
fit.fGarch <- list()
residus <-  matrix(NA,n,d)
colnames(residus) <-  colnames(x)
for(i in 1:d){
  ga <- NULL
  form <- paste("~arma(0,2)+garch(1,1)")
  ga <- garchFit(formula=as.formula(form), data=x[,i], trace=FALSE, cond.dist="norm")
  residus[,i] <- ga@residuals/ga@sigma.t
  fit.fGarch[[i]] <- ga
}
data <- abs(residus)

# calculate the minimum spanning tree as in the Engelke & Volgushev paper
# and use their visual style
p <- 0.95 # same threshold as in their paper
G.est = emp_vario(data=data, p = p)
graph.full <- make_full_graph(d)
MST.est <- igraph::mst(graph=graph.full, weights = 2- Gamma2chi(G.est[ends(graph.full,E(graph.full))]), algorithm = "prim")
vertex_attr(MST.est) <- list(name = colnames(data))
MST.est <- graphicalExtremes:::set_graph_parameters(MST.est)

# save the data as X, plotting coordinates as coords_tree and coords, and minimum spanning tree as g
X <- data
g <- MST.est
save(X, g, file=here("applications/data/currency.Rdata"))

