library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(latex2exp)
library(glmnet)  
library(egg)
library(cowplot)
library(tictoc)
library(here)
library(clusterGeneration)
source(here("simulations/functions_paper.R"))

####################


#################################
##### Plotting graphs      ######
#################################

set.seed(1232)
m <- 1
node_size <- 8
gg <- generate_BA_model(d=100, m=m)$graph
igraph::V(gg)$color <- grDevices::adjustcolor(col = "#4477AA", 
                                                 alpha.f = 0.4)
V(gg)$size <- node_size
V(gg)$label <- ""
pdf(here("figures/tree_graph.pdf"))
plot(gg, layout=layout_with_fr(gg))
dev.off()


m <- 2
gg <- generate_BA_model(d=100, m=m)$graph
igraph::V(gg)$color <- grDevices::adjustcolor(col = "#4477AA", 
                                              alpha.f = 0.4)
V(gg)$size <- node_size
V(gg)$label <- ""
pdf(here("figures/BA_graph.pdf"))
plot(gg, layout=layout_with_fr(gg))
dev.off()


alphad=20
gg <- generate_block_model(ncliques = 6, clique_size = 4, alphad=alphad)$graph #ncliques = 10
igraph::V(gg)$color <- grDevices::adjustcolor(col = "#4477AA", 
                                              alpha.f = 0.4)
V(gg)$size <- node_size
V(gg)$label <- ""

pdf(here("figures/block_graph.pdf"))
plot(gg, layout=layout_with_fr(gg))
dev.off()




#################################
##### Regularization Paths ######
#################################

eglearn_path <- function(d=20, m=2, reg_method, rholist, ...){
  
  BA_model <- generate_BA_model(d=d, m=m)
  g <- BA_model$graph
  G <- BA_model$G
  k <- d * 5
  n <- ceiling(k^{1/0.7})
  p <- 1 - k/n
  X <- rmstable(n=n, d=d, model="HR", par=G)
  
  fit_tmp <- eglearn2(data = X, p=p, rholist = rholist, reg_method = reg_method)
  F1_tmp <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit_tmp$graph[[i]]))
  connected_tmp <- sapply(1:length(rholist), FUN = function(i) ifelse(is_connected(fit_tmp$graph[[i]]), 1, 3))
  sparse_tmp <- sapply(1:length(rholist), FUN = function(i) length(E(fit_tmp$graph[[i]])) / (d*(d-1)/2))
  
  plot(rholist, F1_tmp, type="b", pch = connected_tmp, ylim=c(0, 1), xlab="Penalization parameter rho", ylab = "F score",...)
  par(new=TRUE)
  ## Plot the second plot and put axis scale on right
  plot(rholist, sparse_tmp, pch=15,  xlab="", ylab="", ylim=c(0,1), type="l", lty=2, par=...)
  axis(4, ylim=c(0,1),las=1)
  
  tibble(rholist = rholist, F1 = F1_tmp, sparsity = sparse_tmp, connected = factor(connected_tmp), reg_method = reg_method)
}


set.seed(124124)
rholist <- seq(0.000001, 1, length.out = 30) 
tbl1 <-  eglearn_path(d=20, m=2, reg_method = "ns", rholist=rholist)
rholist <- seq(0.000001, .25, length.out = 30) 
tbl2 <- eglearn_path(d=20, m=2, reg_method = "glasso", rholist=rholist)

tbl3 <- bind_rows(tbl1,tbl2) %>% 
  mutate(reg_method_label = refactor_methods(reg_method, lst_methods))


gg_path <- ggplot(tbl3, aes(col=reg_method_label)) + facet_wrap(~reg_method_label, scales = "free_x") + geom_line(aes(x=rholist, y=F1)) +
  geom_line(aes(x=rholist, y=sparsity), lty = "dashed") +
  geom_point(aes(x=rholist, y=F1, shape = connected), fill="white", size = 2, show.legend = FALSE) +
  ylab("F score / Graph density") +
  xlab(TeX("Tuning parameter $\\rho$")) +
  scale_shape_manual(values = c(21,4)) +  
  scale_fill_manual(values = my_fill[c("gl", "ns")]) +
  scale_colour_manual(values = my_col[c("gl", "ns")]) +
  theme(legend.position = c(0.5, -0.23), legend.direction = "horizontal",  strip.text.x = element_blank()) +
  labs(color="Method:", fill = "Method:")
save_myplot(gg_path, plt_nm = here("figures/ns_gl_path.pdf"), width = 2.5, height = 2.5)
  


#################################
##### Comparison analysis #######
##### trees and BA model  #######
#################################

d <- 20  ## 50, 100

dat_tree <- read_rds(here("figures/data", paste0("sim_study_1_d",d,"_tree.rds"))) %>% 
  unnest(cols = c("perf")) %>%
  mutate(value=if_else(type=="time", log10(abs(value)), value)) %>% 
  pivot_wider(names_from = type, values_from = value) 

dat_BA <- read_rds(here("figures/data", paste0("sim_study_1_d",d,"_BA.rds")))  %>% 
  unnest(cols = c("perf")) %>%
  mutate(value=if_else(type=="time", log10(abs(value)), value)) %>% 
  pivot_wider(names_from = type, values_from = value) 

dat <- bind_rows(dat_tree,dat_BA)
dat <- dat %>% 
  mutate(n_char= factor(n, labels=c("0.5", "1", "2.5", "5"), levels = dat$n %>% sort() %>% unique()))

dat1 <- dat %>% 
  dplyr::select(n_char, reg_method, F1_ns_aic, F1_ns_bic, F1_ns_mbic, m) %>% 
  filter(reg_method=="ns") %>% 
  dplyr::select(-reg_method) %>%  
  pivot_longer(cols = contains("F1"), names_to = "reg_method", values_to = "score") 

dat2 <- dat %>% 
  dplyr::select(n_char,reg_method, F1_max, m) %>% 
  rename(score=F1_max)

dat3 <- bind_rows(dat1,dat2) %>% 
  filter(reg_method %in% c("emst", "MTP2","glasso", "ns", "F1_ns_aic", "F1_ns_bic", "F1_ns_mbic")) %>% 
  mutate(reg_method_label = refactor_methods(reg_method, lst_methods))

gg1 <- ggplot(dat3,  aes(x=n_char, y=score, fill= reg_method_label, col=reg_method_label)) + 
  facet_grid(~ m, scales = "free_x") +
  geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
  xlab("k/d") +
  ylab("F score") +
  scale_fill_manual(values = my_fill) +
  scale_colour_manual(values=my_col) +
  ylim(0,1) +
  theme(legend.position = c(0.5, -0.23), legend.direction = "horizontal",  strip.text.x = element_blank()) +
  labs(color="Method:", fill = "Method:")
save_myplot(gg1, plt_nm = here("figures", paste0("boxplot_d",d,".pdf")), width = 3, height = 3)



#################################
##### Comparison analysis #######
##### block model         #######
#################################

    
  dat <- read_rds(here("figures/data/sim_study_1_d19_block_alpha.rds")) %>%
    unnest(cols = c("perf")) %>%
    mutate(value=if_else(type=="time", log10(abs(value)), value)) %>%
    pivot_wider(names_from = type, values_from = value) %>%
    filter(n %in% c(180, 1800)) %>% # 668
    filter(partial_pos < .5) %>% 
    mutate(n_hi_low = if_else(n==max(n), "2", "1")) %>%
    group_by(reg_method, alphad, n, partial_pos, incoh_pos_GLi, incoh_pos_NSi) %>%
    mutate(reg_method_label = refactor_methods(reg_method, lst_methods)) %>%
    mutate(alphad_char=factor(alphad)) %>%
    mutate(reg_method_n = paste(reg_method_label, n_hi_low, sep = "__"))

  
  g0 <- ggplot(dat %>% filter(reg_method =="ns"),  aes(x=alphad_char, y=F1_max,
                         fill = reg_method_label,
                         col=reg_method_label)) +
    geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
    scale_fill_manual(values = my_fill[c("ns", "gl")]) +
    scale_colour_manual(values = my_col[c("ns", "gl")]) + 
    labs(color="Method:") +
    labs(fill="Method:") 
  
  
  
  g1 <- ggplot(dat %>% filter(reg_method=="glasso"),  aes(x=alphad_char, y=partial_pos)) + 
    geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
    xlab("alphad") +
    ylab(TeX("Positive $\\Theta_{ij}$ proportion")) +
    #ylim(0,.5) +
    xlab("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  g2 <- ggplot(dat,  aes(x=alphad_char, y=F1_max,
                         fill = reg_method_n,
                         col=reg_method_n)) +
    geom_boxplot(alpha=.8, size=.35, outlier.size = .7, show.legend = FALSE) +
    xlab("alphad") +
    ylab("F score") +
    scale_fill_manual(values = create_palette_levels(unique(dat$reg_method_n),
                                                     my_palette_methods)$fills) +
    scale_colour_manual(values=create_palette_levels(unique(dat$reg_method_n),
                                                     my_palette_methods)$cols) +
    ylim(.5,1) +
    xlab("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    labs(color="Method:") +
    labs(fill="Method:") 
  
  
  g3 <- ggplot(dat %>% filter(reg_method=="ns"),  aes(x=alphad_char, y=incoh_pos_NSi, fill= reg_method_label, col=reg_method_label)) + 
    geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
    theme(legend.position = "none") +
    xlab(TeX("$\\alpha$")) +
    scale_fill_manual(values = my_fill) +
    scale_colour_manual(values=my_col) +
    ylab("Incoherence value") +
    labs(color="Method:") +
    labs(fill="Method:") 
    
  
  g4 <- ggplot(dat %>% filter(reg_method=="glasso"),  aes(x=alphad_char, y=incoh_pos_GLi, fill= reg_method_label, col=reg_method_label)) + 
    geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
    theme(legend.position = "none") +
    xlab(TeX("$\\alpha$")) +
    scale_fill_manual(values = my_fill) +
    scale_colour_manual(values=my_col) +
    ylab("Incoherence value") +
    labs(color="Method:") +
    labs(fill="Method:") 
  
  grid1 <- plot_grid(g1 + theme(legend.position = "none"), g2, g3, g4, nrow=2)

  library(ggpubr)
  legend <- get_legend(
    # create some space to the left of the legend
    g0+
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  g5 <- ggarrange(g1, g2, NULL, NULL, g3, g4,
                  nrow = 3, ncol = 2, align = "hv",
                  heights = c(1, -0.1, 1),
                  legend = "bottom", legend.grob = legend)
  
  save_myplot(g5, here("figures/block_alpha.pdf"), height = 5, width = 5)
  
  
  
 
  