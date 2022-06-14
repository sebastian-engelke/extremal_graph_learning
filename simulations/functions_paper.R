sim_study <- function(d = 5, 
                          n = 100,
                          p = NULL,
                          method = c("maxstable", "mpareto"), 
                          m = 2, 
                          gen_model = c("BA", "block"),
                          alphad = 1,
                          reg_method = c("ns", "glasso", "emst", "MTP2"),
                          rhostring = "seq(0.01,1.5,length=15)",
                          incoh = FALSE,
                          rng = NULL){
  ## perform a simulation study to measure performance of the EMTP2 block descent algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - n: number of samples.
  ##    - p: threshold probability.
  ##    - method: data generation method.
  ##    - m: the number of edges added in each setp of the Barabasi-Albert model (m=1 is a tree)
  ##    - reg_method: regression method used to estimate the extremal graphical structure.
  ##    - rholist: the list of penality parameters; must be given a string.
  ##    - noise: "none", "iid" or "tree". The type of noise that is added 
  ##        to the simulation of the model. Only use with method="maxstable" !!!
  ## Returns:
  ##  a tibble with time, KKT condition values and duality gap
  
  # check arguments
  
  rholist <-  eval(parse(text = rhostring))
  
  F1 <- numeric(length(rholist))
  incoh_pos <- numeric(2)
  F1_ic <- NA
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  
  if(gen_model=="BA"){
    BA_model <- generate_BA_model(d = d, m = m)
    G <- BA_model$G
    g <- BA_model$graph
  }
  if(gen_model=="block"){
    block_model <- generate_block_model(ncliques = 6, clique_size = 4, alphad=alphad) 
    g <- block_model$graph
    G <- .05*block_model$G
    d <- nrow(G)
  }
  
  T <- Gamma2Theta(G)
  nb.edges <- length(E(g))
  partial_pos <- (sum(sum(T > 1e-6)) - d)/2 / nb.edges
  partial_neg <- sum(sum(T < -1e-6))/2 / nb.edges
  
  
  if(incoh==TRUE){
    incoh_pos <- sapply(Gamma2Inc(G = G, pessimistic = TRUE), FUN = function(x) mean(x))  
  }
  else incoh_pos <- c(NA,NA)
  
  # perform simulation
  if(method=="maxstable")  X <- rmstable(n=n, d=d, model="HR", par=G)
  if(method=="mpareto")  X <- rmpareto(n=n, d=d, model="HR", par=G)
  
  
  if(reg_method=="emst"){
    ptm <- proc.time()[1]
    fit <- emst(data = X, p=p, method = "vario")
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph))
  }
  else if(reg_method=="MTP2"){
    G_emp <- emp_vario(data = X, p = p)
    ptm <- proc.time()[1]
    G_emtp2 <- emtp2(G_emp, tol=1e-6,verbose = FALSE)$G_emtp2 
    time <- proc.time()[1] - ptm
    adj_emtp2 <- (abs(Gamma2Theta(G_emtp2)) >= 1e-4) 
    graph_emtp2 <- igraph::graph_from_adjacency_matrix(adj_emtp2, mode = "undirected", diag = FALSE)
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=graph_emtp2))
  }
  else{
    ptm <- proc.time()[1]
    fit <- eglearn(data = X, p=p, rholist = rholist, reg_method = reg_method, complete_Gamma = FALSE)
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph[[i]]))
    if(reg_method=="ns") F1_ic <- sapply(1:3, FUN = function(i) F1_score(g=g, gest=fit$graph_ic[[i]])) # 
  }

  F1_max <- max(F1)
  
  tbl <- tibble(type = paste0("time"), 
                value = time) %>% 
    bind_rows(tibble(type = paste0("F1_score", 1:length(rholist)),
                     value =  F1)) %>% 
     bind_rows(tibble(type = c("incoh_pos_GLi", "incoh_pos_NSi"),
                      value =  incoh_pos)) %>% 
    bind_rows(tibble(type = c("partial_pos", "partial_neg"),
                     value =  c(partial_pos,partial_neg))) %>% 
     bind_rows(tibble(type = c("F1_max"),
                      value =  F1_max)) %>% 
    bind_rows(tibble(type = c("F1_ns_aic", "F1_ns_bic", "F1_ns_mbic"),
                     value =  F1_ic))
  
  return(tbl)
}



F1_score <- function(g, gest) {
  (2*ecount(intersection(gest, g)))/( 2*ecount(intersection(gest, g)) + 
                                        ecount(intersection(complementer(g), gest)) +
                                        ecount(intersection(g, complementer(gest))))
}



incoherence_pess <- function(Sig, S, Sc){
  if(length(Sc)!=0)
    1 - max(sapply(Sc, function(e) sum(abs( Sig[e, S] %*% solve(Sig[S, S]) ))))
  else
    1
}
  
incoherence <- function(Sig, S, Sc, signs){
  if(length(Sc)!=0)
    1 - max(abs(sapply(Sc, function(e) sum( Sig[e, S] %*% solve(Sig[S, S]) %*% signs) )))
  else
    1
}
  

GLNSi <- function(Sig, cor.scale=TRUE, pessimistic=FALSE){
  
  d <- ncol(Sig)
  if(cor.scale) Sig <- ( diag(Sig)^(-0.5) * diag(d) ) %*% Sig %*% ( diag(Sig)^(-0.5) * diag(d) )
  
  Th <- solve(Sig)
  Th <- Th * (abs(Th) >= 1e-6)
  
  gr <- graph_from_adjacency_matrix(Th != 0, mode = "undirected")
  S <- which(as.vector(gr[1:d]) == 1)
  Sc <- which(as.vector(gr[1:d]) == 0)
  H <- kronecker(Sig, Sig)
  GLi <- incoherence_pess(H, S, Sc)
  
  NSi <- min(sapply(1:d, function(a){
    if(sum(Th[-a,a] != 0) == 0) return(1)
    else{
      S <- which(Th[,a] != 0)
      S <- S[S != a]
      Sc <- which(Th[,a] == 0)
      if(pessimistic) return(incoherence_pess(Sig, S, Sc))
      else{
        signs <- sign(Th[S,a])
        return(incoherence(Sig, S, Sc, signs))
      }
    }
  }))
  #GLi=NA
  return(c(GLi, NSi))
  
}

# If pessimistic = F (the default), NSi is equal to the pessimistic lower bound on the incoherence,
# with the operator norm
# Otherwise NSi is the exact incoherence depending on the signs of the regression coefficients
# GLi is unchanged by pessimistic
Gamma2Inc <- function(G, cor.scale=TRUE, pessimistic=FALSE){
  
  
  Inc <- t(sapply(1:nrow(G), function(m){
    
    Sig <- Gamma2Sigma(G, m)
    return(GLNSi(Sig, cor.scale, pessimistic))
    
  }))
  
  return(list(GLi=Inc[,1], NSi=Inc[,2]))
  
}


# traditional criteria
# is consistent for a fixed design, fixed p
aic <- function(n, p) 2
bic <- function(n, p) log(n)

# modified BIC of Wang & Leng, JRSSB 2009
# it has a weird justification in the paper
mbic <- function(n, p) log(n) * log(log(p))

# generalized information criterion of Fan & Tang, JRSSB 2013
# seems to be consistent as long as log(p) = O(n^a), a<1
# questions: where is variance in their GLM? is deviance really "sum of square residuals" in Gaussian LR?
gic <- function(n, p) log(log(n)) * log(p)

# extended BIC from Chen & Chen, Biometrika 2008, but in the form of Wang & Zhu, JMA 2011
# only consistent if p = O(n^a), a<1, gam >= 1 - 1/2a
ebic <- function(n, p) log(n) + 2*gam*log(p)

# high-dimensional BIC from Wang & Zhu, JMA 2011
# consistent if p>n, log(p) = O(n^a), a<1, gam > 1
hbic <- function(n, p) 2*gam*log(p)



wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Ahttps://uniart1.wixsite.com/uni-artrgs:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_fn, args = sim_fn_args[i, ]) %>% #ML: fixed the name of the last argument. It used to be fun_args
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}


#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



generate_BA_model <- function(d,m){
  g <- sample_pa(n=d, m=m, zero.appeal=1,directed=F)
  W <- as_adj(g, sparse=F) * matrix(runif(d^2,2, 5), nrow=d) #matrix(2 + rexp(d^2, rate = 1), nrow=d) #matrix(runif(d^2,2, 5), nrow=d) #  # 
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  O <- diag(rowSums(W)) - W
  G <- Theta2Gamma(O)
  return(list(G = G, graph = Gamma2graph(G,to_plot = FALSE)))
}



generate_block_model <- function(ncliques, clique_size, alphad = 1){
  kk <- clique_size
  GG <- matrix(NA, ncliques*(kk-1) + 1, ncliques*(kk-1) + 1)
  for(i in 1:ncliques){
    bigS <- rcorrmatrix(kk, alphad = alphad)
    G1 <- Sigma2Gamma(bigS, full = TRUE)
    if(i==1) GG[1:kk, 1:kk] <- G1
    else GG[(kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1), (kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1)] <- G1
  }
  G <- complete_Gamma(GG)
  round(Gamma2Theta(G),2)
  sum(sum(round(Gamma2Theta(G),2) > 0)) - nrow(G)
  sum(sum(round(Gamma2Theta(G),2) < 0))
  
  return(list(G=G, graph = Gamma2graph(G,to_plot = FALSE)))
}

save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}



lst_methods <- list("emst" = "emst",
                    "MTP2" = "emtp2",
                    "glasso" = "gl",
                    "ns" = "ns",
                    "F1_ns_aic" = "ns_aic",
                    "F1_ns_bic" = "ns_bic",
                    "F1_ns_mbic" = "ns_mbic"
)

my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)



my_palette_methods <- list(
  c("reg_method" = lst_methods$emst, "color" = my_palette$red, "fill" = "white"),
  c("reg_method" = lst_methods$glasso, "color" = my_palette$blue, "fill" = "white"),
  c("reg_method" = lst_methods$MTP2, "color" = my_palette$green, "fill" = "white"),
  c("reg_method" = lst_methods$ns, "color" = my_palette$yellow, "fill" = "white"),
  c("reg_method" = lst_methods$F1_ns_aic, "color" = "black", "fill" = my_palette$pink),
  c("reg_method" = lst_methods$F1_ns_bic, "color" = "black", "fill" = my_palette$light_blue),
  c("reg_method" = lst_methods$F1_ns_mbic, "color" = "black", "fill" = my_palette$grey)
) %>% 
  purrr::transpose() %>% 
  as_tibble() %>% 
  unnest(cols = c(reg_method, color, fill)) 

my_col <-  my_palette_methods %>% 
  dplyr::select(reg_method, color) %>% 
  deframe()

my_fill <-  my_palette_methods %>% 
  dplyr::select(reg_method, fill) %>% 
  deframe()


refactor_methods <- function(methods, lst_methods){
  ## character_vector list with mapping -> factor
  ## refactor column with methods
  
  lst_methods <- lst_methods
  
  
  unique_methods <- unique(methods)
  
  new_levels <- names(lst_methods)
  new_labels <- lst_methods %>% unlist() %>% unname()
  
  factor(methods,
         levels = new_levels,
         labels = new_labels)
}


theme_fct <- function(font_size1=11,  font_size2=7.5){
  theme_set(theme_bw() +
              theme(
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.background = element_blank(),
                strip.background = element_rect(fill = "white"),
                plot.caption=element_text(size=font_size2, hjust=0, 
                                          margin=margin(t=15)),
                text = element_text(size = font_size1),
                axis.ticks = element_blank(),
                axis.text = element_text(size = font_size1),
                panel.grid.major = element_line(size = 0.25)
              ) 
            #+
            # theme_cowplot(font_size = 11)
  )
}

theme_fct()


create_palette_levels <- function(reg_method_n, palette_tbl){
  
  str_spl <- strsplit(reg_method_n, "__")
  
  my_tbl <- tibble(
    reg_method = purrr::map_chr(str_spl, function(el){el[1]}),
    level =  purrr::map_chr(str_spl, function(el){el[2]})
  ) %>%
    left_join(palette_tbl, by = "reg_method") %>%
    mutate(fill = if_else(level == "2", color, fill)) %>%
    mutate(reg_method_lev = paste(reg_method, level, sep = "__"))
  
  my_col <-  my_tbl %>%
    dplyr::select(reg_method_lev, color) %>%
    deframe()
  
  my_fill <-  my_tbl %>%
    dplyr::select(reg_method_lev, fill) %>%
    deframe()
  
  
  list(cols = my_col, fills = my_fill)
  
}
