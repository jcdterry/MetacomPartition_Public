### Load dispersal network

# Network<-make_lattice(  dimvector = c(10, 10),  circular = TRUE) 
# plot.igraph(Network, layout =layout_on_grid)
# Network%>%  as_adjacency_matrix() %>% 
#   as.matrix() -> adjMat
# save(adjMat, file = '../Parameters/Disp10by10loop' )

InitiateSystem <-function( Pars){
  
  set.seed(Pars$seed)
  
  if( Pars$Region_Length== 10){
    load('../Parameters/Disp10by10loop' )
  }else{
    make_lattice(  dimvector = c(Pars$Region_Length, 
                                      Pars$Region_Length),  circular = TRUE) %>%
           as_adjacency_matrix() %>% 
           as.matrix() -> adjMat
    }
  
  N_sp <- Pars$S_global
  N_Nodes <-  Pars$Region_Length^2
  
  ## Generate environmental distributions:
  # ##(uncorrelated)
  # E_nodes   <- runif(N_Nodes)   
  # E_nodes_2 <- runif(N_Nodes)   
  
  ## Directional, but with some Guassian noise 
  # NB Pars$Region_Length = number of nodes, not strictly distance between. However, don't want to double 0 and max
  RL <- Pars$Region_Length
  
  Node_loc<-  data.frame(NodeID = 1:(RL^2),
                         x=rep((1:RL)/RL, each = RL),
                         y=rep((1:RL)/RL, times = RL))
  
  E_nodes   <- rnorm(n=RL^2, mean = Node_loc$x, sd = 0.1)   
  E_nodes_2 <- rnorm(n=RL^2, mean = Node_loc$y, sd = 0.1)
  
  ## Generate Species environmental optima 
  EnvOpt   <- runif(Pars$S_global)
  EnvOpt_2 <- runif(Pars$S_global)
  
  ### becuase 2nd env var doesnt change, it doesn't need to be calculated each time,
  Deviations <-outer( EnvOpt_2,E_nodes_2, `-`) ## species x site matrix of deviations from optimum 
  R_Mat_2 <- cospi(Deviations)^2      ## species x site matrix of population growth rates [not including sensitivity param yet]
  
  ## Generate Totally random site specificty of species (that is fixed through time)
  R_Mat_3 <- matrix(runif( N_sp*N_Nodes, min = 0.2, max =1),
                    nrow = N_sp, ncol = N_Nodes)
  
  ################
  ## generate interactions
  ################
  # 
  # if(is.null(Pars$mutu_conn)){
  #   print('just doing random competiton')
  #  
  # CompMat <- matrix(sample(c(0,1), 
  #                          size = N_sp^2,
  #                          replace= TRUE,
  #                          prob = c(1-Pars$comp_conn,
  #                                   Pars$comp_conn)),
  #                   nrow = N_sp)   
  # }else{
  #     CC <- Pars$comp_conn/2
  #     MC <- Pars$mutu_conn/2
  # 
  #     CompMat <- matrix(sample(c(0,1,-1*Pars$mutu_relstr), 
  #                              size = N_sp^2,
  #                              replace= TRUE,
  #                              prob = c(1-CC-MC,CC, MC)),
  #                       nrow = N_sp)
  #   
  #     CompMat[lower.tri(CompMat)]<-0
  #     CompMat <- CompMat+t(CompMat)
  # }
  ######
  ## Approach using very sparse connectance - just off diagonal and a small sprinkling of randoms set by comp_conn
  
  CompMat <- matrix(sample(c(0,1), 
                           size = N_sp^2,
                           replace= TRUE,
                           prob = c(1-Pars$comp_conn,
                                    Pars$comp_conn)),
                    nrow = N_sp)   
  
  CompMat[row(CompMat) == col(CompMat)-1] <- 1

  
  diag(CompMat) <- 0
  
  ### Starting biomass distribution 
  
  B_MAT <- matrix(0,N_sp,
                  ncol =Pars$Region_Length^2 )
  
  ##### Select species to add
  N_to_add <- ceiling(N_sp/2)
  
  Sp_to_add<-sample.int( N_sp,N_to_add )
  
  Node_to_add <- sample.int(N_Nodes, N_to_add,
                            replace = TRUE)
  
  for(jj in 1:N_to_add){
    B_MAT[Sp_to_add[jj], Node_to_add[jj]] <-1
  }
  
  return(list( DispMat = adjMat, 
               CompMat = CompMat,
               EnvOpt = EnvOpt,
               EnvOpt_2 = EnvOpt_2,
               B_MAT= B_MAT, 
               E_nodes = E_nodes, 
               E_nodes_2= E_nodes_2,
               R_Mat_2 =R_Mat_2,
               R_Mat_3 =R_Mat_3   ))
}

