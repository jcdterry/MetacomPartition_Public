#### Function to calculate R at site 

Calc_RMat <- function(E_nodes, EnvOpt, E_change=0, R_max = 1){


E_nodes<- E_nodes+E_change
Deviations <-outer( EnvOpt,E_nodes, `-`) ## species x site matrix of deviations from optimum 
R_Mat <- R_max*cospi(Deviations)^2      ## species x site matrix of population growth rates [not including sensitivity param yet]

return(R_Mat)
}




Calc_RMat_2 <- function(E_nodes, EnvOpt, E_nodes_2, EnvOpt_2,
                        RelEnvInfl = c(0.5, 0.5), 
                        E_change=0,   R_max = 1){
  
  ### compeletely uncorrelated + muktiplicative envrionmental contributions to R. ()
  ### Both use cos function
  
  E_nodes<- E_nodes+E_change
  Deviations <-outer( EnvOpt,E_nodes, `-`) ## species x site matrix of deviations from optimum 
  R_Mat_1 <-cospi(Deviations)^2      ## species x site matrix of population growth rates [not including sensitivity param yet]
  
  
  ### becuase 2nd env var doesnt change, it doesn't need to be calculated each time, but is also unlikely to take up too much time?
  Deviations <-outer( EnvOpt_2,E_nodes_2, `-`) ## species x site matrix of deviations from optimum 
  R_Mat_2 <- cospi(Deviations)^2      ## species x site matrix of population growth rates [not including sensitivity param yet]

    R_Mat <- R_max * sqrt(((R_Mat_1) * (R_Mat_2) ))
  return(R_Mat)
}

Calc_RMat_3 <- function(E_nodes, EnvOpt, 
                        E_change=0,   R_max = 1,
                        R_Mat_2, R_Mat_3){
  
  ### compeletely uncorrelated + muktiplicative envrionmental contributions to R. ()
  ### Both use cos function
  
  E_nodes<- E_nodes+E_change
  Deviations <-outer( EnvOpt,E_nodes, `-`) ## species x site matrix of deviations from optimum 
  R_Mat_1 <-cospi(Deviations)^2      ## species x site matrix of population growth rates [not including sensitivity param yet]
  
  R_Mat <- R_max * sqrt((R_Mat_1*R_Mat_2))*R_Mat_3
  return(R_Mat)
}




############
## Dispersal
############


TotalDispersal<- function(Disp, B_MAT, DispMat,ImmFrac, N_sp, N_Nodes , ArrivalPop){

## Local Dispersal in
LocalDispersal<-Disp *  ( B_MAT %*% DispMat) ###  species x site matrix of immigration

## Global Dispersal in

N_Arrive = ceiling(N_sp*ImmFrac)

IncomingDispersal <- matrix(0,nrow =N_sp ,ncol =N_Nodes)
Sp_to_add<-sample.int( N_sp,  N_Arrive)

## rather than arriving in just 1 node from outside, now arrive in a random 10% of nodes - 
##  better chance of finding places where it can sustain itself
  
for(jj in 1:N_Arrive){
  ToAddToEachNode <- sample(c(0, ArrivalPop),
                            size= N_Nodes, 
                            prob = c(0.9, 0.1)  ,
                            replace = TRUE)
  IncomingDispersal[Sp_to_add[jj],  ] <-ToAddToEachNode
}

## Update B_MAT with dispersors
return(LocalDispersal + IncomingDispersal)
}

##########
## Core model


Dynamics_Discrete <- function(B_MAT , t_max, timescale,R_Mat,comp_str, CompMat, K){
  
  ##  ### Could first chop B_MAT + co down to size - by removing zeros - May not save that much...
  
  for( t_step in 2:t_max){
    B_MAT <- B_MAT + timescale*B_MAT * (R_Mat   - comp_str * CompMat %*% B_MAT - B_MAT/K)
    B_MAT[B_MAT< 0] <- 0  ## instantly remove negatives (prevents blowups)
  }
  
  if(any(is.na(B_MAT))){break} ## for debugging
  
  
  B_MAT[B_MAT< 10^-12] <- 0   # threshold down any that risk numeric errors

  
  return(B_MAT)
}


Dynamics_Discrete_Sparse <- function(B_MAT , t_max, timescale,R_Mat,comp_str, CompMat, K){
  
  ## converting to a sparse matrix halves the time needed 
  B_MAT <- Matrix::Matrix(B_MAT, sparse = TRUE)
  
  
  for( t_step in 2:t_max){
    B_MAT <- B_MAT + timescale*B_MAT * (R_Mat   - comp_str * CompMat %*% B_MAT - B_MAT/K)
    B_MAT[B_MAT< 0] <- 0  ## instantly remove negatives (prevents blowups)
  }
  
  if(any(is.na(B_MAT))){break} ## for debugging
  
  
  B_MAT[B_MAT< 10^-12] <- 0   # threshold down any that risk numeric errors
  return(B_MAT)
}







