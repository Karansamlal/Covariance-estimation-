##################################################
## Project: Financial case study APG: Covariance estimation methods
## Script purpose: BLOCK-DECO GARCH
## Author: Karan Samlal


# Set up package
library(Hmisc)

# define groups:  metals, energies, grains, livestocks and tropical softs
metals = c("PL","SI","GC","LT","LN","LL","LA","LX","LP") 
energies = c("NG","QS","RB","HO","CO","CL") 
grains = c("W","KW","S","C") 
livestocks = c("LH","LC","FC") 
tropical_softs = c("SB","KC","CT","CC","JO") 



flattenCorrMatrix <- function(cormat) {
  # Gets correlations from a correlation matrix and creates an flat dataframe of it
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}


BLOCK.DECO.GARCH <- function(Cor_mat){
  # convert the correlation matrix to a flat dataframe
  flat.cor = flattenCorrMatrix(Cor_mat)
  
  # convert the correlation values to numeric values
  flat.cor[,3] = as.numeric(flat.cor[,3])
  
  cor.flat = cbind(flat.cor$row, flat.cor$column, flat.cor$cor)
  cor.flat = flat.cor
  cor.flat[,3] = as.numeric(cor.flat[,3])
  cor.flat = data.frame(cor.flat)
  names(cor.flat) <- c("X1", "X2", "X3")
  
  # store the equicorrelations
  rho_metals = c()
  rho_energies = c()
  rho_grains = c()
  rho_livestocks = c()
  rho_tropical_softs = c()
  
  # fill the cross-equicorrelations
  c1 = 1
  c2 = 1
  c3 = 1
  c4 = 1
  c5 = 1
  #print(cor.flat)
  
  for (j in 1:length(cor.flat$X1)){
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% metals){
      rho_metals[c1] = as.numeric(cor.flat$X3[j])
      #print(cor.flat$X3[j])
      c1 = c1 + 1
    }
    if (cor.flat$X1[j] %in% energies  & cor.flat$X2[j]  %in% energies ){
      rho_energies[c2] = as.numeric(cor.flat$X3[j])
      c2 = c2 + 1
    }
    if (cor.flat$X1[j] %in% grains & cor.flat$X2[j]  %in% grains){
      rho_grains[c3] = as.numeric(cor.flat$X3[j])
      c3 = c3 + 1
    }
    if (cor.flat$X1[j] %in% livestocks & cor.flat$X2[j]  %in% livestocks){
      rho_livestocks[c4] = as.numeric(cor.flat$X3[j])
      c4 = c4 + 1
    }
    if (cor.flat$X1[j] %in% tropical_softs & cor.flat$X2[j]  %in% tropical_softs){
      rho_tropical_softs[c5] = as.numeric(cor.flat$X3[j])
      c5 = c5 + 1
    }
    
  }
  
  # get correlation in groups
  mean_rho_metals = mean(rho_metals)
  mean_rho_energies = mean(rho_energies)
  mean_rho_grains = mean(rho_grains)
  mean_rho_livestocks= mean(rho_livestocks)
  mean_rho_tropical_softs= mean(rho_tropical_softs)
  
  ### get correlation between groups
  
  # start with metals
  mean_rho_metals_energies = mean(c(rho_metals, rho_energies))
  mean_rho_metals_grains = mean(c(rho_metals, rho_grains))
  mean_rho_metals_livestocks= mean(c(rho_metals, rho_livestocks))
  mean_rho_metals_tropical_softs = mean(c(rho_metals, rho_tropical_softs))
  
  # then energies
  mean_rho_energies_grains = mean(c(rho_energies, rho_grains))
  mean_rho_energies_livestocks = mean(c(rho_energies, rho_livestocks))
  mean_rho_energies_tropical_softs = mean(c(rho_energies, rho_tropical_softs))
  
  # then grains
  mean_rho_grains_livestocks = mean(c(rho_grains, rho_livestocks))
  mean_rho_grains_tropical_softs = mean(c(rho_grains, rho_tropical_softs))
  
  # then live stocks
  mean_rho_livestocks_tropical_softs = mean(c(rho_livestocks, rho_tropical_softs))
  
  new_cor = c()
  for (j in 1:length(cor.flat$X1)){
    ### first check if we two commodities in the same group
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% metals){
      new_cor[j] = mean_rho_metals
    }
    if (cor.flat$X1[j] %in% energies  & cor.flat$X2[j]  %in% energies ){
      new_cor[j] = mean_rho_energies
    }
    if (cor.flat$X1[j] %in% grains & cor.flat$X2[j]  %in% grains){
      new_cor[j] = mean_rho_grains
    }
    if (cor.flat$X1[j] %in% livestocks & cor.flat$X2[j]  %in% livestocks){
      new_cor[j] = mean_rho_livestocks
    }
    if (cor.flat$X1[j] %in% tropical_softs & cor.flat$X2[j]  %in% tropical_softs){
      new_cor[j] = mean_rho_tropical_softs
    }
    
    ### now check for cross correlations
    
    ### start with metal
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% energies | 
        cor.flat$X2[j] %in% metals & cor.flat$X1[j] %in% energies){
      new_cor[j] = mean_rho_metals_energies
    }
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% grains |
        cor.flat$X2[j] %in% metals & cor.flat$X1[j] %in% grains ){
      new_cor[j] = mean_rho_metals_grains
    }
    
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% livestocks | 
        cor.flat$X2[j] %in% metals & cor.flat$X1[j]  %in% livestocks){
      new_cor[j] = mean_rho_metals_livestocks
    }
    
    if (cor.flat$X1[j] %in% metals & cor.flat$X2[j]  %in% tropical_softs
        | cor.flat$X2[j] %in% metals & cor.flat$X1[j]  %in% tropical_softs){
      new_cor[j] = mean_rho_metals_tropical_softs
    }
    
    ### now energies
    if (cor.flat$X1[j] %in% energies & cor.flat$X2[j]  %in% grains |
        cor.flat$X2[j] %in% energies & cor.flat$X1[j] %in% grains ){
      new_cor[j] = mean_rho_energies_grains
    }
    
    if (cor.flat$X1[j] %in% energies & cor.flat$X2[j]  %in% livestocks | 
        cor.flat$X2[j] %in% energies & cor.flat$X1[j]  %in% livestocks){
      new_cor[j] = mean_rho_energies_livestocks
    }
    
    if (cor.flat$X1[j] %in% energies & cor.flat$X2[j]  %in% tropical_softs
        | cor.flat$X2[j] %in% energies & cor.flat$X1[j]  %in% tropical_softs){
      new_cor[j] = mean_rho_energies_tropical_softs
    }
    
    
    ### then grains
    if (cor.flat$X1[j] %in% grains & cor.flat$X2[j]  %in% livestocks | 
        cor.flat$X2[j] %in% grains & cor.flat$X1[j] %in% livestocks){
      new_cor[j] = mean_rho_grains_livestocks
    }
    if (cor.flat$X1[j] %in% grains & cor.flat$X2[j]  %in%  tropical_softs |
        cor.flat$X2[j] %in% grains & cor.flat$X1[j] %in% tropical_softs ){
      new_cor[j] = mean_rho_grains_tropical_softs
    }
    
    ### then livestocks
    if (cor.flat$X1[j] %in% livestocks & cor.flat$X2[j]  %in% tropical_softs | 
        cor.flat$X2[j] %in% livestocks & cor.flat$X1[j] %in% tropical_softs){
      new_cor[j] = mean_rho_livestocks_tropical_softs
    }
    
  }
  cor.flat$x4 = new_cor
  
  eq_cor =  diag(27)
  # set the equicorrelation to the mean rho at time i
  eq_cor[upper.tri(eq_cor,diag=F)] <-cor.flat$x4
  eq_cor[lower.tri(eq_cor,diag=F)] <- cor.flat$x4
  eq_cor
  return(eq_cor)
  
}
