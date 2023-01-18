
#' @export
nb_lik = function(x,mu,disp){
  #returns negative log likelihood: Var = mu + mu^2/size
  return(-sum(dnbinom(x=x, size = disp, mu = mu,log = TRUE)))
}


#' @export
T_to_p = function(T_stat,alternative = 'two.sided'){
  if(alternative == 'positive'){
    p = 1-pnorm(T_stat)
  }
  if(alternative=='negative'){
    p = 1-pnorm(-T_stat)
  }
  if(alternative=='two.sided'){
    p = 1-pnorm(abs(T_stat))
  }
  return(p)
}

#' @export
gene_level = function(p,w = rep(1,length(p))){
  p = tan((0.5-p)*pi)
  p_total = weighted.mean(p,w,na.rm = T)
  p_total = 1-pcauchy(p_total)
  return (p_total)
}

 
#' @export
celltype_level = function(p,w = rep(1,ncol(p))){
  p = tan((0.5-p)*pi)
  p_total = apply(p,1,function(x){weighted.mean(x,w,na.rm = T)})
  p_total = 1-pcauchy(p_total)
  return (p_total)
}

#' @export
contrast_post = function(betas_all,V_cov_all,index,niche){
  #initialize p value array 
  ngene = dim(betas_all)[3]
  n_type = dim(betas_all)[1]
  #pgt is index type by niche type by gene 
  p = matrix(NA,ngene,1)
  rownames(p) = dimnames(betas_all)[[3]]
  #print('start')
  for(j in c(1:ngene)){
    if(j%%5000 == 0){
      paste0('gene #',j,' out of ', ngene)
    }
    V_cov = V_cov_all[,,j]
    betas = betas_all[,,j]
    #do if  gene is rejected and gene-type has at least 1 rejection
    tryCatch({
      #get sd of contrast
      index_1 = (index-1)*n_type + niche[1]
      index_2 = (index-1)*n_type + niche[2]
      var_contrast = diag(V_cov)[index_1]+diag(V_cov)[index_2]-2*V_cov[index_1,index_2]
      
      #test statistic
      T_stat = betas[index,niche[1]] - betas[index,niche[2]]
      T_stat = T_stat/(sqrt(var_contrast))
      #print(T_stat)
      #print('getting pvalue')
      p_stat = 1-pnorm(T_stat)
      p[j,] = p_stat} #get pval
      , error = function(e) { 
        print(paste0("error",j))
        skip_to_next <<- TRUE}) 
  }
  return(p)
}



check_colloc = function(object,index,niche){
  #initialize list
  colloc = rep(NA, length = length(object@sigma))
  names(colloc) = object@sigma
  #get indices where index is present
  index_ind = which(object@num_cells[,index]>0)
  counter = 1
  for(sig in object@sigma){
    #get indices such that niche is in effective niche
    niche_ind = which(object@effective_niche[[counter]][,niche] > min(object@effective_niche[[counter]][,niche]))
    #get length of overlap
    colloc[counter] = length(niche_ind %in% index_ind)
    counter = counter + 1
  }
  #return lengths of overlap
  return(colloc)
}