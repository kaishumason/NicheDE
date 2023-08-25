
#' @export
nb_lik = function(x,mu,disp){
  #returns negative log likelihood: Var = mu + mu^2/size
  return(-sum(dnbinom(x=x, size = disp, mu = mu,log = TRUE)))
}


#' @export
ultosymmetric=function(m){
  m = m + t(m) - diag(diag(m))
  return (m)
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
gene_level_fisher = function(p,varcov,beta_cov = T){
  #need to take transpose because conversion to vector is columnwise
  p = as.vector(t(p))
  #remove Na's from p
  p = p[is.na(p)==F]
  if(length(p)>0 & is.null(dim(varcov))==F){
    if(beta_cov == T){
      #convert varcov to matrix
      varcov = as.matrix(varcov)
      #make varcov matrix symmetric
      varcov = as.matrix(ultosymmetric(varcov))
      #get diagonal elements
      A = 1/diag(varcov)
      #make into diagonal matrix
      A = diag(A)
      #get variance of ABeta
      varcov = A%*%varcov%*%t(A)
      #convert to correlation matrix
      varcov = cov2cor(varcov)
      #get 1-sided correlation structure for pvalues
      varcov = poolr::mvnconv(varcov, target = "m2lp",side  = 1, cov2cor = F)
    }
    #do Browns test to get combined pvalue
    p_total = poolr::fisher(p,side = 1,R = varcov,adjust = "generalized")$p
    return(p_total)
  }else{
    return(NA)
  }
}

#' @export
celltype_level = function(p,w = rep(1,ncol(p))){
  p = tan((0.5-p)*pi)
  p_total = apply(p,1,function(x){weighted.mean(x,w,na.rm = T)})
  p_total = 1-pcauchy(p_total)
  return (p_total)
}


#' @export
celltype_level_fisher = function(p,varcov){
  #need to take transpose since indices are columnwise when doing vector conversion
  p = t(p)
  #initialize pval vector
  p_total = rep(NA,nrow(p))
  #want to iterate over each row of pvalue matrix
  ind_counter = 0
  if(is.null(dim(varcov))==F){
    #make varcov matrix symmetric
    varcov = as.matrix(ultosymmetric(as.matrix(varcov)))
    #get diagonal elements
    A = 1/diag(varcov)
    #make into diagonal matrix
    A = diag(A)
    #get variance of ABeta
    varcov = A%*%varcov%*%t(A)
    #convert to correlation matrix
    varcov = cov2cor(varcov)
    #get 1-sided correlation structure for pvalues
    varcov = poolr::mvnconv(varcov, target = "m2lp",side  = 1, cov2cor = F)
  }
  for(iter in c(1:ncol(p))){
    #print(iter)
    #get jth column. This corrsponds to pvalue for CT j
    p_CT = p[,iter]
    #make p_CT only those that are non NA
    p_CT = p_CT[is.na(p_CT)==F]
    #get which indices we need in the varcov matrix
    add = length(p_CT)
    #calculate adjusted pvalue if not all NA
    if(add > 1){
      #get range of indices
      range = c((ind_counter+1):(ind_counter + add))
      #get varcov submatrix
      Sigma = varcov[range,range]
      #calculate brown ptotal
      p_total[iter] = gene_level_fisher(p_CT,Sigma,beta_cov = F)
      #update ind_counter
      ind_counter = ind_counter + add
    }else if (add == 1){
      p_total[iter] = p_CT
    }
  }
  ind_counter = 0
  return (p_total)
}




#' @export
contrast_post = function(betas_all,V_cov_all,nulls_all,index,niche){
  #initialize p value array
  ngene = dim(betas_all)[3]
  n_type = dim(betas_all)[1]
  #pgt is index type by niche type by gene
  p = matrix(NA,ngene,1)
  rownames(p) = dimnames(betas_all)[[3]]
  #print('start')
  for(j in c(1:ngene)){
    if(j%%5000 == 0){
      print(paste0('gene #',j,' out of ', ngene))
    }
    #read in null values (remove those interactions)
    null = nulls_all[[j]]
    #reading covariance matrix
    V = as.matrix(V_cov_all[[j]])

    if(is.null(V)){
      V = matrix(NA,n_type^2,n_type^2)
    }
    #make matrix with nulls added
    V_cov = matrix(NA,n_type^2,n_type^2)
    if(length(null)==0){
      V_cov = V
    }else{
      V_cov[-null,-null] = as.matrix(V)}
    betas = betas_all[,,j]


    #do if  gene is rejected and gene-type has at least 1 rejection
    tryCatch({
      #get sd of contrast
      index_1 = (index-1)*n_type + niche[1]
      index_2 = (index-1)*n_type + niche[2]
      #get minimum of the two indices
      id1 = min(index_1,index_2)
      id2 = max(index_1,index_2)
      #compute contrast score
      var_contrast = diag(V_cov)[index_1]+diag(V_cov)[index_2]-2*V_cov[id1,id2]
      #print(diag(V_cov)[index_1])
      #print(diag(V_cov)[index_2])
      #print(V_cov[id1,id2])
      #test statistic
      T_stat = betas[index,niche[1]] - betas[index,niche[2]]
      T_stat = T_stat/(sqrt(var_contrast))
      #print(T_stat)
      #print('getting pvalue')
      p_stat = 1-pnorm(T_stat)
      p[j,] = p_stat} #get pval
      , error = function(e) {
        message(e)
        skip_to_next <<- TRUE})
  }
  return(p)
}

#' @export
contrast_post_new = function(betas_all,V_cov_all,nulls_all,index,niche){
  #initialize p value array
  ngene = length(betas_all)
  if(is.null(dim(V_cov_all[[1]]))){
    n_type = sqrt(length(nulls_all[[1]]))
  }else{
    n_type = dim(betas_all[[1]])[1]
  }
  #pgt is index type by niche type by gene
  p = matrix(NA,ngene,1)
  rownames(p) = dimnames(betas_all)[[3]]
  #print('start')
  for(j in c(1:ngene)){
    if(j%%5000 == 0){
      print(paste0('gene #',j,' out of ', ngene))
    }
    #read in null values (remove those interactions)
    null = nulls_all[[j]]
    #reading covariance matrix
    V = as.matrix(V_cov_all[[j]])

    if(is.null(V)){
      V = matrix(NA,n_type^2,n_type^2)
    }
    #make matrix with nulls added
    V_cov = matrix(NA,n_type^2,n_type^2)
    if(length(null)==0){
      V_cov = V
    }else{
      V_cov[-null,-null] = as.matrix(V)}
    betas = betas_all[[j]]

    if(length(null)!= n_type^2){
      tryCatch({
        #get sd of contrast
        index_1 = (index-1)*n_type + niche[1]
        index_2 = (index-1)*n_type + niche[2]
        #get minimum of the two indices
        id1 = min(index_1,index_2)
        id2 = max(index_1,index_2)
        #compute contrast score
        var_contrast = diag(V_cov)[index_1]+diag(V_cov)[index_2]-2*V_cov[id1,id2]
        #test statistic
        T_stat = betas[index,niche[1]] - betas[index,niche[2]]
        T_stat = T_stat/(sqrt(var_contrast))
        #print(T_stat)
        #print('getting pvalue')
        p_stat = 1-pnorm(T_stat)
        p[j,] = p_stat} #get pval
        , error = function(e) {
          message(e)
          skip_to_next <<- TRUE})
      }
    }
    #do if  gene is rejected and gene-type has at least 1 rejection
  return(p)
}





#' @export
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
    colloc[counter] = sum(niche_ind %in% index_ind == T)
    counter = counter + 1
  }
  #return lengths of overlap
  return(colloc)
}

#need to edit to get postive and negative

#' @export
get_niche_DE_pval = function(object,pos = T){

  #get lok likelihood matrix
  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  #get T_statistic list
  T_stat_list = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_stat_list[[j]] = object@niche_DE[[j]]$T_stat
  }

  if(length(T_stat_list)==1){
    T_stat_list[[2]] = T_stat_list[[1]]
    log_liks = cbind(as.vector(log_liks),as.vector(log_liks))
  }
  #T_stat_list is a list of T_stats from niche-DE
  #log_liks is a matrix with #gene rows and #kernels columns
  #alpha is the desired FDR (we test at level alpha/2 for postive values)
  #convert T_stat to pvalue
  if(pos ==T){
    T_stat_list_sign = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'positive')})
  }
  if(pos==F){
    T_stat_list_sign = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'negative')})
  }
  T_stat_list = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'two.sided')})
  #apply cauchy rule to get gene pvalues and make into a dataframe
  gene_p = lapply(T_stat_list,function(x){apply(x,3,function(y){gene_level(y)})})
  n <- length(gene_p[[1]])
  gene_p = as.matrix(structure(gene_p, row.names = c(NA, -n), class = "data.frame"))
  #get weights for cauchy rule
  suppressWarnings({W = t(apply(log_liks,1,function(x){exp(x-max(x,na.rm = T))}))})
  W[W < 0.1] = 0
  #W = apply(W,1,function(x){x/sum(x)})
  W[is.infinite(W)] = 10e160
  #bind pvalues and weights
  gene_p = cbind(gene_p,W)
  #apply cauchy rule
  gene_p = apply(gene_p,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])})
  #apply BH
  gene_phoch = p.adjust(gene_p,method = "BH")


  #CT level
  CT_p = lapply(T_stat_list,function(x){t(apply(x,3,function(y){celltype_level(y)}))})
  #merge cell type pvalue matrices into an array
  CT_merge = CT_p[[1]]
  if(length(CT_p)>1){
    for(sample in c(2:length(CT_p))){
      CT_merge <- abind::abind(CT_merge,CT_p[[sample]],along=3)
    }
  }
  CT_cauchy = CT_merge[,,1]
  #apply cauchy rule over these
  for(i in c(1:dim(CT_merge)[1])){
    for(j in c(1:dim(CT_merge)[2])){
      CT_cauchy[i,j] = gene_level(p = CT_merge[i,j,],w = W[i,])
    }
  }
  #remove genes that didn't reject before
  CT_hoch = CT_cauchy
  for(j in c(1:nrow(CT_hoch))){
    CT_hoch[j,] = p.adjust(CT_cauchy[j,],method = "BH")
  }

  #interaction Level
  pgt_merge = 0*T_stat_list_sign[[1]]
  valid = 0*T_stat_list_sign[[1]]
  #merge pvalues over all experiemnts
  for(sample in c(1:length(T_stat_list_sign))){
    pgt = tan((0.5-T_stat_list_sign[[sample]])*pi)
    for(j in c(1:dim(pgt)[3])){
      pgt[,,j] = pgt[,,j]*W[j,sample]
      valid[,,j] = valid[,,j] + W[j,sample]*(1-is.na(pgt[,,j]))
    }
    pgt[is.na(pgt)] = 0
    pgt_merge = pgt_merge + pgt
  }


  pgt_merge = pgt_merge/valid
  pgt_cauchy = 1-pcauchy(pgt_merge)

  #do BH on each gene interaction level
  pgt_hoch = pgt_cauchy
  for(j in c(1:dim(pgt_cauchy)[3])){
    for(k in c(1:dim(pgt_cauchy)[1])){
      pgt_hoch[k,,j] = p.adjust(pgt_cauchy[k,,j],method = "BH")
    }
  }


  if(pos==T){
    object@niche_DE_pval_pos = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  if(pos==F){
    object@niche_DE_pval_neg = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  return(object)
}

#' @export
get_niche_DE_pval_fisher = function(object,pos = T){

  #get log likelihood matrix
  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  #get T_statistic list
  T_stat_list = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_stat_list[[j]] = object@niche_DE[[j]]$T_stat
  }

  if(length(T_stat_list)==1){
    T_stat_list[[2]] = T_stat_list[[1]]
    log_liks = cbind(as.vector(log_liks),as.vector(log_liks))
  }
  #T_stat_list is a list of T_stats from niche-DE
  #log_liks is a matrix with #gene rows and #kernels columns
  #alpha is the desired FDR (we test at level alpha/2 for postive values)
  #convert T_stat to pvalue
  if(pos ==T){
    T_stat_list_sign = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'positive')})
    T_stat_list = T_stat_list_sign
  }
  if(pos==F){
    T_stat_list_sign = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'negative')})
    T_stat_list = T_stat_list_sign
  }
  #T_stat_list = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'two.sided')})
  #apply cauchy rule to get gene pvalues and make into a dataframe
  ngene = length(object@gene_names)
  #initialize pvalue matrix
  gene_p = matrix(NA,ngene,length(T_stat_list))
  #iterate over sigmas
  print("Computing Gene Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    #print(k)
    #print("making T_stats")
    if(k > length(object@sigma)){
      gene_p[,k] = gene_p[,1]
    }else{
      #convert array to a list
      T_stats = lapply(seq(dim(T_stat_list[[k]])[3]), function(x) T_stat_list[[k]][ , , x])
      #get varcov matrix
      varcovs = object@niche_DE[[k]]$var_cov
      #iterate over genes
      #calculate brown pvalue
      gene_p[,k] = mapply(gene_level_fisher,T_stats,varcovs,SIMPLIFY = T)
      #gene_p[,k] = gene_level_fisher(T_stat_list[[k]][,,j],varcov[[j]])
    }
  }
  print("Combining Gene Level Pvalues Across Kernel Bandwidths")
  #get weights for cauchy rule
  suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-max(x,na.rm = T))}))})
  W[W < 0.1] = 0
  #bind pvalues and weights
  gene_p = cbind(gene_p,W)
  #apply cauchy rule
  gene_p = apply(gene_p,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])})
  #apply BH
  gene_phoch = p.adjust(gene_p,method = "BH")


  #CT level
  #initialize cell type pvalue matrix
  CT_p = vector(mode='list', length=length(T_stat_list))
  ngene = length(object@gene_names)
  #iterate over sigmas
  print("Computing Cell Type Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    if(k > length(object@sigma)){
      CT_p[[k]] = CT_p[[1]]
    }else{
      #print("making T_stats")
      #convert array to a list
      T_stats = lapply(seq(dim(T_stat_list[[k]])[3]), function(x) T_stat_list[[k]][ , , x])
      #get varcov matrix
      varcovs = object@niche_DE[[k]]$var_cov
      #initialize pvalue matrix
      CT_p_sigma = matrix(NA,ngene,length(object@cell_types))
      #iterate over genes
      #calculate brown pvalue
      CT_pvals = mapply(celltype_level_fisher,T_stats,varcovs,SIMPLIFY = T)
      #rbind values
      for(j in c(1:ngene)){
        CT_p_sigma[j,] = CT_pvals[,j]
      }
    CT_p[[k]] = CT_p_sigma
  }
}
  #merge cell type pvalue matrices into an array
  CT_merge = CT_p[[1]]
  if(length(CT_p)>1){
    for(sample in c(2:length(CT_p))){
      CT_merge <- abind::abind(CT_merge,CT_p[[sample]],along=3)
    }
  }
  CT_cauchy = CT_merge[,,1]
  #apply cauchy rule over these
  print("Combining Cell Type  Level Pvalues Across Kernel Bandwidths")
  for(i in c(1:dim(CT_merge)[1])){
    for(j in c(1:dim(CT_merge)[2])){
      CT_cauchy[i,j] = gene_level(p = CT_merge[i,j,],w = W[i,])
    }
  }
  #remove genes that didn't reject before
  CT_hoch = CT_cauchy
  for(j in c(1:nrow(CT_hoch))){
    CT_hoch[j,] = p.adjust(CT_cauchy[j,],method = "BH")
  }

  print("Computing and Combining interaction Level Pvalues Across Kernel bandwidths")
  #interaction Level
  NG = dim(CT_hoch)[1]
  NCT = dim(CT_hoch)[2]
  pgt_merge = array(0,dim = c(NCT,NCT,NG))
  valid = array(0,dim = c(NCT,NCT,NG))
  #merge pvalues over all experiemnts
  for(sample in c(1:length(T_stat_list_sign))){
    pgt = tan((0.5-T_stat_list_sign[[sample]])*pi)
    for(j in c(1:dim(pgt)[3])){
      pgt[,,j] = pgt[,,j]*W[j,sample]
      valid[,,j] = valid[,,j] + W[j,sample]*(1-is.na(pgt[,,j]))
    }
    pgt[is.na(pgt)] = 0
    pgt_merge = pgt_merge + pgt
  }
  #Make indices with no valid NA
  pgt_merge[valid == 0] = NA
  #normalize
  pgt_merge = pgt_merge/valid
  pgt_cauchy = 1-pcauchy(pgt_merge)

  #do BH on each gene interaction level
  pgt_hoch = pgt_cauchy
  for(j in c(1:dim(pgt_cauchy)[3])){
    for(k in c(1:dim(pgt_cauchy)[1])){
      pgt_hoch[k,,j] = p.adjust(pgt_cauchy[k,,j],method = "BH")
    }
  }
  #make sure column and rownames are set
  rownames(CT_hoch) = object@gene_names
  colnames(CT_hoch) = object@cell_types


  if(pos==T){
    object@niche_DE_pval_pos = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  if(pos==F){
    object@niche_DE_pval_neg = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  return(object)
}

#' @export
get_niche_DE_pval_fisher_new = function(object,pos = T){
  #get log likelihood matrix
  log_liks = unlist(lapply(object@niche_DE[[1]], function(result) result$log_likelihood))
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,unlist(lapply(object@niche_DE[[j]], function(result) result$log_likelihood)))
    }
  }
  #get T_statistic list
  T_stat_list = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_stat_list[[j]] = lapply(object@niche_DE[[j]], function(result) result$T_stat)
  }

  if(length(T_stat_list)==1){
    T_stat_list[[2]] = T_stat_list[[1]]
    log_liks = cbind(as.vector(log_liks),as.vector(log_liks))
  }
  #T_stat_list is a list of T_stats from niche-DE
  #log_liks is a matrix with #gene rows and #kernels columns
  #alpha is the desired FDR (we test at level alpha/2 for postive values)
  #convert T_stat to pvalue
  if(pos ==T){
    T_stat_list_sign = lapply(T_stat_list,function(x){lapply(x,function(x){T_to_p(x,alternative = 'positive')})})
    T_stat_list = T_stat_list_sign
  }
  if(pos==F){
    T_stat_list_sign = lapply(T_stat_list,function(x){lapply(x,function(x){T_to_p(x,alternative = 'negative')})})
    T_stat_list = T_stat_list_sign
  }
  #T_stat_list = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'two.sided')})
  #apply cauchy rule to get gene pvalues and make into a dataframe
  ngene = length(object@gene_names)
  #initialize pvalue matrix
  gene_p = matrix(NA,ngene,length(T_stat_list))
  #iterate over sigmas
  print("Computing Gene Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    #print(k)
    #print("making T_stats")
    if(k > length(object@sigma)){
      gene_p[,k] = gene_p[,1]
    }else{
      #convert array to a list
      #T_stats = lapply(seq(dim(T_stat_list[[k]])[3]), function(x) T_stat_list[[k]][ , , x])
      T_stats = T_stat_list[[k]]
      #print("making pvalues")
      #get varcov list
      varcovs = lapply(object@niche_DE[[k]], function(result) result$Varcov)
      #iterate over genes
      #calculate brown pvalue
      gene_p[,k] = mapply(gene_level_fisher,T_stats,varcovs,SIMPLIFY = T)
    }
  }
  print("Combining Gene Level Pvalues Across Kernel Bandwidths")
  #get weights for cauchy rule
  suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-max(x,na.rm = T))}))})
  W[W < 0.1] = 0
  #bind pvalues and weights
  gene_p = cbind(gene_p,W)
  #apply cauchy rule
  gene_p = apply(gene_p,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])})
  #apply BH
  gene_phoch = p.adjust(gene_p,method = "BH")


  #CT level
  #initialize cell type pvalue matrix
  CT_p = vector(mode='list', length=length(T_stat_list))
  ngene = length(object@gene_names)
  #iterate over sigmas
  print("Computing Cell Type Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    if(k > length(object@sigma)){
      CT_p[[k]] = CT_p[[1]]
    }else{
      #print("making T_stats")
      T_stats = T_stat_list[[k]]
      #get varcov list
      varcovs = lapply(object@niche_DE[[k]], function(result) result$Varcov)
      #initialize pvalue matrix
      CT_p_sigma = matrix(NA,ngene,length(object@cell_types))
      #iterate over genes
      #calculate brown pvalue
      CT_pvals = mapply(celltype_level_fisher,T_stats,varcovs,SIMPLIFY = T)
      #rbind values
      for(j in c(1:ngene)){
        CT_p_sigma[j,] = CT_pvals[[j]]
      }
      CT_p[[k]] = CT_p_sigma
    }
  }
  #merge cell type pvalue matrices into an array
  CT_merge = CT_p[[1]]
  if(length(CT_p)>1){
    for(sample in c(2:length(CT_p))){
      CT_merge <- abind::abind(CT_merge,CT_p[[sample]],along=3)
    }
  }
  CT_cauchy = CT_merge[,,1]
  #apply cauchy rule over these
  print("Combining Cell Type  Level Pvalues Across Kernel Bandwidths")
  for(i in c(1:dim(CT_merge)[1])){
    for(j in c(1:dim(CT_merge)[2])){
      CT_cauchy[i,j] = gene_level(p = CT_merge[i,j,],w = W[i,])
    }
  }
  #remove genes that didn't reject before
  CT_hoch = CT_cauchy
  for(j in c(1:nrow(CT_hoch))){
    CT_hoch[j,] = p.adjust(CT_cauchy[j,],method = "BH")
  }

  print("Computing and Combining interaction Level Pvalues Across Kernel bandwidths")
  #interaction Level
  NG = dim(CT_hoch)[1]
  NCT = dim(CT_hoch)[2]
  pgt_merge = array(0,dim = c(NCT,NCT,NG))
  valid = array(0,dim = c(NCT,NCT,NG))
  #merge pvalues over all experiemnts
  for(sample in c(1:length(T_stat_list_sign))){
    pgt_temp = array(0,dim = c(NCT,NCT,NG))
    for(j in c(1:dim(pgt_temp)[3])){
      x = T_stat_list_sign[[sample]][[j]]
      pgt = tan((0.5-x)*pi)
      pgt_temp[,,j] = pgt*W[j,sample]
      valid[,,j] = valid[,,j] + W[j,sample]*(1-is.na(pgt))
    }
    pgt_temp[is.na(pgt_temp)] = 0
    pgt_merge = pgt_merge + pgt_temp
  }
  #Make indices with no valid NA
  pgt_merge[valid == 0] = NA

  pgt_merge = pgt_merge/valid
  pgt_cauchy = 1-pcauchy(pgt_merge)


  #do BH on each gene interaction level
  pgt_hoch = pgt_cauchy
  for(j in c(1:dim(pgt_cauchy)[3])){
    for(k in c(1:dim(pgt_cauchy)[1])){
      pgt_hoch[k,,j] = p.adjust(pgt_cauchy[k,,j],method = "BH")
    }
  }
  #make sure column and rownames are set
  rownames(CT_hoch) = object@gene_names
  colnames(CT_hoch) = object@cell_types


  if(pos==T){
    object@niche_DE_pval_pos = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  if(pos==F){
    object@niche_DE_pval_neg = list(gene_level = gene_phoch,
                                    cell_type_level = CT_hoch,
                                    interaction_level = pgt_hoch)
  }

  return(object)
}


#' @export
get_niche_DE_pval_raw = function(object,pos = T){
  #get log likelihood matrix
  log_liks = unlist(lapply(object@niche_DE[[1]], function(result) result$log_likelihood))
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,unlist(lapply(object@niche_DE[[j]], function(result) result$log_likelihood)))
    }
  }
  #get T_statistic list
  T_stat_list = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_stat_list[[j]] = lapply(object@niche_DE[[j]], function(result) result$T_stat)
  }

  if(length(T_stat_list)==1){
    T_stat_list[[2]] = T_stat_list[[1]]
    log_liks = cbind(as.vector(log_liks),as.vector(log_liks))
  }
  #T_stat_list is a list of T_stats from niche-DE
  #log_liks is a matrix with #gene rows and #kernels columns
  #alpha is the desired FDR (we test at level alpha/2 for postive values)
  #convert T_stat to pvalue
  if(pos ==T){
    T_stat_list_sign = lapply(T_stat_list,function(x){lapply(x,function(x){T_to_p(x,alternative = 'positive')})})
    T_stat_list = T_stat_list_sign
  }
  if(pos==F){
    T_stat_list_sign = lapply(T_stat_list,function(x){lapply(x,function(x){T_to_p(x,alternative = 'negative')})})
    T_stat_list = T_stat_list_sign
  }
  #T_stat_list = lapply(T_stat_list,function(x){T_to_p(x,alternative = 'two.sided')})
  #apply cauchy rule to get gene pvalues and make into a dataframe
  ngene = length(object@gene_names)
  #initialize pvalue matrix
  gene_p = matrix(NA,ngene,length(T_stat_list))
  #iterate over sigmas
  print("Computing Gene Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    #print(k)
    #print("making T_stats")
    if(k > length(object@sigma)){
      gene_p[,k] = gene_p[,1]
    }else{
      #convert array to a list
      #T_stats = lapply(seq(dim(T_stat_list[[k]])[3]), function(x) T_stat_list[[k]][ , , x])
      T_stats = T_stat_list[[k]]
      #print("making pvalues")
      #get varcov list
      varcovs = lapply(object@niche_DE[[k]], function(result) result$Varcov)
      #iterate over genes
      #calculate brown pvalue
      gene_p[,k] = mapply(gene_level_fisher,T_stats,varcovs,SIMPLIFY = T)
    }
  }
  print("Combining Gene Level Pvalues Across Kernel Bandwidths")
  #get weights for cauchy rule
  suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-max(x,na.rm = T))}))})
  W[W < 0.1] = 0
  #bind pvalues and weights
  gene_p = cbind(gene_p,W)
  #apply cauchy rule
  gene_p = apply(gene_p,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])})
  #apply BH
  gene_phoch = p.adjust(gene_p,method = "BH")


  #CT level
  #initialize cell type pvalue matrix
  CT_p = vector(mode='list', length=length(T_stat_list))
  ngene = length(object@gene_names)
  #iterate over sigmas
  print("Computing Cell Type Level Pvalues")
  for(k in c(1:length(T_stat_list))){
    if(k > length(object@sigma)){
      CT_p[[k]] = CT_p[[1]]
    }else{
      #print("making T_stats")
      T_stats = T_stat_list[[k]]
      #get varcov list
      varcovs = lapply(object@niche_DE[[k]], function(result) result$Varcov)
      #initialize pvalue matrix
      CT_p_sigma = matrix(NA,ngene,length(object@cell_types))
      #iterate over genes
      #calculate brown pvalue
      CT_pvals = mapply(celltype_level_fisher,T_stats,varcovs,SIMPLIFY = T)
      #rbind values
      for(j in c(1:ngene)){
        CT_p_sigma[j,] = CT_pvals[[j]]
      }
      CT_p[[k]] = CT_p_sigma
    }
  }
  #merge cell type pvalue matrices into an array
  CT_merge = CT_p[[1]]
  if(length(CT_p)>1){
    for(sample in c(2:length(CT_p))){
      CT_merge <- abind::abind(CT_merge,CT_p[[sample]],along=3)
    }
  }
  CT_cauchy = CT_merge[,,1]
  #apply cauchy rule over these
  print("Combining Cell Type  Level Pvalues Across Kernel Bandwidths")
  for(i in c(1:dim(CT_merge)[1])){
    for(j in c(1:dim(CT_merge)[2])){
      CT_cauchy[i,j] = gene_level(p = CT_merge[i,j,],w = W[i,])
    }
  }
  #remove genes that didn't reject before
  CT_hoch = CT_cauchy
  for(j in c(1:nrow(CT_hoch))){
    CT_hoch[j,] = p.adjust(CT_cauchy[j,],method = "BH")
  }

  print("Computing and Combining interaction Level Pvalues Across Kernel bandwidths")
  #interaction Level
  NG = dim(CT_hoch)[1]
  NCT = dim(CT_hoch)[2]
  pgt_merge = array(0,dim = c(NCT,NCT,NG))
  valid = array(0,dim = c(NCT,NCT,NG))
  #merge pvalues over all experiemnts
  for(sample in c(1:length(T_stat_list_sign))){
    pgt_temp = array(0,dim = c(NCT,NCT,NG))
    for(j in c(1:dim(pgt_temp)[3])){
      x = T_stat_list_sign[[sample]][[j]]
      pgt = tan((0.5-x)*pi)
      pgt_temp[,,j] = pgt*W[j,sample]
      valid[,,j] = valid[,,j] + W[j,sample]*(1-is.na(pgt))
    }
    pgt_temp[is.na(pgt_temp)] = 0
    pgt_merge = pgt_merge + pgt_temp
  }
  #Make indices with no valid NA
  pgt_merge[valid == 0] = NA

  pgt_merge = pgt_merge/valid
  pgt_cauchy = 1-pcauchy(pgt_merge)


  #do BH on each gene interaction level
  pgt_hoch = pgt_cauchy
  for(j in c(1:dim(pgt_cauchy)[3])){
    for(k in c(1:dim(pgt_cauchy)[1])){
      pgt_hoch[k,,j] = p.adjust(pgt_cauchy[k,,j],method = "BH")
    }
  }
  #make sure column and rownames are set
  rownames(CT_hoch) = object@gene_names
  colnames(CT_hoch) = object@cell_types


  if(pos==T){
    object@niche_DE_pval_pos = list(gene_level = gene_p,
                                    cell_type_level = CT_cauchy,
                                    interaction_level = pgt_cauchy)
  }

  if(pos==F){
    object@niche_DE_pval_neg = list(gene_level = gene_p,
                                    cell_type_level = CT_cauchy,
                                    interaction_level = pgt_cauchy)
  }

  return(object)
}










#' export
niche_DE_core_correct = function(betas_all,V_cov_all,index,niche){
  return("hello world")
}
#
# niche_DE_core_new = function(counts,constants,iter){
#   #get reference expression vector
#   ref_expression = constants$ref_expr[,iter]
#   #get expected expression of gene j
#   EEJ = as.vector(constants$num_cells%*%ref_expression)
#   #get gene number
#   if(iter%%1000 == 0){
#     print(paste0("Gene ",iter," out of ",constants$gene_total))
#   }
#
#   valid = 0
#   liks_val = NA
#   n_type = ncol(constants$effective_niche)
#   #check if we need to do niche-DE
#   if((sum(counts)<constants$C) | (mean(ref_expression<constants$CT_filter)!=1)==F){
#     null = c(1:n_type^2)
#     liks_val = NA
#   }
#   if((sum(counts)>constants$C)&(mean(ref_expression<constants$CT_filter)!=1)){
#     #get pstg matrix
#     #print(j)
#     #t1 = Sys.time()
#     pstg = constants$num_cells%*%as.matrix(diag(ref_expression))/EEJ
#     pstg[,ref_expression<constants$CT_filter] = 0
#     pstg[pstg<0.05]=0
#     #get X
#     #print(1)
#     X = matrix(NA,nrow(pstg),n_type^2)
#     #get counter
#     for(k in c(1:nrow(pstg))){
#       #get spatial covariates
#       X[k,] = as.vector(round(constants$effective_niche[k,],2)%*%t(as.matrix(pstg[k,])))
#     }
#     rm("pstg")
#     #gc()
#     #get index, niche pairs that are non existent
#     null = which(apply(X,2,function(x){sum(x>0,na.rm = T)})<constants$M)
#     #print(length(null))
#     X_partial = X
#     rest = c(1:ncol(X))
#     if(length(null)>0){
#       X_partial = X[,-null]
#       rest = rest[-null]
#     }
#     #get how many important variables there are before batchID is added
#     nvar = ncol(X_partial)
#     rm("X")
#     #add batch
#     if(constants$batch == T &length(unique(constants$batchID)) > 1){
#       batchvar = as.factor(constants$batchID)
#       dummy_col = as.matrix(fastDummies::dummy_cols(batchvar,remove_first_dummy = T, remove_selected_columns = T))
#       #append dummy variables and intercept
#       X_partial = cbind(X_partial,dummy_col)
#     }
#     #gc()
#     #continue if at least one index,niche pair is viable
#     if(length(null)!=n_type^2 & constants$Int == T){
#       tryCatch({
#         #if expected expression for a spot is 0, remove it
#         bad_ind  = which(EEJ==0)
#         #run neg binom regression
#         if(length(bad_ind)>0){
#           full_glm =suppressWarnings({glm(counts[-bad_ind]~X_partial[-bad_ind,] + offset(log(EEJ[-bad_ind])), family = "poisson")}) #do full glm
#         }else{
#           full_glm = suppressWarnings({glm(counts~X_partial + offset(log(EEJ)), family = "poisson")}) #do full glm
#         }
#         mu_hat = as.numeric(exp(predict(full_glm)))#get mean
#         #get dicpersion parameter
#         A = optimize(nb_lik,x = counts,mu = mu_hat, lower = 0.05, upper = 100) #get overdispersion parameter
#         #save dispersion parameter
#         disp = A$minimum
#         #save likelihood
#         liks_val = -A$objective
#         #calculate W matrix for distribution of beta hat
#         W =as.vector(mu_hat/(1 + mu_hat/disp))#get W matrix
#         #W =as.vector(mu_hat)#get W matrix
#         rm("mu_hat")
#         #perform cholesky decomp for finding inverse of X^TWX
#         if(length(bad_ind)>0){
#           X_partial = X_partial[-bad_ind,]
#           #X_partial = Matrix::Matrix(X_partial, sparse=TRUE)
#         }
#         #add intercept
#         X_partial = cbind(rep(1,nrow(X_partial)),X_partial)
#         #make X_partial a matrix
#         X_partial = as.matrix(X_partial)
#         #make X_partial a sparse matrix
#         X_partial = Matrix::Matrix(X_partial, sparse=TRUE)
#         #get variance matrix. Variance is [t(X_partial*W)%*%X_partial]^(-1)
#         var_mat = Matrix::t(X_partial*W)%*%X_partial
#         #subset to get covariance of variables or interest
#         #var_mat = var_mat[1:nvar,1:nvar]
#         rm("X_partial")
#         #gc()
#         #if there are degenerate columns, remove them
#         new_null = c()
#
#         #account for batch ID variables (have n batch ID variables on right and intercept on left)
#         new_null = which(diag(as.matrix(var_mat))[-1]==0)
#         if(length(new_null)>0){
#           #plus one since intercept is removed
#           var_mat = var_mat[-(new_null+1),-(new_null+1)]
#           #no need for the plus one since this is in teh context of the nvar important variables
#           null = sort(c(null,rest[new_null[new_null <= nvar]]))
#           #add to new_null
#           new_null = new_null[new_null <= nvar]
#         }
#
#         if(length(null)!=n_type^2){
#           #get coefficients for important variables and intercept
#           coeff = full_glm$coefficients[1:(nvar+1)]
#           #cholesky decomposition
#           A = Matrix::chol(var_mat,LDL = FALSE,perm = FALSE)
#           #get covaraince matrix
#           V = Matrix::solve(A)%*%Matrix::t(Matrix::solve(A))
#           #get second column to nvar+1 column
#           V = V[2:(n_type^2-length(null)+1),2:(n_type^2-length(null)+1)]
#           #save V as an upper triangular matrix
#           V =  Matrix::Matrix(upper.tri(V,diag = T)*V, sparse=TRUE)
#           #remove large objects
#           rm("A")
#           rm("var_mat")
#           #get sd matrix
#           tau = sqrt(Matrix::diag(V))
#           #print("3")
#           V_ = matrix(NA,n_type,n_type)
#           if(length(null)==0){
#             V_ = matrix(tau,n_type,n_type)
#           }else{
#             V_[c(1:n_type^2)[-null]] = tau}
#           #get beta
#           beta = matrix(NA,n_type,n_type)
#           if(length(new_null)>0){
#             beta[c(1:n_type^2)[-null]] = coeff[-c(1,new_null+1)]
#           }
#
#           if(length(new_null)==0){
#             if(length(null)==0){
#               beta = matrix(coeff[-c(1)],n_type,n_type)
#             }else{
#               beta[c(1:n_type^2)[-null]] = coeff[-c(1)]}
#           }
#           #record test statitistic
#           T_ = beta/V_
#           valid = 1
#
#
#         }
#         #end of if statement
#       }, #get pval
#       error = function(e) {
#         #print(paste0("error,",j))
#         #return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
#         skip_to_next <<- TRUE})
#     }
#     if(length(null)!=n_type^2 & constants$Int == F){
#       tryCatch({
#         #if expected expression for a spot is 0, remove it
#         lm = lm((counts - EEJ) ~ X_partial)
#         rm("X_partial")
#         #gc()
#         sum_lm =  summary(lm)
#         #get log likelihood
#         liks_val = stats::logLik(lm)
#         #get vcov mat
#         V = (sum_lm$cov.unscaled)*(sum_lm$sigma^2)
#         #remove first observation
#         #V = V[-c(1),-c(1)]
#         #see if any bettas have 0 variance
#         new_null = which(diag(as.matrix(V))[-1]==0)
#         if(length(new_null)>0){
#           V = V[-(new_null+1),-(new_null+1)]
#           null = sort(c(null,rest[new_null[new_null <= nvar]]))
#           new_null = new_nul[new_null <= nvar]
#         }
#
#         #remove first observation
#         V = V[-c(1),-c(1)]
#         #remove batch variables
#         V[1:(n_type^2-length(null)+1),1:(n_type^2-length(null)+1)]
#         #save V as an upper triangular matrix
#         V =  Matrix::Matrix(upper.tri(V,diag = T)*V, sparse=TRUE)
#         #get beta coefficients
#         if(length(null)!=n_type^2){
#           #get coefficients
#           coeff = lm$coefficients[1:(nvar+1)]
#           #get number of true coefficents
#           num_coef = n_type^2 - length(null)
#           coeff_T = sum_lm$coefficients[c(1:(num_coef+1)),3]
#           #print('getting beta')
#           beta = matrix(NA,n_type,n_type)
#           T_ = matrix(NA,n_type,n_type)
#           if(length(new_null)>0){
#             beta[c(1:n_type^2)[-null]] = coeff[-c(1,new_null+1)]
#             T_[c(1:n_type^2)[-null]] = coeff_T[-c(1)]
#           }
#
#           if(length(new_null)==0){
#             if(length(null)==0){
#               beta = matrix(coeff[-c(1)],n_type,n_type)
#               T_ = matrix(coeff_T[-c(1)],n_type,n_type)
#             }else{
#               beta[c(1:n_type^2)[-null]] = coeff[-c(1)]
#               T_[c(1:n_type^2)[-null]] = coeff_T[-c(1)]
#             }
#           }
#           #record test statitistic
#           valid = 1
#
#         }
#         #end of if statement
#       }, #get pval
#       error = function(e) {
#         #message(e)
#         #print(paste0("error,"))
#         #return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
#         skip_to_next <<- TRUE})
#     }
#   }
#   #print(pryr::memused())
#   if(valid == 1){
#     A = list(T_ = Matrix::t(T_),betas = Matrix::t(beta), V = V,nulls = null, valid = valid,liks = liks_val)
#     rm(list=ls()[! ls() %in% c("A")])
#     #gc()
#     return (A)
#   }else{
#     A = list(T_ = NA,betas = NA, V = NA,nulls = c(1:n_type^2), valid = 0,liks = liks_val)
#     rm(list=ls()[! ls() %in% c("A")])
#     #gc()
#     return (A)
#   }
#
# }
#
