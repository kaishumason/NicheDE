
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

#' export
#contrast_post = function(betas_all,V_cov_all,index,niche){
  #initialize p value array
  #ngene = dim(betas_all)[3]
  #n_type = dim(betas_all)[1]
  #pgt is index type by niche type by gene
  #p = matrix(NA,ngene,1)
  #rownames(p) = dimnames(betas_all)[[3]]
  #print('start')
  #for(j in c(1:ngene)){
    #if(j%%5000 == 0){
      #print(paste0('gene #',j,' out of ', ngene))
    #}
    #V_cov = V_cov_all[[j]]
    #betas = betas_all[,,j]
    #do if  gene is rejected and gene-type has at least 1 rejection
    #tryCatch({
      #get sd of contrast
      #index_1 = (index-1)*n_type + niche[1]
      #index_2 = (index-1)*n_type + niche[2]
      #var_contrast = diag(V_cov)[index_1]+diag(V_cov)[index_2]-2*V_cov[index_1,index_2]

      #test statistic
      #T_stat = betas[index,niche[1]] - betas[index,niche[2]]
      #T_stat = T_stat/(sqrt(var_contrast))
      #print(T_stat)
      #print('getting pvalue')
      #p_stat = 1-pnorm(T_stat)
      #p[j,] = p_stat} #get pval
      #, error = function(e) {
      #  print(paste0("error",j))
      #  skip_to_next <<- TRUE})
  #}
 # return(p)
#}

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
    V = V_cov_all[[j]]

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
      var_contrast = diag(V_cov)[index_1]+diag(V_cov)[index_2]-2*V_cov[index_1,index_2]

      #test statistic
      T_stat = betas[index,niche[1]] - betas[index,niche[2]]
      T_stat = T_stat/(sqrt(var_contrast))
      #print(T_stat)
      #print('getting pvalue')
      p_stat = 1-pnorm(T_stat)
      p[j,] = p_stat} #get pval
      , error = function(e) {
        skip_to_next <<- TRUE})
  }
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
    colloc[counter] = length(niche_ind %in% index_ind)
    counter = counter + 1
  }
  #return lengths of overlap
  return(colloc)
}



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
  suppressWarnings({W = t(apply(log_liks,1,function(x){exp(x-min(x,na.rm = T))}))})
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
niche_DE_core = function(object,j,sig,CT_filter,C = 150,M = 10,gamma = 0.8,print = T){
  valid = 0
  n_type = ncol(object@num_cells)

  if(j%%1000 == 0 & print == T){
    print(paste0('kernel bandwidth:', sig,' (number ',counter,' out of ',length(object@sigma),' values), ', "Processing Gene #",j,
                 ' out of ',ncol(object@counts)))
  }
  if((sum(object@counts[,j])>C)&(mean(object@ref_expr[,j]<CT_filter)!=1)==F){
    null = c(1:n_type^2)
  }
  if((sum(object@counts[,j])>C)&(mean(object@ref_expr[,j]<CT_filter)!=1)){
    #get pstg matrix
    #print(j)
    #t1 = Sys.time()
    pstg = object@num_cells%*%as.matrix(diag(object@ref_expr[,j]))/object@null_expected_expression[,j]
    pstg[,object@ref_expr[,j]<CT_filter] = 0
    pstg[pstg<0.05]=0
    #get X
    #print(1)
    X = matrix(NA,nrow(pstg),n_type^2)
    for(k in c(1:nrow(pstg))){
      #get feature matrix by multiplying effective niche and pstg vector
      ps = as.matrix(pstg[k,])
      EN_j = round(object@effective_niche[[counter]][k,],2)
      cov_j = ps%*%t(EN_j)
      #make into a vector
      X[k,] = as.vector(t(cov_j))#important to take the transpose
    }
    #get index, niche pairs that are non existent
    null = which(apply(X,2,function(x){sum(x>0,na.rm = T)})<M)
    #print(length(null))
    X_partial = X
    rest = c(1:ncol(X))
    if(length(null)>0){
      X_partial = X[,-null]
      rest = rest[-null]
    }
    #continue if at least one index,niche pair is viable
    if(length(null)!=n_type^2){
      tryCatch({
        #if expected expression for a spot is 0, remove it
        bad_ind  = which(object@null_expected_expression[,j]==0)
        #print('Running GLM')
        #run neg binom regression
        #print()
        if(length(bad_ind)>0){
          full_glm =suppressWarnings({glm(object@counts[-bad_ind,j]~X_partial[-bad_ind,] + offset(log(object@null_expected_expression[-bad_ind,j])), family = "poisson")}) #do full glm
        }else{
          full_glm = suppressWarnings({glm(object@counts[,j]~X_partial + offset(log(object@null_expected_expression[,j])), family = "poisson")}) #do full glm
        }
        mu_hat = exp(predict(full_glm))#get mean
        #print(Sys.time()-t1)
        #get dicpersion parameter
        A = optimize(nb_lik,x = object@counts[,j],mu = mu_hat, lower = 0.05, upper = 100) #get overdispersion parameter
        #save dispersion parameter
        disp = A$minimum
        #save likelihood
        liks_val = -A$objective
        #calculate W matrix for distribution of beta hat
        W =as.vector(mu_hat/(1 + mu_hat/disp))#get W matrix
        #print(3)
        #perform cholesky decomp for finding inverse of X^TWX
        if(length(bad_ind)>0){
          X_partial = as((X_partial[-bad_ind,]),"sparseMatrix")
          #remove bad indices
        }else{
          X_partial = as((X_partial),"sparseMatrix")
        }
        #get variance matrix. Variance is [t(X_partial*W)%*%X_partial]^(-1)
        var_mat = Matrix::t(X_partial*W)%*%X_partial
        #if there are degenerate columns, remove them
        new_null = c()
        if(length(bad_ind)>0){
          new_null = which(diag(as.matrix(var_mat))==0)
          if(length(new_null)>0){
            var_mat = var_mat[-new_null,-new_null]
            null = sort(c(null,rest[new_null]))
          }
        }
        if(length(null)!=n_type^2){
          #cholesky decomposition
          A = Matrix::chol(var_mat,LDL = FALSE,perm = FALSE)
          #get covaraince matrix
          #print("1")
          V = Matrix::solve(A)%*%Matrix::t(Matrix::solve(A))
          #get standard devaition vector
          #print("2")
          #get sd matrix
          tau = sqrt(Matrix::diag(V))
          #print("3")
          V_ = matrix(NA,n_type,n_type)
          if(length(null)==0){
            V_ = matrix(tau,n_type,n_type)
          }else{
            V_[c(1:n_type^2)[-null]] = tau}
          #for full var_cov matrix.
          v_cov = matrix(NA,n_type^2,n_type^2)
          if(length(null)==0){
            v_cov = matrix(V,n_type^2,n_type^2)
          }else{
            v_cov[-null,-null] = as.matrix(V)}
          #print("4")
          #print('getting beta')
          beta = matrix(NA,n_type,n_type)

          if(length(new_null)>0){
            beta[c(1:n_type^2)[-null]] = full_glm$coefficients[-c(1,new_null+1)]
          }

          if(length(new_null)==0){
            if(length(null)==0){
              beta = matrix(full_glm$coefficients[-c(1)],n_type,n_type)
            }else{
              beta[c(1:n_type^2)[-null]] = full_glm$coefficients[-c(1)]}
          }
          #record test statitistic

          T_ = Matrix::t(beta/V_)
          valid = 1
          #
          # T_stat[,,j] = T_
          #
          # betas[,,j] = Matrix::t(beta)
          #
          # var_cov[[j]] = as.matrix(V)
          #
          # nulls[[j]] = null
          #
          # valid[j,counter] = 1

        }
        #end of if statement
      }, #get pval
      error = function(e) {
        #print(paste0("error,",j))
        return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
        skip_to_next <<- TRUE})
    }
  }
  if(valid == 1){
    return (list(T_ = T_,betas = Matrix::t(beta), V = as.matrix(V),nulls = null, valid = valid,liks = liks_val))
  }else{
    return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
  }

}








