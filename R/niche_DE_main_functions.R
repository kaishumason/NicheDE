#' Niche_DE
#'
#' This function performs niche-DE
#' @param object A niche-DE object
#' @param C Minimum total expression of a gene needed for the model to run. Default value is 150.
#' @param M Minimum number of spots containing the index cell type with the
#' niche cell type in its effective niche for (index,niche) niche patterns
#' to be investigated. Default value is 10.
#' @param gamma Percentile a gene needs to be with respect to expression in the
#'  index cell type in order for the model to investigate niche patterns for
#'  that gene in the index cell. Default value is 0.8 (80th percentile)
#' @param print Logical if function should print progress report (kernel, gene #)
#' @return A niche-DE object with niche-DE analysis performed
#' @export
niche_DE = function(object,C = 150,M = 10,gamma = 0.8,print = T){
  #starting Message
  print(paste0('Starting Niche-DE analysis with parameters C = ',C,', M = ',M,', gamma = ', gamma,'.'))
  #initialize list output
  object@niche_DE = vector(mode = "list", length = length(object@sigma))
  names(object@niche_DE) = object@sigma
  counter = 1
  valid = matrix(0,ncol(object@counts),length(object@sigma))
  #iterate over each sigma value
  for(sig in object@sigma){
    print(paste0('Performing Niche-DE analysis with kernel bandwidth:',sig,' (number ',counter,' out of ',length(object@sigma),' values)'))
    #get expression filter (gamma)
    CT_filter = apply(object@ref_expr,1,function(x){quantile(x,gamma)})
    #initialize p value array
    ngene = ncol(object@counts)
    n_type = ncol(object@num_cells)
    dimnames = list(A = colnames(object@num_cells),B  = colnames(object@num_cells), C = colnames(object@counts))
    #pgt is index type by niche type by gene
    T_stat = array(NA,c(n_type,n_type,ngene),dimnames = dimnames)
    #var_cov is too large for even moderately many cell types so will use list instead
    #var_cov = array(NA,c(n_type^2,n_type^2,ngene))
    var_cov = vector(mode='list', length=ngene)
    names(var_cov) = colnames(object@counts)
    nulls = vector(mode='list', length=ngene)
    names(nulls) = colnames(object@counts)
    #initalize betas
    betas = array(NA,c(n_type,n_type,ngene),dimnames = dimnames)
    liks = rep(NA,ngene)
    for(j in c(1:ngene)){
      if(j%%1000 == 0 & print == T){
        print(paste0('kernel bandwidth:', sig,' (number ',counter,' out of ',length(object@sigma),' values), ', "Processing Gene #",j,
                     ' out of ',ncol(object@counts)))
      }
      if((sum(object@counts[,j])>C)&(mean(object@ref_expr[,j]<CT_filter)!=1)==F){
        nulls[[j]] = c(1:n_type^2)
      }
      #do if  gene is rejected and gene-type has at least 1 rejection
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
            liks[j] = -A$objective
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

              T_stat[,,j] = T_

              betas[,,j] = Matrix::t(beta)

              var_cov[[j]] = as.matrix(V)

              nulls[[j]] = null

              valid[j,counter] = 1

            }
            #end of if statement
            }, #get pval
              error = function(e) {
              #print(paste0("error,",j))
              skip_to_next <<- TRUE})
          }
      }
    }
    #save object
    object@niche_DE[[counter]] = list(T_stat = T_stat,beta = betas,var_cov = var_cov,nulls = nulls,log_lik = liks)
    counter = counter + 1
  }
  #get column sums of counts matrix to see how many genes pass filtering
  A = rowSums(valid)
  #get number of genes that pass filtering
  num_pass = sum(A>=1,na.rm = T)
  print('Computing Niche-DE Pvalues')
  object = get_niche_DE_pval(object,pos = T)
  object = get_niche_DE_pval(object,pos = F)
  print(paste0('Niche-DE analysis complete. Number of Genes with niche-DE T-stat equal to ',num_pass))
  if(num_pass < 1000){
    warning('Less than 1000 genes pass. This could be due to insufficient read depth of data or size of C parameter. Consider changing choice of C parameter')
  }
  return(object)
}

#' Niche_DE
#'
#' This function performs niche-DE
#' @param object A niche-DE object
#' @param C Minimum total expression of a gene needed for the model to run. Default value is 150.
#' @param M Minimum number of spots containing the index cell type with the
#' niche cell type in its effective niche for (index,niche) niche patterns
#' to be investigated. Default value is 10.
#' @param gamma Percentile a gene needs to be with respect to expression in the
#'  index cell type in order for the model to investigate niche patterns for
#'  that gene in the index cell. Default value is 0.8 (80th percentile)
#' @param print Logical if function should print progress report (kernel, gene #)
#' @param cores Number of cores to use for paralellization
#' @return A niche-DE object with niche-DE analysis performed
#' @export
#' @importFrom foreach %dopar%
niche_DE_parallel = function(object,C = 150,M = 10,gamma = 0.8,print = T,cores = 4,Int = T){
  #use core niche-DE function and nb_lik function
  nb_lik = function(x,mu,disp){
    #returns negative log likelihood: Var = mu + mu^2/size
    return(-sum(dnbinom(x=x, size = disp, mu = mu,log = TRUE)))
  }

  niche_DE_core = function(object,j,sig,CT_filter,C = 150,M = 10,gamma = 0.8,Int = T){
    valid = 0
    liks_val = NA
    n_type = ncol(object@num_cells)
    #check if we need to do niche-DE
    if((sum(object@counts[,j])<C) | (mean(object@ref_expr[,j]<CT_filter)!=1)==F){
      null = c(1:n_type^2)
      liks_val = NA
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
      #get counter
      counter = which(object@sigma == sig)[1]
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
      if(length(null)!=n_type^2 & Int == T){
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
          mu_hat = as.numeric(exp(predict(full_glm)))#get mean
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
            X_partial = X_partial[-bad_ind,]
            #X_partial = Matrix(X_partial, sparse=TRUE)
            #X_partial = as((X_partial[-bad_ind,]),"sparseMatrix")
            #remove bad indices
          }else{
            #X_partial = as((X_partial),"sparseMatrix")
            #X_partial = Matrix(X_partial, sparse=TRUE)
            #Matrix(regMat, sparse=TRUE)
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
          #return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
          skip_to_next <<- TRUE})
      }
      if(length(null)!=n_type^2 & Int == F){
        tryCatch({
          #if expected expression for a spot is 0, remove it
          lm = lm((object@counts[,j] - object@null_expected_expression[,j]) ~ X_partial)
          sum_lm =  summary(lm)
          #get log likelihood
          liks_val = stats::logLik(lm)
          #get vcov mat
          var_mat = (sum_lm$cov.unscaled)*(sum_lm$sigma^2)
          #remove first observation
          var_mat = var_mat[-c(1),-c(1)]
          #see if any bettas have 0 variance
          new_null = c()
          new_null = which(diag(as.matrix(var_mat))==0)
          if(length(new_null)>0){
            var_mat = var_mat[-new_null,-new_null]
            null = sort(c(null,rest[new_null]))
          }
          #get beta coefficients
          if(length(null)!=n_type^2){

            #print('getting beta')
            beta = matrix(NA,n_type,n_type)
            T_ = matrix(NA,n_type,n_type)
            if(length(new_null)>0){
              beta[c(1:n_type^2)[-null]] = lm$coefficients[-c(1,new_null+1)]
              T_[c(1:n_type^2)[-null]] = sum_lm$coefficients[-c(1,new_null+1),3]
            }

            if(length(new_null)==0){
              if(length(null)==0){
                beta = matrix(lm$coefficients[-c(1)],n_type,n_type)
                T_ = matrix(sum_lm$coefficients[-c(1),3],n_type,n_type)
              }else{
                beta[c(1:n_type^2)[-null]] = lm$coefficients[-c(1)]
                T_[c(1:n_type^2)[-null]] = sum_lm$coefficients[-c(1),3]
                }
            }
            #record test statitistic
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
          #return (list(T_ = matrix(NA,n_type,n_type),betas = matrix(NA,n_type,n_type), V = 0,nulls = c(1:n_type^2), valid = 0,liks = NA))
          skip_to_next <<- TRUE})
      }
    }
    if(valid == 1){
      return (list(T_ = Matrix::t(T_),betas = Matrix::t(beta), V = as.matrix(var_mat),nulls = null, valid = valid,liks = liks_val))
    }else{
      return (list(T_ = 0,betas = 0, V = 0,nulls = c(1:n_type^2), valid = 0,liks = liks_val))
    }

  }


  #starting Message
  print(paste0('Starting Niche-DE analysis with parameters C = ',C,', M = ',M,', gamma = ', gamma,'.'))
  #initialize list output
  object@niche_DE = vector(mode = "list", length = length(object@sigma))
  names(object@niche_DE) = object@sigma
  counter = 1
  valid = matrix(0,ncol(object@counts),length(object@sigma))
  #iterate over each sigma value
  for(sig in object@sigma){
    print(paste0('Performing Niche-DE analysis with kernel bandwidth:',sig,' (number ',counter,' out of ',length(object@sigma),' values)'))
    #get expression filter (gamma)
    CT_filter = apply(object@ref_expr,1,function(x){quantile(x,gamma)})
    #get number of genes
    ngene = ncol(object@counts)
    #save nuber of cell types
    n_type = ncol(object@num_cells)
    #save dimnames for arrays
    dimnames = list(A = colnames(object@num_cells),B  = colnames(object@num_cells), C = colnames(object@counts))
    dimdims = c(ncol(object@num_cells),ncol(object@num_cells),ncol(object@counts))
    #prepare parallelization
    num_cores <- cores
    # Set up the parallel backend using doParallel
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    # Use foreach loop to parallelize "hello" function
    results <- foreach::foreach(i = 1:ngene)%dopar% {
      niche_DE_core(object,i,sig,CT_filter,C,M,gamma,Int)
    }
    #close cluster
    doParallel::stopImplicitCluster()
    #get likelihood list
    liks = unlist(lapply(results, function(result) result$liks))
    #get null entries
    nulls = lapply(results, function(result) result$nulls)
    names(nulls) = colnames(object@counts)
    #get variance covariance list
    var_cov = lapply(results, function(result) result$nulls)
    names(var_cov) = colnames(object@counts)
    #get T_stat_array
    T_stat_list = lapply(results, function(result) result$T_)
    #convert to an array
    T_stat = array(NA, dim = dimdims,dimnames = dimnames)
    dimnames(T_stat) =  dimnames

    #get beta array
    betas_list = lapply(results, function(result) result$betas)
    #convert to an array
    betas <- array(NA, dim = dimdims,dimnames = dimnames)
    dimnames(betas) =  dimnames
    #loop through to fill betas and T_stat
    for(vv in c(1:ngene)){
      T_stat[,,vv] = T_stat_list[[vv]]
      betas[,,vv] = betas_list[[vv]]
    }
    #save valid matrix
    valid[,counter] = unlist(lapply(results, function(result) result$valid))
    #save object
    object@niche_DE[[counter]] = list(T_stat = T_stat,beta = betas,var_cov = var_cov,nulls = nulls,log_lik = liks)
    counter = counter + 1
  }
  #get column sums of counts matrix to see how many genes pass filtering
  A = rowSums(valid)
  #get number of genes that pass filtering
  num_pass = sum(A>=1,na.rm = T)
  print('Computing Niche-DE Pvalues')
  object = get_niche_DE_pval(object,pos = T)
  object = get_niche_DE_pval(object,pos = F)
  print(paste0('Niche-DE analysis complete. Number of Genes with niche-DE T-stat equal to ',num_pass))
  if(num_pass < 1000){
    warning('Less than 1000 genes pass. This could be due to insufficient read depth of data or size of C parameter. Consider changing choice of C parameter')
  }
  return(object)
}

#' Getniche genes for the given index and niche cell types at the desired test.level
#'
#' This function returns genes that show niche patterns at the desired test.level
#'
#' @param object A niche-DE object
#' @param test.level At which test.level to return genes
#' (gene level, cell type level, interaction level)
#' @param index The index cell type
#' @param niche The niche cell type
#' @param direction Character indicating whether to return genes that are (index,niche)+
#' patterns (direction = 'positive') or (index,niche)- (direction = 'negative'). Default value is 'positive'.
#' @param alpha The level at which to perform the Benjamini Hochberg correction. Default value = 0.05
#' @return A vector of genes that are niche significant at the desired FDR,
#' test.level, index cell type, and niche cell type
#' @export
get_niche_DE_genes = function(object,test.level,index,niche,direction = 'positive',alpha = 0.05){
  if((test.level %in% c('gene','cell type','interaction'))==F){
    stop('test.level must be one of gene, cell type, or interaction')
  }

  if((direction %in% c('positive','negative'))==F){
    stop('direction must be one of positive or negative')
  }

  if(direction == 'positive'){
    S = '+'
  }else{
    S = '-'
  }
  if(test.level == 'interaction'){
    print(paste0('Finding Niche-DE',S,' genes at the interaction level between index cell type ',index,' and niche cell type '
           ,niche,'. Performing BH procedure at level ',alpha,'.'))
  }
  if(test.level == 'cell type'){
    print(paste0('Finding Niche-DE',' genes at the cell type level in index cell type ',index,'. Performing BH procedure at level ',alpha,'.'))
  }
  if(test.level == 'gene'){
    print(paste0('Finding Niche-DE',' genes at the gene level','. Performing BH procedure at level ',alpha,'.'))
  }

  #if test.level if gene level
  if(test.level=='gene' & direction == 'positive'){
    #get genes that reject at gene level
    gene_ind = which(object@niche_DE_pval_pos$gene_level<(alpha))
    genes = object@gene_names[gene_ind]
    #get associated pvalues
    pval = object@niche_DE_pval_pos$gene_level[gene_ind]
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Gene')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }

  if(test.level=='gene' & direction == 'negative'){
    #get genes that reject at gene level
    gene_ind = which(object@niche_DE_pval_neg$gene_level<(alpha))
    genes = object@gene_names[gene_ind]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$gene_level[gene_ind]
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Gene')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }



  if((index %in% colnames(object@num_cells))=='negative'){
    stop('Index cell type not found')
  }

  #get index and nice indices
  ct_index = which(colnames(object@num_cells)==index)
  niche_index = which(colnames(object@num_cells)==niche)
  if(length(ct_index)==0){
    stop("Index Cell Type Not Found")
  }
  if(length(niche_index)==0){
    stop("Niche Cell Type Not Found")
  }
  #check to see if they have enough overlap
  colloc = check_colloc(object,ct_index,niche_index)
  for(value in c(1:length(colloc))){
    if(colloc[value]<30){
      warning('Less than 30 observations containing collocalization of ',index,
              ' and ',niche,' at kernel bandwidth ',names(colloc)[value],
              '. Results may be unreliable.')
    }
  }
  #if test.level if cell type level
  if(test.level=='cell type' & direction == 'positive'){
    #get index of index cell type
    ct_index = which(colnames(object@num_cells)==index)
    #get genes that reject at the gene and CT level
    gene_index = which((object@niche_DE_pval_pos$gene_level<(alpha)) & (object@niche_DE_pval_pos$cell_type_level[,ct_index]<(alpha)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_pos$cell_type_level[gene_index,ct_index]
    #save results
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Cell.Type')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }

  if(test.level=='cell type' & direction == 'negative'){
    #get index of index cell type
    ct_index = which(colnames(object@num_cells)==index)
    #get genes that reject at the gene and CT level
    gene_index = which((object@niche_DE_pval_neg$gene_level<(alpha)) & (object@niche_DE_pval_neg$cell_type_level[,ct_index]<(alpha)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$cell_type_level[gene_index,ct_index]
    #save results
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Cell.Type')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }

  if((niche %in% colnames(object@num_cells))==F){
    stop('Niche cell type not found')
  }

  #if test.level if interaction level
  if(test.level =='interaction' & direction=='positive'){
    ct_index = which(colnames(object@num_cells)==index)
    niche_index = which(colnames(object@num_cells)==niche)
    gene_index = which((object@niche_DE_pval_pos$gene_level<(alpha)) &
                         (object@niche_DE_pval_pos$cell_type_level[,ct_index]<(alpha)) &
                         (object@niche_DE_pval_pos$interaction_level[ct_index,niche_index,]<(alpha/2)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_pos$interaction_level[ct_index,niche_index,gene_index]
    #save results
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Interaction')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }
  if(test.level=='interaction' & direction=='negative'){
    ct_index = which(colnames(object@num_cells)==index)
    niche_index = which(colnames(object@num_cells)==niche)
    gene_index = which((object@niche_DE_pval_neg$gene_level<(alpha)) &
                         (object@niche_DE_pval_neg$cell_type_level[,ct_index]<(alpha)) &
                         (object@niche_DE_pval_neg$interaction_level[ct_index,niche_index,]<(alpha/2)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$interaction_level[ct_index,niche_index,gene_index]
    #save results
    result = data.frame(genes,pval)
    if(dim(result)[1]==0){
      print("No Niche-DE Genes at this Resolution")
      return(result)
    }
    colnames(result) = c('Genes','Pvalues.Interaction')
    rownames(result) = c(1:nrow(result))
    print('Returning Niche-DE Genes')
    result = result[order(result[,2]),]
    return(result)
  }
}

#' Get Niche-DE marker genes
#'
#' This function returns genes marker genes in the index cell type when near the first niche cell type realtive to the second one
#'
#' param object A niche-DE object
#' param index The index cell type which we want to find marker genes for
#' param niche1 The niche cell type for the marker genes found
#' param niche2 The niche we wish to compare (index,niche1) patterns to
#' param pos Logical indicating whether to return genes that are (index,niche)+
#' patterns (pos = T) or (index,niche)- (pos = F)
#' param alpha The level at which to perform the Benjamini Hochberg correction. Default value is 0.05.
#' return A vector of genes that are niche marker genes for the index cell type
#'  near the niche1 cell type relative to the niche2 cell type

# #niche_DE_markers_old = function(object,index,niche1,niche2,alpha = 0.05){
#
#
#   if((index %in% colnames(object@num_cells))==F){
#     stop('Index cell type not found')
#   }
#   if((niche1 %in% colnames(object@num_cells))==F){
#     stop('Niche1 cell type not found')
#   }
#   if((niche2 %in% colnames(object@num_cells))==F){
#     stop('Niche2 cell type not found')
#   }
#
#   print(paste0('Finding Niche-DE marker genes in index cell type ',index,' with niche cell type ',niche1,
#                ' relative to niche cell type ',niche2,'. BH procedure performed at level ',alpha,'.'))
#
#   #get beta array
#   betas_all = object@niche_DE[[1]]$beta
#   #get variance covariance array
#   v_cov_all = object@niche_DE[[1]]$var_cov
#   #get index for index and niche cell types
#   index_index = which(colnames(object@num_cells)==index)
#   niche1_index = which(colnames(object@num_cells)==niche1)
#   niche2_index = which(colnames(object@num_cells)==niche2)
#
#   #make sure that collocalization occurs
#   #check to see if they have enough overlap
#   colloc = check_colloc(object,index_index,niche1_index)
#   for(value in c(1:length(colloc))){
#     if(colloc[value]<30){
#       warning('Less than 30 observations containing collocalization of ',index,
#               ' and ',niche1,' at kernel bandwidth ',names(colloc)[value],
#               '. Results may be unreliable.')
#     }
#   }
#
#   #make sure that collocalization occurs
#   #check to see if they have enough overlap
#   colloc = check_colloc(object,index_index,niche2_index)
#   for(value in c(1:length(colloc))){
#     if(colloc[value]<30){
#       warning('Less than 30 observations containing collocalization of ',index,
#               ' and ',niche2,' at kernel bandwidth ',names(colloc)[value],
#               '. Results may be unreliable.')
#     }
#   }
#
#
#
#
#   #get marker pvals
#   pval = contrast_post(betas_all,v_cov_all,index_index,c(niche1_index,niche2_index))
#   #if multiple kernels do this for all kernels
#   if(length(object@sigma)>=2){
#     for(j in c(2:length(object@sigma))){
#       #print(j)
#       betas_all = object@niche_DE[[j]]$beta
#       v_cov_all = object@niche_DE[[j]]$var_cov
#       index_index = which(colnames(object@num_cells)==index)
#       niche1_index = which(colnames(object@num_cells)==niche1)
#       niche2_index = which(colnames(object@num_cells)==niche2)
#       pval = cbind(pval,contrast_post(betas_all,v_cov_all,index_index,c(niche1_index,niche2_index)))
#     }
#   }
#   #apply cauchy combination
#   #record log likelihoods
#   log_liks = object@niche_DE[[1]]$log_lik
#   if(length(object@niche_DE)>=2){
#     for(j in c(2:length(object@niche_DE))){
#       log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
#     }
#   }
#   if(length(object@niche_DE)>=2){
#     log_liks[is.infinite(log_liks)] = 0
#     suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-min(x,na.rm = T))}))})
#   }
#
#   if(length(object@niche_DE)==1){
#     log_liks[is.infinite(log_liks)] = 0
#     suppressWarnings({ W = rep(1,length(log_liks))})
#   }
#
#   #W = apply(W,1,function(x){x/sum(x)})
#   W[is.infinite(W)] = 10e160
#   #bind pvalues and weights
#   contrast = cbind(pval,W)
#   #apply cauchy rule
#   contrast = as.matrix(apply(contrast,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])}))
#   #apply BH
#   contrast_phoch = p.adjust(contrast,method = "BH")
#   #bind genes and their pvalue
#   gene_pval = data.frame(object@gene_names,contrast_phoch)
#   #filter to only those that reject
#   gene_pval = gene_pval[which(gene_pval[,2]<(alpha/2)),]
#   colnames(gene_pval) = c('Genes','Adj.Pvalues')
#   rownames(gene_pval) = c(1:nrow(gene_pval))
#   print('Marker gene analysis complete.')
#   gene_pval = gene_pval[order(gene_pval[,2]),]
#   return(gene_pval)
#
# }


#' Get Niche-DE marker genes(test)
#'
#' This function returns genes marker genes in the index cell type when near the first niche cell type realtive to the second one
#'
#' @param object A niche-DE object
#' @param index The index cell type which we want to find marker genes for
#' @param niche1 The niche cell type for the marker genes found
#' @param niche2 The niche we wish to compare (index,niche1) patterns to
#' @param alpha The level at which to perform the Benjamini Hochberg correction. Default value is 0.05.
#' @return A vector of genes that are niche marker genes for the index cell type
#'  near the niche1 cell type relative to the niche2 cell type
#' @export
niche_DE_markers = function(object,index,niche1,niche2,alpha = 0.05){


  if((index %in% colnames(object@num_cells))==F){
    stop('Index cell type not found')
  }
  if((niche1 %in% colnames(object@num_cells))==F){
    stop('Niche1 cell type not found')
  }
  if((niche2 %in% colnames(object@num_cells))==F){
    stop('Niche2 cell type not found')
  }

  print(paste0('Finding Niche-DE marker genes in index cell type ',index,' with niche cell type ',niche1,
               ' relative to niche cell type ',niche2,'. BH procedure performed at level ',alpha,'.'))

  #get beta array
  betas_all = object@niche_DE[[1]]$beta
  #get variance covariance array
  v_cov_all = object@niche_DE[[1]]$var_cov
  nulls_all = object@niche_DE[[1]]$nulls
  #get index for index and niche cell types
  index_index = which(colnames(object@num_cells)==index)
  niche1_index = which(colnames(object@num_cells)==niche1)
  niche2_index = which(colnames(object@num_cells)==niche2)

  #make sure that collocalization occurs
  #check to see if they have enough overlap
  colloc = check_colloc(object,index_index,niche1_index)
  for(value in c(1:length(colloc))){
    if(colloc[value]<30){
      warning('Less than 30 observations containing collocalization of ',index,
              ' and ',niche1,' at kernel bandwidth ',names(colloc)[value],
              '. Results may be unreliable.')
    }
  }

  #make sure that collocalization occurs
  #check to see if they have enough overlap
  colloc = check_colloc(object,index_index,niche2_index)
  for(value in c(1:length(colloc))){
    if(colloc[value]<30){
      warning('Less than 30 observations containing collocalization of ',index,
              ' and ',niche2,' at kernel bandwidth ',names(colloc)[value],
              '. Results may be unreliable.')
    }
  }




  #get marker pvals
  pval = contrast_post(betas_all,v_cov_all,nulls_all,index_index,c(niche1_index,niche2_index))
  #if multiple kernels do this for all kernels
  if(length(object@sigma)>=2){
    for(j in c(2:length(object@sigma))){
      #print(j)
      betas_all = object@niche_DE[[j]]$beta
      v_cov_all = object@niche_DE[[j]]$var_cov
      nulls_all = object@niche_DE[[j]]$nulls
      index_index = which(colnames(object@num_cells)==index)
      niche1_index = which(colnames(object@num_cells)==niche1)
      niche2_index = which(colnames(object@num_cells)==niche2)
      pval = cbind(pval,contrast_post(betas_all,v_cov_all,nulls_all,index_index,c(niche1_index,niche2_index)))
    }
  }
  #apply cauchy combination
  #record log likelihoods
  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  if(length(object@niche_DE)>=2){
    log_liks[is.infinite(log_liks)] = 0
    suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-min(x,na.rm = T))}))})
  }

  if(length(object@niche_DE)==1){
    log_liks[is.infinite(log_liks)] = 0
    suppressWarnings({ W = rep(1,length(log_liks))})
  }

  #W = apply(W,1,function(x){x/sum(x)})
  W[is.infinite(W)] = 10e160
  #bind pvalues and weights
  contrast = cbind(pval,W)
  #apply cauchy rule
  contrast = as.matrix(apply(contrast,1,function(x){gene_level(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])}))
  #apply BH
  contrast_phoch = p.adjust(contrast,method = "BH")
  #bind genes and their pvalue
  gene_pval = data.frame(object@gene_names,contrast_phoch)
  #filter to only those that reject
  gene_pval_ = gene_pval[which(gene_pval[,2]<(alpha/2)),]
  if(nrow(gene_pval_)==0){
    print("No Significant Genes")
    return(gene_pval)
  }
  colnames(gene_pval_) = c('Genes','Adj.Pvalues')
  rownames(gene_pval_) = c(1:nrow(gene_pval_))
  print('Marker gene analysis complete.')
  gene_pval = gene_pval_[order(gene_pval_[,2]),]
  return(gene_pval_)

}






#' Perform Niche-LR (Ligand receptor analysis) on spot level data
#'
#' This function returns ligands and receptors inferred to be expressed by the given cell types
#'
#' @param object A niche-DE object
#' @param ligand_cell The cell type that expresses the ligand
#' @param receptor_cell The cell type that expresses the receptor
#' @param ligand_target_matrix A matrix that measures the association between
#' ligands and their downstream target genes. Should be target genes by ligands
#' @param lr_mat A matrix that matches ligands with their corresponding receptors.
#' This matrix should have two columns. The first will be ligands and the second
#' will be the corresponding receptors
#' @param K The number of downstream target genes to use when calculating the
#' ligand potential score. Default value is 25.
#' @param M The maximum number of ligands that can pass initial filtering. Default value is 50.
#' @param alpha The level at which to perform the Benjamini Hochberg correction. Default value is 0.05.
#' @param truncation_value The value at which to truncate T statistics. Default value is 3.
#' @return A list of ligand-receptor pairs that are found to be expressed by the
#' specified cell type
#' @export
niche_LR_spot = function(object,ligand_cell,receptor_cell,ligand_target_matrix,lr_mat,K = 25,M = 50,alpha = 0.05,truncation_value = 3){
  #The ligand expressing cell should be the niche cell
  niche = which(colnames(object@num_cells)==ligand_cell)
  #The receptor expressing cell should be the index cell
  index = which(colnames(object@num_cells)==receptor_cell)

  print(paste0('Performing niche-LR with hyperparameters K = ',K,', M = ',M,', alpha = ',alpha,
               ', truncation value = ',truncation_value,'.'))
  print(paste0('Finding ligand receptors between ligand expressing cell type ',ligand_cell,
               ' and receptor expressing cell type ',receptor_cell,'.'))
  #make sure that collocalization occurs
  #check to see if they have enough overlap
  colloc = check_colloc(object,index,niche)
  for(value in c(1:length(colloc))){
    if(colloc[value]<30){
      warning('Less than 30 observations containing collocalization of ',receptor_cell,
              ' and ',ligand_cell,' at kernel bandwidth ',names(colloc)[value],
              '. Results may be unreliable.')
    }
  }


  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }

  if(length(object@niche_DE)>=2){
    #get best kernel for each gene
    top_kernel = apply(log_liks,1,function(x){order(x,decreasing = T)[1]})
  }

  if(length(object@niche_DE)==1){
    #get best kernel for each gene
    top_kernel = rep(1,length(log_liks))
  }

  #get T_statistic list
  T_vector = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_vector[[j]] = object@niche_DE[[j]]$T_stat[index,niche,]
  }
  #truncate niche_DE t-statistic
  for(j in c(1:length(T_vector))){
    T_vector[[j]] = pmin(T_vector[[j]],abs(truncation_value))
    T_vector[[j]] = pmax(T_vector[[j]],-abs(truncation_value))
  }
  #get library reference
  L = object@ref_expr
  #get cutoffs to consider a ligand for each cell type
  CT_filter = apply(L,1,function(x){quantile(x,0.25)})
  #filter ligands
  cand_lig = colnames(L)[which(L[niche,]>CT_filter[niche])]
  filter = which(colnames(ligand_target_matrix)%in% cand_lig)
  ligand_target_matrix = ligand_target_matrix[,filter]

  print('Calculating ligand potential scores')
  #get ligand potential scores
  pear_cor = rep(NA,ncol(ligand_target_matrix))
  score_norm = rep(NA,ncol(ligand_target_matrix))

  #save top niche-DE genes
  top_DE = vector(mode = "list", length = ncol(ligand_target_matrix))
  names(top_DE) = colnames(ligand_target_matrix)
  #iterate over potential ligands
  for(j in c(1:ncol(ligand_target_matrix))){
    ligand = colnames(ligand_target_matrix)[j]
    if(ligand%in% colnames(L)){
      #get which index we need to query for best kernel
      ind = which(colnames(L)==ligand)
      #get best index
      top_index = top_kernel[ind]
      #get scores
      sig = T_vector[[top_index]]
      genes = object@gene_names
      filter = which(is.na(sig)==F)
      genes = genes[filter]
      sig = sig[filter]
      #get ligand_potential scores by summing up T_statistics weighted by scaled niche-net scores
      ligand_vec = ligand_target_matrix
      #filter to only include genes found in data
      filter = which(rownames(ligand_vec)%in% genes)
      #filter
      ligand_vec = ligand_vec[filter,]
      filter = which(genes%in%rownames(ligand_vec))
      sig = sig[filter]
      genes = genes[filter]
      ligand_vec = ligand_vec[genes,]
      ligand_vec[,j] = scale(ligand_vec[,j])
      #get top K downstream genes
      top_cors = order(ligand_vec[,j],decreasing = T)[1:K]
      #get weights based on scaled niche-net scores
      weight = ligand_vec[top_cors,j]/mean(ligand_vec[top_cors,j])
      ##calculate ligand potential scores
      pear_cor[j] = sum(sig[top_cors]*weight)
      top_DE[[j]] = paste((rownames(ligand_vec)[top_cors])[order(sig[top_cors]*weight,decreasing = T)[1:5]],collapse = ',')
      #get normalizing constant
      score_norm[j] = sum(weight^2)
    }
  }
  #get scaled ligand-potential scores
  pear_cor = pear_cor/sqrt(score_norm)
  #get top candidate ligands
  top_index = which(pear_cor>=max(1.64,sort(pear_cor,decreasing = T)[M]))
  top_scores= pear_cor[top_index]
  #get top candidate ligands
  top_genes = colnames(ligand_target_matrix)[top_index]
  if(length(top_genes)==0){
    stop('No candidate ligands')
  }
  #################### test candidate ligands
  print(paste0('Testing candidate ligands for sufficient expression in cell type ',ligand_cell))
  pvalues = c()
  beta = c()
  for(j in c(1:length(top_genes))){
    tryCatch({
      #get candidate ligand
      gene = top_genes[j]
      #get what index it belongs to
      index_gene = which(colnames(L)==gene)
      #get best kernel
      top_index = top_kernel[index_gene]
      #get counts data for that gene
      data = object@counts
      Y = (data[,which(colnames(L)==gene)])
      #which spots to look at
      #niche_index = which(kernel_materials[[top_index]]$EN[,index]>0)
      niche_index = which(object@effective_niche[[top_index]][,index]>min(object@effective_niche[[top_index]][,index]))
      #filter data to only look at these spots
      Y = Y[niche_index]
      #get number of cells per spot
      nst = object@num_cells
      nst = nst[niche_index,]
      nst[,which(L[,index_gene] < CT_filter)] = 0
      #run regression of ligand expression on number of cells per spot
      if(L[niche,index_gene] < CT_filter[niche]){
        beta = c(beta,NA)
        pvalues = c(pvalues,NA)
      }else{
        check = suppressWarnings({glm(Y~nst,family = 'poisson')})
        bad_ind = which(is.na(coef(check)))
        num_bad = sum(bad_ind<(niche+1))
        beta = c(beta,coef(check)[niche+1])
        pvalues = c(pvalues,summary(check)$coefficients[(niche+1-num_bad),4])
      }

    } #get pval
    , error = function(e) {
      skip_to_next <<- TRUE})
  }
  #adjust pvalues and get confirmed ligands
  pvalues = p.adjust(pvalues,method = 'BH')
  ligands = top_genes[which(pvalues<alpha &beta>0)]

  #get candidate receptor
  rec_ind = which(lr_mat[,1]%in% ligands & lr_mat[,2]%in% colnames(L))
  lr_mat = lr_mat[rec_ind,]
  receptors = lr_mat[,2]
  #run regression to see if receptor is expressed by index cell type
  pvalues = c()
  beta = c()
  gene_name = c()
  #iterate over all receptors
  print(paste0('Testing candidate receptors for sufficient expression in cell type ',receptor_cell))
  for(j in c(1:length(receptors))){
    tryCatch({
      #get receptor
      gene = receptors[j]
      #get corresponding ligand
      lig = lr_mat[j,1]
      #get index for the ligand and the optimal kernel
      index_lig = which(colnames(L)==lig)
      top_index = top_kernel[index_lig]
      #get data
      data = object@counts
      Y = data[,which(colnames(L)==gene)]
      #which spots to look at
      #niche_index = which(kernel_materials[[top_index]]$EN[,niche]>0)
      niche_index = which(object@effective_niche[[top_index]][,niche]>min(object@effective_niche[[top_index]][,niche]))
      #filter based on spots to look at
      Y = Y[niche_index]
      nst = object@num_cells[niche_index,]
      nst[,which(L[,index_lig] < CT_filter)] = 0
      #run regression
      if(L[index,index_lig] < CT_filter[index]){
        beta = c(beta,NA)
        pvalues = c(pvalues,NA)
      }else{
        check = suppressWarnings({glm(Y~nst,family = 'poisson')})
        bad_ind = which(is.na(coef(check)))
        num_bad = sum(bad_ind<(index+1))
        beta = c(beta,coef(check)[index+1])
        pvalues = c(pvalues,summary(check)$coefficients[(index+1-num_bad),4])
      }
    }
    , error = function(e) {
      print(paste0("error",j))
      skip_to_next <<- TRUE})
  }
  #adjust pvalues
  pvalues = p.adjust(pvalues,method = "BH")
  rec_mat = cbind(lr_mat,beta,pvalues)
  #look at those with postive beta values
  rec_mat = rec_mat[as.numeric(rec_mat[,3])>0,]
  #name matrix columns
  colnames(rec_mat) = c('ligand','receptor','receptor_beta','receptor_pval')
  LR_pairs = rec_mat[,c(1,2,4)]
  #extract the confirmed ligands and receptors
  if(length(LR_pairs)==0){
    stop('no ligand-receptor pairs to report')
  }
  LR_pairs = LR_pairs[which(as.numeric(LR_pairs[,3])<alpha),c(1:3)]
  for(j in unique(LR_pairs[,1])){
    #get gene/ligand name
    #match to top_DE list
    ID = which(names(top_DE)==j)
    print(ID)
    #get downstream genes
    LR_pairs[which(LR_pairs[,1]==j),3] = top_DE[[ID]]
  }

  colnames(LR_pairs) = c('ligand','receptor','top_downstream_niche_DE_genes')
  rownames(LR_pairs) = c(1:nrow(LR_pairs))

  print(paste0('Returning ligand receptor table between ligand expressing cell type ',ligand_cell,
               ' and receptor expressing cell type ',receptor_cell,'.'))
  return(LR_pairs)

}

#' Perform Niche-LR (Ligand receptor analysis) on single cell level data
#'
#' This function returns ligands and receptors inferred to be expressed by the given cell types
#'
#' @param object A niche-DE object
#' @param ligand_cell The cell type that expresses the ligand
#' @param receptor_cell The cell type that expresses the receptor
#' @param ligand_target_matrix A matrix that measures the assocaition between
#' ligands and their downstream target genes. Should be target genes by ligands
#' @param lr_mat A matrix that matches ligands with their corresponding receptors.
#' This matrix should have two columns. The first will be ligands and the second
#' will be the corresponding receptors
#' @param K The number of downstream target genes to use when calculating the
#' ligand potential score. Default value is 25.
#' @param M The maximum number of ligands that can pass initial filtering. Default value is 50.
#' @param alpha The level at which to perform the Benjamini Hochberg correction. Default value is 0.05.
#' @param alpha_2 The null quantile to compare observed expression to. Default value is 0.5 (50th percentile).
#' @param truncation_value The value at which to truncate T statistics. Default value is 3.
#' @return A list of ligand-receptor pairs that are found to be expressed by the
#' specified cell type
#' @export
niche_LR_cell = function(object,ligand_cell,receptor_cell,ligand_target_matrix,
                         lr_mat,K = 25,M = 50,alpha = 0.05,alpha_2 = 0.5,truncation_value = 3){
  niche = which(colnames(object@num_cells)==ligand_cell)
  index = which(colnames(object@num_cells)==receptor_cell)
  print(paste0('Performing niche-LR with hyperparameters K = ',K,', M = ',M,', alpha = ',alpha,', alpha_2 = ',alpha_2,
               ', truncation value = ',truncation_value,'.'))
  print(paste0('Finding ligand receptors between ligand expressing cell type ',ligand_cell,
               ' and receptor expressing cell type ',receptor_cell,'.'))

  #make sure that collocalization occurs
  #check to see if they have enough overlap
  colloc = check_colloc(object,index,niche)
  for(value in c(1:length(colloc))){
    if(colloc[value]<30){
      warning('Less than 30 observations containing collocalization of ',receptor_cell,
              ' and ',ligand_cell,' at kernel bandwidth ',names(colloc)[value],
              '. Results may be unreliable.')
    }
  }


  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  if(length(object@niche_DE)>=2){
    #get best kernel for each gene
    top_kernel = apply(log_liks,1,function(x){order(x,decreasing = T)[1]})
  }

  if(length(object@niche_DE)==1){
    #get best kernel for each gene
    top_kernel = rep(1,length(log_liks))
  }

  #get T_statistic list
  T_vector = vector(mode = "list", length = length(object@sigma))
  for(j in c(1:length(object@sigma))){
    T_vector[[j]] = object@niche_DE[[j]]$T_stat[index,niche,]
  }
  #truncate niche_DE t-statistic
  for(j in c(1:length(T_vector))){
    T_vector[[j]] = pmin(T_vector[[j]],abs(truncation_value))
    T_vector[[j]] = pmax(T_vector[[j]],-abs(truncation_value))
  }
  #get library reference
  L = object@ref_expr
  #get cutoffs to consider a ligand for each cell type
  CT_filter = apply(L,1,function(x){quantile(x,0.25)})
  #filter ligands
  cand_lig = colnames(L)[which(L[niche,]>CT_filter[niche])]
  filter = which(colnames(ligand_target_matrix)%in% cand_lig)
  ligand_target_matrix = ligand_target_matrix[,filter]

  print('Calculating ligand potential scores')
  #get ligand potential scores
  pear_cor = rep(NA,ncol(ligand_target_matrix))
  score_norm = rep(NA,ncol(ligand_target_matrix))

  #save top niche-DE genes
  top_DE = vector(mode = "list", length = ncol(ligand_target_matrix))
  names(top_DE) = colnames(ligand_target_matrix)

  #iterate over potential ligands
  for(j in c(1:ncol(ligand_target_matrix))){
    ligand = colnames(ligand_target_matrix)[j]
    if(ligand%in% colnames(L)){
      #get which index we need to query for best kernel
      ind = which(colnames(L)==ligand)
      #get best index
      top_index = top_kernel[ind]
      #get scores
      sig = T_vector[[top_index]]
      genes = object@gene_names
      filter = which(is.na(sig)==F)
      genes = genes[filter]
      sig = sig[filter]
      #get ligand_potential scores by summing up T_statistics weighted by scaled niche-net scores
      ligand_vec = ligand_target_matrix
      #filter to only include genes found in data
      filter = which(rownames(ligand_vec)%in% genes)
      #filter
      ligand_vec = ligand_vec[filter,]
      filter = which(genes%in%rownames(ligand_vec))
      sig = sig[filter]
      genes = genes[filter]
      ligand_vec = ligand_vec[genes,]
      ligand_vec[,j] = scale(ligand_vec[,j])
      #get top K downstreamm genes
      top_cors = order(ligand_vec[,j],decreasing = T)[1:K]
      #get weights based on scaled niche-net scors
      weight = ligand_vec[top_cors,j]/mean(ligand_vec[top_cors,j])
      ##calculate ligand potential scores
      pear_cor[j] = sum(sig[top_cors]*weight)
      #save top niche-DE genes of the ligand
      top_DE[[j]] = paste((rownames(ligand_vec)[top_cors])[order(sig[top_cors]*weight,decreasing = T)[1:5]],collapse = ',')
      #get normalizing constant
      score_norm[j] = sum(weight^2)
    }
  }
  #get scaled ligand-potential scores
  pear_cor = pear_cor/sqrt(score_norm)
  #get top candidate ligands
  top_index = which(pear_cor>=max(1.64,sort(pear_cor,decreasing = T)[M]))
  top_scores= pear_cor[top_index]
  #get top candidate ligands
  top_genes = colnames(ligand_target_matrix)[top_index]
  #print(top_genes)

  #################### test candidate ligands
  print(paste0('Testing candidate ligands for sufficient expression in cell type ',ligand_cell))
  pvalues = c()
  beta = c()
  for(j in c(1:length(top_genes))){
    tryCatch({
      #get candidate ligand
      gene = top_genes[j]
      #get what index it belongs to
      index_gene = which(colnames(L)==gene)
      #get best kernel
      top_index = top_kernel[index_gene]
      #get counts data for that gene
      data = object@counts
      Y = (data[,which(colnames(L)==gene)])
      #which spots to look at
      #niche_index = which(kernel_materials[[top_index]]$EN[,index]>0)
      niche_index = which(object@effective_niche[[top_index]][,index]>min(object@effective_niche[[top_index]][,index]))
      #filter data to only look at these spots
      Y = Y[niche_index]
      #get number of cells per spot
      nst = object@num_cells
      nst = nst[niche_index,]
      #print(dim(nst))
      nst[,which(L[,index_gene] < CT_filter)] = 0
      #run regression of ligand expression on number of cells per spot
      #get average expression of ligand
      lambda_niche = mean(Y[nst[,niche]==1])
      #print(niche)
      #print(lambda_niche)
      #get alpha quantile of niche cell type
      lambda_rest = quantile(object@ref_expr[niche,],alpha_2)
      #print(lambda_rest)
      #get test statistic
      Test_stat = lambda_niche-lambda_rest
      #normalize so that null distribution is standard normal
      Test_stat = Test_stat/sqrt(lambda_niche/sum(nst[,niche]==1))
      #append pvalues
      pvalues = c(pvalues,1-pnorm(Test_stat))

    } #get pval
    , error = function(e) {
      #print(paste0("error",j))
      skip_to_next <<- TRUE})
  }
  #adjust pvalues and get confirmed ligands
  pvalues = p.adjust(pvalues,method = 'BH')
  ligands = top_genes[which(pvalues<alpha)]
  #print(ligands)
  #get candidate receptor
  rec_ind = which(lr_mat[,1]%in% ligands & lr_mat[,2]%in% colnames(L))
  lr_mat = lr_mat[rec_ind,]
  receptors = lr_mat[,2]
  #run regression to see if receptor is expressed by index cell type
  pvalues = c()
  beta = c()
  gene_name = c()
  #iterate over all receptors
  print(paste0('Testing candidate receptors for sufficient expression in cell type ',receptor_cell))
  for(j in c(1:length(receptors))){
    tryCatch({
      #get receptor
      gene = receptors[j]
      #get corresponding ligand
      lig = lr_mat[j,1]
      #get index for the ligand and the optimal kernel
      index_lig = which(colnames(L)==lig)
      top_index = top_kernel[index_lig]
      #get data
      data = object@counts
      Y = data[,which(colnames(L)==gene)]
      #which spots to look at
      #niche_index = which(kernel_materials[[top_index]]$EN[,niche]>0)
      niche_index = which(object@effective_niche[[top_index]][,niche]>min(object@effective_niche[[top_index]][,niche]))
      #filter based on spots to look at
      Y = Y[niche_index]
      #get number of cells per spot
      nst = object@num_cells
      nst = nst[niche_index,]
      nst[,which(L[,index_lig] < CT_filter)] = 0
      #run regression
      lambda_niche = mean(Y[nst[,index]==1])
      lambda_rest = quantile(object@ref_expr[index,],alpha_2)
      Test_stat = lambda_niche-lambda_rest
      Test_stat = Test_stat/sqrt(lambda_niche/sum(nst[,index]==1))
      pvalues = c(pvalues,1-pnorm(Test_stat))
    }
    , error = function(e) {
      print(paste0("error",j))
      skip_to_next <<- TRUE})
  }
  #adjust pvalues
  pvalues = p.adjust(pvalues,method = "BH")
  #print(pvalues)
  rec_mat = cbind(lr_mat,pvalues)
  #name matrix columns
  colnames(rec_mat) = c('ligand','receptor','receptor_pval')
  LR_pairs = rec_mat[,c(1,2,3)]
  #print(dim(LR_pairs))
  #extract the confirmed ligands and receptors
  if(length(LR_pairs)==0){
    stop('no ligand-receptor pairs to report')
  }
  LR_pairs = LR_pairs[which(as.numeric(LR_pairs[,3])<alpha),c(1:3)]
  for(j in unique(LR_pairs[,1])){
    #match to top_DE list
    ID = which(names(top_DE)==j)
    #get downstream genes
    LR_pairs[which(LR_pairs[,1]==j),3] = top_DE[[ID]]
  }

  colnames(LR_pairs) = c('ligand','receptor','top_downstream_niche_DE_genes')
  rownames(LR_pairs) = c(1:nrow(LR_pairs))
  print(paste0('Returning ligand receptor table between ligand expressing cell type ',ligand_cell,
               ' and receptor expressing cell type ',receptor_cell,'.'))

  return(LR_pairs)

}

