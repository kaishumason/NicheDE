#' Niche_DE
#'
#' This function performs niche-DE
#'
#' @param C Minimum total expression of a gene needed for the model to run
#' @param M Minimum number of spots containing the index cell type with the 
#' niche cell type in its effective niche for (index,niche) niche patterns
#' to be investigated
#' @param gamma Percentile a gene needs to be with respect to expression in the
#'  index cell type in order for the model to investigate niche patterns for 
#'  that gene in the index cell
#' @return A niche-DE object with niche-DE analysis performed
#' @export
niche_DE = function(object,C = 150,M = 10,gamma = 0.8){
  #intialize list output
  object@niche_DE = vector(mode = "list", length = length(object@sigma))
  names(object@niche_DE) = object@sigma
  counter = 1
  #iterate over each sigma value
  for(sig in object@sigma){
    #get expression filter (gamma)
    CT_filter = apply(object@ref_expr,1,function(x){quantile(x,gamma)})
    #initialize p value array 
    ngene = ncol(object@counts)
    n_type = ncol(object@num_cells)
    dimnames = list(A = colnames(object@num_cells),B  = colnames(object@num_cells), C = colnames(object@counts))
    #pgt is index type by niche type by gene 
    T_stat = array(NA,c(n_type,n_type,ngene),dimnames = dimnames)
    var_cov = array(NA,c(n_type^2,n_type^2,ngene))
    betas = array(NA,c(n_type,n_type,ngene),dimnames = dimnames)
    liks = rep(NA,ngene)
    for(j in c(1:ngene)){
      if(j%%1000 == 0){
        print(paste0('kernel bandwidth:', sig, " Gene #",j, ' out of ',ncol(object@counts)))
      }
      #do if  gene is rejected and gene-type has at least 1 rejection
      if((sum(object@counts[,j])>C)&(mean(object@ref_expr[,j]<CT_filter)!=1)){
        #get pstg matrix 
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
        null = which(apply(X,2,function(x){sum(x>0)})<M)
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
            #print(2)
            if(length(bad_ind)>0){
              full_glm =suppressWarnings({glm(object@counts[-bad_ind,j]~X_partial[-bad_ind,] + offset(log(object@null_expected_expression[-bad_ind,j])), family = "poisson")}) #do full glm
            }else{
              full_glm = suppressWarnings({glm(object@counts[,j]~X_partial + offset(log(object@null_expected_expression[,j])), family = "poisson")}) #do full glm
            }
            mu_hat = exp(predict(full_glm))#get mean
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
            #get variance matrix
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
            #cholesky decomposition
            A = Matrix::chol(var_mat,LDL = FALSE,perm = FALSE)
            #get covaraince matrix
            V = solve(A)%*%Matrix::t(solve(A))
            #get standard devaition vector
            tau = sqrt(diag(V))#get sd matrix
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
            var_cov[,,j] = v_cov}, #get pval
              error = function(e) { 
              print(paste0("error",j))
              skip_to_next <<- TRUE}) 
        }
      }
    }
    #save object
    object@niche_DE[[counter]] = list(T_stat = T_stat,beta = betas,var_cov = var_cov,log_lik = liks)
    counter = counter + 1
  }
  return(object)
}

#' get_niche_DE_pval
#'
#' This function calculates pvalues from a niche-DE analysis
#'
#' @param object A niche-DE object
#' @param pos Logical indicating whether to calculate pvalues for(index,niche)+
#' patterns (pos = T) or (index,niche)- patterns (pos = F)
#' @return A niche-DE object with niche-DE pvalues
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


#' get_niche_DE_genes
#'
#' This function returns genes that show niche patterns at the desired resolution
#'
#' @param object A niche-DE object
#' @param resolution At which resolution to return genes 
#' (gene level, cell type level, interaction level)
#' @param index The index cell type 
#' @param niche The niche cell type 
#' @param pos Logical indicating whether to return genes that are (index,niche)+
#' patterns (pos = T) or (index,niche)- (pos = F)
#' @param alpha The level at which to perform the Benjamini Hochberg correction
#' @return A vector of genes that are niche significant at the desired FDR, 
#' resolution, index cell type, and niche cell type 
#' @export
get_niche_DE_genes = function(object,resolution,index,niche,pos,alpha){
  if((resolution %in% c('gene','cell type','interaction'))==F){
    stop('resolution must be one of gene, cell type, or interaction')
  }
  
  #if resolution if gene level
  if(resolution=='gene' & pos == T){
    #get genes that reject at gene level
    gene_ind = which(object@niche_DE_pval_pos$gene_level<(alpha))
    genes = object@gene_names[gene_ind]
    #get associated pvalues
    pval = object@niche_DE_pval_pos$gene_level[gene_ind]
    result = data.frame(genes,pval)
    colnames(result) = c('Genes','Adj.Gene level Pvalues')
    return(result)
  }
  
  if(resolution=='gene' & pos == F){
    #get genes that reject at gene level
    gene_ind = which(object@niche_DE_pval_neg$gene_level<(alpha))
    genes = object@gene_names[gene_ind]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$gene_level[gene_ind]
    result = data.frame(genes,pval)
    colnames(result) = c('Genes','Adj.Gene level Pvalues')
    return(result)
  }
  
  
  
  if((index %in% colnames(object@num_cells))==F){
    stop('Index cell type not found')
  }
  if((niche %in% colnames(object@num_cells))==F){
    stop('Niche cell type not found')
  }
  
  
  #if resolution if cell type level
  if(resolution=='cell type' & pos == T){
    #get index of index cell type
    ct_index = which(colnames(object@num_cells)==index)
    #get genes that reject at the gene and CT level
    gene_index = which((object@niche_DE_pval_pos$gene_level<(alpha)) & (object@niche_DE_pval_pos$cell_type_level[,ct_index]<(alpha)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_pos$cell_type_level[gene_index,ct_index]
    #save results
    result = data.frame(genes,pval)
    colnames(result) = c('Genes','Adj.Cell Type level Pvalues')
    return(result)
  }
  
  if(resolution=='cell type' & pos == F){
    #get index of index cell type
    ct_index = which(colnames(object@num_cells)==index)
    #get genes that reject at the gene and CT level
    gene_index = which((object@niche_DE_pval_neg$gene_level<(alpha)) & (object@niche_DE_pval_neg$cell_type_level[,ct_index]<(alpha)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$cell_type_level[gene_index,ct_index]
    #save results
    result = data.frame(genes,pval)
    colnames(result) = c('Genes','Adj.Cell Type level Pvalues')
    return(result)
  }
  
  
  #if resolution if interaction level
  if(resolution =='interaction' & pos==T){
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
    colnames(result) = c('Genes','Adj.Interaction level Pvalues')
    return(result)
  }
  if(resolution=='interaction' & pos==F){
    ct_index = which(colnames(object@num_cells)==index)
    niche_index = which(colnames(object@num_cells)==niche)
    gene_index = which((object@niche_DE_pval_pos$gene_level<(alpha)) &
                         (object@niche_DE_pval_neg$cell_type_level[,ct_index]<(alpha)) &
                         (object@niche_DE_pval_neg$interaction_level[ct_index,niche_index,]<(alpha/2)))
    genes = object@gene_names[gene_index]
    #get associated pvalues
    pval = object@niche_DE_pval_neg$interaction_level[ct_index,niche_index,gene_index]
    #save results
    result = data.frame(genes,pval)
    colnames(result) = c('Genes','Adj.Interaction level Pvalues')
    return(result)
  }
}

#' get_niche_DE_genes
#'
#' This function returns genes that show niche patterns at the desired resolution
#'
#' @param object A niche-DE object
#' @param index The index cell type which we want to find marker genes for
#' @param niche1 The niche cell type for the marker genes found
#' @param niche2 The niche we wish to compare (index,niche1) patterns to
#' @param pos Logical indicating whether to return genes that are (index,niche)+
#' patterns (pos = T) or (index,niche)- (pos = F)
#' @param alpha The level at which to perform the Benjamini Hochberg correction
#' @return A vector of genes that are niche marker genes for the index cell type
#'  near the niche1 cell type relative to the niche2 cell type 
#' @export
niche_DE_markers = function(object,index,niche1,niche2,alpha){
  
  if((index %in% colnames(object@num_cells))==F){
    stop('Index cell type not found')
  }
  if((niche1 %in% colnames(object@num_cells))==F){
    stop('Niche1 cell type not found')
  }
  if((niche2 %in% colnames(object@num_cells))==F){
    stop('Niche2 cell type not found')
  }
  
  #get beta array 
  betas_all = object@niche_DE[[1]]$beta
  #get variance covariance array
  v_cov_all = object@niche_DE[[1]]$var_cov
  #get index for index and niche cell types
  index_index = which(colnames(object@num_cells)==index)
  niche1_index = which(colnames(object@num_cells)==niche1)
  niche2_index = which(colnames(object@num_cells)==niche2)
  #get marker pvals
  pval = contrast_post(betas_all,v_cov_all,index_index,c(niche1_index,niche2_index))
  #if multiple kernels do this for all kernels
  if(length(object@sigma)>=2){
    for(j in c(2:length(object@sigma))){
      #print(j)
      betas_all = object@niche_DE[[j]]$beta
      v_cov_all = object@niche_DE[[j]]$var_cov
      index_index = which(colnames(object@num_cells)==index)
      niche1_index = which(colnames(object@num_cells)==niche1)
      niche2_index = which(colnames(object@num_cells)==niche2)
      pval = cbind(pval,contrast_post(betas_all,v_cov_all,index_index,c(niche1_index,niche2_index)))
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
  log_liks[is.infinite(log_liks)] = 0
  suppressWarnings({ W = t(apply(log_liks,1,function(x){exp(x-min(x,na.rm = T))}))})
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
  gene_pval = gene_pval[which(gene_pval[,2]<(alpha/2)),]
  colnames(gene_pval) = c('Genes','Adj.Pvalues')
  return(gene_pval)
  
}


#' get_niche_DE_genes
#'
#' This function returns genes that show niche patterns at the desired resolution
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
#' ligand potential score
#' @param M The maximum number of ligands that can pass initial filtering
#' @param alpha The level at which to perform the Benjamini Hochberg correction
#' @param truncation_value The value at which to truncate T statistics
#' @return A list of ligand-receptor pairs that are found to be expressed by the
#' specified cell type
#' @export
niche_LR_spot = function(object,ligand_cell,receptor_cell,ligand_target_matrix,lr_mat,K,M,alpha,truncation_value = 3){
  #The ligand expressing cell should be the niche cell
  niche = which(colnames(object@num_cells)==ligand_cell)
  #The receptor expressing cell should be the index cell
  index = which(colnames(object@num_cells)==receptor_cell)
  
  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  #get best kernel for each gene
  top_kernel = apply(log_liks,1,function(x){order(x,decreasing = T)[1]})
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
  
  #get ligand potential scores
  pear_cor = rep(NA,ncol(ligand_target_matrix))
  score_norm = rep(NA,ncol(ligand_target_matrix))
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
      #get weights based on scaled niche-net scores
      weight = ligand_vec[top_cors,j]/mean(ligand_vec[top_cors,j])
      ##calculate ligand potential scores
      pear_cor[j] = sum(sig[top_cors]*weight)
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
      print(paste0("error",j))
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
  LR_pairs = LR_pairs[which(as.numeric(LR_pairs[,3])<alpha),c(1:2)]
  return(LR_pairs)
  
}
#' get_niche_DE_genes
#'
#' This function returns genes that show niche patterns at the desired resolution
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
#' ligand potential score
#' @param M The maximum number of ligands that can pass initial filtering
#' @param alpha The null quantile to compare observed epression to
#' @param alpha_2 The level at which to perform the Benjamini Hochberg correction
#' @param truncation_value The value at which to truncate T statistics
#' @return A list of ligand-receptor pairs that are found to be expressed by the
#' specified cell type
#' @export
niche_LR_cell = function(object,ligand_cell,receptor_cell,ligand_target_matrix,
                         lr_mat,K,M,alpha,alpha_2,truncation_value = 3){
  niche = which(colnames(object@num_cells)==ligand_cell)
  index = which(colnames(object@num_cells)==receptor_cell)
  
  log_liks = object@niche_DE[[1]]$log_lik
  if(length(object@niche_DE)>=2){
    for(j in c(2:length(object@niche_DE))){
      log_liks = cbind(log_liks,object@niche_DE[[j]]$log_lik)
    }
  }
  #get best kernel for each gene
  top_kernel = apply(log_liks,1,function(x){order(x,decreasing = T)[1]})
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
  
  #get ligand potential scores
  pear_cor = rep(NA,ncol(ligand_target_matrix))
  score_norm = rep(NA,ncol(ligand_target_matrix))
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
      lambda_rest = quantile(object@ref_expr[niche,],alpha)
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
  ligands = top_genes[which(pvalues<alpha_2)]
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
      lambda_rest = quantile(object@ref_expr[index,],alpha)
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
  LR_pairs = LR_pairs[which(as.numeric(LR_pairs[,3])<alpha_2),c(1:2)]
  return(LR_pairs)
  
}               
                         
                      