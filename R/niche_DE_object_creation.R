library(Matrix)
setClassUnion(name = 'AnyMatrix', members = c("matrix",'data.frame','dgCMatrix'))
#set niche-de
Assay <- setClass(
  Class = 'Niche_DE',
  slots = c(
    counts = 'AnyMatrix',
    coord = 'AnyMatrix',
    sigma = 'vector',
    num_cells = 'AnyMatrix',
    effective_niche = 'list',
    ref_expr = 'AnyMatrix',#null_expected_expression = 'AnyMatrix',
    cell_names = 'vector',
    cell_types = 'vector',
    gene_names = 'vector',
    batch_ID = 'vector',
    spot_distance = 'numeric',
    niche_DE = 'list',
    niche_DE_pval_pos = 'list',
    niche_DE_pval_neg = 'list',
    scale = 'vector',
    Int = 'logical'
  )
)

#' @export
print.Niche_DE = function(object){
  A = paste0('Niche-DE object with ',nrow(object@counts),' observations, ', ncol(object@counts),' genes, ',
         length(unique(object@batch_ID)), ' batch(es), and ', length(object@cell_types), ' cell types.')
  return(A)
}

#' CreateLibraryMatrix
#'
#' This function creates a expression profile matrix for single cell data
#'
#' @param data Single cell rna-seq counts matrix. Dimension should be cells/spots by genes
#' @param cell_type Cell_type assignment matrix. First column is cell names and second column is cell type assignment.
#' @return A library matrix with the average expression profile per cell type
#' @export
CreateLibraryMatrix = function(data,cell_type){
  if( mean(rownames(data)==cell_type[,1])!=1){
    stop('Data rownames and Cell type matrix names do not match')
  }
  print('Computing average expression profile matrix')
  #get unique cell types
  CT = unique(as.vector(cell_type[,2]))
  n_CT = length(CT)
  L = matrix(NA,n_CT,ncol(data))
  rownames(L) = CT
  colnames(L) = colnames(data)
  #iterate over cell types
  for (j in c(1:n_CT)){
    #get cells that belong to this cell type
    cells = which(cell_type[,2]==CT[j])
    #if there are too many cells, downsample
    if(length(cells)>1000){
      print(paste0("Too many cell of type ",CT[j]," downsampling to 1000."))
      cells = sample(cells,1000,replace = F)
    }
    cells = data[cells,]
    L[j,] = apply(cells,2,function(x){mean(x)})
  }
  print('Average expression matrix computed')
  return(L)
}

#' CreateLibraryMatrixFromSeurat
#'
#' This function creates a expression profile matrix for single cell data
#' from a Seurat object
#'
#' @param seurat_object A seurat object
#' @param assay The assay from which we want to extract the counts matrix to
#' calculate the average expression profile.
#' @return A library matrix with the average expression profile per cell type
#' @export
CreateLibraryMatrixFromSeurat = function(seurat_object,assay){
  #get desired assay from seurat object
  sobj_assay = Seurat::GetAssay(seurat_object,assay)
  #get counts matrix
  data = Matrix::t(sobj_assay@counts)
  #get cell type vector
  cell_type = Seurat::Idents(seurat_object)
  if(mean(rownames(data)== names(cell_type))!=1){
    stop('Data rownames and Cell type matrix names do not match')
  }
  print('Computing average expression profile matrix')
  #get unique cell types
  CT = unique(cell_type)
  n_CT = length(CT)
  L = matrix(NA,n_CT,ncol(data))
  rownames(L) = CT
  colnames(L) = colnames(data)
  #iterate over cell types
  for (j in c(1:n_CT)){
    #get cells that belong to this cell type
    cells = which(cell_type==CT[j])
    #if there are too many cells, downsample
    if(length(cells)>1000){
      print(paste0("Too many cell of type ",CT[j],", downsampling to 1000."))
      cells = sample(cells,1000,replace = F)
    }
    cells = data[cells,]
    L[j,] = apply(cells,2,function(x){mean(x)})
  }
  print('Average expression profile matrix computed.')
  return(L)
}


#' CreateNicheDEObject
#'
#' This function creates a niche-DE object
#'
#' @param counts_mat Counts matrix. Dimension should be cells/spots by genes
#' @param coordinate_mat Coordinate matrix. It will be scaled such that the median nearest neighbor distance is 100.
#' @param library_mat Matrix indicating average expression profile for each cell type in the sample
#' @param deconv_mat Deconvolution or cell type assignment matrix of data
#' @param sigma List of kernel bandwidths to use in calculating the effective niche
#' @param Int Boolean of if counts data supplied is integer. Default is true. When performing niche-DE,
#'  Negative binomial regression is performed if True. Linear regression with a gene specific variance is performed if False.
#' @return A niche-DE object
#' @export
CreateNicheDEObject = function(counts_mat,coordinate_mat,library_mat,deconv_mat,sigma,Int = T){
  print('Creating Niche-DE object')
  #make sure that counts matrix is provided
  if (missing(x = counts_mat)) {
    stop("Must provide counts matrix")
  }

  if (Int == T & sum(counts_mat%%1)!=0){
    stop('counts matrix must contain only integers')
  }

  #make sure that cell names (rownames) are not null
  if (is.null(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be non-null')
  }
  if (is.null(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be non-null')
  }
  #make sure that cell names are unique
  if (anyDuplicated(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be unique')
  }
  #make sure that gene names are unique
  if (anyDuplicated(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be unique')
  }


  #make sure that counts_mat and coordinate_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(coordinate_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and coordinate matrix do not match')
  }

  #make sure that counts_mat and deconv_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(deconv_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and deconvolution matrix do not match')
  }

  #make sure that deconv_mat and library_mat have the same cell types in the same order
  if(mean(colnames(deconv_mat)==rownames(library_mat))!=1){
    stop('celltypes of deconvolution matrix and reference expression matrix do not match')
  }

  #make sure that sigma is a vector
  if(length(sigma)>0){
    if((is.vector(sigma) && is.atomic(sigma))==F){
      stop('Sigma must be a vector')
    }
    #make sure that sigma is numeric
    if(is.numeric(sigma)==F){
      stop('sigma must be numeric')
    }
  }
  #make sure that sigma has positive length
  if(length(sigma)==0){
    warning('No sigma(kernel bandwidth) values selected. Default values will be used.
            These default values are only appropriate for data that has a similar resolution to 10X VISIUM
            (55 micrometers in diameter).')
  }

  #calculate number of cells per spot and expected expression
  #get genes that are shared between data and reference expression
  sim_gene = which(colnames(library_mat) %in% colnames(counts_mat))
  #get the gene list
  gene_list = colnames(library_mat)[sim_gene]
  #filter library_mat by removing genes that are not in the counts_mat
  library_mat = library_mat[,sim_gene]
  #get rowsums of the library_matrix (i.e expected library size of a cell type)
  LM = rowSums(library_mat)
  #Get library size of spots
  #get genes that are in the library_mat
  sim_gene = which(colnames(counts_mat)%in% gene_list)
  #filter counts_mat by removing genes that are not in the library_mat
  countsM = counts_mat[,sim_gene]
  #change column and row names
  colnames(countsM) = colnames(counts_mat)[sim_gene]
  rownames(countsM) = rownames(counts_mat)
  #get library size of each spot
  #Lib_spot = rowSums(as.matrix(countsM))
  Lib_spot = apply(countsM,1,function(x){sum(x)})
  #make matrix sparse
  countsM = Matrix::Matrix(counts_mat[,sim_gene], sparse=TRUE)

  #Get expected library size given pi(deconvolution estimate for each spot)
  EL = deconv_mat%*%as.matrix(LM)
  #get expected number of total cells in a spot
  num_cell = Lib_spot/EL
  #get effective niche
  #nst = diag(num_cell[,1])%*%as.matrix(deconv_mat)
  nst = sweep(as.matrix(deconv_mat), MARGIN=1, num_cell[,1], `*`)
  rownames(nst) = rownames(countsM)
  #get expected gene expression given pi
  #EEX = as.matrix(nst)%*%as.matrix(library_mat)
  #rownames(EEX) = rownames(countsM)
  #reorder columns of count data to match that of library reference
  col.order = colnames(library_mat)
  countsM = countsM[,col.order]

  #get min spot distance
  if(dim(coordinate_mat)[1] < 1e4){
    D = as.matrix(dist(coordinate_mat),diag = T)
    min_dist = mean(apply(D,2,function(x){sort(x,decreasing = F)[3]}))
  }else{
    #sample 1000 random points
    inds = sample(c(1:nrow(coordinate_mat)),1000)
    #get distances
    D = Rfast::dista(xnew = as.matrix(coordinate_mat[inds,],ncol = 2,byrow = T),x = coordinate_mat,type = "euclidean",trans = T)
    #get min distance
    mdist = apply(D,1,function(x){sort(x[x>0],decreasing = F)[3]})
    #get median min_dist
    min_dist = median(mdist)
  }

  #scale coordiante matrix so that min_dist = 100
  scale = 100/min_dist
  coordinate_mat = coordinate_mat*scale
  min_dist = 100


  if(length(sigma)==0){
    sigma = c(min_dist*0.001,min_dist,min_dist*2,min_dist*3)
  }

  #make sure that counts_mat and library_mat have the same gene names in the same order
  if(mean(colnames(countsM)==colnames(library_mat))!=1){
    stop('gene names (colnames) of counts matrix and library expression matrix do not match')
  }

  object = new(Class = 'Niche_DE',counts = countsM,coord = coordinate_mat,
               sigma = sigma,num_cells = nst,ref_expr = library_mat,#null_expected_expression = EEX,
               cell_names = rownames(countsM), cell_types = colnames(deconv_mat),
               gene_names = colnames(countsM),batch_ID = rep(1,nrow(countsM)),
               spot_distance = min_dist,scale = scale,Int = Int)
  A = paste0('Niche-DE object created with ',nrow(object@counts),' observations, ', ncol(object@counts),' genes, ',
             length(unique(object@batch_ID)), ' batch(es), and ', length(object@cell_types), ' cell types.')
  print(A)
  return(object)
}

#' CreateNicheDEObjectFromSeurat
#'
#' This function creates a niche-DE object from a seurat object
#'
#' @param seurat_object A spatial seurat object.Coordinate matrix will be extracted via the
#' seurat function 'GetTissueCoordinates'. The coordiantes will be scaled such that the median nearest neighbor distance is 100.
#' @param assay The assay from which to extract the counts matrix from. The counts matrix
#' will be extracted from the counts slot.
#' @param library_mat Matrix indicating average expression profile for each cell type in the sample
#' @param deconv_mat Deconvolution or cell type assignment matrix of data
#' @param sigma List of kernel bandwidths to use in calculating the effective niche
#' @param Int Boolean of if counts data supplied is integer. Default is true. When performing niche-DE,
#'  Negative binomial regression is performed if True. Linear regression with a gene specific variance is performed if False.
#' @return A niche-DE object
#' @export
CreateNicheDEObjectFromSeurat = function(seurat_object,assay,library_mat,deconv_mat,sigma, Int = T){
  print("Creating Niche-DE object")
  #make sure that counts matrix is provided
  if (missing(x = seurat_object)) {
    stop("Must provide seurat object matrix")
  }
  #extract raw counts matrix from seurat object
  sobj_assay = Seurat::GetAssay(seurat_object,assay)
  counts_mat = Matrix::t(sobj_assay@counts)
  #make sure that counts_mat is integers
  if (Int ==T & sum(counts_mat%%1)!=0){
    stop('counts matrix must contain only integers')
  }

  #extract coordinate matrix from seurat object
  #coordinate_mat = Seurat::GetTissueCoordinates(seurat_object,image = names(seurat_object@images)[1])
  #slice = names(seurat_object@images)[1]
  x = as.numeric(seurat_object@images[[1]]@coordinates$imagerow)
  y = as.numeric(seurat_object@images[[1]]@coordinates$imagecol)
  coordinate_mat = cbind(x,y)
  rownames(coordinate_mat) = rownames(seurat_object@images[[1]]@coordinates)
  #make sure that cell names (rownames) are not null
  if (is.null(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be non-null')
  }
  if (is.null(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be non-null')
  }
  #make sure that cell names are unique
  if (anyDuplicated(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be unique')
  }
  #make sure that gene names are unique
  if (anyDuplicated(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be unique')
  }


  #make sure that counts_mat and coordinate_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(coordinate_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and coordinate matrix do not match')
  }

  #make sure that counts_mat and deconv_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(deconv_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and deconvolution matrix do not match')
  }

  #make sure that deconv_mat and library_mat have the same cell types in the same order
  if(mean(colnames(deconv_mat)==rownames(library_mat))!=1){
    stop('celltypes of deconvolution matrix and reference expression matrix do not match')
  }



  #make sure that sigma is a vector
  if(length(sigma)>0){
    if((is.vector(sigma) && is.atomic(sigma))==F){
      stop('Sigma must be a vector')
    }
    #make sure that sigma is numeric
    if(is.numeric(sigma)==F){
      stop('sigma must be numeric')
    }
  }
  #make sure that sigma has positive length
  if(length(sigma)==0){
    warning('No sigma(kernel bandwidth) values selected. Default values will be used.
            These default values are only appropriate for data that has a similar resolution to 10X VISIUM
            (55 micrometers in diameter).')
  }


  #calculate number of cells per spot and expected expression
  #get genes that are shared between data and reference expression
  sim_gene = which(colnames(library_mat) %in% colnames(counts_mat))
  #get the gene list
  gene_list = colnames(library_mat)[sim_gene]
  #filter library_mat by removing genes that are not in the counts_mat
  library_mat = library_mat[,sim_gene]
  #get rowsums of the library_matrix (i.e expected library size of a cell type)
  LM = rowSums(library_mat)
  #Get library size of spots
  #get genes that are in the library_mat
  sim_gene = which(colnames(counts_mat)%in% gene_list)
  #filter counts_mat by removing genes that are not in the library_mat
  countsM = as.matrix(counts_mat[,sim_gene])
  colnames(countsM) = colnames(counts_mat)[sim_gene]
  rownames(countsM) = rownames(counts_mat)
  #get library size of each spot
  Lib_spot = rowSums(countsM)
  #make matrix sparse
  countsM = Matrix::Matrix(counts_mat[,sim_gene], sparse=TRUE)

  #Get expected library size given pi(deconvolution estimate for each spot)
  EL = deconv_mat%*%as.matrix(LM)
  #get expected number of total cells in a spot
  num_cell = Lib_spot/EL
  #get effective niche
  nst = diag(num_cell[,1])%*%as.matrix(deconv_mat)
  rownames(nst) = rownames(countsM)
  #get expected gene expression given pi
  #EEX = as.matrix(nst)%*%as.matrix(library_mat)
  #rownames(EEX) = rownames(countsM)
  #reorder columns of count data to match that of library reference
  col.order = colnames(library_mat)
  countsM = countsM[,col.order]

  #get min spot distance
  D = as.matrix(dist(coordinate_mat),diag = T)
  min_dist = mean(apply(D,2,function(x){sort(x,decreasing = F)[3]}))

  #scale coordiante matrix so that min_dist = 100
  scale = 100/min_dist
  coordinate_mat = coordinate_mat*scale
  min_dist = 100

  #make sure that counts_mat and library_mat have the same gene names in the same order
  if(mean(colnames(countsM)==colnames(library_mat))!=1){
    stop('gene names (colnames) of counts matrix and library expression matrix do not match')
  }

  object = new(Class = 'Niche_DE',counts = countsM,coord = coordinate_mat,
               sigma = sigma,num_cells = nst,ref_expr = library_mat,#null_expected_expression = EEX,
               cell_names = rownames(countsM),cell_types = colnames(deconv_mat),
               gene_names = colnames(countsM),batch_ID = rep(1,nrow(countsM)),
               spot_distance = min_dist,scale = scale,Int = Int)
  #make sure that counts_mat and

  A = paste0('Niche-DE object created with ',nrow(object@counts),' observations, ', ncol(object@counts),' genes, ',
             length(unique(object@batch_ID)), ' batch(es), and ', length(object@cell_types), ' cell types.')
  print(A)
  return(object)
}


#' MergeObjects
#'
#' This function merges niche-DE objects
#'
#' @param objects A list of niche-DE objects
#' @return A niche-DE object with each niche-DE object that concatenates
#' each niche-DE object used in the input list. Coordinates are scaled
#' to ensure that kernel bandwidths are consistent across datasets
#' @export
MergeObjects = function(objects){
  #main idea: iterate through batches and merge them
  if(is.list(objects)==F){
    stop('Input must be a list of niche_DE objects')
  }
  #get reference object
  reference_obj = objects[[1]]
  coord_merge = reference_obj@coord
  counts_merge = reference_obj@counts
  num_cells_merge = reference_obj@num_cells
  #refr_spot_distance = reference_obj@min_distance
  refr_spot_distance = 100/reference_obj@scale
  batch_ID_merge = reference_obj@batch_ID
  scales = c(1)
  for(j in c(2:length(objects))){
    #check if sigma and L match
    if(setequal(objects[[j]]@sigma,reference_obj@sigma)==F){
      stop('all objects being merged must consider the same kernel bandwidths')
    }
    if(setequal(objects[[j]]@Int,reference_obj@Int)==F){
      stop('all counts must all be integer valued(object@Int = T) or all be continuous(object@Int = F)')
    }
    if(mean(objects[[j]]@ref_expr==reference_obj@ref_expr)==F){
      stop('all objects being merged must consider the same kernel bandwidths')
    }


    #merge coordinates so that they scale the same
    #scale = (objects[[j]]@min_distance)/refr_spot_distance
    scale = (100/objects[[j]]@scale)/refr_spot_distance
    scales = c(scales,scale)
    coord_merge = rbind(coord_merge,objects[[j]]@coord*1/scale)
    #merge num cells
    if(mean(colnames(num_cells_merge) == colnames(objects[[j]]@num_cells))!=1){
      stop('cell types must be the same in merged objects')
    }
    num_cells_merge = rbind(num_cells_merge,objects[[j]]@num_cells)

    #merge actual expresssion
    if(mean(colnames(counts_merge) == colnames(objects[[j]]@counts))!=1){
      stop('genes must be the same and in the same order in merged objects')
    }
    counts_merge = rbind(counts_merge,objects[[j]]@counts)



    #merge batch ID
    batch_ID_new = objects[[j]]@batch_ID + max(batch_ID_merge)
    batch_ID_merge = c(batch_ID_merge,batch_ID_new)

    #rename cells/spots
    rownames(counts_merge) = c(1:nrow(counts_merge))
    rownames(coord_merge) = c(1:nrow(coord_merge))
    rownames(num_cells_merge) = c(1:nrow(num_cells_merge))

    cell_names = rownames(counts_merge)
    gene_names = colnames(counts_merge)

  }
  #make merged object
  object = new(Class = 'Niche_DE',counts = counts_merge,coord = coord_merge,
               sigma = reference_obj@sigma,num_cells = num_cells_merge,
               ref_expr = reference_obj@ref_expr,
               cell_names = cell_names, cell_types = colnames(num_cells_merge),
               gene_names = gene_names,batch_ID = batch_ID_merge,
               spot_distance = 100,scale = scales,Int = reference_obj@Int)
  return(object)
}



#' CalculateEffectiveNiche
#'
#' This function calculates the effective niche of a niche-DE object
#'
#' @param object A niche-DE object
#' @param batch_size Number of cells to calculate the effective niche for at a time
#' @param cutoff Minimum kernel similarity. Similarities below this value get truncated to 0
#' @param standardize A logical indicator for if the effective niche should be z-normalized by column
#' @return A niche-DE object the effective niche calculated.
#' The effective niche is a list with each entry corresponding to a kernel bandwidth
#' @export
CalculateEffectiveNicheLargeScale = function(object,batch_size = 1000,cutoff = 0.05,standardize = TRUE){
  EN = matrix(NA,nrow(object@coord),length(object@cell_types))
  rownames(EN) = rownames(object@coord)
  colnames(EN) = object@cell_types
  #initialize effective niche list
  object@effective_niche = vector(mode = "list", length = length(object@sigma))
  #name each item of the list after the kernel bandwidth being used
  names(object@effective_niche) = object@sigma
  #ge counter for kernel bandwidth
  counter_sig = 1
  for(sig in object@sigma){
    print(paste0('Calculating effective niche for kernel bandwith ', sig,
                 '(',counter_sig,' out of ',length(object@sigma),' values).'))
    #calculate effective niche for each individual dataset
    counter = 0
    for(ID in c(1:length(unique(object@batch_ID)))){
      print(paste0("Calculating effective niche for batch ",ID))
      #get coordinates for dataset
      coord_ID = object@coord[which(object@batch_ID == ID),]
      #initialize effective niche dataset
      EN_dataset = matrix(NA,nrow(coord_ID),length(object@cell_types))
      #get number of cells per spot for dataset
      num_cell_ID = object@num_cells[object@batch_ID == ID,]
      #get number of cells per batch
      C = nrow(coord_ID)
      num_iter = ceiling(C/batch_size)

      print("Separating Tissue Into Spatially Local Regions")
      square_side = ceiling(sqrt(num_iter))

      #aggregate all cells that belong to same radius
      x_range = range(coord_ID[,1])
      y_range = range(coord_ID[,2])

      x_length = x_range[2] - x_range[1]
      y_length = y_range[2] - y_range[1]

      x_bins = seq(x_range[1],x_range[2]  + x_length/(2*square_side),length.out = square_side)
      y_bins = seq(y_range[1],y_range[2] + y_length/(2*square_side),length.out = square_side)
      #x_bins = seq(x_range[1],x_range[2],length.out = square_side)
      #y_bins = seq(y_range[1],y_range[2],length.out = square_side)

      x_bin_index <- findInterval(coord_ID[,1], x_bins)
      y_bin_index <- findInterval(coord_ID[,2], y_bins)
      cell_counter = 0
      for (j in c(1:square_side)) {
        for(k in c(1:square_side)){
          inds = which(x_bin_index == j & y_bin_index == k)
          Ncell = length(inds)
          if(Ncell == 0){
            next
          }
          print(paste0("Calculating effective niche for cells ",
                       cell_counter + 1, " to ", cell_counter + Ncell))

          x_range = range(coord_ID[inds, 1])
          y_range = range(coord_ID[inds, 2])
          x_range = c(x_range[1] - sig * sqrt(-log(cutoff)),
                      x_range[2] + sig * sqrt(-log(cutoff)))
          y_range = c(y_range[1] - sig * sqrt(-log(cutoff)),
                      y_range[2] + sig * sqrt(-log(cutoff)))
          cand = which(spatstat.utils::inside.range(coord_ID[,
                                                             1], x_range) & spatstat.utils::inside.range(coord_ID[,
                                                                                                                  2], y_range))
          D = Rfast::dista(xnew = coord_ID[inds,
                                           ,drop = FALSE], x = coord_ID[cand,
                                                                        ,drop = FALSE], type = "euclidean", trans = T)
          D = exp(-D^2/sig^2)
          D[D<cutoff] = 0
          EN_dataset[inds, ] = D %*% num_cell_ID[cand,
          ]
          cell_counter = cell_counter + Ncell
        }
      }
      #bind EN of this dataset to EN of other datasets
      if(counter == 0){
        EN[which(object@batch_ID == ID),] = EN_dataset
        ref_size = mean(rowSums(num_cell_ID))
      } else{
        #scale EN
        scale = ref_size/mean(rowSums(num_cell_ID))
        EN[which(object@batch_ID == ID),] = EN_dataset*scale
      }
      #make counter bigger
      counter = counter + 1
    }
    #normalize columns of EN
    if(standardize == T){
      EN = apply(EN,2,function(x){(x-mean(x,na.rm = T))/sd(x,na.rm =T)})
    }
    EN[is.na(EN)] = 0
    #add to list of effective niches (one for each sigma)
    object@effective_niche[[counter_sig]] = EN
    counter_sig = counter_sig + 1
  }
  print('Effective niche calculated')
  return(object)
}








#' Filter_NDE
#'
#' This function filters a niche-DE object to only include specific observations
#'
#' @param object A niche-DE object
#' @param cell_names Cell names of observations that should be kept
#' @return A niche-DE object that only includes the specified cells.
#' We recommend using this function after calculating the effective niche on the whole dataset
#' @export
Filter_NDE = function(object,cell_names){
  if(mean(cell_names%in%object@cell_names)!=1){
    stop('Some cell names not present in niche-DE object')
  }
  object@counts = object@counts[rownames(object@counts) %in% cell_names,]
  object@coord = object@coord[rownames(object@coord) %in% cell_names,]
  object@num_cells = object@num_cells[rownames(object@num_cells) %in% cell_names,]
  if(length(object@effective_niche) > 0){
    for(j in c(1:length(object@effective_niche))){
      object@effective_niche[[j]] = object@effective_niche[[j]][rownames(object@effective_niche[[j]]) %in% cell_names,]
    }
  }
  #object@null_expected_expression = object@null_expected_expression[rownames(object@null_expected_expression) %in% cell_names,]
  object@batch_ID = object@batch_ID[which(object@cell_names%in% cell_names)]
  object@cell_names = object@cell_names[which(object@cell_names%in% cell_names)]
  return(object)
}









#' CalculateEffectiveNiche
#'
#' This function calculates the effective niche of a niche-DE object
#'
#' @param object A niche-DE object
#' @param cutoff Minimum kernel similarity. Similarities below this value get truncated to 0
#' @return A niche-DE object the effective niche calculated.
#' The effective niche is a list with each entry corresponding to a kernel bandwidth
#' @export
CalculateEffectiveNiche = function(object,cutoff = 0.05){
  object@effective_niche = vector(mode = "list", length = length(object@sigma))
  names(object@effective_niche) = object@sigma
  counter_sig = 1
  for(sig in object@sigma){
    print(paste0('Calculating effective niche for kernel bandwith ', sig,
                 '(',counter_sig,' out of ',length(object@sigma),' values).'))
    #calculate effective niche for each individual dataset
    counter = 0
    for(ID in c(1:length(unique(object@batch_ID)))){
      #get coordinates for dataset
      coord_ID = object@coord[object@batch_ID == ID,]
      #get number of cells per spot for dataset
      num_cell_ID = object@num_cells[object@batch_ID == ID,]
      #calculate distance matrix for dataset
      K = as.matrix(dist(coord_ID,method = "euclidean",diag = TRUE))
      #transform to get kernel matrix for dataset
      K = exp(-K^2/sig^2)
      #truncate values less than cutoff to be 0
      K[K<cutoff] = 0
      #calculate effective niche
      EN_dataset = K%*%num_cell_ID
      #bind EN of this dataset to EN of other datasets
      if(counter == 0){
        EN = EN_dataset
        ref_size = mean(rowSums(num_cell_ID))
      } else{
        #scale EN
        scale = ref_size/mean(rowSums(num_cell_ID))
        EN = rbind(EN,EN_dataset*scale)
      }
      #make counter bigger
      counter = counter + 1
    }
    #normalize columns of EN
    EN = apply(EN,2,function(x){(x-mean(x))/sd(x)})
    EN[is.na(EN)] = 0
    #EN = apply(EN,2,function(x){x-mean(x)})
    #add to list of effective niches (one for each sigma)
    object@effective_niche[[counter_sig]] = EN
    counter_sig = counter_sig + 1
  }
  print('Effective niche calculated')
  return(object)
}



