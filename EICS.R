#' library(EPIC)
#' load("D:/guoman/cibersort/EPIC_SVR/data/BRef_guman_533.Rdata")
#' mixture <- read.table("D:/guoman/cibersort/EPIC_SVR/data/example2-百分-mix2.txt",header=T,sep="\t",row.names=1,check.names=F)
#' mixture <- data.matrix(mixture)
#' res <- EICS(bulk=mixture, reference=BRef_guoman_533)#较差
#' res <- EICS(bulk=mixture, reference=BRef.rda)#reference取EPIC包中data数据BRef.rda
library(e1071)
library(parallel)
library(preprocessCore)

CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)     #radial##支持向量机，其中参数nu的值用来改变分类算法的准确性，在本算法中，这三个数值都要进行计算
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') 
    out <- mclapply(1:svn_itor, res, mc.cores=1) else ###mclapply是并行计算函数，在windows系统下无法开启多核并行计算，如果不是windows就可以
      out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)  ##### 
  ##### 这里输出的是三个参数下的支持向量模型
  nusvm <- rep(0,svn_itor) ###初始化。重复函数,第一个参数重复的内容，第二个参数重复的次数
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV #提取模型的coefs和SV相乘作为权重
    weights[which(weights<0)]<-0  ###把小于0的权重赋值为0
    w<-weights/sum(weights)  ###
    u <- sweep(X,MARGIN=2,w,'*') ###sweep函数对矩阵进行计算，MARGIN=2代表列，MARGIN=1代表行，此处为在x矩阵按列乘以W
    k <- apply(u, 1, sum) ##对u按行求和
    nusvm[t] <- sqrt((mean((k - y)^2))) #平均值：mean(x)
    corrv[t] <- cor(k, y)#相关系数
    t <- t + 1
  }  #最后得到 nusvm和corrv，里面分别对应三个模型下的数值
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]] ###以最小的nusvm值挑选对应的模型
  
  #get and normalize coefficients #获得并且归一化系数
  q <- t(model$coefs) %*% model$SV  ##提取最优模型的coefs和SV相乘
  q[which(q<0)]<-0
  w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#' @keywords internal
scaleCounts <- function(counts, sigGenes=NULL, renormGenes=NULL, normFact=NULL){
  if (is.null(sigGenes))
    sigGenes <- 1:nrow(counts)
  
  if (is.null(normFact)){
    if (is.null(renormGenes))
      renormGenes <- 1:nrow(counts)
    normFact <- colSums(counts[renormGenes,,drop=FALSE], na.rm=TRUE)
  }
  counts <- t( t(counts[sigGenes,,drop=FALSE]) / normFact) * 1e6
  # Need to take the transpose so that the division is made on the correct
  # elements
  return(list(counts=counts, normFact=normFact))
}


EICS <- function(bulk, reference=NULL, mRNA_cell=NULL, mRNA_cell_sub=NULL,
                 sigGenes=NULL, scaleExprs=TRUE, withOtherCells=TRUE,
                 constrainedSum=TRUE, rangeBasedOptim=FALSE){
  # First get the value of the reference profiles depending on the input
  # 'reference'.
  with_w <- TRUE
  if (is.null(reference)){
    reference <- EPIC::TRef
  } else if (is.character(reference)){
    if (reference %in% prebuiltRefNames){
      reference <- get(reference, pos="package:EPIC")
      # Replace the char defining the reference name by the corresponding
      # pre-built reference values.
    } else
      stop("The reference, '", reference, "' is not part of the allowed ",
           "references:", paste(prebuiltRefNames, collapse=", "))
  } else if (is.list(reference)){
    refListNames <- names(reference)
    if ( (!all(c("refProfiles", "sigGenes") %in% refListNames)) ||
         (("refProfiles" %in% refListNames) && !is.null(sigGenes)) )
      stop("Reference, when given as a list needs to contain at least the ",
           "fields 'refProfiles' and 'sigGenes' (sigGenes could also be ",
           "given as input to EPIC instead)")
    if (!("refProfiles.var" %in% refListNames)){
      warning("'refProfiles.var' not defined; using identical weights ",
              "for all genes")
      reference$refProfiles.var <- 0
      with_w <- FALSE
    }
    if ((length(reference$refProfiles.var) > 1) &&
        (!identical(dim(reference$refProfiles.var), dim(reference$refProfiles))
         || !identical(dimnames(reference$refProfiles.var),
                       dimnames(reference$refProfiles))))
      stop("The dimensions and dimnames of 'reference$refProfiles' and ",
           "'reference$refProfiles.var' need to be the same")
  } else {
    stop("Unknown format for 'reference'")
  }
  
  refProfiles <- reference$refProfiles
  refProfiles.var <- reference$refProfiles.var
  
  nSamples <- NCOL(bulk); samplesNames <- colnames(bulk)
  if (is.null(samplesNames)){
    samplesNames <- 1:nSamples
    colnames(bulk) <- samplesNames
  }
  nRefCells <- NCOL(refProfiles); refCellsNames <- colnames(refProfiles)
  
  # Checking the correct format of the input variables
  if (!is.matrix(bulk) && !is.data.frame(bulk))
    stop("'bulk' needs to be given as a matrix or data.frame")
  if (!is.matrix(refProfiles) && !is.data.frame(refProfiles))
    stop("'reference$refProfiles' needs to be given as a matrix or data.frame")
  if (with_w && (!is.matrix(refProfiles.var) && !is.data.frame(refProfiles.var)))
    stop("'reference$refProfiles.var' needs to be given as a matrix or ",
         "data.frame when present.")
  
  # Keeping only common genes and normalizing the counts based on these common
  # genes
  bulk_NA <- apply(is.na(bulk), MARGIN=1, FUN=all)
  if (any(bulk_NA)){
    warning(sum(bulk_NA), " genes are NA in all bulk samples, removing these.")
    bulk <- bulk[!bulk_NA,]
  }
  bulkGenes <- rownames(bulk)
  if (anyDuplicated(bulkGenes))
    stop("There are some duplicated gene names in 'bulk'")
  refGenes <- rownames(refProfiles)
  if (anyDuplicated(refGenes))
    stop("There are some duplicated gene names in 'refGenes'")
  commonGenes <- intersect(bulkGenes, refGenes)
  
  if (is.null(sigGenes))
    sigGenes <- reference$sigGenes
  sigGenes <- sigGenes[sigGenes %in% commonGenes]
  nSigGenes <- length(sigGenes)
  if (nSigGenes < nRefCells)
    stop("There are only ", nSigGenes, " signature genes",
         " matching common genes between bulk and reference profiles,",
         " but there should be more signature genes than reference cells")
  
  if (scaleExprs){
    if (length(commonGenes) < 2e3)
      warning("there are few genes in common between the bulk samples and ",
              "reference cells:", length(commonGenes), ", so the data scaling ",
              "might be an issue")
    # The value of 2e3 is arbitrary, but should be a respectable number for the
    # data renormalization.
    bulk <- scaleCounts(bulk, sigGenes, commonGenes)$counts
    temp <- scaleCounts(refProfiles, sigGenes, commonGenes)
    refProfiles <- temp$counts
    if (with_w)
      refProfiles.var <- scaleCounts(refProfiles.var, sigGenes,
                                     normFact=temp$normFact)$counts
    # the refProfiles.var is normalized by the same factors as refProfiles.
  } else {
    bulk <- bulk[sigGenes,]
    refProfiles <- refProfiles[sigGenes,]
    if (with_w)
      refProfiles.var <- refProfiles.var[sigGenes,]
  }
  
  if (is.null(mRNA_cell))
    mRNA_cell <- EPIC::mRNA_cell_default
  
  if (!is.null(mRNA_cell_sub)){
    if (is.null(names(mRNA_cell_sub)) || !is.numeric(mRNA_cell_sub))
      stop("When mRNA_cell_sub is given, it needs to be a named numeric vector")
    mRNA_cell[names(mRNA_cell_sub)] <- mRNA_cell_sub
  }
  
  
  if (with_w && !rangeBasedOptim){
    # Computing the weight to give to each gene
    w <- rowSums(refProfiles / (refProfiles.var + 1e-12), na.rm=TRUE)
    # 1e-12 to avoid divisions by 0: like this if refProfiles and refProfiles.var
    # both are 0 for a given element, it will result in a weight of 0.
    med_w <- stats::median(w[w>0], na.rm=TRUE)
    w[w> 100*med_w] <- 100*med_w
    # Set a limit for the big w to still not give too much weights to these genes
  } else
    w <- 1
  
  # Defining the constraints for the fit of the proportions.
  if (withOtherCells){
    cMin <- 0
  } else {
    cMin <- 0.99
    # So that when we assume no cells without known reference profile are
    # present, we ask the sum of the mRNA proportions of known cells to be at
    # least 0.99.
  }
  cMax <- 1
  ui <- diag(nRefCells)
  ci <- rep(0,nRefCells)
  if (constrainedSum){
    ui <- rbind(ui, rep(1, nRefCells), rep(-1, nRefCells))
    ci <- c(ci, cMin, -cMax)
  }
  # ui and ci define the constraints, in the form "ui %*% x - ci >= 0".
  # The first constraints are that the proportion of each cell must be
  # positive; then if constrainedSum is true, we want the sum of all proportions
  # to be bigger than cMin and this sum must also be smaller or equal to cMax.
  cInitProp <- (min(1, cMax)-1e-5) / nRefCells
  # used as an initial guess for the optimized - start with all cell fractions
  # equal. We added the -1e-5 in cInitProp because the optimizer needs to have
  # the initial guess inside the admissible region and not on its boundary
  
  # Estimating for each sample the proportion of the mRNA per cell type.
  
  absolute=TRUE
  abs_method='sig.score'
  output <- matrix()
  b <- bulk[,1]
  result <- CoreAlg(refProfiles, b, absolute, abs_method)
  output<-result$w
  otherCells=1-sum(output)
  output <- cbind(output,otherCells)
  
  for(cSample in 2:nSamples){
    b <- bulk[,cSample]
    print(cSample)
    result <- CoreAlg(refProfiles, b, absolute, abs_method)
    #get results
    w <- result$w
    otherCells=1-sum(w)
    out <- cbind(w,otherCells)
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    output <- rbind(output, out)
  }
  
  #save results
  samplename <- colnames(bulk)
  out0 <- cbind(samplename,output)
  header <- c('Mixture',colnames(refProfiles),"otherCells")
  write.table(rbind(header,out0), file="Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- output
  rownames(obj) <- colnames(bulk)
  obj
  
}


