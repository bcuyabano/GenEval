#' Simulate genotypes
#'
#' This function simulates numeric SNP-genotypes.
#'
#' @param n The number of individuals to simulate.
#' @param AlleleFreq A vector with the required allele frequencies for each SNP, with \code{length(AlleleFreq)} being the number of SNPs.
#' @param LD Logical indicating if SNP-genotypes should be simulated in linkage disequilibrium.
#' @param R The correlation matrix for the linkage disequilibrium between SNPs. If \code{R=NULL}, a default matrix with 0.95^(lag between SNPs) will be generated.
#' @param phased Lagical indicating if phased genotypes should be returned.
#' @return The simulated SNPs, as a single matrix with allele counts (\code{phased=FALSE}), or as a list of matrices with t0,1 at each allele (\code{phased=TRUE}).
#' @examples simGeno(n=10,AlleleFreq=rep(0.5,5))
#' @examples simGeno(n=10,AlleleFreq=rep(0.5,5),phased=TRUE)
#' @examples simGeno(n=10,AlleleFreq=rep(0.5,5),LD=TRUE)
#' @examples R_LD <- diag(1,5)
#' @examples R_LD[abs(col(R_LD)-row(R_LD)) == 1] <- 0.95
#' @examples R_LD[abs(col(R_LD)-row(R_LD)) == 2] <- 0.85
#' @examples R_LD[abs(col(R_LD)-row(R_LD)) == 3] <- 0.75
#' @examples R_LD[abs(col(R_LD)-row(R_LD)) == 4] <- 0.65
#' @examples simGeno(n=10,AlleleFreq=rep(0.5,5),LD=TRUE,R=R_LD)
#' @examples geno <- simGeno(n=1000,AlleleFreq=rep(0.5,5),LD=TRUE,R=R_LD,phased=TRUE)
#' @examples str(geno)
#' @examples cor(rbind(geno[[1]],geno[[2]]))
#' @export
simGeno <- function(n,AlleleFreq,LD=FALSE,R=NULL,phased=FALSE){
  # internal function
  mvrbinom <- function(p){
    return(rbinom(n,1,p))
  }
  # create list to store genotypes
  M <- vector("list",2)
  if(!LD){
    #--------------------------------------#
    # simulate SNPs in linkage equilibrium #
    #--------------------------------------#
    if(!is.null(R)){
      warning("Unused argument: 'R' (due to 'LD=FALSE').")
    }
    M[[1]] <- apply(matrix(AlleleFreq,nrow=1),2,mvrbinom)
    M[[2]] <- apply(matrix(AlleleFreq,nrow=1),2,mvrbinom)
  } else{
    #-----------------------------------------#
    # simulate SNPs in linkage disequilibrium #
    #-----------------------------------------#
    if(is.null(R)){
      # generate R if R=NULL
      tmp <- c(1,0.95^(1:(length(AlleleFreq)-1)))
      R <- matrix(NA,length(AlleleFreq),length(AlleleFreq))
      aux <- as.numeric(abs(row(R)-col(R)))+1
      R <- matrix(tmp[aux],length(AlleleFreq),length(AlleleFreq))
      gc()
    } else{
      # verify that R is correctly defined
      if(class(R)[1] == "matrix"){
        if(!isSymmetric(R)){
          stop("'R' must be a symmetric matrix.")
        }
        if(nrow(R) != length(AlleleFreq)){
          stop("Dimensions of 'R' must match the length of 'AlleleFreq'.")
        }
        if(any(diag(R) != 1)){
          stop("Elements in 'diag(R)' must be equal to 1.")
        }
        if(any(R < 0 | R > 1)){
          stop("Elements in 'R' must be within the (0,1) interval.")
        }
      } else stop("'R' must be a matrix.")
    }
    M[[1]] <- matrix(0,n,length(AlleleFreq))
    M[[1]][mvrnorm(n,qnorm(AlleleFreq),R) > 0] <- 1
    M[[2]] <- matrix(0,n,length(AlleleFreq))
    M[[2]][mvrnorm(n,qnorm(AlleleFreq),R) > 0] <- 1
  }
  #-------------#
  # return data #
  #-------------#
  if(!phased){
    M <- M[[1]] + M[[2]]
    gc()
  }
  return(M)
}

#' Simulate phenotypes
#'
#' This function simulates phenotypes (continuous or binary) based on the quantitative trait loci (QTL), only with the genetic values and residuals (\code{y=g+e}).
#'
#' @param M A matrix or a list of matrices of the (unphased) SNP-genotypes set to be the QTL.
#' @param h2 The heritability.
#' @param Vy The phenotypic variance (for continuous traits). Default is \code{Vy=1}.
#' @param binary Logical indicating whether simulated phenotypes should be binary. Default is \code{binary=FALSE}.
#' @param p If \code{binary=TRUE}, the probability of observing 1's in the phenotypes, with default \code{p=0.5}.
#' @param sim.error An error threshold to be considered on the orthogonalization of genetic values and residuals.
#' @return
#' \describe{
#' \item{y}{The simulated phenotypes.}
#' \item{a}{The simulated QTL effects.}
#' \item{g}{The simulated genetic values.}
#' \item{e}{The simulated residuals.}
#' \item{h2}{The simulated heritability.}
#' }
#' @examples ## general examples ##
#' @examples AF <- runif(500,0.05,1)
#' @examples geno <- simGeno(n=1000,AlleleFreq=AF,LD=TRUE)
#' @examples pheno <- simPheno(M=geno,h2=0.6,Vy=10)
#' @examples str(pheno)
#' @examples pheno <- simPheno(M=geno,h2=0.6,binary=TRUE)
#' @examples str(pheno)
#' @examples mean(pheno$y) # this is the realized probability of observing y = 1
#' @examples pheno <- simPheno(M=geno,h2=0.6,binary=TRUE,p=0.2)
#' @examples str(pheno)
#' @examples mean(pheno$y) # this is the realized probability of observing y = 1
#'
#' @examples ## example with multiple components with different heritabilities ##
#' @examples AF <- runif(500,0.05,1)
#' @examples geno <- simGeno(n=1000,AlleleFreq=AF,LD=TRUE)
#' @examples pheno <- simPheno(M=list(geno[,1:350],geno[,351:500]),h2=c(0.2,0.4),Vy=10)
#' @examples str(pheno)
#' @examples pheno <- simPheno(M=list(geno[,1:350],geno[,351:500]),h2=c(0.2,0.4),Vy=10,binary=TRUE,p=0.3)
#' @examples str(pheno)
#' @examples mean(pheno$y) # this is the realized probability of observing y = 1
#' @seealso \code{\link{simGeno}}
#' @export
simPheno <- function(M,h2,Vy=1,binary=FALSE,p=NULL,sim.error=1e-4){
  # center genotypes
  if(class(M)[1] == "matrix"){
    M <- scale(M,scale=FALSE)
  } else M <- lapply(M,scale,scale=FALSE)
  # simulate genetic effects/values
  if(class(M)[1] == "matrix"){
    # preliminary QTL effects
    a <- as.numeric(scale(rnorm(ncol(M),0,1)))
    # preliminary genetic values
    g <- as.numeric(M%*%a)
    # rescale QTL effects to achieve var(g) = h2
    a <- sqrt(h2*Vy/var(g))*a
    # recalculate genetic values
    g <- as.numeric(M%*%a)
    # start phenotypes with genetic values
    y <- g
  } else{
    a <- vector("list",length(M))
    g <- vector("list",length(M))
    for(i in 1:length(M)){
      # preliminary QTL effects
      a[[i]] <- as.numeric(scale(rnorm(ncol(M[[i]]),0,1)))
      # preliminary genetic values
      g[[i]] <- as.numeric(M[[i]]%*%a[[i]])
      if(i == 1){
        # rescale QTL effects to achieve var(g) = h2
        a[[i]] <- sqrt(h2[[i]]*Vy/var(g[[i]]))*a[[i]]
        # recalculate genetic values
        g[[i]] <- as.numeric(M[[i]]%*%a[[i]])
        # start phenotypes with genetic values
        y <- g[[i]]
      } else{
        for(j in 1:(i-1)){
          # orthogonalize genetic values
          g[[i]] <- sqrt(Vy*h2[[i]]/var(g[[i]]))*(g[[i]] - as.numeric(crossprod(g[[j]],g[[i]])/crossprod(g[[j]]))*g[[j]])
          # repeat orthogonalization until genetic variance is perfect
          while(abs(var(g[[i]]) - Vy*h2[[i]]) > sim.error){
            g[[i]] <- sqrt(Vy*h2[[i]]/var(g[[i]]))*(g[[i]] - as.numeric(crossprod(g[[j]],g[[i]])/crossprod(g[[j]]))*g[[j]])
            # rescale QTL effects to achieve var(g) = h2
            a[[i]] <- as.numeric(ginv(M[[i]])%*%g[[i]])
            a[[i]] <- sqrt(h2[[i]]*Vy/var(g[[i]]))*a[[i]]
            # recalculate genetic values
            g[[i]] <- as.numeric(M[[i]]%*%a[[i]])
          }
        }
        # add genetic values to phenotypes
        y <- y + g[[i]]
      }
    }
  }
  # simulate residuals
  e <- as.numeric(scale(rnorm(ifelse(class(M)[1] == "list",nrow(M[[1]]),nrow(M)),0,1)))
  # orthogonalize residuals to genetic values
  s2e <- Vy*(1-sum(h2))
  e <- sqrt(s2e/var(e))*(e - as.numeric(crossprod(y,e)/crossprod(y))*y)
  # repeat orthogonalization until residual variance is perfect
  while(abs(var(e) - s2e) > sim.error){
    e <- sqrt(s2e/var(e))*(e - as.numeric(crossprod(y,e)/crossprod(y))*y)
  }
  # add residuals to phenotypes
  y <- y + e
  # get heritability
  if(class(M)[1] == "matrix"){
    h2sim <- var(g)/var(y)
  } else{
    h2sim <- unlist(lapply(g,var))/var(y)
  }
  # discretize phenotypes if binary=TRUE
  if(binary){
    if(is.null(p)){
      p <- 0.5
    }
    findk <- function(k){
      t <- exp(y+k)/(1+exp(y+k))
      return(mean(t)-p)
    }
    k <- uniroot(findk,c(-50,50))$root
    y <- as.numeric(runif(length(y),0,1) < exp(y+k)/(1+exp(y+k)))
  }
  #-------------#
  # return data #
  #-------------#
  return(list(y=y,a=a,g=g,e=e,h2=h2sim))
}

#' Simulate a population info
#'
#' This function simulates a population info, with individuals' ID, sex, generation and parent information.
#'
#' @param n The size of each generation simulated. If a single number is provided, all generations will be of the same size.
#' @param n.gen The number of generations to simulate after generation zero. The generation counter starts at zero, being that the founder population. By default, \code{n.gen=0}, and if \code{n.gen} is set to any positive number, e.g. \code{n.gen=1}, that means that the info on the founder generation and on one generation after that will be generated.
#' @param p.male The proportion of male individuals in the population.
#' @param pop.info The information on a base population on which the code should build new generations.
#' @param pedigree Should the pedigree be generated? Default is \code{pedigree=FALSE}.
#' @param pre.id An ID label to paste at the begining of the default ID generated by the function.
#' @return A data frame with the population info, and the pedigree when required.
#' \describe{
#' \item{info}{The population info.}
#' \item{ped}{The pedigree matrix.}
#' }
#' @examples simPopInfo(10,0)
#' @examples simPopInfo(10,2)
#' @examples simPopInfo(3,1,pedigree=TRUE)
#' @examples simPopInfo(10,2,0.25)
#' @examples simPopInfo(c(5,10,15),2)
#' @export
simPopInfo <- function(n,n.gen=0,p.male=0.5,pop.info=NULL,pedigree=FALSE,pre.id=""){
  #-----------------------#
  # check values/elements #
  #-----------------------#
  if(any(n <= 0)){
    stop("Values in 'n' must be positive.")
  }
  if(any(n > 999999)){
    stop("I cannot simulate a generation with more than 999,999 individuals. Don't be too ambitious!")
  }
  if(n.gen < 0){
    stop("'n.gen' cannot be negative.")
  }
  if(n.gen > 999){
    stop("I cannot simulate more than 999 generations, sorry.")
  }
  if(p.male < 0 | p.male > 1){
    stop("'p.male' must be within the zero-one interval.")
  }
  if(n.gen > 0){
    if(length(n) == 1){
      n <- rep(n,n.gen+1)
    } else if(length(n) != n.gen + 1){
      stop("length(n) must be equal to 'n.gen + 1'")
    }
  }
  #-------------------#
  # internal function #
  #-------------------#
  char0 <- function(n){
    return(paste(rep(0,6-n),collapse=""))
  }
  #--------------------------#
  # simulate population info #
  #--------------------------#
  # generation zero
  if(is.null(pop.info)){
    pop.info <- data.frame(ID=1:n[1],parent1=rep(NA,n[1]),parent2=rep(NA,n[1]),sex=rep(0,n[1]),generation=rep(0,n[1]))
    pop.info$sex[sample(n[1],round(p.male*n[1]))] <- 1
    nID <- matrix(nchar(as.character(pop.info$ID[pop.info$generation == 0])),ncol=1)
    pop.info$ID <- paste0(pre.id,"G",paste(rep(0,2),collapse=""),"0",c("M","F")[2-pop.info$sex[pop.info$generation == 0]],apply(nID,1,char0),pop.info$ID[pop.info$generation == 0])
  }
  if(n.gen > 0){
    for(g in 1:n.gen){
      tmp <- data.frame(ID=1:n[g+1],parent1=rep(NA,n[g+1]),parent2=rep(NA,n[g+1]),sex=rep(0,n[g+1]),generation=rep(g,n[g+1]))
      tmp$sex[sample(n[g+1],round(p.male*n[g+1]))] <- 1
      nID <- matrix(nchar(as.character(tmp$ID[tmp$generation == g])),ncol=1)
      tmp$ID <- paste0(pre.id,"G",paste(rep(0,3-nchar(as.character(g))),collapse=""),g,c("M","F")[2-tmp$sex[tmp$generation == g]],apply(nID,1,char0),tmp$ID[tmp$generation == g])
      tmp$parent1 <- sample(pop.info$ID[pop.info$sex == 1 & pop.info$generation >= g-4],n[g+1],replace=TRUE,prob=(pop.info$generation[pop.info$sex == 1 & pop.info$generation >= g-4]+1)/sum(pop.info$generation[pop.info$sex == 1 & pop.info$generation >= g-4]+1))
      tmp$parent2 <- sample(pop.info$ID[pop.info$sex == 0 & pop.info$generation >= g-2],n[g+1],replace=TRUE,prob=(pop.info$generation[pop.info$sex == 0 & pop.info$generation >= g-2]+1)/sum(pop.info$generation[pop.info$sex == 0 & pop.info$generation >= g-2]+1))
      pop.info <- rbind(pop.info,tmp)
      rm(tmp)
    }
  }
  #-------------#
  # return data #
  #-------------#
  if(pedigree)
  {
    # generate pedigree matrix, if required
    A <- 2*kinship(pop.info$ID,pop.info$parent1,pop.info$parent2,2-pop.info$sex)
    return(list(info=pop.info,ped=A))
  } else return(pop.info)
}

#' Simulates genotypes for generations of individuals
#'
#' This function outputs SNP-genotypes for generations of individuals, based on the given SNPs for generation zero (founder population) and on the population info in the same format generated by \code{simPopInfo}.
#'
#' @param M A phased SNP matrix for generation zero.
#' @param pop.info The population info, in the same format generated by \code{simPopInfo}.
#' @return A matrix with the SNPs for all individuals.
#' @examples M <- simGeno(10,rep(0.5,3),phased=TRUE)
#' @examples pop.info <- simPopInfo(c(10,15),1)
#' @examples popBreed(M,pop.info)
#' @seealso \code{\link{simGeno}} \code{\link{simPopInfo}}
#' @export
popBreed <- function(M,pop.info){
  #--------------------#
  # internal functions #
  #--------------------#
  ParentGeno <- function(M,k){
    ListGeno <- function(M,k){
      genoTMP <- function(M){
        return(M[k,])
      }
      return(lapply(M,genoTMP))
    }
    return(matrix(unlist(ListGeno(M,k)),nrow=length(M),byrow=TRUE))
  }
  geno.offspring <- function(i){
    geno.parent <- list(ParentGeno(M,which(pop1$ID == pop2[i,]$parent1)),ParentGeno(M,which(pop1$ID == pop2[i,]$parent2)))
    cross.over <- rpois(2,0.001*ncol(M[[1]]))
    geno.i <- numeric(0)
    for(k in 1:2)
    {
      if(cross.over[k] == 0){
        geno.i <- rbind(geno.i,geno.parent[[k]][sample(nrow(geno.parent[[k]]),nrow(geno.parent[[k]])/2),])
      } else{
        if(cross.over[k] > ncol(M[[1]])-1){
          cross.over[k] <- ncol(M[[1]])-1
        }
        tmp <- sort(sample(ncol(M[[1]])-1,cross.over[k]))
        aux1 <- c(1,tmp+1)
        aux2 <- c(tmp,ncol(M[[1]]))
        rm(tmp)
        tmp <- numeric(0)
        for(j in 1:length(aux1)){
          if(length(M) > 2){
            tmp <- cbind(tmp,geno.parent[[k]][sample(nrow(geno.parent[[k]]),nrow(geno.parent[[k]])/2),aux1[j]:aux2[j]])
          } else tmp <- c(tmp,geno.parent[[k]][sample(nrow(geno.parent[[k]]),nrow(geno.parent[[k]])/2),aux1[j]:aux2[j]])
        }
        geno.i <- rbind(geno.i,tmp)
        rm(tmp)
      }
    }
    return(geno.i)
  }
  getRow <- function(A,n){
    return(A[n,])
  }
  #----------------------#
  # execute the breeding #
  #----------------------#
  for(g in 1:max(pop.info$generation)){
    pop1 <- pop.info[1:nrow(M[[1]]),]
    pop2 <- pop.info[pop.info$generation == g,]
    Mtmp <- lapply(lapply(lapply(strsplit(as.character(1:nrow(pop2)),""),paste,collapse=""),as.numeric),geno.offspring)
    for(k in 1:length(M)){
      M[[k]] <- rbind(M[[k]],matrix(unlist(lapply(Mtmp,getRow,n=k)),length(Mtmp),ncol(M[[1]]),byrow=TRUE))
    }
  }
  #-------------#
  # return data #
  #-------------#
  return(M)
}

#' Genomic relationship matrix
#'
#' This function builds the genomic relationship matrix based on SNP-genotypes.
#'
#' @param M A SNP-genotype matrix (unphased).
#' @param method The method to use for standardizing the SNP-genotypes. Deafult is \code{method="VanRaden"}, and alternatively can be set \code{method="Individual"}, which will standardize each SNP individually.
#' @return The genomic relationship matrix.
#' @examples M <- simGeno(5,rep(0.5,1000))
#' @examples M <- M[,apply(M,2,var) != 0]
#' @examples mkGRM(M)
#' @examples mkGRM(M,method="Individual")
#' @export
mkGRM <- function(M,method="VanRaden")
{
  if(method == "VanRaden"){
    M <- M-matrix(colMeans(M),nrow(M),ncol(M),byrow=TRUE)
    vM <- sum(colVars(M))
  } else if(method == "Individual"){
    M <- scale(M)
    vM <- ncol(M)
  } else stop("Please choose a valid method.")
  return(tcrossprod(M)/vM)
}

#' Pedigree relationship matrix
#'
#' This function builds the pedigree relationship matrix based on a population info in the same format as the generated by function \code{simPopInfo}. If you are generating the population info using \code{simPopInfo}, the pedigree relationship matrix can be created directly from \code{simPopInfo} by setting the parameter \code{pedigree=TRUE}.
#'
#' @param pop.info The population info, in the same format generated by \code{simPopInfo}.
#' @return The pedigree relationship matrix.
#' @examples mkPED(simPopInfo(10,0))
#' @export
mkPED <- function(pop.info)
{
  return(2*kinship(pop.info$ID,pop.info$parent1,pop.info$parent2,2-pop.info$sex))
}

#' Inverse of a matrix
#'
#' This function inverts any matrix, using the best method for that specific matrix.
#'
#' @param A A matrix.
#' @param use.eigen A logical indicating if the eigen-decomposition should be used to obtain the inverse. Valid only for symmetric matrices.
#' @param U A matrix of the eigenvectors (optional).
#' @param L A vector with the eigenvalues (optional).
#' @return The inverse of the matrix.
#' @export
MatrixInv <- function(A,use.eigen=FALSE,U=NULL,L=NULL)
{
  if(!use.eigen){
    Ainv <- try(solve(A),silent=TRUE)
    if(class(Ainv)[1] == "try-error")
    {
      Ainv <- try(solve(A+diag(1e-6,nrow(A))),silent=TRUE)
      if(class(Ainv)[1] == "try-error")
      {
        Ainv <- ginv(A)
      }
    }
  } else{
    if(!isSymmetric(A)){
      stop("'use.eigen=TRUE' is only valid for square and symmetric matrices.")
    }
    if(is.null(U) || is.null(L)){
      tmp <- eigen(A)
      U <- tmp$vectors
      L <- tmp$values
      rm(tmp)
    } else{
      if(nrow(U) != ncol(U)){
        stop("Check the dimensions of 'U'.")
      }
      if(length(L) != ncol(U)){
        stop("The number of eigenvalues must match the number of eigenvectors.")
      }
    }
  }
  return(Ainv)
}

#' REML estimates of variance components
#'
#' This function estimates the variance components and the fixed effects in a mixed model, using the restricted maximum likelihood (REML).
#'
#' @param y The response variable.
#' @param X A matrix of fixed effects.
#' @param Z The design matrix (or list of matrices) of the random effects. If all design matrices are an identity, this argument does not need to be provided.
#' @param VarMat The variance-covariance matrix (or list of matrices) of the random effects.
#' @param R A variance-covariance matrix for the residuals, if they are not independent and identically distributed (i.e. R != In)
#' @return
#' \describe{
#' \item{fixed}{The estimated fixed effects. The first element is always the intercept.}
#' \item{VarComp}{The estimated variance components for all random effects and for the residual.}
#' }
#' @examples ## SINGLE COMPONENT WITH FIXED EFFECTS ##
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=runif(100,0.05,0.5))
#' @examples # simulate phenotypes
#' @examples pheno <- simPheno(M,0.6,Vy=10)
#' @examples # add fixed effects
#' @examples X <- cbind(rep(1,nrow(M)),runif(nrow(M)),rbinom(nrow(M),1,0.5))
#' @examples b <- c(5,-1,2)
#' @examples y <- pheno$y + as.numeric(X%*%b)
#' @examples # REML
#' @examples G <- mkGRM(M)
#' @examples REMLsolve(y=y,X=X,VarMat=G)
#'
#' @examples ## TWO VARIANCE COMPONENTS ##
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=runif(100,0.05,0.5))
#' @examples # simulate phenotypes
#' @examples pheno <- simPheno(list(M[,1:70],M[,71:100]),c(0.4,0.2),Vy=10)
#' @examples # add fixed effects
#' @examples X <- cbind(rep(1,nrow(M)),runif(nrow(M)),rbinom(nrow(M),1,0.5))
#' @examples b <- c(5,-1,2)
#' @examples y <- pheno$y + as.numeric(X%*%b)
#' @examples # REML
#' @examples G <- list(mkGRM(M[,1:70]),mkGRM(M[,71:100]))
#' @examples REMLsolve(y=y,X=X,VarMat=G)
#' @examples REMLsolve(y=y,X=X,VarMat=list(diag(1,70),diag(1,30)),Z=list(M[,1:70],M[,71:100]))
#' @examples tmp <- REMLsolve(y=y,X=X,VarMat=list(diag(1,70),diag(1,30)),Z=list(M[,1:70],M[,71:100]))$VarComp
#' @examples c(sum(apply(M[,1:70],2,var)),sum(apply(M[,71:100],2,var)),1)*tmp
#' @export
REMLsolve <- function(y,X=NULL,Z=NULL,VarMat=NULL,R=NULL){
  #--------------------#
  # internal functions #
  #--------------------#
  unlistA <- function(A){
    if(class(A)[1] == "list"){
      Atmp <- numeric(0)
      for(j in 1:length(A)){
        Atmp <- cbind(Atmp,A[[j]])
      };rm(j)
      return(Atmp)
    } else return(A)
  }
  ImatList <- function(k,n){
    if(k > 1){
      A <- vector("list",k)
      for(j in 1:k){
        A[[j]] <- diag(1,n)
      };rm(j)
      return(A)
    } else return(diag(1,n))
  }
  weightA <- function(A,x){
    if(class(A)[1] == "list"){
      for(j in 1:length(A)){
        A[[j]] <- A[[j]]*x[j]
      };rm(j)
      return(A)
    } else return(A*x)
  }
  sumA <- function(A){
    if(class(A)[1] == "list"){
      for(j in 2:length(A)){
        A[[1]] <- A[[1]] + A[[j]]
      };rm(j)
      return(A[[1]])
    } else return(A)
  }
  #------------------------------------------------------------#
  # identify number of random components (excluding residuals) #
  #------------------------------------------------------------#
  M <- ifelse(class(VarMat)[1] == "list",length(VarMat),1)
  #-------------------------------------------#
  # define variances of Zu (if Z != identity) #
  #-------------------------------------------#
  s2z <- NULL
  if(!is.null(Z)){
    if(class(VarMat)[1] == "list"){
      s2z <- unlist(lapply(lapply(Z,apply,2,var),sum))
      tmpG <- function(x){
        return(scale(Z[[x]],scale=FALSE)%*%tcrossprod(VarMat[[x]],scale(Z[[x]],scale=FALSE))/s2z[x])
      }
      tmp <- list(0)
      for(j in 1:M){
        tmp[[j]] <- j
      }
      G <- lapply(tmp,tmpG)
    } else{
      s2z <- sum(apply(Z,2,var))
      G <- scale(Z,scale=FALSE)%*%tcrossprod(VarMat,scale(Z,scale=FALSE))/s2z
    }
  } else G <- VarMat
  if(is.null(R)){
    R <- diag(1,length(y))
  }
  if(is.null(X)){
    X <- matrix(1,length(y),1)
  }
  if(class(G)[1] == "list"){
    G[[length(G)+1]] <- R
  } else G <- list(G,R)
  names(G) <- c(paste0("Var",1:M),"Residuals")
  #------#
  # REML #
  #------#
  tmp <- MMEst(y,Cofactor=X,VarList=G)[[1]]
  s2 <- tmp$Sigma2
  if(!is.null(s2z)){
    s2[-length(s2)] <- s2[-length(s2)]/s2z
  }
  if(!is.null(Z)){
    Ztmp <- Z
  } else Ztmp <- ImatList(ifelse(class(VarMat)[1] == "list",length(VarMat),1),length(y))
  # variance(s) of the random effects weighted by variance component(s)
  # and aligned to phenotypes by design matrix(ces) Z
  Vtmp <- weightA(VarMat,s2[-length(s2)])
  if(class(Vtmp)[1] == "list"){
    for(j in 1:length(Vtmp)){
      Vtmp[[j]] <- tcrossprod(tcrossprod(Ztmp[[j]],Vtmp[[j]]),Ztmp[[j]])
    };rm(j)
  } else Vtmp <- tcrossprod(tcrossprod(Ztmp,Vtmp),Ztmp)
  Ztmp <- unlistA(Ztmp)
  # total phenotypic variance
  if(is.null(R)){
    Vy <- diag(1,length(y))*s2[length(s2)]
  } else Vy <- R*s2[length(s2)]
  Vy <- Vy + sumA(Vtmp)
  Vy_inv <- MatrixInv(Vy)
  rm(Vy);gc()
  if(!is.null(X)){
    b <- as.numeric((MatrixInv(crossprod(X,Vy_inv)%*%X)%*%crossprod(X,Vy_inv)%*%y))
  } else b <- as.numeric(apply(Vy_inv,2,sum)%*%y)/sum(Vy_inv)
  #------------------#
  # return solutions #
  #------------------#
  return(list(FixedEff=b,VarComp=s2))
}

#' BLUP
#'
#' This function performs BLUP in a mixed model. If not provided, the function will estimate the variance components using REML.
#'
#' @param y The response variable.
#' @param X A matrix of fixed effects.
#' @param FixedEff The effects of the fixed effects, if previously estimated or provided.
#' @param Z The design matrix (or list of matrices) of the random effects. If all design matrices are an identity, this argument does not need to be provided.
#' @param VarMat The variance-covariance matrix (or list of matrices) of the random effects.
#' @param R A variance-covariance matrix for the residuals, if they are not independent and identically distributed (i.e. if R != In).
#' @param VarComp The variance components of each random effect, including the residuals. If \code{VarComp=NULL}, the function will estimate these variances using REML.
#' @param TestGroup The index of observations from a test population. This will make the function perform prediction over these observations.
#' @param verbose Logical: should steps of the co,putation be printed on screen?
#' @return
#' \describe{
#' \item{fixed}{The estimated fixed effects. The first element is always the intercept.}
#' \item{VarComp}{The estimated variance components for all random effects, and for the residual.}
#' \item{random}{The predicted random effects.}
#' \item{Sigma2_BLUP}{The variance-covariance matrix of the parameters estimated/predict with BLUP.}
#' \item{Rpred}{The accuracy of the prediction performed for \code{TestGroup}.}
#' }
#' @examples ## SINGLE COMPONENT WITH FIXED EFFECTS ##
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=runif(100,0.05,0.5))
#' @examples # simulate phenotypes
#' @examples pheno <- simPheno(M,0.6,Vy=10)
#' @examples # add fixed effects
#' @examples X <- cbind(rep(1,nrow(M)),runif(nrow(M)),rbinom(nrow(M),1,0.5))
#' @examples b <- c(5,-1,2)
#' @examples y <- pheno$y + as.numeric(X%*%b)
#' @examples # REML
#' @examples G <- mkGRM(M)
#' @examples tmp <- REMLsolve(y=y,X=X,VarMat=G)
#' @examples tmp
#' @examples # BLUP
#' @examples aux <- BLUPsolve(y=y,X=X,VarMat=G,VarComp=tmp$VarComp)
#' @examples str(aux)
#' @examples # BLUP estimating variance components
#' @examples aux <- BLUPsolve(y=y,X=X,VarMat=G)
#' @examples str(aux)
#'
#' @examples ## TWO VARIANCE COMPONENTS ##
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=runif(100,0.05,0.5))
#' @examples # simulate phenotypes
#' @examples pheno <- simPheno(list(M[,1:70],M[,71:100]),c(0.4,0.2),Vy=10)
#' @examples # BLUP using GRMs
#' @examples G <- list(mkGRM(M[,1:70]),mkGRM(M[,71:100]))
#' @examples aux <- BLUPsolve(y=pheno$y,VarMat=G)
#' @examples str(aux)
#' @examples # prediction
#' @examples aux <- BLUPsolve(y=pheno$y,VarMat=G,TestGroup=401:500)
#' @examples str(aux)
#' @examples # BLUP using genotype matrices
#' @examples aux <- BLUPsolve(y=pheno$y,VarMat=list(diag(1,70),diag(1,30)),Z=list(M[,1:70],M[,71:100]))
#' @examples str(aux)
#' @examples # prediction
#' @examples aux <- BLUPsolve(y=pheno$y,VarMat=list(diag(1,70),diag(1,30)),Z=list(M[,1:70],M[,71:100]),TestGroup=401:500)
#' @examples str(aux)
#' @seealso \code{\link{REMLsolve}}
#' @export
BLUPsolve <- function(y,X=NULL,FixedEff=NULL,Z=NULL,VarMat=NULL,R=NULL,VarComp=NULL,TestGroup=NULL,verbose=TRUE){
  #---------------#
  # safety checks #
  #---------------#
  if(!is.null(FixedEff)){
    if(is.null(X)){
      if(length(FixedEff) > 1){
        stop("If more than the overall mean is provided in 'FixedEff', you must also provide 'X'.")
      }
    } else if(ncol(X) != length(FixedEff)){
      stop("Number of columns in 'X' does not match the length of 'FixedEff'.")
    }
  }
  #--------------------#
  # internal functions #
  #--------------------#
  unlistA <- function(A){
    if(class(A)[1] == "list"){
      Atmp <- numeric(0)
      for(j in 1:length(A)){
        Atmp <- cbind(Atmp,A[[j]])
      };rm(j)
      return(Atmp)
    } else return(A)
  }
  ImatList <- function(k,n){
    if(k > 1){
      A <- vector("list",k)
      for(j in 1:k){
        A[[j]] <- diag(1,n)
      };rm(j)
      return(A)
    } else return(diag(1,n))
  }
  weightA <- function(A,x){
    if(class(A)[1] == "list"){
      for(j in 1:length(A)){
        A[[j]] <- A[[j]]*x[j]
      };rm(j)
      return(A)
    } else return(A*x)
  }
  sumA <- function(A){
    if(class(A)[1] == "list"){
      for(j in 2:length(A)){
        A[[1]] <- A[[1]] + A[[j]]
      };rm(j)
      return(A[[1]])
    } else return(A)
  }
  subA <- function(A,x=NULL,row=TRUE,col=FALSE){
    if(!is.null(A) & !is.null(x)){
      if(class(A)[1] == "list"){
        if(row){
          # A <- lapply(A,subset,subset=x)
          for(j in 1:length(A)){
            A[[j]] <- subset(A[[j]],subset=x)
          }
        }
        if(col){
          # A <- lapply(A,subset,select=x)
          for(j in 1:length(A)){
            A[[j]] <- subset(A[[j]],select=x)
          }
        }
      } else{
        if(row){
          A <- subset(A,subset=x)
        }
        if(col){
          A <- subset(A,select=x)
        }
      }
    }
    return(A)
  }
  subY <- function(y,x=NULL){
    if(!is.null(x)){
      y <- y[x]
    }
    return(y)
  }
  #-----------------------------------------------------------#
  # create logical subset TRN vector if test group is defined #
  #-----------------------------------------------------------#
  if(!is.null(TestGroup)){
    TRN <- !(1:length(y) %in% TestGroup)
  } else TRN <- NULL
  #---------------------------------------------#
  # run REML if variance components are missing #
  #---------------------------------------------#
  if(is.null(VarComp)){
    if(verbose){
      cat("Estimating variance components with REML... ",sep="")
    }
    if(is.null(Z)){
      tmp <- REMLsolve(subY(y,TRN),X=subA(X,TRN),Z=subA(Z,TRN),VarMat=subA(A=VarMat,x=TRN,col=TRUE),R=subA(R,TRN,col=TRUE))
    } else tmp <- REMLsolve(subY(y,TRN),X=subA(X,TRN),Z=subA(Z,TRN),VarMat=VarMat,R=subA(R,TRN,col=TRUE))
    if(is.null(FixedEff)){
      FixedEff <- tmp$FixedEff
    } else{
      if(any(abs(FixedEff - tmp$FixedEff) > 1e-4)){
        FixedEff <- tmp$FixedEff
        warning("'FixedEff' provided differed too much (> 1e-4) from those weighted by the variance components (obtained with REML), and were therefore updated.")
      }
    }
    VarComp <- tmp$VarComp
    if(verbose){
      cat("DONE\n",sep="")
    }
  }
  #--------------------------------#
  # prepare elements to solve BLUP #
  #--------------------------------#
  if(verbose){
    cat("Solving Henderson's MME to obtain BLUPs... ",sep="")
  }
  # variance of random effects into a single matrix
  S <- as.matrix(bdiag(weightA(VarMat,VarComp[-length(VarComp)])))
  # design matrix(ces) for the random effects subset if test group is defined
  if(!is.null(Z)){
    Ztmp <- subA(Z,x=TRN)
  } else Ztmp <- subA(ImatList(ifelse(class(VarMat)[1] == "list",length(VarMat),1),length(y)),TRN)
  # variance(s) of the random effects weighted by variance component(s)
  # and aligned to phenotypes by design matrix(ces) Z
  Vtmp <- weightA(VarMat,VarComp[-length(VarComp)])
  if(class(Vtmp)[1] == "list"){
    for(j in 1:length(Vtmp)){
      Vtmp[[j]] <- tcrossprod(tcrossprod(Ztmp[[j]],Vtmp[[j]]),Ztmp[[j]])
    };rm(j)
  } else Vtmp <- tcrossprod(tcrossprod(Ztmp,Vtmp),Ztmp)
  Ztmp <- unlistA(Ztmp)
  # total phenotypic variance
  if(is.null(R)){
    Vy <- subA(diag(1,length(y)),TRN,col=TRUE)*VarComp[length(VarComp)]
  } else Vy <- subA(R,TRN,col=TRUE)*VarComp[length(VarComp)]
  Vy <- Vy + sumA(Vtmp)
  Vy_inv <- MatrixInv(Vy)
  rm(Vy);gc()
  #-------------------------------------------#
  # solve BLUP for fixed effects if necessary #
  #-------------------------------------------#
  if(is.null(FixedEff)){
    if(!is.null(X)){
      FixedEff <- as.numeric((MatrixInv(crossprod(subA(X,TRN),Vy_inv)%*%subA(X,TRN))%*%crossprod(subA(X,TRN),Vy_inv)%*%subY(y,TRN)))
    } else FixedEff <- as.numeric(apply(Vy_inv,2,sum)%*%subY(y,TRN))/sum(Vy_inv)
  }
  #-------------------------------------------------------#
  # correct phenotypes for the fixed effects              #
  # start ypred on the fixed effects (will be used later) #
  #-------------------------------------------------------#
  if(!is.null(X)){
    ypred <- as.numeric(X%*%FixedEff)
  } else ypred <- FixedEff
  yorig <- y
  y <- y - ypred
  #------------------------------#
  # solve BLUP on random effects #
  #------------------------------#
  if(!is.null(X)){
    XpVy_inv <- crossprod(subA(X,TRN),Vy_inv)
    Pr2 <- Vy_inv - crossprod(XpVy_inv,MatrixInv(XpVy_inv%*%subA(X,TRN)))%*%XpVy_inv
  } else Pr2 <- Vy_inv - matrix(apply(Vy_inv,1,sum),ncol=1)%*%matrix(apply(Vy_inv,2,sum),nrow=1)/sum(Vy_inv)
  if(ifelse(class(VarMat)[1] == "list",length(VarMat),1) == 1 & is.null(Z) & is.null(TRN)){
    # random effects
    aux <- as.numeric(S%*%Vy_inv%*%subY(y,TRN))
    # reliability of random effects
    r2aux <- diag(S%*%Pr2)
  } else{
    # random effects
    aux <- as.numeric(S%*%t(Ztmp)%*%Vy_inv%*%subY(y,TRN))
    # reliability of random effects
    r2aux <- diag(S%*%t(Ztmp)%*%Pr2%*%Ztmp)
  }
  if(class(VarMat)[1] == "list"){
    tmp <- aux
    r2tmp <- r2aux
    aux <- vector("list",length(VarMat))
    r2aux <- vector("list",length(VarMat))
    for(j in 1:length(VarMat)){
      aux[[j]] <- tmp[1:ncol(VarMat[[j]])]
      tmp <- tmp[-(1:ncol(VarMat[[j]]))]
      r2aux[[j]] <- r2tmp[1:ncol(VarMat[[j]])]
      r2tmp <- r2tmp[-(1:ncol(VarMat[[j]]))]
    };rm(j)
  }
  if(verbose){
    cat("DONE\n",sep="")
  }
  #--------------------------------------------#
  # print extra information if verbose == TRUE #
  #--------------------------------------------#
  if(verbose){
    if(!is.null(Z)){
      if(class(Z)[1] == "list"){
        ZI <- rep(FALSE,length(Z))
        for(j in 1:length(Z)){
          if(isSymmetric(Z[[j]])){
            ZI[j] <- ifelse(sum(diag(Z[[j]] == 1)) == nrow(Z[[j]]),TRUE,FALSE)
          }
        }
        if(sum(ZI) == length(ZI)){
          cat("\nNOTICE: The model assumed that Z's are all identity matrices.",sep="")
        } else if(sum(ZI) > 0){
          cat("\nNOTICE: The model assumed that Z is an identity matrix for random component(s): ",sep="")
          cat(which(ZI),".",sep="")
        }
      } else{
        ZI <- FALSE
        if(isSymmetric(Z)){
          ZI <- ifelse(sum(diag(Z == 1)) == nrow(Z),TRUE,FALSE)
        }
        if(ZI){
          cat("\nNOTICE: The model assumed that Z is an identity matrix.",sep="")
        }
      }
      if(sum(ZI) > 0){
        cat("\n        Therefore Zu is equivalent to RandomEff for such component(s).\n",sep="")
      }
    } else {
      if(class(VarComp)[1] == "list"){
        cat("\nNOTICE: The model assumed that Z's are all identity matrices.",sep="")
        cat("\n        Therefore Zu is equivalent to RandomEff for all components.\n",sep="")
      } else{
        cat("\nNOTICE: The model assumed that Z is an identity matrix.",sep="")
        cat("\n        Therefore Zu is equivalent to RandomEff.\n",sep="")
      }
    }
  }
  #------------------------------------#
  # predicted phenotypes and residuals #
  #------------------------------------#
  g <- numeric(0)
  if(class(VarMat)[1] == "list"){
    for(j in 1:length(VarMat)){
      if(!is.null(Z)){
        g <- cbind(g,as.numeric(Z[[j]]%*%aux[[j]]))
      } else g <- cbind(g,aux[[j]])
    };rm(j)
    ypred <- ypred + apply(g,1,sum)
  } else{
    if(!is.null(Z)){
      g <- as.numeric(Z%*%aux)
    } else g <- aux
    ypred <- ypred + g
  }
  e <- yorig - ypred
  #-----------#
  # model fit #
  #-----------#
  Zu <- as.data.frame(g)
  names(Zu) <- names(VarComp[-length(VarComp)])
  if(class(VarMat)[1] == "list"){
    if(is.null(TestGroup)){
      FitAcc <- c(as.numeric(cor(y,g)),cor(y,apply(g,1,sum)))
    } else     FitAcc <- c(as.numeric(cor(y[TRN],g[TRN,])),cor(y[TRN],apply(g[TRN,],1,sum)))
    names(FitAcc) <- c(names(VarComp[-length(VarComp)]),"CorOverSum")
  } else FitAcc <- cor(subY(y,TRN),subY(g,TRN))
  FitAccY <- cor(subY(yorig,TRN),subY(ypred,TRN))
  if(!is.null(TRN)){
    Zufit <- Zu[TRN,]
  } else Zufit <- Zu
  modelFit <- list(yCorrectedFixed=subY(y,TRN),yfit=subY(ypred,TRN),Zu=Zufit,residuals=subY(e,TRN),FitAccRand=FitAcc,FitAccY=FitAccY)
  #------------------#
  # model's accuracy #
  #------------------#
  if(!is.null(TestGroup)){
    if(class(VarMat)[1] == "list"){
      PredAcc <- c(as.numeric(cor(y[!TRN],g[!TRN,])),cor(y[!TRN],apply(g[!TRN,],1,sum)))
      names(PredAcc) <- c(names(VarComp[-length(VarComp)]),"CorOverSum")
    } else PredAcc <- cor(y[!TRN],g[!TRN])
    PredAccY <- cor(yorig[!TRN],ypred[!TRN])
    Zupred <- Zu[!TRN,]
    modelPred <- list(yCorrectedFixed=subY(y,!TRN),ypred=subY(ypred,!TRN),Zu=Zupred,residuals=subY(e,!TRN),PredAccRand=PredAcc,PredAccY=PredAccY)
    return(list(FixedEff=FixedEff,VarComp=VarComp,RandomEff=list(eff=aux,r2=r2aux),modelFit=modelFit,modelPred=modelPred))
  } else return(list(FixedEff=FixedEff,VarComp=VarComp,RandomEff=list(eff=aux,r2=r2aux),modelFit=modelFit))
}

#' Kappa statistics
#'
#' This function calculates the kappa statistics, to compare two variance-covariance matrices (G1 and G2) defined for the same group of individuals.
#' Kappa is computed using the eigen-decomposition of the matrices, \code{G=ULU'} in which U and L are the eigen-vectors and eigen-values, respectively. \code{Kmat = U2'(G1)U2}, and \code{kappa=diag(Kmat)}.
#' Finally, to compare G1 and G2, it is enough to plot \code{lambda} against \code{kappa}. The more linear their relationship, the more equivalent matrices G1 and G2 are, in terms of variance structure.
#'
#' @param G1 The Reference matrix.
#' @param G2 The matrix which we want to compare with G1.
#' @param U2 The eigen-vectors of G1. If not provided, the function will obtain them internally.
#' @param L2 The eigen-values of G2. If not provided, the function will obtain them internally.
#' @param Kmat Should the entire Kappa matrix be returned? Default is \code{Kmat=FALSE}.
#' @return
#' \describe{
#' \item{lambda}{The eigen-values of G2.}
#' \item{kappa}{The kappa statistics that compare G2 to G1.}
#' \item{Kmat}{The whole Kappa matrix.}
#' }
#' @examples # simulate uncorrelated genotypes
#' @examples M <- simGeno(n=100,AlleleFreq=rep(0.5,200))
#' @examples # three GRMs from the simulated genotypes
#' @examples # the first with 100 SNPs
#' @examples G1 <- mkGRM(M[,1:100])
#' @examples # the second with 100 SNPs that don't overlap those from the first
#' @examples G2 <- mkGRM(M[,101:200])
#' @examples # the third with all SNPs (100 SNPs overlap with those from the first matrix)
#' @examples G3 <- mkGRM(M)
#' @examples #-------------------------------------------------------------------#
#' @examples # compare GRMs using kappa
#' @examples # verify how that is different from directly comparing the values
#' @examples # verify how that is also different from comparing the eigen-values
#' @examples #-------------------------------------------------------------------#
#' @examples ###########
#' @examples # G1 x G2 #
#' @examples ###########
#' @examples par(mar=c(4.5,4.5,1,1))
#' @examples # plot G1 x G2
#' @examples plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G2),G2[upper.tri(G2)]),xlab="G1",ylab="G2")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G1) x eigen-values(G2)
#' @examples plot(eigen(G1)$values,eigen(G2)$values,xlab="eigen-values(G1)",ylab="eigen-values(G2)")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G2)xkappa
#' @examples tmp <- Gkappa(G1,G2)
#' @examples plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
#' @examples abline(0,1,lty=2,col=2)
#' @examples ###########
#' @examples # G1 x G3 #
#' @examples ###########
#' @examples par(mar=c(4.5,4.5,1,1))
#' @examples # plot G1xG3
#' @examples plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G3),G3[upper.tri(G3)]),xlab="G1",ylab="G3")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G1)xeigen-values(G3)
#' @examples plot(eigen(G1)$values,eigen(G3)$values,xlab="eigen-values(G1)",ylab="eigen-values(G3)")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G3)xkappa
#' @examples tmp <- Gkappa(G1,G3)
#' @examples plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
#' @examples abline(0,1,lty=2,col=2)
#'
#' @examples # simulate correlated genotypes
#' @examples M <- simGeno(n=100,AlleleFreq=rep(0.5,1000),LD=TRUE,phased=FALSE)
#' @examples # three GRMs from the simulated genotypes
#' @examples # the first with 500 SNPs
#' @examples G1 <- mkGRM(M[,1:500])
#' @examples # the second with 500 SNPs that don't overlap those from the first
#' @examples G2 <- mkGRM(M[,501:1000])
#' @examples # the third with all SNPs (500 SNPs overlap with those from the first matrix)
#' @examples G3 <- mkGRM(M)
#' @examples #-------------------------------------------------------------------#
#' @examples # compare GRMs using kappa
#' @examples # verify how that is different from directly comparing the values
#' @examples # verify how that is also different from comparing the eigen-values
#' @examples #-------------------------------------------------------------------#
#' @examples ###########
#' @examples # G1 x G2 #
#' @examples ###########
#' @examples par(mfrow=c(2,2),mar=c(4.5,4.5,1,1))
#' @examples # plot G1xG2
#' @examples plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G2),G2[upper.tri(G2)]),xlab="G1",ylab="G2")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G1)xeigen-values(G2)
#' @examples plot(eigen(G1)$values,eigen(G2)$values,xlab="lambda(G1)",ylab="lambda(G2)")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G2)xkappa
#' @examples tmp <- Gkappa(G1,G2)
#' @examples plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
#' @examples abline(0,1,lty=2,col=2)
#' @examples plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
#' @examples ###########
#' @examples # G1 x G3 #
#' @examples ###########
#' @examples par(mfrow=c(2,2),mar=c(4.5,4.5,1,1))
#' @examples # plot G1xG3
#' @examples plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G3),G3[upper.tri(G3)]),xlab="G1",ylab="G3")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G1)xeigen-values(G3)
#' @examples plot(eigen(G1)$values,eigen(G3)$values,xlab="lambda(G1)",ylab="lambda(G3)")
#' @examples abline(0,1,lty=2,col=2)
#' @examples # plot eigen-values(G3)xkappa
#' @examples tmp <- Gkappa(G1,G3)
#' @examples plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
#' @examples abline(0,1,lty=2,col=2)
#' @examples plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
#' @export
Gkappa <- function(G1,G2=NULL,U2=NULL,L2=NULL,Kmat=FALSE)
{
  if(is.null(U2) | is.null(L2)){
    tmp <- eigen(G2)
    U2 <- tmp$vectors
    L2 <- tmp$values
  }
  K <- crossprod(U2,G1)%*%U2
  if(Kmat){
    return(list(lambda=L2,kappa=diag(K),Kmat=K))
  } else return(list(lambda=L2,kappa=diag(K)))
}

#' Calculates the relationship between two populations
#'
#' This function calculates the relationship between two populations based on the singular-value decomposition of their genotype matrices. The relationship is obtained as in Cuyabano and Gondro, 2020 (in prep.), and the two populations must have the SNP-genotypes on the same set of SNPs.
#'
#' @param M1 The SNP-genotypes of population one.
#' @param M2 The SNP-genotypes of population two.
#' @return The parameters of the linear regression (intercept and slope) on the components built from the SVDs.
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=rep(0.5,100),LD=TRUE)
#' @examples # separate two populations of individuals
#' @examples # note that the populations do not need to have the same size
#' @examples M1 <- M[1:300,]
#' @examples M2 <- M[301:500,]
#' @examples popCor(M1,M2)
#' @export
popCor <- function(M1,M2)
{
  scaleM <- function(M)
  {
    centerM <- function(M)
    {
      return(scale(M,center=TRUE,scale=FALSE))
    }
    return(centerM(M)/sqrt(sum(apply(M,2,var))))
  }
  M1svd <- svd(scaleM(M1))
  M2svd <- svd(scaleM(M2))
  Rsvd <- svd(sqrt(nrow(M2)/nrow(M1))*crossprod(M2svd$v,M1svd$v)%*%diag(M1svd$d))
  tmp <- lm(Rsvd$d ~ M2svd$d)$coef
  tmp <- data.frame(a=-abs(tmp[1]),b=tmp[2],r=tmp[2]-abs(tmp[1]))
  return(tmp)
}

#' Calculates the expected limit of prediction accuracy
#'
#' This function calculates the expected limit of prediction accuracy, using the relationship between two populations based on the singular-value decomposition of their genotype matrices. This limit is obtained as in Cuyabano and Gondro, 2020 (in prep), and the two populations must have the SNP-genotypes on the same set of SNPs.
#'
#' @param a The intercept obtained with \code{popCor}.
#' @param b The slope obtained with \code{popCor}.
#' @param n The size of the test population (the population for which predictions are to be made).
#' @param h2 The heritability of the trait.
#' @param p The confidence level of the interval for the expected limit of prediction accuracy.
#' @return A data frame with the mean expected limit of prediction accuracy, as well as the boundaries of its confidence interval.
#' @examples # simulate genotypes
#' @examples M <- simGeno(n=500,AlleleFreq=rep(0.5,100))
#' @examples # separate two populations of individuals
#' @examples # note that the populations do not need to have the same size
#' @examples M1 <- M[1:300,]
#' @examples M2 <- M[301:500,]
#' @examples tmp <- popCor(M1,M2)
#' @examples LimPredAcc(a=tmp$a,b=tmp$b,n=200,h2=0.5)
#' @examples LimPredAcc(a=tmp$a,b=tmp$b,n=200,h2=seq(0.1,0.9,0.1))
#' @examples # compare the expected limits of prediction accuracy to the classic sqrt(h2)
#' @examples par(mar=c(4.5,4.5,1,1))
#' @examples curve(sqrt(x),xlim=c(0,1),ylim=c(0,1),xlab=expression(h^2),ylab="Limit of prediction accuracy")
#' @examples aux <- LimPredAcc(a=tmp$a,b=tmp$b,n=200,h2=seq(0.01,0.99,0.01))
#' @examples lines(seq(0.01,0.99,0.01),aux$mean,lty=2)
#' @examples lines(seq(0.01,0.99,0.01),aux$CI_lower,lty=3,col=2)
#' @examples lines(seq(0.01,0.99,0.01),aux$CI_upper,lty=3,col=2)
#' @examples legend("topleft",lty=c(1,2,3),col=c(1,1,2),legend=c("sqrt(h2)","R_expected","R_CI"))
#' @export
LimPredAcc <- function(a,b,n,h2,p=0.95){
  g <- a + b
  alpha <- 1-p
  Rmu <- ((1+sqrt(h2))^g-(1-sqrt(h2))^g)/((1+sqrt(h2))^g+(1-sqrt(h2))^g)
  Rlb <- (exp(-2*qnorm(1-alpha/2)/sqrt(n-3))*(1+sqrt(h2))^g-(1-sqrt(h2))^g)/(exp(-2*qnorm(1-alpha/2)/sqrt(n-3))*(1+sqrt(h2))^g+(1-sqrt(h2))^g)
  Rup <- (exp(2*qnorm(1-alpha/2)/sqrt(n-3))*(1+sqrt(h2))^g-(1-sqrt(h2))^g)/(exp(2*qnorm(1-alpha/2)/sqrt(n-3))*(1+sqrt(h2))^g+(1-sqrt(h2))^g)
  return(data.frame(mean=Rmu,CI_lower=Rlb,CI_upper=Rup))
}

#' Confidence interval for normal distributed data
#'
#' This function calculates the confidence interval for the mean of a set of normal distributed data.
#'
#' @param x A vector containing values to calculate confidence interval.
#' @param levels A vector containing the levels (factor) from the data, to calculate separate confidence intervals.
#' @param p The confidence level of the interval.
#' @return A data frame with the confidence interals.
#' @export
ConfIntMeanNorm <- function(x,levels=NULL,p=0.95)
{
  if(is.null(levels))
  {
    tmp <- t.test(x,p.level=p)$conf.int
    return(data.frame(min=tmp[1],max=tmp[2]))
  } else
  {
    tmp <- numeric(0)
    group <- levels(factor(levels))
    for(k in group)
    {
      tmp <- rbind(tmp,t.test(x[levels == k],p.level=p)$conf.int)
    }
    return(data.frame(levels=group,min=tmp[,1],max=tmp[,2]))
  }
}

#' Bootstrap confidence interval
#'
#' This function calculates the bootstrap confidence interval for the mean of a set of data.
#'
#' @param x A vector containing values to calculate confidence interval.
#' @param levels A vector containing the levels (factor) from the data, to calculate separate confidence intervals.
#' @param p The confidence level of the interval.
#' @param n.rep The number of bootstrap replicates.
#' @return A data frame with the confidence interals.
#' @export
ConfIntMeanBoot <- function(x,levels=NULL,conf=0.95,n.rep=10^4)
{
  if(length(x) == 1)
  {
    return(rep(x,2))
    stop("Warning: only one observation found, no confidence interval could be calculated")
  }
  levels <- levels[!is.na(x)]
  x <- x[!is.na(x)]
  if(is.null(levels))
  {
    u.boot <- one.boot(x, mean, R=n.rep)
    if(conf == 0.95)
    {
      conflim <- boot.ci(u.boot, type="bca")
    } else conflim <- boot.ci(u.boot,conf=conf, type="bca")
    return(conflim$bca[4:5])
    # return(data.frame(min=conflim$bca[4],max=conflim$bca[5]))
  } else
  {
    tmp <- numeric(0)
    group <- levels(factor(levels))
    for(k in group)
    {
      u.boot <- one.boot(x[levels == k], mean, R=n.rep)
      if(conf == 0.95)
      {
        conflim <- boot.ci(u.boot, type="bca")
      } else conflim <- boot.ci(u.boot,conf=conf, type="bca")
      tmp <- rbind(tmp,c(conflim$bca[4],conflim$bca[5]))
    }
    return(data.frame(levels=group,min=tmp[,1],max=tmp[,2]))
  }
}

#' The Marchenko-Pastur density
#'
#' Density function of the Marchenko-Pastur distribution with rate parameter equal to \code{a} and variance equal to \code{s2}.
#'
#' @param x A vectors of quantiles.
#' @param a Rate parameter (n/m = #observations/#variables).
#' @param s2 Variance parameter.
#' @return The density.
#' @export
dMP <- function(x,a,s2)
{
  a.minus <- s2*(1-sqrt(a))^2
  a.plus <- s2*(1+sqrt(a))^2
  dens <- rep(0,length(x))
  i <- which(x >= a.minus & x <= a.plus)
  dens[i] <- sqrt((a.plus-x[i])*(x[i]-a.minus))/(2*pi*s2*a*x[i])
  return(dens)
}

#' The Marchenko-Pastur cumulative distribution
#'
#' Cumulative distribution function of the Marchenko-Pastur distribution with rate parameter equal to \code{a} and variance equal to \code{s2}.
#'
#' @param x A vectors of quantiles.
#' @param a Rate parameter (n/m = #observations/#variables).
#' @param s2 Variance parameter.
#' @return The cumulative probability distribution.
#' @export
pMP <- function(x,a,s2)
{
  a.minus <- s2*(1-sqrt(a))^2
  a.plus <- s2*(1+sqrt(a))^2
  cdf <- rep(0,length(x))
  i <- which(x > a.plus)
  cdf[i] <- 1
  dMP.tmp <- function(t)
  {
    return(dMP(t,a,s2))
  }
  i <- which(x >= 0 & x <= a.plus)
  if(a > 1)
  {
    cdf[i] <- 1-1/a
  }
  for(i in which(x >= a.minus & x <= a.plus))
  {
    cdf[i] <- cdf[i] + integrate(dMP.tmp,lower=a.minus,upper=x[i])$value
  }
  return(cdf)
}

#' Plot eigen-values comparatively with The Marchenko-Pastur distribution
#'
#' @param L A vectors of eigen-values.
#' @param a Rate parameter (n/m = #observations/#variables).
#' @param s2 Variance parameter.
#' @param f The type of function to be plotted, \code{"density"} or \code{"cumulative"} (default).
#' @export
plotMP <- function(L,a,s2,f="cumulative"){
  par(mar=c(4.5,4.5,1,1))
  if(f == "density"){
    hist(L,breaks=length(L)/5,col=8,freq=FALSE,axes=FALSE,main="",xlab=expression(lambda * " (eigen-values)"),ylab="density")
    abline(h=0)
    axis(1,pos=0)
    axis(2)
    curve(dMP(x,nrow(M)/ncol(M),1),lwd=2,col=2,add=TRUE)
    legend("topright",bty="n",pch=c(22,NA),lty=c(NA,1),lwd=c(1,2),col=c(1,2),pt.bg=c(8,NA),legend=c("real data density","MP density function"))
  } else if("cumulative"){
    plot(quantile(L,seq(0,1,length.out=500)),seq(0,1,length.out=500),ylim=c(0,1),axes=FALSE,main="",xlab=expression(lambda * " (eigen-values)"),ylab="cumulative probability distribution")
    abline(h=0)
    axis(1,pos=0)
    axis(2)
    curve(pMP(x,nrow(M)/ncol(M),1),lwd=2,col=2,add=TRUE)
    legend("topleft",bty="n",pch=c(1,NA),lty=c(NA,1),lwd=c(1,2),col=c(1,2),legend=c("real data","MP cdf"))
  }
}
