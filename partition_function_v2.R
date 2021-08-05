Cut_plane_standard <- function(V,V.ID,V.label,tolerance,dist.max=NA,filter.ind=FALSE,method){
  d <- dim(V)[2]
  ess <- NULL
  xy.crv <- NULL
  xy.crv.af <- NULL
  sample.cut.count <- 1
  V.left.ID <-NULL
  V.right.ID <- NULL

  
  #if(sum(!is.na(V.label))>0){
  # cut if there is at least one NA label in the obs & at least two valid labels & at least from two groups.
  if(!filter.ind ){  
    if(is.na(dist.max)){
      V.unique <- unique(V[!is.na(V.label),])
      d1 <- dim(V.unique)[1]
      if(d1==1){
        dist.max <- 0
      }else{
        dist.max <- max(dist(V.unique))
      }
    }
    ##skip the cut if all obs are from the same group(pausing condition)
    ## or if labled vertices are same  
    skip.index <- (dist.max==0)
  }else{
    skip.index =1 # skip if all vertices are unlabled
  }
  
  
  if(!skip.index){
    #set.seed(Sys.time())

    x.min.star  <-  -sqrt(2)/2
    x.max.star <-    sqrt(2)/2
    y.min.star <- -sqrt(2)/2 
    y.max.star <-  sqrt(2)/2
  
    if(n.p>0){
      xy.crv <- cbind(rep(NA,n.p),rep(NA,n.p))
      xy.crv.af <- cbind(rep(NA,n.p),rep(NA,n.p))
      curves.bf <- cbind(rep(NA,n.p),rep(NA,n.p))
      curves.af <- cbind(rep(NA,n.p),rep(NA,n.p))
    } 
    

    #set.seed(seed)

    run.index <- 1
    
    while(run.index){
      
      deg <- sample(1:3,size=1)
      Pc     <- matrix(NA,nrow=(deg+1),ncol=2)
      x2 <- NULL
      y2 <- NULL
      if(deg>1){
        x2 <- sort(runif((deg-1),min=x.min.star,max=x.max.star))
      } 
      
      y2 <- runif((deg+1),min=y.min.star,max=y.max.star)
      Pc[,1]<-c(x.min.star,x2,x.max.star)
      Pc[,2]<- c(y2) 
      
      Pc0 <- Pc
      
      r   <- runif(1,-pi,pi)
      p.new <- mapf(V,r)
      x.temp <- as.matrix(p.new[,1])
      minmax <- f_bezier_minmax(Pc0,deg)
      
      minP0 <- min(p.new[,2])
      maxP0 <- max(p.new[,2])
      
      minPc <- minmax[1] #min(Pc[,2])
      maxPc <- minmax[2] #max(Pc[,2])
      a1 <- minP0 - maxPc
      a2 <- maxP0 - minPc
      #set.seed(689)
      # method 1: peice-wise linear; method 2:  Beizer curve; method 3: b-spline ; 
        if(method==2){
        lift <- runif(1,a1,a2)
        Pc[,2] <- Pc0[,2]+lift
        y.temp <- apply(x.temp,1,f_bezier,Pc=Pc,deg=deg) #comment this line, when you run consistency checking; 
      }  
      
      #set.seed(Sys.time())
      #x.temp <- as.vector(x.temp)
      
      left <- (sign(y.temp-p.new[,2])==1)
       any(left==TRUE) && any(left==FALSE)
      
      if(any(left==TRUE) && any(left==FALSE)){
      #cat(Pc,"\n")
      #if(TRUE){ # Use this for consistency checking 
        run.index <- 0
        
        V.left.ID <- V.ID[left]
        V.right.ID <- V.ID[!left]
        
        if(n.p>0){

          
          if(method==2){
            points <- f_bezier_arc_length(deg,Pc,n.p)
            xy.crv <- points$curves
            ess <- points$ess
          }          
          xy.crv.af <- rev.mapf(xy.crv,r)
        }

      }else{
        sample.cut.count <- sample.cut.count+1
        run.index <- (1/sample.cut.count>tolerance)*1
        skip.index <- 1-run.index
      }
      }
    
  }
  
  if(skip.index){
    dist.max=0
  }
  
  return(list(curves.bf= xy.crv,
              curves.af= xy.crv.af,
              dist.max=dist.max,
              sample.cut.count=sample.cut.count,
              skip.index=skip.index,
              V.left.ID=V.left.ID,
              V.right.ID=V.right.ID,
              ess = ess))
}

Generative_Process <- function(partition,V.all,V.label,tau,group.level,group.len,tolerance,method,save.cut=FALSE){
  l <- partition$l #l <- length(partition$Polytopes)
  d <- dim(V.all)[2]
  tau.v <- partition$tau
  
  cut.order <- partition$cut.order
 
  partition$logw     <-    log(1/N) # if resampling = TRUE; if  not; revise this sentence.
 
  Cut <- list()
  cut.idx <- 0
  skip.all <- 0  #only set it to 1 when  (1) neither cond1 or cond2 is  true OR (2) cost exceeds the budget .
  group.level2 <- as.factor(group.level)
  
  if(tau.v[l]>=tau){
    skip.all <- 1
  }else{
    Lambdas <- sapply(partition$Polytopes, "[[", 2)
    filters.all <- sapply(partition$Polytopes, "[[", 5)
    
    ##shrink the candidate space, only choose from polytopes whose
    #  (i) Lambda>0
    #       AND
    #  (ii) obs are not from same group (the obs with known labels)
    cond1 <- (Lambdas>0) & (!filters.all)
    cond2 <- (apply(is.infinite(sapply(partition$Polytopes,"[[",3)/0),2,sum)>1)
    cond12 <- sum(cond1&cond2)
    if(cond12==0){
      skip.all <- 1
    }else{
      sample.space <- which(cond1&cond2)
      if(length(sample.space)==1){
        j <- sample.space
      }else{
        j <- sample(x=c(sample.space),size=1,prob = Lambdas[sample.space])
      }
      
      V.temp.ID <- partition$Polytopes[[j]]$V.ID
      V.temp <- V.all[V.temp.ID,]
      V.temp.label <- V.label[V.temp.ID]
      dist.max <- partition$Polytopes[[j]]$dist.max
      filter.ind <- partition$Polytopes[[j]]$filter.ind
      Cut <- Cut_plane_standard(V.temp,V.temp.ID,V.temp.label,tolerance,dist.max,filter.ind,method)
      
      if(Cut$skip.index!=1){
        
        cut.order <- c(cut.order,j)
        partition$cut.order <- cut.order
        
        V.temp.ID.left <- Cut$V.left.ID
        V.temp.ID.right <- Cut$V.right.ID
        
        if(length(V.temp.ID.left)==1){
          V.temp.left <- t(as.matrix(V.all[V.temp.ID.left,]))
        }else{
          V.temp.left <- V.all[V.temp.ID.left,]
        }
        
        if(length(V.temp.ID.right)==1){
          V.temp.right <- t(as.matrix(V.all[V.temp.ID.right,]))
        }else{
          V.temp.right <- V.all[V.temp.ID.right,]
        }
 
        if(save.cut){
          partition$Cut[[l]] <- Cut$curves.af  #change setting here to save more/or other Cut related info.
        }
        
        cut.idx <- 1
        l <- l+1
        
        
        
        V.lbl.left <- V.label[V.temp.ID.left]
        V.lbl.right <- V.label[V.temp.ID.right]
        
      
        #if(sum(!is.na(V.lbl.left))>0 & sum(is.na(V.lbl.left))>0 & (sum(table(V.lbl.left)!=0)>1)){
         if(sum(!is.na(V.lbl.left))>0 & (sum(table(V.lbl.left)!=0)>1)){
          
          id.train.left <-  which(id.train %in% V.temp.ID.left )
          dist.max.left <- max(dist.all[id.train.left,id.train.left])
          
          filter.ind.left<-FALSE
        }else{
          dist.max.left <- 0
          filter.ind.left<- TRUE
        }
        
        
        
        #if(sum(!is.na(V.lbl.right))>0 & sum(is.na(V.lbl.right))>0 & (sum(table(V.lbl.right)!=0)>1)){
        if(sum(!is.na(V.lbl.right))>0 & (sum(table(V.lbl.right)!=0)>1)){
          #dist.max.right <- max(dist.all[V.temp.ID.right,V.temp.ID.right])
          id.train.right <-  which(id.train %in% V.temp.ID.right  )
          dist.max.right <- max(dist.all[id.train.right,id.train.right])
          filter.ind.right<-FALSE
        }else{
          dist.max.right <- 0
          filter.ind.right <- TRUE
        }
        
        
        Lambda.left <- dist.max.left/2
        Lambda.right <- dist.max.right/2

 
        countbygroup.left <- as.numeric(table(c(V.label[V.temp.ID.left],group.level2))-1)
        countbygroup.right <- as.numeric(table(c(V.label[V.temp.ID.right],group.level2))-1)
        
        cbg.lr <- countbygroup.left + countbygroup.right + alpha0
        cbg.l <- countbygroup.left +  alpha0
        cbg.r <- countbygroup.right +  alpha0
      
        logwt <-  sum(lgamma(cbg.l))-sum(lgamma(sum(cbg.l)))+sum(lgamma(cbg.r))-sum(lgamma(sum(cbg.r))) -sum(lgamma(cbg.lr))+sum(lgamma(sum(cbg.lr)))-sum(lgamma(alpha0))+lgamma(sum(alpha0)) 
       
        partition$logw     <-   logwt+ partition$logw
        
        partition$Polytopes[[j]] <- list(V.ID=V.temp.ID.left,Lambda=Lambda.left,CountByGroup=countbygroup.left,dist.max=dist.max.left,filter.ind=filter.ind.left)
        partition$Polytopes[[l]] <- list(V.ID=V.temp.ID.right,Lambda=Lambda.right,CountByGroup=countbygroup.right,dist.max=dist.max.right,filter.ind=filter.ind.right)
        partition$ess <- c(partition$ess,Cut$ess)
          
        Lambdas <- sapply(partition$Polytopes, "[[", 2)
        sum_Lbds <- sum(Lambdas)
        if(sum_Lbds>0){
          tau.v.plus1 <- rexp(n=1,sum_Lbds)+tau.v[l-1]
          tau.v <- c(tau.v,tau.v.plus1)
          partition$tau <- tau.v
        }else{
          tau.v.plus1 <- Inf
          tau.v <- c(tau.v,tau.v.plus1)
          partition$tau <- tau.v
          skip.all <- 1
        }
        
        
      }
      
      if(Cut$skip.index==1){
        partition$Polytopes[[j]]$dist.max <- 0
        partition$Polytopes[[j]]$filter.ind<- TRUE
        partition$Polytopes[[j]]$Lambda <- 0
      }
      
      
    }
    
  }
  
  partition$cut.indx   <- cut.idx 
  partition$skip.index <- skip.all
  partition$l          <- l


  return(partition)
}

partition_initial <- function(group.level,group.len,V.all,V.label,l.max,tau,tolerance){
  group.level2 <- as.factor(group.level)
  n <- length(V.label)
  V.temp.ID <- c(1:n)
  V.temp <- V.all[V.temp.ID,]
  d <- dim(V.all)[2]
  Polytopes <- list()
  #dist.all <- as.matrix(dist(V.all))
  
  dist.all <- as.matrix(dist(V.all[id.train,]))
  
  countbygroup <- as.numeric(table(c(V.label,group.level2))-1)
  
  dist.max <- max(dist.all)
  Polytopes[[1]] <- list(V.ID=V.temp.ID,  Lambda=dist.max/2, CountByGroup=countbygroup,dist.max=dist.max, filter.ind=FALSE)
  
  
  Partition.Inital <- list(Polytopes=Polytopes)
  
  output <- list(Partition.Inital=Partition.Inital,dist.all=dist.all)
  return(output)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


weight_f <-function(Partition.t,logl.tminus1,logl.t,logWt,lent1){
  lent2 <- sapply(Partition.t,"[[","l")
  diff_indx <- which((lent2-lent1)!=0) 
  if(length(diff_indx)>0){
    for(j in diff_indx){
      cbg.new <- sapply(Partition.t[[j]]$Polytopes, "[[", 3) # count by group for the new partion
      multi.beta.a <- apply(cbg.new,2, (function(x){x+alpha0}))
      logl.t[j] <- sum(lgamma(multi.beta.a))-sum(lgamma(apply(multi.beta.a,2,sum)))-dim(cbg.new)[2]*(sum(lgamma(alpha0))-lgamma(sum(alpha0)))
      logWt[j] <- logWt[j]+logl.t[j]-logl.tminus1[j]
      logl.tminus1[j] <- logl.t[j]
    }
    
  }
  return(out=list(logl.tminus1=logl.tminus1,logl.t=logl.t,logWt=logWt))
}
