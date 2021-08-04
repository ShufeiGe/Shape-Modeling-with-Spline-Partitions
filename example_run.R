

tau <- 5
file.out=paste0("./result/out","_cuts_tau",tau,"_v1.RData")
file.out2=paste0("./result/out","_cuts_tau",tau,"_v2.RData")
library(purrr)
library(parallel)
source("curves_f.R")
source("partition_function_v2.R")

method <-2
no.cores <-1
l.max <- Inf  # no. of cuts
N <- 5 # no. particles in each repetition
exp.rep <- 1  # number of repetitions of the experiments

tolerance <- 0.01

load("yinyang_test.RData")

 
table(data.train[,3])


d     <- dim(data.train)[2]-1
V.all <-  as.matrix(data.train[,1:d])



ALPHA <- 1/1000

split.seed <-888
set.seed(split.seed)
seed_v <- sample(x=1:10000000,size=exp.rep)

rep.max <- 1  # number of trees; fix it to 1; 
Wt.m.all <-rep(NA,length = exp.rep)
acc.all <-rep(NA,length=exp.rep)
p.selected <- NULL
n.p <- 200       # length.out <- n.p # number of points on the proposed line; #if possible output curves.


if(l.max==Inf){
  l.max=dim(data.train)[1]
}

min <- matrix(apply(V.all,2,min),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
max <- matrix(apply(V.all,2,max),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
V.all <- (V.all-min)/(max-min)
rm("min","max")

apply(V.all,2,range)
V.all <- (V.all-0.5)
apply(V.all,2,range)

V.label <- data.train[,(d+1)]
table(V.label)
id.test <- which(is.na(V.label))

if(length(id.test)==0){
  id.train <- c(1:length(V.label))
}else{
  id.train <- c(1:length(V.label))[-id.test]
}

#V.label[id.test]<- NA #set the unknow labels as NA; otherwise, they will be grouped as a new group.
group <- as.factor(V.label)

group.level <- levels(group)
group.len <- length(group.level)



part_init <- partition_initial(group.level,group.len,V.all,V.label,l.max,tau,tolerance)
Partition.Inital <- part_init$Partition.Inital
dist.all <- part_init$dist.all
rm("part_init")
pb <- txtProgressBar(min = 0, max = exp.rep*rep.max, style = 3)

t0 <- proc.time()
for(exp.rep.i in 1:exp.rep){
  set.seed(seed_v[exp.rep.i])
  #----inference -----start from here--------
  alpha0 <- matrix((round(table(group)/table(group)[1],2))*ALPHA,ncol=1) #try math paramter
  
  table.out <- c()
  
  Partition.t0 <- list()
  Lambda0 <- Partition.Inital$Polytopes[[1]]$Lambda
  
  cbg0 <- sapply(Partition.Inital$Polytopes, "[[", 3)
  multi.beta.0 <- apply(cbg0,2, (function(x){x+alpha0}))
  logl.tminus0 <- rep(sum(lgamma(multi.beta.0))-sum(lgamma(apply(multi.beta.0,2,sum)))-dim(cbg0)[2]*(sum(lgamma(alpha0))-lgamma(sum(alpha0))),N)
  
  
  Cut <- list()
  Cut[[1]] <- cbind(rep(NA,n.p),rep(NA,n.p))
  for(i in 1:N){
    tau.temp <- rexp(1,rate=Lambda0)
     Partition.t0[[i]] <- list(Polytopes=Partition.Inital$Polytopes,tau=tau.temp,Cut=Cut,cut.indx=0,skip.index=0,l=1,cut.order=NULL,logw=logl.tminus0[1],ess=NULL)
    
  }
  #rm("Cut")
  
  Wt0 <- rep(1/N,N)   # Weights vector for N particles at time t=0
  logWt0 <- log(Wt0)   # log Weights vector for N particles at time t=0
  
  label.predict.all <- NULL
  l.cut.all <- NULL
  tau.li0 <- min(map_dbl(Partition.t0,~rev(.$tau)[1]))
  for(rep in 1:rep.max){
    setTxtProgressBar(pb, ((exp.rep.i-1)*rep))
    
    Partition.t <- Partition.t0
    Wt <- Wt0        # Weights vector for N particles at time t
    logWt <- logWt0  # log Weights vector for N particles at time t
    logl.tminus1 <- logl.tminus0 #log likelihood for N particles at time t-1
    logl.t <- logl.tminus1 #log likelihood for N particles at time t
    # create progress bar
    tau.li <- tau.li0
    end.cond <- FALSE
    
    skip.index.N <- rep(0,N)
    
    t <- 0  # here t is the number of cuts; we use vector tau to denote the cost at each cut.
    t.idx <- rep(0,l.max)
    while(tau.li<tau){
      #resampling-------------start
      if(t>0){
        
        jC <- sample(x=1:N,size=N,prob =Wt,replace = TRUE )
        
        Partition.t.new <- list()
        for(i in 1:N){
          Partition.t.new[[i]] <- Partition.t[[jC[i]]]
        }
        
        Partition.t <- Partition.t.new
        rm("Partition.t.new")
      }
      
      Partition.t <-mclapply(Partition.t,function(partition.temp){
        if(!partition.temp$skip.index){
          partition.temp <- Generative_Process(partition.temp,V.all,V.label,tau,group.level,group.len,tolerance,method,save.cut=TRUE)
        }
        return(partition.temp)
      }, mc.cores = no.cores)
      
      logWt        <-  sapply(Partition.t,"[[","logw")
      Wt <- exp((logWt-max(logWt)))/sum(exp((logWt-max(logWt))))
      t <- max(length(Partition.t[[which.max(Wt)]]$tau)-1,1)
      tau.li <- min(map_dbl(Partition.t,~rev(.$tau)[1]))
      end.cond <- min(skip.index.N)
      
      t.idx[t] <- t.idx[t]+1
      if(end.cond | (t+1)>l.max | tau.li>tau|(t.idx[t]>100) ){
        tau.li <- Inf
      }
      
      count.by.group.temp <- sapply(Partition.t[[which.max(Wt)]]$Polytopes, "[[", 3)
      count.by.group.temp.idx <- apply(count.by.group.temp,2,which.max)
      ID.temp <- sapply(Partition.t[[which.max(Wt)]]$Polytopes, "[[", 1,simplify = FALSE)
      label.predic <- rep(NA,length(V.label))
      for(id.j in 1:length(ID.temp)){
        label.predic[ID.temp[[id.j]]] <- as.numeric(group.level[count.by.group.temp.idx[id.j]])
      }
      label.predict.all<- cbind(label.predict.all,label.predic)
      l.cut.all <- rbind(l.cut.all,t)
      #}
      if(t%%100==0){
        base::save(Partition.t,l.cut.all,label.predict.all,V.all,V.label,file=file.out)
      }
      
    }
  }
  
  
  Wt.m <- which.max(Wt)
  Wt.m.all[exp.rep.i]  <- Wt.m
  
  zz0 <- (proc.time()-t0)
  if((exp.rep.i%%exp.rep)==0){
    zz = round(zz0[3]/60/60,4)
    if (zz < 1) {
      zz = round(zz * 60)
      if (zz <= 1) {
        zz = "1 minute"
      } else {
        zz = sprintf("%d minutes", zz)
      }
    } else if (zz > 24) {
      zz = round(zz/24)
      if (zz <= 1) {
        zz = "1 day"
      } else {
        zz = sprintf("%d days", zz)
      }
    } else {
      zz = round(zz)
      if (zz <= 1) {
        zz = "1 hour"
      } else {
        zz = sprintf("%d hours", zz)
      }
    }
    cat(sprintf(paste("\n"," Running time: ", zz, sep="" )))
  }
}

rm("dist.all")
base::save.image(file=file.out)
zz0 <- (proc.time()-t0)

j <- which.max(Wt)
Partition.t.Wt.max=Partition.t[[j]]
base::save(Partition.t.Wt.max,V.all,V.label,file=file.out2)
 
