# function1: map points to the new coordinate system.
# the new coordinate system is defined by centre (c),rotation (r),translation (t).
# Input:
#   P: points to map, a n by 2 matrix; per row per point;
#   c: centre, a 1 by 2 matrix or a vector of length 2.
#   r: angle of rotation takes value from 0 to 2*pi.
#   t: translation, a scalar, takes value from (-infinity, infinity)
# Output:
#   a 2 by n matrix reprents the points in the new coordinate system.
#rm(list=ls())

runif_polar_rho <- function(n,R=1){
  u <- runif(n)
  r <- sqrt(u)*R
  return(r)
}


mapf <- function(P,c, r, t){
  n <- nrow(P)
  P.new <- P - matrix(c,nrow = n, ncol=2,byrow=TRUE)
  a <- -r 
  # M =[cos(.), -sin();sin(.),cos(.)]; 
  # rotating coordinate system by an angle r is equivalently to rotating points by -r
  r.M <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P.new%*%t(r.M)
  P.new[,2] <- P.new[,2]-t
 return(P.new)
}

mapf_v2 <- function(P,c, r, t1,t2){
  n <- nrow(P)
  P.new <- P - matrix(c,nrow = n, ncol=2,byrow=TRUE)
  a <- -r 
  # M =[cos(.), -sin();sin(.),cos(.)]; 
  # rotating coordinate system by an angle r is equivalently to rotating points by -r
  r.M <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P.new%*%t(r.M)
  P.new[,2] <- P.new[,2]-t2
  P.new[,1] <- P.new[,1]-t1
  return(P.new)
}

mapf_v3 <- function(P,r){
  n <- nrow(P)
  a <- -r 
  # M =[cos(.), -sin();sin(.),cos(.)]; 
  # rotating coordinate system by an angle r is equivalently to rotating points by -r
  r.M <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P%*%t(r.M)
  return(P.new)
}
rev.mapf_v3 <- function(P,r){
  n <- nrow(P)
  a <- -r 
  r.M.inv <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P%*%t(r.M.inv)
  return(P.new)
}

# rotate points back; reverse of mapf()
# M^{-1}=[cos(.),sin(.); -sin(.),cos(.)]
rev.mapf <- function(P,c,r,t){
  n <- nrow(P)
  a <- -r 
  P.new <- P
  P.new[,2] <- P.new[,2]+t
  r.M.inv <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P.new%*%t(r.M.inv)
  P.new <- P.new + matrix(c,nrow = n, ncol=2,byrow=TRUE)
  return(P.new)
}

rev.mapf_v2 <- function(P,c,r,t1,t2){
  n <- nrow(P)
  a <- -r 
  P.new <- P
  P.new[,1] <- P.new[,1]+t1
  P.new[,2] <- P.new[,2]+t2
  r.M.inv <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P.new%*%t(r.M.inv)
  P.new <- P.new + matrix(c,nrow = n, ncol=2,byrow=TRUE)
  return(P.new)
}



#set.seed(123)
#n <- 100
#P <- cbind(9,9)
##P <- cbind(runif(n,-10,10),runif(n,-10,10))
#c <- c(10,10)
#r <- pi/4
#t <- 2


#P.new <- mapf(P,c,r,t)
#par(mfrow=c(2,2))
#plot(P,col=1,ylim=c(-20,20),xlim=c(-20,20))
#lines(P.new,col=2,type="p")
#(P.new)

# function2: peicewise linear; return f(x)
# Input: Pc, knots on the peicewise linear,a m by 2 matrix, per row per point, 
#        and the values in column 1 are sorted in ascending order.
#        x, predictor.
# Output: f(x)
pl <- function(x,Pc){

  cc <- which(Pc[,1]-x==0)
  
  if(length(cc)!=0){
    y <- Pc[cc,2]
  }else{
    m <- nrow(Pc)
    order <- order(Pc[,1],decreasing = FALSE)
    Pc <- Pc[order,]
    
    if(x<Pc[1,1] | x>Pc[m,1]){
      stop("x is beyond the user-defined domain!!!")
    }else{
      Int <-cbind(Pc[1:(m-1),1],Pc[2:m,1])
      sign.1 <- sign(Int[,1]-x)+sign(Int[,2]-x)
      cc <- which(sign.1==0)
      
      p1 <- Pc[cc,]
      p2 <- Pc[cc+1,]
      y <- (p2[2]-p1[2])*(x-p1[1])/(p2[1]-p1[1])+p1[2]
    }
      
    
  }
  
  return(y)
  
}
#Pc <- matrix(runif(8),4,2)
#x <-matrix(runif(5,min=min(Pc[,1]),max=max(Pc[1,])),ncol=1)
#
#pl(x[1],Pc)
#pl(x[2],Pc)
#pl(x[3],Pc)
#pl(x[4],Pc)
#pl(x[4],Pc)
#apply(x,1,pl,Pc=Pc)

# function3: bezier curve (order =1,2,3)
# Input : 
#       x, predictor
#       Pc, control points, (deg+1) by 2 matrix, per row per point
#       deg = 1, 2, 3; 1 for linear; 2 for quadratic, 3 for cubic)
# Output: y; (x,y) represents a point on the bezier curve.
#
#f_bezier <- function(x,Pc,deg){
# dg2 <- nrow(Pc) -1
# if( dg2 !=deg){
#   stop("number of rows of Pc doesn't match the degree")
# }
# 
# if(deg==1){
#   M <- matrix(c(1,0,-1,1),2,2,byrow=TRUE)
#   Mtx<- M%*%Pc[1:2,1]
#   Mtx[1]<-Mtx[1]-x
#   t <- Re(polyroot(Mtx))
#   y <- c(1,t)%*%M%*%Pc[1:2,2]
# } 
#
# if(deg==2){
#   M <- matrix(c(1,0,0,-2,2,0,1,-2,1),3,3,byrow=TRUE)
#   Mtx<- M%*%Pc[1:3,1]
#   Mtx[1]<-Mtx[1]-x
# 
#   roots <- polyroot(Mtx)
#   r.r   <- Re(roots)
#   r1    <- (r.r>(1+10e-10))
#   r2    <- Mod(roots)-r.r
#   rr    <- r1+r2
#   ii    <- which(rr==0)
#   
#   if(length(ii)!=1){
#     stop("out of scope")
#   }else{
#     t <- r.r[ii]
#     y <- c(1,t,t^2)%*%M%*%Pc[1:3,2] 
#   }
#   
#   
# }
# 
# if(deg==3){
#   M <- matrix(c(1,0,0,0,-3,3,0,0,3,-6,3,0,-1,3,-3,1),4,4,byrow=TRUE)
#   Mtx<- M%*%Pc[1:4,1]
#   Mtx[1]<-Mtx[1]-x
#
#
#   roots <- polyroot(Mtx)
#   r.r   <- Re(roots)
#   r1    <- (r.r>(1+10e-10))
#   r2    <- Mod(roots)-r.r
#   rr    <- r1+r2
#   ii    <- which(rr==0)
#   
#   if(length(ii)!=1){
#     stop("out of scope")
#   }else{
#     t <- r.r[ii]
#     y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2] 
#   }
#  
#   
#
#   
# }
# 
# return(y)
#}



f_bezier <- function(x,Pc,deg){
  dg2 <- nrow(Pc) -1
  if( dg2 !=deg){
    stop("number of rows of Pc doesn't match the degree")
  }
  
  if(deg==1){
    M <- matrix(c(1,0,-1,1),2,2,byrow=TRUE)
    Mtx<- M%*%Pc[1:2,1]
    Mtx[1]<-Mtx[1]-x
    t <- Re(polyroot(Mtx))
    y <- c(1,t)%*%M%*%Pc[1:2,2]
  }
  
  if(deg==2){
    M <- matrix(c(1,0,0,-2,2,0,1,-2,1),3,3,byrow=TRUE)
    Mtx<- M%*%Pc[1:3,1]
    Mtx[1]<-Mtx[1]-x
    
    roots <- polyroot(Mtx)
    r.r <- Re(roots)
    r.r1 <- (-10e-3 < r.r) & (r.r< 1+10e-3)
    if(sum(r.r1)==1){
      t <- r.r[r.r1]
      y <- c(1,t,t^2)%*%M%*%Pc[1:3,2]
    }else if(sum(r.r1)==0){
      stop("a bug in f_bezier()-1!")
    }else{
      r.i   <- (abs(Im(roots))<=10e-3) & r.r1
      if(sum(r.i)==0){
        stop("a bug in f_bezier()-2!")
      }else{
        r.r2 <- r.r[r.i]
        t <- min(max(r.r2,0),1)
        y <- c(1,t,t^2)%*%M%*%Pc[1:3,2]
      }
    }
  }
  
  
  if(deg==3){
    M <- matrix(c(1,0,0,0,-3,3,0,0,3,-6,3,0,-1,3,-3,1),4,4,byrow=TRUE)
    Mtx<- M%*%Pc[1:4,1]
    Mtx[1]<-Mtx[1]-x

    roots <- polyroot(Mtx)
    r.r <- Re(roots)
    r.r1 <- (-10e-3 < r.r) & (r.r< 1+10e-3)
    if(sum(r.r1)==1){
      t <- r.r[r.r1]
      y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2]
    }else if(sum(r.r1)==0){
      stop("a bug in f_bezier()-1!")
    }else{
      r.i   <- (abs(Im(roots))<=10e-3) & r.r1
      if(sum(r.i)==0){
        stop("a bug in f_bezier()-2!")
      }else{
        r.r2 <- r.r[r.i]
        t <- min(max(r.r2,0),1)
        y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2]
      }
    }
    
  }
  
  return(y)
}

#f_bezier(x,Pc,deg)
#roots <- polyroot(Mtx)
#r.r <- Re(roots)
#r.r1 <- (-10e-3 < r.r) & (r.r< 1+10e-3)
#if(sum(r.r1)==1){
#  t <- r.r[r.r1]
#  y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2]
#}else if(sum(r.r1)==0){
#  stop("bug in f_bezier()-1!")
#}else{
#  r.i   <- (abs(Im(roots))<=10e-3) & r.r1
#  if(sum(r.i)==0){
#    stop("bug in f_bezier()-2!")
#  }else{
#    r.r2 <- r.r[r.i]
#    t <- min(max(r.r2,0),1)
#    y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2]
#  }
#}


#roots <- polyroot(Mtx)
#r.r <- Re(roots)
#r.i   <- which((abs(Im(roots))<=10e-3))
#if(length(r.i)==0){
#  stop("bug in f_bezier()-1!")
#}else{
#  r.r1 <- (-10e-3 < r.r[r.i]) & (r.r[r.i] < 1+10e-3)
#  if(length(r.r1)==0){
#    stop("bug in f_bezier()-2!")
#  }else{
#    r.r2 <- r.r[r.r1]
#    t <- min(max(r.r2,0),1)
#    y <- c(1,t,t^2,t^3)%*%M%*%Pc[1:4,2]
#  }
#}



#Pc <- matrix(runif(8),3,2)
#Pc.ord <-order(Pc[,1])
#Pc <-Pc[Pc.ord,]
#deg <- 2
#x <- as.matrix(runif(100, min=min(Pc[,1]),max=max(Pc[,1])))
#f_bezier(x[1],Pc,deg)
#f_bezier(x[2],Pc,deg)
#f_bezier(x[3],Pc,deg)
#f_bezier(x[4],Pc,deg)
#y<-apply(x,1,f_bezier,Pc=Pc,deg=deg)
#x <- as.vector(x)
#ox <- order(x)
#par(mfrow=c(1,2))
#plot(x[ox],y[ox],type="l")
# t <- seq(0, 1, length=100)
#plot(bezier(t=t, p=Pc, deg=deg),type="l",xlab = "",ylab="")



f_bezier_minmax <- function(Pc,deg){
  dg2 <- nrow(Pc) -1
  minmax <- c(NA,NA)
  cc <- NULL
  if( dg2 !=deg){
    stop("number of rows of Pc doesn't match the degree")
  }
  
  if(deg==1){
    M <- matrix(c(1,0,-1,1),2,2,byrow=TRUE)
    minmax <- sort(c(c(1,0)%*%M%*%Pc[1:2,2],c(1,1)%*%M%*%Pc[1:2,2]))
    #Pc <- cbind(sort(c(0,1)),runif(2))
    #M <- matrix(c(1,0,-1,1),2,2,byrow=TRUE)
    #t <- seq(0,1,length.out = 20)
    #Mtx <- M%*%Pc[,1]
    #Mty <- M%*%Pc[,2]
    #x <- apply(as.matrix(t),1,function(x){sum((c(1,x))*Mtx)})
    #y <- apply(as.matrix(t),1,function(x){sum((c(1,x))*Mty)})
    #plot(x,y,type="l")
    #f_bezier_minmax(Pc,deg=1)
  }
  
  if(deg==2){
    M <- matrix(c(1,0,0,-2,2,0,1,-2,1),3,3,byrow=TRUE)
    Mp <- matrix(c(0,1,0,0,0,2),2,3,byrow=TRUE)
    Mty<- Mp%*%M%*%Pc[1:3,2]
    z <- Re(polyroot(Mty))
    if(z)
      if(0<z & z<1){
        cc <- c(1,z,z^2)%*%M%*%Pc[1:3,2]
      }
    
    cc0 <- c(1,0,0)%*%M%*%Pc[1:3,2]
    cc1 <- c(1,1,1)%*%M%*%Pc[1:3,2]
    cc01 <- c(cc,cc0,cc1)
    minmax <- c(min(cc01),max(cc01))
    #Pc <- cbind(sort(c(0,1,runif(1))),runif(3))
    #t <- seq(0,1,length.out = 20)
    # Mtx <- M%*%Pc[1:3,1]
    # Mty <- M%*%Pc[1:3,2]
    #x <- apply(as.matrix(t),1,function(x){sum((c(1,x,x^2))*Mtx)})
    #y <- apply(as.matrix(t),1,function(x){sum((c(1,x,x^2))*Mty)})
    #plot(x,y,type="l")
    
  }
  
  
  if(deg==3){
    cc <- NULL
    Mp <- matrix(c(0,0,0,1,0,0,0,2,0,0,0,3),3,4)
    M <- matrix(c(1,0,0,0,-3,3,0,0,3,-6,3,0,-1,3,-3,1),4,4,byrow=TRUE)
    
    #Pc <- cbind(sort(c(0,1,runif(2))),runif(4))
    #t <- seq(0,1,length.out = 20)
    #Mtx <- M%*%Pc[,1]
    #Mty <- M%*%Pc[,2]
    #x <- apply(as.matrix(t),1,function(x){sum((c(1,x,x^2,x^3))*Mtx)})
    #y <- apply(as.matrix(t),1,function(x){sum((c(1,x,x^2,x^3))*Mty)})
    #plot(x,y,type="l")
    
    Mty<- Mp%*%M%*%Pc[,2]
    z <- polyroot(Mty)
    re.z <- Re(z)
    z1 <- re.z[1]
    z2 <- re.z[2]
    if(sum(Mod(z)==re.z)==2){
      if(z1>0 & z1<1){
        cc <- c(cc,(c(1,z1,z1^2,z1^3)%*%M%*%Pc[,2]))
      }
      if(z2>0 & z2<1){
        cc <- c(cc,(c(1,z2,z2^2,z2^3)%*%M%*%Pc[,2]))
      }
    }
    cc0 <- c(1,0,0,0)%*%M%*%Pc[,2]
    cc1 <- c(1,1,1,1)%*%M%*%Pc[,2]
    cc01 <- c(cc,cc0,cc1)
    minmax <- c(min(cc01),max(cc01))
  }
  return(minmax)
}

#function4: distance of point P=(x,y) to a line segment 
#Input: p2ls(P0,P1,P2)
#      point P=(x,y)
#      line segment (defined by two points), P1 and P2
p2ls <- function(p0,p1,p2){
  p10 <- p1-p0
  p12 <- p1-p2
  p20 <- p2-p0
  
  p10_p12 <- sum(p10*p12)
  p20_p21 <- sum(p20*(-p12))
  
  if(p10_p12<=0){
    d <- sqrt(sum(p10^2))
  }else if(p20_p21<=0){
    d <- sqrt(sum(p20^2))
  }else{
    cp <- sum(p20*(-p12))/(sum(p12^2))
    vp <- p20 - cp*(-p12)
    d  <- sqrt(sum(vp^2))
  }
  return(d)
}

#p1=c(0,0)
#p2=c(3,0)
#p0=c(4,4)
#p2ls(p0,p1,p2)
#p2ls(c(3,3),p1,p2)
#p2ls(c(2,2),p1,p2)
#p2ls(c(-4,-4),p1,p2)
#p2ls(c(2,-4),p1,p2)

#generate a random b-spline function
bsf <- function(order,nknots,knots,sigma){
  require(fda)
  #knots    <- seq(xlim[1], xlim[2], length.out = nknots)
  #nbasis   = nknots + order - 2
  #basis = create.bspline.basis(xlim,nbasis,order,knots)
  basis  <-  create.bspline.basis(norder=order,breaks=knots)
  nbasis <-  nknots + order - 2
  coeffs <-  rnorm(nbasis,0,sd=sigma)
  return(list(basis=basis,coeffs=coeffs))
}

f_bezier_arc_length <- function(deg,Pc,n){
  n2 <- n
 # n <- 1000
  if(deg==1){
    M <- matrix(c(1,0,-1,1),2,2,byrow=TRUE)
    Mtx <- M%*%Pc[,1]
    Mty <- M%*%Pc[,2]
    tt <- runif(n)
    tt <- sort(tt)
    w_n <- rep(1/n,n) 
    (Ess <- 1/sum(w_n^2))
    #id0<- sample(1:n,size=n,prob = w_n,replace = TRUE)
    #id<- sample(id0,size=n2)
    
    x1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x))*Mtx)})
    y1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x))*Mty)})
  }else if(deg==2){
    M <- matrix(c(1,0,0,-2,2,0,1,-2,1),3,3,byrow=TRUE)
    Mtx <- M%*%Pc[,1]
    Mty <- M%*%Pc[,2]
    tt <- runif(n)#seq(0,1,length.out = 50)
    tt <- sort(tt)
    w <- apply(as.matrix(tt),1,function(x){sqrt(sum((c(0,1,2*x)%*%M%*%Pc)^2))})
    w_n <- w/sum(w)
    (Ess <- 1/sum(w_n^2))
     
    
    #id0<- sample(1:n,size=n,prob = w_n,replace = TRUE)
    #id<- sample(id0,size=n2)
    x1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x,x^2))*Mtx)})
    y1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x,x^2))*Mty)})
  }else if(deg==3){
    M <- matrix(c(1,0,0,0,-3,3,0,0,3,-6,3,0,-1,3,-3,1),4,4,byrow=TRUE)
    Mtx <- M%*%Pc[,1]
    Mty <- M%*%Pc[,2]
    tt <- runif(n)#seq(0,1,length.out = 50)
    tt <- sort(tt)
    w <- apply(as.matrix(tt),1,function(x){sqrt(sum((c(0,1,2*x,3*x^2)%*%M%*%Pc)^2))})
    w_n <- w/sum(w)
    (Ess <- 1/sum(w_n^2))
 
    #id0<- sample(1:n,size=n,prob = w_n,replace = TRUE)
    #id<- sample(id0,size=n2)
  
    x1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x,x^2,x^3))*Mtx)})
    y1 <- apply(as.matrix(tt),1,function(x){sum((c(1,x,x^2,x^3))*Mty)})
  }
  #return(list(curves = cbind(x1[id],y1[id]), ess = Ess))
  return(list(curves = cbind(x1,y1), ess = Ess))
}

#n <- 50
#deg <- 1
#Pc <- cbind(c(0,1),runif(2))
#points <- f_bezier_arc_length(deg,Pc,n)
#points$ess
#plot(points$curves,type="l")  
#lines(points$curves,type="p",col=2)  
#
#deg <- 2
#Pc <- cbind(c(0,0.3,1),runif(3))
#points <- f_bezier_arc_length(deg,Pc,n)
#points$ess
#id2 <- order(points$curves[,1])
#plot(points$curves[id2,],type="l")  
#lines(points$curves,type="p",col=2)  
#
#
#deg <-3
#Pc <- cbind(c(0,0.3,0.5,1),runif(4))
#
#points <- f_bezier_arc_length(deg,Pc,n)
#points$ess
#id2 <- order(points$curves[,1])
#plot(points$curves[id2,],type="l")  
#lines(points$curves,type="p",col=2)  
#

#order=4
#nknots=5
#knots <- seq(0,1,length.out = nknots)
#xlim=c(0,1)
#sigma=100
#bsf.test <- bsf(order,nknots,knots,sigma)
#x <- seq(0,1,length.out = 100)
#y <- apply((eval.basis(x,bsf.test$basis)%*%bsf.test$coeffs),1,sum)
#plot(x,y,type="l",col=2)

