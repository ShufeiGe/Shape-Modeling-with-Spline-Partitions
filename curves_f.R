# functions : (1) mapf() map points to the new coordinate system; (2) rev.mapf() : map points to the original coordinate system.
# the new coordinate system is defined by centre (c),rotation (r),translation (t).
# Input:
#   P: points to map, a n by 2 matrix; per row per point;
#   r: angle of rotation takes value from 0 to 2*pi.
# Output:
#   a 2 by n matrix reprents the points in the new coordinate system.

mapf <- function(P,r){
  n <- nrow(P)
  a <- -r 
  # M =[cos(.), -sin();sin(.),cos(.)]; 
  # rotating coordinate system by an angle r is equivalently to rotating points by -r
  r.M <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P%*%t(r.M)
  return(P.new)
}
rev.mapf <- function(P,r){
  n <- nrow(P)
  a <- -r 
  r.M.inv <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),nrow=2,ncol=2,byrow = T)
  P.new <- P%*%t(r.M.inv)
  return(P.new)
}

 
# f_bezier():  given control points and x-coordinate(s) of interest, return the crossponding y coordinate(s).
# Input: Pc, control points of the bezier curve, a m by 2 matrix, per row per point, 
#        and the values in column 1 are sorted in ascending order.
#        deg, degree of the bezier curve, deg=m-1
#        x, predictor.
# Output: f(x)
 
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
 
# f_bezier_minmax() :  given control points, return the lower and upper bounds of the  bezier curve.
# Input: Pc, control points of the bezier curve, a m by 2 matrix, per row per point, 
#        and the values in column 1 are sorted in ascending order.
#       deg: degree of the bezier curve, deg=m-1
# Output: lower and upper bounds.


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
    
    
  }
  
  
  if(deg==3){
    cc <- NULL
    Mp <- matrix(c(0,0,0,1,0,0,0,2,0,0,0,3),3,4)
    M <- matrix(c(1,0,0,0,-3,3,0,0,3,-6,3,0,-1,3,-3,1),4,4,byrow=TRUE)
    
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
