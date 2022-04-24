setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 
library(tidyr)

dirname <- paste0(getwd(),"/stn80.csv")
stn80 <- fread(dirname)
available_stations <- fread(paste0(getwd(),"/", "available_stations_3.csv"))


av <- merge(available_stations,stn80[,1:4], all=FALSE)
av <- av[,c(1,5,6,7)]
#plot(unlist(aba[,5]),unlist(aba[,6]))


eucl_dist <- function(x1,x2){
  as.numeric(sqrt((x2[1]-x1[1])^2+(x2[2] - x1[2])^2))
}


# Haversine formula (https://en.wikipedia.org/wiki/Haversine_formula): 
# it can be implemented also with:   
# library(geosphere)
# distHaversine(c(53.36575, 7.348687), c(53.36507, 7.348940))
# I left it because we are going to do an estimation of the distance.
# First of all, the haversine formula compute the exact distance called "crow-flight"
# which is lat and long without altitude.
# Then, we will apply Pitagora between the difference in the altitude and this 
# distance.
crow_flight_dist <- function(loc1,loc2) {
  rad <- pi/180
  a1 <- loc1[2]*rad
  a2 <- loc1[1]*rad
  b1 <- loc2[2]*rad
  b2 <- loc2[1]*rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1)*cos(b1)*(sin(dlon/2))^2
  c <- 2*atan2(sqrt(a), sqrt(1 - a))
  R <- 6378137   # metri
  d <- R*c
  return(d)
}



# this is for testing. The result should be  80.18433
# crow_flight_dist(c(53.36575, 7.348687), c(53.36507, 7.348940))
est_dist <- function(loc1,loc2){
  d <-  sqrt(crow_flight_dist(loc1, loc2)**2 + (loc1[3] - loc2[3])**2)
  return(d)
}

# testing for the estimated distance:
# locc1 <- unlist(c(39.7772,-98.7783,554.7))
# locc2 <- unlist(c(32.1300,-81.2100,14.0))
# crow_flight_dist(locc1,locc2)
# est_dist(locc1,locc2)




# matrix with euclidian distances
# you need to pass only three columns
matrix_dist <- function(set){
  
  a <- matrix(0, nrow = dim(set)[1], ncol = dim(set)[1])
  for (j in 1:dim(set)[1]) { 
    for (i in 1:dim(set)[1]) { 
      a[i,j] <- eucl_dist(unlist(set[j,]),unlist(set[i,]))
    }
  }
  return(a)
}

set1 <- matrix_dist(av[,c(2,3)])



# lower triangular matrix with estimated distances
# you need to pass only three columns
matrix_est_dist <- function(set){
  
  a <- matrix(0, nrow = dim(set)[1], ncol = dim(set)[1])
  for (j in 1:(dim(set)[1] - 1)) { 
    for (i in (j+1):dim(set)[1]) { 
      a[i,j] <- est_dist(unlist(set[j,]),unlist(set[i,]))
    }
  }
  return(a)
}
set2 <- matrix_est_dist(av[,c(2,3,4)])
set2 <- set2 / 1000   #in Km
a <- set2
a[a <= 0] <- Inf




# Agglomerative clustering 
# single-linkage -> When a new cluster is formed, 
# the (dis)similarities between it and the other clusters and/or 
# individual entities resent are computed based on the (dis)similarity
# between the nearest two members of each group
#
# problem: it soffers the "chain effect"
#
# how to use it: c is a lower triangular matrix(Inf instead of zeros),
# n the maximum number of iterations, clusters the number of desired cluster 
# and dist_max the distance maximum allowed in the formation of a new cluster.
# So, there are three way to stop the clustering.
#
# output: a list with: a string vector of the clusters (notice that a function is needed
# to link them to the right stations), the final matrix of the distances and the 
# minimum distance in the last matrix
hier_single_clustering <- function(c, n, clusters, dist_max){
  
  b <- 1:dim(c)[2]
  
  for( k in 1:n){
    
    d <- c
    minn <- which(c == min(c[,1:dim(c)[2]]), arr.ind = TRUE)
    minn_1 <- min(minn)
    minn_2 <- max(minn)
    l_1 <- dim(c)[1]
    l_2 <- dim(c)[2]
    l1 <- l_1 - 2
    b_clus <- b[minn]
    
    if(length(b) == 1 || length(b) == clusters || min(c[,1:dim(c)[2]]) >= dist_max){
      break
    }
    
    if( !(minn_1 == 1 || minn_2 == l_1 || l_1 == 2 || (minn_2 - minn_1) == 1) ){
      
      c <- rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(minn_2 - 1),],d[(minn_2+1):l_1,])
      c <- cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(minn_2 - 1)],c[,(minn_2+1):l_2])
      b <- c(b[1:(minn_1-1)],b[(minn_1+1):(minn_2-1)],b[(minn_2+1):l_1])
      aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_1]
      bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_2]
      
    }else{
      
      if(l_1 == 2){
        
        c <- c[2,1]
        b <- paste(min(b_clus),max(b_clus),sep=",")
        break
        
      }else if( (minn_2 - minn_1) == 1){
        
        if(minn_1 == 1){
          
          c <- matrix(d[(minn_2+1):l_1,],l_1 -2,l_1)
          c <- matrix(c[,(minn_2+1):l_2],l_1 -2,l_1 -2)
          b <- b[(minn_2+1):l_1]
          aaa1 <- d[minn_1,c((minn_2+1):l_1)]
          aaa2 <- d[minn_2,c((minn_2+1):l_1)]
          bbb1 <- d[c((minn_2+1):l_1),minn_1]
          bbb2 <- d[c((minn_2+1):l_1),minn_2]
          
        }else if(minn_2 == l_1){
          
          c <- matrix(d[1:(minn_1 - 1),],l_1 -2,l_1)
          c <- matrix(c[,1:(minn_1 - 1)],l_1 -2,l_1 -2)
          b <- b[1:(minn_1-1)]
          aaa1 <- d[minn_1,c(1:(minn_1-1))]
          aaa2 <- d[minn_2,c(1:(minn_1-1))]
          bbb1 <- d[c(1:(minn_1-1)),minn_1]
          bbb2 <- d[c(1:(minn_1-1)),minn_2]
          
        }else{
          
          c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
          c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2 )
          b <- c(b[1:(minn_1-1)],b[(minn_2+1):l_1])
          aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_2+1):l_1)]
          aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_2+1):l_1)]
          bbb1 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_1]
          bbb2 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_2]
        }
        
      }else if(minn_1 == 1){ 
        
        c <- matrix(rbind(d[2:(minn_2 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
        c <- matrix(cbind(c[,2:(minn_2 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2)
        b <- c(b[2:(minn_2-1)],b[(minn_2+1):l_1])
        aaa1 <- d[minn_1,c(2:(minn_2-1),(minn_2+1):l_1)]
        aaa2 <- d[minn_2,c(2:(minn_2-1),(minn_2+1):l_1)]
        bbb1 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_1]
        bbb2 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_2]
        
      }else{
        
        c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(l_1-1),]),l_1 -2,l_1)
        c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(l_2-1)]),l_1 -2,l_1 -2)
        b <- c(b[1:(minn_1-1)],b[(minn_1+1):(l_1-1)])
        aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_1]
        bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_2]
      }
    }
    
    c <- rbind(c,t(rep(Inf,l1)))
    c[l1 + 1,1:l1] <- pmin(pmin(aaa1,aaa2),pmin(bbb1,bbb2))
    c <- cbind(c,rep(Inf,l1+1))
    b <- c(b,paste(min(b_clus),max(b_clus),sep=","))
  }
  return(list(b,c,min(c[,1:dim(c)[2]])))
}

aaaaa1 <- hier_single_clustering(a,140,10,1000)



pmax2 <- function(a,b){   # pmax doesn't work well with our Inf
  a[is.infinite(a)] <- NA
  b[is.infinite(b)] <- NA
  e <- pmax(a,b, na.rm=TRUE)
  e[is.na(e)] <- Inf
  e
}

# Agglomerative clustering
# complete-likeage -> When a new cluster is formed, the (dis)similarities 
# between it and the other clusters and/or individual entities present
# are computed based on the (dis)similarity between the farthest two members 
# of each group

hier_complete_clustering <- function(c, n, clusters, dist_max){
  
  b <- 1:dim(c)[2]
  
  for( k in 1:n){
    
    d <- c
    minn <- which(c == min(c[,1:dim(c)[2]]), arr.ind = TRUE)
    minn_1 <- min(minn)
    minn_2 <- max(minn)
    l_1 <- dim(c)[1]
    l_2 <- dim(c)[2] 
    l1 <- l_1 - 2
    b_clus <- b[minn]
    
    if(length(b) == 1 || length(b) == clusters || min(c[,1:dim(c)[2]]) >= dist_max){
      break
    }
    
    if( !(minn_1 == 1 || minn_2 == l_1 || l_1 == 2 || (minn_2 - minn_1) == 1) ){
      
      c <- rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(minn_2 - 1),],d[(minn_2+1):l_1,])
      c <- cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(minn_2 - 1)],c[,(minn_2+1):l_2])
      b <- c(b[1:(minn_1-1)],b[(minn_1+1):(minn_2-1)],b[(minn_2+1):l_1])
      aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_1]
      bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_2]
      
    }else{
      
      if(l_1 == 2){
        
        c <- c[2,1]
        b <- paste(min(b_clus),max(b_clus),sep=",")
        break
        
      }else if( (minn_2 - minn_1) == 1){
        
        if(minn_1 == 1){
          
          c <- matrix(d[(minn_2+1):l_1,],l_1 -2,l_1)
          c <- matrix(c[,(minn_2+1):l_2],l_1 -2,l_1 -2)
          b <- b[(minn_2+1):l_1]
          aaa1 <- d[minn_1,c((minn_2+1):l_1)]
          aaa2 <- d[minn_2,c((minn_2+1):l_1)]
          bbb1 <- d[c((minn_2+1):l_1),minn_1]
          bbb2 <- d[c((minn_2+1):l_1),minn_2]
          
        }else if(minn_2 == l_1){
          
          c <- matrix(d[1:(minn_1 - 1),],l_1 -2,l_1)
          c <- matrix(c[,1:(minn_1 - 1)],l_1 -2,l_1 -2)
          b <- b[1:(minn_1-1)]
          aaa1 <- d[minn_1,c(1:(minn_1-1))]
          aaa2 <- d[minn_2,c(1:(minn_1-1))]
          bbb1 <- d[c(1:(minn_1-1)),minn_1]
          bbb2 <- d[c(1:(minn_1-1)),minn_2]
          
        }else{
          
          c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
          c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2 )
          b <- c(b[1:(minn_1-1)],b[(minn_2+1):l_1])
          aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_2+1):l_1)]
          aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_2+1):l_1)]
          bbb1 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_1]
          bbb2 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_2]
        }
        
      }else if(minn_1 == 1){ 
        
        c <- matrix(rbind(d[2:(minn_2 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
        c <- matrix(cbind(c[,2:(minn_2 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2)
        b <- c(b[2:(minn_2-1)],b[(minn_2+1):l_1])
        aaa1 <- d[minn_1,c(2:(minn_2-1),(minn_2+1):l_1)]
        aaa2 <- d[minn_2,c(2:(minn_2-1),(minn_2+1):l_1)]
        bbb1 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_1]
        bbb2 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_2]
        
      }else{
        
        c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(l_1-1),]),l_1 -2,l_1)
        c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(l_2-1)]),l_1 -2,l_1 -2)
        b <- c(b[1:(minn_1-1)],b[(minn_1+1):(l_1-1)])
        aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_1]
        bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_2]
      }
    }
    
    c <- rbind(c,t(rep(Inf,l1)))
    c[l1 + 1,1:l1] <- pmax2(pmax2(aaa1,aaa2),pmax2(bbb1,bbb2))
    c <- cbind(c,rep(Inf,l1+1))
    b <- c(b,paste(min(b_clus),max(b_clus),sep=","))
    
  }
  return(list(b,c,min(c[,1:dim(c)[2]])))
}

aaaaa2 <- hier_complete_clustering(a,140,10,100)





pmean <- function(a,b){
  a[is.infinite(a)] <- NA
  b[is.infinite(b)] <- NA
  e <- colMeans(rbind(a,b), na.rm=TRUE)
  e[is.na(e)] <- Inf
  e
}

# Agglomerative clustering
# average-likeage -> When a new cluster is formed, the (dis)similarities 
# between it and the other clusters and/or individual entities present are
# computed based on the average (dis)similarity between all members in each group 

hier_average_clustering <- function(c, n, clusters, dist_max){
  
  b <- 1:dim(c)[2]
  
  for( k in 1:n){
    
    d <- c
    minn <- which(c == min(c[,1:dim(c)[2]]), arr.ind = TRUE)
    minn_1 <- min(minn)
    minn_2 <- max(minn)
    l_1 <- dim(c)[1]
    l_2 <- dim(c)[2] 
    l1 <- l_1 - 2
    b_clus <- b[minn]
    
    if(length(b) == 1 || length(b) == clusters || min(c[,1:dim(c)[2]]) >= dist_max){
      break
    }
    
    if( !(minn_1 == 1 || minn_2 == l_1 || l_1 == 2 || (minn_2 - minn_1) == 1) ){
      
      c <- rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(minn_2 - 1),],d[(minn_2+1):l_1,])
      c <- cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(minn_2 - 1)],c[,(minn_2+1):l_2])
      b <- c(b[1:(minn_1-1)],b[(minn_1+1):(minn_2-1)],b[(minn_2+1):l_1])
      aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1)]
      bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_1]
      bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(minn_2-1),(minn_2+1):l_1),minn_2]
      
    }else{
      
      if(l_1 == 2){
        
        c <- c[2,1]
        b <- paste(min(b_clus),max(b_clus),sep=",")
        break
        
      }else if( (minn_2 - minn_1) == 1){
        
        if(minn_1 == 1){
          
          c <- matrix(d[(minn_2+1):l_1,],l_1 -2,l_1)
          c <- matrix(c[,(minn_2+1):l_2],l_1 -2,l_1 -2)
          b <- b[(minn_2+1):l_1]
          aaa1 <- d[minn_1,c((minn_2+1):l_1)]
          aaa2 <- d[minn_2,c((minn_2+1):l_1)]
          bbb1 <- d[c((minn_2+1):l_1),minn_1]
          bbb2 <- d[c((minn_2+1):l_1),minn_2]
          
        }else if(minn_2 == l_1){
          
          c <- matrix(d[1:(minn_1 - 1),],l_1 -2,l_1)
          c <- matrix(c[,1:(minn_1 - 1)],l_1 -2,l_1 -2)
          b <- b[1:(minn_1-1)]
          aaa1 <- d[minn_1,c(1:(minn_1-1))]
          aaa2 <- d[minn_2,c(1:(minn_1-1))]
          bbb1 <- d[c(1:(minn_1-1)),minn_1]
          bbb2 <- d[c(1:(minn_1-1)),minn_2]
          
        }else{
          
          c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
          c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2 )
          b <- c(b[1:(minn_1-1)],b[(minn_2+1):l_1])
          aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_2+1):l_1)]
          aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_2+1):l_1)]
          bbb1 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_1]
          bbb2 <- d[c(1:(minn_1-1),(minn_2+1):l_1),minn_2]
        }
        
      }else if(minn_1 == 1){ 
        
        c <- matrix(rbind(d[2:(minn_2 - 1),],d[(minn_2+1):l_1,]),l_1 -2,l_1)
        c <- matrix(cbind(c[,2:(minn_2 - 1)],c[,(minn_2+1):l_2]),l_1 -2,l_1 -2)
        b <- c(b[2:(minn_2-1)],b[(minn_2+1):l_1])
        aaa1 <- d[minn_1,c(2:(minn_2-1),(minn_2+1):l_1)]
        aaa2 <- d[minn_2,c(2:(minn_2-1),(minn_2+1):l_1)]
        bbb1 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_1]
        bbb2 <- d[c(2:(minn_2-1),(minn_2+1):l_1),minn_2]
        
      }else{
        
        c <- matrix(rbind(d[1:(minn_1 - 1),],d[(minn_1+1):(l_1-1),]),l_1 -2,l_1)
        c <- matrix(cbind(c[,1:(minn_1 - 1)],c[,(minn_1+1):(l_2-1)]),l_1 -2,l_1 -2)
        b <- c(b[1:(minn_1-1)],b[(minn_1+1):(l_1-1)])
        aaa1 <- d[minn_1,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        aaa2 <- d[minn_2,c(1:(minn_1-1),(minn_1+1):(l_2-1))]
        bbb1 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_1]
        bbb2 <- d[c(1:(minn_1-1),(minn_1+1):(l_2-1)),minn_2]
      }
    }
    
    c <- rbind(c,t(rep(Inf,l1)))
    c[l1 + 1,1:l1] <- pmean(pmean(aaa1,aaa2),pmean(bbb1,bbb2))
    c <- cbind(c,rep(Inf,l1+1))
    b <- c(b,paste(min(b_clus),max(b_clus),sep=","))
  }
  #out <<- c
  return(list(b,c,min(c[,1:dim(c)[2]])))
}

aaaaa3 <- hier_average_clustering(a,180,10,100)

aaaaa4 <- hier_average_clustering(a[1:30,1:30],180,10,100)





################################################################################

#library(stringi)

# These are all functions to check that there are no repetations in the clustering 
# and they extract the ID of the stations 
char_num <- function(x){
  return(as.numeric(unlist(x)))
}

char_num_l <- function(x,sep=","){
  stat <- strsplit(x,sep)
  return(lapply(stat, char_num))
}

check_no_rep <- function(x,av){
  a <- sort(unlist(x)) # unroll 
  numb_av <- dim(av)[1]
  if(numb_av != length(a) || any(sort(unique(a)) != (1:numb_av)) ){
    # print("something is wrong!!")
    return(0)
  }else{
    # print("ok,no rep and same length")
    return(1)
  }
}

retrieve_id <- function(s,y){ # clusters, stations' table
  x <- char_num_l(s)
  if(check_no_rep(x,y) != 0){
    stat2 <- lapply(x,function(x,y){y[x,]$ID},y)  # DON'T CHANGE THE field ID
    return(stat2)
  }else{
    return(0)
  } 
}

aaaaa3_clus_id <- retrieve_id(aaaaa3[[1]],av)
