setwd("/home/user/Desktop/method")
setwd("./data")

#example data
dataFilePaths <- list("workloada.csv", "workloadb.csv", "workloadc.csv")
firstLines    <- list( 37,              45,              36)
lastLines     <- list( 306,             314,             305)

getPreprocedData <- function(dataFilePaths, firstLines, lastLines){
  data <- list()
  for(i in 1:length(dataFilePaths)){
    data[[i]] <- read.table(dataFilePaths[[i]], header=TRUE, sep=",", dec=".", fileEncoding="windows-1252") 
    data[[i]] <- data[[i]][firstLines[[i]]:lastLines[[i]],]
    data[[i]][,1] <- NULL
    naColNames <- colnames(data[[i]])[colSums(is.na(data[[i]])) > 0]
    data[[i]][ ,naColNames] <- NULL
    variances <- data.frame(data[[i]][1,], row.names=c("variance"))
    variances[,] <- 0
    for(j in 1:length(data[[i]])){
      variances[j] <- var(x=data[[i]][,j], na.rm=TRUE) 
    }
    names <- colnames(variances)[variances <= 0]
    data[[i]][,names] <- NULL
  }
  return (data)
}

getGroupsFromCluster <- function(data){
  p_mtx <- matrix(nrow=ncol(data), ncol=ncol(data))
  rownames(p_mtx) <- colnames(data)
  colnames(p_mtx) <- colnames(data)
  for(i in 1:ncol(data)){
    for(j in i:ncol(data)){
      p_var <- cor(x = data[,i], y = data[,j], method = "pearson")
      if(p_var >= 0){
        p_diss <- 1-p_var
      }else{
        p_diss <- p_var+1
      }
      p_mtx[i,j] <- p_diss
      p_mtx[j,i] <- p_diss
    }
  }
  diss_mtx <- as.dist(p_mtx)
  cluster_h <- hclust(diss_mtx)
  cluster <- cutree(cluster_h, h=0.1)
  groups <- list()
  for(i in 1:max(cluster)){
    groups[[length(groups)+1]] <- list()
  }
  for(i in 1:length(cluster)){
    groups[[cluster[i]]][length(groups[[cluster[i]]])+1] <- names(cluster[i])
  }
  return(list(groups,p_mtx))
}

getSlicedGroups <- function(exp1Groups, exp2Groups){
  groups <- list()
  for(i in 1:length(exp1Groups)){
    for(j in 1:length(exp2Groups)){
      if(length(exp1Groups[[i]]) != 0 && length(exp2Groups[[j]]) != 0){
        expGroupIntersect <- intersect(exp1Groups[[i]], exp2Groups[[j]])
        if(length(expGroupIntersect) != 0){
          exp1Groups[[i]] <- setdiff(exp1Groups[[i]], expGroupIntersect)
          exp2Groups[[j]] <- setdiff(exp2Groups[[j]], expGroupIntersect)
          groups[[length(groups)+1]] <- expGroupIntersect
        }
      }
    }
  }
  for(i in 1:length(exp1Groups)){
    if(length(exp1Groups[[i]]) != 0 ){
      groups[[length(groups)+1]]<-exp1Groups[[i]]
    }
  }
  for(j in 1:length(exp2Groups)){
    if(length(exp2Groups[[j]]) != 0 ){
      groups[[length(groups)+1]] <- exp2Groups[[j]]
    }
  }
  return (groups)
}

getRepresentativeCounters <- function(groups, p_mtx){
  rCounters <- list()
  for(i in 1:length(groups)){
    lastSum <- Inf
    lastCounter <- ""
    for(j in 1:length(groups[[i]])){
      counter <- groups[[i]][[j]]
      sum <- 0
      for(k in 1:length(groups[[i]])){
        if(j!=k){
          for(z in 1:length(p_mtx)){
            if(!is.na(match(groups[[i]][[j]],rownames(p_mtx[[z]]))) & !is.na(match(groups[[i]][[k]],rownames(p_mtx[[z]])))){
              sum <- sum + p_mtx[[z]][groups[[i]][[j]],groups[[i]][[k]]]
            }
            else{
              sum <- sum + 1
            }
          }
        }
      }
      if(sum < lastSum){
        lastSum <- sum
        lastCounter <- counter
      }
    }
    rCounters[[i]] <- lastCounter
  }
  return (rCounters)
}

data <- getPreprocedData(dataFilePaths, firstLines, lastLines)

groups <- list()
p_mtx <- list()
for(i in 1:length(data)){
  result <- getGroupsFromCluster(data[[i]])
  groups[[i]]<- result[[1]]
  p_mtx[[i]] <- result[[2]]
}

newGroups <- list()
newGroups[[1]] <- getSlicedGroups(groups[[1]], groups[[2]])
#length(newGroups[[1]])
for(i in 3:length(groups)){
  j <- (i-2)
  newGroups[[j+1]] <- getSlicedGroups(newGroups[[j]], groups[[i]])    
  #print(length(newGroups[[j+1]]))
}

counters <- getRepresentativeCounters(newGroups[[length(newGroups)]], p_mtx)
counters

