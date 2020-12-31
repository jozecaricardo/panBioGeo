####### Hovenkamp procedure to detect allopatry ########
########################################################
is_allopatry <- function(results, phylogeny, coordin){
  
  ##### packages ######
  library(stringr)
  library(stringi)
  library(ctv)
  library(phangorn)
  library(raster)
  library(viridis)
  library(dismo)
  library(phytools)
  library(geiger)
  
  #### raster files in .tif and shape files in .shp ####
  # mapas.p <- list.files(path = 'out', pattern = 'tif', full.names = T)
  # presenceRaster <- stack(mapas.p)
  
  # presenceShape <- rgdal::readOGR(dsn = 'out/pointshape_E_L_cribrarius.shp', verbose = FALSE)
  
  
  ############# preparing the phylogeny ###################
  
  # print('1) preparing the coordinates...')
  coordin <- matrix(as.matrix(coordin[, c(2, 3)]), nrow(coordin),
           2, dimnames = list(coordin[,1], colnames(coordin)[c(2, 3)]))
  
  # colores <- setNames(sample(viridis(n = length(unique(rownames(coordin))))),
  #                unique(rownames(coordin)))
  
  if(!is.null(phylogeny)){
    print('1) preparing the phylogeny...')
    print('...checking names from dataset and names in the species phylogeny.')
    if(length(name.check(phylogeny, coordin)) != 1){
      chk <- name.check(phylogeny, coordin)
      phylogeny <- drop.tip(phylogeny, chk$tree_not_data)
    }
    
    # phylogenetic trees and nodes:
    if(is.null(phylogeny$edge.length) == TRUE){
      print('...creating the branch lengths of a phylogeny equal to one...')
      phylogeny <- compute.brlen(phylogeny, 1)
    }
  }
  
  #############################################
  ######### looping the nodes #################
  
  arv <- reorder(phylogeny, "postorder") #reordenando os levels
  
  e1 <- arv$edge[, 1]
  e2 <- arv$edge[, 2]
  # EL <- arv$edge.length
  # e1 #nós internos
  # e2 #nós terminais e raiz
  # EL #ramos
  
  
  nb.tip <- length(arv$tip.label)
  nb.node <- arv$Nnode
  
  # x <- reorder(x = as.factor(ranges), X = c(1:length(ranges)))
  
  #### grid numbers in which that species occurs ######
  ### grid numnber list of all species:
  gridList <- list()
  
  results <- results[-dim(results)[1],] # without root
  
  for(j in colnames(results)){
    gridList[[j]] <- names(which(results[,j] == 1))  # grid numbers of all species
  }
  
  nll <- length(gridList)
  
  # ranges <- LETTERS[1:nll]
  
  # x <- reorder(x = as.factor(ranges), X = c(1:length(ranges)))
  
  # suma1 <- matrix(NA, nb.tip + nb.node, nll)
  
  # cor <- cores[as.integer(x)]
  # cols <- setNames(cor, x)
  # ARE <- to.matrix(ranges, x)
  
  # ARE <- ARE[-c(dim(ARE)[1]:length(arv$tip.label) + 1),]
  
  # lvls <- levels(x)
  # x <- as.integer(x)
  
  # suma2 <- matrix(NA, nb.tip + nb.node, nll)
  
  # TIPS <- 1:nb.tip
  
  # suma2[cbind(TIPS, x[1:nb.tip])] <- lvls
  
  
  ## open a file
  logfile <- "out/vicarius.txt"
  
  cat(c("species", "pattern", "\n"), file = logfile, sep="\t")
  
  ##### help variables ####
  tabela <- mrca(phy = arv, full = FALSE)
  tabela.v <- as.vector(tabela)
  names(tabela.v) <- rep(colnames(tabela), dim(tabela)[1])
  
  for(i in seq(from = 1, by = 2, length.out = nb.node)){
    j <- i + 1L
    ance <- e1[i]
    des1 <- e2[i]
    des2 <- e2[j]
    
    lis <- arv$tip.label[arv$edge[which(arv$edge[,1] == ance), 2]]
    
    if(any(is.na(lis)) == FALSE){
      comparison <- intersect(gridList[[lis[1]]], gridList[[lis[2]]])
      
      # vicariance or sympatry?
      if(!is.null(comparison)){
        cat(c(paste0(c(lis), collapse = ', '), 'SYMPATRY', '\n'), sep = '\t',
            file = logfile, append = TRUE)
      } else {
        cat(c(paste0(c(lis), collapse = ', '), 'ALLOPATRY', '\n'), sep = '\t',
            file = logfile, append = TRUE)
      }
    } else if(any(is.na(lis)) == TRUE){
      # internal and terminal nodes
      node_des1 <- which(tabela %in% des1)		
      cladis_des1 <- unique(names(tabela.v)[node_des1])
      
      node_des2 <- which(tabela %in% des2)		
      cladis_des2 <- unique(names(tabela.v)[node_des2])
      
      
      comparison <- intersect(unlist(gridList[cladis_des1]), unlist(gridList[cladis_des2]))
      
      # vicariance or sympatry?
      if(!is.null(comparison)){
        cat(c(paste0(c(cladis_des1, cladis_des2), collapse = ', '), 'SYMPATRY', '\n'), sep = '\t',
            file = logfile, append = TRUE)
      } else {
        cat(c(paste0(c(cladis_des1, cladis_des2), collapse = ', '), 'ALLOPATRY', '\n'), sep = '\t',
            file = logfile, append = TRUE)
      }
      
    }
  }
  print('Done.')
  ##########################################################
  ##########################################################
}




