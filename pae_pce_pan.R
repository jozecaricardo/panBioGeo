#### PAE-PCE analysis ####
# matric = the presence-absence matrix from terminal_node function
# rou = number of iterations in the Ratchet parsimony
pae_pce <- function(matric, k = 25, minit = 5, maxit = 10000, oneGrid = FALSE, shapeFile, resolut,
                    gridView = FALSE, labelGrid = FALSE, legendSpecies = TRUE, sobrepo = FALSE,
                    nVoltas = NULL, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, N = 1){
	# install.packages('ctv')
  # install.packages("stringr")
  # install.packages('stringi')
  # install.packages('TreeSearch')
  library(stringr)
  library(stringi)
  library(ctv)
  library(phangorn)
  library(raster)
  library(viridis)
  library(TreeSearch) # TBR search
	# install.views('Phylogenetics')

	#Convert the character matrix to a phyDat object 
	# in order to infer a tree. By assigning the matrix 
	# type as 'user' we can specify the components of the 
	# matrix (in this case, binary 1s and 0s).

	tempMatrix <- as.phyDat(matric, type = 'USER', levels = c(0, 1))

	colu <- ncol(matric)
	if(oneGrid == FALSE){
	  # creating a restriction:
	  riVec <- rep(1.0, colu)
	} else {
	  # creating a restriction:
	  ciVec <- rep(1.0, colu)
	  riVec <- rep(0, colu)
	}

	contagem <- 0
	
	lista <- list()
	
	print('Please, remember that the resulting rasters will be kept in your setted directory!')
	
	if(sobrepo == TRUE){
	  print('Please do not forgive for indicating the number of iterations in the nVoltas argument!')
	  if(is.null(nVoltas)){
	    stop('You should to indicate the number of iterations in the nVoltas argument!')
	  }
	}
	
	######################
	##### shape file #####
	grid <- raster(extent(shapeFile), resolution = resolut, crs = CRS("+proj=longlat +datum=WGS84"))
	grid <- raster::extend(grid, c(1, 1))
	gridPolygon <- rasterToPolygons(grid)
	suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
	#proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
	# projection(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")
	
	
	# clipping the intersected cells:
	suppressWarnings(cropped_map <- raster::intersect(gridPolygon, shapeFile))
	
	# producing a raster of the shapefile
	mask.raster <- raster(extent(shapeFile), resolution = resolut,
	                      crs = CRS("+proj=longlat +datum=WGS84"))
	suppressWarnings(r <- rasterize(shapeFile, mask.raster))
	proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
	# mask.raster[is.na(mask.raster)] <- 0
	r <- merge(r, mask.raster)
	
	colores <- viridis(n = length(r[!is.na(r)]))
	
	if(sobrepo == FALSE){
	  plot(shapeFile, axes = TRUE, main = 'Generalized tracks converted to grids',
	       sub = paste0(c('resolution: ', resolut[1], ' x ', resolut[2]), collapse = ''), cex.main = 0.9,
	       cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
	  abline(h = 0, col = 'red', lty = 2) # equator
	  
	  iterat <- 0
	} else if(sobrepo == TRUE){
	  iterat <- 0
	}
	
	# set the cells associated with the shapfile to the specified value
	r[r == 1] <- NA
	
	
	####### Detecting only those grids being supported by autaphomorphies ########
	if(oneGrid == TRUE){
	  
	  #### PAE-PCE ###
	  while(length(which((riVec == 0) == TRUE)) > 1 && length(which(ciVec > 0) == TRUE) > 1 || length(ciVec) > 1){
	    contagem <- contagem + 1
	    
	    # We can get a random starting tree for a parsimony search using the
	    #  random.addition function:
	    # ra.tre <- njs(dist.hamming(tempMatrix))
	    
	    if(dim(as.data.frame(tempMatrix))[1] == 0){
	      print('This iteration has to be stopped because there is no more autapomorphies!')
	      break
	    }
	    
	    ### TBR + RATCHET searches
	    dm <- dist.hamming(tempMatrix)
	    ra.tre <- NJ(dm)
	    
	    # Use the pratchet function (parsimony ratchet) to find the most
	    # parsimonious tree. Specifying k means the algorithm will search   
	    # through k possible trees to find the most parsimonious solution.
	    treeIt1 <- list()
	    
	    for(rou in 1:N){
	      treeIt1[[rou]] <- pratchet(tempMatrix, k = k, trace = 0, minit = minit, maxit = maxit,
	                                 start = ra.tre, all = TRUE)
	      
	      if(length(treeIt1[[rou]]) > 1){
	        treeIt1[[rou]] <- consensus(treeIt1[[rou]], p = 0.50)
	      }
	      
	    }
	    
	    if(is.null(treeIt1)){
	      stop('This iteration has to be stopped because there is no more autapomorphies!')
	    }
	    
	    treeT <- consensus(treeIt1, p = 0.50) # 'majority rule': change p to 0.5
	    # plot(treeT)
	    
	    # treeComp <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
	    #                               rearrangements = 'SPR')
	    # treeComp <- multi2di(treeComp)
	    
	    # scoreTreeInit <- parsimony(treeComp, tempMatrix, method = 'fitch')
	    
	    
	    # treeIt <- list()
	    # treeIt <- TBR(treeComp)
	    # treeIt <- consensus(treeIt1, p = 1.0) # 'majority rule': change p to 0.5
	    
	    # treeIt2[[times]] <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
	    #                            rearrangements = 'SPR')
	    # parsimony(c(treeIt1, treeIt2), tempMatrix, method = 'fitch')
	    
	    # Root the tree by the designated outgroup 
	    # (write the species name as it appears in the tree, 
	    # but add an underscore to fill the space).
	    strictTree <- root(treeT, outgroup = "ROOT", resolve.root = TRUE)
	    # assign edge length
	    # treeRatchet <- acctran(chartree, skullchars)
	    
	    # CI of each character (i.e., taxon):
	    ciVec <- CI(strictTree, data = tempMatrix, sitewise = TRUE)
	    ciVec[which(is.na(ciVec))] <- 0
	    
	    # RI of each character (i.e., taxon):
	    riVec <- RI(strictTree, data = tempMatrix, sitewise = TRUE)
	    riVec[which(is.na(riVec))] <- 0
	    
	    # if(any(ciVec == 1.0) == TRUE){
	    #	colu <- colu - 1
	    #}
	    
	    ## first load packages
	    require(phytools)
	    tree_temp <- strictTree
	    if(is.null(tree_temp$edge.length) == TRUE){
	      print('...creating the branch lengths of a tree equal to one...')
	      tree_temp <- compute.brlen(tree_temp, 1)
	    }
	    
	    #Discrete mapping
	    inter1 <- which(ciVec > 0); inter2 <- which(riVec == 0)
	    interT <- intersect(inter1, inter2)
	    matTemp <- matric[, interT] # temporary matrix with grids supported by autaphomorphies
	    
	    lista_T <- list()
	    
	    nomesCOL <- colnames(matTemp) # species which are synapomorphies
	    if(is.null(nomesCOL) && contagem != 1){
	      print('This analysis had to be stopped because there is no more synapomorphies!')
	      
	      break
	      
	    } else if(is.null(nomesCOL) && contagem == 1){
	      stop('This iteration had to be stopped because there is no more synapomorphies!')
	    }
	    
	    for(j in nomesCOL){
	      discre <- matTemp[,j]
	      for(a in 1:length(discre)){
	        ifelse(discre[a] == 1, discre[a] <- j, discre[a] <- 'absent')
	      } # vector with the presence and absence of one species
	      #This maps the character (i.e., taxon).
	      # x11()
	      # plotBranchbyTrait(chartree, tempMatrix[,which(ciVec == 1)])
	      ## simulate single stochastic character map using empirical Bayes method
	      discre <- as.factor(discre)  ### vetor de presença da espécie
	      
	      lista_T[[j]] <- names(discre)[which(discre == j)]
	    }
	    
	    n_occur <- data.frame(table(unlist(lista_T))) #frequency of each grid
	    # n_occur
	    
	    
	    syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
	    # syn_grids
	    if(length(syn_grids) < 2){
	      print('This iteration had to be stopped because there is no more synapomorphies!')
	      break
	    }
	    
	    speciesNames <- str_replace(names(syn_grids), "[0123456789]", "")
	    
	    speciesNames <- str_replace(speciesNames, "[0123456789]", "")
	    
	    # data frame with the results...
	    frameTemp <- data.frame(spp = speciesNames, grid_n = syn_grids, row.names = 1:length(syn_grids))
	    
	    # matrix with species names and grid numbers:
	    resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
	                             1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))
	    
	    if(length(unique(rownames(resulPaeRaster))) == 0){
	      stop('This iteration had to be stopped because there is no more synapomorphies!')
	    }
	    
	    posic <- 0
	    for(pos in 1:length(unique(rownames(resulPaeRaster)))){
	      posic[pos] <- which(colnames(matric) == unique(rownames(resulPaeRaster))[pos])
	    }
	    
	    # changing the matrix for the next iteration:	
	    tempMatrix <- subset(tempMatrix, select=-posic, site.pattern = FALSE)
	    
	    ###########################
	    
	    ## the number of species supporting the grid number
	    speciesNumber <- data.frame(table(resulPaeRaster))
	    # speciesNumber
	    
	    similarNames <- list()
	    for(volta in 1:dim(speciesNumber)[1]){
	      similarNames[[volta]] <- resulPaeRaster[which(resulPaeRaster[,1] == speciesNumber$resulPaeRaster[volta]), ]
	    }
	    # similarNames
	    
	    ## grid numbers that are generalized tracks in this iteration
	    gridIt <- unique(as.numeric(unlist(similarNames)))
	    # gridIt
	    
	    
	    # set the cells associated with the shapfile to the specified value
	    r[r > 0] <- NA
	    
	    for(j in 1:dim(speciesNumber)[1]){
	      gTrack <- as.numeric(as.character(speciesNumber[j,1]))
	      values(r)[gTrack] <- 1
	    }
	    
	    if(sobrepo == TRUE){
	      plot(r, axes = FALSE, legend = FALSE, add = TRUE, col = colores[(contagem + nVoltas)],
	         alpha = 0.60)
	    } else if(sobrepo == FALSE){
	      plot(r, axes = FALSE, legend = FALSE, add = TRUE, col = colores[contagem],
	           alpha = 0.50)
	    }
	    
	    # producing a raster:
	    #convert the raster to points for plotting the number of a grid
	    map.r <- as.data.frame(rasterToPoints(r))
	    pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1) # raster com as presenças
	    suppressWarnings(writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
	                                 format = "GTiff", overwrite = TRUE))
	    
	    if(gridView == TRUE){
	      plot(cropped_map, add = TRUE)
	      
	      map.r$gridNumber <- speciesNumber$resulPaeRaster
	      
	      if(labelGrid == TRUE){
	        text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = (2.5 * resolut[1]) / 10, col = 'red')
	      }
	    }
	    
	    if(legendSpecies == TRUE){
	      print('Please, check the legendSpeciesB.txt file in the /out directory...')
	      #local <- locator()
	      logfile <- "out/legendSpeciesB.txt"
	      cat(c("grid_number", "species", "\n"), file = logfile, sep="\t")
	      for(i in as.numeric(as.character(speciesNumber$resulPaeRaster))){
	        cat(c(paste0(c('grid_', i), collapse = ''), rownames(subset(resulPaeRaster,
	                       resulPaeRaster[,'grid_n'] == i)), '\n'), sep = '\t',
	            file = logfile, append = TRUE)
	      }
	    }
	    
	    # legend
	    if(sobrepo == FALSE){
	      iterat[contagem] <- contagem
	    } else if(sobrepo == TRUE){
	      iterat[contagem] <- contagem + nVoltas
	    }
	    
	    
	    
	    lista[[contagem]] <- frameTemp
	    
	    print(paste0(c('The iteration number ', contagem, ' has finished.'), collapse = ''))
	    
	  } # closing 'while' looping
	  
	  if(is.null(nomesCOL) || dim(as.data.frame(tempMatrix))[1] == 0){
	    iterat <- iterat[contagem - 1]
	  } else if(!is.null(nomesCOL) || dim(as.data.frame(tempMatrix))[1] != 0){
	    iterat <- iterat
	  }
	  
	  if(sobrepo == TRUE){
	    if(is.null(iterat)){
	      break
	    } else if(!is.null(iterat)){
	      legend(x = 'bottomright', legend = iterat, pch = 15, col = colores[(nVoltas + 1):(length(iterat) + nVoltas)],
	           title = 'Adding generalized track of the iterations...',
	           title.col = 'blue', pt.cex = 1.5, cex = 0.8)
	    }
	  } else if (sobrepo == FALSE){
	    if(is.null(iterat)){
	      break
	    } else if(!is.null(iterat)){
	      legend(x = 'topright', legend = iterat, pch = 15, col = colores[1:length(iterat)],
	           title = 'Adding generalized track of the iterations...',
	           title.col = 'red', pt.cex = 1.5, cex = 0.8)
	    }
	  }
	  
	  if(contagem == 1){
	    print('There is no result!')
	  } else if(contagem > 1){
	    return(lista)
	  }
	  # closing if... else... condition
	} else {
	  
	  ####### Detecting those grids being supported by synapomorphies ########
	  
	  #### PAE-PCE ###
	  
	  while(length(which((riVec > 0) == TRUE)) > 1 || length(riVec) > 1){
	    contagem <- contagem + 1
	    
	    # We can get a random starting tree for a parsimony search using the
	    #  random.addition function:
	    # ra.tre <- njs(dist.hamming(tempMatrix))
	    
	    if(dim(as.data.frame(tempMatrix))[1] == 0){
	      print('This iteration had to be stopped because there is no more synapomorphies!')
	      break
	    }
	    
	    ### TBR + RATCHET searches
	    dm <- dist.hamming(tempMatrix)
	    ra.tre <- NJ(dm)
	    
	    # Use the pratchet function (parsimony ratchet) to find the most
	    # parsimonious tree. Specifying k means the algorithm will search   
	    # through k possible trees to find the most parsimonious solution.
	    treeIt1 <- list()
	    
	    for(rou in 1:N){
	      treeIt1[[rou]] <- pratchet(tempMatrix, k = k, trace = 0, minit = minit, maxit = maxit,
	                         start = ra.tre, all = TRUE)
	      
	      if(length(treeIt1[[rou]]) > 1){
	        treeIt1[[rou]] <- consensus(treeIt1[[rou]], p = 0.50)
	      }
	    }
	    
	    if(is.null(treeIt1)){
	      stop('This iteration had to be stopped because there is no more autapomorphies!')
	    }
	    
	    
	    
	    treeT <- consensus(treeIt1, p = 0.50) # 'majority rule': change p to 0.5
	    
	    # treeComp <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
	    #                               rearrangements = 'SPR')
	    # treeComp <- multi2di(treeComp)
	    
	    # scoreTreeInit <- parsimony(treeComp, tempMatrix, method = 'fitch')
	    
	    
	    # treeIt <- list()
	    # treeIt <- TBR(treeComp)
	    # treeIt <- consensus(treeIt1, p = 1.0) # 'majority rule': change p to 0.5
	    
	    # treeIt2[[times]] <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
	    #                            rearrangements = 'SPR')
	    # parsimony(c(treeIt1, treeIt2), tempMatrix, method = 'fitch')
	    
	    # Root the tree by the designated outgroup 
	    # (write the species name as it appears in the tree, 
	    # but add an underscore to fill the space).
	    strictTree <- root(treeT, outgroup = "ROOT", resolve.root = TRUE)
	    # x11()
	    # plot(strictTree)
	    # assign edge length
	    # treeRatchet <- acctran(tree = treeIt2, data = tempMatrix)
	    
	    # ancestral states:
	    # ancResul <- ancestral.pars(tree = strictTree, data = tempMatrix)
	    # plotAnc(tree = strictTree, data = ancResul, attr(ancResul, "index")[2])
	    
	    # CI of each character (i.e., taxon):
	    #ciVec <- CI(strictTree, data = tempMatrix, sitewise = TRUE)
	    # ciVec[which(is.na(ciVec))] <- 0
	    
	    # RI of each character (i.e., taxon):
	    riVec <- RI(strictTree, data = tempMatrix, sitewise = TRUE)
	    riVec[which(is.na(riVec))] <- 0
	    
	    # if(any(ciVec == 1.0) == TRUE){
	    #	colu <- colu - 1
	    #}
	    
	    ## first load packages
	    require(phytools)
	    tree_temp <- strictTree
	    if(is.null(tree_temp$edge.length) == TRUE){
	      print('...creating the branch lengths of a tree equal to one...')
	      tree_temp <- compute.brlen(tree_temp, 1)
	    }
	    
	    #Discrete mapping
	    matTemp <- matric[,which(riVec > 0)] # temporary matrix with synapomorphies
	    
	    lista_T <- list()
	    
	    nomesCOL <- colnames(matTemp) # species which are synapomorphies
	    if(is.null(nomesCOL) && contagem != 1){
	      print('This analysis had to be stopped because there is no more synapomorphies!')
	      break
	    } else if(is.null(nomesCOL) && contagem == 1){
	      print('This iteration had to be stopped because there is no more synapomorphies!')
	      next
	    }
	    
	    for(j in nomesCOL){
	      discre <- matTemp[,j]
	      for(a in 1:length(discre)){
	        ifelse(discre[a] == 1, discre[a] <- j, discre[a] <- 'absent')
	      } # vector with the presence and absence of one species
	      #This maps the character (i.e., taxon).
	      # x11()
	      # plotBranchbyTrait(chartree, tempMatrix[,which(ciVec == 1)])
	      ## simulate single stochastic character map using empirical Bayes method
	      discre <- as.factor(discre)  ### vetor de presença da espécie
	      
	      lista_T[[j]] <- names(discre)[which(discre == j)]
	    }
	    
	    n_occur <- data.frame(table(unlist(lista_T))) #frequency of each grid
	    # n_occur
	    
	    
	    syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
	    # syn_grids
	    if(length(syn_grids) < 2){
	      print('This iteration has to be stopped because there is no more synapomorphies!')
	      next
	    }
	    
	    speciesNames <- str_replace(names(syn_grids), "[0123456789]", "")
	    
	    speciesNames <- str_replace(speciesNames, "[0123456789]", "")
	    
	    # data frame with the results...
	    frameTemp <- data.frame(spp = speciesNames, grid_n = syn_grids, row.names = 1:length(syn_grids))
	    
	    # matrix with species names and grid numbers:
	    resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
	                             1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))
	    
	    if(length(unique(rownames(resulPaeRaster))) == 0){
	      stop('This iteration had to be stopped because there is no more synapomorphies!')
	    }
	    
	    posic <- 0
	    for(pos in 1:length(unique(rownames(resulPaeRaster)))){
	      posic[pos] <- which(colnames(matric) == unique(rownames(resulPaeRaster))[pos])
	    }
	    
	    # changing the matrix for the next iteration:	
	    tempMatrix <- subset(tempMatrix, select=-posic, site.pattern = FALSE)
	    
	    ###########################
	    
	    ## the number of species supporting the grid number
	    speciesNumber <- data.frame(table(resulPaeRaster))
	    # speciesNumber
	    
	    similarNames <- list()
	    for(volta in 1:dim(speciesNumber)[1]){
	      similarNames[[volta]] <- resulPaeRaster[which(resulPaeRaster[,1] == speciesNumber$resulPaeRaster[volta]), ]
	    }
	    # similarNames
	    
	    ## grid numbers that are generalized tracks in this iteration
	    gridIt <- unique(as.numeric(unlist(similarNames)))
	    # gridIt
	    
	    # set the cells associated with the shapfile to the specified value
	    r[r > 0] <- NA
	    
	    for(j in 1:dim(speciesNumber)[1]){
	      gTrack <- as.numeric(as.character(speciesNumber[j,1]))
	      values(r)[gTrack] <- 1
	    }
	    
	    plot(r, axes = FALSE, legend = FALSE, add = TRUE, col = colores[contagem],
	         alpha = 0.50)
	    
	    # producing a raster:
	    #convert the raster to points for plotting the number of a grid
	    map.r <- as.data.frame(rasterToPoints(r))
	    pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1) # raster com as presenças
	    suppressWarnings(writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
	                format = "GTiff", overwrite = TRUE))
	    
	    if(gridView == TRUE){
	      plot(cropped_map, add = TRUE)
	      
	      map.r$gridNumber <- speciesNumber$resulPaeRaster
	      
	      if(labelGrid == TRUE){
	        text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = (2.5 * resolut[1]) / 10, col = 'red')
	      }
	    }
	    
	    if(legendSpecies == TRUE){
	      print('Please, check the legendSpeciesA.txt file in the /out directory...')
	      #local <- locator()
	      logfile <- "out/legendSpeciesA.txt"
	      cat(c("grid_number", "species", "\n"), file = logfile, sep="\t")
	      for(i in as.numeric(as.character(speciesNumber$resulPaeRaster))){
	        cat(c(paste0(c('grid_', i), collapse = ''), rownames(subset(resulPaeRaster,
	           resulPaeRaster[,'grid_n'] == i)), '\n'), sep = '\t',
	            file = logfile, append = TRUE)
	      }
	    }
	      
	    
	    lista[[contagem]] <- frameTemp
	    
	    print(paste0(c('The iteration number ', contagem, ' has finished.'), collapse = ''))
	    
	  } # close while looping
	  if(is.null(nomesCOL) || dim(as.data.frame(tempMatrix))[1] == 0){
	    x <- seq(1:(contagem - 1))
	  } else {
	    x <- seq(1:contagem)
	  }
	  legend(x = 'topright', legend = x, pch = 15, col = colores[x],
	         title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
	  
	  if(contagem == 1){
	    print('There is no result!')
	  } else if(contagem != 1){
	    return(lista)
	  }
	} # close if... else... condition
}