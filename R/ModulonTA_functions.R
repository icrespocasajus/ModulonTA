suppressMessages(require(igraph))

#' @title Jaccard Distance
#' @description Calculate Jaccard distance between two vectors.
#' @param a A vector.
#' @param b A vector.
#' @return Jaccard distance calculated as the ratio between the intersection and the union of the two vectors.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname jaccard
#' @export 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' @title Redundancy
#' Calculate redundancy between two vectors
#' @param a A vector
#' @param b A vector
#' @return Redundancy calculated as the ratio between the intersection and the minimum length of the two vectors.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' redundancy(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname redundancy
#' @export 
redundancy <- function(a, b) {
  intersection = length(intersect(a, b))
  min = min(length(a),length(b))
  return (intersection/min)
}


#' @title Redundancy with respect to a Reference
#' Calculate redundancy of a vector with respect to another vector.
#' @param a A query vector
#' @param b A reference vector
#' @return Redundancy calculated as the ratio between the intersection and the length of the reference vector
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' redundancy(c('A','B','C'),c('A','B','C','D','E'))
#'  }
#' }
#' @rdname redundancy.wrt
#' @export 
redundancy.wrt <- function(a, b) {
  intersection = length(intersect(a, b))
  ref = length(b)
  return (intersection/ref)
}

#' @title Overlap between two vectors
#' Calculate the number of common elements between two vectors
#' @param a A vector
#' @param b A vector
#' @return The integer with the number of common elements
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' overlap(c('A','B','C'),c('A','B','C','D','E'))
#'  }
#' }
#' @rdname overlap
#' @export 
overlap <- function(a, b) {
  intersection = length(intersect(a, b))
  return (intersection)
}


#' @title Range between 0 and 1
#' @description Normalize absolute values of a numeric vector
#' @param x A numeric vector
#' @return Normalized values (x-min(x))/(max(x)-min(x))
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' range01(c('-50','-5','5','50','100')
#'  }
#' }
#' @rdname range01
#' @export 
range01 <- function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))}




#' @title Split Composite Names
#' @description Split a vector of composite names into a matrix of simple names.
#' @param x character vector
#' @param split character to split each element of vector on, see strsplit
#' @param ... other arguments are passed to strsplit
#' @return Normalized values (x-min(x))/(max(x)-min(x))
#' @details This function is the same as strsplit except that the output value is a matrix instead of a list. 
#' The first column of the matrix contains the first component from each element of x, the second column contains the second components etc.
#'  The number of columns is equal to the maximum number of components for any element of x.
#'   The motivation for this function in the limma package is handle input columns which are composites of two or more annotation fields.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' x = c("2__cc.1","2__cc.2"); strsplit2(x,split="__")
#'  }
#' }
#' @rdname strsplit2
#' @export 
strsplit2 = function (x, split, ...) 
{
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i]) 
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}




#' @title Target Analysis
#' @description Find common targets between a given modulon connected component and modulon constituent elements out of the connected component.
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param modulons A list with as many elements as modulons/clusters containing the constituent elements.
#' @param cc Connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @return List of the modulon constituent elements sharing targets with a given modulon connected component.
#' @details This function is searches for the modulon constituent elements sharing targets with a given modulon connected component, for each connected component and modulon in a given list of connected components split by their source modulon.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs
#' )
#' }
#' }
#' @rdname target.analysis
#' @export 
target.analysis = function(net,modulons,cc){
  network = net
  modulons = modulons
  cc = cc
  regulons = split(network$Target,network$Source)
  target.analysis.results = list()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      query.cc = names(cc[[modulon.tmp]])[j]
      core.tmp = cc[[modulon.tmp]][[query.cc]]
      
      core.targets.tmp = as.character(unlist(regulons[core.tmp]))
      regulons.subset.core = regulons[core.tmp]
      regulons.subset.no.core = regulons[setdiff(modulons[[modulon.tmp]],core.tmp)]
      
      no.core.redundancy.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        redundancy.tmp = redundancy.wrt(tf.targets.tmp,core.targets.tmp)
        return(redundancy.tmp)
      })
      no.core.redundancy.df = data.frame(TF=names(as.data.frame(no.core.redundancy.wrt.core)),Redundancy = as.numeric(as.character(as.data.frame(no.core.redundancy.wrt.core))))
      rownames(no.core.redundancy.df)=no.core.redundancy.df$TF
      no.core.redundancy.df = no.core.redundancy.df[order(no.core.redundancy.df$Redundancy,decreasing = T),]
      
      no.core.similarity.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        redundancy.tmp = jaccard(tf.targets.tmp,core.targets.tmp)
        return(redundancy.tmp)
      })
      no.core.similarity.df = data.frame(TF=names(as.data.frame(no.core.similarity.wrt.core)),Similarity = as.numeric(as.character(as.data.frame(no.core.similarity.wrt.core))))
      rownames(no.core.similarity.df)=no.core.similarity.df$TF
      no.core.similarity.df = no.core.similarity.df[order(no.core.similarity.df$Similarity,decreasing = T),]
      
      no.core.overlap.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        overlap.tmp = overlap(tf.targets.tmp,core.targets.tmp)
        return(overlap.tmp)
      })
      no.core.overlap.df = data.frame(TF=names(as.data.frame(no.core.overlap.wrt.core)),Overlap = as.numeric(as.character(as.data.frame(no.core.overlap.wrt.core))))
      rownames(no.core.overlap.df)=no.core.overlap.df$TF
      no.core.overlap.df = no.core.overlap.df[order(no.core.overlap.df$Overlap,decreasing = T),]
      
      target.analysis.results[[paste(modulon.tmp,query.cc,sep = '__')]]=list(Redundancy=no.core.redundancy.df,Similarity = no.core.similarity.df,Overlap=no.core.overlap.df)
      
    }  
    
  }
  return(target.analysis.results)
}



#' @title Target Analysis manual query
#' @description Find common targets between a given modulon connected component and modulon constituent elements out of the connected component; the specific modulon and connected component have to be specified.
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param modulons A list with as many elements as modulons/clusters containing the constituent elements.
#' @param cc Connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @param query.mod Name of the query modulon
#' @param query.cc Name of the connected component
#' @return List with as many elements as the modulon constituent elements out of a given modulon connected component. Each element includes a dataframe with 4 columns: TF, Redundancy, Similarity and Overlap
#' @details This function is searches for the modulon constituent elements sharing targets with a given modulon connected component.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.manual.query(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs,
#' query.mod = '3',
#' query.cc= 'cc.3')
#' }
#' }
#' @rdname target.analysis.manual.query
#' @export 
target.analysis.manual.query = function(net,modulon,cc,query.mod,query.cc){
  network = net
  modulons = modulon
  cc = cc
  modulon.tmp = query.mod
  core.tmp = cc[[modulon.tmp]][[query.cc]]
  
  regulons = split(network$Target,network$Source)
  core.targets.tmp = as.character(unlist(regulons[core.tmp]))
  regulons.subset.core = regulons[core.tmp]
  regulons.subset.no.core = regulons[setdiff(modulons[[modulon.tmp]],core.tmp)]
  
  no.core.redundancy.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = redundancy.wrt(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.redundancy.df = data.frame(TF=names(as.data.frame(no.core.redundancy.wrt.core)),Redundancy = as.numeric(as.character(as.data.frame(no.core.redundancy.wrt.core))))
  rownames(no.core.redundancy.df)=no.core.redundancy.df$TF
  no.core.redundancy.df = no.core.redundancy.df[order(no.core.redundancy.df$Redundancy,decreasing = T),]
  
  no.core.similarity.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = jaccard(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.similarity.df = data.frame(TF=names(as.data.frame(no.core.similarity.wrt.core)),Similarity = as.numeric(as.character(as.data.frame(no.core.similarity.wrt.core))))
  rownames(no.core.similarity.df)=no.core.similarity.df$TF
  no.core.similarity.df = no.core.similarity.df[order(no.core.similarity.df$Similarity,decreasing = T),]
  
  no.core.overlap.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    overlap.tmp = overlap(tf.targets.tmp,core.targets.tmp)
    return(overlap.tmp)
  })
  no.core.overlap.df = data.frame(TF=names(as.data.frame(no.core.overlap.wrt.core)),Overlap = as.numeric(as.character(as.data.frame(no.core.overlap.wrt.core))))
  rownames(no.core.overlap.df)=no.core.overlap.df$TF
  no.core.overlap.df = no.core.overlap.df[order(no.core.overlap.df$Overlap,decreasing = T),]
  
  target.analysis.results = list()
  target.analysis.results[[paste(query.mod,query.cc,sep = '__')]]=list(Redundancy=no.core.redundancy.df,Similarity = no.core.similarity.df,Overlap=no.core.overlap.df)
  return(target.analysis.results)
}






#' @title Find satellite transcription factors
#' @description Find, among the constituent elements of a given modulon, transcription factors sharing targets with a modulon connected component above a given threshold of similarity, redundancy or overlap (see jaccard(),redundancy.wrt() and overlap() functions)
#' @param data A list object where each element contains a dataframe  including the 'Similarity','Redundancy','Overlap' values of the modulon constituent elements not included within a given modulon connected component.
#' @param feature Feature to be considered for the satellite detection. Possible values: c('Similarity','Redundancy','Overlap').
#' @param threshold Cutoff value for either similarity, redundancy or overlap (number of common targets).
#' @return List object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with the transcription factors above the threshold for the feature considered.
#' @details This function collects the modulon constituent elements sharing targets with a given modulon above a given threshold of similarity, redundancy or overlap.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon.wrt.cc.manual.query(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs,
#' query.mod = '3',
#' query.cc= 'cc.3') %>% Find.Sat()
#' }
#' }
#' @rdname Find.Sat
#' @export 
Find.Sat = function(data,feature = 'Redundancy',threshold = 0,quant.prob = NULL){
  Satellites = lapply(data,function(x){
    x[[feature]][x[[feature]][feature]>threshold,'TF']
  })
  return(Satellites)
}


#' @title Filter satellite transcription factors by discriminant power
#' @description Filter satellite transcription factors by discriminant power values derived from one or more discriminant analysis.
#' @param sat.data List object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param DA.data List object with the discriminant score of all the regulons derived from one or more discriminant analysis.
#' @param DA OPLS-DA to be considered. This argument can take as value any of the names of DA.data list or several of them. If DA = 'Any', all the names of DA.data are considered.
#' @param top.percent Quantile (in %) to be considered from the top of the discriminant score ranking.
#' @return List object with as many elements as the modulon connected components provided in the input. Each element of the list contains a character vector with the corresponding satellite transcription factors that are among the top discriminant regulons from one or more discriminant analysis.
#' @details This function remove satellite transcription factors with low discriminant power.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon.wrt.cc.manual.query(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs,
#' query.mod = '3',
#' query.cc= 'cc.3') %>% Find.Sat() %>% Filter.Sat(DA.data = DA.TILs)
#' }
#' }
#' @rdname Filter.Sat
#' @export 
Filter.Sat = function(sat.data,DA.data,DA=c('Any'),top.percent=10){
  DA.only.top = lapply(DA.data,function(x){
    y = x[x$Weights > quantile(x$Weights,prob = 1-(top.percent/100)), ,drop=FALSE]
    return(y)
  })
  DA.only.top.names = lapply(DA.data,function(x){
    y = rownames(x)[x$Weights > quantile(x$Weights,prob = 1-(top.percent/100))]
    return(y)
  })
  
  if(DA == 'Any'){
    selection = names(DA.data)
  }else{selection = DA}
  
  DA.only.top.names.selected = DA.only.top.names[selection]
  
  DA.selection = unique(as.character(unlist(DA.only.top.names.selected)))
  
  Satellites.Filtered = lapply(sat.data,function(x){
    return(intersect(x,DA.selection))
  })
  return(Satellites.Filtered)
}





#' @title Calculate modulon membership
#' @description Calculation of the correlation (Spearman) between the first principal component (PC1) of a subset of the regulon activity matrix with only modulon constituent elements and each of these elements.
#' @param data Matrix of dataframe with regulon activity.
#' @param mod List object with each modulon constituent elements.
#' @return List object with as many elements as modulons. Each element of the list contains a character vector with the corresponding modulon constituent elements.
#' @details This function subsets the regulon activity matrix to include only modulon constituent elements. After centering and scaling the resulting matrix, a PCA is performed; the PC1 is then correlated with the regulon activity of each modulon constituent element. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#' Modulon.Membership(data = RegAUC.TILs, mod = modulons.TILs)
#' }
#' }
#' @rdname Modulon.Membership
#' @export 
Modulon.Membership = function(data,mod){
  membership.results = lapply(mod,function(mod.tmp){
    prcomp.res = prcomp(data[,mod.tmp],center = T,scale. = T,rank. = 1)
    prcomp.res.df = as.data.frame(prcomp.res[["x"]])
    cor.results = as.data.frame(abs(t(cor(x = prcomp.res.df[rownames(data),'PC1'],y=data[rownames(data),mod.tmp],method = 'spearman'))))
    colnames(cor.results)='R.PC1'
    cor.results = cor.results[order(cor.results[,'R.PC1'],decreasing = T),,drop=F]
    return(cor.results)
  })
  return(membership.results)
}



#' @title Calculate modulon regulatory core membership
#' @description Calculation of the correlation (Spearman) between the first principal component (PC1) of a subset of the regulon activity matrix with only modulon regulatory core constituent elements and each of these elements.
#' @param data Matrix of dataframe with regulon activity.
#' @param mod A character vector with a given modulon constituent elements.
#' @param core A character vector with a given modulon regulatory core constituent elements
#' @return List object with as many elements as modulons. Each element of the list contains a character vector with the corresponding modulon constituent elements.
#' @details This function subsets the regulon activity matrix to include only a given modulon regulatory core constituent elements. After centering and scaling the resulting matrix, a PCA is performed; the PC1 is then correlated with the regulon activity of each modulon constituent element. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#' Core.Membership.manual(data = RegAUC.TILs, mod = modulons.TILs[['3']],core=cc.TILs[['3']][['cc.3']])
#' }
#' }
#' @rdname Core.Membership.manual
#' @export 
Core.Membership.manual = function(data,mod,core){
  prcomp.res = prcomp(data[,core],center = T,scale. = T,rank. = 1)
  prcomp.res.df = as.data.frame(prcomp.res[["x"]])
  cor.results = as.data.frame(abs(t(cor(x = prcomp.res.df[rownames(data),'PC1'],y=data[rownames(data),mod],method = 'spearman'))))
  colnames(cor.results)='R.PC1'
  cor.results = cor.results[order(cor.results[,'R.PC1'],decreasing = T),,drop=F]
  return(cor.results)
}



#' @title Plot modulon target similarity, redundancy or overlap.
#' @description Plot a heatmap displaying modulon target similarity, redundancy or overlap.
#' @param net List R object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param mod List object with each modulon constituent elements.
#' @param cc List with the connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @param regulatory.core Regulatory core names.
#' @param feature Target analysis feature to be displayed; one of c("Redundancy","Similarity","Overlap")
#' @param RegAUC Regulon activity matrix.
#' @return List object with as many elements as modulons; each element contain a heatmap.
#' @details This function generates a heatmap to explore the results of the modulon target analysis
#' @examples 
#' \dontrun{
#' if(interactive()){
#' plots = Modulon.heatmap(
#'  net = network.TILs,
#'  mod = modulons.TILs,
#'  cc = cc.TILs,
#'  regulatory.core = Modulon.Cores.TILs,
#'  feature = 'Redundancy',
#'  RegAUC = RegAUC.TILs)
#'  print(plots[['3']])
#' }
#' }
#' @rdname Modulon.heatmap
#' @export
Modulon.heatmap = function(net,mod,cc,regulatory.core,feature='Redundancy',RegAUC){
  # Libraries
  library(stringr)
  library(pheatmap)
  library(operators)
  
  regulons = split(net$Target,net$Source)
  
  annotation = c()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      cc.tmp = names(cc[[modulon.tmp]])[j]
      for(k in 1:length(cc[[modulon.tmp]][[cc.tmp]])){
        TF.tmp = cc[[modulon.tmp]][[cc.tmp]][k]
        annotation = c(annotation,paste(modulon.tmp,cc.tmp,TF.tmp,sep = '__'))
      }
    } 
  }  
  annotation.df = data.frame(Modulon=ModulonCore::strsplit2(annotation,'__')[,1],
                             cc=ModulonCore::strsplit2(annotation,'__')[,2],
                             TF = ModulonCore::strsplit2(annotation,'__')[,3])
  
  rownames(annotation.df)=annotation.df$TF
  
  annotation.df$id = paste(annotation.df$Modulon,annotation.df$cc,sep = '__')
  annotation.df$Regulatory.Core = ifelse(annotation.df$id %in% regulatory.core,'Yes','Not')
  annotation.df$Regulatory.Core.Annotation = annotation.df$cc
  annotation.df$Regulatory.Core.Annotation[annotation.df$id %!in% regulatory.core]=NA
  
  tab.tmp = table(annotation.df$id) > 1
  true.cc = names(tab.tmp)[tab.tmp ]
  annotation.df[annotation.df$id %!in% true.cc,'cc']=NA
  
  # Add modulon membership
  Modulon.Membership.results = Modulon.Membership(data=RegAUC,mod=mod)
  
  # Generate plots
  plots = list()
  for(h in 1:length(names(mod))){
    name.tmp = names(mod)[h]
    
    regulons.tmp = regulons[mod[[name.tmp]]]
    
    feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
    rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
    colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
    
    for(i in 1:nrow(feature.df)){
      row.tmp = rownames(feature.df)[i]
      for(j in 1:ncol(feature.df)){
        col.tmp = colnames(feature.df)[j]
        if(feature == 'Redundancy'){feature.df[row.tmp,col.tmp]=redundancy(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
        if(feature == 'Similarity'){feature.df[row.tmp,col.tmp]=jaccard(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
        if(feature == 'Overlap'){feature.df[row.tmp,col.tmp]=overlap(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
      }
      
    }
    feature.df=as.matrix(feature.df)
    diag(feature.df)=NA
    
    # Annotation
    
    annotation.c = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Connected.Component=annotation.df[rownames(feature.df),'cc'],
                              Modulon.Membership=Modulon.Membership.results[[name.tmp]][rownames(feature.df),'R.PC1'])
    rownames(annotation.c)=rownames(feature.df)
    
    annotation.r = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Connected.Component=annotation.df[rownames(feature.df),'cc'])
    rownames(annotation.r)=rownames(feature.df)
    
    ann_colors = list(
      Modulon.Membership = c("white", "darkgreen"),
      Regulatory.Core = c(Yes = 'black', Not='white')
    )
    
    # Heatmap input
    phm.input = feature.df[order(annotation.c$Connected.Component,decreasing = T),order(annotation.c$Connected.Component,decreasing = T)]
    
    
    # Gaps connected components
    generate.gaps = function(data,col){
      character.tmp = data[,col,drop=T]
      gaps = c()
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&!(is.na(character.tmp[i+1]))&!(character.tmp[i] == character.tmp[i+1])){gaps = c(gaps,i)}
      }
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&(is.na(character.tmp[i+1]))){gaps = c(gaps,i)}
      }
      return(gaps)
    }
    
    gaps = generate.gaps(data=annotation.c[rownames(phm.input),],col='Connected.Component')
    
    phm.tmp = pheatmap::pheatmap(
      phm.input ,
      main = paste(feature,' Modulon ',name.tmp,sep = ' '),
      display_numbers = F,
      cluster_rows = F,
      cluster_cols = F,
      gaps_row = gaps,
      gaps_col = gaps,
      annotation_col = annotation.c,
      annotation_row = annotation.r,
      annotation_color = ann_colors,
      fontsize_row = 8,
      show_rownames = 8,
      cellwidth = 8,
      cellheight = 8,
      fontsize_col = 8,
      scale='none'
    )
    dev.off()
    gc()
    plots[[name.tmp]]=phm.tmp
  }
  return(plots)
}







#' @title Plot modulon target similarity, redundancy or overlap wrt specific connected components.
#' @description Plot a heatmap displaying modulon target similarity, redundancy or overlap wrt specific connected components.
#' @param net List R object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param mod List object with each modulon constituent elements.
#' @param cc List with the connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @param regulatory.core Regulatory core names.
#' @param feature Target analysis feature to be displayed; one of c("Redundancy","Similarity","Overlap")
#' @param RegAUC Regulon activity matrix.
#' @return List object with as many elements as modulons; each element contain a heatmap.
#' @details This function generates a heatmap to explore the results of the modulon target analysis
#' @examples 
#' \dontrun{
#' if(interactive()){
#' plots = Modulon.heatmap(
#'  net = network.TILs,
#'  mod = modulons.TILs,
#'  cc = cc.TILs,
#'  regulatory.core = Modulon.Cores.TILs,
#'  feature = 'Redundancy',
#'  RegAUC = RegAUC.TILs)
#'  print(plots[['3']])
#' }
#' }
#' @rdname Modulon.heatmap.sat
#' @export
Modulon.heatmap.sat = function(net,mod,cc,regulatory.core,feature='Redundancy',sat,DA.data,DA='Any',RegAUC){
  # Libraries
  library(stringr)
  library(pheatmap)
  library(operators)
  
  regulons = split(net$Target,net$Source)
  
  annotation = c()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      cc.tmp = names(cc[[modulon.tmp]])[j]
      for(k in 1:length(cc[[modulon.tmp]][[cc.tmp]])){
        TF.tmp = cc[[modulon.tmp]][[cc.tmp]][k]
        annotation = c(annotation,paste(modulon.tmp,cc.tmp,TF.tmp,sep = '__'))
      }
    } 
  }  
  annotation.df = data.frame(Modulon=ModulonCore::strsplit2(annotation,'__')[,1],
                             cc=ModulonCore::strsplit2(annotation,'__')[,2],
                             TF = ModulonCore::strsplit2(annotation,'__')[,3])
  
  rownames(annotation.df)=annotation.df$TF
  
  annotation.df$id = paste(annotation.df$Modulon,annotation.df$cc,sep = '__')
  annotation.df$Regulatory.Core = ifelse(annotation.df$id %in% regulatory.core,'Yes','Not')
  annotation.df$Regulatory.Core.Annotation = annotation.df$cc
  annotation.df$Regulatory.Core.Annotation[annotation.df$id %!in% regulatory.core]=NA
  
  tab.tmp = table(annotation.df$id) > 1
  true.cc = names(tab.tmp)[tab.tmp ]
  annotation.df[annotation.df$id %!in% true.cc,'cc']=NA
  
  # Add modulon membership
  Modulon.Membership.results = Modulon.Membership(data=RegAUC,mod=mod)
  
  # Generate plots
  plots = list()
  for(h in 1:length(regulatory.core)){
      mod.tmp = ModulonCore::strsplit2(regulatory.core[h],'__')[,1]
      cc.tmp = ModulonCore::strsplit2(regulatory.core[h],'__')[,2]
      regulons.tmp = regulons[mod[[mod.tmp]]]
      regulons.cc.tmp = regulons[cc[[mod.tmp]][[cc.tmp]]]
      cc.targets.tmp = as.character(unlist(regulons.cc.tmp))
      
      feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
      rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
      colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
      
      for(i in 1:nrow(feature.df)){
        row.tmp = rownames(feature.df)[i]
        for(j in 1:ncol(feature.df)){
          col.tmp = colnames(feature.df)[j]
          if(feature == 'Redundancy'){feature.df[row.tmp,col.tmp]=redundancy(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
          if(feature == 'Similarity'){feature.df[row.tmp,col.tmp]=jaccard(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
          if(feature == 'Overlap'){feature.df[row.tmp,col.tmp]=overlap(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
        }
      }  
      feature.df=as.matrix(feature.df)
      diag(feature.df)=NA
      
      # Highlight only the selected regulatory core
      annotation.df$Regulatory.Core = ifelse(rownames(annotation.df) %in% names(regulons.cc.tmp),'Yes','Not')
      
      # Add core membership to the annotation
      Core.Membership.results = Core.Membership.manual(data=RegAUC,mod=mod[[mod.tmp]],core=cc[[mod.tmp]][[cc.tmp]] )
      

      
      # Annotation
      
      annotation.c = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                                Regulatory.Core.Satellite = ifelse(rownames(feature.df) %in% sat[[paste(mod.tmp,cc.tmp,sep = '__')]],'Yes','Not'),
                                Connected.Component=annotation.df[rownames(feature.df),'cc'],
                                Regulatory.Core.Membership = Core.Membership.results[rownames(feature.df),'R.PC1'],
                                Modulon.Membership=Modulon.Membership.results[[mod.tmp]][rownames(feature.df),'R.PC1'])
      rownames(annotation.c)=rownames(feature.df)
      
      annotation.r = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                                Regulatory.Core.Satellite = ifelse(rownames(feature.df) %in% sat[[paste(mod.tmp,cc.tmp,sep = '__')]],'Yes','Not'),
                                Connected.Component=annotation.df[rownames(feature.df),'cc'])
      rownames(annotation.r)=rownames(feature.df)
      
      ann_colors = list(
        Modulon.Membership = c("white", "darkgreen"),
        #Discriminant.Power = c("white", "black"),
        Regulatory.Core.Membership = c("white", "firebrick"),
        Regulatory.Core = c(Yes = 'black', Not='white'),
        Regulatory.Core.Satellite = c(Yes = 'red', Not='white')
      )
      
      # Add DA
      
      for(da in 1:length(names(DA.data))){
        name.tmp = names(DA.data)[da]
        annotation.c[,name.tmp] = DA.data[[name.tmp]][rownames(feature.df),'Weights']
      }
      

      my.breaks <- c(seq(-1, 1, by=0.001)) 
      my.colors <- colorRampPalette(colors = c("white", "black"))(length(my.breaks))
      my.colors.df = data.frame(colors = my.colors,breaks=as.numeric(round(my.breaks,digits = 3)))

      for(da in 1:length(names(DA.data))){
        name.tmp = names(DA.data)[da]
        max.tmp=as.numeric(round(max(annotation.c[,name.tmp]),digits = 3))
        min.tmp=as.numeric(round(min(annotation.c[,name.tmp]),digits = 3))
        max.tmp.color=my.colors.df[(my.colors.df$breaks == max.tmp),'colors']
        min.tmp.color=my.colors.df[(my.colors.df$breaks == min.tmp),'colors']
        ann_colors[[name.tmp]] = c(min.tmp.color, max.tmp.color)
      }
      # Heatmap input
      phm.input = feature.df[order(annotation.c$Connected.Component,decreasing = T),order(annotation.c$Connected.Component,decreasing = T)]
      
      
      # Gaps connected components
      generate.gaps = function(data,col){
        character.tmp = data[,col,drop=T]
        gaps = c()
        for(i in 1:length(character.tmp)){
          if(!(is.na(character.tmp[i]))&!(is.na(character.tmp[i+1]))&!(character.tmp[i] == character.tmp[i+1])){gaps = c(gaps,i)}
        }
        for(i in 1:length(character.tmp)){
          if(!(is.na(character.tmp[i]))&(is.na(character.tmp[i+1]))){gaps = c(gaps,i)}
        }
        return(gaps)
      }
      
      gaps = generate.gaps(data=annotation.c[rownames(phm.input),],col='Connected.Component')
      
      phm.tmp = pheatmap::pheatmap(
        phm.input ,
        main = paste(feature,' Modulon ',mod.tmp,'core ',cc.tmp,sep = ' '),
        display_numbers = F,
        cluster_rows = F,
        cluster_cols = F,
        gaps_row = gaps,
        gaps_col = gaps,
        annotation_col = annotation.c,
        annotation_row = annotation.r,
        annotation_color = ann_colors,
        fontsize_row = 8,
        show_rownames = 8,
        cellwidth = 8,
        cellheight = 8,
        fontsize_col = 8,
        scale='none'
      )
      dev.off()
      gc()
      plots[[paste(mod.tmp,cc.tmp,sep = '__')]]=phm.tmp
  }
  return(plots)
}










#' @title Plot modulon target similarity, redundancy or overlap.
#' @description Plot a heatmap displaying modulon target similarity, redundancy or overlap.
#' @param net List R object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param mod List object with each modulon constituent elements.
#' @param cc List with the connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @param regulatory.core Regulatory core names.
#' @param feature Target analysis feature to be displayed; one of c("Redundancy","Similarity","Overlap")
#' @param RegAUC Regulon activity matrix.
#' @return List object with as many elements as modulons; each element contain a heatmap.
#' @details This function generates a heatmap to explore the results of the modulon target analysis
#' @examples 
#' \dontrun{
#' if(interactive()){
#' plots = Modulon.corrplot(
#'  net = network.TILs,
#'  mod = modulons.TILs,
#'  cc = cc.TILs,
#'  regulatory.core = Modulon.Cores.TILs,
#'  feature = 'Redundancy',
#'  RegAUC = RegAUC.TILs)
#'  print(plots[['3']])
#' }
#' }
#' @rdname Modulon.corrplot
#' @export
Modulon.corrplot = function(net,mod,cc,regulatory.core,feature='Redundancy',RegAUC){
  # Libraries
  library(stringr)
  library(corrplot)
  library(RColorBrewer)
  library(operators)
  
  regulons = split(net$Target,net$Source)
  
  annotation = c()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      cc.tmp = names(cc[[modulon.tmp]])[j]
      for(k in 1:length(cc[[modulon.tmp]][[cc.tmp]])){
        TF.tmp = cc[[modulon.tmp]][[cc.tmp]][k]
        annotation = c(annotation,paste(modulon.tmp,cc.tmp,TF.tmp,sep = '__'))
      }
    } 
  }  
  annotation.df = data.frame(Modulon=ModulonCore::strsplit2(annotation,'__')[,1],
                             cc=ModulonCore::strsplit2(annotation,'__')[,2],
                             TF = ModulonCore::strsplit2(annotation,'__')[,3])
  
  rownames(annotation.df)=annotation.df$TF
  
  annotation.df$id = paste(annotation.df$Modulon,annotation.df$cc,sep = '__')
  annotation.df$Regulatory.Core = ifelse(annotation.df$id %in% regulatory.core,'Yes','Not')
  annotation.df$Regulatory.Core.Annotation = annotation.df$cc
  annotation.df$Regulatory.Core.Annotation[annotation.df$id %!in% regulatory.core]=NA
  
  tab.tmp = table(annotation.df$id) > 1
  true.cc = names(tab.tmp)[tab.tmp ]
  annotation.df[annotation.df$id %!in% true.cc,'cc']=NA
  
  # Add modulon membership
  Modulon.Membership.results = Modulon.Membership(data=RegAUC,mod=mod)
  
  # Generate plots
  plots = list()
  for(h in 1:length(names(mod))){
    name.tmp = names(mod)[h]
    
    regulons.tmp = regulons[mod[[name.tmp]]]
    
    feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
    rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
    colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
    
    for(i in 1:nrow(feature.df)){
      row.tmp = rownames(feature.df)[i]
      for(j in 1:ncol(feature.df)){
        col.tmp = colnames(feature.df)[j]
        if(feature == 'Redundancy'){feature.df[row.tmp,col.tmp]=redundancy(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
        if(feature == 'Similarity'){feature.df[row.tmp,col.tmp]=jaccard(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
        if(feature == 'Overlap'){feature.df[row.tmp,col.tmp]=overlap(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])}
      }
      
    }
    feature.df=as.matrix(feature.df)
    diag(feature.df)=NA
    
    # Annotation
    
    annotation.c = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Connected.Component=annotation.df[rownames(feature.df),'cc'],
                              Modulon.Membership=Modulon.Membership.results[[name.tmp]][rownames(feature.df),'R.PC1'])
    rownames(annotation.c)=rownames(feature.df)
    
    annotation.r = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Connected.Component=annotation.df[rownames(feature.df),'cc'])
    rownames(annotation.r)=rownames(feature.df)
    
    ann_colors = list(
      Modulon.Membership = c("white", "darkgreen"),
      Regulatory.Core = c(Yes = 'black', Not='white')
    )
    
    # Heatmap input
    phm.input = feature.df[order(annotation.c$Connected.Component,decreasing = T),order(annotation.c$Connected.Component,decreasing = T)]
    
    
    # Gaps connected components
    generate.gaps = function(data,col){
      character.tmp = data[,col,drop=T]
      gaps = c()
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&!(is.na(character.tmp[i+1]))&!(character.tmp[i] == character.tmp[i+1])){gaps = c(gaps,i)}
      }
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&(is.na(character.tmp[i+1]))){gaps = c(gaps,i)}
      }
      return(gaps)
    }
    
    gaps = generate.gaps(data=annotation.c[rownames(phm.input),],col='Connected.Component')
    
    phm.tmp = corrplot::corrplot(
      phm.input,
      type="upper",
      order="hclust",
      col=RColorBrewer::brewer.pal(n=8, name="RdYlBu")
    )  
    dev.off()
    gc()
    plots[[name.tmp]]=phm.tmp
  }
  return(plots)
}




#' @title Modulon target analysis: similarity, redundancy and overlap
#' @description Generate the matrix with all pair-wise comparisons within a given modulon for each of the following features:"Redundancy","Similarity" and "Overlap".
#' @param net List R object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param mod List object with each modulon constituent elements.
#' @param mod.query The name of the modulon to be analysed.
#' @return List object with three objects ("Redundancy","Similarity" and "Overlap") with the corresponding matrices.
#' @details This function perform the modulon target analysis.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon(net=network.TILs,mod=modulons.TILs,mod.query = '3')
#' }
#' }
#' @rdname target.analysis.modulon
#' @export
target.analysis.modulon = function(net,mod,mod.query){
  library(stringr)

  regulons = split(net$Target,net$Source)
  name.tmp = mod.query
  regulons.tmp = regulons[mod[[name.tmp]]]
  
  feature.df.list = list()
  
  # Redundancy
  feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
  rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
  colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
  for(i in 1:nrow(feature.df)){
    row.tmp = rownames(feature.df)[i]
    for(j in 1:ncol(feature.df)){
      col.tmp = colnames(feature.df)[j]
      feature.df[row.tmp,col.tmp]=redundancy(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])
    }
  }
  feature.df=as.matrix(feature.df)
  diag(feature.df)=0
  feature.df.list[['Redundancy']]=feature.df 
  
  # Similarity
  feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
  rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
  colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
  for(i in 1:nrow(feature.df)){
    row.tmp = rownames(feature.df)[i]
    for(j in 1:ncol(feature.df)){
      col.tmp = colnames(feature.df)[j]
      feature.df[row.tmp,col.tmp]=jaccard(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])
    }
  }
  feature.df=as.matrix(feature.df)
  diag(feature.df)=0
  feature.df.list[['Similarity']]=feature.df 
  
  
  # Overlap
  feature.df = as.data.frame(matrix(NA,length(names(regulons.tmp)),length(names(regulons.tmp))))
  rownames(feature.df)=  stringr::str_sort( names(regulons.tmp),numeric = T)
  colnames(feature.df)= stringr::str_sort( names(regulons.tmp),numeric = T)
  for(i in 1:nrow(feature.df)){
    row.tmp = rownames(feature.df)[i]
    for(j in 1:ncol(feature.df)){
      col.tmp = colnames(feature.df)[j]
      feature.df[row.tmp,col.tmp]=overlap(regulons.tmp[[row.tmp]],regulons.tmp[[col.tmp]])
    }
  }
  feature.df=as.matrix(feature.df)
  diag(feature.df)=0
  feature.df.list[['Overlap']]=feature.df 
  
  return(feature.df.list)
}

#' @title Plot results of the modulon target analysis
#' @description Generate the plot for the matrix with all pair-wise comparisons within a given modulon for one of the following features:"Redundancy","Similarity" and "Overlap".
#' @param data List object with three elements ("Redundancy","Similarity" and "Overlap") with the corresponding matrices.
#' @param feature Target analysis feature to be displayed; one of c("Redundancy","Similarity","Overlap").
#' @param core Character vector with modulon regulatory core constituent elements.
#' @param core.sat Character vector with modulon regulatory core detected satellites.
#' @return A corrplot() object.
#' @details This function generates a plot showing the results of the modulon target analysis of a given modulon.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' results.target.analysis.modulon=target.analysis.modulon(net=network.TILs,mod=modulons.TILs,mod.query = '3')
#' target.analysis.modulon.plot(data=results.target.analysis.modulon,feature = 'Redundancy')
#' }
#' }
#' @rdname target.analysis.modulon.plot
#' @export
target.analysis.modulon.plot = function(data,feature='Redundancy',color = 'YlGn',core=c(),core.sat=c()){
  library(corrplot)
  feature.df = data[[feature]]
  if(feature == 'Redundancy'){color = 'YlGn'}
  if(feature == 'Similarity'){color = 'Purples'}
  if(feature == 'Overlap'){color = 'YlOrBr'}
  return(corrplot::corrplot(
    feature.df,
    type = 'lower',
    method = 'pie' , 
    order = 'alphabet', 
    tl.col = ifelse(rownames(feature.df)%in%core,'black',ifelse(rownames(feature.df)%in%core.sat,'red','darkgrey')),
    cl.ratio = 0.2,
    tl.srt = 45,
    col = corrplot::COL1(color, 10),
    is.corr = F,
    title = paste("\n\n\n",'Modulon ',feature,sep = ""),tl.cex=0.75))
}


















#' @title Target Analysis wrt modulon connected components
#' @description Find common targets between a given modulon connected component and modulon constituent elements out of the connected component.
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param modulons A list with as many elements as modulons/clusters containing the constituent elements.
#' @param cc Connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @return List of the modulon constituent elements sharing targets with a given modulon connected component.
#' @details This function is searches for the modulon constituent elements sharing targets with a given modulon connected component, for each connected component and modulon in a given list of connected components split by their source modulon.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon.wrt.cc(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs
#' )
#' }
#' }
#' @rdname target.analysis.modulon.wrt.cc
#' @export 
target.analysis.modulon.wrt.cc = function(net,modulons,cc){
  network = net
  modulons = modulons
  cc = cc
  regulons = split(network$Target,network$Source)
  target.analysis.results = list()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      query.cc = names(cc[[modulon.tmp]])[j]
      core.tmp = cc[[modulon.tmp]][[query.cc]]
      
      core.targets.tmp = as.character(unlist(regulons[core.tmp]))
      regulons.subset.core = regulons[core.tmp]
      regulons.subset.no.core = regulons[setdiff(modulons[[modulon.tmp]],core.tmp)]
      
      no.core.redundancy.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        redundancy.tmp = redundancy.wrt(tf.targets.tmp,core.targets.tmp)
        return(redundancy.tmp)
      })
      no.core.redundancy.df = data.frame(TF=names(as.data.frame(no.core.redundancy.wrt.core)),Redundancy = as.numeric(as.character(as.data.frame(no.core.redundancy.wrt.core))))
      rownames(no.core.redundancy.df)=no.core.redundancy.df$TF
      no.core.redundancy.df = no.core.redundancy.df[order(no.core.redundancy.df$Redundancy,decreasing = T),]
      
      no.core.similarity.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        redundancy.tmp = jaccard(tf.targets.tmp,core.targets.tmp)
        return(redundancy.tmp)
      })
      no.core.similarity.df = data.frame(TF=names(as.data.frame(no.core.similarity.wrt.core)),Similarity = as.numeric(as.character(as.data.frame(no.core.similarity.wrt.core))))
      rownames(no.core.similarity.df)=no.core.similarity.df$TF
      no.core.similarity.df = no.core.similarity.df[order(no.core.similarity.df$Similarity,decreasing = T),]
      
      no.core.overlap.wrt.core = lapply(regulons.subset.no.core,function(x){
        tf.targets.tmp = as.character(unlist(x))
        overlap.tmp = overlap(tf.targets.tmp,core.targets.tmp)
        return(overlap.tmp)
      })
      no.core.overlap.df = data.frame(TF=names(as.data.frame(no.core.overlap.wrt.core)),Overlap = as.numeric(as.character(as.data.frame(no.core.overlap.wrt.core))))
      rownames(no.core.overlap.df)=no.core.overlap.df$TF
      no.core.overlap.df = no.core.overlap.df[order(no.core.overlap.df$Overlap,decreasing = T),]
      
      target.analysis.results[[paste(modulon.tmp,query.cc,sep = '__')]]=list(Redundancy=no.core.redundancy.df,Similarity = no.core.similarity.df,Overlap=no.core.overlap.df)
      
    }  
    
  }
  return(target.analysis.results)
}



#' @title Target Analysis with respect to modulon connected components: manual query
#' @description Find common targets between a given modulon connected component and modulon constituent elements out of the connected component; the specific modulon and connected component have to be specified.
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param modulons A list with as many elements as modulons/clusters containing the constituent elements.
#' @param cc Connected components generated with find.connected.components() or regulatory cores as the output of core().
#' @param query.mod Name of the query modulon.
#' @param query.cc Name of the connected component.
#' @return List with as many elements as the modulon constituent elements out of a given modulon connected component. Each element includes a dataframe with 4 columns: TF, Redundancy, Similarity and Overlap
#' @details This function is searches for the modulon constituent elements sharing targets with a given modulon connected component.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon.wrt.cc.manual.query(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs,
#' query.mod = '3',
#' query.cc= 'cc.3')
#' }
#' }
#' @rdname target.analysis.modulon.wrt.cc.manual.query
#' @export 
target.analysis.modulon.wrt.cc.manual.query = function(net,modulon,cc,mod.query,cc.query){
  network = net
  modulons = modulon
  cc = cc
  modulon.tmp = mod.query
  core.tmp = cc[[modulon.tmp]][[cc.query]]
  
  regulons = split(network$Target,network$Source)
  core.targets.tmp = as.character(unlist(regulons[core.tmp]))
  regulons.subset.core = regulons[core.tmp]
  regulons.subset.no.core = regulons[setdiff(modulons[[modulon.tmp]],core.tmp)]
  
  no.core.redundancy.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = redundancy.wrt(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.redundancy.df = data.frame(TF=names(as.data.frame(no.core.redundancy.wrt.core)),Redundancy = as.numeric(as.character(as.data.frame(no.core.redundancy.wrt.core))))
  rownames(no.core.redundancy.df)=no.core.redundancy.df$TF
  no.core.redundancy.df = no.core.redundancy.df[order(no.core.redundancy.df$Redundancy,decreasing = T),]
  
  no.core.similarity.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = jaccard(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.similarity.df = data.frame(TF=names(as.data.frame(no.core.similarity.wrt.core)),Similarity = as.numeric(as.character(as.data.frame(no.core.similarity.wrt.core))))
  rownames(no.core.similarity.df)=no.core.similarity.df$TF
  no.core.similarity.df = no.core.similarity.df[order(no.core.similarity.df$Similarity,decreasing = T),]
  
  no.core.overlap.wrt.core = lapply(regulons.subset.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    overlap.tmp = overlap(tf.targets.tmp,core.targets.tmp)
    return(overlap.tmp)
  })
  no.core.overlap.df = data.frame(TF=names(as.data.frame(no.core.overlap.wrt.core)),Overlap = as.numeric(as.character(as.data.frame(no.core.overlap.wrt.core))))
  rownames(no.core.overlap.df)=no.core.overlap.df$TF
  no.core.overlap.df = no.core.overlap.df[order(no.core.overlap.df$Overlap,decreasing = T),]
  
  target.analysis.results = list()
  target.analysis.results[[paste(modulon.tmp,cc.query,sep = '__')]]=list(Redundancy=no.core.redundancy.df,Similarity = no.core.similarity.df,Overlap=no.core.overlap.df)
  return(target.analysis.results)
}














#' @title Plot modulon target similarity, redundancy or overlap wrt specific connected components.
#' @description Plot a heatmap displaying modulon target similarity, redundancy or overlap wrt specific connected components.
#' @param net List R object with as many elements as the modulon connected components provided as the input. Each element of the list contains a character vector with satellite transcription factors.
#' @param mod List object with each modulon constituent elements.
#' @param cc List with the connected components generated with find.connected.components() or regulatory cores as the output of core()
#' @param mod.query Name of the query modulon.
#' @param cc.query Name of the connected component considered as modulon regulatory core.
#' @param feature Target analysis feature to be displayed; one of c("Redundancy","Similarity","Overlap")
#' @param RegAUC Regulon activity matrix.
#' @return List object with as many elements as modulons; each element contain a heatmap.
#' @details This function generates a heatmap to explore the results of the modulon target analysis
#' @examples 
#' \dontrun{
#' if(interactive()){
#' plots = Modulon.heatmap(
#'  net = network.TILs,
#'  mod = modulons.TILs,
#'  cc = cc.TILs,
#'  regulatory.core = Modulon.Cores.TILs,
#'  feature = 'Redundancy',
#'  sat = satellites.filtered,
#'  RegAUC = RegAUC.TILs,
#'  color= 'YlGn')
#'  print(plots[['3']])
#' }
#' }
#' @rdname Modulon.complexheatmap
#' @export
Modulon.complexheatmap = function(net,mod,cc,mod.query,cc.query,feature='Redundancy',sat,DA.data,DA='Any',RegAUC,color= 'YlGn'){
  # Libraries
  library(stringr)
  library(operators)
  library(ComplexHeatmap)
  library(corrplot)
  library(randomcoloR)
  
  if(feature == 'Redundancy'){color = 'YlGn'}
  if(feature == 'Similarity'){color = 'Purples'}
  if(feature == 'Overlap'){color = 'YlOrBr'}
  
  
  regulons = split(net$Target,net$Source)
  regulatory.core = cc[[mod.query]][[cc.query]]
  core.sat = satellites.filtered[[paste(mod.query,cc.query,sep = '__')]]
  
  annotation = c()
  for(i in 1:length(names(cc))){
    modulon.tmp = names(cc)[i]
    for(j in 1:length(names(cc[[modulon.tmp]]))){
      cc.tmp = names(cc[[modulon.tmp]])[j]
      for(k in 1:length(cc[[modulon.tmp]][[cc.tmp]])){
        TF.tmp = cc[[modulon.tmp]][[cc.tmp]][k]
        annotation = c(annotation,paste(modulon.tmp,cc.tmp,TF.tmp,sep = '__'))
      }
    } 
  }  
  annotation.df = data.frame(Modulon=ModulonCore::strsplit2(annotation,'__')[,1],
                             cc=ModulonCore::strsplit2(annotation,'__')[,2],
                             TF = ModulonCore::strsplit2(annotation,'__')[,3])
  
  rownames(annotation.df)=annotation.df$TF
  
  # Select modulon of interest
  
  annotation.df = annotation.df[annotation.df$Modulon==mod.query,]
  
  annotation.df$id = paste(annotation.df$Modulon,annotation.df$cc,sep = '__')
  annotation.df$Regulatory.Core = ifelse(annotation.df$TF %in% regulatory.core,'Yes','Not')
  annotation.df$Regulatory.Core.Annotation = annotation.df$cc
  annotation.df$Regulatory.Core.Annotation[annotation.df$TF %!in% regulatory.core]=NA
  
  annotation.df$Satellites = ifelse(annotation.df$TF %in% core.sat,'Yes','Not')
  
  
  tab.tmp = table(annotation.df$id) > 1
  true.cc = names(tab.tmp)[tab.tmp ]
  annotation.df[annotation.df$id %!in% true.cc,'cc']=NA
  
  # Add modulon membership
  Modulon.Membership.results = Modulon.Membership(data=RegAUC,mod=mod)
  
  # Generate plots
  #plots = list()
  #for(h in 1:length(regulatory.core)){
    mod.tmp = mod.query
    cc.tmp = cc.query
    regulons.tmp = regulons[mod[[mod.tmp]]]
    regulons.cc.tmp = regulons[cc[[mod.tmp]][[cc.tmp]]]
    cc.targets.tmp = as.character(unlist(regulons.cc.tmp))
    
    feature.df = results.target.analysis.modulon[[feature]]
    feature.df=as.matrix(feature.df)
    diag(feature.df)=NA
    
    
    # Add core membership to the annotation
    Core.Membership.results = Core.Membership.manual(data=RegAUC,mod=mod[[mod.tmp]],core=cc[[mod.tmp]][[cc.tmp]] )
    
    
    
    # Annotation
    
    annotation.c = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Regulatory.Core.Satellite = ifelse(rownames(feature.df) %in% sat[[paste(mod.tmp,cc.tmp,sep = '__')]],'Yes','Not'),
                              Connected.Component=annotation.df[rownames(feature.df),'cc'],
                              Regulatory.Core.Membership = Core.Membership.results[rownames(feature.df),'R.PC1']#,
                              #Modulon.Membership=Modulon.Membership.results[[mod.tmp]][rownames(feature.df),'R.PC1']
                              )
    rownames(annotation.c)=rownames(feature.df)
    
    annotation.r = data.frame(Regulatory.Core = annotation.df[rownames(feature.df),'Regulatory.Core'],
                              Regulatory.Core.Satellite = ifelse(rownames(feature.df) %in% sat[[paste(mod.tmp,cc.tmp,sep = '__')]],'Yes','Not')#,
                              #Connected.Component=annotation.df[rownames(feature.df),'cc']
                              )
    rownames(annotation.r)=rownames(feature.df)
    
  cc.color = distinctColorPalette(k = length(unique(annotation.c$Connected.Component)))
    ann_colors = list(
      Modulon.Membership = c("white", "darkgreen"),
      Connnected.Component = cc.color,
      #Discriminant.Power = c("white", "black"),
      Regulatory.Core.Membership = c("white", "firebrick"),
      Regulatory.Core = c(Yes = 'black', Not='grey'),
      Regulatory.Core.Satellite = c(Yes = 'red', Not='grey')
    )
    
    my.breaks <- c(seq(-1, 1, by=0.001)) 
    my.colors <- colorRampPalette(colors = c("white", "black"))(length(my.breaks))
    my.colors.df = data.frame(colors = my.colors,breaks=as.numeric(round(my.breaks,digits = 3)))
    
    
    # Heatmap input
    phm.input = feature.df[order(annotation.c$Connected.Component,decreasing = T),order(annotation.c$Connected.Component,decreasing = T)]
    annotation.c=annotation.c[rownames( phm.input),]
    annotation.r=annotation.r[rownames( phm.input),]
    
    # Gaps connected components
    generate.gaps = function(data,col){
      character.tmp = data[,col,drop=T]
      gaps = c()
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&!(is.na(character.tmp[i+1]))&!(character.tmp[i] == character.tmp[i+1])){gaps = c(gaps,i)}
      }
      for(i in 1:length(character.tmp)){
        if(!(is.na(character.tmp[i]))&(is.na(character.tmp[i+1]))){gaps = c(gaps,i)}
      }
      return(gaps)
    }
    
    gaps = generate.gaps(data=annotation.c[rownames(phm.input),],col='Connected.Component')
    
    #Edit for ComplexHeatmap
    annotation.c$Regulatory.Core[is.na(annotation.c$Regulatory.Core)] = 'Not'
    annotation.r$Regulatory.Core[is.na(annotation.c$Regulatory.Core)] = 'Not'
    
    phm.tmp = ComplexHeatmap::pheatmap(
      name = feature,
      legend = T,
     # annotation_legend = c(T,T,T,T),
      color = corrplot::COL1(color, 10),
      phm.input ,
      main = paste(feature,' Modulon ',mod.tmp,'core ',cc.tmp,sep = ' '),
      display_numbers = F,
      cluster_rows = F,
      cluster_cols = F,
      gaps_row = gaps,
      gaps_col = gaps,
      annotation_col = annotation.c,
      annotation_row = annotation.r,
      annotation_color = ann_colors,
      fontsize_row = 8,
      show_rownames = 8,
      cellwidth = 8,
      cellheight = 8,
      fontsize_col = 8,
      scale='none'
    )
    
    mat = results.target.analysis.modulon[[feature]]
    mat=mat[rownames(phm.input),colnames(phm.input)]
    feature.wrt.core = results.target.analysis.modulon.wrt.cc.w.core[[paste(mod.query,cc.query,sep = '__')]][[feature]][rownames(mat),feature]
    row_ha = rowAnnotation(Feature = anno_barplot(feature.wrt.core),annotation_legend_param = list(just = c("left", "bottom")))
    
    
    # DA heatmap
    annotation.c.DA = annotation.c
    for(da in 1:length(names(DA.data))){
      name.tmp = names(DA.data)[da]
      annotation.c.DA[,name.tmp] = DA.data[[name.tmp]][rownames(feature.df),'Weights']
    }
    heatmap.DA = as.matrix(annotation.c.DA[,names(DA.data),drop=F])
    colnames(heatmap.DA)=gsub('_Vs_background','',colnames(heatmap.DA))
    
    phm.DA.tmp = ComplexHeatmap::pheatmap(
      name = 'OPLS-DA',
      #legend = F,
      main = 'OPLS-DA',
      color = corrplot::COL1('Greys', 10),
      heatmap.DA ,
      #main = paste(feature,' Modulon ',mod.tmp,'Discriminant Power',sep = ' '),
      display_numbers = F,
      cluster_rows = F,
      cluster_cols = F,
      gaps_row = gaps,
      #gaps_col = gaps,
      #annotation_col = annotation.c,
      #annotation_row = annotation.r,
      #annotation_color = ann_colors,
      fontsize_row = 8,
      show_rownames = 8,
      cellwidth = 8,
      cellheight = 8,
      fontsize_col = 8,
      scale='none'
    )
    
   # phm.input2 = phm.input
  #  phm.input2[is.na(phm.input2)]=0
    
    Heatmap3D(phm.input2, name = "mat", column_title = "This is a 3D heatmap",col=c('white','red'),cluster_columns = F, cluster_rows=F)
    
    
    
    
    plot.tmp = phm.tmp + row_ha
    plot.listt = phm.tmp + row_ha + phm.DA.tmp
    plot.list = phm.DA.tmp + row_ha + phm.tmp
    
    plot.list = phm.DA.tmp + phm.tmp
    
    #draw(final.plot, auto_adjust = c(F))
   plot= draw(plot.list, auto_adjust = F, padding = unit(c(2, 2, 10, 2), "mm"),heatmap_legend_side = "bottom", annotation_legend_side = "right",merge_legend = TRUE)
   plot= draw(plot.list, auto_adjust = F, padding = unit(c(2, 2, 10, 2), "mm"),heatmap_legend_side = "bottom", annotation_legend_side = "right")
   
    decorate_heatmap_body("OPLS-DA", {grid.text("OPLS-DA", y = unit(1, "npc") + unit(2, "mm"), just = "bottom")})
   
    
    dev.off()
    gc()
   
  return(plot)
}









#' @title Target Analysis with respect to modulon connected components: manual query
#' @description Find common targets between a given modulon connected component and modulon constituent elements out of the connected component; the specific modulon and connected component have to be specified.
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param mod A list with as many elements as modulons/clusters containing the constituent elements.
#' @param cc Connected components generated with find.connected.components() or regulatory cores as the output of core().
#' @param mod.query Name of the query modulon.
#' @param cc.query Name of the connected component.
#' @return List with as many elements as the modulon constituent elements out of a given modulon connected component. Each element includes a dataframe with 4 columns: TF, Redundancy, Similarity and Overlap
#' @details This function is searches for the modulon constituent elements sharing targets with a given modulon connected component.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' target.analysis.modulon.wrt.cc.manual.query(net = network.TILs,
#' mod = modulons.TILs,
#' cc = cc.TILs,
#' mod.query = '3',
#' cc.query = 'cc.3')
#' }
#' }
#' @rdname target.analysis.modulon.wrt.cc.manual.query.2
#' @export 
target.analysis.modulon.wrt.cc.manual.query.2 = function(net,mod,cc,mod.query,cc.query){
  network = net
  modulons = mod
  cc = cc
  modulon.tmp = mod.query
  core.tmp = cc[[modulon.tmp]][[cc.query]]
  
  regulons = split(network$Target,network$Source)
  core.targets.tmp = as.character(unlist(regulons[core.tmp]))
  regulons.subset.core = regulons[core.tmp]
  regulons.subset.no.core = regulons[setdiff(modulons[[modulon.tmp]],core.tmp)]
  regulons.subset.core.and.no.core = c(regulons.subset.core,regulons.subset.no.core)
  
  no.core.redundancy.wrt.core = lapply(regulons.subset.core.and.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = redundancy.wrt(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.redundancy.df = data.frame(TF=names(as.data.frame(no.core.redundancy.wrt.core)),Redundancy = as.numeric(as.character(as.data.frame(no.core.redundancy.wrt.core))))
  rownames(no.core.redundancy.df)=no.core.redundancy.df$TF
  no.core.redundancy.df = no.core.redundancy.df[order(no.core.redundancy.df$Redundancy,decreasing = T),]
  
  no.core.similarity.wrt.core = lapply(regulons.subset.core.and.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    redundancy.tmp = jaccard(tf.targets.tmp,core.targets.tmp)
    return(redundancy.tmp)
  })
  no.core.similarity.df = data.frame(TF=names(as.data.frame(no.core.similarity.wrt.core)),Similarity = as.numeric(as.character(as.data.frame(no.core.similarity.wrt.core))))
  rownames(no.core.similarity.df)=no.core.similarity.df$TF
  no.core.similarity.df = no.core.similarity.df[order(no.core.similarity.df$Similarity,decreasing = T),]
  
  no.core.overlap.wrt.core = lapply(regulons.subset.core.and.no.core,function(x){
    tf.targets.tmp = as.character(unlist(x))
    overlap.tmp = overlap(tf.targets.tmp,core.targets.tmp)
    return(overlap.tmp)
  })
  no.core.overlap.df = data.frame(TF=names(as.data.frame(no.core.overlap.wrt.core)),Overlap = as.numeric(as.character(as.data.frame(no.core.overlap.wrt.core))))
  rownames(no.core.overlap.df)=no.core.overlap.df$TF
  no.core.overlap.df = no.core.overlap.df[order(no.core.overlap.df$Overlap,decreasing = T),]
  
  target.analysis.results = list()
  target.analysis.results[[paste(modulon.tmp,cc.query,sep = '__')]]=list(Redundancy=no.core.redundancy.df,Similarity = no.core.similarity.df,Overlap=no.core.overlap.df)
  return(target.analysis.results)
}
