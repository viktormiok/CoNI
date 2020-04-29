#CoNI Functions
############


####################################
#Functions
####################################
CoNI<- function(data1, data2,data1name="data1",data2name="data2", padjustData2=TRUE, correlateDFs=TRUE,splitData1=TRUE,old_split=NULL,split_number=2,outputDir="./CoNIOutput/",iteration_start=1,wait_iteration=0, numCores=6) { #It could be something else

  #Check if input objects are defined
  do_objectsExist(data1,data2)

  #Start measuring time
  start_time <- Sys.time()

  #Check if output directory exists
  check_outputDir(outputDir)

  #Test if sample names are the same in both data sets
  compare_sampleNames(data1,data2)

  #Make sure column names are appropiate
  colnames(data1)<-make.names(colnames(data1),unique=TRUE)
  colnames(data2)<-make.names(colnames(data2),unique=TRUE)

  if(!file.exists(paste(data2name,"_Tablesignificant.csv",sep=""))){
    #Get significant correlations between metabolites
    print("Calculating correlations of data2")
    normMetabo_Tablesignificant<-sig_correlation2(data2,padjustData2)
    #Get indexes for the rows and columns for the metabo data
    normMetabo_Tablesignificant$RowIndex<-apply(normMetabo_Tablesignificant,1,function(x){return(which(colnames(data2)[1:ncol(data2)]==x[1]))})
    normMetabo_Tablesignificant$ColIndex<-apply(normMetabo_Tablesignificant,1,function(x){return(which(colnames(data2)[1:ncol(data2)]==x[2]))})
    data.table::fwrite(normMetabo_Tablesignificant,paste(data2name,"_Tablesignificant.csv",sep=""))
    normMetabo_Tablesignificant<-data.table::fread(paste("./",data2name,"_Tablesignificant.csv",sep=""))
  }else{
    normMetabo_Tablesignificant<-data.table::fread(paste("./",data2name,"_Tablesignificant.csv",sep=""))
  }

  #Get low variance genes
  data1<-get_lowvarFeatures(data1) #This step was criticised
  #Remove those with too many zeros
  data1 <- data1[, which(as.numeric(colSums(data1 != 0)) > ceiling(nrow(data1)/2))] #At least two samples have a value higher than zero
  #data1<-as.data.frame(data1)

  #Get only those genes that correlate with the metabolites
  if(correlateDFs & !file.exists(paste("./",data1name,"_",data2name,".csv",sep=""))){
    print("Calculating correlations between data2 and data1")

    #Get Column indices of all metabolites
    metabo_indices<-unique(c(normMetabo_Tablesignificant$RowIndex,normMetabo_Tablesignificant$ColIndex))
    #Subset metabolites to correlate with genes
    SubSetdata2<-data2[,metabo_indices]

    ResultsCorDfs <- sig_correlation2Dfs(SubSetdata2,data1)
    genesSig <- unique(ResultsCorDfs$gene)

    #####################
    #Correct split number
    #missing
    #Add function missing
    ####################

    print(paste(length(genesSig),"features were kept from 'big' df",sep=" "))
    data1 <- data1[,genesSig]
    data.table::fwrite(data1,paste("./",data1name,"_",data2name,".csv",sep=""))
    data1<-data.table::fread(paste("./",data1name,"_",data2name,".csv",sep=""))
  }else if(file.exists("./data1.csv")){
    data1<-data.table::fread(paste("./",data1name,"_",data2name,".csv",sep=""))
  }

  if(ncol(data1)*nrow(normMetabo_Tablesignificant)>500){
    splitData1<-TRUE
  }

  #Split Data Frame
  if(splitData1){

    ls_dfs<-split_df(data1,old_split,split_number)
    print(paste("Data1 will be split into",length(ls_dfs),"parts",sep=" "))

    for (i in iteration_start:length(ls_dfs)){
      df_iter<-ls_dfs[[i]]

      #Wait 30 seconds to avoid overheating and errors in the parallelization...
      #If usign macs fan control, less pauses can be made...
      #if(i%%10==0){ #Every ten iterations, a pause will be made
      #  Sys.sleep(30)
      #}
      Sys.sleep(wait_iteration)


      #Convert to data.frames
      df_iter<-as.data.frame(df_iter)
      normMetabo_Tablesignificant<-as.data.frame(normMetabo_Tablesignificant)


      #Register parallel backend
      #library(doSNOW)
      if(is.null(numCores)){
        numCores<-detectCores()-2
        cat("Running parallelization with ",numCores," cores\n",sep="")
      }else{
        cat("Running parallelization with ",numCores," cores\n",sep="")
      }

      print(paste('Running CoNI Split Number',i,sep=" "))

      cl<-makeCluster(numCores)
      registerDoSNOW(cl)
      #registerDoParallel(numCores)

      df_results = foreach(j = 1:ncol(df_iter), .export =c("df_iter","normMetabo_Tablesignificant","data2","isNA_NAN"),.packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.inorder = FALSE) %dopar% { #Loop table significant metabolites %dopar%
        results2 =foreach(i = 1:nrow(normMetabo_Tablesignificant),.export =c("df_iter","normMetabo_Tablesignificant","data2","isNA_NAN"),.packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.verbose = FALSE,.inorder = FALSE) %dopar% { #Loop genes
          index1<-normMetabo_Tablesignificant[i,6]#Index column of first metabolite
          index2<-normMetabo_Tablesignificant[i,7]#Index column of second metabolite

          #Get metabolites names and gene name
          metabolite_name1<-normMetabo_Tablesignificant[i,1]
          metabolite_name2<-normMetabo_Tablesignificant[i,2]
          gene_name<-colnames(df_iter)[j]

          #Get correlation between metabolites
          cor_coefficient<-normMetabo_Tablesignificant[i,3]
          cor_pvalue<-normMetabo_Tablesignificant[i,4]

          #Calculate partial correlation between metabolites partialing out gene
          pcor_result<-pcor.test(data2[,index1],data2[,index2],df_iter[,j],method="p")
          pcor_pvalue<-pcor_result[[2]]
          pcor_coefficient<-pcor_result[[1]]

          #Sometimes the computer is not precise in float representation...
          if(pcor_coefficient > 1){
            pcor_coefficient<-0.999
          }else if(pcor_coefficient < -1){
            pcor_coefficient<- -0.999
          }

          #Correlation metabolites vs gene
          cor_m1_vs_g <- cor(data2[,index1],df_iter[,j])
          cor_m2_vs_g <- cor(data2[,index2],df_iter[,j])

          #Test if partialcorrelation coefficient differs from correlation coefficient
          cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(data2),
                                           alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

          cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(data2),
                                            alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

          cdgo_pvalue <- cdgo@steiger1980$p.value
          cdgo2_pvalue<- cdgo2@steiger1980$p.value

          rowtoprint<-list(metabolite_name1,metabolite_name2,gene_name, #change to list instead of cbind.data.frame
                           pcor_coefficient,pcor_pvalue,cor_coefficient,
                           cor_pvalue,cdgo_pvalue,cdgo2_pvalue)
        }
      }
      stopCluster(cl)

      #Get rid of rows with only NAs --> Not working...
      #df_results<-df_results[rowSums(is.na(df_results)) != ncol(df_results), ]
      df_results<-as.data.frame(df_results)
      colnames(df_results)<-c("metabolite_name1",	"metabolite_name2",	"gene_name",	"pcor_coefficient",	"pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue")
      df_results<-as.matrix(df_results)

      oldw <- getOption("warn")
      options(warn = -1)

      #Save result to memory
      writeT<-writeTable(df_results,num_cores = numCores,outputDir = outputDir,iteration = i) #Try to writ using fwrite
      if(!length(writeT)==0){write.csv(df_results,paste(outputDir,"CoNIOutputSplit",i,".csv",sep=""))}#If fwrite fails it is written with write.csv

      options(warn = oldw)

      #Remove results
      rm(df_results)

      #Print times
      iteration_endtime <- Sys.time()
      if(i==iteration_start){
        iteration_time<-difftime(iteration_endtime,start_time,units='mins')
        cat("Iteration time:",iteration_time,"minutes","\n",sep=" ")
        iteration_time_between<-iteration_endtime
      }else{
        iteration_time<-difftime(iteration_endtime,iteration_time_between,units='mins')
        cat("Iteration time:",iteration_time,"minutes","\n",sep=" ")
        iteration_time_between<-iteration_endtime
      }

    }
    #Merge output results CoNI
    CoNIOutput <- merge_outpuSplitFiles(outputDir)
    print('CoNI ran successfully')

    #Output processing time
    end_time <- Sys.time()
    total_time<-difftime(end_time,start_time,units='hours')
    cat(total_time,"hours", "\n",sep=" ")
    return(CoNIOutput)
  }else{
    print('Split is not necessary. Running standard CoNI...')
    splitData1<-FALSE
  }



  if (splitData1==FALSE){
    print('Running CoNI...')
    df_results = foreach(j = 1:ncol(data1), .combine=rbind, .inorder=FALSE) %dopar% {#Loop genes
      results2 = foreach(i = 1:nrow(normMetabo_Tablesignificant), .combine=rbind,.inorder=FALSE) %dopar% {#Loop table significant metabolites

        index1<-normMetabo_Tablesignificant[i,6]#Index column of first metabolite
        index2<-normMetabo_Tablesignificant[i,7]#Index column of second metabolite

        #Get metabolites names and gene name
        metabolite_name1<-normMetabo_Tablesignificant[i,1]
        metabolite_name2<-normMetabo_Tablesignificant[i,2]
        gene_name<-colnames(data1)[j]

        #Get correlation between metabolites
        cor_coefficient<-normMetabo_Tablesignificant[i,3]
        cor_pvalue<-normMetabo_Tablesignificant[i,4]

        #Calculate partial correlation between metabolites partialing out gene
        pcor_result<-pcor.test(data2[,index1],data2[,index2],data1[,j],method="p")
        pcor_pvalue<-pcor_result[[2]]
        pcor_coefficient<-pcor_result[[1]]

        if(pcor_coefficient > 1){
          pcor_coefficient<-0.999
        }else if(pcor_coefficient < -1){
          pcor_coefficient<- -0.999
        }

        #Correlation metabolites vs gene
        cor_m1_vs_g <- cor(data2[,index1],data1[,j])
        cor_m2_vs_g <- cor(data2[,index2],data1[,j])

        #Test if partialcorrelation coefficient differs from correlation coefficient
        cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(data2),
                                         alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(data2),
                                          alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo_pvalue <- cdgo@steiger1980$p.value
        cdgo2_pvalue<- cdgo2@steiger1980$p.value
        rowtoprint<-cbind.data.frame(metabolite_name1,metabolite_name2,gene_name,
                                     pcor_coefficient,pcor_pvalue,cor_coefficient,
                                     cor_pvalue,cdgo_pvalue,cdgo2_pvalue)

      }
    }
    #Output processing time
    print('CoNI ran successfully')
    end_time <- Sys.time()
    print(end_time - start_time)
    return(df_results)
  }

}


#This function tries to write a table with fread
writeTable <- function(results_write,num_cores,outputDir,iteration) {
  out <- tryCatch(
    {
      data.table::fwrite(results_write, paste(outputDir,"CoNIOutputSplit",iteration,".csv",sep=""),nThread=num_cores)
    },
    error=function(cond) {
      message('fwrite failed')
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return('NA')
    }
  )
  return(out)
}



#This function tests if the result of a calculation is NA or NAN
isNA_NAN<-function(result){
  if(is.na(result)|is.nan(result)){
    an<- TRUE
  }else{
    an <- FALSE
  }
  an
}

#This function tests if the input files exist, if they do not it will output an error and
#end CoNI
do_objectsExist<-function(gene_exp,norm_metabo_dat){
  if(missing(gene_exp) | missing(norm_metabo_dat)){
    message("Input objects are missing")
    stop("CoNI end")
  }else{
    print("Input objects are defined")
  }
}


#This function wills split the 'big' omics data into smaller data frames to improve computations, avoid run out of memory
#The second option is probably better... need to review..
#In the second otion the data is splitted according to the number given + 1 that are the remaining of the division (modulo)
split_df<-function(AbundantDF,numberSplitDF=5,numberSplitDF_2=NULL){
  dt_list<-list()
  if(is.null(numberSplitDF_2)){
    if(ncol(AbundantDF) %% numberSplitDF !=0){
      SplitFirstn<-floor(ncol(AbundantDF)/numberSplitDF)
      remainingDF<-floor(ncol(AbundantDF)/numberSplitDF) + ncol(AbundantDF) %% numberSplitDF
      for (i in 1:numberSplitDF){
        if (i==1){
          dt_list[[i]]<-AbundantDF[,1:SplitFirstn]
        }else if (i<numberSplitDF){
          start<-SplitFirstn*(i-1)+1
          end<-SplitFirstn*i
          dt_list[[i]]<-AbundantDF[,start:end]
        }else{
          start<-SplitFirstn*(numberSplitDF-1)+1
          end<-length(AbundantDF)
          dt_list[[i]]<-AbundantDF[,start:end]
        }
      }
    }else{
      Splitn<-ncol(AbundantDF)/numberSplitDF
      for (i in 1:numberSplitDF){
        if (i==1){
          dt_list[[i]]<-AbundantDF[,1:Splitn]
        }else{
          start<-Splitn*(i-1)+1
          end<-Splitn*i
          dt_list[[i]]<-AbundantDF[,start:end]
        }
      }
    }
  }else{
    if(ncol(AbundantDF) %% numberSplitDF_2 !=0){
      SplitFirstParts<-floor(ncol(AbundantDF)/numberSplitDF_2)
      start<-1 #Start to slice is position 1 of data frame
      end<-SplitFirstParts #
      i=1
      while(end <= SplitFirstParts*numberSplitDF_2){
        dt_list[[i]]<-AbundantDF[,start:end]
        start<-start+SplitFirstParts
        end<-start+SplitFirstParts-1
        i<-i+1
      }
      start<-SplitFirstParts*numberSplitDF_2+1
      dt_list[[i]]<-AbundantDF[,start:ncol(AbundantDF)]
    }else{
      split_size<-ncol(AbundantDF)/numberSplitDF_2
      start<-1 #Start to slice is position 1 of data frame
      end<-split_size
      i=1
      while(end <= ncol(AbundantDF)){
        dt_list[[i]]<-AbundantDF[,start:end]
        start<-start+split_size
        end<-start+split_size-1
        i<-i+1
      }
    }

  }
  dt_list
}

#This function compares the sample names between the two omics provided
#If they do not match CoNI fails to run
compare_sampleNames<-function(df1,df2){
  Rowsdf1<-rownames(df1)[order(rownames(df1))]
  Rowsdf2<-rownames(df2)[order(rownames(df2))]
  if(!identical(Rowsdf1,Rowsdf2)){
    print('Sample names between datasets do not match')
    print('Make sure omics data comes from the same samples and that sample names are consistent across datasets')
    stop('CoNI end')
  }else{
    print('Samples match')
  }
}


#This function gets rid of columns and/or rows with at least one NA
get_ridNA<-function(df, type='column'){
  `%notin%` <- Negate(`%in%`)
  if(type %notin% c('column','row','both') || !is.data.frame(df)){
    print('Choose column, row or both, and provide a data frame')
    stop('Exit function')
  }

  if(type=='both'){
    df<-df[rowSums(is.na(df)) == 0 , colSums(is.na(df)) == 0]
  }else if(type=='column'){
    df<-df[ , colSums(is.na(df)) == 0]
  }else{
    df<-df[rowSums(is.na(df)) == 0 , ]
  }
  df
}

##This function gets the upper part of the matrix after calculating the
#correlation coefficients between all pairs of elements
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]], #Get rownames
    column = rownames(cormat)[col(cormat)[ut]], #Get colnames
    cor  =(cormat)[ut], #Get correlation coefficients
    p = pmat[ut] #Get p values
  )
}

#NewFunction more similar to CONI
#Get significant correlations
sig_correlation2<-function(input_data1,padj=TRUE){
  corr<-Hmisc::rcorr(as.matrix(input_data1),type='p')
  corr_table<-flattenCorrMatrix(corr$r,corr$P)
  corr_table$adj.p<-p.adjust(corr_table$p)

  if(padj){
    corr_tableSig <- corr_table %>% filter(adj.p<0.05)
    if(nrow(corr_tableSig) == 0){
      print('No features significantly correlate after padjustment for the reduced df')
      print('Using non adjusted pvalues')
      corr_tableSig<-corr_table %>% filter(p<0.05)}
  }else{
    print("Ajustment for multiple testing was set to FALSE for 'small' df")
    corr_tableSig<-corr_table %>% filter(p<0.05)
  }
  print(paste('Significant correlations',nrow(corr_tableSig),sep=" "))
  corr_tableSig
}

#This functions input are two data frames (e.g. metabolites and genes)
#it calculates the correlation matrix and creates a table with only significant pairs
#No correction for multiple testing is done
sig_correlation2Dfs<-function(metabolite_data,gene_expression){
  n <- t(!is.na(metabolite_data)) %*% (!is.na(gene_expression)) # same as count.pairwise(x,y) from psych package/ Matches number of samples
  r <- cor(metabolite_data, gene_expression, use = "pairwise.complete.obs") # MUCH MUCH faster than corr.test()
  cor2pvalue = function(r, n) {
    t <- (r*sqrt(n-2))/sqrt(1-r^2)
    p <- 2*(1 - pt(abs(t),(n-2)))
    se <- sqrt((1-r*r)/(n-2))
    out <- list(r, n, t, p, se)
    names(out) <- c("r", "n", "t", "p", "se")
    return(out)
  }

  # Get a list with matrices of correlation, pvalues, standard error, etc.
  result = cor2pvalue(r,n)
  rcoeffMatrix<-result$r
  pvalueMatrix<-result$p

  rows<-rownames(rcoeffMatrix)
  cols<-colnames(pvalueMatrix)

  df <- data.frame(metabolite=character(),
                   gene=character(),
                   cor=double(),
                   pvalue=double(),
                   stringsAsFactors=FALSE)

  for(i in rows){
    for(j in cols){
      if (pvalueMatrix[i,j]>0.05){
        next
      }else{
        cor<-rcoeffMatrix[i,j]
        pvalue<-pvalueMatrix[i,j]
        df<- df %>% add_row(metabolite = i, gene = j, cor = cor, pvalue = pvalue)
      }

    }
  }
  df
}




#Get low variance features
get_lowvarFeatures<-function(df){
  df<-as.data.frame(t(df))
  df$Var<-genefilter::rowVars(as.matrix(df))
  df<-subset(df,Var<0.5)
  df<-df[,-ncol(df)] #get rid of the column with the variances
  df<-as.data.frame(t(df))
}


#This function reads the output files generated by CoNIParallel
#The folder should not contain any other files except CoNI output
merge_outpuSplitFiles<-function(outputDir){
  #outputDir<-gsub('\\.','',outputDir)
  #outputDir<-gsub('\\/','',outputDir)
  file_list <- list.files(outputDir)
  file_list<-file_list[grep("CoNI",file_list)]
  for (file in file_list){
    # if the merged dataset doesn't exist, create it
    if (!exists("datasetResultsCoNI")){
      datasetResultsCoNI <- data.table::fread(paste(outputDir,file,sep=""), header=TRUE, sep=",")
    }else{
      temp_dataset <-data.table::fread(paste(outputDir,file,sep=""), header=TRUE, sep=",")
      datasetResultsCoNI<-rbind(datasetResultsCoNI, temp_dataset)
      rm(temp_dataset)
    }
  }
  datasetResultsCoNI
}


#Check if output directory exists and if it does not it will create it
check_outputDir<-function(outputDir){
  if (file.exists(paste(outputDir,sep=""))) {
    print("Output directory exists")
  } else {
    print("Output directory does not exist - creating directory ... ")
    dir.create(file.path(paste(outputDir,sep="")))
  }
}

#Add if gene is a low expressed gene
#Expression table contains in one column the genes and a second column that specifies if the gene is a low expressed gene = LowExpression
#or not = ""
expression_level<-function(CoNI_results, expressionTable){
  expressionLevel<-expressionLevelTable[match(CoNI_results$gene_name,expressionTable$Gene),"ExpressionLevel"]
  expressionLevel
}


#Create network using as input the output of CoNI
generate_network<-function(ResultsCoNI, colorNodesTable){
  results_SteigerAdjust <- ResultsCoNI[,1:3] #Get pair metabolites and gene

  #Summarize results for network construction
  df <- plyr::ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
              Genes=length(gene_name),
              GenesString=paste0(unique(gene_name),collapse=";"))
  colnames(df) <- c("from","to","weight","Genes")
  clinksd <- df
  clinksd$type <- "hyperlink"
  clinksd$width <- clinksd$weight/max(clinksd$weight) #Calculate a width based on the maximum number of genes per connection
  cnodes <- data.frame("Name"=unique(c(as.character(df$from),as.character(df$to))),stringsAsFactors=F)#Get the nodes (metabolites)
  #Assign colors to nodes
  m <- merge(cnodes,colorNodesTable,by.x="Name",by.y=colnames(colorNodesTable)[1],all.x=T)

  cnodesd <- m
  #Change column names
  colnames(clinksd)[3] <- "weightreal"
  colnames(clinksd)[6] <- "weight"

  #Create graph
  netd <- igraph::graph_from_data_frame(d=clinksd, vertices=cnodesd, directed=F)
  netdhfd <- igraph::simplify(netd,remove.multiple=F)

  netdhfd
}

get_variableName <- function(variable) {
  deparse(substitute(variable))
}


generate_network_2<-function(ResultsCoNI, colorNodesTable,outputDir="./",outputFileName="ResultsCoNI"){
  results_SteigerAdjust <- ResultsCoNI[,c(1:7)] #Get pair metabolites, gene and pcor and cor information... change to add more information
  #Summarize results for network construction
  df<-plyr::ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
            weightreal=length(gene_name),
            Genes=paste0(unique(gene_name),collapse=";"),
            #ActualGeneNames=paste0(unique(ActualGeneName),collapse=";"),
            PcorValues=paste0(pcor_coefficient,collapse=";"),
            CorValues=paste0(cor_coefficient,collapse=";"),
            PcorAverage=mean(pcor_coefficient),
            CorAverage=mean(cor_coefficient),
            DirectionPcor=ifelse(sum(pcor_coefficient<0)>sum(pcor_coefficient>0),'Negative',ifelse(sum(pcor_coefficient<0) == sum(pcor_coefficient>0),'Balanced','Positive')))
  colnames(df)[1:2] <- c("from","to")
  clinksd <- df
  clinksd$type <- "hyperlink"
  clinksd$weight <- clinksd$weightreal/max(clinksd$weightreal) #Calculate a width based on the maximum number of genes per connection
  #Save table
  write.csv(clinksd,paste(outputDir,"TableForNetwork_",outputFileName,".csv",sep=""))

  cnodes <- data.frame("Name"=unique(c(as.character(df$from),as.character(df$to))),stringsAsFactors=F)#Get the nodes (metabolites)
  #Assign colors to nodes
  m <- merge(cnodes,colorNodesTable,by.x="Name",by.y=colnames(colorNodesTable)[1],all.x=T)
  cnodesd <- m
  #Change column names
  #colnames(clinksd)[10] <- "weight"

  #Create graph
  netd <- igraph::graph_from_data_frame(d=clinksd, vertices=cnodesd, directed=F)
  netdhfd <- igraph::simplify(netd,remove.multiple=F)

  netdhfd
}


#Find local regulated genes
find_localRegulatedFeatures<-function(ResultsCoNI,network){
  ls2 <- length(unique(ResultsCoNI$gene_name)) #get number of genes affecting metabolites
  #Distance = 2 -> Second level neighborhood?
  df <- list()

  for(i in names(igraph::V(network))){ #loop nodes of graph
    l <- igraph::V(network)$name[neighbors(network, i)] #Get first level neighbors of node in iteration
    l1 <- list()
    for(j in l){
      l1[[j]] <- igraph::V(network)$name[neighbors(network, j)]#Loop neighbors of the node in the iteration and get their neighbors (Second level neighborhood)
    }
    l1 <- unique(unlist(l1)) #Get unique 2nd level neighbors
    s <- subset(ResultsCoNI, ((metabolite_name1==i & metabolite_name2 %in% l) | (metabolite_name2==i & metabolite_name1 %in% l)) |
                  ((metabolite_name1 %in% l & metabolite_name2 %in% l1) | (metabolite_name2 %in% l & metabolite_name1 %in% l1)))
    a <- length(unique(paste0(s$metabolite_name1,"_",s$metabolite_name2)))
    d <- length(unique(s$gene_name))
    e <- nrow(s)
    s <- droplevels(s)
    b <- table(s$gene_name)
    df[[i]] <- data.frame("Node1"=rep(i,length(b)),"Edges"=rep(a,length(b)),"Draws"=rep(e,length(b)),"GenesInTotal"=rep(d,length(b)),as.data.frame(b))
  }
  #Generate result data frame
  res2 <- do.call(rbind.data.frame, df)

  #we use the binomial distribution to test if the enrichment is significant as we can draw a gene for an area more often
  res2$Pval <- apply(res2,1,function(x){dbinom(as.numeric(x[[6]]),as.numeric(x[[3]]),1/ls2)})
  res2$Padj <- p.adjust(res2$Pval)
  res2 <- res2[order(res2$Padj),]
  res2Sig <- subset(res2,res2$Padj<0.05)

  res2Sig
  #Here we also use the Pvalue as cuotff for this tiny example
  #res2Sig <- subset(res2,res2$Pval<0.05)
}

tableLRGs_Metabolites<-function(CoNIResults,LRGenes){
  CoNIResults_LRGs<-CoNIResults[CoNIResults$gene_name %in% LRGenes,]
  Gene_TableLRGs<- plyr::ddply(CoNIResults_LRGs, .(metabolite_name1,metabolite_name2), plyr::summarize,
                         Genes=paste(gene_name,collapse=","))
  #Join Metabolite pairs
  CoNIResults_LRGs_MetaboliteJoined<-unite(CoNIResults_LRGs,MetabolitePair,metabolite_name1,metabolite_name2,sep="-")
  CoNIResults_LRGs_MetaboliteJoined<-CoNIResults_LRGs_MetaboliteJoined[,c(1,2)]

  #Chowate table Genes and their corresponding Metabolite pairs
  LRGs_and_MPairs <- plyr::ddply(CoNIResults_LRGs_MetaboliteJoined, .(gene_name), plyr::summarize,
                           MetabolitePairs=paste(MetabolitePair,collapse=","))
  #Temporary table
  temp<-as.data.frame(CoNIResults_LRGs[,c(1:3)])

  #Add to the LRGs and Metabolites pairs the unique individual metabolites
  LRGs_and_MPairs$Metabolites<-ddply(temp, .(gene_name), plyr::summarize,
                                     Metabolites=paste(unique(c(as.character(metabolite_name1),as.character(metabolite_name2))),collapse=","))[,2]
  colnames(LRGs_and_MPairs)<-c("Local Regulated Gene","Metabolite pairs","Metabolites")
  LRGs_and_MPairs
}
