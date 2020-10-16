#CoNI Functions
############

####################################
#Functions
####################################
CoNI<- function(driverD, linkedD,outputName="CoNIOutput",driverDname="driverD",linkedDname="linkedD",padjustlinkedD=TRUE, correlateDFs=TRUE,splitdriverD=TRUE,split_number=2,outputDir="./CoNIOutput/",delPrevious=FALSE,delIntermediaryFiles=TRUE,iteration_start=1,numCores=NULL,verbose=TRUE,inhouseSteigerPcor=FALSE) {

  #Check if input objects are defined
  do_objectsExist(driverD,linkedD,verbose)

  #Start measuring time
  start_time <- Sys.time()

  #Check if output directory exists
  check_outputDir(outputDir,verbose)

  #Check if previous files are present and delete accordingly
  check_previous(delPrevious,iteration=iteration_start,outDir=outputDir,verb=verbose)

  #Test if sample names are the same in both data sets
  compare_sampleNames(driverD,linkedD)

  #Make sure column names are appropiate
  colnames(driverD)<-make.names(colnames(driverD),unique=TRUE)
  colnames(linkedD)<-make.names(colnames(linkedD),unique=TRUE)

  if(!file.exists(paste(outputDir,"KeptFeatures_",linkedDname,".csv",sep=""))){
    #Get significant correlations between metabolites
    if(verbose){print("Calculating correlations of linked Data")}
    normMetabo_Tablesignificant<-sig_correlation2(input_driverD = linkedD,padj = padjustlinkedD,verb = verbose)
    #Get indexes for the rows and columns for the metabo data
    normMetabo_Tablesignificant$RowIndex<-apply(normMetabo_Tablesignificant,1,function(x){return(which(colnames(linkedD)[1:ncol(linkedD)]==x[1]))})
    normMetabo_Tablesignificant$ColIndex<-apply(normMetabo_Tablesignificant,1,function(x){return(which(colnames(linkedD)[1:ncol(linkedD)]==x[2]))})
    data.table::fwrite(normMetabo_Tablesignificant,paste(outputDir,"KeptFeatures_",linkedDname,".csv",sep=""))
    normMetabo_Tablesignificant<-data.table::fread(paste(outputDir,"KeptFeatures_",linkedDname,".csv",sep=""))
  }else{
    normMetabo_Tablesignificant<-data.table::fread(paste(outputDir,"KeptFeatures_",linkedDname,".csv",sep=""))
  }

  #Get low variance genes
  driverD<-get_lowvarFeatures(driverD) #This step was criticised
  #Remove those with too many zeros
  driverD <- driverD[, which(as.numeric(colSums(driverD != 0)) > ceiling(nrow(driverD)/2))] #At least two samples have a value higher than zero
  #driverD<-as.data.frame(driverD)

  #Get only those genes that correlate with the metabolites
  if(correlateDFs & !file.exists(paste(outputDir,"KeptFeatures_",driverDname,".csv",sep=""))){
    if(verbose){print("Calculating correlations between linked Data and driver Data")}

    #Get Column indices of all metabolites
    metabo_indices<-unique(c(normMetabo_Tablesignificant$RowIndex,normMetabo_Tablesignificant$ColIndex))
    #Subset metabolites to correlate with genes
    SubSetlinkedD<-linkedD[,metabo_indices]

    ResultsCorDfs <- sig_correlation2Dfs(SubSetlinkedD,driverD)
    genesSig <- unique(ResultsCorDfs$gene)
    if(verbose){print(paste(length(genesSig),"features were kept from driver Data",sep=" "))}
    driverD <- driverD[,genesSig]
    data.table::fwrite(driverD,paste(outputDir,"KeptFeatures_",driverDname,".csv",sep=""))
    driverD<-data.table::fread(paste(outputDir,"KeptFeatures_",driverDname,".csv",sep=""))
  }else if(file.exists(paste(outputDir,"KeptFeatures_",driverDname,".csv",sep=""))){
    driverD<-data.table::fread(paste(outputDir,"KeptFeatures_",driverDname,".csv",sep=""))
  }

  if(ncol(driverD)*nrow(normMetabo_Tablesignificant)>10000){
    if(splitdriverD==FALSE){
      cat('For computational purposes a split will be performed\n')
      splitdriverD<-TRUE
    # }else if(ncol(driverD)/10 > split_number){
    #   cat('Provided split number is too small, split number was adjusted')
    #   split_number<-ncol(driverD)/10
     }
  }

  #Split Data Frame
  if(splitdriverD){

    ls_dfs<-split_df(driverD,split_number)
    print(paste("Driver Data was split into",length(ls_dfs),"parts",sep=" "))

    for (i in iteration_start:length(ls_dfs)){
      df_iter<-ls_dfs[[i]]

      # #Wait 30 seconds to avoid overheating and errors in the parallelization...
      # #If usign macs fan control, less pauses can be made...
      # #if(i%%10==0){ #Every ten iterations, a pause will be made
      # #  Sys.sleep(30)
      # #}
      # Sys.sleep(wait_iteration)


      #Convert to data.frames
      df_iter<-as.data.frame(df_iter)
      normMetabo_Tablesignificant<-as.data.frame(normMetabo_Tablesignificant)


      #Register parallel backend
      #library(doSNOW)
      if(is.null(numCores)){
        numCores<-detectCores()-2
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }else{
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }

      if(verbose){print(paste('Running CoNI Split Number',i,sep=" "))}

      cl<-makeCluster(numCores)
      registerDoSNOW(cl)
      #registerDoParallel(numCores)

      pb<-tkProgressBar(max=ncol(driverD))
      df_results = foreach(j = 1:ncol(df_iter), .packages = c("ppcor", "doParallel","cocor","progress") , .combine=rbind,.inorder = FALSE,.options.snow=taskBar1(pb,ncol(driverD))) %dopar% { #Loop table significant metabolites %dopar%
        results2 =foreach(i = 1:nrow(normMetabo_Tablesignificant),.packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.verbose = FALSE,.inorder = FALSE) %dopar% { #Loop genes
          index1<-normMetabo_Tablesignificant[i,6]#Index column of first metabolite
          index2<-normMetabo_Tablesignificant[i,7]#Index column of second metabolite

          #Get metabolites names and gene name
          Feature_1_linkedD<-normMetabo_Tablesignificant[i,1]
          Feature_2_linkedD<-normMetabo_Tablesignificant[i,2]
          Feature_driverD<-colnames(df_iter)[j]

          #Get correlation between metabolites
          cor_coefficient<-normMetabo_Tablesignificant[i,3]
          cor_pvalue<-normMetabo_Tablesignificant[i,4]

          #############################
          #Calculate partial correlation between metabolites partialing out gene
          #j=metabolite1,k=metabolite2,h=gene
          #r.jk_h
          if(inhouseSteigerPcor){#

            #Part of original output
            pcor_result<-pcor.test(linkedD[,index1],linkedD[,index2],df_iter[,j],method="p")
            pcor_pvalue<-pcor_result[[2]]
            pcor_coefficient<-pcor_result[[1]]

            #########################################
            #New Test using partial correlation values
            #r.kh_j or r.jh_k
            #Metabolite 2 and Gene | Metabolite 1
            pcor_res_kh_j<-pcor.test(linkedD[,index2],df_iter[,j],linkedD[,index2],method="p")
            pcor_res_kh_jCoef<-pcor_res_kh_j[[1]]
            pcor_res_kh_jpval<-pcor_res_kh_j[[2]]

            #Metabolite 1 and Gene | Metabolite 2
            pcor_res_jh_k<-pcor.test(linkedD[,index1],df_iter[,j],linkedD[,index2],method="p")
            pcor_res_jh_kCoef<-pcor_res_jh_k[[1]]
            pcor_res_jh_kpval<-pcor_res_jh_k[[2]]

            #Metabolite 1 and Metabolite 2 | Gene
            pcor_res_jk_h<-pcor.test(linkedD[,index1],linkedD[,index2],df_iter[,j],method="p")
            pcor_res_jk_hCoef<-pcor_res_jk_h[[1]]
            pcor_res_jk_hpval<-pcor_res_jk_h[[2]]


            TSteigerPcorInput<- cocor.dep.groups.overlap(r.jk=pcor_res_jh_kCoef, r.jh=pcor_res_kh_jCoef, r.kh=pcor_res_jk_hCoef, n=nrow(linkedD),
                                             alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

            TSteigerPcorInput_pvalue<-TSteigerPcorInput@steiger1980$p.value

            direction<-ifelse(cor_coefficient==pcor_res_jk_hCoef,"equal",ifelse(cor_coefficient<pcor_res_jk_hCoef,"increase","decrease"))

            #Original output
            pcor_result<-pcor.test(linkedD[,index1],linkedD[,index2],df_iter[,j],method="p")
            pcor_pvalue<-pcor_result[[2]]
            pcor_coefficient<-pcor_result[[1]]

            #Sometimes the computer is not precise in float representation...
            if(pcor_coefficient > 1){
              pcor_coefficient<-0.999
            }else if(pcor_coefficient < -1){
              pcor_coefficient<- -0.999
            }

            #Correlation metabolites vs gene
            cor_m1_vs_g <- cor(linkedD[,index1],df_iter[,j])
            cor_m2_vs_g <- cor(linkedD[,index2],df_iter[,j])

            #Test if partialcorrelation coefficient differs from correlation coefficient
            cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(linkedD),
                                             alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

            cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(linkedD),
                                              alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

            cdgo_pvalue <- cdgo@steiger1980$p.value
            cdgo2_pvalue<- cdgo2@steiger1980$p.value

            #############
            #D suggestion
            #New Test
            SteigerM1M2_M1G<- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=cor_m1_vs_g, r.kh=cor_m2_vs_g, n=nrow(linkedD),
                                                       alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
            SteigerM1M2_M1G_pvalue<-SteigerM1M2_M1G@steiger1980$p.value

            SteigerM1M2_M2G<- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=cor_m2_vs_g, r.kh=cor_m1_vs_g, n=nrow(linkedD),
                                                       alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
            SteigerM1M2_M2G_pvalue<-SteigerM1M2_M2G@steiger1980$p.value




            rowtoprint<-list(Feature_1_linkedD,Feature_2_linkedD,Feature_driverD, #change to list instead of cbind.data.frame
                             pcor_coefficient,pcor_pvalue,cor_coefficient,
                             cor_pvalue,cdgo_pvalue,cdgo2_pvalue,pcor_res_kh_jCoef,pcor_res_jh_kCoef,cor_coefficient,pcor_res_jk_hCoef,direction,TSteigerPcorInput_pvalue,SteigerM1M2_M1G_pvalue,SteigerM1M2_M2G_pvalue)

            #New test
            # cdgo <- cocor.dep.groups.overlap_inHouse(r.jk=cor_coefficient[[1]], r.jk_h = pcor_coefficient,r.kh_j = pcor_res_kh_jCoef,
            #                                          n = nrow(linkedD),alternative="two.sided", alpha=0.05, conf.level=0.95)
            # cdgo_pvalue<-cdgo$p.value
            # #
            # # cdgo2 <- cocor.dep.groups.overlap_inHouse(r.jk=cor_coefficient[[1]], r.jk_h = pcor_coefficient,r.kh_j = pcor_res_jh_kCoef,
            # #                                          n = nrow(linkedD),alternative="two.sided", alpha=0.05, conf.level=0.95)
            # # cdgo2_pvalue<-cdgo2$p.value
            #
            # rowtoprint<-list(Feature_1_linkedD,Feature_2_linkedD,Feature_driverD, #change to list instead of cbind.data.frame
            #                  pcor_coefficient,pcor_pvalue,cor_coefficient,
            #                  cor_pvalue,cdgo_pvalue,cdgo2_pvalue)
          }else{
            #Calculate partial correlation between metabolites partialing out gene
            pcor_result<-pcor.test(linkedD[,index1],linkedD[,index2],df_iter[,j],method="p")
            pcor_pvalue<-pcor_result[[2]]
            pcor_coefficient<-pcor_result[[1]]

            #Sometimes the computer is not precise in float representation...
            if(pcor_coefficient > 1){
              pcor_coefficient<-0.999
            }else if(pcor_coefficient < -1){
              pcor_coefficient<- -0.999
            }

            #Correlation metabolites vs gene
            cor_m1_vs_g <- cor(linkedD[,index1],df_iter[,j])
            cor_m2_vs_g <- cor(linkedD[,index2],df_iter[,j])

            #Test if partialcorrelation coefficient differs from correlation coefficient
            cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(linkedD),
                                             alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

            cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(linkedD),
                                              alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

            cdgo_pvalue <- cdgo@steiger1980$p.value
            cdgo2_pvalue<- cdgo2@steiger1980$p.value

            rowtoprint<-list(Feature_1_linkedD,Feature_2_linkedD,Feature_driverD, #change to list instead of cbind.data.frame
                             pcor_coefficient,pcor_pvalue,cor_coefficient,
                             cor_pvalue,cdgo_pvalue,cdgo2_pvalue)

          }
          #############################

        }

      }
      close(pb)
      stopCluster(cl)

      #Get rid of rows with only NAs --> Not working...
      #df_results<-df_results[rowSums(is.na(df_results)) != ncol(df_results), ]
      df_results<-as.data.frame(df_results)
      #Note cor_coefficient is the same as Cor_M1M2
      if(inhouseSteigerPcor){
        colnames(df_results)<-c("Feature_1_linkedD",	"Feature_2_linkedD",	"Feature_driverD",	"pcor_coefficient",	"pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue","Pcor_M2G_M1","Pcor_M1G_M2","Cor_M1M2","Pcor_M1M2_G","Direction","NewTESTPcor","SteigerM1M2_M1G","SteigerM1M2_M2G")
      }else{
        colnames(df_results)<-c("Feature_1_linkedD",	"Feature_2_linkedD",	"Feature_driverD",	"pcor_coefficient",	"pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue")
      }
      df_results<-as.matrix(df_results)

      #oldw <- getOption("warn")
      #options(warn = -1)

      #Save result to memory
      writeT<-writeTable(df_results,num_cores = numCores,outputDir = outputDir,iteration = i) #Try to writ using fwrite
      if(!length(writeT)==0){write.csv(df_results,paste(outputDir,"CoNIOutputSplit",i,".csv",sep=""))}#If fwrite fails it is written with write.csv

      #options(warn = oldw)

      #Remove results
      rm(df_results)

      #Print times
      iteration_endtime <- Sys.time()
      if(i==iteration_start){
        iteration_time<-difftime(iteration_endtime,start_time,units='mins')
        if(verbose){cat("Iteration time:",iteration_time,"minutes","\n",sep=" ")}
        iteration_time_between<-iteration_endtime
      }else{
        iteration_time<-difftime(iteration_endtime,iteration_time_between,units='mins')
        if(verbose){cat("Iteration time:",iteration_time,"minutes","\n",sep=" ")}
        iteration_time_between<-iteration_endtime
      }

    }
    #Merge output results CoNI
    CoNIOutput <- merge_outpuSplitFiles(outputDir)
    #Change column names
    #colnames(CoNIOutput)<-c(paste0("Feature_1_",linkedDname),	paste0("Feature_2_",linkedDname),	paste0("Feature_",driverDname),	"pcor_coefficient",	"pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue")

    #################
    #This step might be problematic with very big data... I need to add an alternative way
    #To do... add alternative method for adjusted steiger
    #Add adjusted steiger pvalue
    #CoNIOutput$cdgo_adjusted<-p.adjust(CoNIOutput$cdgo_pvalue)
    #CoNIOutput$cdgo2_adjusted<-p.adjust(CoNIOutput$cdgo2_pvalue)

    #Save results
    suppressMessages(data.table::fwrite(CoNIOutput, paste(outputDir,outputName,".csv",sep=""),nThread=numCores))
    #Delete intermediary files
    delIntFiles(delIntermediaryFiles,outputDir)

    #Output processing time
    end_time <- Sys.time()
    total_time<-difftime(end_time,start_time,units='hours')
    cat(total_time,"hours", "\n",sep=" ")
    print('CoNI ran successfully')
    return(CoNIOutput)
  }else{
    print('Split was set to FALSE')
    splitdriverD<-FALSE
  }

  if (splitdriverD==FALSE){
    #Register parallel backend
    #library(doSNOW)
    if(is.null(numCores)){
      numCores<-detectCores()-2
      if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
    }else{
      if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
    }

    cl<-makeCluster(numCores)
    registerDoSNOW(cl)

    print('Running CoNI...')

    df_results = foreach(j = 1:ncol(driverD), .combine=rbind,.packages = c("ppcor", "doParallel","cocor"), .inorder=FALSE) %dopar% {#Loop genes
      results2 = foreach(i = 1:nrow(normMetabo_Tablesignificant), .combine=rbind,.packages = c("ppcor", "doParallel","cocor") ,.inorder=FALSE) %dopar% {#Loop table significant metabolites

        index1<-normMetabo_Tablesignificant[i,6]#Index column of first metabolite
        index2<-normMetabo_Tablesignificant[i,7]#Index column of second metabolite

        #Get metabolites names and gene name
        Feature_1_linkedD<-normMetabo_Tablesignificant[i,1]
        Feature_2_linkedD<-normMetabo_Tablesignificant[i,2]
        Feature_driverD<-colnames(driverD)[j]

        #Get correlation between metabolites
        cor_coefficient<-normMetabo_Tablesignificant[i,3]
        cor_pvalue<-normMetabo_Tablesignificant[i,4]

        #Calculate partial correlation between metabolites partialing out gene
        pcor_result<-pcor.test(linkedD[,index1],linkedD[,index2],driverD[,j],method="p")
        pcor_pvalue<-pcor_result[[2]]
        pcor_coefficient<-pcor_result[[1]]

        if(pcor_coefficient > 1){
          pcor_coefficient<-0.999
        }else if(pcor_coefficient < -1){
          pcor_coefficient<- -0.999
        }

        #Correlation metabolites vs gene
        cor_m1_vs_g <- cor(linkedD[,index1],driverD[,j])
        cor_m2_vs_g <- cor(linkedD[,index2],driverD[,j])

        #Test if partialcorrelation coefficient differs from correlation coefficient
        cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(linkedD),
                                         alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(linkedD),
                                          alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo_pvalue <- cdgo@steiger1980$p.value
        cdgo2_pvalue<- cdgo2@steiger1980$p.value
        rowtoprint<-cbind.data.frame(Feature_1_linkedD,Feature_2_linkedD,Feature_driverD,
                                     pcor_coefficient,pcor_pvalue,cor_coefficient,
                                     cor_pvalue,cdgo_pvalue,cdgo2_pvalue)

      }
    }
    #colnames(df_results)<-c(paste0("Feature_1_",linkedDname),	paste0("Feature_2_",linkedDname),	paste0("Feature_",driverDname),	"pcor_coefficient",	"pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue")
    #Save results
    suppressMessages(data.table::fwrite(df_results, paste(outputDir,outputName,".csv",sep=""),nThread=numCores))

    #Output processing time
    end_time <- Sys.time()
    print(end_time - start_time)
    print('CoNI ran successfully')
    return(df_results)
  }

}


taskBar1<-function(pb,ntasks){
  progress<-function(n) setTkProgressBar(pb,n, label=paste(round(n/ntasks*100,0), "%"))
  opts<-list(progress=progress)
  return(opts)
}

#################
#Ratios approach

#Calculate log2FoldChanges
calculateLog2FoldChange<-function(matx){
  matx<-as.matrix(matx)
  combinations<-combn(ncol(matx),2) #Get all possible combination
  mat <- matrix(, nrow = nrow(matx), ncol = ncol(combinations)) #create empty matrix to assign all possible combinations
  namesResults<-c()
  for (col in 1:ncol(combinations)){
    index1<-combinations[1,col]
    index2<-combinations[2,col]

    #Names
    met1_name<-colnames(matx)[index1]
    met2_name<-colnames(matx)[index2]
    nameComb<-paste(met1_name,"__",met2_name,sep = "")
    namesResults<-c(namesResults,nameComb)

    #Log2FoldChange
    foldChange<- matx[,index1]/matx[,index2]
    log2FC<- log(foldChange,2)
    mat[,col]<-log2FC
  }
  colnames(mat)<-namesResults
  mat
}

#gene_exp = data1
#log2FCData = data2
RaReNI <- function(data1, data2,
                   splitLog2FC=TRUE,
                   split_number=2,
                   outputDir="./RaReNIResults/",
                   iteration_start=1,
                   numCores=NULL,
                   verbose=FALSE) { #It could be something else


  #Check if input objects are defined
  do_objectsExist(data1,data2,verbose)

  #Check if output directory exists
  check_outputDir(outputDir,verbose)

  #Start measuring time
  start_time <- Sys.time()

  log2FC_data2<-as.data.frame(calculateLog2FoldChange(data2))

  #Make sure column names are appropiate
  colnames(data1)<-make.names(colnames(data1),unique=TRUE)
  colnames(log2FC_data2)<-make.names(colnames(log2FC_data2),unique=TRUE)

  #Split Data Frame
  if(splitLog2FC){

    #Split Log2FC data in parts
    ls_dfs<-split_df(log2FC_data2,split_number)
    print(paste("Log2FC Data will be split into",length(ls_dfs),"parts",sep=" "))

    #Loop Log2FC parts
    for (i in iteration_start:length(ls_dfs)){
      df_iter<-ls_dfs[[i]]
      print(paste('Running RaReNI Split Number',i,sep=" "))

      #Register parallel backend
      #library(doSNOW)
      if(is.null(numCores)){
        numCores<-detectCores()-2
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }else{
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }

      #Register parallel backend
      cl<-makeCluster(numCores)
      registerDoSNOW(cl)

      df_results = foreach(i = 1:ncol(df_iter),.packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.inorder = FALSE) %dopar% { #Loop log2FC data
        results2 =foreach(j = 1:ncol(data1),.packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.verbose = FALSE,.inorder = FALSE) %dopar% { #Loop genes

          #Get names
          GeneName<-colnames(data1)[j] #gene name
          LFC2Name<-colnames(df_iter)[i] #log2Fold change between metabolite pair...

          #Vectors
          geneVec<-data1[,j]
          LFCMetVec<-df_iter[,i]

          #Linear model
          #Linear model
          lmIter <- lm(LFCMetVec~geneVec)
          modelSummary <- summary(lmIter)  # capture model summary as an object
          modelCoeffs <- modelSummary$coefficients  # model coefficients
          beta.estimate <- modelCoeffs[get_variableName(geneVec), "Estimate"]  # get beta estimate for speed
          std.error <- modelCoeffs[get_variableName(geneVec), "Std. Error"]  # get std.error for speed
          t_value <- beta.estimate/std.error  # calc t statistic
          p_value <- modelCoeffs[get_variableName(geneVec), "Pr(>|t|)"] # calc p Value
          #p_adjust<-p.adjust(p_value)
          f_statistic <- modelSummary$fstatistic[[1]]  # fstatistic
          f <- summary(lmIter)$fstatistic  # parameters for model p-value calc
          model_p <- pf(f[1], f[2], f[3], lower=FALSE)[[1]] #model p-value
          #model_p_adjust<-p.adjust(model_p)
          rsquared <- summary(lmIter)$r.squared
          rsquare_adj <- summary(lmIter)$adj.r.squared

          rowtoprint<-list(LFC2Name,GeneName, #change to list instead of cbind.data.frame
                           beta.estimate,std.error,t_value,p_value,
                           f_statistic,model_p,rsquared,rsquare_adj)
        }
      }

      stopCluster(cl)

      df_results<-as.data.frame(df_results)
      colnames(df_results)<-c("MetabolitePair",	"gene_name",	"beta.estimate",	"std.error",	"t_value",	"p_value",	"f_statistic",	"model_pvalue","r_squared","adjusted_r_squared")
      df_results<-as.matrix(df_results)


      #Save result to memory
      writeT<-writeTable(df_results,num_cores = numCores,outputDir = outputDir,iteration = i) #Try to write using fwrite
      if(!length(writeT)==0){write.csv(df_results,paste(outputDir,"RaReNIResultsSplit",i,".csv",sep=""))}#If fwrite fails it is written with write.csv

      #Remove results
      rm(df_results)
    }
    #Merge output results CoNI
    RaReNIResults <- merge_outpuSplitFiles(outputDir)
    print('RaReNI ran successfully')
    return(RaReNIResults)
  }
}





#This function tries to write a table with fread
writeTable <- function(results_write,num_cores,outputDir,iteration) {
  out <- tryCatch(
    {
      suppressMessages(data.table::fwrite(results_write, paste(outputDir,"CoNIOutputSplit",iteration,".csv",sep=""),nThread=num_cores))
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

#This function checks previous files and deletes according to the User input option
#under construction... modify...
check_previous<-function(del,iteration,outDir,verb=verbose){
  if(del){
    filesDel<-list.files(outDir,pattern = "CoNIOutput")
    if(length(filesDel)>0){
      sapply(filesDel, function(f){file.remove(paste0(outDir,f))})
    }
    filesDel2<-list.files(outDir,pattern = "KeptFeatures_")
    if(length(filesDel2)>0){
      sapply(filesDel2, function(f){file.remove(paste0(outDir,f))})
    }
  }else if(iteration>1){
    #check if there are files... If there are none... Stop
    if(verb){
      cat("Previous files are present\n")
      cat("Starting from iteration ",iteration,"\n",sep="")
    }
  }else{
    filesDel<-list.files(outDir,pattern = "CoNIOutputSplit")
    if(length(filesDel)>0){
      sapply(filesDel, function(f){file.remove(paste0(outDir,f))})
    }
  }
}

delIntFiles<-function(del,outDir){
  if(del){
    filesDel<-list.files(outDir,pattern = "CoNIOutputSplit")
    sapply(filesDel, function(f){file.remove(paste0(outDir,f))})
  }

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
do_objectsExist<-function(gene_exp,norm_metabo_dat,verb=verbose){
  if(missing(gene_exp) | missing(norm_metabo_dat)){
    message("Input objects are missing")
    stop("CoNI end")
  }else{
    if(verb){
      print("Input objects are defined")
    }
  }
}


#This function wills split the 'big' omics data into smaller data frames to improve computations, avoid run out of memory
#The second option is probably better... need to review..
#In the second otion the data is splitted according to the number given + 1 that are the remaining of the division (modulo)
split_df<-function(AbundantDF,numberSplitDF_2=2){
  dt_list<-list()
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
    dt_list[[i]]<-AbundantDF[,start:ncol(AbundantDF),drop=FALSE] #Avoid losing column name when is a single column that is left
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
sig_correlation2<-function(input_driverD,padj=TRUE,method="BH", verb=verbose){
  corr<-Hmisc::rcorr(as.matrix(input_driverD),type='p')
  corr_table<-flattenCorrMatrix(corr$r,corr$P)
  corr_table$adj.p<-p.adjust(corr_table$p,method = method)

  if(padj){
    corr_tableSig <- corr_table %>% filter(adj.p<0.05)
    if(nrow(corr_tableSig) == 0){
      print('No features significantly correlate after padjustment for linkedD')
      print('Using non adjusted pvalues')
      corr_tableSig<-corr_table %>% filter(p<0.05)}
  }else{
    print("Ajustment for multiple testing was set to FALSE for correlations in linked Data")
    corr_tableSig<-corr_table %>% filter(p<0.05)
  }
  if(verb){print(paste('Significant correlations',nrow(corr_tableSig),sep=" "))}
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
#This part might be too slow... there are probably more efficient ways
  for(i in rows){
    for(j in cols){
      if (pvalueMatrix[i,j]>0.05){
        next
      }else{
        cor<-rcoeffMatrix[i,j]
        pvalue<-pvalueMatrix[i,j]
        df<- df %>% add_row(metabolite = i, gene = j, cor = cor, pvalue = pvalue) #This part might be inefficient and slow
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


#This function reads the output split files generated by CoNI and generates a single result object
merge_outpuSplitFiles<-function(outputDir){
  #outputDir<-gsub('\\.','',outputDir)
  #outputDir<-gsub('\\/','',outputDir)
  file_list <- list.files(outputDir)
  file_list<-file_list[grep("CoNIOutputSplit",file_list)]
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
check_outputDir<-function(outputDir,verb=verbose){
  if (file.exists(paste(outputDir,sep=""))) {
    if(verb){print("Output directory exists")}
  } else {
    print("Output directory does not exist - creating directory ... ")
    dir.create(file.path(paste(outputDir,sep="")))
  }
}

#Add if gene is a low expressed gene
#Expression table contains in one column the genes and a second column that specifies if the gene is a low expressed gene = LowExpression
#or not = ""
expression_level<-function(CoNI_results, expressionTable){
  expressionLevel<-expressionLevelTable[match(CoNI_results$Feature_driverD,expressionTable$Gene),"ExpressionLevel"]
  expressionLevel
}


#Create network using as input the output of CoNI
generate_network<-function(ResultsCoNI, colorNodesTable){
  results_SteigerAdjust <- ResultsCoNI[,1:3] #Get pair metabolites and gene

  #Summarize results for network construction
  df <- plyr::ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
              Genes=length(Feature_driverD),
              GenesString=paste0(unique(Feature_driverD),collapse=";"))
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
  results_SteigerAdjust <- ResultsCoNI[,c(1:7,10:11)] #Get pair metabolites, gene and pcor and cor information... change to add more information
  #Summarize results for network construction
  df<-plyr::ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
            weightreal=length(Feature_driverD),
            Genes=paste0(unique(Feature_driverD),collapse=";"),
            #ActualGeneNames=paste0(unique(ActualGeneName),collapse=";"),
            PcorValues=paste0(pcor_coefficient,collapse=";"),
            CorValues=paste0(cor_coefficient,collapse=";"),
            PcorAverage=mean(pcor_coefficient),
            CorAverage=mean(cor_coefficient),
            # PcorLink1Driver=paste0(Pcor_M1G_M2,collapse=";"),
            # PcorLink2Driver=paste0(Pcor_M2G_M1,collapse=";"),
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
  ls2 <- length(unique(ResultsCoNI$Feature_driverD)) #get number of genes affecting metabolites
  #Distance = 2 -> Second level neighborhood?
  df <- list()

  for(i in names(igraph::V(network))){ #loop nodes of graph
    l <- igraph::V(network)$name[neighbors(network, i)] #Get first level neighbors of node in iteration
    l1 <- list()
    for(j in l){
      l1[[j]] <- igraph::V(network)$name[neighbors(network, j)]#Loop neighbors of the node in the iteration and get their neighbors (Second level neighborhood)
    }
    l1 <- unique(unlist(l1)) #Get unique 2nd level neighbors
    s <- subset(ResultsCoNI, ((Feature_1_linkedD==i & Feature_2_linkedD %in% l) | (Feature_2_linkedD==i & Feature_1_linkedD %in% l)) |
                  ((Feature_1_linkedD %in% l & Feature_2_linkedD %in% l1) | (Feature_2_linkedD %in% l & Feature_1_linkedD %in% l1)))
    a <- length(unique(paste0(s$Feature_1_linkedD,"_",s$Feature_2_linkedD)))
    d <- length(unique(s$Feature_driverD))
    e <- nrow(s)
    s <- droplevels(s)
    b <- table(s$Feature_driverD)
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

tableLRDFs_LFs<-function(CoNIResults,LRDFs){
  CoNIResults_LRDFs<-CoNIResults[CoNIResults$Feature_driverD %in% LRDFs,]
  Gene_TableLRDFs<- plyr::ddply(CoNIResults_LRDFs, plyr::.(Feature_1_linkedD,Feature_2_linkedD), plyr::summarize,
                         Genes=paste(Feature_driverD,collapse=","))
  #Join Metabolite pairs
  CoNIResults_LRDFs_MetaboliteJoined<-tidyr::unite(CoNIResults_LRDFs,MetabolitePair,Feature_1_linkedD,Feature_2_linkedD,sep="-")
  CoNIResults_LRDFs_MetaboliteJoined<-CoNIResults_LRDFs_MetaboliteJoined[,c(1,2)]

  #Chowate table Genes and their corresponding Metabolite pairs
  LRDFs_and_MPairs <- plyr::ddply(CoNIResults_LRDFs_MetaboliteJoined, plyr::.(Feature_driverD), plyr::summarize,
                           MetabolitePairs=paste(MetabolitePair,collapse=","))
  #Temporary table
  temp<-as.data.frame(CoNIResults_LRDFs[,c(1:3)])

  #Add to the LRDFs and Metabolites pairs the unique individual metabolites
  LRDFs_and_MPairs$Metabolites<-plyr::ddply(temp, plyr::.(Feature_driverD), plyr::summarize,
                                     Metabolites=paste(unique(c(as.character(Feature_1_linkedD),as.character(Feature_2_linkedD))),collapse=","))[,2]
  colnames(LRDFs_and_MPairs)<-c("Local Regulated Driver Feature","Linked Feature pairs","Linked Features")
  LRDFs_and_MPairs
}
