#CoNI Functions

#' Correlation guided Network Integration
#' @description CoNI is the main function of Correlation guided Network Integration (CoNI). Input data should come from two sources (e.g., gene expression and metabolite expression), and it should come from the same samples. It calculates all pairwise correlations of the second data input elements and the partial correlations of these pairwise elements with respect to the elements of the first data input. Both data inputs can be prefiltered to include only those elements that significantly correlate. The first data input can be prefiltered to keep just low variance elements (var<0.5). A Steiger test is used to identify significant changes between the correlation and partial correlation values. Results can be visually represented in a Network.
#' @param edgeD Object to use as first data input (e.g., protein expression)
#' @param vertexD Object to use as second data input (e.g., metabolite expression)
#' @param outputDir Output Directory where results are stored
#' @param saveRaw logical. If TRUE the raw output of CoNI is saved in the output directory (outputDir)
#' @param outputNameRaw Name for the raw output file if saved
#' @param onlySgRes logical. If TRUE CoNI output is filtered and only significant results are kept
#' @param multipleTAdj logical. If TRUE it will filter results after adjustment of multiple testing
#' @param padjustvertexD logical. If TRUE vertexD is filtered according to the significant adjusted p-value of its pairwise correlations
#' @param correlateDFs logical. If TRUE the elements that significantly correlate of vertexD are correlated with the elements of edgeD. Only the elements that significantly correlate are kept
#' @param filter_highVarianceEdge logical. If TRUE features of edgeD with high variance are filtered out
#' @param splitedgeD logical. If TRUE edgeD will be split in n subsets for the computation (some instances n+1). Keep as TRUE unless the data input is small
#' @param split_number Number of parts to split the elements of edgeD
#' @param delPrevious logical. If TRUE previous files of a previous run are deleted
#' @param delIntermediaryFiles logical. If TRUE the output file of every iteration is deleted and only a single file with all results is kept
#' @param iteration_start Iteration start for CoNI. Useful if run is interrupted as one can restart from the last iteration
#' @param numCores Cores assigned for parallelization
#' @param verbose logical. If TRUE output in the console is more verbose
#' @param more_coef logical. If TRUE it will include the partial correlation of edge and vertex Features
#' @param edgeDname File name extension for the edge features that significantly correlate with at least one vertex feature. This file will be read if the function is called again with the same input and with delPrevious=FALSE
#' @param vertexDname File name extension for the vertex features that are involved in at least one significant correlation. This file will be read if the function is called again with the same input and with delPrevious=FALSE
#' @param saveFiles logical. If FALSE CoNI function will not save any file to disk
#' @return CoNI returns a data.frame with the correlation coefficients of the vertex-pairs, the partial correlation coefficients for every triplet, and the pvalue of the Steiger tests
#' @examples
#' #Run CoNI
#'
#' #Load gene expression - Toy dataset of two treatments
#' data(GeneExpToy)
#' #Samples in rows and genes in columns
#' GeneExp <- as.data.frame(t(GeneExpToy))
#' hfd_gene <- GeneExp[1:8,] #high fat diet
#' chow_gene<- GeneExp[9:nrow(GeneExp),] #chow diet
#' #Load metabolite expression - Toy dataset of two treatments
#' data(MetaboExpToy)
#' MetaboExp <- MetaboExpToy
#' hfd_metabo <- MetaboExp[11:18,] #high fat diet
#' chow_metabo <- MetaboExp[1:10,] #chow diet
#' #Match row names both data sets
#' rownames(hfd_metabo)<-rownames(hfd_gene)
#' rownames(chow_metabo)<-rownames(chow_gene)
#'
#' #Run CoNI with tiny example and no significance testing
#' #High fat diet
#' #For big datasets it is recommended to set splitedgeD to TRUE
#' #and split_number should be adjusted accordingly
#' #See vignette for an example
#' #Running CoNI with only a tiny dataset
#' \donttest{
#'  CoNIResultsHFD <- CoNI(hfd_gene,hfd_metabo,
#'                         numCores = 2,
#'                         onlySgRes = FALSE,
#'                         filter_highVarianceEdge=FALSE,
#'                         padjustvertexD = FALSE,
#'                         correlateDFs = FALSE,
#'                         edgeDname="HFD_genes",
#'                         vertexDname = "HFD_metabolites",
#'                         saveFiles = FALSE,
#'                         splitedgeD = FALSE,
#'                         outputDir = "./")
#'}
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import dplyr
#' @import ppcor
#' @import cocor
#' @importFrom data.table fwrite fread
#' @importFrom stats cor p.adjust
#' @importFrom utils write.csv
#' @export
CoNI<- function(edgeD, vertexD,
                outputDir = "./CoNIOutput/",
                saveRaw = TRUE, outputNameRaw = "CoNIOutput",
                onlySgRes = FALSE, multipleTAdj = TRUE,
                padjustvertexD = TRUE, correlateDFs=TRUE,
                filter_highVarianceEdge = TRUE,
                splitedgeD = TRUE, split_number = 2,
                delPrevious = FALSE, delIntermediaryFiles = TRUE,
                iteration_start = 1, numCores = NULL,
                verbose = TRUE,
                more_coef = FALSE,
                edgeDname = "edgeD",vertexDname = "vertexD",
                saveFiles = TRUE) {

  j <- NULL
  #Set delPrevious to FALSE if iteration start > 1
  if(iteration_start>1){
    if(verbose){cat('Iteration start > 1')}
    delPrevious<-FALSE
    file_list <- list.files(outputDir)
    file_list <- file_list[grep("CoNIOutputSplit", file_list)]
    if(length(file_list) < (iteration_start-1)){
      stop("Previous files were not found in the output directory")
    }
  }

  #Check parameters
  ParaList<-as.list(match.call())
  checkInputParameters(ParaList)

  #Check if input objects are defined
  do_objectsExist(edgeD, vertexD,verbose)

  #Start measuring time
  start_time <- Sys.time()

  #Check if output directory exists
  check_outputDir(outputDir, verbose)

  #Check if previous files are present and delete accordingly
  check_previous(delPrevious, iteration = iteration_start, outDir = outputDir, verb = verbose)

  #Split number cannot be less than 2
  if(split_number < 2){
    if(verbose){print("Split number cannot be less than two")}
    split_number<-2
    }

  #Test if sample names are the same in both data sets
  compare_sampleNames(edgeD, vertexD)

  #Make sure column names are appropiate
  colnames(edgeD) <- make.names(colnames(edgeD), unique = TRUE)
  colnames(vertexD) <- make.names(colnames(vertexD), unique = TRUE)

  if(!file.exists(paste(outputDir, "KeptFeatures_", vertexDname, ".csv", sep = ""))){
    #Get significant correlations between metabolites
    if(verbose){print("Calculating correlations of vertex Data")}
    normvertexD_Tablesignificant <- sig_correlation2(input_edgeD = vertexD,padj = padjustvertexD,verb = verbose)
    #Get indexes for the rows and columns for the metabo data
    normvertexD_Tablesignificant$RowIndex <- apply(normvertexD_Tablesignificant, 1, function(x){return(which(colnames(vertexD)[1:ncol(vertexD)]==x[1]))})
    normvertexD_Tablesignificant$ColIndex <- apply(normvertexD_Tablesignificant, 1, function(x){return(which(colnames(vertexD)[1:ncol(vertexD)]==x[2]))})
    if(saveFiles){
      fwrite(normvertexD_Tablesignificant, paste(outputDir, "KeptFeatures_", vertexDname, ".csv", sep=""))
      normvertexD_Tablesignificant <- fread(paste(outputDir, "KeptFeatures_", vertexDname, ".csv", sep=""))
    }
  }else{
    normvertexD_Tablesignificant <- fread(paste(outputDir, "KeptFeatures_", vertexDname, ".csv", sep=""))
  }

  #Get low variance edge features (e.g. genes)
  if(filter_highVarianceEdge){
    edgeD <- get_lowvarFeatures(edgeD) #This step was criticised
    if(!nrow(edgeD)>0){
      stop("After filtering high variance edge features none remained")
    }
  }
  #Remove those with too many zeros
  edgeD <- edgeD[, which(as.numeric(colSums(edgeD != 0)) > ceiling(nrow(edgeD)/2))] #At least two samples have a value higher than zero
  #edgeD<-as.data.frame(edgeD)

  #Get only those genes that correlate with the metabolites
  if(correlateDFs & !file.exists(paste(outputDir, "KeptFeatures_", edgeDname, ".csv", sep=""))){
    if(verbose){print("Calculating correlations between vertex Data and edge Data")}

    #Get Column indices of all metabolites
    metabo_indices<-unique(c(normvertexD_Tablesignificant$RowIndex,normvertexD_Tablesignificant$ColIndex))
    #Subset metabolites to correlate with genes
    SubSetvertexD<-vertexD[,metabo_indices]

    ResultsCorDfs <- sig_correlation2Dfs(SubSetvertexD,edgeD)
    genesSig <- unique(ResultsCorDfs$gene)
    if(verbose){print(paste(length(genesSig),"features were kept from edge Data",sep=" "))}
    edgeD <- edgeD[,genesSig]
    if(saveFiles){
      fwrite(edgeD,paste(outputDir,"KeptFeatures_",edgeDname,".csv",sep=""))
      edgeD <- fread(paste(outputDir,"KeptFeatures_",edgeDname,".csv",sep=""))
    }
  }else if(file.exists(paste(outputDir,"KeptFeatures_",edgeDname,".csv",sep=""))){
    edgeD <- fread(paste(outputDir,"KeptFeatures_",edgeDname,".csv",sep=""))
  }

  if(ncol(edgeD)>2000){
    if(splitedgeD==FALSE){
      cat('For computational purposes a split will be performed\n')
      splitedgeD<-TRUE
      split_number<-round(ncol(edgeD)*0.02)
    }
    if(split_number > ncol(edgeD)){
      if(verbose){print("Cannot split edgeD by the number provided, it exceeds edgeD dimensions")}
      split_number<-round(ncol(edgeD)*0.02)
      if(split_number<1){
        print("Cannot split less than 2")
        splitedgeD<-TRUE
        split_number<-round(ncol(edgeD)*0.02)
      }
    }
  }



  #Split Data Frame
  if(splitedgeD){

    ls_dfs<-split_df(edgeD,split_number)
    print(paste("Edge Data was split into",length(ls_dfs),"parts",sep=" "))

    for (i in iteration_start:length(ls_dfs)){
      df_iter<-ls_dfs[[i]]


      #Convert to data.frames
      df_iter<-as.data.frame(df_iter)
      normvertexD_Tablesignificant<-as.data.frame(normvertexD_Tablesignificant)


      #Register parallel backend
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        # use 2 cores in CRAN/Travis/AppVeyor
        numCores <- 2
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }else if(is.null(numCores)){
        numCores<-detectCores()-2
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }else{
        if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
      }


      if(verbose){print(paste('Running CoNI Split Number',i,sep=" "))}

      cl<-makeCluster(numCores,setup_timeout = 0.5)
      #registerDoSNOW(cl)
      registerDoParallel(cl)

      #Set progress bar
      # pb<-tkProgressBar(max=ncol(edgeD))

      #Run operations in parallel
      df_results = foreach(j = 1:ncol(df_iter), .packages = c("ppcor", "doParallel","cocor") , .combine=rbind,.inorder = FALSE) %dopar% { #Loop table significant metabolites %dopar% .options.snow=taskBar1(pb,ncol(edgeD))
        results2 =foreach(i = 1:nrow(normvertexD_Tablesignificant),.packages = c("ppcor", "doParallel","cocor"), .combine=rbind,.inorder = FALSE) %dopar% { #Loop genes

          index1<-normvertexD_Tablesignificant[i,6]#Index column of first vertex feature (e.g. metabolite)
          index2<-normvertexD_Tablesignificant[i,7]#Index column of second vertex feature (e.g. metabolite)

          #Get vertex features names and edge feature name (e.g. names for metabolites and gene)
          Feature_1_vertexD<-normvertexD_Tablesignificant[i,1]
          Feature_2_vertexD<-normvertexD_Tablesignificant[i,2]
          Feature_edgeD<-colnames(df_iter)[j]

          #Get correlation between vertex features (e.g. metabolites)
          cor_coefficient<-normvertexD_Tablesignificant[i,3]
          cor_pvalue<-normvertexD_Tablesignificant[i,4]

          #Calculate partial correlation between vertex features partialling out edge feature (e.g. metabolites and gene)
          pcor_result<-pcor.test(vertexD[,index1],vertexD[,index2],df_iter[,j],method="p")
          pcor_pvalue<-pcor_result[[2]]
          pcor_coefficient<-pcor_result[[1]]

          #Sometimes the computer is not precise in float representation...
          #For numbers very close to 1 and -1 it is problematic
          if(pcor_coefficient > 1){
            pcor_coefficient<-0.999
          }else if(pcor_coefficient < -1){
            pcor_coefficient<- -0.999
          }

          #Correlation vertex feature and edge feature (e.g metabolites vs gene)
          cor_m1_vs_g <- cor(vertexD[,index1],df_iter[,j])
          cor_m2_vs_g <- cor(vertexD[,index2],df_iter[,j])

          #Steiger test
          cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(vertexD),
                                           alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

          cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(vertexD),
                                            alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')

          cdgo_pvalue <- cdgo@steiger1980$p.value
          cdgo2_pvalue<- cdgo2@steiger1980$p.value




          #vertex Feature 1 and edge Feature partialling out vertex Feature 2 (e.g. Metabolite 1 and Gene | Metabolite 2)
          pcor_res_jh_k <- tryCatch({pcor.test(vertexD[,index1],df_iter[,j],vertexD[,index2],method="p")},
            error=function(cond) {
              message('Partial correlation LF1_edge|LF2 failed')
              message("Here's the original error message:")
              message(cond)
              # Choose a return value in case of error
              return('NA')
            }
          )
          if(is.na(pcor_res_jh_k[[1]])){
            pcor_res_jh_kCoef<-"NA"
            pcor_res_jh_kpval<-"NA"
          }else{
            pcor_res_jh_kCoef<-pcor_res_jh_k[[1]]
            pcor_res_jh_kpval<-pcor_res_jh_k[[2]]
          }


          #vertex Feature 2 and edge Feature partialling out vertex Feature 1 (e.g. Metabolite 2 and Gene | Metabolite 1)
          pcor_res_kh_j<-tryCatch({pcor.test(vertexD[,index2],df_iter[,j],vertexD[,index1],method="p")},
                   error=function(cond) {
                     message('Partial correlation LF2_edge|LF1 failed')
                     message("Here's the original error message:")
                     message(cond)
                     # Choose a return value in case of error
                     return('NA')
                   }
          )
          if(is.na(pcor_res_kh_j[[1]])){
            pcor_res_kh_jCoef<-"NA"
            pcor_res_kh_jpval<-"NA"
          }else{
            pcor_res_kh_jCoef<-pcor_res_kh_j[[1]]
            pcor_res_kh_jpval<-pcor_res_kh_j[[2]]
          }

          rowtoprint<-list(Feature_1_vertexD,Feature_2_vertexD,Feature_edgeD,
                           pcor_coefficient,pcor_pvalue,cor_coefficient,
                           cor_pvalue,cdgo_pvalue,cdgo2_pvalue,
                           pcor_res_jh_kCoef,pcor_res_jh_kpval,
                           pcor_res_kh_jCoef,pcor_res_kh_jpval)
        }

      }
      # close(pb)
      stopCluster(cl)

      df_results<-as.data.frame(df_results)

      #Set column names
      colnames(df_results)<-c("Feature_1_vertexD",	"Feature_2_vertexD",	"Feature_edgeD",	"pcor_coefficient",
                              "pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue",
                              "pcor_LF1_edge__LF2", "pcor_pvalue_LF1_edge__LF2",
                              "pcor_LF2_edge__LF1", "pcor_pvalue_LF2_edge__LF1")

      df_results<-as.matrix(df_results)

      #oldw <- getOption("warn")
      #options(warn = -1)

      #Save result to memory
      writeT<-writeTable(df_results,num_cores = numCores,outputDir = outputDir,iteration = i) #Try to write using fwrite
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

    #################
    #This step might be problematic with very big data, enough RAM is needed to avoid errors
    #Add adjusted steiger pvalue
    CoNIOutput$cdgo_adjusted<-p.adjust(CoNIOutput$cdgo_pvalue)
    CoNIOutput$cdgo2_adjusted<-p.adjust(CoNIOutput$cdgo2_pvalue)
    CoNIOutput<-as.data.frame(CoNIOutput)

    if(!more_coef){
      CoNIOutput<-CoNIOutput[,c(1:9,14:length(CoNIOutput))]
    }

    #Save raw results
    if(saveFiles & saveRaw){
      suppressMessages(fwrite(CoNIOutput, paste(outputDir,outputNameRaw,"_Raw",".gz",sep=""),nThread=numCores))
    }

    #Keep only significant results
    if(onlySgRes && multipleTAdj){ #adjustment for multiple testing
      CoNIOutput<-CoNIOutput %>% filter(cor_pvalue<=0.05) %>% filter(.data$cdgo_adjusted<=0.05 & .data$cdgo2_adjusted<=0.05)
    }else if(onlySgRes){ #without adjustment for multiple testing
      CoNIOutput<-CoNIOutput %>% filter(cor_pvalue<=0.05) %>% filter(.data$cdgo_pvalue<=0.05 & .data$cdgo2_pvalue<=0.05)
    }

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
    splitedgeD<-FALSE
  }

  if (splitedgeD==FALSE){
    #Register parallel backend
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      numCores <- 2
      if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
    }else if(is.null(numCores)){
      numCores<-detectCores()-2
      if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
    }else{
      if(verbose){cat("Running parallelization with ",numCores," cores\n",sep="")}
    }

    cl<-makeCluster(numCores)
    #registerDoSNOW(cl)
    registerDoParallel(cl)

    print('Running CoNI...')

    #Set progress bar
    # pb<-tkProgressBar(max=ncol(edgeD))

    #Convert to data.frames
    edgeD<-as.data.frame(edgeD)
    normvertexD_Tablesignificant<-as.data.frame(normvertexD_Tablesignificant)

    df_results = foreach(j = 1:ncol(edgeD),.packages = c("ppcor", "doParallel","cocor"), .combine=rbind,.inorder = FALSE) %dopar% {#Loop genes .options.snow=taskBar1(pb,ncol(edgeD))
      results2 = foreach(i = 1:nrow(normvertexD_Tablesignificant),.packages = c("ppcor", "doParallel","cocor") ,.combine=rbind,.inorder = FALSE) %dopar% {#Loop table significant metabolites

        index1<-normvertexD_Tablesignificant[i,6]#Index column of first metabolite
        index2<-normvertexD_Tablesignificant[i,7]#Index column of second metabolite

        #Get vertex features names and edge feature name (e.g. names for metabolites and gene)
        Feature_1_vertexD<-normvertexD_Tablesignificant[i,1]
        Feature_2_vertexD<-normvertexD_Tablesignificant[i,2]
        Feature_edgeD<-colnames(edgeD)[j]

        #Get correlation between vertex features (e.g. metabolites)
        cor_coefficient<-normvertexD_Tablesignificant[i,3]
        cor_pvalue<-normvertexD_Tablesignificant[i,4]

        #############################
        #Calculate partial correlation between vertex features partialing out edge feature (e.g. metabolites and gene)
        pcor_result<-pcor.test(vertexD[,index1],vertexD[,index2],edgeD[,j],method="p")
        pcor_pvalue<-pcor_result[[2]]
        pcor_coefficient<-pcor_result[[1]]

        #Sometimes the computer is not precise in float representation...
        #For numbers very close to 1 and -1 it is problematic
        if(pcor_coefficient > 1){
          pcor_coefficient<-0.999
        }else if(pcor_coefficient < -1){
          pcor_coefficient<- -0.999
        }

        #Correlation vertex feature and edge feature (e.g metabolites vs gene)
        cor_m1_vs_g <- cor(vertexD[,index1],edgeD[,j])
        cor_m2_vs_g <- cor(vertexD[,index2],edgeD[,j])

        #Test if partial correlation coefficient differs from correlation coefficient
        #j=vertex feature 1 (e.g. metabolite1)
        #k=vertex feature 1 (e.g.metabolite2)
        #h=edge feature (e.g. gene)

        #Steiger test
        cdgo <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m1_vs_g[1], n=nrow(vertexD),
                                         alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo2 <- cocor.dep.groups.overlap(r.jk=cor_coefficient[[1]], r.jh=pcor_coefficient, r.kh=cor_m2_vs_g[1], n=nrow(vertexD),
                                          alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0, test='steiger1980')
        cdgo_pvalue <- cdgo@steiger1980$p.value
        cdgo2_pvalue<- cdgo2@steiger1980$p.value


        #vertex Feature 1 and edge Feature partialling out vertex Feature 2 (e.g. Metabolite 1 and Gene | Metabolite 2)
        pcor_res_jh_k <- tryCatch({pcor.test(vertexD[,index1],edgeD[,j],vertexD[,index2],method="p")},
                                  error=function(cond) {
                                    message('Partial correlation LF1_edge|LF2 failed')
                                    message("Here's the original error message:")
                                    message(cond)
                                    # Choose a return value in case of error
                                    return('NA')
                                  }
        )
        if(is.na(pcor_res_jh_k[[1]])){
          pcor_res_jh_kCoef<-"NA"
          pcor_res_jh_kpval<-"NA"
        }else{
          pcor_res_jh_kCoef<-pcor_res_jh_k[[1]]
          pcor_res_jh_kpval<-pcor_res_jh_k[[2]]
        }

        #vertex Feature 2 and edge Feature partialling out vertex Feature 1 (e.g. Metabolite 2 and Gene | Metabolite 1)
        pcor_res_kh_j<-tryCatch({pcor.test(vertexD[,index2],edgeD[,j],vertexD[,index1],method="p")},
                                error=function(cond) {
                                  message('Partial correlation LF2_edge|LF1 failed')
                                  message("Here's the original error message:")
                                  message(cond)
                                  # Choose a return value in case of error
                                  return('NA')
                                }
        )
        if(is.na(pcor_res_kh_j[[1]])){
          pcor_res_kh_jCoef<-"NA"
          pcor_res_kh_jpval<-"NA"
        }else{
          pcor_res_kh_jCoef<-pcor_res_kh_j[[1]]
          pcor_res_kh_jpval<-pcor_res_kh_j[[2]]
        }


        rowtoprint<-list(Feature_1_vertexD,Feature_2_vertexD,Feature_edgeD,
                                     pcor_coefficient,pcor_pvalue,cor_coefficient,
                                     cor_pvalue,cdgo_pvalue,cdgo2_pvalue,
                                     pcor_res_jh_kCoef,pcor_res_jh_kpval,
                                     pcor_res_kh_jCoef,pcor_res_kh_jpval)

      }
    }

    # close(pb)
    stopCluster(cl)
    #Remove weird format
    df_results<-as.data.frame(df_results) #create data frame
    #Remove nested list format
    results<-list()
    for(i in 1:ncol(df_results)){#loop columns
      resultC<-unlist(df_results[i]) #unlist column
      results[[i]]<-unname(resultC) #save result to list
    }
    df_results<-as.data.frame(do.call(cbind,results)) #cbind list to create data frame


    # df_results<-sapply(df_results[,1:ncol(df_results)],function(x){
    #     x<-unlist(x)})
    # df_results<-as.data.frame(df_results)
    #Set numeric columns as numeric
    df_results[,4:13]<-sapply(df_results[, 4:13], function(x){
    as.numeric(as.character(x))
    })


    #Set column names
    colnames(df_results)<-c("Feature_1_vertexD",	"Feature_2_vertexD",	"Feature_edgeD",	"pcor_coefficient",
                            "pcor_pvalue",	"cor_coefficient",	"cor_pvalue",	"cdgo_pvalue",	"cdgo2_pvalue",
                            "pcor_LF1_edge__LF2", "pcor_pvalue_LF1_edge__LF2",
                            "pcor_LF2_edge__LF1", "pcor_pvalue_LF2_edge__LF1")


    #Add adjusted steiger pvalue
    df_results$cdgo_adjusted<-p.adjust(df_results$cdgo_pvalue)
    df_results$cdgo2_adjusted<-p.adjust(df_results$cdgo2_pvalue)

    #Save result to memory
    if(saveFiles & saveRaw){
      df_results_raw<-as.matrix(df_results)
      suppressMessages(fwrite(df_results_raw, paste(outputDir,outputNameRaw,"_Raw",".gz",sep=""),nThread=numCores,quote = TRUE))
    }

    if(!more_coef){
      df_results<-df_results[,c(1:9,14:length(df_results))]
    }

    #Filter significance
    if(onlySgRes){
      df_results<-df_results %>% filter(cor_pvalue<=0.05) %>% filter(.data$cdgo_adjusted<=0.05 & .data$cdgo2_adjusted<=0.05)
    }


    #Output processing time
    end_time <- Sys.time()
    total_time<-difftime(end_time,start_time,units='hours')
    cat(total_time,"hours", "\n",sep=" ")
    print('CoNI ran successfully')
    return(df_results)
  }

}

#' Check input parameters
#' @description Internal use. Function  to check if input parameters are of the right class
#' @keywords internal
#' @importFrom methods is
#' @return No return value, called for side effects
checkInputParameters<-function(ParaList){
  #Functions used
  matchLs<-function(L1,L2){
    Idx<-match(L1,L2)
    IdxOut<-Idx[!is.na(Idx)]
    IdxOut
  }

  #Check path
  ParamPathName<-c("outputDir")
  ParamPathL<-ParaList[matchLs(ParamPathName,names(ParaList))]

  if(length(ParamPathL)>0){
    param<-eval(ParamPathL[[1]])
    if(param=="./"){
      print("Output in working directory")
    }else{
      LPathParts<-strsplit(param,split = "/")[[1]]
      LPathParts<-LPathParts[-length(LPathParts)]
      Path<-paste(LPathParts,collapse="/")
      if(!dir.exists(Path)){
        stop("Path provided for the new directoy does not exist")
      }
    }
  }



  #Check parameters that should be characters
  ParamChNames<-  c("outputName","edgeDname","vertexDname")
  ChParaList<-ParaList[matchLs(ParamChNames,names(ParaList))]   #Obtain the parameter values given by the user
  ChNameList <-ParamChNames[matchLs(names(ParaList),ParamChNames)] #and names

  if(length(ChParaList)>0){
    for(i in 1:length(ChParaList)){#Loop parameters to make sure they are of class character
      param<-eval(ChParaList[[i]])
      if (!is(param, "character")) {
        stop(paste0("Wrong class for input '",ChNameList[i],"'. Character is expected and ",class(param)," was given"))
      }
    }
  }


  #Check parameters that should be nummeric
  ParamNumNames<-  c("split_number","numCores")
  NumParaList<-ParaList[matchLs(ParamNumNames,names(ParaList))]
  NumNameList <-ParamNumNames[matchLs(names(ParaList),ParamNumNames)]

  if(length(NumParaList)>0){
    for(i in 1:length(NumParaList)){#Loop parameters to make sure they are of class numeric
      param<-eval(NumParaList[[i]])
      if(is.null(param) && NumNameList[i] == "numCores"){
        next
      }else if(is.null(param)){
        stop(paste0("Input ",NumNameList[i]," cannot be null"))
      }

      if (!is(param, "numeric")) {
        stop(paste0("Wrong class for input '",NumNameList[i],"'. Numeric is expected and ",class(param)," was given"))
      }
      if(param<0){
        stop(paste0("Input '",NumNameList[i],"' is not a positive integer"))
      }
    }

  }

  #Check parameters that should be logical
  ParamLogNames<-  c("padjustvertexD","onlySgRes","correlateDFs","splitedgeD",
                     "delPrevious","delIntermediaryFiles","verbose",
                     "filter_highVarianceEdge","more_coef","saveRaw")
  LogParaList<-ParaList[matchLs(ParamLogNames,names(ParaList))]
  LogNameList <-ParamLogNames[matchLs(names(ParaList),ParamLogNames)]

  if(length(LogParaList)>0){
    for(i in 1:length(LogParaList)){#Loop parameters to make sure they are of class logical
      param<-eval(LogParaList[[i]])
      if (!is(param, "logical")) {
        stop(paste0("Wrong class for input '",LogNameList[i],"'. Logical is expected and ",class(param)," was given"))
      }
    }
  }
}

#' Write table
#' @description Internal use. This function tries to write a table with fread, if it fails it returns "NA"
#' @keywords internal
#' @importFrom data.table fwrite
#' @return Returns NA if it fails to write file to memory, otherwise no return value
writeTable <- function(results_write,num_cores,outputDir,iteration) {
  out <- tryCatch(
    {
      suppressMessages(fwrite(results_write, paste(outputDir,"CoNIOutputSplit",iteration,".csv",sep=""),nThread=num_cores))
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

#' Check previous files
#' @description Internal use. This function checks previous files and deletes according to the User input option
#' Requires some modifications...
#' @keywords internal
#' @return No return value, called for side effects
check_previous<-function(del,iteration,outDir,verb){
  if(del){
    filesDel<-unique(c(list.files(outDir,pattern = "^CoNIOutput"),list.files(outDir,pattern = "^.*_Raw")))
    if(length(filesDel)>0){
      sapply(filesDel, function(fi){file.remove(paste0(outDir,fi))})
    }
    filesDel2<-list.files(outDir,pattern = "KeptFeatures_")
    if(length(filesDel2)>0){
      sapply(filesDel2, function(fi){file.remove(paste0(outDir,fi))})
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
      sapply(filesDel, function(fi){file.remove(paste0(outDir,fi))})
    }
  }
}

#' Delete intermediary files
#' @description Internal use. This function deletes intermediary files.
#' @keywords internal
#' @return No return value, called for side effects
delIntFiles<-function(del,outDir){
  if(del){
    filesDel<-list.files(outDir,pattern = "CoNIOutputSplit")
    sapply(filesDel, function(fi){file.remove(paste0(outDir,fi))})
  }

}

#' Check if files exist
#' @description Internal use. This function tests if the input files exist, if they do not it will output an error and
#' end CoNI
#' @keywords internal
#' @return No return value, called for side effects
do_objectsExist<-function(gene_exp,norm_metabo_dat,verb){
  if(missing(gene_exp) | missing(norm_metabo_dat)){
    message("Input objects are missing")
    stop("CoNI end")
  }else{
    if(verb){
      print("Input objects are defined")
    }
  }
}

#' Split dataset
#' @description Internal use. This function wills split the 'big' omics data into smaller data frames to improve computations, to avoid running out of memory.
#' @keywords internal
#' @return A list of data.frames
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

#' Compare sample names
#' @description Internal use. This function compares the sample names between the two datasets provided. If names do not match CoNI stops and outputs an error message.
#' @keywords internal
#' @return No return value, called for side effects
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

#' Flatten
#' @description Internal use. This function gets the upper part of the matrix after calculating the correlation coefficients between all pairs of elements
#' cormat : matrix of the correlation coefficients
#' pmat : matrix of the correlation p-values
#' @keywords internal
#' @return A data.frame with all pairwise correlations of the input elements and their respective p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]], #Get rownames
    column = rownames(cormat)[col(cormat)[ut]], #Get colnames
    cor  =(cormat)[ut], #Get correlation coefficients
    p = pmat[ut] #Get p values
  )
}

#' Pairwise correlations
#' @description Internal use. This function calculates the pairwise correlations of a matrix (it uses Hmisc::rcorr function) and gets the significant correlations
#' @keywords internal
#' @importFrom Hmisc rcorr
#' @return A data.frame with the significant correlations of a correlation matrix
sig_correlation2<-function(input_edgeD,padj=TRUE,method="BH", verb){
  corr<-rcorr(as.matrix(input_edgeD),type='p')
  corr_table<-flattenCorrMatrix(corr$r,corr$P)
  corr_table$adj.p<-p.adjust(corr_table$p,method = method)

  if(padj){
    corr_tableSig <- corr_table %>% filter(.data$adj.p<0.05)
    if(nrow(corr_tableSig) == 0){
      print('No features significantly correlate after padjustment for vertexD')
      print('Using non adjusted pvalues')
      corr_tableSig<-corr_table %>% filter(.data$p<0.05)}
  }else{
    print("Ajustment for multiple testing was set to FALSE for correlations in vertex Data")
    corr_tableSig<-corr_table %>% filter(.data$p<0.05)
  }
  if(verb){print(paste('Significant correlations',nrow(corr_tableSig),sep=" "))}
  corr_tableSig
}

#' Significant correlations 2 Df
#' @description Internal use. This function input are two data frames (e.g. metabolites and genes). It calculates the correlation matrix and creates a table with only significant pairs. No correction for multiple testing is done
#' @keywords internal
#' @importFrom stats pt
#' @return A data.frame with vertex-edge significant correlation coefficients and their p-values
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
#This part is slow... needs to be improved
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

#' Low variance features
#' @description Internal use. Get low variance features.
#' @import genefilter
#' @importFrom genefilter rowVars
#' @keywords internal
#' @return A data.frame where low variance features were removed
get_lowvarFeatures<-function(df){
  Var<-NULL
  df<-as.data.frame(t(df))
  df$Var<-genefilter::rowVars(as.matrix(df))
  df<-subset(df,Var<0.5)
  df<-df[,-ncol(df)] #get rid of the column with the variances
  df<-as.data.frame(t(df))
}

#' Merge Files.
#' @description Internal use. This function reads the output split files generated by CoNI and generates a single result object. It is slow. Probably there are other faster alternatives.
#' @keywords internal
#' @importFrom data.table fread
#' @return A single data.table with the results of CoNI
merge_outpuSplitFiles<-function(outputDir){
  #outputDir<-gsub('\\.','',outputDir)
  #outputDir<-gsub('\\/','',outputDir)
  file_list <- list.files(outputDir)
  file_list<-file_list[grep("CoNIOutputSplit",file_list)]
  for (file in file_list){
    # if the merged dataset doesn't exist, create it
    if (!exists("datasetResultsCoNI")){
      datasetResultsCoNI <- fread(paste(outputDir,file,sep=""), header=TRUE, sep=",")
    }else{
      temp_dataset <- fread(paste(outputDir,file,sep=""), header=TRUE, sep=",")
      datasetResultsCoNI<-rbind(datasetResultsCoNI, temp_dataset)
      rm(temp_dataset)
    }
  }
  datasetResultsCoNI
}

#' Output directory
#' @description Internal use. This function checks if the output directory exists and if it does not, it will create it
#' @keywords internal
#' @return No return value, called for side effects
check_outputDir<-function(outputDir,verb){
  if (file.exists(paste(outputDir,sep=""))) {
    if(verb){print("Output directory exists")}
  } else {
    print("Output directory does not exist - creating directory ... ")
    dir.create(file.path(paste(outputDir,sep="")))
  }
}

#' Create network
#' @description This function creates a network using as input the output of CoNI and a table specifying the colors for the nodes.
#' @param ResultsCoNI The input of the function are the results of CoNI.
#' @param colorVertexNetwork logical. If TRUE, the table colorVertexTable has to be provided to specify vertex colors
#' @param colorVertexTable Table specifying the colors for the nodes (vertex features). The first column should contain the names matching the features of the vertex Data and the colors or other data can be specified in the rest of the columns
#' @param outputDir Output directory where the network is saved as well as the file that was used to generate the network.
#' @param outputFileName The name of the file used to create the network.
#' @param Class Optional data frame with at least two columns, first column contains all vertex features and another column the vertex feature class (column named "Class"). Necessary for treatment comparisons based on class
#' @param saveFiles logical. If FALSE TableForNetwork_`outputFileName`.csv and Network_`outputFileName`.graphml are not saved to disk
#' @importFrom plyr ddply
#' @importFrom igraph graph_from_data_frame simplify V E degree hub_score transitivity closeness betweenness eigen_centrality centr_betw centr_clo centr_degree edge_betweenness write_graph
#' @return Returns an igraph object (network) constructed from ResultsCoNI. The network includes the following network statistics
#' \itemize{
##'  \item{"degree"}{The number of the vertex adjacent edges}
##'  \item{"hub_score"}{The principal eigenvector of A*t(A), where A is the adjacency matrix of the graph}
##'  \item{"transitivity"}{Probability that the adjacent vertices of a vertex are connected}
##'  \item{"closeness"}{Steps required to access every other vertex from a given vertex}
##'  \item{"betweenness"}{(roughly) The number of geodesics (shortest paths) going through a vertex or an edge}
##'  \item{"eigen_centrality"}{Takes a graph (graph) and returns the eigenvector centralities of positions v within it}
##'  \item{"centralized_betweenness"}{The vertice-level centrality score according to the betweenness of vertices}
##'  \item{"centralized_closeness"}{The vertice-level centrality score according to the closeness of vertices}
##'  \item{"centralized_degree"}{The vertice-level centrality score according to the degrees of vertices}
##' }
##' For more details see igraph package
#' @examples
#' #Generate Network
#'
#' #Load color nodes table
#' data(MetColorTable)
#' #Assign colors according to "Class" column
#' MetColorTable<-assign_colorsAnnotation(MetColorTable)
#' #Load CoNI results
#' data(CoNIResultsHFDToy)
#'
#' #Generate Network
#' HFDNetwork<-generate_network(ResultsCoNI = CoNIResultsHFDToy,
#'                              colorVertexNetwork = TRUE,
#'                              colorVertexTable = MetColorTable,
#'                              outputDir = "./",
#'                              outputFileName = "HFD",
#'                              saveFiles = FALSE)
#' @export
generate_network<-function(ResultsCoNI,
                           colorVertexNetwork = TRUE,
                           colorVertexTable,
                           outputDir = "./", outputFileName="ResultsCoNI",
                           Class=NULL,
                           saveFiles = TRUE){
  Feature_edgeD <-pcor_coefficient <- cor_coefficient <- NULL
  pcor_LF1_edge__LF2 <-pcor_pvalue_LF1_edge__LF2 <-pcor_LF2_edge__LF1 <- pcor_pvalue_LF2_edge__LF1 <- NULL

  saveNetwork=TRUE
  #results_SteigerAdjust<-ResultsCoNI
  if(ncol(ResultsCoNI)>11){
    results_SteigerAdjust <- ResultsCoNI[,c(1:7,10:13)]
  }else{
    results_SteigerAdjust <- ResultsCoNI[,c(1:7)]
  }

   #Get pair metabolites, gene and pcor and cor information... change to add more information
  #Summarize results for network construction
  df<-ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
            weightreal=length(Feature_edgeD),
            Genes=paste0(unique(Feature_edgeD),collapse=";"),
            #ActualGeneNames=paste0(unique(ActualGeneName),collapse=";"),
            PcorValues=paste0(pcor_coefficient,collapse=";"),
            CorValues=paste0(cor_coefficient,collapse=";"),
            PcorAverage=mean(pcor_coefficient),
            CorAverage=mean(cor_coefficient)
            )
  if(ncol(ResultsCoNI)>11){
    df_2<-ddply(results_SteigerAdjust,c(1,2),plyr::summarize,
                      PcorLink1edge=paste0(pcor_LF1_edge__LF2,collapse=";"),
                      PcorLink1edge_pvalue=paste0(pcor_pvalue_LF1_edge__LF2,collapse=";"),
                      PcorLink2edge=paste0(pcor_LF2_edge__LF1,collapse=";"),
                      PcorLink2edge_pvalue=paste0(pcor_pvalue_LF2_edge__LF1,collapse=";"))
    df<-cbind(df,df_2[,3:ncol(df_2)])
  }



  colnames(df)[1:2] <- c("from","to")
  clinksd <- df
  clinksd$type <- "hyperlink"
  clinksd$weight <- clinksd$weightreal/max(clinksd$weightreal) #Calculate a width based on the maximum number of genes per connection
  #Save table
  if(!is.null(Class)){
    if(ncol(Class)<2){
      stop("The 'Class' data frame provided is not correct")
    }else if(!any(grepl("Class$",colnames(Class),ignore.case=TRUE))){
      stop("The 'Class' data frame does not contain a 'Class' column")
    }
    idxClass<-grep("Class$",colnames(Class),ignore.case = TRUE)
    clinksd$Vertex1_Class<-Class[match(clinksd$from,make.names(Class[,1])),idxClass] #make.names is necessary as the results files returns make.names feature names
    clinksd$Vertex2_Class<-Class[match(clinksd$to,make.names(Class[,1])),idxClass]
  }
  if(saveFiles){
    write.csv(clinksd,paste(outputDir,"TableForNetwork_",outputFileName,".csv",sep=""),row.names=FALSE)
  }

  cnodes <- data.frame("Name"=unique(c(as.character(df$from),as.character(df$to))),stringsAsFactors=FALSE)#Get the nodes (metabolites)

  if(colorVertexNetwork){
    #Assign colors to nodes
    m <- merge(cnodes,colorVertexTable,by.x="Name",by.y=colnames(colorVertexTable)[1],all.x=TRUE)
    cnodesd <- m
  }else{
    cnodesd <- cnodes
  }

  #Change column names
  #colnames(clinksd)[10] <- "weight"

  #Create graph
  netd <- igraph::graph_from_data_frame(d=clinksd, vertices=cnodesd, directed=FALSE)
  netd_simple <- igraph::simplify(netd,remove.multiple=FALSE)

  #Add network statistics
  igraph::V(netd_simple)$degree<-degree(netd_simple, mode="all")
  igraph::V(netd_simple)$hub_score<-hub_score(netd_simple, weights=NA)$vector
  igraph::V(netd_simple)$transitivity<-transitivity(netd_simple, type="local")
  igraph::V(netd_simple)$closeness<-closeness(netd_simple, mode="all", weights=NA)
  igraph::V(netd_simple)$betweenness<-betweenness(netd_simple, directed=FALSE, weights=NA)
  igraph::V(netd_simple)$eigen_centrality<-eigen_centrality(netd_simple, directed=FALSE, weights=NA)$vector
  igraph::V(netd_simple)$centralized_betweenness<-centr_betw(netd_simple, directed=FALSE, normalized=TRUE)$res
  igraph::V(netd_simple)$centralized_closeness<-centr_clo(netd_simple, mode="all", normalized=TRUE)$res
  igraph::V(netd_simple)$centralized_degree<-centr_degree(netd_simple, mode="all", normalized=TRUE)$res
  #V(netd_simple)$membership_community_edgeBetweenes<-cluster_edge_betweenness(netd_simple,directed = F)$membership

  #Add edge betweeness
  igraph::E(netd_simple)$betweeness <- edge_betweenness(netd_simple, directed=FALSE, weights=NA)

  if(saveFiles & saveNetwork){
    write_graph(netd_simple,file=paste0(outputDir,"Network_",outputFileName,".graphml"),format="graphml")
  }
  return(netd_simple)
}

#' Find local controlling features
#' @description This function applies for a selected subnetwork a binomial test using the frequency of appearance of an edge feature and the total number of edge features. The probability corresponds to 1/n_df, where n_df corresponds to the total number of edge features in the network.
#' The selected subnetwork corresponds to the second level neighborhood of a specific node. The test is applied to all possible second level neighborhoods in the network.
#' @param ResultsCoNI The output of CoNI (after p-adjustment)
#' @param network Network created with the function generate_network
#' @param padjust logical. Filter output based on adjusted p values
#' @return Returns a data.frame with the results of the binomial tests. Significant results correspond to local controlling features
#' @importFrom igraph V neighbors
#' @importFrom stats dbinom
#' @import dplyr
#' @examples
#' #Load color nodes table
#' data(MetColorTable)
#'
#' #Assign colors according to "Class" column
#' MetColorTable<-assign_colorsAnnotation(MetColorTable)
#'
#' #Load CoNI results
#' data(CoNIResultsHFDToy)
#'
#' #Generate Network
#' #Note: Colors not visible when ploting in Igraph
#' HFDNetwork<-generate_network(ResultsCoNI = CoNIResultsHFDToy,
#'                              colorVertexNetwork = TRUE,
#'                              colorVertexTable = MetColorTable,
#'                              Class = MetColorTable,
#'                              outputDir = "./",
#'                              outputFileName = "HFD",
#'                              saveFiles = FALSE)
#'
#'#Note: For this tiny example nothing is significant
#' LCG_BinomialTestTableHFD<- find_localControllingFeatures(ResultsCoNI = CoNIResultsHFDToy,
#'                                                          network = HFDNetwork )
#' LCGenes_HFD<-as.character(unique(LCG_BinomialTestTableHFD$edgeFeatures))
#'
#' @export
find_localControllingFeatures<-function(ResultsCoNI,network,padjust=TRUE){
  Feature_1_vertexD <- Feature_2_vertexD <- NULL

  ls2 <- length(unique(ResultsCoNI$Feature_edgeD)) #get number of genes affecting metabolites
  #Distance = 2 -> Second level neighborhood?
  df <- list()

  for(i in names(igraph::V(network))){ #loop nodes of graph
    l <- igraph::V(network)$name[neighbors(network, i)] #Get first level neighbors of node in iteration
    l1 <- list()
    for(j in l){ #loop the first neighbors and get their neighbors (Second level neighborhood)
      l1[[j]] <- igraph::V(network)$name[neighbors(network, j)]
    }
    l1 <- unique(unlist(l1)) #Get unique 2nd level neighbors
    #Subset the CoNI Results table to include only the second level neighborhood
    s <- subset(ResultsCoNI, ((Feature_1_vertexD==i & Feature_2_vertexD %in% l) | (Feature_2_vertexD==i & Feature_1_vertexD %in% l)) |
                  ((Feature_1_vertexD %in% l & Feature_2_vertexD %in% l1) | (Feature_2_vertexD %in% l & Feature_1_vertexD %in% l1)))
    #Get the total number of edges in the neighborhood
    EdgesNo <- length(unique(paste0(s$Feature_1_vertexD,"_",s$Feature_2_vertexD)))
    #Get the unique total number of edge features (e.g., genes) in the neighborhood
    DrF_totalNo <- length(unique(s$Feature_edgeD))
    #The amount of edge features (e.g., genes) (with repetitions) found in the second level neighborhood
    DrF_wRepNo <- nrow(s)
    s <- droplevels(s)
    #The number of times every edge feature (e.g., genes) appears in the neighborhood. It is a table.
    b <- table(s$Feature_edgeD)
    TotalNumberGenes<-length(unique(s$Feature_edgeD))
    df[[i]] <- data.frame("Node1"=rep(i,DrF_totalNo),"Edges"=rep(EdgesNo,DrF_totalNo),"Draws"=rep(DrF_wRepNo,DrF_totalNo),"GenesInTotal"=rep(DrF_totalNo,DrF_totalNo),as.data.frame(b))
  }
  #Generate result data frame
  res2 <- do.call(rbind.data.frame, df)

  #we use the binomial distribution to test if the enrichment is significant as we can draw a gene for an area more often
  res2$Pval <- apply(res2,1,function(x){dbinom(as.numeric(x[[6]]),as.numeric(x[[3]]),1/ls2)})
  res2$Padj <- p.adjust(res2$Pval)
  res2 <- res2[order(res2$Padj),]
  res2 <- res2 %>% rename(edgeFeatures = .data$Var1)

  if(padjust){
    res2 <- subset(res2,res2$Padj<0.05)
  }else{
    res2 <- subset(res2,res2$Pval<0.05)
  }
  res2
}

#'Linker Features by magnitude of effect
#'@description This function outputs the linker features with the strongest effect on the correlation of the vertex features
#'@param ResultsCoNI The output of CoNI
#'@param topn Top n number of features to output
#'@return Returns a data.frame, a filtered version of ResultsCoNI, showing the top n features
#'with the strongest effect, that is, the highest difference between the partial correlation and correlation coefficient.
#'@importFrom rlang .data
#'@examples
#' data(CoNIResultsHFDToy)
#' Top10HFD<-top_n_LF_byMagnitude(CoNIResults_HFD,topn = 10)
#'@export
top_n_LF_byMagnitude<-function(ResultsCoNI, topn=10){
  ResultsCoNI<-ResultsCoNI %>% mutate(difference=abs(.data$cor_coefficient - .data$pcor_coefficient)) %>% arrange(desc(.data$difference))
  lEdgeFeatures<-unique(ResultsCoNI$Feature_edgeD)
  if(length(lEdgeFeatures)>=topn){
    selectedEdgeFeatures<-lEdgeFeatures[1:topn]
  }else{
    selectedEdgeFeatures<-lEdgeFeatures[1:length(selectedEdgeFeatures)]
  }
  Out<-as.data.frame(ResultsCoNI[ResultsCoNI$Feature_edgeD %in% selectedEdgeFeatures,])
  return(Out)
}

#' Table local controlling edge features and vertex pairs
#' @description This function creates a table of the local controlling edge features
#' @param CoNIResults The output of CoNI (after p-adjustment)
#' @param LCFs Local controlling edge features as a vector
#' @return A data.frane of local controlling edge features and their respective vertex pairs, and unique vertexes.
#' @examples
#' #Load CoNI results
#' data(CoNIResultsHFDToy)
#' #Note: arbitrary list of genes, not Local controlling features
#' tableLCFs_VFs(CoNIResultsHFDToy, c("Lilr4b","Rps12"))
#' @importFrom plyr ddply
#' @importFrom tidyr unite
#' @export
tableLCFs_VFs<-function(CoNIResults,LCFs){
  Feature_1_vertexD<-Feature_2_vertexD<-Feature_edgeD<-MetabolitePair<-NULL

  CoNIResults_LCFs<-CoNIResults[CoNIResults$Feature_edgeD %in% LCFs,]
  Gene_TableLCFs<- plyr::ddply(CoNIResults_LCFs, plyr::.(Feature_1_vertexD,Feature_2_vertexD), plyr::summarize,
                         Genes=paste(Feature_edgeD,collapse=","))
  #Join Metabolite pairs
  CoNIResults_LCFs_MetaboliteJoined<-tidyr::unite(CoNIResults_LCFs,MetabolitePair,Feature_1_vertexD,Feature_2_vertexD,sep="-")
  CoNIResults_LCFs_MetaboliteJoined<-CoNIResults_LCFs_MetaboliteJoined[,c(1,2)]

  #Chowate table Genes and their corresponding Metabolite pairs
  LCFs_and_MPairs <- plyr::ddply(CoNIResults_LCFs_MetaboliteJoined, plyr::.(Feature_edgeD), plyr::summarize,
                           MetabolitePairs=paste(MetabolitePair,collapse=","))
  #Temporary table
  temp<-as.data.frame(CoNIResults_LCFs[,c(1:3)])

  #Add to the LCFs and Metabolites pairs the unique individual metabolites
  LCFs_and_MPairs$Metabolites<-plyr::ddply(temp, plyr::.(Feature_edgeD), plyr::summarize,
                                     Metabolites=paste(unique(c(as.character(Feature_1_vertexD),as.character(Feature_2_vertexD))),collapse=","))[,2]
  colnames(LCFs_and_MPairs)<-c("Local Controlling edge Feature","Vertex Feature pairs","Vertex Features")
  LCFs_and_MPairs
}

#' Compare triplets
#' @description Compare vertexFeature-vertexFeature-edgeFeature between two treatments, that is, find the shared triplets between two different CoNI runs.
#' @param Treat1_path TableForNetwork_file1 (file generated by CoNI) with path for Treatment 1
#' @param Treat2_path TableForNetwork_file2 (file generated by CoNI) with path for Treatment 2
#' @param OutputName Output file name with path
#' @return A data.frame with the shared triplets (vertex1 vertex2 edge_feature) between two CoNI runs
#' @examples
#' #For an example see the vignette
#' @importFrom utils write.csv
#' @export
Compare_Triplets<-function(Treat1_path,Treat2_path,
                           OutputName="Shared_Genes_and_Edges_Treat1vsTreat2.csv"){
  path_C<-file.path(find.package("CoNI"),"python")
  runPython<-tryCatch({system(paste0('python3 ',path_C,'/Compare_Triplets.py ',Treat1_path," ",Treat2_path," ",OutputName))},
           error=function(cond) {
             # Choose a return value in case of error
             return('Error')
           }
  )
  if(runPython==9009 || runPython == 127 ||  runPython == 2){
    runPython<-tryCatch({system(paste0('python ',path_C,'/Compare_Triplets.py ',Treat1_path," ",Treat2_path," ",OutputName))},
                        error=function(cond) {
                                   # Choose a return value in case of error
                                   return('Error')
                                 }
    )
    if(runPython==9009 || runPython == 127 ||  runPython == 2){stop("Make sure python3 is installed and in your path")}
  }
  # system(paste0('python3 ',path_C,'/Compare_Triplets.py ',Treat1_path," ",Treat2_path," ",OutputName))
  Output<-read.csv(OutputName,sep="\t")
  return(Output)
}

#' Table VertexClass pairs of shared Edge Features
#' @description Compare VertexClass pairs of the shared Edge Features of two treatments (e.g., lipid-class-pairs per shared gene)
#' @param Treat1_path TableForNetwork_file (file generated by CoNI) with path of Treatment 1
#' @param Treat2_path TableForNetwork_file (file generated by CoNI) with path of Treatment 2
#' @param OutputName Output file name with path
#' @param Treat1Name Name of treatment one, default Treat1
#' @param Treat2Name Name of treatment one, default Treat2
#' @return A data.frame with all possible vertex-class pairs and their numbers per edge-feature and treatment.
#' @examples
#' #For an example see the vignette
#' @importFrom utils read.csv
#' @export
Compare_VertexClasses_sharedEdgeFeatures<-function(Treat1_path,Treat2_path,OutputName="Shared_Genes_and_Edges_Treat1vsTreat2.csv",Treat1Name="Treat1",Treat2Name="Treat2"){
  are_ClassColumnsPresent<-function(DFTreatment){
    boolVec<-c("Vertex1_Class","Vertex2_Class") %in% colnames(DFTreatment)
    if(!all(boolVec)){
      stop("Error: Make sure to add Class when you create your networks")
    }
  }

  DFTreat1<-read.csv(Treat1_path,nrows = 2,header=TRUE)
  DFTreat2<-read.csv(Treat1_path,nrows = 2,header=TRUE)
  are_ClassColumnsPresent(DFTreat1)
  are_ClassColumnsPresent(DFTreat2)


  path_C<-file.path(find.package("CoNI"),"python")
  runPython<-tryCatch({system(paste0('python3 ',path_C,'/ComparisonClasses.py ',Treat1_path," ",Treat2_path," ",OutputName," ",Treat1Name," ",Treat2Name))},
                      error=function(cond) {
                        # Choose a return value in case of error
                        return('Error')
                      }
  )
  if(runPython==9009 || runPython == 127 ||  runPython == 2){
    runPython<-tryCatch({system(paste0('python ',path_C,'/ComparisonClasses.py ',Treat1_path," ",Treat2_path," ",OutputName," ",Treat1Name," ",Treat2Name))},
                        error=function(cond) {
                          # Choose a return value in case of error
                          return('Error')
                        }
    )
    if(runPython==9009 || runPython == 127 ||  runPython == 2){stop("Make sure python3 is installed and in your path")}
  }
  # system(paste0('python3 ',path_C,'/ComparisonClasses.py ',Treat1_path," ",Treat2_path," ",OutputName," ",Treat1Name," ",Treat2Name))
  Output<-read.csv(OutputName,sep="\t")
  Output[,Treat1Name]<-as.numeric(gsub("Pair Class Missing",0,Output[,Treat1Name]))
  Output[,Treat2Name]<-as.numeric(gsub("Pair Class Missing",0,Output[,Treat2Name]))
  return(Output)
}

#' Vertex-class pairs profile of one shared edge feature
#' @description This function will create a barplot from the output of Compare_VertexClasses_sharedEdgeFeatures for a specific shared Edge Feature (e.g., a shared gene).
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param edgeF Edge feature present in output of Compare_VertexClasses_sharedEdgeFeatures
#' @param treat1 Name of treatment one, default Treatment1. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param treat2 Name of treatment one, default Treatment2. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param factorOrder A list specifying the order of the treatments.
#' @param col1 Color for Treatment 1
#' @param col2 Color for Treatment 2
#' @param EdgeFeatureType Type of Edge Feature (e.g., Gene)
#' @param xlb Name for x-axis
#' @param ylb Name for the y-axis
#' @param szaxisTxt Size axis text
#' @param szaxisTitle Size axis titles
#' @export
#' @return A ggplot object for a barplot. The barplot shows the vertex-class pairs profile of a single shared edge feature between two treatments
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' create_edgeFBarplot(CompTreatTable = VertexClassesSharedGenes_HFDvsChow,
#'                     edgeF = "Lilr4b",
#'                     treat1 = "HFD",
#'                     treat2 = "Chow",
#'                     factorOrder = c("HFD","Chow"),
#'                     EdgeFeatureType = "Gene")
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom forcats fct_relevel
#' @importFrom tidyselect vars_select_helpers
create_edgeFBarplot<-function(CompTreatTable,edgeF,treat1="Treatment1",treat2="Treatment2",
                              factorOrder=NULL,col1="red",col2="blue",EdgeFeatureType="Edge Feature",
                              xlb="Vertex-Class Pairs",
                              ylb="Number of pairs",
                              szaxisTxt=12,szaxisTitle=12){

  treatment<-number_pairs<-VertexClassPair<-NULL

  CompTreatTable[,treat1]<-gsub("Pair Class Missing",0,CompTreatTable[,treat1])
  CompTreatTable[,treat2]<-gsub("Pair Class Missing",0,CompTreatTable[,treat2])

  CompTreatTableF<-CompTreatTable %>% filter(.data$EdgeFeature==edgeF)
  CompTreatTableF<-CompTreatTableF[,c(2:4)]
  #Make sure columns are numeric
  CompTreatTableF[,2]<-as.numeric(CompTreatTableF[,2]) #Treat1
  CompTreatTableF[,3]<-as.numeric(CompTreatTableF[,3]) #Treat2
  #Relationship of edge features with lipid classes for specific gene
  CompTreatTableF_VertexClasses<-CompTreatTableF %>% group_by(.data$VertexClassPair) %>% summarise(across(tidyselect::vars_select_helpers$where(is.numeric),sum))
  ResultsFBarplot <- gather(CompTreatTableF_VertexClasses, treatment, number_pairs,-VertexClassPair)
  #Reorder factors
  if(!is.null(factorOrder)){
    ResultsFBarplot <- ResultsFBarplot %>% mutate(treatment = fct_relevel(.data$treatment,factorOrder))
  }

  #Create bar plot
  p<-ggplot(ResultsFBarplot, aes(x=.data$VertexClassPair, y=.data$number_pairs, fill=.data$treatment)) +
    geom_bar(width = 0.4,stat="identity",position="dodge") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5,size=20),
          axis.title=element_text(size = szaxisTitle,face="bold", colour = "black"),
          axis.text = element_text(size = szaxisTxt))+
    scale_fill_manual("treatment", values = c(col1,col2))+
    geom_text(aes(label=.data$number_pairs),size=3, position=position_dodge(width=0.9), vjust=-0.25)+
    ggtitle(paste0(EdgeFeatureType,": ",edgeF))+
    xlab(xlb)+
    ylab(ylb)
  return(p)
}

#' Vertex-class pairs profile of shared features
#' @description This function will create a barplot from the output of Compare_VertexClasses_sharedEdgeFeatures using all shared Edge Features (e.g., genes).
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param treat1 Name of treatment one, default Treatment1. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param treat2 Name of treatment one, default Treatment2. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param factorOrder A list specifying the order of the treatments.
#' @param col1 Color for Treatment 1
#' @param col2 Color for Treatment 2
#' @param maxpairs If number of class-vertex-pairs > maxpairs, display number pairs on top of bar
#' @param xlb Name for x-axis
#' @param ylb Name for the y-axis
#' @param szggrepel Size ggrepel labels
#' @param nudgey Nudge y ggrepel
#' @param nudgex Nudge x ggrepel
#' @param szaxisTxt Size axis text
#' @param szaxisTitle Size axis title
#' @export
#' @return A ggplot object for a barplot. The barplot shows the vertex-class pairs profile of all shared edge features between treatments
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' create_GlobalBarplot(CompTreatTable = VertexClassesSharedGenes_HFDvsChow,
#'                      treat1 = "HFD",
#'                      treat2 = "Chow",
#'                      factorOrder = c("HFD","Chow"),
#'                      col1="red",
#'                      col2 ="blue",
#'                      maxpairs = 1,
#'                      szggrepel = 6,
#'                      szaxisTxt = 15,
#'                      szaxisTitle = 15,
#'                      xlb = "Metabolite-pair classes")
#' @import ggplot2
#' @import ggrepel
#' @importFrom tidyr gather
#' @importFrom forcats fct_relevel
#' @importFrom rlang .data
#' @importFrom tidyselect vars_select_helpers
create_GlobalBarplot<-function(CompTreatTable,
                               treat1="Treatment1",
                               treat2="Treatment2",
                               factorOrder=NULL,
                               col1="red",
                               col2="blue",
                               maxpairs=1,
                               xlb="Vertex-Class Pairs",
                               ylb="Number of pairs",
                               szggrepel =3.5,
                               nudgey=0.5,
                               nudgex=0.5,
                               szaxisTxt=12,
                               szaxisTitle=12){

  treatment <- number_pairs <- VertexClassPair <- NULL

  CompTreatTable[,treat1]<-gsub("Pair Class Missing",0,CompTreatTable[,treat1])
  CompTreatTable[,treat2]<-gsub("Pair Class Missing",0,CompTreatTable[,treat2])

  #Get rid of edge features, as we want a global view
  CompTreatTable_NoEdgeFeatures<-CompTreatTable[,c(2:4)]
  #Make sure columns are numeric
  CompTreatTable_NoEdgeFeatures[,2]<-as.numeric(CompTreatTable_NoEdgeFeatures[,2]) #Treat1
  CompTreatTable_NoEdgeFeatures[,3]<-as.numeric(CompTreatTable_NoEdgeFeatures[,3]) #Treat2
  #Global view of the relationship of edge features with lipid classes
  CompTreatTable_VertexClasses<-CompTreatTable_NoEdgeFeatures %>% group_by(.data$VertexClassPair) %>% summarise(across(tidyselect::vars_select_helpers$where(is.numeric),sum))
  GlobalResultsFBarplot <- tidyr::gather(CompTreatTable_VertexClasses, treatment, number_pairs,-VertexClassPair)
  #Reorder factors
  if(!is.null(factorOrder)){
    GlobalResultsFBarplot <- GlobalResultsFBarplot %>% mutate(treatment = forcats::fct_relevel(.data$treatment,factorOrder))
  }

  #Create bar plot
  p<-ggplot(GlobalResultsFBarplot, aes(x=.data$VertexClassPair, y=.data$number_pairs, fill=.data$treatment)) +
    geom_bar(width = 0.4,stat="identity",position="dodge") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.title=element_text(size = szaxisTitle,face="bold", colour = "black"),
          axis.text = element_text(size = szaxisTxt))+
    geom_text(data = subset(GlobalResultsFBarplot, number_pairs >= maxpairs & treatment==treat1),
              aes(label=.data$number_pairs),
              show.legend  = FALSE ,
              size=szggrepel,
              position=position_dodge(width=0.9),
              vjust=-0.25)+
    geom_text_repel(data = subset(GlobalResultsFBarplot, number_pairs >=  maxpairs & treatment==treat2),
                    aes(label=.data$number_pairs),
                    show.legend = FALSE,
                    size = szggrepel,
                    min.segment.length = unit(0, 'lines'),
                    hjust=0,
                    nudge_y = nudgey,
                    nudge_x = nudgex,
                    direction="y",
                    point.padding =NA,
                    force = 0.1,
                    segment.alpha=0.3,
                    max.overlaps=Inf)+
    scale_fill_manual("treatment", values = c(col1,col2))+
    scale_colour_manual(values=c(col2, col1))+
    xlab(xlb)+
    ylab(ylb)+
    coord_cartesian(clip = "off")
  return(p)
}

#'Labels to colors
#' @description Internal use. This function is modified version of labels2colors of WGCNA and the internet
#' @keywords internal
#' @return A character vector with colors
labels2colors_2<-function (labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey", commonColorCode = TRUE) {
  standardColors<-function (n = NULL) {
    if (is.null(n))
      return(.GlobalStandardColors)
    if ((n > 0) && (n <= length(.GlobalStandardColors))) {
      return(.GlobalStandardColors[c(1:n)])
    }
    else {
      stop("Invalid number of standard colors requested.")
    }
  }

  # This code forms a vector of color names in which the first entries are given by BaseColors and the rest
  # is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".
  BaseColors = c("turquoise","blue","brown","yellow","green","red","pink","magenta",
                 "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
                 "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
                 "darkturquoise", "darkgrey",
                 "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue",
                 "paleturquoise", "violet", "darkolivegreen", "darkmagenta" );

  RColors = colors()[-grep("grey", colors())];
  RColors = RColors[-grep("gray", RColors)];
  RColors = RColors[-grep("black", RColors)];
  InBase = match(BaseColors, RColors);
  ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
  nExtras = length(ExtraColors);

  # Here is the vector of colors that should be used by all functions:
  .GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:nExtras) +sin(13*c(1:nExtras))) )] );

  if (is.null(colorSeq))
    colorSeq = standardColors()
  if (is.numeric(labels)) {
    if (zeroIsGrey)
      minLabel = 0
    else minLabel = 1
    if (any(labels < 0, na.rm = TRUE))
      minLabel = min(c(labels), na.rm = TRUE)
    nLabels = labels
  }
  else {
    if (commonColorCode) {
      factors = factor(c(as.matrix(as.data.frame(labels))))
      nLabels = as.numeric(factors)
      dim(nLabels) = dim(labels)
    }
    else {
      labels = as.matrix(as.data.frame(labels))
      factors = list()
      for (c in 1:ncol(labels)) factors[[c]] = factor(labels[, c])
      nLabels = sapply(factors, as.numeric)
    }
  }
  if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
    nRepeats = as.integer((max(labels) - 1)/length(colorSeq)) + 1
    warning(paste("labels2colors: Number of labels exceeds number of avilable colors.",
                  "Some colors will be repeated", nRepeats, "times."))
    extColorSeq = colorSeq
    for (rep in 1:nRepeats) extColorSeq = c(extColorSeq, paste(colorSeq, ".", rep, sep = ""))
  }
  else {
    nRepeats = 1
    extColorSeq = colorSeq
  }
  colors = rep("grey", length(nLabels))
  fin = !is.na(nLabels)
  colors[!fin] = naColor
  finLabels = nLabels[fin]
  colors[fin][finLabels != 0] = extColorSeq[finLabels[finLabels != 0]]
  if (!is.null(dim(labels)))
    dim(colors) = dim(labels)
  colors
}

#' Get class rgb color
#' @description Internal use. This function gets the class rgb color of a specific (lipid) class
#' @keywords internal
#' @return A character object, that corresponds to a color in hexadecimal format
getcolor<-function(ClassM,tableColor){
  IDxclass<-grep("class$",colnames(tableColor),ignore.case = TRUE)
  IDxrgb<-grep("rgb",colnames(tableColor),ignore.case = TRUE)
  clhex<-unique(tableColor[which(tableColor[,IDxclass]==ClassM),IDxrgb])
  return(clhex)
}

#' Assing Colors to Class
#'@description This function assigns two color columns (color name and rgb) to an annotation data frame according to a column named 'Class' or 'class'
#'@param AnnotationDf Annotation data frame that contains a factor variable to use to assign colors
#'@param col  Column with factor variable that will be used to assign colors
#'@export
#'@return The input data.frame with two extra columns specifying the colors
#'        for all vertexes according to their respective vertex-class
#'@importFrom gplots col2hex
assign_colorsAnnotation<-function(AnnotationDf,col="Class"){
  IDxclass<-grep(paste0(col,"$"),colnames(AnnotationDf),ignore.case = TRUE)
  AnnotationDf$Color<-labels2colors_2(as.numeric(as.factor(AnnotationDf[,IDxclass])))
  AnnotationDf$ColorRgb<-col2hex(AnnotationDf$Color)
  return(AnnotationDf)
}

#' Get colors
#' @description Internal use. This function gets the rgb colors of every (lipid) class and names them according to the class
#' @keywords internal
#' @return A character vector with the rgb colors named after the vertex-class (e.g. lipid class)
obtain_groupcolors<-function(Annotation){
  group.colors<-c()
  IDxclass<-grep("class$",colnames(Annotation),ignore.case = TRUE)
  Classes<-unique(Annotation[,IDxclass])
  for(class in Classes){
    RgbColor<-getcolor(class,Annotation)
    group.colors<-c(group.colors,RgbColor)
  }
  names(group.colors)<-unique(Annotation[,IDxclass])
  return(group.colors)
}

#' Number lipid features per class
#' @description Internal use. This function counts for every edge feature the number of vertex features (e.g. lipids) per class. It disregards the number of vertex pairs...
#' @keywords internal
#' @importFrom tidyr separate
#' @import dplyr
#' @return A data.frame with the number of vertex features per class and edge feature
countClassPerEdgeFeature<-function(ResTable,treatment="Chow"){
  EdgeFeatures<-unique(ResTable$EdgeFeature)
  ResCountVertexClass<-data.frame(
    EdgeFeature=character(),
    VertexClass=character(),
    Count=numeric(),
    stringsAsFactors = FALSE
  )

  for (edgeFeature in EdgeFeatures){
    FractionEdgeFeature<-ResTable %>% filter(.data$EdgeFeature==edgeFeature) #Count for every EdgeFeature
    FractionEdgeFeature <- FractionEdgeFeature %>% tidyr::separate(.data$VertexClassPair, c("Vertex1", "Vertex2"), "_") #Split the Vertex pairs, get two columns
    TrColumn<-which(colnames(FractionEdgeFeature) == treatment) #Get index of the desired treatment to count Vertexs
    FractionEdgeFeature<-FractionEdgeFeature[FractionEdgeFeature[,TrColumn]>0,] #Get the instances that are not zero
    VertexsPresent<-unique(c(FractionEdgeFeature$Vertex1, FractionEdgeFeature$Vertex2)) #Get the unique Vertexs the gene/EdgeFeature is connected to
    for (vertex in VertexsPresent){ #Loop the Vertexs present
      IdxVertex1<-grep(vertex,FractionEdgeFeature$Vertex1) #From the filtered table (no zero values), get the row indexes of the first Vertex of the Vertex pairs (first Vertex column)
      IdxVertex2<-grep(vertex,FractionEdgeFeature$Vertex2) #From the filtered table (no zero values), get the row indexes of the second Vertex of the Vertex pairs (second Vertex column)
      SumVertex1<-sum(FractionEdgeFeature[IdxVertex1,TrColumn]) #Using the indexes sum the number of times first column
      SumVertex2<-sum(FractionEdgeFeature[IdxVertex2,TrColumn]) #Using the indexes sum the number of times second column
      TotalVertexClassEdgeFeature<-SumVertex1+SumVertex2 #Get the total, number of times that EdgeFeature is involved with that specific Vertex class
      ResCountVertexClass<-rbind(ResCountVertexClass,c(EdgeFeature=edgeFeature,VertexClass=vertex,Count=TotalVertexClassEdgeFeature)) #Add to result table
    }
  }
  colnames(ResCountVertexClass)<-c("EdgeFeature","VertexClass","Count")
  ResCountVertexClass <- ResCountVertexClass %>% filter(.data$Count>0)
  ResCountVertexClass$Count <- as.numeric(ResCountVertexClass$Count) #make sure is numeric
  ResCountVertexClass<-ResCountVertexClass %>% arrange(group_by = .data$EdgeFeature, dplyr::desc(.data$Count))
  return(ResCountVertexClass)
}

#' Split function
#' @description Internal use. Function to split the EdgeFeatures in smaller groups
#' @keywords internal
#' @return A list of character vectors, each vector contains n edge features
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

#' Number Vertex features per class for every shared edge feature
#' @description Internal use. This function creates a barplot depicting the number of vertex features per class for every edge feature. To use this function one has to split first the file (or not if it is small) with the funciton chunk2
#' @keywords internal
#' @import ggrepel
#' @import ggplot2
#' @return ggplot object for a barplot depicting the number of vertex features per class for every edge feature
barplot_VertexsPerEdgeFeature<-function(SplitFile,title="Vertex Features per class",AnnotationWithColors,ggrepelL=TRUE,xlb="Gene",szggrepel=2.5,szTitle=12,szaxisTxt=12,szaxisTitle=12,szLegendText=10,szLegendKey=1){

  EdgeFeature <- Count <- VertexClass <- NULL

  group.colors<-obtain_groupcolors(AnnotationWithColors)
  #Get table that counts how many times the EdgeFeature appears,
  #every time corresponds to a number of a specific Vertex class
  TimesEdgeFeatureTable<-table(SplitFile$EdgeFeature)
  #Create a variable called Id to set order in group barplot
  Id<-c()
  for (i in names(TimesEdgeFeatureTable)){
    Id<-c(Id,1:TimesEdgeFeatureTable[i])
  }
  SplitFile$Id<-as.factor(Id)


  g<-ggplot(SplitFile, aes(x=EdgeFeature, y=as.numeric(Count), fill=VertexClass,group=Id))+
    geom_bar(width = 1,stat="identity",position=position_dodge(0.7)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.title = element_text(size = szLegendText+2),
          legend.text = element_text(size = szLegendText),
          legend.key.size = unit(szLegendKey,"cm"),
          plot.title = element_text(hjust = 0.5,size=szTitle,face="bold"),
          axis.title.y=element_text(angle=90),
          axis.title=element_text(size = szaxisTitle,face="bold", colour = "black"),
          axis.text = element_text(size = szaxisTxt)) +
    scale_fill_manual(values=group.colors,
                      name="Vertex Class",) +
    xlab(xlb)+
    ylab("Count")
  if(ggrepelL){
    g<-g+geom_text_repel(aes(label = VertexClass),
                         size=szggrepel,color="black",
                         min.segment.length = unit(0, 'lines'),
                         nudge_y=0.1,
                         nudge_x=0.1,
                         vjust = -1,
                         force = 0.1,
                         segment.alpha=0.3,
                         max.overlaps=Inf
    )
  }
  g<-g + ggtitle(title)

  return(g)
}

#' Vertex Class profile per edge feature (one treatment)
#' @description This function creates a barplot or barplots showing the number of vertex features per class for every shared edge feature between two treatments
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param Annotation Data frame that includes the rgb colors for every class. The column 'class' (or 'Class') has to be present and also the column 'ColorRgb'
#' @param chunks To avoid a non readable dense plot the results can be spitted in multiple plots
#' @param treat Specify the treatment for which the plot will be created. It should be one of the two treatments in the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param small logical. If only a few edge features are in the input set as TRUE. A single plot will be created
#' @param ggrep logical. If TRUE includes ggrepel labels for every bar
#' @param xlb x-axis label
#' @param onlyTable logical. If TRUE a table is returned instead of a plot
#' @param szTitle Size title
#' @param szaxisTxt Size axis text
#' @param szaxisTitle Size axis title
#' @param ... Other parameters for inner functions, mainly ggplot2 visual parameters
#' @return A list of ggplot objects to create different barplots. The barplots show the number of vertex features per class for every shared
#' edge feature between two treatments. The barplots restrict to one of the compared treatments. An alternative output
#'is a data.frame with the number of vertex features per class and edge feature (onlyTable=TRUE)
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' data(MetColorTable)
#' #Note: No differences in example as all the Output of CoNI was kept
#' getVertexsPerEdgeFeature(CompTreatTable = VertexClassesSharedGenes_HFDvsChow,
#'                          Annotation = MetColorTable,
#'                          chunks = 2,
#'                          treat = "HFD")
#' @export
#' @import ggplot2
#' @import ggrepel
getVertexsPerEdgeFeature<-function(CompTreatTable, Annotation,
                                   chunks = 5, treat=NULL,
                                   small = FALSE,
                                   ggrep = TRUE,
                                   xlb = "Gene",
                                   onlyTable = FALSE,
                                   szTitle = 12,
                                   szaxisTxt = 12,
                                   szaxisTitle = 12, ...){
  if(is.null(treat)){
    stop("Specify treatment")
  }
  EdgeFeatureVertex<-Table<-countClassPerEdgeFeature(CompTreatTable,treatment = treat) #First count per EdgeFeature the number of Vertexs for each class, function above
  # EdgeFeatureVertex$Count <- as.numeric(EdgeFeatureVertex$Count)
  # EdgeFeatureVertex<-EdgeFeatureVertex %>% arrange(group_by = EdgeFeature, dplyr::desc(Count)) #Order the table from high to low per EdgeFeature

  if(onlyTable){
    return(EdgeFeatureVertex)
  }

  #Get the EdgeFeatures of the table
  EdgeFeatures<-unique(EdgeFeatureVertex$EdgeFeature)
  barplots<-list()
  if(!small){
    #Split the results in n pieces so the plots are readable
    SplitIntoPieces<-chunk2(x = EdgeFeatures,chunks)
    for(i in 1:chunks){
      Split<-EdgeFeatureVertex[EdgeFeatureVertex$EdgeFeature %in% SplitIntoPieces[[i]],]
      barplots[[i]]<-barplot_VertexsPerEdgeFeature(Split,title = treat,AnnotationWithColors = Annotation,ggrepelL = ggrep,xlb = xlb,szaxisTxt = szaxisTxt,szaxisTitle=szaxisTitle,...)
      cat(max(Split$Count),"\n")
    }
  }else{
    barplots[[1]]<-barplot_VertexsPerEdgeFeature(EdgeFeatureVertex,title=treat,AnnotationWithColors = Annotation,ggrepelL=ggrep,xlb = xlb,szaxisTxt = szaxisTxt,szaxisTitle=szaxisTitle,szTitle=szTitle,...)
    cat(max(EdgeFeatureVertex$Count),"\n")
  }
  return(barplots)

}

#' Vertex-Class profile per edge feature Side-by-Side (two treatments)
#' @description This function creates a grid of barplots. The barplot of one side depicts the number of class vertex features per edge feature for treatment 1 and the other side the same barplot for treatment 2. Results of both Treatments are side by side for better comparison.
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param Treat1 Name treatment 1 as in table CompTreatTable
#' @param Treat2 Name treatment 2 as in table CompTreatTable
#' @param Annotation Data frame that includes the rgb colors for every class. The column 'class' (or 'Class') has to be present and also the column 'ColorRgb'
#' @param chunks To avoid a non readable dense plot the results can be spitted in multiple plots
#' @param ggrep logical. If TRUE includes ggrepel labels for every bar
#' @param xlb Change the x-axis label
#' @param onlyT logical. If TRUE a table is returned instead of a grid of plots
#' @param small logical. If only a few edge features are in the input set as TRUE. A single plot will be created
#' @param ... Other parameters for inner functions, mainly ggplot2 visual parameters
#' @return A gtable containing side-by-side barplots, one for each treatment, showing the number of vertex features per class for every shared edge feature
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' VCSGs<-VertexClassesSharedGenes_HFDvsChow
#' data(MetColorTable)
#' HFD_vs_Chow_LCP_Gene<-getVertexsPerEdgeFeature_and_Grid(VCSGs,
#'                                                         "HFD","Chow",
#'                                                         Annotation=MetColorTable,
#'                                                         ggrep=FALSE,
#'                                                         small = FALSE,
#'                                                         chunks = 3,
#'                                                         szLegendKey=0.2)
#' plot(HFD_vs_Chow_LCP_Gene)
#' @export
#' @import ggplot2
#' @import ggrepel
#' @importFrom  gridExtra arrangeGrob
#' @importFrom utils capture.output
getVertexsPerEdgeFeature_and_Grid<-function(CompTreatTable,
                                            Treat1, Treat2, Annotation,
                                            chunks = 3, ggrep = TRUE,
                                            xlb = "Edge Feature",
                                            onlyT = FALSE,
                                            small = FALSE,...){
  if(small){
    ylimTreat1<-capture.output(Treat1Plot<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat1,Annotation=Annotation,ggrep = ggrep,xlb = xlb,...))
    ylimTreat2<-capture.output(Treat2Plot<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat2,Annotation=Annotation,ggrep = ggrep, xlb =xlb,...))

    #Get tables
    Treat1Table<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat1,Annotation=Annotation,ggrep = ggrep,xlb = xlb,onlyTable = TRUE)
    Treat2Table<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat2,Annotation=Annotation,ggrep = ggrep, xlb =xlb,onlyTable = TRUE)

    if(onlyT){
      Table_VertexProfile<-Treat1Table %>% full_join(Treat2Table,by= c(colnames(Treat1Table)[1],colnames(Treat1Table)[2]),suffix=c(paste0("_",Treat1),paste0("_",Treat2)))
      Table_VertexProfile <-  Table_VertexProfile %>% mutate_at(vars( starts_with("Count_") ),
                                                                ~if_else( is.na(.), 0, .) )
      return(Table_VertexProfile)
    }

    ylim<-cbind(as.numeric(ylimTreat1),as.numeric(ylimTreat2))
    ylim_max<-apply(ylim,1,max)

    #Assign limits
    plots<-c()
    Tr1<-Treat1Plot[[1]]+ylim(0,ylim_max)
    Tr2<-Treat2Plot[[1]]+ylim(0,ylim_max)
    plots[[1]]<-Tr1
    plots[[2]]<-Tr2
    arrangeGrob(grobs=plots,ncol=2)
  }else{
    #Get ylimits and plots
    ylimTreat1<-capture.output(Treat1Plots<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat1,Annotation=Annotation,ggrep = ggrep,xlb = xlb,...))
    ylimTreat2<-capture.output(Treat2Plots<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,small=small,treat = Treat2,Annotation=Annotation,ggrep = ggrep, xlb =xlb,...))

    #Get tables
    Treat1Table<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,treat = Treat1,Annotation=Annotation,ggrep = ggrep,xlb = xlb,onlyTable = TRUE)
    Treat2Table<-getVertexsPerEdgeFeature(CompTreatTable,chunks = chunks,treat = Treat2,Annotation=Annotation,ggrep = ggrep, xlb =xlb,onlyTable = TRUE)

    if(onlyT){
      Table_VertexProfile<-Treat1Table %>% full_join(Treat2Table,by= c(colnames(Treat1Table)[1],colnames(Treat1Table)[2]),suffix=c(paste0("_",Treat1),paste0("_",Treat2)))
      Table_VertexProfile <-  Table_VertexProfile %>% mutate_at(vars( starts_with("Count_") ),
                                                                ~if_else( is.na(.), 0, .) )
      return(Table_VertexProfile)
    }


    ylim<-cbind(as.numeric(ylimTreat1),as.numeric(ylimTreat2))
    ylim_max<-apply(ylim,1,max)

    #Assign limits
    plots<-c()
    i<-0 #As there are 2 plots per chunk, I need to loop two numbers simultaneously. One to assign the limits and one to add to the list
    for(j in seq(1,chunks*2,2)){
      i<-i+1
      Tr1<-Treat1Plots[[i]]+ylim(0,ylim_max[i])
      Tr2<-Treat2Plots[[i]]+ylim(0,ylim_max[i])
      plots[[j]]<-Tr1
      plots[[j+1]]<-Tr2
    }
    return(arrangeGrob(grobs=plots,ncol=2)) #instead of grid.arrange that plots to console to avoid error related to size of screen
  }
}

#'Stacked Global Barplot (One treatment)
#' @description This function will create a stacked barplot from the output of Compare_VertexClasses_sharedEdgeFeatures using all shared Edge Features (e.g., genes) between two treatments.
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param treat Name of treatment to display. It should match the column name in the output of Compare_VertexClasses_sharedEdgeFeatures
#' @param xlb Name for x-axis
#' @param ylb Name for y-axis
#' @param max_pairsLegend If number of Edge Features >= max_pairsLegend, display number of Edge Features as label with ggrepel
#' @param mx.overlaps Max number of overlaps ggrepel
#' @param szggrepel Size ggrepel labels
#' @param force Repelling force for ggrepel labels
#' @param szTitle Size Title
#' @param szaxisTxt Size axis text
#' @param szaxisTitle Size axis titles
#' @param ylim Optional y-limits of the plot
#' @import ggplot2
#' @import ggrepel
#' @return A ggplot object to create a stacked barplot. The stacked barplot shows the vertex-class pairs profile of all shared edge features but restricted to a single treatment. Every bar consists of multiple edge features (stacked) that are represented with different colors
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' create_stackedGlobalBarplot_perTreatment(CompTreatTable = VertexClassesSharedGenes_HFDvsChow,
#'                                          treat = "HFD",
#'                                          max_pairsLegend = 9,
#'                                          xlb = "Metabolite-class-pairs")
#' @export
#' @importFrom tidyr gather
#' @importFrom forcats fct_relevel
#' @importFrom gplots col2hex
#' @importFrom rlang .data
#' @importFrom tidyselect vars_select_helpers
create_stackedGlobalBarplot_perTreatment<-function(CompTreatTable,
                                                   treat=NULL,
                                                   xlb="Vertex-Class Pairs",
                                                   ylb="Number of pairs",
                                                   max_pairsLegend = 2,
                                                   mx.overlaps = Inf,
                                                   szggrepel=6,
                                                   force=0.1,
                                                   szTitle=12, szaxisTxt=12, szaxisTitle=12,
                                                   ylim=NULL){

  treatment <- number_pairs <- VertexClassPair <- EdgeFeature <- label <- NULL

  CompTreatTable[,ncol(CompTreatTable)]<-as.numeric(gsub("Pair Class Missing",0,CompTreatTable[,ncol(CompTreatTable)]))
  CompTreatTable[,ncol(CompTreatTable)-1]<-as.numeric(gsub("Pair Class Missing",0,CompTreatTable[,ncol(CompTreatTable)-1]))

  if(is.null(treat)){
    print("Provide treatment to filter data e.g., treat='HFD'")
  }else{
    #Global view of the relationship of edge features with metabolite classes
    Stacked<-CompTreatTable %>% group_by(.data$VertexClassPair,.data$EdgeFeature) %>% summarise(across(tidyselect::vars_select_helpers$where(is.numeric),sum))
    StackedFBarplot <- tidyr::gather(Stacked, treatment, number_pairs,c(-VertexClassPair,-EdgeFeature))
    StackedFBarplotTreat<-StackedFBarplot %>% filter(.data$treatment==treat)

    ColorTable<-data.frame(
      EdgeFeature=StackedFBarplotTreat$EdgeFeature,
      ColorGroup=labels2colors_2(as.numeric(as.factor(StackedFBarplotTreat$EdgeFeature))))
    ColorGroupRgb<-col2hex(ColorTable$ColorGroup)
    names(ColorGroupRgb) = ColorTable$EdgeFeature
    #StackedFBarplotTreat$EdgeFeature <- factor(StackedFBarplotTreat$EdgeFeature, levels = ColorTable$EdgeFeature)

    countmaxperLipidClass<-StackedFBarplotTreat %>% group_by(.data$VertexClassPair) %>% summarise(number_pairs = sum(.data$number_pairs))
    cat(max(countmaxperLipidClass$number_pairs))



    TreatStacked <- ggplot(StackedFBarplotTreat, aes(x =.data$VertexClassPair , y = .data$number_pairs, fill=.data$EdgeFeature))+
      geom_bar(stat="identity") +
      theme(legend.position = "none") + #Remove legend
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.title = element_text(hjust = 0.5,size=szTitle,face="bold"),
            axis.title=element_text(size = szaxisTitle,face="bold", colour = "black"),
            axis.text = element_text(size = szaxisTxt))+
      geom_text_repel(data = StackedFBarplotTreat %>%
                      mutate(label = ifelse(number_pairs>=max_pairsLegend,paste0(EdgeFeature,": ",number_pairs),"")),
                      aes(label=label),
                      size = szggrepel,
                      position = position_stack(vjust = 0.5),
                      color="black",
                      alpha=0.6,
                      force= force,
                      max.overlaps=mx.overlaps,
                      min.segment.length = unit(0, 'lines')) +
      scale_colour_manual(values = ColorGroupRgb) +
      scale_fill_manual(values=ColorGroupRgb)+
      xlab(xlb)+
      ylab(ylb)
    if(!is.null(ylim)){
      TreatStacked<-TreatStacked+ ylim(c(0,ylim))
    }
    TreatStacked<-TreatStacked + ggtitle(treat)
    return(TreatStacked)
  }

}

#' Stacked Global Barplot Side-by-side (two treatments)
#' @description This function will create a stacked barplot from the output of Compare_VertexClasses_sharedEdgeFeatures using all shared Edge Features (e.g., genes) between two treatments. Results of both Treatments are side by side for better comparison.
#' @param CompTreatTable Output of Compare_VertexClasses_sharedEdgeFeatures
#' @param Treat1 Name treatment 1 as in table CompTreatTable
#' @param Treat2 Name treatment 2 as in table CompTreatTable
#' @param ggrep logical. If TRUE includes ggrepel labels for every bar
#' @param max_pairsLegend If number of Edge Features >= max_pairsLegend, display number of Edge Features as ggrepel label
#' @param mx.overlaps Max number of overlaps ggrepel
#' @param szggrepel Size ggrepel labels
#' @param force Repelling force for ggrepel labels
#' @param xlb Name for x-axis
#' @param ... Other parameters for inner functions, mainly ggplot2 visual parameters
#' @return A gtable containing stacked barplots. The barplots show the vertex-class pairs profile of all shared edge features between two treatments (one bar plot per treatment). Every bar consists of multiple edge features that are depicted with different colors
#' @examples
#' data(VertexClassesSharedGenes_HFDvsChow)
#' VCSGs<-VertexClassesSharedGenes_HFDvsChow
#' HFD_vs_Chow_stackedBarplot<-getstackedGlobalBarplot_and_Grid(VCSGs,
#'                                                              Treat1 = "HFD",
#'                                                              Treat2 = "Chow",
#'                                                              xlb = "Metabolite-class-pairs",
#'                                                              max_pairsLegend=9)
#' plot(HFD_vs_Chow_stackedBarplot)
#' @importFrom utils capture.output
#' @export
getstackedGlobalBarplot_and_Grid<-function(CompTreatTable, Treat1, Treat2,
                                           ggrep=TRUE, max_pairsLegend = 1, force = 0.1, mx.overlaps=Inf, szggrepel=6,
                                           xlb = "Vertex-Class Pairs",...){
  #Get ylimits and plots
  ylimTreat1<-capture.output(Treat1Plot<-create_stackedGlobalBarplot_perTreatment(CompTreatTable = CompTreatTable,
                                                                                  treat = Treat1,
                                                                                  max_pairsLegend = max_pairsLegend,
                                                                                  force = force,
                                                                                  xlb = xlb,
                                                                                  mx.overlaps = mx.overlaps,
                                                                                  szggrepel = szggrepel, ...))
  ylimTreat2<-capture.output(Treat2Plot<-create_stackedGlobalBarplot_perTreatment(CompTreatTable = CompTreatTable,
                                                                                  treat = Treat2,
                                                                                  max_pairsLegend = max_pairsLegend,
                                                                                  force = force,
                                                                                  xlb = xlb,
                                                                                  mx.overlaps = mx.overlaps,
                                                                                  szggrepel = szggrepel, ...))

  ylim<-cbind(as.numeric(ylimTreat1),as.numeric(ylimTreat2))
  ylim_max<-apply(ylim,1,max)

  #Assign limits
  plots<-c()
  Tr1<-Treat1Plot+ylim(0,ylim_max[1])
  Tr2<-Treat2Plot+ylim(0,ylim_max[1])
  plots[[1]]<-Tr1
  plots[[2]]<-Tr2
  arrangeGrob(grobs=plots,ncol=2) #instead of grid.arrange that plots to console to avoid error related to size of screen
}

#' Bipartite Table
#' @description Internal use. This function creates a table that is used to create a simple bipartite network
#' @keywords internal
#' @return A matrix containing two columns, one for vertexes and one for edge features. The matching happens if they are adjacent to one another
createBipartiteTable<-function(CoNINetworkTable){
  resultFinal<-list()
  for(n in 1:nrow(CoNINetworkTable)){
    vertex1<-CoNINetworkTable[n,1] #Get vertex 1
    vertex2<-CoNINetworkTable[n,2] #Get vertex 2
    Links<-CoNINetworkTable[n,4] #Get Links
    LinkL<-unlist(strsplit(Links,";")) #Split to get list
    # PcorValuesL_Link1vsDriver<-unlist(strsplit(CoNINetworkTable[n,9],";")) #Get pcor values and split to get list
    # PcorValuesL_Link2vsDriver<-unlist(strsplit(CoNINetworkTable[n,10],";")) #Get pcor values and split to get list
    rows_complete<-list()
    #pcor=numeric()
    for (i in 1:length(LinkL)){
      rowvertex1<-c(from=vertex1,to=LinkL[i]) #pcor=PcorValuesL_Link1vsDriver[i]
      rowvertex2<-c(from=vertex2,to=LinkL[i]) #pcor=PcorValuesL_Link2vsDriver[i]
      rows<-rbind(rowvertex1,rowvertex2)
      rows_complete[[i]]<-rows
    }
    rows_complete<-do.call(rbind,rows_complete)
    # result = foreach (n = 1:nrow(CoNINetworkTable), .combine=rbind) %do% {
    # }
    resultFinal[[n]]<-rows_complete
  }
  resultFinal<-do.call(rbind,resultFinal)
  return(resultFinal)
}

#' Bipartite Network
#' @description This function creates a simple bipartite graph, it shows the driver and linker features as nodes.
#' @param TableNetwork TableForNetwork_file (file generated by CoNI) with path
#' @param colorVertexTable Table specifying the colors for the vertex features. The first column should contain the names matching the features of the vertex Data and another column should specify the colors (column name: Colors).
#' @param incidenceMatrix logical. If TRUE it returns a hypergraph incidence matrix instead of a bipartite graph
#' @return An igraph object for a bipartite graph or a hypergraph incidence matrix to represent ResultsCoNI. Basic network statistics are included in the bipartite graph. See generate_network function for details or consult the igraph package
#' @examples
#' #See vignette for an example
#' @export
#' @importFrom gplots col2hex
#' @importFrom igraph graph_from_data_frame simplify V E degree hub_score transitivity closeness betweenness eigen_centrality centr_betw centr_clo centr_degree edge_betweenness write_graph as_ids get.incidence
createBipartiteGraph<-function(TableNetwork,colorVertexTable,incidenceMatrix=FALSE){
  TableNetwork<-read.csv(TableNetwork)
  #Create bipartite table
  bipartiteTable<-createBipartiteTable(TableNetwork)
  #Remove redundancy
  bipartiteTable<-unique(bipartiteTable)
  bipartiteTable<-as.data.frame(bipartiteTable)

  #Check nodes if there are identical names in vertex and linked features
  LinkedFeaturesIdentical<-bipartiteTable$to[bipartiteTable$to %in% bipartiteTable$from]
  if(length(LinkedFeaturesIdentical)>0){
    LinkedFeaturesIdentical<-unique(LinkedFeaturesIdentical)
    for(feature in LinkedFeaturesIdentical){
      bipartiteTable$to[bipartiteTable$to==feature]<-paste0(feature,"_linkedF")
    }
  }

  #Create graph
  cnodes <- data.frame("Name"=unique(c(as.character(bipartiteTable$from),as.character(bipartiteTable$to))),stringsAsFactors=FALSE)#Get the nodes vertexFeature-EdgeFeature
  #Assign colors to nodes
  m <- merge(cnodes,colorVertexTable,by.x="Name",by.y=colnames(colorVertexTable)[1],all.x=TRUE)
  #Assign grey color to Edges
  m<- m %>% mutate(type=ifelse(is.na(m[,2]),"EdgeFeature","VertexFeature"))#Might be problematic but minimum annotation file should contain three columns, vertex-feature,color and colorRgb
  #m[is.na(m[,2]),2]<-"EdgeFeature"
  idx_colorColumn<-grep("color$",colnames(m),ignore.case = TRUE)
  m[is.na(m[,idx_colorColumn]),idx_colorColumn]<-"grey"
  m$ColorRgb<-col2hex(m[,idx_colorColumn])
  #Create graph
  netd <- graph_from_data_frame(d=bipartiteTable, vertices=m, directed=FALSE)
  netd <- simplify(netd,remove.multiple=FALSE)

  #Bipartite option
  #bipartite_mapping(netd)$type this function is giving me problems
  igraph::V(netd)$type <- ifelse(igraph::as_ids(V(netd)) %in% bipartiteTable$from,TRUE,FALSE)
  igraph::V(netd)$shape <- ifelse(V(netd)$type, "circle", "square")

  #Add network stats to bipartite graph... at the moment ignoring that it is a bipartite graph
  #stats are as as if it were a simple network
  igraph::V(netd)$degree <- degree(netd, mode="all")
  igraph::V(netd)$hub_score <- hub_score(netd, weights=NA)$vector
  igraph::V(netd)$transitivity <- transitivity(netd, type="local")
  igraph::V(netd)$closeness <- closeness(netd, mode="all", weights=NA)
  igraph::V(netd)$betweenness <- betweenness(netd, directed=FALSE, weights=NA)
  igraph::V(netd)$eigen_centrality <- eigen_centrality(netd, directed=FALSE, weights=NA)$vector
  igraph::V(netd)$centralized_betweenness <- centr_betw(netd, directed=FALSE, normalized=TRUE)$res
  igraph::V(netd)$centralized_closeness <- centr_clo(netd, mode="all", normalized=TRUE)$res
  igraph::V(netd)$centralized_degree <- centr_degree(netd, mode="all", normalized=TRUE)$res
  #V(netd)$membership_community_edgeBetweenes<-cluster_edge_betweenness(netd,directed = F)$membership

  #Add edge betweeness
  igraph::E(netd)$betweeness <- edge_betweenness(netd, directed=FALSE, weights=NA)

  if(incidenceMatrix){
    incidenceM<-get.incidence(netd)
    return(as.matrix(incidenceM))
  }else{
    return(netd)
  }
}

#' Network Statistics
#' @description This function calculates simple network statistics and returns them as a dataframe
#' @param Network An Igraph network
#' @return Returns a data.frame with nine rows with the following network statistics:
#'  \itemize{
##'  \item{"net_avg_pathL"}{Shortest paths between vertices}
##'  \item{"net_edge_density"}{Graph density, ratio of the number of edges and the number of possible edges}
##'  \item{"net_transitivity"}{Probability that the adjacent vertices of a vertex are connected}
##'  \item{"net_diameter"}{Length of the longest geodesic}
##'  \item{"net_nodes_first_path_diameter"}{The nodes along the first found path with the length of diameter}
##'  \item{"net_eigenvalue"}{The eigenvalue corresponding to the centrality scores.}
##'  \item{"net_centralized_betweenessIdx"}{The graph level centrality index after centralizing the graph according to the betweenness of vertices}
##'  \item{"net_centralized_closenessIdx"}{The graph level centrality index after centralizing the graph according to the closeness of vertices}
##'  \item{"net_centralized_degreeIdx"}{The graph level centrality index after centralizing the graph according to the degrees of vertices}
##' }
#' For more information on the statistics consult the igraph package.
#' @examples
#' #Load color nodes table
#' data(MetColorTable)
#' #Assign colors according to "Class" column
#' MetColorTable<-assign_colorsAnnotation(MetColorTable)
#' #Load CoNI results
#' data(CoNIResultsHFDToy)
#'
#' #Generate Network
#' HFDNetwork<-generate_network(ResultsCoNI = CoNIResultsHFDToy,
#'                              colorVertexNetwork = TRUE,
#'                              colorVertexTable = MetColorTable,
#'                              outputDir = "./",
#'                              outputFileName = "HFD",
#'                              saveFiles = FALSE)
#' NetStats(HFDNetwork)
#' @importFrom  tibble rownames_to_column
#' @importFrom igraph mean_distance edge_density transitivity diameter get_diameter eigen_centrality centr_betw centr_clo centr_degree
#' @import dplyr
#' @export
NetStats<-function(Network){
  NetworkStatsTable<-data.frame(Value=t(data.frame(
    net_avg_pathL=mean_distance(Network, directed=F),
    net_edge_density=edge_density(Network, loops=F),
    net_transitivity=transitivity(Network, type="global"),
    net_diameter=diameter(Network, directed=F, weights=NA),
    net_nodes_first_path_diameter= paste(names(get_diameter(Network, directed=TRUE)),collapse=","),#returns the nodes along the first found path of that distance
    net_eigenvalue=eigen_centrality(Network, directed=FALSE, weights=NA)$value,
    net_centralized_betweenessIdx=centr_betw(Network, directed=F, normalized=TRUE)$centralization,
    net_centralized_closenessIdx=centr_clo(Network, mode="all", normalized=TRUE)$centralization,
    net_centralized_degreeIdx=centr_degree(Network, mode="all", normalized=TRUE)$centralization
    #net_community__modularity_edgeBetweenes=modularity(cluster_edge_betweenness(Network,directed = F))
  )))
  NetworkStatsTable<-NetworkStatsTable %>% rownames_to_column("Network_statistic")
  return(NetworkStatsTable)
}

#' Get vertexes for edge feature
#' @keywords internal
#' @return A character vector with the vertexes connected to a given edge feature
getvertexes_edgeFeature<-function(edgeFeature,CoNIResults){
  Tvertexes<-CoNIResults[CoNIResults$Feature_edgeD==edgeFeature,c(1,2),drop=FALSE]
  vertexes<-unique(unlist(c(Tvertexes[,1],Tvertexes[,2])))
  return(vertexes)
}

#' Correlation vs Partial correlation
#' @description This function fits two linear models on standardize data and plots the results. It generates a scatter plot with two regression lines, where the slopes correspond to the correlation and partial correlation coefficients (blue for cor and red for pcor)
#' @param ResultsCoNI The significant results generated by CoNI
#' @param edgeFeature The edge feature to explore e.g. Fabp2 (for a gene)
#' @param vertexD Vertex data that was given as input to CoNI
#' @param edgeD Edge data that was given as input to CoNI
#' @param vertexFeatures The vertex features to include as a list. If not specified all metabolites available in combination with the edgeFeature will be used
#' @param outputDir Output directory with path
#' @param label_edgeFeature Name for plot title e.g. Gene or Protein
#' @param plot_to_screen logical. If TRUE plots will be outputted to the plotting screen
#' @param fname File name to save the plots
#' @param height height of the plotting area for the saved file
#' @param width width of the plotting are for the saved file
#' @param saveFiles logical. If FALSE plot is not saved to disk
#' @import ggrepel
#' @import ggplot2
#' @importFrom stats lm as.formula resid
#' @importFrom rlang .data
#' @examples
#' #Load gene expression - Toy dataset of two treatments
#' data(GeneExpToy)
#' #Samples in rows and genes in columns
#' GeneExp <- as.data.frame(t(GeneExpToy))
#' hfd_gene <- GeneExp[1:8,] #high fat diet
#' chow_gene<- GeneExp[9:nrow(GeneExp),] #chow diet
#'
#' #Load metabolite expression - Toy dataset of two treatments
#' data(MetaboExpToy)
#' MetaboExp <- MetaboExpToy
#' hfd_metabo <- MetaboExp[11:18,] #high fat diet
#' chow_metabo <- MetaboExp[1:10,] #chow diet
#'
#' #Match row names both data sets
#' rownames(hfd_metabo)<-rownames(hfd_gene)
#' rownames(chow_metabo)<-rownames(chow_gene)
#'
#' #Load CoNI results
#' data(CoNIResultsHFDToy)
#'
#' plotPcorvsCor(ResultsCoNI = CoNIResultsHFDToy,
#'               edgeFeature = "Arfrp1",
#'               vertexFeatures = c("PC.ae.C40.2", "SM..OH..C22.1"),
#'               vertexD = hfd_metabo,
#'               edgeD = hfd_gene,
#'               label_edgeFeature = "Gene",
#'               plot_to_screen = TRUE,
#'               height = 10,
#'               saveFiles = FALSE)
#' @return Returns a ggplot object for a scatter plot with two regression lines.
#' The blue line is the regression of the vertex features, and the red line is the regression
#' of the resulting residuals after regressing each vertex feature with the edge feature.
#' The slope of the blue line corresponds to the pearson correlation coefficient and the slope of the red line
#' to the partial correlation coefficient
#' @export
plotPcorvsCor<-function(ResultsCoNI,
                        edgeFeature,
                        vertexD, edgeD,
                        vertexFeatures=NULL,
                        outputDir="./",
                        fname,
                        label_edgeFeature="Edge Feature",
                        plot_to_screen=TRUE,
                        height=10,width=8,
                        saveFiles=FALSE){

  fname<-paste0(outputDir,edgeFeature,".pdf")
  ResultsCoNIfull<-ResultsCoNI %>% filter(.data$Feature_edgeD==edgeFeature)
  if(!is.null(vertexFeatures)){
    llF1<-sapply(ResultsCoNIfull$Feature_1_vertexD,function(feature){feature %in% vertexFeatures})
    ResultsCoNIfull<-ResultsCoNIfull[llF1,]
    llF2<-sapply(ResultsCoNIfull$Feature_2_vertexD,function(feature){feature %in% vertexFeatures})
    ResultsCoNIfull<-ResultsCoNIfull[llF2,]
  }

  plots<-list()
  for(i in 1:nrow(ResultsCoNIfull)){
    ResultsCoNI<-ResultsCoNIfull[i, ,drop=FALSE]
    #Example Fabp2
    vertexes_edgeFeature<-getvertexes_edgeFeature(edgeFeature = edgeFeature,CoNIResults = ResultsCoNI)
    M1<-vertexes_edgeFeature[1]
    M2<-vertexes_edgeFeature[2]

    edgeFeature_vertex_Expression<-as.data.frame(cbind(edgeD[,edgeFeature,drop=FALSE],vertexD[,vertexes_edgeFeature,drop=FALSE]))
    fN<-ncol(edgeFeature_vertex_Expression)-1
    edgeFeature_vertex_Expression[,1:fN]<-apply(edgeFeature_vertex_Expression[,1:fN],2,as.numeric)

    #Linear model vertex1 and vertex2
    fitM1M2<-lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[2], "~",colnames(edgeFeature_vertex_Expression)[3])), data=edgeFeature_vertex_Expression)
    summary(fitM1M2)

    #Residuals vertex 1 and vertex 2
    eM1M2<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[2], "~",colnames(edgeFeature_vertex_Expression)[3])), data=edgeFeature_vertex_Expression))
    eM2M1<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[3], "~",colnames(edgeFeature_vertex_Expression)[2])), data=edgeFeature_vertex_Expression))

    #Residuals vertex 1 and edgeFeature 1
    eM1G1<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[2], "~",colnames(edgeFeature_vertex_Expression)[1])), data=edgeFeature_vertex_Expression))
    eG1M1<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[1], "~",colnames(edgeFeature_vertex_Expression)[2])), data=edgeFeature_vertex_Expression))

    #Residuals vertex 2 and edgeFeature 1
    eM2G1<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[3], "~",colnames(edgeFeature_vertex_Expression)[1])), data=edgeFeature_vertex_Expression))
    eG1M2<-resid(lm(as.formula(paste0(colnames(edgeFeature_vertex_Expression)[1], "~",colnames(edgeFeature_vertex_Expression)[3])), data=edgeFeature_vertex_Expression))

    ResMatrix<-as.data.frame(cbind(eM1M2=eM1M2,eM1G1=eM1G1,eM2G1=eM2G1,eG1M1=eG1M1,eM2M1=eM2M1,eG1M2=eG1M2))

    #Scale
    NewDF<-as.data.frame(cbind(vertex1=edgeFeature_vertex_Expression[[2]],vertex2=edgeFeature_vertex_Expression[[3]],eM1G1=ResMatrix$eM1G1,eM2G1=ResMatrix$eM2G1))
    NewDF<-as.data.frame(scale(NewDF))

    plots[[i]]<-ggplot(NewDF,aes(.data$vertex1,.data$vertex2)) +
      geom_point()+
      stat_smooth(method="lm",se=FALSE)+
      geom_point(data = NewDF,aes(eM1G1,eM2G1),color="red")+
      stat_smooth(data = NewDF,aes(eM1G1,eM2G1),color="red",method="lm",se=FALSE)+
      xlab(M1)+
      ylab(M2)+
      ggtitle(paste0(label_edgeFeature,": ",edgeFeature))+
      theme(plot.title = element_text(size=14,color="red", hjust = 0.5,face="bold.italic"))
  }
  if(plot_to_screen){
    sapply(plots,plot)
  }
  if(saveFiles){
    if(length(plots)>1){
      plots_arrange<-arrangeGrob(grobs=plots,ncol=2)
      ggsave(filename=fname, plot=plots_arrange, width=width, height=height)
    }else{
      ggsave(filename=fname, plot=plots[[1]], width=6, height=4)
    }
  }
}

