\name{CoNI}
\alias{CoNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Correlation guided Network Integration
%%  ~~function to do ... ~~
}
\description{This is the main function of Correlation guided Network Integration (CoNI). It takes as input data from two sources (e.g. gene expression and metabolite expression). It calculates all pariwise correlations of the elements of the second data input, and the partial correlations of these pairwise elements with respect to the elements of the first data input. Both data inputs are pre-filtered to include only those elements that significantly correlate. The first data input (typically the largest data input) is also prefiltered to keep only elements with low variance (var<0.5). A modified steiger test is used to identify significant changes between the correaltion and partial correlation values. Results can be visually represented in a Network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CoNI(data1,data2,data2Name="data2",padjustdata2=TRUE,
      correlateDFs=TRUE,splitData1=TRUE,split_number=2,outputDir="./OutputDir/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{Object to use as first data input.
}
  \item{data2}{Object to use as second data input.
}
  \item{data2Name}{File name to save significant correlation results of object data2. This file will be read if function is called again with the same file names.
}
  \item{padjustdata2}{logical. If TRUE adjusted p-values will be used to filter second data input.
}
  \item{correlateDFs}{logical. If TRUE the elements that significantly correlate of data2 will be correlated with the elements of data1. Only the elements that significantly correlate are kept.
}
  \item{splitData1}{logical. If TRUE data1 (tipically the largest data input) will be split into n subsets for the computation. Keep as TRUE unless data input is small (FALSE option under construction).
}
  \item{split_number}{Number of parts to split the elements of data1.
}
  \item{outputDir}{Output Directory where (split) results files are stored.
}
  \item{iteration_start}{Iteration start for CoNI. Useful when run is interrrupted.
}
  \item{wait_iteration}{Time lag between iterations.
}
  \item{numCores}{Cores assigned for parallelization.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Load gene expression
data(GeneExp)
#Samples in rows and genes in columns
GeneExp <- as.data.frame(t(GeneExp))
hfd_gene <- GeneExp[1:8,]

#Load metabolite expression
data(MetaboExp)
hfd_metabo <- MetaboExp[11:18,]

#Match row names both data sets
rownames(hfd_metabo)<-rownames(hfd_gene)

#Run CoNI
CoNIResults <- CoNI(hfd_gene,hfd_metabo,
                    data1name="HFD_genes",data2name = "HFD_metabolites",
                    padjustData2 = FALSE,
                    split_number = 2)

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
