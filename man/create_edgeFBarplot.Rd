% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoNI_functions.R
\name{create_edgeFBarplot}
\alias{create_edgeFBarplot}
\title{Vertex-class pairs profile of one shared feature between two treatments}
\usage{
create_edgeFBarplot(
  CompTreatTable,
  edgeF,
  treat1 = "Treatment1",
  treat2 = "Treatment2",
  factorOrder = NULL,
  col1 = "red",
  col2 = "blue",
  EdgeFeatureType = "Edge Feature",
  xlb = "Vertex-Class Pairs",
  ylb = "Number of pairs",
  szaxisTxt = 12,
  szaxisTitle = 12
)
}
\arguments{
\item{CompTreatTable}{Output of Compare_VertexClasses_sharedEdgeFeatures}

\item{edgeF}{Edge feature present in output of Compare_VertexClasses_sharedEdgeFeatures}

\item{treat1}{Name of treatment one, default Treatment1. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures}

\item{treat2}{Name of treatment one, default Treatment2. It should match the column names of the output of Compare_VertexClasses_sharedEdgeFeatures}

\item{factorOrder}{A list specifying the order of the treatments.}

\item{col1}{Color for Treatment 1}

\item{col2}{Color for Treatment 2}

\item{EdgeFeatureType}{Type of Edge Feature (e.g., Gene)}

\item{xlb}{Name for x-axis}

\item{ylb}{Name for the y-axis}

\item{szaxisTxt}{Size axis text}

\item{szaxisTitle}{Size axis titles}
}
\value{
A barplot showing the vertex-class pairs profile of the shared edge feature between treatments
}
\description{
This function will create a barplot from the output of Compare_VertexClasses_sharedEdgeFeatures for a specific shared Edge Feature (e.g., a shared gene).
}
\examples{
create_edgeFBarplot(CompTreatTable = LCP_HFDvsChow,
                    edgeF = "Otulin",
                    treat1 = "HFD",
                    treat2 = "Chow",
                    factorOrder = c("HFD","Chow"),
                    EdgeFeatureType = "Gene")
}