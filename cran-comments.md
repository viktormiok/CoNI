## Test enviroments
* local x86_64-apple-darwin17.0 (64-bit), R 4.0.3
* win-builder release x86_64-w64-mingw32 (64-bit), R 4.1.1
* win-builder devel x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2021-09-23 r80951)


## R CMD check results

24. September 2021
There were no ERRORs or WARNINGs.

There was 1 NOTE in local:
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘José Manuel Monroy Kuhn <nolozz@gmail.com>’

New submission

There were 3 NOTEs in win-builder release:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'José Manuel Monroy Kuhn <nolozz@gmail.com>'

New submission

** running examples for arch 'i386' ... [47s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
     user system elapsed
CoNI 0.15   0.05   12.09
** running examples for arch 'x64' ... [50s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
     user system elapsed
CoNI 0.19   0.04   12.05


There were 2 NOTEs in win-builder devel:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'JosÈ Manuel Monroy Kuhn <nolozz@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  CoNI (3:48)
  hypergraph (17:170)
  omics (17:39)

* checking examples ... [51s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
     user system elapsed
CoNI 0.17   0.05   12.99


25. September 2021
I wrapped CoNI example in \donttest{} to pass the test
Only one note remains

Comments from CRAN:

If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or 
'F' as vector names.), e.g.:
   man/createBipartiteGraph.Rd
   man/getstackedGlobalBarplot_and_Grid.Rd
   man/getVertexsPerEdgeFeature_and_Grid.Rd
   man/getVertexsPerEdgeFeature.Rd
   man/plotPcorvsCor.Rd

Please add \value to .Rd files regarding exported methods and explain 
the functions results in the documentation. Please write about the 
structure of the output (class) and also what the output means. (If a 
function does not return a value, please document that too, e.g. 
\value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 23 .Rd files, e.g.:
      assign_colorsAnnotation.Rd: \value
      barplot_VertexsPerEdgeFeature.Rd: \value
      check_outputDir.Rd: \value
      check_previous.Rd: \value
      checkInputParameters.Rd: \value
      chunk2.Rd: \value
      
Please ensure that you do not use more than 2 cores in your examples, 
vignettes, etc.

28. September 2021
There was 1 NOTE in local:
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘José Manuel Monroy Kuhn <nolozz@gmail.com>’

New submission

There was 1 NOTE in win-builder release:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jos� Manuel Monroy Kuhn <nolozz@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  CoNI (3:49)
  al (17:242)
  et (17:239)
  hypergraph (17:170)
  omics (17:39)
  
There was 1 NOTE in win-builder devel:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jos� Manuel Monroy Kuhn <nolozz@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  CoNI (3:49)
  al (17:242)
  et (17:239)
  hypergraph (17:170)
  omics (17:39)



