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


