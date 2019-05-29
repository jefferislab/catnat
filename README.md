[![Travis-CI Build Status](https://api.travis-ci.org/alexanderbates/catnat.svg?branch=master)](https://travis-ci.org/jefferislab/catnat)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](http://jefferislab.github.io/catnat/reference/)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

# catnat
R Package for use with [rcatmaid](https://github.com/jefferis/rcatmaid) and [nat](https://github.com/jefferis/rcatmaid). **catnat** provides some higher level analysis function for, for example, clustering synapses within a neuron's tree structure, clustering together neurons by synapse position in 3D space and splitting a neuron into different compartments (e.g. axon-dendrite-primary neurite) and visualising these splits and clusters. The packages remains heavily in development.

## What's in the package currently?
```r
# install
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferislab/catnat")

# use
library(catnat)

# some useful functions
?seesplit3d
?flow.centrality
?get.synapses
?cluster.by.synapses
```

Note: Installation depends on the [devtools](http://CRAN.R-project.org/package=devtools) package. 
Windows users may need to install [Rtools](http://www.murdoch-sutherland.com/Rtools/).

## Acknowledgements

**catnat** depends on [rcatmaid](https://github.com/jefferis/rcatmaid) and [nat](https://github.com/jefferis/nat),
which were developed principally by Greg Jefferis.

rcatmaid is based on python code presently visible at:

* https://github.com/catmaid/CATMAID/blob/master/scripts/remote/access.py
* https://github.com/catmaid/CATMAID/blob/master/django/applications/catmaid/urls.py
* https://github.com/schlegelp/CATMAID-to-Blender/blob/master/CATMAIDImport.py

by Albert Cardona and Philipp Schlegel. Released under the GPL-3 license.
