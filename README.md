[![Travis-CI Build Status](https://api.travis-ci.org/alexanderbates/catnat.svg?branch=master)](https://api.travis-ci.org/alexanderbates/catnat)

# catnat
R Package for use with [rcatmaid](https://github.com/jefferis/rcatmaid) and [nat](https://github.com/jefferis/rcatmaid). Proivdes some higher level analysis function for, for example, clustering synapses within a neuron's tree structure, clustering together neurons by synapse position in 3D space and splitting a neuron into different compartments (e.g. axon-dendrite-primary neurite) and visualising these splits and clusters. In development.

## What's in the package currently?
```r
# install
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferis/nat")
devtools::install_github("jefferis/rcatmaid")
devtools::install_github("alexanderbates/catnat")
library(catnat)
?seesplit3d
?flow.centrality
?get.synapses
?cluster.by.synapses
```

Note: Windows users need [Rtools](http://www.murdoch-sutherland.com/Rtools/) and
[devtools](http://CRAN.R-project.org/package=devtools) to install this way.

## Acknowledgements

rcatmaid and nat were developed principally by Greg Jefferis.

rcatmaid is based on python code presently visible at:

* https://github.com/catmaid/CATMAID/blob/master/scripts/remote/access.py
* https://github.com/catmaid/CATMAID/blob/master/django/applications/catmaid/urls.py
* https://github.com/schlegelp/CATMAID-to-Blender/blob/master/CATMAIDImport.py

by Albert Cardona and Philipp Schlegel. Released under the GPL-3 license.
