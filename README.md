# mixsqp: Fast maximum-likelihood estimation of mixture proportions using sequential quadratic programming

[![CRAN status badge](https://www.r-pkg.org/badges/version/mixsqp)](https://cran.r-project.org/package=mixsqp)
[![Travis Build Status](https://travis-ci.org/stephenslab/mixsqp.svg?branch=master)](https://travis-ci.org/stephenslab/mixsqp)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/i8744qet66w5uhe2?svg=true)](https://ci.appveyor.com/project/pcarbo/mixsqp)
[![codecov](https://codecov.io/gh/stephenslab/mixsqp/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/mixsqp)

The mixsqp R package provides algorithms based on [sequential 
quadratic programming][sqp] for maximum likelihood estimation of the
mixture proportions in a finite mixture model where the component
densities are known. The SQP algorithm is expected to obtain solutions
that are at least as accurate as the state-of-the-art MOSEK
interior-point solver (called via the "KWDual" function in the
[REBayes package][rebayes]), and is expected to compute these
solutions much more quickly in large data sets.

For more details on the methods, please see the [journal
paper][jcgs-paper] or the [arXiv preprint][arxiv-paper].

The methods were originally implemented in [Julia][julia];
please see [here][mixsqp-julia] for the Julia implementation.

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find the mixsqp package or any of the source code in this
repository useful for your work, please cite:

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai
> Anitescu. [A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.][jcgs-paper]
> To appear in the *Journal of Computational and Graphical Statistics.*

## License

Copyright (c) 2017-2020, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Quick Start

Install mixsqp from [CRAN](http://www.r-pkg.org/pkg/varbvs):

```R
install.packages("mixsqp")
```

For more detailed installation instructions, see the "Setup" section
below.

Once you have installed the package, load the package in R:

```R
library(mixsqp)
```

Next, run the small example provided with the mixsqp function:

```R
example("mixsqp")
```

For a more detailed illustration of the SQP algorithm applied to the
problem of computing maximum-likelihood estimates for a mixture model,
read through the [introductory vignette][mixsqp-vignette].

To learn more, visit the [package website][mixsqp-website], or view
the "mixsqp" help page:

```R
help("mixsqp")
```

## Setup

To install mixsqp from [CRAN](http://www.r-pkg.org/pkg/varbvs), in R
run:

```R
install.packages("mixsqp")
```

Alternatively, to install the latest version of the mixsqp package
from GitHub, use [devtools][devtools]:

```R
install.packages("devtools")
library(devtools)
install_github("stephenslab/mixsqp",build_vignettes = TRUE)
```

This command should automatically install all required packages if
they are not installed already.

If you have cloned the repository locally, you can install the package
with the `install_local` function from devtools. Assuming your working
directory contains the mixsqp repository, run this code to install the
package:

```R
library(devtools)
install_local("mixsqp")
```

### Additional setup notes

+ Compiling the mixsqp package from source will require a C++ compiler
setup that is appropriate for the the R installed on your
computer. For details, refer to the [CRAN documentation][cran]. For
Mac computers, see [these notes][compiling-macos].

+ To use the (optional) alternative solver, "mixkwdual", which is
mostly useful for comparisons of the different optimization methods,
you will need to install the [REBayes][rebayes] package. The REBayes
package, in turn, requires the [Rmosek][mosek] package. Refer to the
[MOSEK documentation][mosek-docs] for instructions on installing the
Rmosek package. Once you have followed these steps, you can run
[this example](inst/code/test.rmosek.R) to verify that Rmosek is
correctly installed. Installation of the REBayes package also allows
you to build [the vignette][mixsqp-vignette] and view it locally:

```R
devtools::install_github("stephenslab/mixsqp",build_vignettes = TRUE)
library(mixsqp)
vignette("mixsqp-intro")
```

## Developer notes

### Testing the package

To install and test the mixsqp R package, run the following commands
in the shell:

```bash
R CMD build mixsqp
R CMD INSTALL mixsqp_0.3-31.tar.gz
R CMD check --as-cran mixsqp_0.3-31.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

### Updating the C++ source and documentation

When any changes are made to [roxygen2][roxygen2] markup or to the C++
code in the `src` directory, simply run `devtools::document()` to 
update the [RcppExports.cpp](src/RcppExports.cpp), the NAMESPACE file,
and the package documentation files in the `man` directory.

### Updating the pkgdown site

Run this line of R code to build the website (make sure you have an
Internet connection while running the code):

```R
pkgdown::build_site()
```

Version 1.4.1 of pkgdown was used.

### Other developer notes

Add file `pre-commit` to `.git/hooks` in the git repository to
prevent commits that don't include a change to the package version.

## Credits

The mixsqp R package was developed by [Youngseok Kim][youngseok] and
[Peter Carbonetto][peter] at the [University of Chicago][uchicago],
with contributions from [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai].

[mixsqp-website]: https://stephenslab.github.io/mixsqp
[mixsqp-vignette]: https://stephenslab.github.io/mixsqp/articles/mixsqp-intro.html
[sqp]: https://neos-guide.org/content/sequential-quadratic-programming
[arxiv-paper]: https://arxiv.org/abs/1806.01412
[jcgs-paper]: https://doi.org/10.1080/10618600.2019.1689985
[julia]: https://julialang.org
[mixsqp-julia]: https://github.com/stephenslab/mixsqp-paper
[issues]: https://github.com/stephenslab/mixsqp/issues
[rebayes]: https://cran.r-project.org/package=REBayes
[mosek]: https://www.mosek.com
[mosek-docs]: https://www.mosek.com/documentation
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
[mit-license]: https://opensource.org/licenses/mit-license.html
[devtools]: https://github.com/r-lib/devtools
[roxygen2]: https://cran.r-project.org/package=roxygen2
[cran]: https://cran.r-project.org
[compiling-macos]: https://pcarbo.github.io/pcarbo/r-macos.html
