# mixSQP: R package for fast maximum-likelihood estimation of mixture proportions using sequential quadratic programming

[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/b7cp8eo5e7ikid0i?svg=true)](https://ci.appveyor.com/project/pcarbo/mixsqp)

The mixSQP R package provides algorithms based on
[sequential quadratic programming][sqp] for maximum likelihood
estimation of the mixture proportions in a finite mixture model where
the component densities are known. The SQP algorithm is expected to
obtain solutions that are at least as accurate as the state-of-the-art
MOSEK interior-point solver (called via the "KWDual" function in
the [REBayes package][rebayes]), and is expected to compute these
solutions much more quickly in large data sets.

For more details on the mixSQP algorithm, please see
[our paper on arXiv][arxiv-paper].

See also the [Julia implementation of mixSQP][mixsqp-julia], which is
currently faster than the R implementation, particularly for large
data sets.

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find that the mixSQP package, or any of the source code in this
repository, is useful for your work, please cite our paper:

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai
> Anitescu. *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*
> [arXiv:1806.01412][arxiv-paper].

## License

Copyright (c) 2017-2018, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the
[MIT license][mit-license]. See
file [LICENSE](LICENSE) for the full text of the license.

## Setup

This command should automatically retrieve and install the ashr
package from Github. If it does not, install ashr separately using
devtools:

To install the latest version of the mixSQP package from GitHub,
use [devtools][devtools]:

```R
install.packages("devtools")
library(devtools)
install_github("youngseok-kim/mixsqp",build_vignettes = TRUE,
               upgrade_dependencies = FALSE)
```

This command should automatically install all required packages if
they are not installed already.

Alternatively, if you have cloned the repository locally, you can
install the package with the `install_local` function from
devtools. Assuming your working directory contains the mixSQP
repository, run this code to install the package:

```R
library(devtools)
list.files(pattern = "mixSQP") # Should output "mixSQP".
nstall_local("mixSQP",build_vignettes = TRUE,upgrade_dependencies = FALSE)
```

**Important note:** Compiling the mixSQP package from source will
require a C++ compiler setup that is appropriate for the the R
installed on your computer. For details, refer to the
[CRAN documentation][cran]. For Mac computers, see
[these notes][compiling-macos].

## Developer notes

### Testing the package

To install and test the mixSQP R package, run the following commands
in the shell:

```bash
R CMD build mixSQP
R CMD INSTALL mixSQP_0.1-41.tar.gz
R CMD check --as-cran mixSQP_0.1-41.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

### Updating the C++ source and documentation

When any changes are made to [roxygen2][roxygen2] markup or to the C++
code in the src directory, simply run `devtools::document()` to
update the [RcppExports.cpp](src/RcppExports.cpp), the NAMESPACE file,
and the package documentation files in the man directory.

### Updating the pkgdown site

Run this line of R code to build the website (make sure you have an
Internet connection while running the code):

```R
pkgdown::build_site(mathjax = FALSE)
```

### Other developer notes

+ Add file `pre-commit` to `.git/hooks` in the git repository to
prevent commits that don't include a change to the package version.

## Credits

This project was developed by [Youngseok Kim][youngseok] at the
[University of Chicago][uchicago], with contributions from
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai].

[sqp]: https://neos-guide.org/content/sequential-quadratic-programming
[arxiv-paper]: https://arxiv.org/abs/1806.01412
[mixsqp-julia]: https://github.com/stephenslab/mixsqp-paper
[issues]: https://github.com/youngseok-kim/mixsqp/issues
[rebayes]: https://cran.r-project.org/package=REBayes
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
