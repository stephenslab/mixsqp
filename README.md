# mixSQP: R package for fast maximum-likelihood estimation of mixture proportions using sequential quadratic programming

[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/b7cp8eo5e7ikid0i?svg=true)](https://ci.appveyor.com/project/pcarbo/mixsqp)

The mixSQP R package provides algorithms based on
[sequential quadratic programming][sqp] for maximum likelihood
estimation of the mixture proportions in a finite mixture model where
the component densities are known. For large data sets (large sample
sizes), the SQP algorithm can be orders of magnitude faster than the
state-of-the-art interior-point solver implemented in the
[REBayes package][rebayes]. For more details on the mixSQP algorithm,
please see [our paper on arXiv][arxiv-paper].

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
[MIT license](https://opensource.org/licenses/mit-license.html). See
file [LICENSE](LICENSE) for the full text of the license.

## Setup

*Add setup instructions here.*

NOTE: Proper C++ compiler setup is required.

This command should automatically retrieve and install the ashr
package from Github. If it does not, install ashr separately using
devtools:

library(devtools)
install_github("stephens999/ashr")

Alternatively, if you have cloned the repository locally, you can
install the package with the `install_local` function:

```R
```

## Developer notes

+ To install and test the mixSQP R package, run the following commands
in the shell:

```bash
R CMD build mixSQP
R CMD INSTALL mixSQP_0.1-11.tar.gz
R CMD check --as-cran mixSQP_0.1-11.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

## Credits

This project was developed by [Youngseok Kim][youngseok] at the
[University of Chicago][uchicago], with contributions from
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai].

[sqp]: https://neos-guide.org/content/sequential-quadratic-programming
[arxiv-paper]: https://arxiv.org/abs/1806.01412
[issues]: https://github.com/youngseok-kim/mixsqp/issues
[rebayes]: https://cran.r-project.org/package=REBayes
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
