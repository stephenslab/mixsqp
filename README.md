# mixSQP: R package for fast maximum-likelihood estimation of mixture proportions using sequential quadratic programming

[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/b7cp8eo5e7ikid0i?svg=true)](https://ci.appveyor.com/project/pcarbo/mixsqp)

The mixSQP R package provides algorithms based on sequential quadratic
programming for maximum likelihood estimation of the mixture
proportions in a finite mixture model where the component densities
are known. For large data sets (large sample sizes), the SQP algorithm
can be orders of magnitude faster than a state-of-the-art
interior-point solver implemented by the "KWDual" method in the
[REBayes package][rebayes].

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our manuscript, Kim *et al* (2018). The full
citation is given above.

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai Anitescu
> (2018). *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*
> (Submitted for review.)


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

## Credits

This project was developed by [Youngseok Kim][youngseok] at the
[University of Chicago][uchicago], with contributions from
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai].

[rebayes]: https://cran.r-project.org/package=REBayes
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
