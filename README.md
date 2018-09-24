# mixSQP: Sequential Quadratic Programming for fast maximum-likelihood estimation of mixture proportions

[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/b7cp8eo5e7ikid0i?svg=true)](https://ci.appveyor.com/project/pcarbo/mixsqp)

The mixSQP package provides algorithms based on
[sequential quadratic programming][sqp] (SQP) for maximum likelihood
estimation of the mixture proportions in a finite mixture model
    where the component densities are known. For large data sets
    (large sample size), the SQP algorithm can be orders of magnitude
    faster than a state-of-the-art interior-point solver implemented
    by function 'KWDual' in the 'REBayes' package.

mix-SQP: Sequential Quadratic Programming for fast maximum-likelihood
estimation of mixture proportions

We provide code implementing optimization method for
maximum-likelihood estimation of mixture proportions, including a fast
algorithm based on sequential quadratic programming, which we call
"mix-SQP".

This repository contains code resources to accompany our research
paper,

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai Anitescu
> (2018). *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*
> (Submitted for review.)

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our manuscript, Kim *et al* (2018). The full
citation is given above.

## License

Copyright (c) 2017-2018, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html). See
file [LICENSE](LICENSE) for the full text of the license.

## Setup

*Add setup instructions here.*

This command should automatically retrieve and install the ashr package from Github. If it does not, install ashr separately using devtools:

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

[sqp]:
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
