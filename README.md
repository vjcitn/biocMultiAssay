# MultiAssayExperiment <a href='https://waldronlab.io/MultiAssayExperiment'><img src='https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/MultiAssayExperiment/MultiAssayExperiment.png' align="right" height="100" /></a>

## Software For The Integration Of Multi-Omics Experiments In Bioconductor

https://doi.org/10.1158/0008-5472.CAN-17-0344

[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/MultiAssayExperiment.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MultiAssayExperiment)
[![Platforms](http://www.bioconductor.org/shields/availability/release/MultiAssayExperiment.svg)](https://www.bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html#archives)
[![wercker status](https://app.wercker.com/status/2aa523f23142715771256b85187d7bcb/s/devel "wercker status")](https://app.wercker.com/project/byKey/2aa523f23142715771256b85187d7bcb)
[![Coverage Status](https://codecov.io/github/waldronlab/MultiAssayExperiment/coverage.svg?branch=devel)](https://codecov.io/github/waldronlab/MultiAssayExperiment?branch=devel)
[![Downloads](http://www.bioconductor.org/shields/downloads/release/MultiAssayExperiment.svg)](https://bioconductor.org/packages/stats/bioc/MultiAssayExperiment)

## Installation

We recommend installing the stable release version of MultiAssayExperiment in
Bioconductor. This can be done using `BiocManager`:

```
if (!require("BiocManager"))
    install.packages("BiocManager")

library(BiocManager)

install("MultiAssayExperiment")
```

## Schematic

Here is a visual overview of the `MultiAssayExperiment` class.
The three main components are:

* `colData` - phenotype data
* `ExperimentList` - a list of tables from experimental results
* `sampleMap` - a graph representation of 'samples' and 'participants'

_Note_. For simplicity, we use the terms 'samples' and 'participants' to
elucidate the relationship although `MultiAssayExperiment` can work on other
biological specimens.

<img src="vignettes/MultiAssayExperiment.png" alt="MultiAssayExperiment schematic"/>

## Cheatsheet

<a href="https://github.com/waldronlab/cheatsheets/blob/main/MultiAssayExperiment_QuickRef.pdf"><img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/main/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/></a>

## Ready-to-use `MultiAssayExperiment` objects

For easy-to-use and ready-made MultiAssayExperiment objects, use the
`curatedTCGAData` experiment data package.

```
install("curatedTCGAData")
```

## Companion package for working with TCGA data

TCGAutils is a handy package for working with `MultiAssayExperiment` data
objects from `curatedTCGAData`. It is highly recommended to use `TCGAutils` for
identifier manipulation, sample identification and more.

```
install("TCGAutils")
```

## Documentation

The `MultiAssayExperiment` API is available by browsing to the
[API wiki](https://github.com/waldronlab/MultiAssayExperiment/wiki/MultiAssayExperiment-API).

## The `MultiAssayExperiment` Bioconductor Special Interest Group

This group meets remotely to discuss this project approximately every 3 weeks.
If you are interested, please join the
[MultiAssayExperiment Google Group](https://groups.google.com/forum/#!forum/biocmultiassay)
and see the
[calendar](https://www.google.com/calendar/embed?src=9ar0qc8mpkv6b9intgmdcdf0ss%40group.calendar.google.com&ctz=America/New_York)
of upcoming meetings.

## Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect
all people who contribute through reporting issues, posting feature
requests, updating documentation, submitting pull requests or patches,
and other activities.

We are committed to making participation in this project a
harassment-free experience for everyone, regardless of level of
experience, gender, gender identity and expression, sexual
orientation, disability, personal appearance, body size, race, age, or
religion.

Examples of unacceptable behavior by participants include the use of
sexual language or imagery, derogatory comments or personal attacks,
trolling, public or private harassment, insults, or other
disrespectful conduct.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct. Project
maintainers who do not follow the Code of Conduct may be removed from
the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior
may be reported by opening an issue or contacting one or more of the
project maintainers.

This Code of Conduct is adapted from the [Contributor
Covenant](http://contributor-covenant.org), version 1.0.0, available
at
[http://contributor-covenant.org/version/1/0/0/](http://contributor-covenant.org/version/1/0/0/)
