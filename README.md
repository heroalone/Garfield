
<div style="display: flex; align-items: center;">
    <img src="images/Garfield_logo_new.png" alt="Garfield Logo" width="222" height="222" style="margin-right: 20px;">
    <h1 style="flex: 1;">Garfield: Genetic Association by Random Forest and InterpretivE Logic Decisions</h1>
</div>

### Full documentation at: https://Garfield.readthedocs.io/ 

## Table of Contents
=================

- [Overview <a name="user\-content\-workflow"></a>](#overview)
- [Installation <a name="user\-content\-install"></a>](#install)
- [Usage](#usage)
- [Contributing](#contribute)
- [License](#license)
- [Contact](#contact)

## Project Overview <a name="overview"></a>

Garfield - Genetic Association by Random Forest and InterpretivE Logic Decisions, helps to identify genetic heterogeneity and various interaction effects.


## Installation <a name="install"></a>

Provide instructions on how to install and set up the project locally. Include any dependencies or prerequisites that need to be installed.
### Installation via `conda` [recommended]


### Installation of dependencies step by step

Download the lastest version:

```
```

1. Perl packages:
  
    ```
2. R packages
    ```

## Preparation of input files <a name="usage"></a>

## Usage by examples
```bash
```
# Case1: gene-based
$ command 1


Here is an example:

```
Garfield 
```

Options
    -help | -h
            Prints the help message and exits.

    --input [required]
           - RDS files. <fig1.rds,fig2.rds...>

    --labels [optional]
           - Labesl for each figure. Default: <A,B,C,D...>

    --output [optional]
           - Output files for the graph. Default: cowplot_mer_fig.pdf

    --ncol [optional]
           - Number of columns on the graph.

    --base_height [optional]
           - The height (in inches) of each sub-plot

    --base_aspect_ratio [optional]
           -  The aspect ratio of each sub-plot. Default: 1.6
```
```

# Case2: paralog-based
$ command 2

# Case3: potential synthetic associations
$ command 3

Please see details in [Document](http://xxx) websites.


Trait:
missing: NA,
negative: default is considered as NA, please specify --trait include_negative to keep these as normal trait values;

Best suggestion: removing all NAs, and keep the id of trait and plink genotypes in the same order.


plink --indiv-sort <mode name> [filename]
# This allows you to specify how samples should be sorted when generating new datasets. The four modes are:


genotype：no-missing 
/groups/nordborg/user/haijun.liu/software/plink2 --bfile cubicV3.maf0.02.chr1 --genotyping-rate should be 1.




## How to cite <a name="cite"></a>
A paper describing the current study is under preparation.


## Warranty and License  <a name="license"></a>
You acknowledge and agree that the software is provided to you on an "AS IS" basis, without warranty of any kind, express or implied.

Garfield is released under the GPLv3, which is provided in a separate license file. Briefly, it allows users to legally copy, distribute, and modify. However, any distribution of derivative code must also be open source and abide by the same GPLv3 agreement..

## Contributing <a name="contribute"></a>
If Garfield really helps you with new findings, or you have further comments or suggestions, please contact:<a name="contact"></a><br>
Dr. Haijun Liu (haijun.liu[at]gmi.oeaw.ac.at).

If you have any good ideas to improve (either programmatically or algorithmically) and want to contribute, please feel free to make a pull request or let me know to invite you as collaborators; YOUR HELP IS ALWAYS VERY WELCOME!

You could report any problems and bugs directly here as well: https://github.com/heroalone/Garfield/issues.
