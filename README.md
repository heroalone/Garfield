
<table>
<tr>
<td width="20%">
<img src="images/Garfield_logo_new.png" alt="Garfield Logo" width="147" height="147">
</td>
<td width="80%">
<h1> Garfield: Genetic Association by Random Forest and InterpretivE Logic Decisions</h1>
</td>
</tr>
</table>


## Overview

The conventional genetic association analysis often assumes a single causal factor underlying each locus in an additive manner. While this simplifies modeling and aids in findings, it may not capture the true complexity of genetic architecture. Garfield integrates variable selection through Random Forest and interpretation via logic gates. Its goal is to identify heterogeneity or other interactions involving multiple variants.

## Table of Contents

- [Installation](#install)
- [Inputs preparation](#input)
- [Options](#option)
- [Usage & examples](#usage)
- [License](#license)
- [Citation](#cite)
- [Contributing/Contact](#contribute)
- [Q & A](#Q_A)


## Installation <a name="install"></a>

#### Installation via `conda` (recommended)

#### Or download and install locally
```bash
git clone git@github.com:heroalone/Garfield.git
cd Garfield
perl INSTALL.pl
```

#### Dependencies

The following dependencies are required and will be installed:

1. Perl modules (perl ≥ 5 should always work, while v5.32.1 is successfully tested):
    - [Pod::Usage](https://metacpan.org/dist/Pod-Usage)
    - [Getopt::Long::Subcommand](https://metacpan.org/pod/Getopt::Long::Subcommand)
    - [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

2. R packages (R ≥ 3.1 should work, while v3.5.1 is successfully tested):
    - [genio](https://cran.r-project.org/web/packages/genio/index.html)
    - [ranger](https://cran.r-project.org/web/packages/ranger/index.html)
    - [logicFS](https://www.bioconductor.org/packages/release/bioc/html/logicFS.html)

---


## Preparation of Input Files

#### - Phenotype [--trait|-t \<file\>]

The phenotype input is a 3-column tab-delimited file with the following structure (FAM_ID, IND_ID, and values):
```
111    111    0
222    222    1
333    333    0
...
```

**Note:** The program considers negative values as "NA" by default and removes corresponding samples. Use "--trait_include_negative" to retain negative values as normal trait values.

It is recommended to remove all samples with missing phenotype values in advance and maintain the same order of samples in both phenotype and genotype files.



#### - Genotype [--genotype|-g \<file\>]

The genotype should be in [PLINK binary format](https://www.cog-genomics.org/plink/1.9/formats#bed): `file.bed`, `file.bim` and `file.fam`. You can easily convert other formats into it, for example from VCF:

```bash
plink --vcf input.vcf.gz --keep sample_overlap_with_trait.txt --indiv-sort file test_trait.txt --maf 0.02 --make-bed --out output
```

  - `--vcf` : filename of input vcf
  - `--keep` : only keep listed samples in the given file, with family IDs in the first column and within-family IDs in the second column
  - `--indiv-sort`: specify how samples should be sorted in the genotype, here the 'file' mode is used to take the order in the specify file "test_trait.txt". For other modes please see [here](https://www.cog-genomics.org/plink/1.9/data#indiv_sort).
  - `--maf` : filters out all variants with minor allele frequency below this provided threshold, default 0.01
  - `--out` : specify the prefix for the output file. Here the output will be `output.bed`, `output.bim`, and `output.fam`.
 see more other parameters in the [plink documentation](https://www.cog-genomics.org/plink/1.9/).

**Note:** Missing genotypes are not allowed (in both random forest and logic gates analyses), please use various imputation methods in advance to make 100% genotyping rates. See [Beagle](http://faculty.washington.edu/browning/beagle/beagle.html) used in our study.


### Other input files required by specific mode(s)

#### BED (Browser Extensible Data) file

```bash
--bed <file.bed>
```
Not to be confused with plink genotype bed file. The BED file here contains coordinates and gene names, with the first four columns as `[chrom start end gene_id]`. Note that the chrom IDs (1 or chr1, etc) should be consistent with that in genotype file.

#### fai index or genome file

```bash
--faidx <file.fai>
```

The `.fai` (fasta index) or `.genome` file contains coordinates for each chromosome, with the first two columns as `[chrom length]`.

#### Gene set file

```bash
--geneset <gene.set>
```

Tab-delimited text file with each row listing two or more gene IDs that will be analyzed together.

#### variant list (e.g. GWAS peaks)

```bash
--INDpeak <gwaspeak.list> 
```

A list of markers, with only one column and the name should be present in the genotype .bim file.



## Options <a name="option"></a>

Additional parameters in addition to above inputs:

```bash
--extension|-e <int>
```
- Extension of intervals for each flanking of gene. \[Default: 50000\]

```bash
--outdir|-o <string>
```
- Specify the directory of outputs. \[default: ./\]

```bash
--prefix|-p <string>
```

- Specify the prefix of output pseudo-genotypes. \[default: genotype_phenotype_subcommand\]

```bash
--temporary|-tmp <string>
```

- Specify the temporary directory for intermediate processing. \[default: ./tmp\]

```bash
--threads|-@ <int>
```

- Specify the threads that can be used. \[default: 1\]

```bash
--window|-w <int>  
```

- Window size of each sliding window. \[Default: 50000\]

```bash
--step|-s <int>
```

- Step size of each sliding window. \[Default: 25000\]

```bash
--INDpeak <gwaspeak.list>
```

- A list of markers, only one column and the name should be present in the genotype .bim file.

```bash
--rmLD2peak|-rm <float>
```

- Variants that show LD r2 above this level with that of the provided variant list are excluded \[Default: 0.3\]

```bash
--LDprune_rsq|-prune <float>
```

- Variant pruning is applied based on the given LD r2 threshold here (and the "--indep-pairwise 10 3 \<float\>" in plink with be applied); set to 1 to cancel this pruning process \[Default: 0.9\]

#### For detailed documentation locally, you can type:

```bash
Garfield --help|-h
```

or

```bash
Garfield <subcomand> --help|-h
```



## Usage & Examples <a name="usage"></a>

Garfield includes 4 subcommands:

### Gene

```bash
Garfield Gene --genotype input --trait trait --output output --temporary ./tmp --prefix gene2trait.test --threads 2 --bed gene.bed --extension 20000
```

### Window

```bash
Garfield Window --genotype input --trait trait --output output --temporary ./tmp --prefix window.trait.test --threads 5 --faidx fasta.fai --window 50000 --step 20000
```

### GeneSet

```bash
Garfield GeneSet --genotype input --trait trait --output output --temporary ./tmp --prefix geneSet.trait.test --threads 2 --bed gene.bed --geneset genepairs.txt --extension 20000
```

### Ghost

```bash
Garfield Ghost --genotype input --INDpeak gwas.peaks --extension 100000 --rmLD2peak 0.3 --LDprune_rsq 0.9 --output output --temporary ./tmp --prefix test2peak --threads 5
```

## How to cite <a name="cite"></a>

A paper describing the current study is under preparation.

## Warranty and License <a name="license"></a>

You acknowledge and agree that the software is provided to you on an **AS IS** basis, without warranty of any kind, express or implied.

Garfield is released under the GPLv3 license. See the separate license file for details. Briefly, it allows users to legally copy, distribute, and modify. However, any distribution of derivative code must also be open source and abide by the same GPLv3 agreement.

## Contributing <a name="contribute"></a>

Any feedback, comments, suggestions or contributions to help improve Garfield is highly welcome. If you have encountered any issues or have ideas for new features, please let me know. There are a few ways to contribute:

#### Contact the Developers

For questions or comments, please contact the lead developer Dr. Haijun Liu: haijun.liu@gmi.oeaw.ac.at.  
I'm also excited to hear about any new findings you may have using this software, DO let me know!

#### Submit Issues

You can report any problems, bugs or feature requests directly on our [GitHub issues page](https://github.com/heroalone/Garfield/issues).

#### Contribute Code

If you would like to contribute code, please fork the repo and submit a pull request. Any contributions, large or small, are greatly appreciated. You can also reach out if you would like to be added as a collaborator.

#### Spread the Word

Another way to support the project is to tell others about Garfield. If you find it useful in your own research, please consider citing our upcoming paper or mentioning Garfield in your publications and presentations.


## Q & A <a name="Q_A"></a>

- 
- 
- 

