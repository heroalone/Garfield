
<div style="display: flex; align-items: center;">
    <img src="images/Garfield_logo_new.png" alt="Garfield Logo" width="222" height="222" style="margin-right: 20px;">
    <h1 style="flex: 1;">Garfield: Genetic Association by Random Forest and InterpretivE Logic Decisions</h1>
</div>


## Project Overview <a name="overview"></a>

The current genetic association analysis mainly assumes that single causal factor underlying each locus, and usually in additive manner. This simplifies modeling and helps in a great of findings, however, it may not always reflect the true complexity of genetic architecture.
Garfield integrates variable selection by Random Forest and interpretation by logic gates, aims to identify the presence of heterogeneity or other interactions involving multiple variants.


## Table of Contents

- [Installation <a name="doc\-content\-install"></a>](#install)
- [Inputs <a name="doc\-content\-input"></a>](#input)
- [Options <a name="doc\-content\-option"></a>](#option)
- [Usage <a name="doc\-content\-usage"></a>](#usage)
- [Citation <a name="doc\-content\-cite"></a>](#cite)
- [License <a name="doc\-content\-license"></a>](#license)
- [Contributing/Contact <a name="doc\-content\-contribute"></a>](#contribute)
- [Q & A <a name="doc\-content\-Q_A"></a>](#Q_A)




## Installation <a name="install"></a>

### Installation via `conda` [recommended]


### Or download and install locally
```bash
git clone git@github.com:heroalone/Garfield.git
cd Garfield
perl INSTALL.pl
```

#### These dependencies will be required and installed
1. Perl modules (perl ≥ 5 should always work, while v5.32.1 is successfully tested):
    -  <a href="https://metacpan.org/dist/Pod-Usage" target="_blank">"Pod::Usage"</a>>
    -  ["Getopt::Long::Subcommand"](https://metacpan.org/pod/Getopt::Long::Subcommand)
    -  ["Parallel::ForkManager"](https://metacpan.org/pod/Parallel::ForkManager)

2. R packages (R ≥ 3.1 should work, while v3.5.1 is successfully tested):
    -  ["genio"](https://cran.r-project.org/web/packages/genio/index.html)
    -  ["ranger"](https://cran.r-project.org/web/packages/ranger/index.html)
    -  ["logicFS"](https://www.bioconductor.org/packages/release/bioc/html/logicFS.html)



=================
## Preparation of input files <a name="input"></a>
#### Phenotype [--trait|-t \<file\>]:
3rd column of tab delimited file with each line of FAM_ID, IND_ID, and values. For example:
```
111    111    0
222    222    1
333    333    0
...
```

Note: 
The program would consider negative values as "NA" in default, and remove corresponding samples; please specify "--trait_include_negative" to keep negative as normal trait values.

It's always suggested to remove all the samples with missing phenotype values in advance, and keep the sample of both phenotype and genotype in the same order.


#### Genotype [--genotype|-g <file>]:
in <a href="https://www.cog-genomics.org/plink/1.9/formats#bed" target="_blank">plink binary format</a>: file.bed, file.bim and file.fam. 
You could simply transform other formats into it, for example from vcf:
```bash
plink --vcf input.vcf.gz --keep [sample_overlap_with_trait.txt] --indiv-sort file test_trait.txt --maf 0.02 --make-bed --out output
```
       --vcf : filename of input vcf.
       --keep : only keep listed samples in the two columns file, with family IDs in the first column and within-family IDs in the second column.
       --indiv-sort: specify how samples should be sorted in the genotype, here the 'file' mode is used to take the order in the specify file "test_trait.txt". For other modes please see [here](https://www.cog-genomics.org/plink/1.9/data#indiv_sort).
       --maf : filters out all variants with minor allele frequency below this provided threshold, default 0.01.
       --out : specify the prefix for the output file. Here the output will be output.bed, output.bim, and output.fam.
       see more other parameters in [plink documentation](https://www.cog-genomics.org/plink/1.9/).


Note:
Missing genotypes is not allowed (in both random forest and logic gates analyses), please use various imputation methods in advance to make 100% genotypying rates. <a href="http://faculty.washington.edu/browning/beagle/beagle.html" target="_blank">Beagle</a> is used in our study.


### Other input files required by specific mode(s)

```bash
--bed <file.bed>
```
       file with coordinates and gene names, with the first four columns as [chrom start end gene_id].
       Note that the chrom IDs (1 or chr1, etc) should be consistent with that in genotype file
```bash
--faidx <file.fai>
```
       .fai (fasta index) or .genome file with coordinates for each chromosome, with the first two columns as [chrom length]
```bash
--geneset <gene.set>
```
       Tab-delimited text file with each row listing two or more gene_ids that will be analized together [gene_id1 --tab-- gene_id2]
```bash
--INDpeak <gwas.peak.list>
```
       a list of markers, only one column and the name should be present in the genotype .bim file.


## Options <a name="option"></a>
General parameters in addition to above inputs:

--extension|-e <int>
       - Extension of intervals for each flanking of gene. [Default: 50000]

--outdir|-o <string>

       - Specify the directory of outputs. [default: ./]

--prefix|-p <string>

       - Specify the prefix of output pseudo-genotypes. [default: genotype_phenotype_subcommand]

--temporary|-tmp <string>

       - Specify the temporary directory for intermediate processing. [default: ./tmp]

--threads|-@ <int>

       - Specify the threads that can be used. [default: 1]

--window|-w <int>

       - Window size of each sliding window. [Default: 50000]

--step|-s <int>

       - Step size of each sliding window. [Default: 25000]

--INDpeak <gwaspeak.list>

       - a list of markers, only one column and the name should be present in the genotype .bim file.

--rmLD2peak|-rm <float>

       - Variants that show LD r2 above this level with that of the provided variant list (--INDpeak) are excluded [Default: 0.3]

--LDprune_rsq|-prune <float>

       - Variant pruning is applied based on the given LD r2 threshold here (and the "--indep-pairwise 10 3 <float>" in plink with be applied); set to 1 to cancel this additional pruning process [Default: 0.9]

For detailed documentation locally, you can run:
```bash
Garfield --help|-h
```
or
```bash
Garfield \<subcomand\> --help|-h
```



## Usage & Examples <a name="usage"></a>
### There're 4 modes (subcommands) are currently included:
```bash
- Gene
# Run Garfield based on genes and an extended flanking interval:
  Garfield Gene --genotype input --trait trait.txt --output output --temporary ./tmp --prefix gene2trait.test --threads 2 --bed gene.bed --extension 20000
```

```bash
- Window
# Run Garfield based on sliding windows of the genome:
  Garfield Window --genotype input --trait trait --output output --temporary ./tmp --prefix window.trait.test --threads 2 --faidx fasta.fai --window 50000 --step 20000
```

```bash
- GeneSet
# Run Garfield based on given gene sets:
  Garfield GeneSet --genotype input --trait trait --output output --temporary ./tmp --prefix geneSet.trait.test --threads 2 --bed gene.bed --geneset genepairs.txt --extension 20000
```

```bash
- Ghost
# Run Garfield to search potential synthetic association to individual peak genotypes, within an specific flanking interval:
  Garfield Ghost --genotype input --INDpeak peak.list --extension 100000 --rmLD2peak 0.3 --LDprune_rsq 0.9 --output output --temporary ./tmp --prefix test2peak --threads 2

```


## How to cite <a name="cite"></a>
A paper describing the current study is under preparation.


## Warranty and License  <a name="license"></a>
You acknowledge and agree that the software is provided to you on an "AS IS" basis, without warranty of any kind, express or implied.

Garfield is released under the GPLv3, which is provided in a separate license file. Briefly, it allows users to legally copy, distribute, and modify. However, any distribution of derivative code must also be open source and abide by the same GPLv3 agreement..

## Contributing <a name="contribute"></a>
If Garfield really helps you with new findings, or you have further comments or suggestions, please contact: 
<a name="contact">Dr. Haijun Liu (haijun.liu[at]gmi.oeaw.ac.at).</a>

If you have any good ideas to improve (either programmatically or algorithmically) and want to contribute, please feel free to make a pull request or let me know to invite you as collaborators; EVERY LITTLE BIT HELPS AND IS ALWAYS WELCOME!

You could report any problems and bugs directly here as well: https://github.com/heroalone/Garfield/issues.


## Q & A <a name="Q_A"></a>

