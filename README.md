
<table>
<tr>
<td width="21%">
<img src="images/Garfield_logo_new.png" alt="Garfield Logo" width="147" height="147">
</td>
<td width="79%">
<h1> Garfield: Genetic Association by Random Forest and InterpretivE Logic Decisions</h1>
</td>
</tr>
</table>


## Overview
Conventional genetic association analysis often assumes a single causal factor underlying each locus in an additive manner. While this simplifies modeling and aids in findings (a lot), it may not capture the true complexity of genetic architecture. Garfield integrates variable selection through random forest and interpretation via logic gates, with the goal to identify heterogeneity or other interactions involving multiple variants.

## Table of Contents

- [Installation](#install)
- [Inputs preparation](#input)
- [Options](#option)
- [Usage & examples](#usage)
- [Outputs](#output)
- [License](#license)
- [Citation](#cite)
- [Contributing/Contact](#contribute)
- [Q & A](#Q_A)


## Installation <a name="install"></a>

#### download and install locally
```bash
git clone git@github.com:heroalone/Garfield.git
cd Garfield
perl INSTALL.pl
```


#### Dependencies
The following dependencies are required and will be installed:
```
Perl modules (perl ≥ 5 should always work, while v5.32.1 is used in the present study):
  - [Pod::Usage](https://metacpan.org/dist/Pod-Usage)
  - [Getopt::Long::Subcommand](https://metacpan.org/pod/Getopt::Long::Subcommand)
  - [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

R packages (R ≥ 3.1 should work, while v3.5.1 is used in the present study):
  - [genio](https://cran.r-project.org/web/packages/genio/index.html)
  - [ranger](https://cran.r-project.org/web/packages/ranger/index.html)
  - [logicFS](https://www.bioconductor.org/packages/release/bioc/html/logicFS.html)
```

#### Test the success of installation
```bash
Garfield Gene --genotype ./example/test.genotype \
--trait ./example/test.trait.txt \
--bed ./example/test.geneAnno.bed \
--extension 20000 \
--outdir ./test \
--temporary ./tmp \
--prefix test \
--threads 1
```

You should got the message below and two files in the `test` folder: `"Garfield.bestDNF.test.txt"` and `bash"Garfield.Geno.test.tped"`.
```bash
DONE! running time:  x wallclock secs ( x usr x sys +  x cusr  x csys =  x CPU)
  /\_/\    \~~~ MEOW ~~~
 ( o.o )
  > ^ <
```

and to see detailed help info at any time: `Garfield -h` or `Garfield Gene -h`.

---

## Preparation of Input Files <a name="input"></a>
####  Phenotype [--trait|-t \<file\>]
The phenotype input is a 3-column space-delimited file: `FAM_ID`, `IND_ID`, and `phenotype_values`
```
111 111 0
222 222 1
333 333 0
...
```

**Note:** The program by default treats negative values as missing ("NA"） and corresponding samples are removed in further analysis, which can be retaind with `--keep_negative 0`.

It is highly recommended to remove all missing phenotypes in advance and it is necessarily to maintain sample order between genotype and phenotype files.


####  Genotype [--genotype|-g \<file\>]
The genotype should be in [PLINK binary format](https://www.cog-genomics.org/plink/1.9/formats#bed): `file.bed`, `file.bim` and `file.fam`. You can easily convert other formats into it, for example from VCF:

```bash
plink --vcf input.vcf.gz \
--keep sample_overlap_with_trait.txt \
--indiv-sort file test_trait.txt \
--maf 0.02 --make-bed \
--out prefix_output
```

  - `--vcf` : filename of input vcf
  - `--keep` : only keep listed samples in the given file, with family IDs in the first column and within-family IDs in the second column
  - `--indiv-sort`: specify how samples should be sorted in the genotype, here the 'file' mode is used to take the order in the specify file "test_trait.txt". For other modes please see [here](https://www.cog-genomics.org/plink/1.9/data#indiv_sort).
  - `--maf` : filters out all variants with minor allele frequency below this provided threshold, default 0.01
  - `--out` : specify the prefix for the output file. Here the output will be `prefix_output.bed`, `prefix_output.bim`, and `prefix_output.fam`.
 see more options in [plink docs](https://www.cog-genomics.org/plink/1.9/).

**Note:** Missing genotypes are disallowed (in both random forest and logic regression analyses), please use various imputation methods in advance to make 100% genotyping rates. See [Beagle](http://faculty.washington.edu/browning/beagle/beagle.html) used in our study. It's also suggested to filter out genotypes with high initial missing rates after the imputation, since the accuracy rates of which are generally not very high.


### Other input files required by specific mode(s)

#### BED (Browser Extensible Data) file
`--bed <file.bed>` : Not to be confused with plink genotype bed file. The BED file here contains coordinates of each gene, with the first four columns as `[chrom start end gene_id]`. Chrom IDs (1 or chr1, etc) must match genotype file.

#### fai index or genome file
`--faidx <file.fai>` : The `.fai` (fasta index) or `.genome` file contains coordinates for each chromosome, with the first two columns as `[chrom length]`.

#### Gene set file
`--geneset <gene.set>` : Tab-delimited text file listing two or more gene IDs per row that will be analyzed together.

#### variant list (e.g. GWAS peaks)
`--INDpeak <gwaspeak.list>` : A list of one-column markers, which should be present in the genotype `.bim` file.


## Options <a name="option"></a>
Additional parameters in addition to above inputs:

  - `--extension|-e <int>` : Extension of intervals for each flanking of gene. \[`Default: 50000`\]
  - `--outdir|-o <string>` : Specify the directory of outputs. \[`Default: ./`\]
  - `--prefix|-p <string>` : Specify the prefix of output pseudo-genotypes. \[`Default: genotype_phenotype_subcommand`\]
  - `--temporary|-tmp <string>` : Specify the temporary directory for intermediate processing. \[`Default: ./tmp`\]
  - `--threads|-@ <int>` : Specify the threads that can be used. \[`Default: 1`\]
  - `--window|-w <int>` : Window size of each sliding window. \[`Default: 50000`\]
  - `--step|-s <int>` : Step size of each sliding window. \[`Default: 25000`\]
  - `--INDpeak <gwaspeak.list>` : A list of markers, only one column and the name should be present in the genotype ``.bim`` file.

  - `--rmLD2peak|-rm <float>` : Variants that show LD r2 above this level with that of the provided variant list are excluded. \[`Default: 0.3`\]
  - `--LDprune_rsq|-prune <float>` : Variant pruning is applied based on the given LD r2 threshold here (and the `"--indep-pairwise 10 3 \<float\>"` in plink with be applied); set to `1` to cancel this pruning process. \[`Default: 0.9`\]
  - `--keep_negative|-keep` : Taking negative trait values as missing, `o` ==> no, `1` ==> yes [`Default: 1`].

  - `--help|-h` : Show detailed documentation locally, which can be run with `Garfield --help` or `Garfield <subcomand> --help`.



## Usage & Examples <a name="usage"></a>

Garfield includes 4 subcommands:
#### Gene (based on genes, each with extending 20Kb flanking regions)
```bash
Garfield Gene --genotype input \
--trait phenotype.txt \
--outdir output \
--temporary ./tmp \
--prefix gene2trait.test \
--threads 2 \
--bed gene.bed \
--extension 20000
```

#### Window (based on sliding windows, with window size of 50Kb and step 20Kb)
```bash
Garfield Window --genotype input \
--trait phenotype.txt \
--outdir output \
--temporary ./tmp \
--prefix window.trait.test \
--threads 5 \
--faidx fasta.fai \
--window 50000 \
--step 20000
```

#### GeneSet (based on gene sets, each gene extended flanking 20Kb)
```bash
Garfield GeneSet --genotype input \
--trait phenotype.txt \
--outdir output \
--temporary ./tmp \
--prefix geneSet.trait.test \
--threads 2 \
--bed gene.bed \
--geneset genepairs.txt \
--extension 20000
```

#### Ghost (analyze potential synthetic association using genotype as phenotype for given variants)
```bash
Garfield Ghost --genotype input \
--INDpeak gwas.peaks \
--extension 100000 \
--rmLD2peak 0.3 \
--LDprune_rsq 0.9 \
--outdir output \
--temporary ./tmp \
--prefix test2peak \
--threads 5
```


## Understanding of Outputs <a name="output"></a>
Two outcomes produced from Garfield: the pseudo-genotypes `**Garfield.Geno.\*.tped**` (plink transposed PED format) and `**Garfield.bestDNF.\*.txt**` descrbing the disjunctive normal form (DNF) of each pseudo-genotype.

#### The DNF output
is a 3-column tab-delimited file: `chrom`, `marker ID`, and `DNF`. For example:
```bash
10  trait.10_1700000_1800000_10.17 rs22 \& !rs66
```
This indicates the likely presence of heterogeneity between variants of rs22 and rs66. In this expression, the `allele 1` of rs22 and the `allele 0` of rs66 would lead to `allele 1` in the pseudo-genotype, while all the other allelic combinations of rs22 and rs66 consist of `allele 0` of the pseudo-genotype.

#### The .tped genotype file 
2N+4 space-delimited, where N represents the sample size. Each row is for one variant, with the first four fields `[chrom, marker_ID, genetic position (cM, mark 0 if unknown), genomic/physical coordinate]`, and the following each two fields are genotypes for each sample listed in the .tfam file. The genomic positions are all 1 by default, you can modify it by one of the marker or the start of each gene/window. Here's an example for two variants in four samples:
```bash
10 trait.10_1700000_1800000_10.17 0 1 1 1 2 2 1 1 2 2
10 trait.10_1800000_1900000_10.18 0 1 2 2 2 2 1 1 1 1
```

The `.tfam` file must be accompanied by this `.tped`, please copy and replace the raw `.fam/.tfam` file with the same base name as .tped, this can only happen when the sample order is identical between genotype and phenotype files.
Then this genotype can be loaded into plink with `--tfile`, or used directly in association mapping.

#### association mapping
```bash
emmax -t BASE_of_tped \
-k file_of_kinship \
-c file_of_pca \
-p file_of_trait \
-o prefix_of_output \
-d 3
```

---

## How to cite <a name="cite"></a>
A paper describing the current study is under preparation.

## Warranty and License <a name="license"></a>
Please acknowledge and agree that the software is provided **AS IS**  without warranty.

Garfield is released under the `GPLv3` license. See the separate license file for details. Briefly, it allows copying and distribution with attribution and requiring derivative works also be open source.

## Contributing <a name="contribute"></a>
Feedbacks, comments, issues or ideas to improve Garfield are highly welcome. Several ways to contribute:

#### Contact
For questions or comments, please contact Dr. Haijun Liu: heroalone@qq.com.
I'm also excited to hear about any new findings you may have using this software, `**DO**` let me know!

#### Submit Issues
You can report any problems, bugs or feature requests on [GitHub issues](https://github.com/heroalone/Garfield/issues).

#### Contribute Code
Contribute code by forkingand submitting a pull request. Any contributions, large or small, are greatly appreciated. You may also reach out to be added as a collaborator.

#### Spread the Word
Support the project by telling others, and consider citing our upcoming paper or mentioning Garfield in your publications and presentations.


## Q & A <a name="Q_A"></a>

- I'm sure all the perl modules are installed sucessfully, but still got `Can't locate *.pm in @INC` error.
We tried to add the module PATH into your system, obvioursly it's failed, please try to:
 1) add `use lib '/path/to/your/module/directory';` into Garfield file under `use strict;`;
 2) add the path into @INC with `echo 'export PERL5LIB=$PERL5LIB:/path/to/Garfield/lib' >> ~/.bashrc` and `source ~/.bashrc `;
 3) or simply to run it by specifying the library path with `-I` parameter, try: ```perl -I /path/to/Garfield/lib /path/to/Garfield/Garfiled Gene -h ```

- 
