[![Release](https://img.shields.io/github/release/birollab/ntRoot.svg)](https://github.com/birollab/ntRoot/releases)
[![link](https://img.shields.io/badge/ntRoot-manuscript-brightgreen)](https://doi.org/10.1093/bioadv/vbaf287)
[![link](https://img.shields.io/badge/HostSeq-manuscript-yellow)](https://doi.org/10.1139/gen-2026-0022)
[![Zenodo](https://img.shields.io/badge/ntRoot-zenodo-red)](https://zenodo.org/doi/10.5281/zenodo.10869033)
[![Conda](https://img.shields.io/conda/dn/bioconda/ntroot?label=Conda)](https://anaconda.org/bioconda/ntroot)

![Logo](https://github.com/birollab/ntRoot/blob/main/ntroot-logo_colors.png)

# ntRoot

Scalable Ancestry Predictions from Genomic Data

## Contents

1. [Credit](#credit)
2. [Description](#description)
3. [Installation](#install)
4. [Dependencies](#dependencies)
5. [Usage](#usage)
6. [Interactive LAI HTML visualization](#lai)
7. [Human ancestry predictions](#data)
8. [Demo](#demo)
9. [Documentation](#docs)	
10. [Citing](#citing)
11. [License](#license)

## Credit  <a name=credit></a>
Written by René L Warren and Lauren Coombe

## Description <a name=description></a>
ntRoot is a framework for ancestry inference from genomic data, offering both Local Ancestry Inference (LAI) and Global Ancestry Inference (GAI). Leveraging integrated variant call sets from the 1000 Genomes Project (1kGP), ntRoot provides accurate predictions(1) of human super-population ancestry with speed and efficiency from Whole Genome Sequencing (WGS) datasets and complete or draft-stage Whole Genome Assemblies (WGA). Through streamlined processing, flexible genomic input, and integrated local ancestry visualization, ntRoot holds promises for human ancestry inference of small-to-large patient/individual cohorts, enabling association studies with demographics and facilitating deeper insights into population genetics and disease risk factors.

(1) Tested on base-accurate quality data, including Illumina short read and newer nanopore (ONT, KitV14) & PacBio CCS HiFi long read datasets, complete reference genomes and polished, Oxford Nanopore Technology long read GoldRush, Flye and Shasta draft genome assemblies 

## Installation <a name=install></a>

Installing ntRoot using conda (recommended):
```
conda install -c bioconda -c conda-forge ntroot
```

Installing ntRoot from the source code:
```
git clone https://github.com/birollab/ntRoot.git
cd ntRoot
```
No compilation is required for ntRoot (only the dependencies), so simply add the ntRoot repository to your PATH.

### Dependencies <a name=dependencies></a>
- python 3.9+
- perl
- [ntEdit 2.0.2+](https://github.com/birollab/ntEdit)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [samtools](https://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [plotly](https://plotly.com/)
- [pandas](https://pandas.pydata.org/)


## Usage <a name=usage></a>
```
usage: ntroot [-h] [-r REFERENCE] [--reads READS] [--genome GENOME [GENOME ...]] -l L [--exome] [-k K] [--tile TILE] [--lai] [-t T] [-z Z] [-j J] [--cutoff CUTOFF] [-Y Y]
              [--custom_vcf CUSTOM_VCF] [--masked] [--exome_bed EXOME_BED] [--strip_info] [-v] [-V] [-n] [-f]

ntRoot: Ancestry inference from genomic data

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Reference genome (FASTA, Multi-FASTA, and/or gzipped compatible)
  --reads READS         Prefix of input reads file(s) for detecting SNVs. All files in the working directory with the specified prefix will be used. (fastq, fasta, gz, bz, zip)
  --genome GENOME [GENOME ...]
                        Genome assembly file(s) for detecting SNVs compared to --reference
  -l L                  input IVC VCF file with annotated variants (e.g., 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz, clinvar.vcf, etc.)
  --exome               Input reads for detecting SNVs are from whole exome sequencing. If provided, must also specify either --exome_bed or --masked. --cutoff 2 is implied unless otherwise specified.
  -k K                  k-mer size
  --tile TILE           Tile size for ancestry fraction inference (bp) [default=5000000]
  --lai                 Output ancestry predictons per tile in a separate output file
  --recomb_bed RECOMB_BED
                        Recombination-tile BED file, tab-separated with columns: chrom, start, end, name, and an optional cM width. Chromosome IDs must match the reference/VCF. When provided, ntRoot runs in recombination mode and outputs recombination-resolution LAI and a recombination-based global ancestry (instead of the fixed-tile outputs).
  -t T                  Number of threads [default=4]
  -z Z                  Minimum contig length [default=100]
  -j J                  controls size of k-mer subset. When checking subset of k-mers, check every jth k-mer [default=3]
  --cutoff CUTOFF       Minimum coverage of k-mers in ntEdit Bloom filter. Solid k-mers are used if set to 0 [0]
  -Y Y                  Ratio of number of k-mers in the k subset that should be present to accept an edit (higher=stringent) [default=0.55]
  --custom_vcf CUSTOM_VCF
                        Input VCF for computing ancestry. When specified, ntRoot will skip the ntEdit step, and predict ancestry from the provided VCF.
  --masked              Exome Mode (--exome) only: Indicates that the reference genome provided with --reference has all NON-targeted exonic regions hard-masked. 
  --exome_bed EXOME_BED
                        Exome Mode (--exome) only: BED file of exome targeted regions. 
  --strip_info          When using --custom_vcf, strip the existing INFO field from the input VCF.
  -v, --verbose         Verbose mode [default=False]
  -V, --version         show program's version number and exit
  -n, --dry-run         Print out the commands that will be executed
  -f, --force           Run all ntRoot steps, regardless of existing output files

Note: please specify --reads OR --genome (not both)
If you have any questions about ntRoot, please open an issue at https://github.com/birollab/ntRoot
```

## Interactive Local Ancestry (LAI) visualization <a name=lai></a>

When `--lai` is specified, ntRoot automatically generates an interactive HTML visualization of chromosomal local ancestry assignments as part of the Snakemake workflow.

The visualization enables interactive exploration of chromosome-specific and genome-wide ancestry patterns through linked chromosome, ancestry, and summary views.

The resulting visualization provides:

- Chromosome-length-scaled local ancestry tracks
- Global ancestry composition summaries
- Interactive exploration of chromosome-specific ancestry patterns
- Interactive exploration of genome-wide ancestry patterns
- Hover and click highlighting across linked views
- Standalone HTML output viewable in any modern web browser

Users can interactively explore ancestry patterns within individual chromosomes or across the entire genome using the generated standalone HTML report. 

Example visualizations:
<br>
[Interactive LAI demo1 (HTML)](https://www.birollab.ca/viz/ERR3988858_ntedit_k55_variants.vcf_ntroot-lai-interactive_tile5000000.html)
<br>
[Interactive LAI demo2 (HTML)](https://www.birollab.ca/viz/ERR3242010_ntedit_k55_variants.vcf_ntroot-lai-interactive_tile5000000.html)



## Human ancestry predictions <a name=data></a>

Using the 1kGP integrated variant call set.

Download this archive:
<pre>
wget https://zenodo.org/records/22151063/files/ntroot_supplementary_zenodo.tar.gz
</pre>
  
from:
<pre>
https://doi.org/10.5281/zenodo.10869033
</pre>

unzip and untar:
<pre>
tar xvzf ntroot_supplementary_zenodo.tar.gz
</pre>

access the files:
<pre>
cd ./ntroot_supplementary_zenodo/data
ls

1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz
genetic_map_to_bed.py
GRCh38.fa.gz
readme
recomb_tiles_0.5cM.bed
</pre>

**Running ntRoot with whole genome sequencing reads or genome assemblies**

Users will specify:
<pre>
ntroot --reference GRCh38.fa.gz (--reads FILE_PREFIX OR --genome FILE) -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz -k 55
</pre>

Example command:
<pre>
ntroot -k 55 --reference GRCh38.fa.gz --reads ERR3242308_ -t 48 -Y 0.55 -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz
</pre>

To use a BED file to specify tile coordinates (instead of using the default fixed-size tiles), run:
<pre>
ntroot -k 55 --reference GRCh38.fa.gz --reads ERR3242308_ -t 48 -Y 0.55 -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz --recomb_bed recomb_tiles_0.5cM.bed
</pre>
The `recomb_tiles_0.5cM.bed` file is provided in the ntRoot supplementary archive as an example.

**Running ntRoot with whole exome sequencing**

If your input reads are from whole exome sequencing, the regions of your reference genome that are NOT targeted exonic regions should be hard-masked (converted to Ns):
<pre>
ntroot -k 55 --exome --masked --reference masked_GRCh38.fa.gz --reads test_exome -t 48 -Y 0.55 -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz
</pre>
ntRoot can perform the masking automatically if you do not already have a masked reference file. In that case, provide a BED file with all the targeted regions, and ntRoot will use bedtools to mask the reference regions that are NOT targeted regions:
<pre>
ntroot -k 55 --exome --exome_bed exome_seq_regions.bed --reference GRCh38.fa.gz --reads test_exome -t 48 -Y 0.55 -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz
</pre>

**Running ntRoot using a pre-existing VCF file:**

If you would like to infer ancestry from a pre-existing VCF file:
<pre>
ntroot -r GRCh38.fa.gz --custom_vcf third_party.vcf -l 1000GP_integrated_snv_v2a_27022019.GRCh38.phased_gt1.vcf.gz
</pre>

Note: For more advanced users, and for ancestry predictions on organisms other than human, please contact us.


## Demo <a name=demo></a>
To test your installation:
```
cd demo
./run_ntroot_demo.sh
```
Ensure that the ntRoot installation is available on your PATH.

For example interactive Local Ancestry Inference (LAI) visualizations generated by ntRoot, see:
<br>
[Interactive LAI demo1 (HTML)](https://www.birollab.ca/viz/ERR3988858_ntedit_k55_variants.vcf_ntroot-lai-interactive_tile5000000.html)
<br>
[Interactive LAI demo2 (HTML)](https://www.birollab.ca/viz/ERR3242010_ntedit_k55_variants.vcf_ntroot-lai-interactive_tile5000000.html)


## Documentation <a name=docs></a>

Refer to the README.md file on how to install and run ntRoot.
Our [preprint](https://doi.org/10.1101/2024.03.26.586646) contains information about the software and its performance.
![ntRoot PSB poster](https://github.com/birollab/ntRoot/blob/main/ntRootPSB2025.png)
This [PSB2025 poster](https://doi.org/10.7490/f1000research.1119849.1), also available on [zenodo](https://doi.org/10.5281/zenodo.13844277), contains additional information, benchmarks and results.


## Citing <a name=citing></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/birollab/ntRoot.svg)](https://github.com/birollab/ntRoot/stargazers) and for using and promoting this free software! We hope that ntRoot is useful to you and your research.

If you use ntRoot, please cite:

[ntRoot: human ancestry inference at scale, from genomic data](https://doi.org/10.1093/bioadv/vbaf287)
<pre>
Warren RL, Coombe L, Wong J, Kazemi P, Birol I.
ntRoot: computational inference of human ancestry at scale from genomic data. 
Bioinformatics Advances, vbaf287, https://doi.org/10.1093/bioadv/vbaf287
</pre>

Citing the poster:
<pre>
Warren R, Coombe L, Wong J et al. ntRoot: Scalable ancestry predictions from genome sequencing data [version 1]. F1000Research 2024, 13:1021 (poster) (https://doi.org/10.7490/f1000research.1119849.1)
</pre>

See how we applied `ntRoot` to [infer ancestry for 10,250 Canadian participants from the HostSeq cohort](https://doi.org/10.1101/2025.06.10.658783):
<pre>
Concordance and dissonance: A genome-wide analysis of self-declared versus inferred ancestry in 10,250 participants from the HostSeq cohort
René L Warren, Inanc Birol, CGEn HostSeq Initiative
bioRxiv 2025.06.10.658783; doi: https://doi.org/10.1101/2025.06.10.658783
</pre>

## License <a name=license></a>

ntRoot Copyright (c) 2024-present British Columbia Cancer Agency Branch.  All rights reserved.

ntRoot is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
