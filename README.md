# SkewX

  [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
  [![DOI]()](https://doi.org/10.1101/gr.279396.124)

## Introduction

**SkewX** is a Nextflow pipeline to measure skewed X inactivation from long-read sequencing of native DNA, with either PacBio or Nanopore technologies.

The pipeline accepts haplotagged BAM files (with `HP` and `PS` tags and 5mCG methylation information) as the default input. It clusters reads based on their methylation profile over CpG islands and uses the haplotype information already embedded in the BAM to measure skew in X inactivation.

Alternatively, the full pipeline can be run from unphased modbam files, performing variant calling with DeepVariant and phasing with WhatsHap before the methylation analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.

## Pipeline summary

### Default mode: haplotagged BAM input

The required input is a haplotagged BAM file with `HP`/`PS` tags and 5mCG methylation information. Then:

1. Index the BAM file (skipped if a `.bai` index already exists alongside the BAM)
2. Calculate coverage with [Mosdepth](https://github.com/brentp/mosdepth)
3. Extract haplotype-tagged reads over CpG islands with [Samtools](https://www.htslib.org/)
4. Cluster reads based on methylation profile with [NanoMethViz](https://www.bioconductor.org/packages/release/bioc/html/NanoMethViz.html)
5. Derive phase block statistics from `HP`/`PS` BAM tags (chrX only)
6. Measure skew in X inactivation and generate a report for each individual

### Full pipeline mode (`--phased_bam false`)

The required input is modbam files with 5mCG information. Then:

1. If the reads are not already aligned, align to the reference genome with [Minimap2](https://github.com/lh3/minimap2)
2. If multiple samples per individual are present, merge them into a single BAM file
3. Call variants with [DeepVariant](https://github.com/google/deepvariant)
4. Phase SNPs with [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html)
5. Haplotype-tag reads with [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html)
6. Calculate coverage with [Mosdepth](https://github.com/brentp/mosdepth)
7. Extract haplotype-tagged reads over CpG islands with [Samtools](https://www.htslib.org/)
8. Cluster reads based on methylation profile with [NanoMethViz](https://www.bioconductor.org/packages/release/bioc/html/NanoMethViz.html)
9. Measure skew in X inactivation and generate a report for each individual

## Quick Start

1. Install or module load [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility.

3. If using Singularity, ensure it is mounted to your home directory:
   ```bash
   export NXF_SINGULARITY_HOME_MOUNT=true
   ```

4. Prepare a samplesheet (see [Input](#input) below).

5. Start running your analysis:

   ```bash
   nextflow main.nf --input samplesheet.csv --outdir skew_results --cgi_bedfile additional_files/CGIs_CHM13v2.0_chrX.bed -profile singularity
   ```

## Input

### Default mode: haplotagged BAM samplesheet

The samplesheet must be a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `individual` | Unique identifier for the individual |
| `sample` | Sample/tissue name (e.g. `blood`, `NSC`) |
| `haplotagged_bam` | Path to the haplotagged BAM file (must contain `HP` and `PS` tags and 5mCG methylation) |

Example:
```csv
individual,sample,haplotagged_bam
sample1,blood,/path/to/sample1_blood.hp.bam
sample1,NSC,/path/to/sample1_NSC.hp.bam
sample2,blood,/path/to/sample2_blood.hp.bam
```

Multiple samples (tissues) per individual are supported — they will be grouped for reporting.

### Full pipeline mode samplesheet (`--phased_bam false`)

| Column | Description |
|--------|-------------|
| `individual` | Unique identifier for the individual |
| `sample` | Sample/tissue name |
| `modbam_5mCG` | Path to the modbam file with 5mCG methylation calls |

Example:
```csv
individual,sample,modbam_5mCG
sample1,blood,/path/to/sample1_blood.bam
sample1,NSC,/path/to/sample1_NSC.bam
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to samplesheet CSV |
| `--outdir` | required | Output directory |
| `--cgi_bedfile` | required | BED file of CpG islands on chrX |
| `--phased_bam` | `true` | Start from haplotagged BAM (default). Set to `false` to run the full pipeline from unphased modbam |
| `--reference` | `null` | Reference FASTA (required only when `--phased_bam false`) |
| `--ubam` | `false` | Input BAMs are unmapped; align with Minimap2 first (full pipeline mode only) |
| `--deepvariant_region` | `chrX` | Genomic region for DeepVariant and phase block extraction |
| `--deepvariant_model` | `ONT_R104` | DeepVariant model (`WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA`) |
| `--deepvariant_num_shards` | `24` | Number of CPUs for DeepVariant make_examples step |
| `--lrs` | `ont` | Long-read sequencing platform (`ont` or `pacbio`) |

## Credits

SkewX was originally written by Quentin Gouil, James Lancaster and Ed Yang.

We thank the following people for their extensive assistance in the development of this pipeline:

- Kathleen Zeglinski for her superior Nextflow expertise
- Shian Su for implementing new features in NanoMethViz

## Citations

If you use **SkewX** for your analysis, please cite it using the following doi: [10.1101/gr.279396.124](https://doi.org/10.1101/gr.279396.124)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
