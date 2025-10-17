# ampliseq: Documentation

## Quick Start

To run the pipeline with a test dataset, execute the following command:

```bash
nextflow run ampliseq -profile quick_test,docker --outdir output_test
```

This will run the pipeline with a minimal test dataset and the docker profile.

## Pipeline Steps

The ampliseq pipeline processes raw sequencing data to produce an Amplicon Sequence Variant (ASV) table and taxonomic assignments. The main steps are as follows:

1.  **Input Validation:** The pipeline validates the input samplesheet and read files.
2.  **Adapter and Primer Trimming:** `Cutadapt` is used to remove adapter sequences and PCR primers from the raw reads.
3.  **Read Filtering and Denoising:** The pipeline uses `DADA2` for:
    *   Filtering low-quality reads.
    *   Learning sequencing error rates.
    *   Denoising reads into Amplicon Sequence Variants (ASVs).
    *   Merging paired-end reads.
    *   Removing chimeric sequences.
4.  **ASV Processing:** Optional steps include clustering similar ASVs with `VSEARCH` and filtering ASVs by length.
5.  **Taxonomic Classification:** Assigns taxonomy to the final ASVs (details in the [Usage](usage.md) documentation).
6.  **Reporting:** The pipeline generates a `MultiQC` report summarizing the results from all steps.

For more detailed information, please refer to the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)