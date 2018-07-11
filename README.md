# virus-detection-nf
Detection of presence of viral genome integrations using RNA-seq data.

## Installation

### Nextflow
Install `nextflow` following the [instructions](https://www.nextflow.io/docs/latest/getstarted.html).

Be sure to run at least Nextflow version 0.30.2.

### Singularity
Install `singularity` following the instructions at
https://singularity.lbl.gov/install-linux

### virus-detection-nf pipeline

The most convenient way is to install `virus-detection-nf` is to use `nextflow`'s built-in `pull` command
```bash
nextflow pull obenauflab/virus-detection-nf
```

## Documentation

The workflow consists of 3 steps:

* Centrifuge: Detection and quantification of presence of viral sequences in the overall dataset
* BWA + Manta: Alignment of the data to the human genome + detected viral genomes to determine integration sites
* Sailfish: Collect and quantify splice-variants in the dataset. 

```bash
nextflow run obenauflab/virus-detection-nf --help
```

## Credits
[Nextflow](https://github.com/nextflow-io/nextflow):  Paolo Di Tommaso

[Singularity](https://singularity.lbl.gov): Singularityware

[Centrifuge](https://ccb.jhu.edu/software/centrifuge/): Daehwan Kim, Salzberg Lab

[bwa](http://bio-bwa.sourceforge.net/): Hl3

[Sailfish](https://github.com/kingsfordgroup/sailfish) - Kingsford Group
