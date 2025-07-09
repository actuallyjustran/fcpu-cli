
# fcpu-cli

This repository provides a command-line interface for running GWAS (Genome-Wide Association Studies) using the high-performance [`FarmCPUpp`](https://github.com/ArtemZhou/FarmCPUpp) R package. It allows researchers to analyze genetic association using genotype and phenotype input files, supporting both `.csv` and `.vcf` formats.

## Functions
- Runs GWAS via `FarmCPUpp::farmcpu()`
- Supports the following inputs:
    - VCF (`--vcf`)
    - CSV Genotype and Map files (`--geno`, `--map`)
- Automatically converts data to `big.matrix`, a required format for FarmCPUpp
- Auto checks chromosome format
- Saves trait-specific GWAS results to CSV

## Installation

This wrapper requires the installation of a Docker program, such as [`Docker Desktop`](https://www.docker.com/products/docker-desktop/). Run the following script in CMD or WSL:

```docker build -t farmcpupp-cli .```

The required dependencies are baked into the Dockerfile.

An example of a command line for usage after building the container is shown below:

```docker run --rm -v $(pwd):/app farmcpupp-cli --pheno test-pheno.csv --vcf test.vcf```

Note: ```$(pwd)``` changes to ```"$PWD"``` when going from CMD to WSL


## Input Requirements

### Phenotype (`--pheno`)
A phenotype CSV file that:
- Has a `Taxa` column
- Has Trait columns
- Is comma delimited

Example:
```csv
Taxa,TestPheno
IRIS_313-8285,2.3
IRIS_313-8349,3.1
...
```

### Genotype inputs

#### Option 1: VCF

```bash
--vcf input.vcf
```

- Must contain valid genotype calls in GT format (0/0, 0/1, 1/1)
- Chromosome field should be numeric

#### Option 2: CSV Genotype and Map files

```bash
--geno input-geno.csv --map input-map.csv
```

- Genotype: rows = individuals, columns = SNPs
- Map file must have: `SNP`, `Chromosome`, and `Position`


## Running the software

### From VCF:
```bash
Rscript run_farmcpupp.R --pheno test-pheno.csv --vcf test.vcf
```

### From CSV genotype and map:
```bash
Rscript run_farmcpupp.R --pheno test-pheno.csv --geno test-gd.csv --map test-gm.csv
```

## 📤 Output

One result CSV per trait:
```
FarmCPUpp_<trait>.csv
```

Includes:
- SNP
- Chromosome
- Position
- p-value
- Effect estimates
- Standard error
- t-statistics


## License
This software is free and open source under the terms of the MIT License.
