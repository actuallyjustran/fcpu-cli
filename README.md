
# fcpu-cli

This repository provides a command-line interface for running Genome-Wide Association Studies (GWAS) using the [`FarmCPUpp`](https://github.com/amkusmec/FarmCPUpp) package. It allows researchers to analyze genetic association using genotype and phenotype input files, supporting both `.csv`, `.tsv`, `.vcf`, and `.vcf.gz` formats.

## Functions
- Runs GWAS via `FarmCPUpp::farmcpu()`
- Supports the following inputs:
    - VCF (`--vcf`)
    - CSV Genotype and Map files (`--geno`, `--map`)
- Automatically converts data to `big.matrix`, a required format for FarmCPUpp
- Auto checks chromosome format
- Saves trait-specific GWAS results to CSV
- Generates Manhattan and QQ plots

## Installation

This wrapper requires the installation of a Docker program, such as [`Docker Desktop`](https://www.docker.com/products/docker-desktop/). Run the following script in CMD or WSL:

```docker build -t fcpu-cli .```

The required dependencies are baked into the Dockerfile.

An example of a command line for usage after building the container is shown below:

```CLI
docker run --rm -v "$PWD":/usr/src/app -w /usr/src/app farmcpupp Rscript run_farmcpu.R --vcf data/test-150.vcf.gz --pheno data/test-pheno.csv
```

Note: ```$(pwd)``` changes to ```"$PWD"``` when going from CMD to WSL


## Input

### Phenotype (`--pheno`)
A phenotype CSV file that:
- Has a `Taxa` column
- Has Trait columns
- Is tab-delimited (`.tsv`) or comma-delimited (`.csv`)

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
- .vcf.gz is also an accepted file

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

## Output

- `FarmCPUpp_<trait>.csv` – GWAS results per trait
- `FarmCPUpp_<trait>_manhattan.png` – Manhattan plot of p-values
- `FarmCPUpp_<trait>_qq.png` – QQ plot of observed vs expected p-values
<p align="center">
  <img src="https://github.com/actuallyjustran/fcpu-cli/blob/cd41c80f3e5cb6238d22c67c5c2e910b5c0b8831/results/FarmCPUpp_TestPheno_manhattan.png?raw=true" alt="Manhattan plot" width="600"/>
  <br>
  <em>Generated Manhattan plot</em>
</p>

<p align="center">
  <img src="https://github.com/actuallyjustran/fcpu-cli/blob/cd41c80f3e5cb6238d22c67c5c2e910b5c0b8831/results/FarmCPUpp_TestPheno_qq.png?raw=true" alt="QQ plot" width="600"/>
  <br>
  <em>Generated QQ plot</em>
</p>


If the GWAS model produces valid p-values, plots will be generated automatically.

GWAS Results include:
- SNP
- Chromosome
- Position
- p-value
- Effect estimates
- Standard error
- t-statistics


## License
This software is free and open source under the terms of the MIT License.
