# GRO-seq Pipeline

The GRO-seq analysis pipeline processes raw FASTQ data through a series of steps including adapter trimming, quality control, genome mapping, peak calling, mRNA contamination, pause index, and TSS enrichment analysis. Using Singularity for reproducibility and supporting batch analysis of multiple samples.

# Part I Workflow

![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/GRO-seq_pipeline.png)

# Part II Requirements

1. **Recommended System Configuration**:

     * 8-core CPU
     * 80 GB RAM

2. **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

     * **Step 1: Install System Dependencies**

       ```bash
       # Update package lists and install dependencies
       sudo apt-get update
       sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
            libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
       curl wget git
       ```
     
     * **Step 2: Install Go Language**
     
       ```bash
       # Download and install Go
       wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
       sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
       rm go1.21.3.linux-amd64.tar.gz
     
       # Configure Go environment variables and apply them
       echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
       echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
       source ~/.bashrc
       ```
     
     * **Step 3: Download, Build, and Install Singularity**
     
       ```bash
       # Note: The script navigates to /mnt/share/software. 
       # You can change this to your preferred directory for source code.
       cd /mnt/share/software
     
       # Download the Singularity CE source code
       wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz
     
       # Extract the archive and clean up
       tar -xvzf singularity-ce-4.0.1.tar.gz
       rm singularity-ce-4.0.1.tar.gz
       cd singularity-ce-4.0.1
     
       # Configure the build
       ./mconfig
     
       # Build Singularity (this can be time-consuming)
       cd builddir
       make
     
       # Install Singularity to the system
       sudo make install
       ```
     
     * **Step 4: Verify the Installation**
     
       ```bash
       # Check the installed version
       singularity --version
       
       # Display help information
       singularity -h
       ```
3. **snakemake**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. 

```bash
pip install snakemake
```

4. **Pipeline Files**:

     * `GRO-seq.smk`
     * `groseq.sif` (The Singularity container)

5. **Reference Data**: A directory containing all necessary reference files (e.g., bowtie2 indices and GTF annotation, etc.).

**Note on Sequencing Type:**
This pipeline supports paired-end (PE)  sequencing data. The example below shows the format for paired-end data.

### 1. Prepare the Reference Genome

The pipeline requires several pre-built reference files. Below are the steps to generate them for the human hg38 genome using GENCODE annotations.

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p data
cd data
```

#### Common Reference Files for Both Modes

We require the following base files:

**Download Reference Files:**

```bash
# Download Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz

# Download GTF Annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz

# Unzip the files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v48.primary_assembly.annotation.gtf.gz
```

##### Create rRNA Reference:

```bash
# Create a list of rRNA intervals from the GTF for bowtie2's rRNA index
awk '$3=="gene" && (/gene_type "rRNA_pseudogene"/ || /gene_type "rRNA"/) {print $1":"$4"-"$5}' gencode.v48.primary_assembly.annotation.gtf > gencode.v48.primary_assembly.rDNA.annotation.txt
singularity exec groseq.sif samtools faidx GRCh38.primary_assembly.genome.fa \
    -r gencode.v48.primary_assembly.rDNA.annotation.txt > GRCh38.primary_assembly.rDNA.fa
```

##### Build bowtie2 Indices:

```bash
# Build the main Genome Index
mkdir -p GRCh38.primary_assembly.genome.bowtie2_index
singularity exec groseq.sif bowtie2-build --threads 16 \
	GRCh38.primary_assembly.genome.fa \
	GRCh38.primary_assembly.genome.bowtie2_index/GRCh38.primary_assembly.genome

# Build the rRNA Index
mkdir -p GRCh38.primary_assembly.rDNA.bowtie2_index
singularity exec groseq.sif bowtie2-build --threads 16 \
	GRCh38.primary_assembly.rDNA.fa \
	GRCh38.primary_assembly.rDNA.bowtie2_index/GRCh38.primary_assembly.rDNA
```
### 2. Prepare the test fastq data

The data run by this pipeline is from SRR25723299 in the SRA database.The specific processing method is as follows

#### Create a dedicated directory for the test sra data:

```bash
mkdir -p data/samples
cd data/samples
```
#### Download the test sra data
```bash
prefetch SRR25723299
```
#### Convert sra data to fastq data
```bash
fastq-dump --split-files SRR25723299\SRR25723299.sra
```

#### Randomly sample fastq data:

```bash
singularity exec groseq.sif seqtk sample -s100 SRR25723299_1.fastq 27526 > SRR25723299_R1.fastq
singularity exec groseq.sif seqtk sample -s100 SRR25723299_2.fastq 27526 > SRR25723299_R2.fastq
pigz -p 8 SRR25723299_R1.fastq
pigz -p 8 SRR25723299_R2.fastq
```
### 3. Check the GRO-seq snakemake workflow

The specific snakemake workflow is as follows

#### Dry-run:

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -np \
          -s GRO-seq.smk \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/test/singularity/:/root/"
```
#### Draw a detailed DAG picture

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -s GRO-seq.smk \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/test/singularity/:/root/" \
          --dag  | \
dot -Tsvg > dag.svg
```

![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/dag.png)
### 4. Check the current working directory

The initial working structure for snakemake in /mnt/liuq/test/singularity/:
## Input Structure and Interpretation

Before the pipeline, the input directory contain several files and directories. Below is a detailed explanation of what each file you should supply and how it can be used.
```bash
├── config.yaml
├── dag.svg
├── data
│   ├── gencode.v48.primary_assembly.annotation.gtf
│   ├── gencode.v48.primary_assembly_exons_no_first_exon.bed
│   ├── gencode.v48.primary_assembly_introns.bed
│   ├── GRCh38.primary_assembly.genome.bowtie2_index
│   │   ├── GRCh38.primary_assembly.genome.1.bt2
│   │   ├── GRCh38.primary_assembly.genome.2.bt2
│   │   ├── GRCh38.primary_assembly.genome.3.bt2
│   │   ├── GRCh38.primary_assembly.genome.4.bt2
│   │   ├── GRCh38.primary_assembly.genome.rev.1.bt2
│   │   └── GRCh38.primary_assembly.genome.rev.2.bt2
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── GRCh38.primary_assembly.rDNA.bowtie2_index
│   │   ├── GRCh38.primary_assembly.rDNA.1.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.2.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.3.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.4.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.rev.1.bt2
│   │   └── GRCh38.primary_assembly.rDNA.rev.2.bt2
│   ├── minus_TSS.tsv
│   ├── plus_TSS.tsv
│   └── samples
│       ├── SRR25723299_R1.fastq.gz
│       └── SRR25723299_R2.fastq.gz
├── GRO-seq.smk
├── groseq.sif
└── scripts
    ├── filter_paired_fq.pl
    ├── pause_index.py
    ├── PEPPRO.R
    ├── plot_pause_index.r
    ├── plot_tss_enrichment.R
    └── pyTssEnrichment.py

```
- **GRO-seq.smk** — The main Snakemake workflow script.
- **config.yaml** — Configuration file containing paths, parameters, and sample information.
  ⚠️ Must be located in the same directory as `GRO-seq.smk`.
- **proseq.sif** — Singularity container image with all required software and dependencies pre-installed.
- **dag.svg**— Detailed pipeline DAG picture.
- **scripts** — Additional scripts required for program execution.
- **data** — Reference genome index and rDNA index for bowtie2,samples which were sequenced and other required files which are listed as follows.
#### Per-Sample Files (`data/samples/`)

- **`plus_TSS.tsv`**
  - **Content**: Browser Extensible Data (BED) format file used to describe the TSS position of the **coding** strand of the chromosome.Its fourth column is the transcript ID, and the fifth column is occupied by 0.
  
  - **Application**: It's the TSS position of the **coding** strand of the chromosome and can be used for TssEnrichment rule in snakemake workflow.
  
  - **Example**
  
```bash
chr1    11121   11122   ENST00000832824.1       0       +
chr1    11125   11126   ENST00000832825.1       0       +
chr1    11410   11411   ENST00000832826.1       0       +
chr1    11411   11412   ENST00000832827.1       0       +
chr1    11426   11427   ENST00000832828.1       0       +
```

- **`minus_TSS.tsv`**
  
  - **Content**: Browser Extensible Data (BED) format file used to describe the TSS position of the **template** strand of the chromosome.Its fourth column is the transcript ID, and the fifth column is occupied by 0.
  
  - **Application**: It's the TSS position of the **template** strand of the chromosome and can be used for TssEnrichment rule in snakemake workflow.
  
  - **Example**
  
```bash
chr1    15104   15105   ENST00000831746.1       0       -
chr1    29344   29345   ENST00000831201.1       0       -
chr1    29359   29360   ENST00000831165.1       0       -
chr1    29365   29366   ENST00000831154.1       0       -
chr1    29365   29366   ENST00000831157.1       0       -
```

- **`rDNA.bowtie2_index`**
  
  - **Content**: This is the rDNA index created by bowtie2,The building method is as above.
  - **Application**: It's the primary reference for rDNA read alignment.
- **`genome.fa`**
  - **Content**: This is the genome fasta file.
  - **Application**: It is used for pause_index rule in snakemake workflow.
- **`genome.bowtie2_index`**
  - **Content**:  This is the whole genome index created by bowtie2,The building method is as above.
  - **Application**:  It's the primary reference for whole genome read alignment.
- **`introns.bed`**
  - **Content**: Browser Extensible Data (BED) format file used to describe the introns position of the chromosome.Its fourth column is the transcript ID, and the fifth column is the gene ID.
  - **Application**: It's the introns position of the chromosome and can be used for exons_introns_sort rule in snakemake workflow.
```bash
chr1    11211   12010   ENST00000832824.1       ENSG00000290825.2       +
chr1    12227   12613   ENST00000832824.1       ENSG00000290825.2       +
chr1    12721   13453   ENST00000832824.1       ENSG00000290825.2       +
chr1    11211   12010   ENST00000832825.1       ENSG00000290825.2       +
chr1    12227   12613   ENST00000832825.1       ENSG00000290825.2       +
```
- **`exons_no_first_exon.bed`**
  - **Content**: Browser Extensible Data (BED) format file used to describe the exons (do not include the first exon of each transcript) position of the chromosome.Its fourth column is the transcript ID, and the fifth column is the gene ID.
  - **Application**:  It's the exons (do not include the first exon of each transcript) position of the chromosome and can be used for exons_introns_sort rule in snakemake workflow.
```bash
chr1    12010   12227   ENST00000832824.1       ENSG00000290825.2       +
chr1    12613   12721   ENST00000832824.1       ENSG00000290825.2       +
chr1    13453   14413   ENST00000832824.1       ENSG00000290825.2       +
chr1    12010   12227   ENST00000832825.1       ENSG00000290825.2       +
chr1    12613   12721   ENST00000832825.1       ENSG00000290825.2       +
```
- **`annotation.gtf`**
  - **Content**: This is the genome annotation gtf file.
  - **Application**: It is used for pause_index rule in snakemake workflow.
# Part III Running

After adding the corresponding parameters in config.yaml,you can execute the pipeline using a single command, the only important parameter is thread (**-j**).
**Example code**

* **Step 1: Edit `config.yaml`**

```shell
samples:
    SRR25723299: 
        R1 : "data/samples/SRR25723299_R1.fastq.gz"
        R2 : "data/samples/SRR25723299_R2.fastq.gz"
fn_gtf : "data/gencode.v48.primary_assembly.annotation.gtf"
fa_in : 'data/GRCh38.primary_assembly.genome.fa'
exon_name : "data/gencode.v48.primary_assembly_exons_no_first_exon.bed"
intron_name : "data/gencode.v48.primary_assembly_introns.bed"
BOWTIE2_RDNAIDX : "data/GRCh38.primary_assembly.rDNA.bowtie2_index/GRCh38.primary_assembly.rDNA"
BOWTIE2IDX : "data/GRCh38.primary_assembly.genome.bowtie2_index/GRCh38.primary_assembly.genome"
home : "/mnt/liuq/test/singularity"
container : "proseq.sif"
species : "hsa"
```

* **Step 2: run snakemake**

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -s GRO-seq.smk \
		  -j 25 \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/test/singularity/:/root/"
```

**Command Parameters**

**edit `config.yaml`**

- `samples`: Path to the FASTQ file. For paired-end data (required)
- `BOWTIE2_RDNAIDX`: Path to the directory where bowtie2 reference build with prefix (required)
- `BOWTIE2IDX`: Path to the directory where bowtie2 rDNA reference build with prefix (required)
- `container`: Path to the singularity environment file (required)
- `exon_name`: Browser Extensible Data (BED) format file used to describe the exons (do not include the first exon of each transcript) position of the chromosome.Its fourth column is the transcript ID, and the fifth column is the gene ID (required)
- `intron_name`: Browser Extensible Data (BED) format file used to describe the introns position of the chromosome.Its fourth column is the transcript ID, and the fifth column is the gene ID (required)
- `fa_in` : Reference Genome FASTA (required)
- `fn_gtf` : Reference Genome GTF Annotation (required)
- `home`: The root directory where running program (required)
- `species`: species which was sequenced (required)

**run snakemake**

- `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
- `--singularity-args`: Allows passing additional arguments to the Singularity runtime (e.g., `--bind`, `--nv`, or custom options).
- `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel when executing workflow rules.
- `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references. The format is `/project_directory:/project_directory`. Multiple directories can be mounted by separating them with commas, for example: `/path1:/path1,/path2:/path2` (required)

Then delete the intermediate files and folders

```bash
rm -rf bed/ qc/ temp/ intermediate/*/HOMER_tag data/multiqc_report/*/multiqc_data
rm intermediate/count_tts.filtered.txt intermediate/count_tts.raw.txt intermediate/tts_down_50k.txt \
intermediate/*tss_tts.txt intermediate/*tss.txt \
data/*_fail_qc_dups.bam data/prealignments_*.bam \
data/*_temp.bam data/*_sort.bam data/*_noMT.bam intermediate/*_chr_order.txt data/*_noMT.bam.bai data/*_PE1.bam data/*_PE1.bam.bai
```
## Output Structure and Interpretation

After the pipeline completes, the output directory will contain several files and directories. Below is a detailed explanation of what each file is and how it can be used.
```bash
├── config.yaml
├── dag.svg
├── data
│   ├── gencode.v48.primary_assembly.annotation.gtf
│   ├── gencode.v48.primary_assembly_exons_no_first_exon.bed
│   ├── gencode.v48.primary_assembly_introns.bed
│   ├── GRCh38.primary_assembly.genome.bowtie2_index
│   │   ├── GRCh38.primary_assembly.genome.1.bt2
│   │   ├── GRCh38.primary_assembly.genome.2.bt2
│   │   ├── GRCh38.primary_assembly.genome.3.bt2
│   │   ├── GRCh38.primary_assembly.genome.4.bt2
│   │   ├── GRCh38.primary_assembly.genome.rev.1.bt2
│   │   └── GRCh38.primary_assembly.genome.rev.2.bt2
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── GRCh38.primary_assembly.rDNA.bowtie2_index
│   │   ├── GRCh38.primary_assembly.rDNA.1.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.2.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.3.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.4.bt2
│   │   ├── GRCh38.primary_assembly.rDNA.rev.1.bt2
│   │   └── GRCh38.primary_assembly.rDNA.rev.2.bt2
│   ├── minus_TSS.tsv
│   ├── multiqc_report
│   │   └── SRR25723299
│   │       └── multiqc_report.html
│   ├── plus_TSS.tsv
│   ├── samples
│   │   ├── SRR25723299_R1.fastq.gz
│   │   └── SRR25723299_R2.fastq.gz
│   ├── SRR25723299
│   │   ├── SRR25723299_R1_unmap_filter_paired_fastqc.html
│   │   ├── SRR25723299_R1_unmap_filter_paired_fastqc.zip
│   │   ├── SRR25723299_R2_unmap_filter_paired_fastqc.html
│   │   └── SRR25723299_R2_unmap_filter_paired_fastqc.zip
│   ├── SRR25723299_minus.bw
│   └── SRR25723299_plus.bw
├── figure
│   ├── SRR25723299_mRNA_contamination.pdf
│   ├── SRR25723299_mRNA_contamination.png
│   ├── SRR25723299_pause_index.pdf
│   ├── SRR25723299_pause_index.png
│   ├── SRR25723299_TSSenrichment.pdf
│   └── SRR25723299_TSSenrichment.png
├── GRO-seq.smk
├── intermediate
│   └── SRR25723299
│       ├── transcript.gtf
│       └── transcript.txt
├── groseq.sif
├── scripts
│   ├── filter_paired_fq.pl
│   ├── pause_index.py
│   ├── PEPPRO.R
│   ├── plot_pause_index.r
│   ├── plot_tss_enrichment.R
│   └── pyTssEnrichment.py
└── stat
    └── SRR25723299_align_with_bt2.txt
```

- **`multiqc_report.html`**
  Open multiqc_report.html in a web browser to explore all sections interactively.
  
   - **Application**: This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .
  
       - **General Statistics**: A combined table summarizing important metrics for each sample:
       
       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/general_statistic.png)
       - **FastQC**: Quality-control metrics on raw and trimmed reads, including  
        'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
           'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
           'Sequence Length Distribution', 'Sequence Duplication Levels',  
           'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.
         
       - **Sequence Quality Histograms**: The mean quality value across each base position in the read. 
       
       
       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/fastqc_per_base_sequence_quality_plot.png)
           
      - **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.  
      
       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/fastqc_per_sequence_quality_scores_plot.png)
      
  
- **`fastqc.html(zip)`**
  
  - **Content**: Open fastqc.html in a web browser to explore all sections interactively.Similar to the above multiqc results.
  
  - **Application**:This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .Similar to the above multiqc results.
  
- **`minus.bw & plus.bw`**
  - **Content**: Two BigWig files that represent the PRO-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format.

  - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes. Highly expressed genes will appear as peaks.

       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/UCSC_bw.png)
  
- **`mRNA_contamination.pdf(png)`**
  - **Content**: For Assessment of nascent RNA purity,calculates the exon to intron read density ratio. And PRO-seq libraries have an increased promoter emphasis and higher mRNA contamination indicated by an increase in reads in promoters and exons at the cost of reads in introns and promoter flanking regions ..

  - **Application**: It is used to check the contamination of mRNA, the peak value is required to be less than 2.
  
       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/SRR25723299_mRNA_contamination.png)


- **`pause_index.pdf(png)`**
  - **Content**: Define the pause index as the ratio of the density of reads in the *pausing* region versus the density in the corresponding gene body,then plots the frequency distribution of the pause index across genes. A greater pause index indicates a more efficient run-on, as a higher value indicates that paused polymerases efficiently incorporate the modified NTPs. An efficient run-on process has a median pause index greater than 10.
  
  - **Application**: It is used to check the transcription pause, the peak value is required to be more than 10.
    
       ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/SRR25723299_pause_index.png)
  
- **`TSSenrichment.pdf(png)`**
  
  - **Content**: For Assessment of run-on efficiency, aggregates sequencing reads at 2000 bases upstream and downstream of a reference set of TSSs to plot and calculate a TSS enrichment score. The normalized TSS enrichment score is the ratio of the average coverage in 100-bp windows, with the numerator centered at the TSS peak summit and the denominator in the background at the edge of the 2000-bp window.
  
  - **Application**:  It is used to check the TSS enrichment,requires the negative strand to precede the positive strand.
  
    ![](https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/picture/SRR25723299_TSSenrichment.png)
  
  
  
- **`transcript.gtf(txt)`**
  - **Content**:  Indicates the location information of transcripts that meet PRO-seq.
  - **Application**:  It is used to find transcriptional regulatory sequences


- **`align_with_bt2.txt`**
  - **Content**: Bowtie2 log of alignment to genome.
  - **Application**:  It used to check the alignment rate, unique alignment rate and other information.

```bash
20746 reads; of these:
  20746 (100.00%) were paired; of these:
    9951 (47.97%) aligned concordantly 0 times
    7148 (34.45%) aligned concordantly exactly 1 time
    3647 (17.58%) aligned concordantly >1 times
    ----
    9951 pairs aligned concordantly 0 times; of these:
      708 (7.11%) aligned discordantly 1 time
    ----
    9243 pairs aligned 0 times concordantly or discordantly; of these:
      18486 mates make up the pairs; of these:
        14014 (75.81%) aligned 0 times
        2449 (13.25%) aligned exactly 1 time
        2023 (10.94%) aligned >1 times
66.22% overall alignment rate
```

Here, the alignment rate is 66.22% and the unique alignment rate is 34.45%.

## Video Tutorials

Watch this video tutorial to see a complete walk through of running the pipeline:

https://github.com/Haolab-BIG/GRO-seq-Processing-Pipeline/raw/main/snakemake.mp4

