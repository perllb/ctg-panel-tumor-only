# ctg-panel-tumor-only
Nextflow pipeline to run panels in tumor-only (somatic) mode in Dragen

- **NB!** This pipeline is not supported for automation (yet)

- Dragen Tumor-Only analysis 
- NEW: Must run ctg-parse-samplesheet (Davids script) before running. This will create one samplesheet pr project, which contain fastq1, fastq2 and bam columns in the [Data] section. 
  - You then have two alternatives: Manually join these two samplesheet into one that will be run together, or run separately for each samplesheet. (See how to specify samplesheet below). 
- The pipeline still includes demultiplexing. 
 - It supports running several projects on the same samplesheet/pipeline run (but can of course run just one project). 
 - It will do a common demux over all projects 
  - And then move fastq files to respective projects output dir. 
  - Note: In `ctg-delivery/panel-tumor-only/` there will be created a "metaid" folder, in which the demux output is written. 
    - Then, for each project, a new output dir will be created with project-id as name.
    - The metaid-folder will remain, even after then project-fastqs are moved to their respective project folders - so the Logs and Reports will be kept from the common demux.
- The pipeline can be initiated using the `panel-tumor-only-driver` - it will generate a `ctg-projects/panel-tumor-only/<metaid>` folder, with nf pipeline, samplesheet, config and bin - and start pipeline from here. 


## The following steps are performed by the pipeline:

* `Demultiplexing` (dragen bcl-conversion): Converts raw basecalls to fastq, and demultiplex samples based on index. Adapters are trimmed if added to samplesheet [Settings].
* `Alignment` (dragen): Map and align fastqs to reference genome
* `Variant calling` (dragen): Variants (SNV and SV) are called in the target region (specified with bed). CNV not yet available - will be available with Dragen 3.9, but will require panel of normal.  
* `Annotation` (nirvana): Clinical grade annotation of all variants passing basic filtering in Dragen. 
* `Dragen metrics`: Compiling Dragen alignment and coverage metrics to table.
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
* `MultiQC`: Summarizes FastQC and Dragen metrics into one document (https://multiqc.info/).


## Input files
  
- Start with the labsheet. 
1. Samplesheet (`CTG_SampleSheet.panel-tumor-only.labsheet.csv`)
- **NB!** Note that this is the samplesheet structure you need BEFORE running ctg-parse-samplesheet (Davids script).

- This will then be parse by ctg-parse-samplesheet (Davids) and create one samplesheet pr project, e.g. `CTG_SampleSheet.panel-tumor-only.2022_023.csv`
- You need to specify this when running the driver. E.g:
`panel-tumor-only-driver -s CTG_SampleSheet.panel-tumor-only.2022_023.csv`

### Samplesheet requirements

- The samplesheet (before parsing) format is standard IEM generated sheet, with the following modifications:

#### I. Added under [Header]:
```
[Header]
metaid,2022_023
PipelineName,panel-tumor-only,,,,,,,,,,,,
PipelineProfile,gmck,,,,,,,,,,,,
Paired,true,,,,,,,,,,,,
..
..
```

- Must be specified as **metaid,[project-id]** (can in theory use any name, but good to use project id. ctg-project folder will be named by the metaid)
- Comma separated, no spaces.
- Only **metaid** will work (not metaID or MetaID etc.)


#### II. Additional columns added after [Data]:

| Column | Supported values |
| ------ | -------- |
| Sample_Ref | hg38 / hg19 / mm10 : hg38, hg19 and mm10 are currently set up for dragen |
| panel | comprehensive / core / gmck : Twist-Comprehensive-Exome or Twist-Core-Exome panel bed files exist for both hg38 and hg19. GMCK ONLY with hg19!   |
| annotate | y / n : set 'y' for nirvana annotation |

- Also note that **no** Sample_Name should be added. Leave that column blank!

### Samplesheet template (.csv)

Samplesheet name: `CTG_SampleSheet.panel-tumor-only.labsheet.csv`

```
[Header]
metaid,2021_088,
PipelineName,panel-tumor-only,,,,,,,,,,,,
PipelineProfile,gmck,,,,,,,,,,,,
Paired,true,,,,,,,,,,,,
IEMFileVersion,5
Date,2021-04-29
Workflow,GenerateFASTQ
Application,NovaSeq FASTQ Only
Instrument Type,NovaSeq
Assay,TWIST
"Index Adapters,""IDT- UD Indexes (96 Indexes)"""
Chemistry,Amplicon

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index2,I5_Index_ID,index,Sample_Project,Sample_Ref,panel,annotate
OL_01_06_B_150ng,,,,3_A10,IDT_10nt_UDI_i5_202,GCAGCGTTAG,IDT_10nt_UDI_i7_202,TGTGACGTAT,2022_023,hg19,gmck,y
OL_01_08_C,,,,3_A12,IDT_10nt_UDI_i5_204,AGAGCTCTAC,IDT_10nt_UDI_i7_204,ATCCAGCATG,2022_023,hg19,gmck,y
L13_01_06_I,,,,3_A2,IDT_10nt_UDI_i5_194,GCGCTAGAGA,IDT_10nt_UDI_i7_194,CCAACTGTAT,2022_023,hg19,gmck,y
M12_02_01_C,,,,3_A3,IDT_10nt_UDI_i5_195,GCGATGCCAA,IDT_10nt_UDI_i7_195,GAGTAGTGAC,2022_023,hg19,gmck,y
L30_01_01,,,,3_A4,IDT_10nt_UDI_i5_196,ACGAAGCGAT,IDT_10nt_UDI_i7_196,CCAGTATTCT,2022_023,hg19,gmck,y 
```


## USAGE (manual run with nextflow)
Alternative is to run with driver (see below). 

1. Clone and build the Singularity container for this pipeline: (https://github.com/perllb/ctg-exome/tree/master/container). 
2. Edit the nextflow.config file to fit your project and system. Set directory from where you want to run nextflow (containing config and samplesheet) as `basedir`. (from where you execute the pipeline).
4. Edit your samplesheet to match the example samplesheet
5. Run pipeline 
```
nohup nextflow run pipe-panel-tumor-only.nf > log.pipe-panel-tumor-only.nf.txt &
```

## USAGE with driver 
For automated execution of pipeline and workflow.

- Must be started from within runfolder root directory.
- Needs:
 1. Runfolder (from where it is started)
 2. CTG_Samplesheet.. (with format as specified in ***Samplesheet requirements*** above). If not specified, driver will try to run with `runfolder/CTG_SampleSheet.panel.tumor-only.csv` from runfolder. If this does not exists, it will crash. 
   - To specify a specific samplesheet to use, run with the -s flag to specify samplesheet: `panel-tumor-only-driver -s CTG_SampleSheet.panel-tumor-only.2022_023.csv`
 3. Genome reference on dragen. Must be in /staging/$species/reference/$ref. (Species is extracted from `ref` field in SampleSheet (hg* = human, mm* = mouse)- This already exists for hg19, hg38 and mm10.
 4. Nirvana annotation. Must be set up before run, and path to nirvanadir must be added to config. Currently in `/projects/fs1/shared/references/Nirvana/`
 5. Target bed files. Defined path in nextflow.config. If using `gmck`, `core` or `comprehensive` in samplesheet "panel", the pipeline will point to the correct bed file (via nextflow.config). So as long as `Sample_Ref` is hg38, hg19 or mm10, and `panel` is one of these 3, you do not have to add any new bed files..
 
```
####################################
# PANEL SOMATIC TUMOR ONLY  driver #
####################################

Usage: panel-tumor-only-driver [ -i META_ID ] [ -s SAMPLESHEET ] [ -a INDEX-TYPE ] [ -b BCL2FASTQ-ARG ] [ -f FFPE ] [ -g DRAGEN_ARG ] [ -r RESUME ] [ -c CUSTOM-GENOME ] [ -t CUSTOM-TARGET ] [ -p PADDING ] [ -d DEMUX-OFF ] [ -h HELP ]


Optional arguments:
META-ID    -i : Set 'meta-id' for runfolder (e.g. 210330-10x). Default: Takes date of runfolder + run ID in runfolder name and adds panel-tumor-only as suffix. E.g. '210330_A00681_0334_AHWFKTDMXX' becomes 210330_0334-panel-tumor-only
SAMPLESHEET   -s : Set samplesheet used for run (Default: CTG_SampleSheet.csv)
INDEX-TYPE    -a : Set -a if single index uses. (Default: dual)
BCL2FASTQ-ARG -b : String with bcl2fastq argument. e.g. '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20'
FFPE          -f : Set -f if FFPE (then '--vc-enable-orientation-bias-filter true' will be set)
DRAGEN-ARG    -g : String with Dragen arguments. e.g. -g '--vc-enable-gatk-acceleration --vc-decoy-contigs'
RESUME        -r : Set if to resume nf-pipeline
CUSTOM-GENOME -c : Path to custom reference genome if needed. Skip if human/mouse defined in samplesheet
CUSTOM-TARGET -t : Path to custom target bed file if TWIST Comprehensive or Core is not used. Skip if core/comprehensive is defined in samplesheet
PADDING       -p : Set bp padding around target coordinates (default: 20)
DEMUX-OFF     -d : Set flag to skip mkfastq (then fastq must be in FQDIR)
HELP          -h : print help message
```

***Run driver with default settings***
This requires the current files and directories to be in correct name and location:
- `CTG_SampleSheet.panel-tumor-only.csv` in runfolder

```
cd runfolder 
panel-tumor-only-driver
```

***Run driver with single-index samples and bcl2fastq-arguments location***
```
cd runfolder 
exome-driver -a -b '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20' 
```

***Run driver with specific samplesheet (not CTG_SampleSheet.csv in runfolder)***
```
cd runfolder 
panel-tumor-only-driver -s /path/to/sheet.csv
```

## Functions of panel-tumor-only-driver

1. Creates project folder, containing:
   - nextflow.config (copied from ctg-exome pipeline dir, and edited based on default driver params and user-specified parameters)
   - pipe-exome.nf (copied from ctg-exome pipeline dir)
   - samplesheet (copied from the one specified in driver)
   - bin (scripts needed in the pipeline)
2. Creates pipeline output directory
   - default is specified in driver script (/projets/fs1/nas-sync/ctg-delivery/panel-tumor-only/<metaid>)
3. Creates QC log output directory
   - in which qc output of pipeline is copied 
4. Starts pipe-exome.nf

**Currently supports the following panels**
- `Twist Core Exome` (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- `Twist Comprehensive Exome`: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)

**Padding**
Default padding is 20bp. Can be altered in nextflow.config.

**Filters for Coverage Reports**
By default, the following filters are applied for coverage reports:
- `--remove-duplicates true`: Duplicate reads are ignored
- `--qc-coverage-ignore-overlaps true`: Resolve all of the alignments for each fragment and avoid double-counting any overlapping bases
- `--qc-coverage-filters-1 "mapq<20,bq<20"`: Reads with MAPQ<20 and BQ<20 are ignored

**Currently supports the following references**
- hg38
- mm10
- Custom references can be added in config file (set "panel" column to "custom" in samplesheet, and add customgenome=/path/to/custom/reference in nextflow.config)


## Output:
* CTG-output
    * `fastq`: Contains raw fastq files from demultiplexing.
    * `dragen`: 
      * `metrics`: Output from dragen alignment + variant calling. Contains metrics from alignment and calling, raw .vcfs and .bam files. 
         * `sv`: Output metrics and candidate vcf from sv calling.       
      * `vcf`: 
         * `*hard-filtered.vcf` contain the SNV variants (and some SVs) that passed filters. 
         * `*DiploidSV.vcf` contain SVs that passes the candidate filters. 
         * Annotated and filtered vcf are also written here.
    * `qc`: Quality control output. 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * dragen metrics: Summarized metrics for each sample.
        * multiqc output: Summarizing FastQC, dragen metrics and demultiplexing (https://multiqc.info/)



## Container
- `ngs-tools` Singularity container contain NGS-related tools, embedded in the repo: 
https://github.com/perllb/ctg-wgs/tree/master/container 

## References
- Dragen version: 3.9
    - User guide: (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf). 
    - Relase notes: (https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/Illumina-DRAGEN-Bio-IT-Platform-3.7-Release-Notes-1000000142362-v00.pdf)
- Nirvana annotator: https://illumina.github.io/NirvanaDocumentation/
- Twist Core Exome (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- Twist Comprehensive Exome: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)




