# EBAME 2019: Longread / Nanopore sequencing Workshop
# Metagenomics of a gueuze type lambic beer from the county of Sussex, United Kingdom.


## Introduction

Today we aim to investigate the community composition of a blended (gueuze) Lambic beer brewed in the South Downs, UK. Traditional styles of beer and wine are commonly fermented using domesticated strains of yeasts (_Saccharomyces_ _cerevisiae_). While such strains show phenotypic diversity within and between different beers, wines, geographic regions and even breweries, most fermentations use a defined monoculture of _S. cervevisiae_ to give reproducible results. To maintain reproducibility, fermentation is undertaken in controlled anaerobic environments that limit exposure to contaminants from the environment which can lead to spoilage. 

Lambic beers are unlike traditional ales from the UK as the fermentation process is reliant on the natural seeding of the prepared fermentable (termed wort) by yeasts and bacteria from the local environment. Following the heating and boiling of grain to produce a range of fermentable sugars, the hot wort needs to be cooled before it is suitable for fermentation. There are many ways to achieve this under sterile conditions; however, one traditional method is to use "Coolships" which consisted of large open metal vats with a high surface to volume ratio. While highly effective in rapidly cooling the wort to fermentable temperatures, the open nature of coolships provides an opportunity for the natural microflora of the local area to enter the wort before it is sealed in barrels to undergo fermentation. The beer we are examining today has been seeded and fermented with wild microflora originating from the Firle region of the South Downs, UK. The beer itself is a limited edition run of 1500 bottles and is a blend on Lambic beer brewed in 2016 and 2017. I would like to thank brewer Gary Brandon for getting hold of the beer.  


![alt text](https://github.com/BadgerRob/Staging/blob/master/AllagashCoolship_1200.jpg "Coolship")  
_Cooling wort in a coolship_  

![alt text](https://github.com/BadgerRob/Staging/blob/master/image3.jpeg "Firle Microflora")  
_The ingredients of the beer_  

![alt text](https://github.com/BadgerRob/Staging/blob/master/image4.jpeg "beer")
_The beer_  

DNA was extracted from 250 ml of beer following centrifugation to collect a pellet consisting of cells, cellular debris and breakdown material (termed Trub). Some of this beer has been fermenting in barrels since 2016.     

![alt text](https://github.com/BadgerRob/Staging/blob/master/Trubb.jpg "Cloudy beer")
_A 50 ml falcon tube of cloudy beer_  

![alt text](https://github.com/BadgerRob/Staging/blob/master/Trubb2.jpg "Trub pellet")
_Pellet of cells and trub from 250 ml of beer. Spun at 4000 rpm @ 4 degrees C._  

Lysis was undertaken at 55 degrees C with SDS, beta mercaptoethanol and protinaseK. DNA extraction was then undertaken gently with phenol:chloroform:isoamyl alcohol, salted out with Sodium Acetate and 2.5 x 100% ethanol, then suspended in water and treated with RNAseA. DNA underwent library preparation with ONT kit LSK-109 and run on an r9.4.1 flowcell with 1562 active pores detected. DNA sequencing was undertaken using MinKNOW and preliminary basecalling was undertaken using guppy (fast mode). No barcoding/indexing was used in this sequencing run.


## Setup  

Organise into groups of 3 - 4.  

## Data  

Sample fast5 files:  
Sample quick fastq files:  
Sample high accuracy fastq files:  
minikraken database:  
workshop database:  


## Basecalling

Nanopore sequencing results in fast5 files that contain raw signal data termed "squiggles". This signal needs to be processed into the `.fastq` format for onward analysis. This is undertaken through a process called 'basecalling'. The current program released for use by Oxford Nanopore is called `Guppy` and can be implemented in both GPU and CPU modes. Two forms of basecalling are available, 'fast' and 'high-accuracy' (HAC). HAC basecalling implements a 'flipflop' basecalling algorithm which is highly computationally intensive and thus slower than the fast basecalling method. Compare the two methods on the subset of fast5 files.  

Guppy fast basecalling:
```
guppy_basecaller -r --input_path path/to/fast5/ --save_path /path/to/fastq/ --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 4 --num_callers 2

```

Guppy high_accuracy basecalling:
```
guppy_basecaller -r --input_path path/to/fast5/ --save_path /path/to/fastq/ --config dna_r9.4.1_450bps_hac.cfg  --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 4 --num_callers 2

```

|Flag                      | Description               | 
| -------------------------|:-------------------------:| 
| `guppy_basecaller`       |calls Guppy                | 
| `-r`                     |recursive mode             | 
| `--input-path`           |path to fast5 dir/         |
| `--save-path`            |path to output fastq files |
| `--qscore_filtering`     |filter for quality score   |
| `--min_qscore`           |minimum qscore for pass    |
| `cpu_threads_per_caller` |number of threads          |
| `--num_callers`          |number of basecallers      |
| `--config`               |configuration file         |



Reads are output as `.fastq` files containing 4000 reads and quality data per file. Sequences contain a `@` followed by a header then sequence. This is separated from quality data by `+`. 

(_Optional_) If you to watch in real time how many sequences are being written you can change to the directory where your fastq files are being written (/pass) and enter the bash one-liner:

```

watch -n 10 'find . -name "*.fastq" -exec grep 'read=' -c {} \; | paste -sd+ | bc'

```

|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `watch`                     |invoke 'watch' command: tells the computer to watch for something      | 
| `- n 10`                    |every 10 seconds                                                        | 
| `find .`                    |look in every directory below where you are                             |
| `-name "*.fastq"`           |target of find which will find every file with 'fastq' in its name      |
| `-exec`                     |execute the following command on find .                                 |
| `grep`                      |command to look for stuff within each file found by 'find . -name'      |
| `'read='`                   |the thing grep is looking for (1 per sequence header)                   |
| `-c`                        |count the number of 'read=' from grep                                   |
| `{} \; \| paste -sd+ \| bc'`|paste the output from grep to basic calculator and display on screen    |  

### Observations

How long did the different basecalling methods take to run?  
How do the identities differ at the individual read level when using a simple blast search of NCBI databases?  

## Read QC  

Before starting any analysis it is often advised to check the number of reads and quality of your run. You can start by using a simple bash one liner to count all reads in `pass/`. This can be done as reads are being basecalled.

Count the number of fastq reads in the Guppy pass dir.

```

cat pass/*.fastq | grep 'read=' - -c

```

|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `cat`                       |display content                                                         | 
| `pass/*.fastq`              |of all files in pass ending in .fastq                                   | 
| `\|`                        |run grep                                                                |
| `grep`                      |call grep search                                                        |
| `"read="`                   |look for lines with "read=" in                                          |
| `-`                         |target the output from `cat`                                            |
| `-c`                        |count                                                                   |


Once basecalling has completed you can create a single `.fastq` file for onward analysis. Piping outputs from `cat pass/*.fastq` can also be used if storage is limited.  

```

cat path/to/pass/*.fastq > workshop.reads.fastq

```

### Resample reads (optional extra)

If required, you can resample reads using fastqSample command from the program Canu.
To resample 15,000 reads with the same length distribution but no less than 1000bp:

```
cp path/to/workshop.reads.fastq path/to/reads.fastq.u.fastq

fastqSample -U -p 150000 -m 1000 -I /path/to/reads.fastq -O /path/to/reads.downsample.fastq.

```
Note: the file name must be `FILENAME.fastq.u.fastq` but the path must show `FILENAME.fastq`.  

|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `-U`                        |unpaired reads used for nanopore sequencing                              | 
| `-p`                        |total number of random reads to resample                                | 
| `-m`                        |minimum read length to include                                          |
| `-I`                        |path/to/input/reads.fastq                                               |
| `-O`                        |output sampled reads                                                    |
| `-max`                      |optional flag to sample longest reads                                   |



### Nanoplot


## Fixing broken fastq files with Seqkit sana (optional)

Sometimes errors can occur when preparing a `.fastq` file for analysis. This can cause problems in down stream processing. [Seqkit](https://github.com/shenwei356/seqkit) is designed to help identify errors and salvage broken `.fastq` files.

If not installed, Seqkit can be installed in your conda environment by:

```

conda install -c bioconda seqkit

```

Run [Seqkit sana](https://bioinf.shenwei.me/seqkit/usage/#sana) to sanitize a `.fastq` file.  

```

seqkit sana  workshop.reads.fastq -o rescued.workshop.reads.fastq

```

## Taxonomic identification using Kraken2.

Kraken and Kraken2 provide a means to assign taxonomic identification to reads using a k-mer based indexing against a reference database. We provide a small reference database compiled for this workshop as well as the minikraken2 database. (ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz) Other databases such as the Loman labs [microbial-fat-free](https://lomanlab.github.io/mockcommunity/mc_databases.html) and [maxikraken](https://lomanlab.github.io/mockcommunity/mc_databases.html) are also available. 

### Optional extra information

Custom reference databases can be created using `kraken2-build --download-library`, `--download-taxonomy` and `--build` [commands](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases). Mick Wattson has written [Perl scripts](https://github.com/mw55309/Kraken_db_install_scripts) to aid in customisation. An example of the creation of custom databases by the Loman lab can be found [here](http://porecamp.github.io/2017/metagenomics.html).

Run `kraken2` on the sanitized `workshop.reads.fastq` file provided in this tutorial using the `kraken2_workshop_db`. 

```

kraken2 --db path/to/kraken2_workshop_db/ --threads 8 --report path/to/output/report.txt path/to/workshop.reads.fastq > path/to/output

```

|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `kraken2`                   |call kraken2                                                            | 
| `--db`                      |database name flag                                                      | 
| `--threads`                 |number of threads to use                                                |
| `--report`                  |generate a user friendly report.txt of taxonomy                         |


### Observations

Have a look at both the output and `report.txt` files using head and more to get a first look at the sample. Use `head` and `more` bash commands.  

Does this look correct for a lambic beer?  
Anything odd in the sample?  
How much odd stuff is in the sample?  
Why do you think this has had positives hits in the kraken2 databases?  
What industry do you think the [brewers](https://www.burningskybeer.com/beers/coolship-release-no-2/) sourced the coolships from?   

Try the assembly again using the minikraken2 database and see how your database can affect your results.  

## Visualization of output

While scrolling through the kraken2 outputs can be fun and somewhat alarming, it is not the most efficient or user friendly way of eyeballing your sample. Here we present two methods to display rough community composition of the sample in a user friendly and interactive way.

### Krona

Krona produces an interactive `.html` file based on your `--report` file. While not fully integrated with kraken2, the use of the report file gives an overall approximation of your sample diversity based on individual reads. Try this on the kraken outputs from the different databases and/or basecalling modes. 

```

ktImportTaxonomy -q 2 -t 3 report.txt -o kraken_krona_report.html

```
|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `ktImportTaxonomy`          |call  KronaTools Import taxonomy                                        | 
| `q 1 -t 3`                  |for compatibility with kraken2 output                                   | 
| `report.txt.`               |Kraken2 report.txt file                                                 |
| `-o`                        |HTML output                                                             |


Copy the html files to your local machine and open in your preferred browser (tested in firefox).

```

scp USERNAME@IP:/path/to/report.txt ~/Desktop

```

It should look something like this:  
![alt text](https://github.com/BadgerRob/Staging/blob/master/Krona.png "Krona report")  


### Pavian  

[Pavian](https://github.com/fbreitwieser/pavian) is an R based program that is useful to produce Sankey plots and much more. It can be run on your local machine if you have R [installed](https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/installr.html). You may need to install `r-base-dev`. To set up Pavian in R on your local machine, open an R terminal and enter the following.  

```
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

```

To run Pavian enter into R terminal:  

```

pavian::runApp(port=5000)

```
You can now access Pavian at http://127.0.0.1:5000 in a web browser if it does not load automatically.  

Alternatively a shiny route is available.

```

shiny::runGitHub("fbreitwieser/pavian", subdir = "inst/shinyapp")

```

You should now be presented with a user interface to which you can browse for your report files on your local machine.  

![alt text](https://github.com/BadgerRob/Staging/blob/master/pavin_snap.png "Pavian input")

Once you have loaded your file, navigate to the "sample" tab and try interacting with the plot. It should look something like this:  

![alt text](https://github.com/BadgerRob/Staging/blob/master/kraken.png)

### Observations
How to the different databases affect your results?  
How does read depth affect your results?  
How does basecalling mode affect your results?  

## Assembly minimap2/miniasm/racon

### Minimap2  

[Minimap2](https://github.com/lh3/minimap2) is a program that has been developed to deal with mapping long and noisy raw nanopore reads. Two modes are used in the following assembly, `minimap2 -x ava-ont` and `minimap2 -ax map-ont`. The former performs an exhaustive "All v ALL" pairwise alignments on the read sets to find and map overlaps between reads. The latter maps long noisy read to a reference sequence. Minimap2 was developed to replace BWA for mapping long noisy reads from both nanopore and Pac-Bio sequencing runs. 


```
minimap2 -x ava-ont -t 8 workshop.reads.fastq workshop.reads.fastq | gzip -1 > workshop.paf.gz

```
Note: The output of this file is in a compressed pairwise alignment format (.paf.gz).  


### Miniasm  

[Miniasm](https://github.com/lh3/miniasm) is then used to perform an assembly using the identified overlaps in the `.paf` file. No error correction is undertaken thus the assembled contigs will have the approximate error structure of raw reads.  

```
miniasm -f workshop.reads.fastq qorkshop.paf.gz > workshop.contigs.gfa

```

Miniasm produces a graphical fragment assembly (`.gfa`) file containing assembled contigs and contig overlaps. This file type can be viewed in the program `bandage` to give a visual representation of a metagenome. However, due to the error associated with this approach, read polishing can be undertaken to increase sequence accuracy. Contigs can be extracted from a `.gfa` file and stored in a `.fasta` format using the following awk command.

```

awk '/^S/{print ">"$2"\n"$3}' workshop.contigs.gfa | fold > workshop.contigs.fasta

```

## Polishing with racon

Polishing a sequence refers to the process of identifying and correcting errors in a sequence based on a consensus or raw reads. Some methods use raw signal data from the fast5 files to aid in the identification and correction. A number of polishing programs are in circulation which include the gold standard [Nanopolish](https://github.com/jts/nanopolish), the ONT release [Medaka](https://github.com/nanoporetech/medaka) and the ultra-fast [Racon](https://github.com/isovic/racon). Each approach has advantages and disadvantages to their use. Nanopolish is computationally intensive but uses raw signal data contained in fast5 files to aid error correction. This also relies on the retention of the large `.fast5` files from a sequencing run. Medaka is reliable and relatively fast and racon is ultra fast but does not use raw squiggle data.  

The first step in polishing an assembly is to remap the raw reads back to the assembled contigs. This is done using `minimap2 -ax map-ont`.  

```
minimap2 -t 8 -ax map-ont workshop.contigs.fasta workshop.reads.fastq | gzip -1 > workshop.reads_to_assembly.paf.gz

```

Racon is then used to polish the assembled contigs using the mapped raw reads.  

```

racon -t 12 workshop.reads.fastq workshop.reads_to_assembly.paf.gz workshop.contigs.fasta > workshop.contigs.racon.fasta

```

## Kraken2 contig identification

kraken2 can be run on the assembled contigs in the same way as before, also using multiple databases for comparison.

```

kraken2 --db path/to/kraken2_workshop_db/ --threads 8 --report workshop.contigs.racon.txt workshop.contigs.racon.fasta > workshop.contigs.kraken

```

### Observations

What has happened to the number of taxa in the kraken2 report?  
How do you explain this effect?  
Try blast searching some contigs from your report against the NCBI database.  


## Flye assembly

The assemblers [Flye](https://github.com/fenderglass/Flye) and [Canu](https://github.com/marbl/canu) are available to perform assemblies which include error correction steps. Canu was primarily designed to assemble whole genomes from sequenced isolates and is more computationally intensive that Flye. Flye has a --meta flag with designed parameters to assemble long read metagenomes. Here we will run Flye on our raw reads.


```

flye --nano-raw path/to/workshop.reads.fastq -g 1g -o flye_workshop/ -t 8 --meta

```

|Flag                         | Description                                                            | 
| ----------------------------|:----------------------------------------------------------------------:| 
| `flye`                      |call  Flye                                                              | 
| `--nano-raw`                |using nanopore uncorrected raw reads as input                           | 
| `-g`                        |estimated genome size for primary coverage estimate                     |
| `-o`                        |output dir                                                              |
| `-t`                        |number of threads                                                       |
| `--meta`                    |metagenome assembly rules                                               |

Note: The assembly can now be polished with one of the formentioned programs.

### Observations

How does the fly assembly differ from the minimap2/miniasm assembly?  
How does it differ from the raw read Kraken2 report?  

## Summary

Long read sequencing provides a means to assemble metagenomes. Due to the length of reads, assemblies of larger complete contigs are possible relative to short, high accuracy reads. This often permits a greater understanding of community composition, genetic diversity as well as a greater resolution of the genetic context of genes of interest.  

So, what is _"Firle Microflora"_?
