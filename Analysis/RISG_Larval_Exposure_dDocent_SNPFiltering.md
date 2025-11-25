# RISG Diel-Cycling CADO Larval Exposure Experiment 
Author: A. Zyck
Date: August 2024

Bioinformatic analysis of EecSeq data from RISG larval stress exposure experiment that was conducted in summer 2022. Samples were processed for sequencing using the EecSeq protocol from summer 2023 - winter 2024. 

In the exposure experiment, there were three replicate spawning blocks, labeled B2, B3, B4. Before exposure, four samples were collected for genomic processing (T0). Larvae were reared in a ambient (Con) treatment or a stress (Hi) treatment. Ambient conditions were filtered seawater from Narragansett Bay. The stress treatment was a diel-cycling of pH and dissolved oxygen (DO) with ambient conditions during the day, cycling down to a pH of 7.0 and DO of 1.5 mg/L over night. There were three replicates per treatment. At the end of each exposure period, larval samples were collected for each bucket. 

All bioinformatic analyses will be performed in KITT (PLOMEE Lab server).

Data uploaded and analyzed on KITT. User logged in before following steps are completed.

In home directory, make new project directory. 

```
$ mkdir RISG_Larval

$ cd RISG_Larval
```

**Data location: `PATH /home/azyck/RISG_Larval`**

## In Terminal:

### 1. Setup: Downloading software programs, creating environments anaconda folders

**Downloading Bioconda (skip this step if Bioconda has already been downloaded).**

```javascript
# downloading Miniconda software
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh

# restarting window with source command
$ source ~/.bashrc

# adding different channels
$ conda config --add channels defaults
	# should get: "warning: 'defaults' already in 'channels' list, moving to the top"
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

# to see if this worked with cat command
$ cat .condarc
```
> Bioconda is a bioinformatics software package manager. See more at [https://bioconda.github.io](https:///bioconda.github.io).

**Create and activate a dDocent conda environment:**

```
$ conda create -n risg_larv ddocent
$ conda activate risg_larv

# the beginning of the line should then look like: (risg_larv) [username@KITT ~]$
```

Copy sequencing files from storage in KITT to this working directory. Code shown for B2 files, repeated for B3 and B4. 

```
$ cp /RAID_STORAGE4/Shared/RISG_CADO/demultiplexed/larval/B2_T0* .

$ cp /RAID_STORAGE4/Shared/RISG_CADO/demultiplexed/larval/B2_TCE* .
```

**Checking the quality of data post-demultiplexing**

```
$ mkdir demultiplexed_fastqc_results
$ cd demultiplexed_fastqc_results
$ fastqc ../*fq.gz
$ mv ../*fastqc.* .
```

**Multiqc Analysis**

Note: Mutliqc won't run for me in my conda environment, but runs in base. 

```
$ multiqc .

$ mv multiqc_report.html demultiplexed_multiqc_report.html
```

I pushed the demultiplexed multiqc report to Github, in Output directory. The report can be viewed [here](https://github.com/amyzyck/RISG_LarvalCADO/blob/main/Output/demultiplexed_multiqc_report.html).

Overall the samoles look pretty good. A few samples have a lower number of sequences that I would like, but we'll move on. 

### Read Trimming, Mapping, and SNP Calling

Using [**dDocent**](http://www.ddocent.com/)

> dDocent is a bioinformatics program created by Dr. Jon Purtiz that is specifically designed for different types of RAD sequencing.

**Jon downloaded a version of dDocent on my KITT account `dDocent_ngs` that can be used for Expressed Exome Capture Sequencing (EecSeq) and pooled samples (larval pools). It is located in the `RISG_Larval` directory**

If this is your first time running dDocent, I recommend going through the [Quick Start Guide](http://www.ddocent.com/quick/). I also recommend:

1. Reading through the [User Guide](http://www.ddocent.com/UserGuide/).
2. Completing the [Assembly Tutorial](http://www.ddocent.com/assembly/), using the simulated dataset.

**Create and activate a dDocent conda environment (if you did not do so previously):**

```
$ conda create -n risg_larv ddocent
$ conda activate risg_larv

# the beginning of the line should then look like: (risg_larv) [azyck@KITT ~]$
```

**Make directory for dDocent and link files and `dDocent_ngs` into directory**

Starting in `RISG_Larval` directory:

```
$ mkdir RISG_ddocent
$ cd RISG_ddocent/

$ ln -s ../*.fq.gz .
$ ln -s ../dDocent_ngs .
```


**Running dDocent for trimming, mapping, and SNP calling**

I previously ran each step at a time with checkpoints in between. The read trimming has worked well. I determined optimal values for match score, mismatch score, and gap penalty for read mapping (A = 2, B = 4, and O = 6). Now I'm running dDocent will all steps at once.  

```
$ bash dDocent_ngs
```

```
dDocent 3.1.0 

Contact jpuritz@uri.edu with any problems 

 
Checking for required software

All required software is installed!

dDocent version 3.1.0 started Fri Aug 16 21:17:12 EDT 2024 

30 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 30 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
2
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
4
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
6
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
yes
Is this a pooled sequencing data set?  Please type yes or no and press [ENTER]
yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.

dDocent started Fri Aug 16 21:17:12 EDT 2024

dDocent finished Mon Aug 19 21:08:47 EDT 2024

After filtering, kept 33920239 out of a possible 48286413 Sites
```

```
$ mv TotalRawSNPs.vcf.gz raw.total.vcf.gz 
```

### Variant Filtering 

In the working directory, make a new directory for filtering

In `PATH:/home/azyck/RISG_Larval/RISG_ddocent`

```
$ mkdir RISG_Larval_SNPFiltering
$ cd RISG_Larval_SNPFiltering
```

Link vcf file to this directory:

```
$ ln -s ../raw.total.vcf.gz .
```

Because the larval samples are pooled, any filtering steps applied are not based on genotype information. 

I'll be following steps from a scipt that Jon shared with me for filtering pooled larval samples. I couldn't get the original script to run properly, so I ran the steps from the script one-by-one.  

This first step comes directly from the script Jon shared. It's set up to work with multiple vcf files, but can work with just one. Here are the specfic steps it's doing:

- Starts with bcftools view to read the gzipped VCF file raw.total.vcf.gz.
- Applies the first filter to keep variants where less than 75% of samples have missing data.
- Uses bcftools +setGT to set genotypes to missing (.) for samples where the depth is less than 20
- Applies the final set of filters:
   - Keeps variants where less than 25% of samples have missing data
   - Keeps variants where the ratio of alternate allele observations (AO) to reference allele observations (RO) is greater than 0.015
- Outputs the result as a compressed VCF file
- Redirects error messages to Larv.filter.errors.

```
$ bcftools view -i 'F_MISSING<0.75' raw.total.vcf.gz | \
> bcftools +setGT -- -t q -n . -i "FORMAT/DP<20" 2> RISG.filter.errors | \
> bcftools view -i 'F_MISSING<0.25 && INFO/AO/INFO/RO > 0.015' -O z -o RISG.TRSdp.20.g75.total.recode.vcf.gz 2>> RISG.filter.errors 

$ zcat RISG.TRSdp.20.g75.total.recode.vcf.gz | grep -v '^#' | wc -l

output: 
2035486
```

Next , I'm splitting up the variants from mitochondrial and nuclear DNA. 

```
$ zcat RISG.TRSdp.20.g75.total.recode.vcf.gz | mawk '!/NC_007175.2/' > RISG.TRSdp.20.g75.nDNA.vcf

$ cat RISG.TRSdp.20.g75.nDNA.vcf | grep -v '^#' | wc -l

output:
2035199
```

The rest of the filtering steps from the script are saved in a separate script [`dDocent_ngs_filters2`](https://github.com/amyzyck/RISG_LarvJuvCADO/blob/main/Scripts/larval_scripts/dDocent_ngs_filters2). This is modified from another original script Jon made [`dDocent_ngs_filters`](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Scripts/dDocent_ngs_filters). This script filters the remaining sites based on properly paired status and quality depth ratio. 

```
$ curl -L -O https://raw.githubusercontent.com/amyzyck/RISG_LarvJuvCADO/refs/heads/main/Scripts/larval_scripts/dDocent_ngs_filters2
$ chmod +x dDocent_ngs_filters2

$ ./dDocent_ngs_filters2 RISG.TRSdp.20.g75.nDNA.vcf RISGTRSdp20g75nDNA

output:
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed. 

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters 

Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of additional sites filtered based on properly paired status
 7620 of 2035199 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 298905 of 2035199 


                                                                                                                        
                                                                                                                        
                                                Histogram of mean depth per site                                        
     100000 +-------------------------------------------------------------------------------------------------------+   
            |     +     +    +    ***    +     +    +     +     +     +     +    +     +     +     +    +     +     |   
      90000 |-+                  *****                    'meandepthpersite' using (bin($1,binwidth)):(1.0) *******-|   
            |                    *****                                                                              |   
            |                    *******                                                                            |   
      80000 |-+                 ****** **                                                                         +-|   
            |                   ****** ***                                                                          |   
      70000 |-+                 ****** ***                                                                        +-|   
            |                   ****** ****                                                                         |   
      60000 |-+               ******** *****                                                                      +-|   
            |                 * ****** ******                                                                       |   
      50000 |-+               * ****** ********                                                                   +-|   
            |                 * ****** ****** **                                                                    |   
            |                 * ****** ****** ***                                                                   |   
      40000 |-+               * ****** ****** *****                                                               +-|   
            |                ** ****** ****** *******                                                               |   
      30000 |-+              ** ****** ****** **********                                                          +-|   
            |                ** ****** ****** ******* ****                                                          |   
      20000 |-+              ** ****** ****** ******* *********                                                   +-|   
            |               *** ****** ****** ******* ****** *******                                                |   
            |               *** ****** ****** ******* ****** *********                                              |   
      10000 |-+             *** ****** ****** ******* ****** ******* *                                            +-|   
            |     +     +  **** ****** ****** ******* ****** ******* *******************     +     +    +     +     |   
          0 +-------------------------------------------------------------------------------------------------------+   
            10    15    20   25    30    35    40   45    50    55    60    65   70    75    80    85   90    95   100  
                                                           Mean Depth                                                   
                                                                                                                        
The 95% cutoff would be 80
Would you like to use a different maximum mean depth cutoff than 80, yes or no
yes
Please enter new cutoff
70
Number of sites filtered based on maximum mean depth
 109826 of 1728798 

Total number of sites filtered
 416227 of 2035199 

Remaining sites
 1618972 

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in RISGTRSdp20g75nDNA.filterstats
```


Next, I'm going to split this VCF up into a separate vcf for each Block. In Jon's CASE analysis, he found distinct allele frequency shift patterns specific to Block. 


First I'm going to pull the sample names out of the vcf file

If the Output_prefix.FIL.recode.vcf is gzipped, make sure it's gzipped with `bgzip`

```
$ tabix -p vcf RISGTRSdp20g75nDNA.FIL.recode.vcf.gz
$ bcftools query -l RISGTRSdp20g75nDNA.FIL.recode.vcf.gz > samples
```

I'm going to apply a few additional filtering steps, following what Jon did. This additional filtering will allow up to 10% missing data, keeps only bi-allelic sites, calculates allele frequencies, and then filters based on allele frequency (0.015 < AAF < 0.985)

First Block 2: 

```
$ bcftools view -S <(grep '^B2' samples) RISGTRSdp20g75nDNA.FIL.recode.vcf.gz | bcftools view -i 'F_MISSING<0.1'  | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B2.RISG.vcf.gz

$ bcftools index B2.RISG.vcf.gz

$ bcftools query -f '%CHROM\t%POS\n' B2.RISG.vcf.gz > B2.pos

$ cat B2.pos | wc -l
702163
```

Block 3: 

```
$ bcftools view -S <(grep '^B3' samples) RISGTRSdp20g75nDNA.FIL.recode.vcf.gz | bcftools view -i 'F_MISSING<0.1'  | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B3.RISG.vcf.gz

$ bcftools index B3.RISG.vcf.gz

$ bcftools query -f '%CHROM\t%POS\n' B3.RISG.vcf.gz > B3.pos

$ cat B3.pos | wc -l
144949
```

Block 4: 

```
$ bcftools view -S <(grep '^B4' samples) RISGTRSdp20g75nDNA.FIL.recode.vcf.gz | bcftools view -i 'F_MISSING<0.1'  | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B4.RISG.vcf.gz

$ bcftools index B4.RISG.vcf.gz

$ bcftools query -f '%CHROM\t%POS\n' B4.RISG.vcf.gz > B4.pos

$ cat B4.pos | wc -l
732518
```

Combine all .pos files, sorts them, and removes duplicates.

```
$ cat B*.pos | sort | uniq > fil.pos

$ cat fil.pos | wc -l
1053839
```

Create a new VCF file (SNP.RISG.TRSdp.20.B90.2a.perp.vcf.gz) containing only the SNPs present in fil.pos.

```
$ bcftools view -R fil.pos --threads 40 -m2 -M 4 RISGTRSdp20g75nDNA.FIL.recode.vcf.gz -O z -o SNP.RISG.TRSdp.20.B90.2a.perp.vcf.gz
```

I'm going to an extra step where I extract SNP positions from the final combined VCF file and convert them to BED format. I'll then use bedtools to intersect these positions with a gene annotation file (I'll ask Jon for this).
From here, I can extract gene names (LOC IDs) and create a unique, sorted list in `RISG.study.background.LOC`.

```
$ bcftools view --threads 40 SNP.RISG.TRSdp.20.B90.2a.perp.vcf.gz | mawk '!/#/' | cut -f 1,2 | mawk '{print $1"\t"$2-1"\t"$2}' > total.snp.bed

# I will have to do this step later 
$ bedtools intersect -wb -a total.snp.bed -b sorted.ref3.0.gene.sc.hmask.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > RISG.study.background.LOC
```

#### Convert VCF to PoPoolation2 (SYNC files)

The final VCF file needs to converted to a different format in order for it to work in PoPoolation2. Jon gave me a script to convert the file, called `VCFtoPopPool.py`. 

```
$ bcftools view --threads 40 B2.RISG.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
$ python2 VCFtoPopPool.py temp.vcf RISG.Block2.sync 
$ rm temp.vcf

$ bcftools view --threads 40 B3.RISG.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
$ python2 VCFtoPopPool.py temp.vcf RISG.Block3.sync 
$ rm temp.vcf

$ bcftools view --threads 40 B4.RISG.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
$ python2 VCFtoPopPool.py temp.vcf RISG.Block4.sync 
$ rm temp.vcf

$ bcftools view --threads 40 SNP.RISG.TRSdp.20.B90.2a.perp.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
$ python2 VCFtoPopPool.py temp.vcf RISG.All.Blocks.sync 
$ rm temp.vcf
```

#### Add coverage stats to sync file and filter by minimum coverage

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/Puritz_etal_CASE/refs/heads/main/scripts/add_cov_sync
$ chmod +x add_cov_sync
```

First, I added the coverage stats to the sync file, then filtered out sites with a minimum coverage less than 10 and a mean coverage less than 24

**Block 2**

```
$ mawk -f scripts/add_cov_sync RISG.Block2.sync | mawk '$14 > 10 && $16 > 24'> RISG.Block2.cov.sync

$ cat RISG.Block2.cov.sync | wc -l
702068
```

95 sites filtered out.

**Block 3**

```
$ mawk -f scripts/add_cov_sync RISG.Block3.sync | mawk '$14 > 4 && $16 > 24'> RISG.Block3.cov.sync

$ cat RISG.Block3.cov.sync | wc -l
144939
```

10 sites filtered out. I tried different minimum coverage values from 0-10 and the same number of loci were filtered out each time. 

**Block 4**

```
$ mawk -f scripts/add_cov_sync RISG.Block4.sync | mawk '$14 > 10 && $16 > 24'> RISG.Block4.cov.sync

$ cat RISG.Block4.cov.sync | wc -l
732115
```

403 sites filtered out.

**Full dataset**

```
$ mawk -f scripts/add_cov_sync RISG.All.Blocks.sync | mawk '$34 > 4 && $36 > 24'> RISG.All.Blocks.cov.sync

$ cat RISG.All.Blocks.cov.sync | wc -l
1011194
```

~42,081 sites filtered out

The remainder of the analysis will be conducted in R.