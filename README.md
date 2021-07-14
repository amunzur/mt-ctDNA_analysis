# mt-ctDNA_analysis

## Scripts: 
* `convert_bamtosam.sh`: Takes the bam files with mito chromosomes only and converts them to sam files and saves to a different location. 
The file also does some filtering, such as removing reads with MQ less than 20 and length less than 40. 

* `extract_MT_reads.sh`: From the WGS bam files, extract the mitochondrial chromosomes and save the bam files to a different location. 

* `fragment_size_plot.R`: Reads the sam files as a text file and makes a plot of the fragment sizes using R. In progress. 

* `get_read_numbers.sh`: Compare the read numbers as we filter the unaltered WGS bam files. 

* `make_bamList.R`: This file makes a list of the bam files to generate tnvstats and to run Matti's variant calling pipeline. 
