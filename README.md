[logo]:https://imgur.com/a/Wl8eyDz

# Qmatey
A taxonomic profiler capable of robust microbiome analysis 
## Features
* Automated pipeline starting with QC-filtered sequencing data and ending with host-associated microbiome quantification
* Supports all Next-Generation Sequencing platforms
* Identifies strain-level microbial interactions 
* Capable of species-level, genus-level, and family-level taxonomic analysis 
* QC-plots to directly interact with quantified metagenomic data
## Installation 
Clone or download the git repository to a desired location 

```
$ git clone https://github.com/B-Kristy/Qmatey.git
```

## Dependencies
* R version compatible with the following dependencies: 
   * ggplot2 
   * plotly
   * plyr
   * dplyr
   * htmlwidgets 
   * car 
* Java 
  * Try 'sudo apt install default-jre'
* Datamash
  * Try 'sudo apt install datamash' 
* If you are performing a local BLAST, you will require an NCBI sequencing database compiled into one directory

## Setting Up a Project Directory 
A project directory should contain the following sub-directories:
* Input Sequences
  * This is where your QC-filtered sequencing data will go.
* Reference Genome
  * This is where your host associated reference genome(s) will go. Reference genomes must be in **FastQ format**.
* Configuration file
  * The format of the configuration file can be taken from the tools directory of the Qmatey Repository. 
## Preparing A Database Directory for a Local BLAST
If necessary, install the lftp tool to navigate NCBI's FTP site:
```
sudo apt-get install lftp
```
To obtain an NCBI sequencing database, go to the FTP site ftp://ftp.ncbi.nlm.nih.gov/blast/db/.
Create a database directory and select the sequencing database from the FTP site you wish to obtain using the following commands: 
```
mkdir database_directory
cd database_directory
wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.00.tar.gz
```
In the above code, I am extracting the nt.00.tar.gz nucleotide database file into my database directory. A complete sequencing database will require an extensive amount of space. 

Next, uncompress the database files and remove them with the following commands:
```
tar -xzvf nt00.tar.gz
rm nt00.tar.gz
```

Your database directory should now have the desired, uncompressed database files.

## Configuration

Variable | Usage | Input
-------------- | ------------------------------------------------------------------- | -----
input_dir      | the path to QC-filtered sequencing data                             | e.g. /path/to/dir/
ref_dir        | the path to the host reference genomes                              | e.g. /path/to/dir
db_dir         | the path to a local NCBI sequencing database                            | e.g. /path/to/dir/nt or NA
threads        | the maximum number of subprocesses that can run simultaneously      | integer 
tool_dir       | the path to QMatey's tools                                          | e.g. /path_to_github_repository/tools
strain_level   | An option for strain-level taxonomic analysis                       | TRUE or FALSE
species_level  | An option for species-level taxonomic analysis                      | TRUE or FALSE
genus_level    | An option for genus-level taxonomic analysis                        | TRUE or FALSE
family_level   | An option for family-level taxonomic analysis                       | TRUE or FALSE
blast_location | An option to perform BLAST locally or remotely                      | LOCAL or REMOTE
remote_db_dir  | the NCBI database for remote BLAST performance                    | e.g. nt, 16s, nr, etc. or NA

## Usage 
Before running QMatey, make sure you have:
* Created a project directory with all the required subdirectories
* Have QC-filtered sequencing data in the appropriate input directory
* Obtained a host-reference genome in .fastq format in the appropriate reference genome directory
* **for a local BLAST**: Obtained an NCBI database directory and have the uncompressed files in one directory
* **for a remote BLAST**: identified an NCBI database directory online 
* Have a correctly edited configuration file within the project directory 

From the command line, type: 
```
$ bash <path to github repository>/Qmatey_v0.1.sh <path to project directory>/Qmatey.config
```
# License 
<a href="https://github.com/tararickman/metagenome/blob/add-license-1/LICENSE"> GNU General Public License v3.0
