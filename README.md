# QMatey
A taxonomic profiler capable of microbiome analysis designed by biologists.
## Features
* Automated pipeline starting with raw sequencing data and ending with interactive visualizations 
* Supports all Next-Generation Sequencing platforms
* Identifies strain-level microbial interactions 
* Capable of both species-level and genus-level clustering 
* QC-plots to directly interact with metagenomic data 
* Easy to learn
## Installation 
Clone or download the git repository to a desired location 

```
$ git clone https://tararickman/QMatey.git
```

## Dependencies
* R version compatible with the following dependencies: ggplot2 plotly
* Java 
  * Try 'sudo apt install default-jre'
* Datamash
  * Try 'sudo apt install datamash' 
* Any NCBI sequence database 
 
## Setting Up a Project Directory 
A project directory should contain the following sub-directories:
* Input Sequences
  * This is where your QC-filtered sequencing data will go.
* Reference Genome
  * This is where your reference genome(s) will go. Reference genomes must be in **FastQ format**.
* Configuration file
  * The format of the configuration file can be taken from the tools directory of the QMatey Repository. 
## Preparing A Database Directory 
IF necessary, install the lftp tool to navigate NCBI's FTP site:
```
sudo apt-get install lftp
```
To obtain an NCBI sequencing database, go to the FTP site ftp://ftp.ncbi.nlm.nih.gov/blast/db/.
Create a database directory and select the sequencing database from the FTP site you wish to obtain using the following commands: 
```
cd $database_directory
wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.tar.gz
```
In the above code, I am extracting all compressed files of the **nucleotide (nt) database** into my database directory. A complete sequencing database will require an extensive amount of space. 

Next, uncompress the database files amd remove it with the following commands:
```
tar -xzvf *nt.tar.gz
rm *nt.tar.gz
```

Your database directory should now have the desired, uncompressed database files.

## Configuration

Variable | Usage | Input
-------------- | ------------------------------------------------------------------- | -----
input_dir      | the path to QC-filtered sequencing data                             | e.g. /path/to/dir/
ref_dir        | the path to the host reference genomes                              | e.g. /path/to/dir
db_dir         | the path to the NCBI sequencing database                            | e.g. /path/to/dir/nt
threads        | the maximum number of subprocesses that can run simultaneously      | integer 
tool_dir       | the path to QMatey's tools                                          | e.g. /path_to_github_repository/tools

## Usage 
Before running QMatey, make sure you have:
* Created a project directory with all the required subdirectories
* Have QC-filtered sequencing data in the appropriate input directory
* Obtained a host-reference genome in .fastq format in the appropriate reference genome directory
* Obtained an NCBI database directory and have the uncompressed files in one directory
* Have a correctly edited configuration file in the project directory 

From the command line, type: 
```
$ bash <path to github repository>/QMatey.sh <path to project directory>/QMatey.config
```
# License 
<a href="https://github.com/tararickman/metagenome/blob/add-license-1/LICENSE"> GNU General Public License v3.0
