# **Composite Transposon Finder (TnComp_finder)**

Composite Transposon Finder (TnComp_finder) is a program for the prediction of putative composite transposons in bacterial and archaeal genomes based on insertion sequence replicas in a relatively short span. It works by comparing nucleotide sequences from bacterial and archaeal genomes to a custom transposon database. Results are writen in report files and pre-anotated GenBank files to help in subsequent manual curation.

---

## RELEASE

Version 1.0.0 - November 11, 2019.

Available from <https://github.com/danillo-alvarenga/tncomp_finder>.

---

## REQUIREMENTS

TnComp_finder runs on Python 3.4 or newer and has been tested on Debian 9, Ubuntu 16.04 and CentOS 6.7. It should work on any modern GNU/Linux distro. TnComp_finder depends on the Biopython library version 1.66 (for Python 3), the BLAST+ suite version 2.2.28 and Prodigal version 2.6.1. Newer versions of these packages are likely to work, but if you encounter problems you may try downgrading to the approved versions.

---

### DATABASES

TnComp_finder depends on a transposon database. As an alternative to the database included with this program, it is possible to use custom databases by retrieving complete nucleotide sequences for transposable elements and writing them into a multifasta file called `transposons.fna`. Headers should follow the pattern `$ID|$FAMILY` (where `$ID` is the database ID number for the sequence and `$FAMILY` is the family the transposon is classified in). This file should be moved to a subdirectory called `db` in the directory where the `TnComp_finder.py` script is stored. We recomend getting sequences from ISfinder (<https://www-is.biotoul.fr/>).

---

## INSTALLATION

For TnComp_finder to work the blastn and prodigal executables need to be installed and available in your path. Follow installation instructions from distributors of these packages for your particular distro. In Debian and derivative distros, these programs are usually available in the official repositories and can be installed by issuing the following commands in a terminal in superuser mode:  
`apt install ncbi-blast+`  
`apt install prodigal`  
`apt install python3-biopython`  

Biopython can also be installed as superuser by pip3:  
`pip3 install biopython`  

For convenience, after downloading and extracting the latest TnComp_finder release from GitHub you can move the TnComp_finder directory to the desired destination and add it to your path:  
`mv TnComp_finder-1.0.0/ bioinformatics/`  
`cd bioinformatics/TnComp_finder-1.0.0/`  
`echo 'export PATH=$PATH:'$(pwd) >> ~/.bashrc`  
`source ~/.bashrc`  

>**Note:** CentOS 6 only ships with the legacy Python 2 version. In case you are using a Red Hat-based distribution in which Python 3 is still unavailable, the latest Python version can be installed from source by issuing the following commands to the terminal as superuser:  
>1. Install the OpenSSL development package: `yum install openssl-devel`  
>2. Download the latest stable Python version (Python 3.7.0 at this moment): `wget https://www.python.org/ftp/python/3.7.0/Python-3.7.0.tgz`  
>3. Extract the package: `tar -xzf Python-3.7.0.tgz`  
>4. Change to the decompressed directory: `cd Python-3.7.0/`  
>5. Configure the build scripts: `./configure`  
>6. Build the executables: `make`  
>7. Install the built binaries: `make install`  

>**Note:** Older or newer versions of these programs as well as other operating systems might also work, but since they have not been tested this is not guaranteed.
>**Note:** The main scripts looks for the `db` directory in the same directory as the main script is stored. Moving it will cause the script to fail.

---

## USAGE

`TnComp_finder.py [-h] [-v] -f sequences.fasta [sequences.fasta ...] [-o directory] [-p threads] [-g] [-i %] [-c %] [-d bp] [-e bp] [-s | -t]`  

optional arguments:

`-h`, `--help` | show a help message and exit  
`-v`, `--version` | show version and exit  

`-f sequences.fasta [sequences.fasta ...]`, `--file sequences.fasta [sequences.fasta ...]` | target sequences  
`-o directory`, `--out directory` | output directory  
`-g`, `--gbk` | write a genbank file with predictions  

`-p threads`, `--processors threads` | processor threads available for analyses  

`-i percentage`, `--identity percentage` | minimum identity with the database  
`-c percentage`, `--coverage percentage` | minimum coverage with the database  
`-d base pairs`, `--distance base pairs` | maximum distance between transposons  

`-e base pairs`, `--extend base pairs` | retrieve sequence with extended borders  

`-k `, `--flanking` | ignore candidates with similar flanking regions  
`-s`, `--same` | consider only copies on the same strand  
`-t`, `--different` | consider only copies on different strands  

In order to run TnComp_finder, indicate the target genome file in the nucleotide fasta format with the `-f`|`--file` argument and add flags adjusting the default parameters to your preferences. Optionally, indicate an output directory for the analyses results with the `-o`|`--out` argument. TnComp_finder will output an annotated genbank file if you add the `-g`|`--genbank` flag to the command line.

Complete genomes can only be analyzed in a single processor core, but multiple files and diferent contigs from draft genomes in multifasta files can be analyzed concurrently in the number of cores you determine by adding the `-p`|`--processors` parameter followed by the number of threads to be used.

The `-i`|`--identity` parameter determines the minimum percentage identical nucleotides between query and database sequences and defaults to 90 %. Likewise, the `-c`|`--coverage` parameter indicates the minimum allowed percentage of alignment size compared to the database sequence size, defaulting to 95 %. Finally, `-d`|`--distance` establishes the farthest the detected features can be located in the query sequence and defaults to 20000 bp. Borders for the candidates found can be extended by adding a number of base pairs with `-e`|`--extend`. Estimation of replicas may be restricted to either those in the same strand with `-s`|`--same` or in different strands with `-t`|`--different`.

The `-k`|`--flanking` flag should be used in combination with the `-e`|`--extend` parameter. This mode allows to compare the regions flanking the putative composite transposon to each other and exclude it from the list of candidates when these regions are too similar. More specifically, it disregards candidates presenting alignments in at least half the length of the extended regions that are above the provided identity and coverage thresholds.

**Examples**:  
`TnComp_finder -f Escherichia\_coli\_A123\_complete\_genome.fasta`  
`TnComp_finder -f Escherichia\_coli\_B456\_draft\_genome.fasta Escherichia\_coli\_C789\_draft\_genome.fasta -o results -gk -p 4 -i 80 -c 70 -d 10000`  

---

## RESULTS

TnComp_finder generates a report file indicating the query sequence followed by potential candidates, detailing A and B transposon length, distance between replicas, their general characteristics, sequences and predicted ORF coordinates in the contemplated genome region. If you add the `-g`|`--genbank` flag, the program will also output a .gbk file with annotations corresponding to the analysis results which can be manually curated in programs such as the Artemis genome browser.

The report and gbk files include sequences for each of the predicted transposon copies and for the whole sequence, from first to last ORF. If the `-e`|`--extend` parameter was included, an additional sequence that is longer upstream and downstream than the whole sequence is provided.

Information about run is stored in the file `info.txt`, which is also written in the working directory.

---

## CITATION

If you find this program useful in your research, please cite the following references:
- Alvarenga DO, Varani AM (2019). TnComp_finder: prokaryotic composite transposon finder. Available from <https://github.com/danillo-alvarenga/tncomp_finder>.
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421
- Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Freidberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423
- Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11:119
- Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M (2006) ISfinder: the reference centre for bacterial insertion sequences. Nucleic Acids Res 34:D32-D36

---

## ISSUES AND REQUESTS

If you experience any issues or would like to see support for an additional feature implemented in TnComp_finder, please file a request in the GitHub issue tracker or email it to the developer. Please feel free to contact the developer with any further questions regarding this software. You can reach him at <mailto:danillo.alvarenga@gmail.com>.

---

## LICENSE

Copyright © 2019 Danillo Oliveira Alvarenga

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/agpl-3.0.html>.
