# LncRNA-Snakemake
Full LncRNA data analysis workflow from sample.fastq to IsoMat.txt

**NOTE: Internal use only, do not publish this code without removing non-original files from /lib**

## Dependencies
* [Python 2.7](https://www.python.org/downloads/release/python-2715/) & [Python 3.8](https://www.python.org/downloads/release/python-380/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [STAR](https://github.com/alexdobin/STAR/releases)
* [RSEM](https://github.com/deweylab/RSEM)
* [stringtie](https://github.com/gpertea/stringtie)
* [slncky](https://github.com/slncky/slncky)
* [GTFSharp](https://github.com/acesnik/LncRNADiscovery)

## Using this workflow
**1.** Download this repository

**2.** Place your data in the project folder
   
&emsp;&emsp;Each replicate should be placed in its own folder, and the two sample names should have the format \[replicate name\]\_R\[1\\2\].fastq
   
**3.** Download genome files

&emsp;&emsp;Three files are required: \[genome\].fa, \[genome\].\[ensemble\].gtf, and \[genome\].\[ensemble\].gff3.

&emsp;&emsp;Place these files within the \\lib\\Genome folder.

**4.** Download Slncky annotation files

&emsp;&emsp;Place these within \\lib\\annotations.
   
**5.** Rename parameters within the snakemake file

&emsp;&emsp;`GENOME_VERSION` - set to your genome version (what you named the files in step 3).
   
&emsp;&emsp;`ENSEMBL_VERSION` - set to your ensemble version (what you named the files in step 3).
   
&emsp;&emsp;`SLNCKY_ANNOTATIONS` - set to your genome version (hg38 for human, mm10 for mouse)

**6.** Run `snakemake --use-conda` from within the project folder in a command line

**7.** Check replicate folders and output folder for results

&emsp;&emsp;Any steps pertaining to an individual replicate will be in that replicate's folder. All combined steps will be within \\output.

&emsp;&emsp;\\output\\logs will give run information for each step and also any errors that occured.

## Acknowledgements

  * Some code within the snakemake file was adapted from previous work by [Anthony Cesnik](https://github.com/acesnik/LncRNADiscovery), credit is also denoted within the file.
  * `GTFSharp`: Cesnik, A. J.; Yang, B.; Truong, A.; Spiniello, M.; Steinbrink, M.; Shortreed, M. R.; Frey, B. L.; Jarrard, D. F.; Smith, L. M. “Long Noncoding RNAs AC009014.3 and Newly Discovered XPLAID Differentiate Aggressive and Indolent Prostate Cancers.” Translational Oncology, 2018, 11, 808–814.
  * `RSEM`: Li B and Dewey CN (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12 (1), 323.
  * `STAR`: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, and Gingeras TR (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1), 15–21.
  * `stringtie`: Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT & Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122.
  * `slncky`: Chen J, Shishkin AA, Zhu X, Kadri S, Maza I, Guttman M, Hanna JH, Regev A, and Garber M (2016). Evolutionary analysis across mammals reveals distinct classes of long non-coding RNAs. Genome Biol 17(1), 19.

