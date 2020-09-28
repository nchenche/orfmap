ORFMap
===========

ORFMap - A tool aimed at scanning a genome for stop-codons delimited sequences (ORFs) and annotating them.


| test | test |
| --- | --- |
| test | test |

Description
-----------

From a genomic fasta file and its associated GFF, the program first scans the genome to retrieve all sequences
delimited by stop codons. Only sequences of at least 60 nucleotides long are kept by default.

Those so-called ORF sequences are then annotated depending upon GFF element type(s) used as a reference. The CDS element type is always used as a reference but others can be added. By default an ORF sequence has 5 possible annotations whose 3 are listed below:
* c_CDS                if the ORF overlap with a CDS in the same phase
  - nc_5-CDS           if the 5' extremity of the c_CDS is at least 60 nucleotides long
  - nc_3-CDS           if the 3' extremity of the c_CDS is at least 60 nucleotides long	
* nc_ovp-CDS           if the ORF overlap with a CDS in a different phase
* nc_intergenic        if the ORF do not overlap with anything	
 
Note that if an ORF sequence is tagged as 'c_CDS', this sequence is further processed to be cut at its 5' and 3' extremities that do not overlap with the CDS. If their length is above or equal to 60 nucleotides, then those subsequences can be assigned as nc_5-CDS and/or nc_3-CDS.
 

The user can also specify what GFF element type(s) can be used as reference(s) to annotate ORF sequences in addition to the CDS type. For instance, if the user adds the tRNA element type, ORF sequences could now be assigned as nc_ovp-tRNA if they overlap with a tRNA. Thus 6 assignments would now be possible for an ORF sequence:
* c_CDS                if the ORF overlap with a CDS in the same phase
  - nc_5-CDS           if the 5' extremity of the c_CDS is at least 60 nucleotides long
  - nc_3-CDS           if the 3' extremity of the c_CDS is at least 60 nucleotides long
* nc_ovp-CDS           if the ORF overlap with a CDS in a different phase
* nc_ovp-tRNA          if the ORF overlap with a tRNA
* nc_intergenic        if the ORF do not overlap with anything


**Note on default parameters**:
* CDS in the only element type used as a reference to annotate ORF sequences.
* the minimum nucleotide number required to consider an ORF sequence is set at 60 nucleotides
* an ORF sequence is considered as overlapping with an element (e.g. CDS) if at least 70 % of its sequence overlap with the element or if this element is totally included within the ORF sequence


----------------------------------------
Installation procedure from distribution
----------------------------------------


I. First steps
------------

1. Uncompress and untar the package:

```bash
tar -xzvf orfmap-0.0.tgz
```

2. Go to the ORFMap directory
 
```bash
cd ORFMap-0.0
```


II. Install the package in a virtual environment (the recommended way to avoid dependencies conflict)
------------------------------------------------------------------

1. Install virtualenv (if necessary)

```bash
pip install virtualenv
```

2. Create a virtual environment (for python3)

```bash
virtualenv -p python3 env_orfmap
```

3. Activate your virtual environment

```bash
source env_orfmap/bin/activate
```

**Note that once activated:**:
* you should see the name of your virtual environment in brackets on your terminal line
* any python commands will now work within your virtual environment
	


4. Install ORFMap in your virtual environment

```bash
python setup.py install
```

**Note**: once installed, you should be able to run orfmap (see below). Once you don't need to use it, you can deactivate or exit your virtual environment by executing in the terminal:

```bash
deactivate
```

From this installation, everytime you'll want to use orfmap you'll need to activate your dedicated virtual environment.



II_bis. Install the package in the standard python libraries of your system (not recommended)
---------------------------------------------------------------------------------------------

In ORFMap-0.0/:

```bash
python setup.py install
```


----------------------------------------
Running ORFMap
----------------------------------------

Basic run
---------

ORFMap requires two input files: a genomic fasta file (-fna) and its associated GFF file (-gff). The most basic run can be executed by typing:

```
run_orfmap -fna mygenome.fna -gff mygenome.gff
```

All of the ORF sequences are annotated relative to the CDS element type only. Thus 5 possible annotations are possible:
* c_CDS                if the ORF overlap with a CDS in the same phase
  - nc_5-CDS           if the 5' extremity of the c_CDS is at least 60 nucleotides long
  - nc_3-CDS           if the 3' extremity of the c_CDS is at least 60 nucleotides long	
* nc_ovp-CDS           if the ORF overlap with a CDS in a different phase
* nc_intergenic        if the ORF do not overlap with anything	


The output will be two separated files with the prefix "mapping_orf_":
* mapping_orf_mygenome.fa: 	a proteic fasta file of all the ORFs sequences found
* mapping_orf_mygenome.gff:	A GFF file describing all the ORFs sequences found
  
By default, the two output files will contain all possible 5 annotations mentionned above.


Usage description
-----------------

To see all options available:

```
run_orfmap -h
```

This command will show:

usage: run_orfmap [-h] -fna [FNA] -gff [GFF] [-type TYPE [TYPE ...]] [-o_include O_INCLUDE [O_INCLUDE ...]] [-o_exclude O_EXCLUDE [O_EXCLUDE ...]] [-orf_len [ORF_LEN]]
                  [-co_ovp [CO_OVP]] [-out [OUT]]

Genomic mapping of pseudo-ORF

optional arguments:
- -h, --help            show this help message and exit
- -fna [FNA]            Genomic fasta file (.fna)
- -gff [GFF]            GFF annotation file (.gff)
- -type TYPE [TYPE ...] Type feature(s) a flag is desired for ('CDS' in included by default).
- -o_include O_INCLUDE [O_INCLUDE ...]  Type feature(s) and/or Status attribute(s) desired to be written in the output (all by default).
- -o_exclude O_EXCLUDE [O_EXCLUDE ...]  Type feature(s) and/or Status attribute(s) desired to be excluded (None by default).
- -orf_len [ORF_LEN]    Minimum number of nucleotides required to define a sequence between two consecutive stop codons as an ORF sequence (60 nucleotides by default).
-co_ovp [CO_OVP]      Cutoff defining the minimum CDS overlapping ORF fraction required to label on ORF as a CDS. By default, an ORF sequence will be tagged as a CDS if at least 70 per cent of its sequence overlap with the CDS sequence.
- -out [OUT]            Output directory


Except -fna and -gff arguments that are mandatory, all others are optional.

By default:
- -type 	CDS
- -o_include 	'all'
- -o_exclude	None
- -orf_len	60
- -co_ovp	0.7
- -out		'./'


Usage examples
-----------------

This command will use tRNA and snRNA element types as a reference to annotate ORF sequences and the outputs will be created in myResults/
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -out myResults
```

7 annotations will be possible for an ORF sequence:
 - c_CDS 		if the ORF overlap with a CDS in the same phase
   - nc_5-CDS		if the 5' extremity of the c_CDS is at least 60 nucleotides long
   - nc_3-CDS		if the 3' extremity of the c_CDS is at least 60 nucleotides long
 - nc_ovp-CDS		if the ORF overlap with a CDS in a different phase
 - nc_ovp-tRNA		if the ORF overlap with a tRNA
 - nc_ovp-snRNA	if the ORF overlap with a snRNA
 - nc_intergenic	if the ORF do not overlap with anything
 
To make the output files having only ORF sequences mapped as nc_ovp-tRNA and nc_ovp-snRNA:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -o_include nc_ovp-tRNA nc_ovp-snRNA -out myResults
```

To make the output files having all ORF sequences except those mapped as c_CDS:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -o_exclude c_CDS -out myResults
```

or:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -o_exclude coding -out myResults
```

**Note**: -o_include and -o_exclude take either feature types or a status attribute as arguments. Feature types have to be amongst the possible annotations for ORF sequences (e.g. c_CDS, nc_5-CDS, nc_intergenic...) while status attribute is either 'coding' or 'non-coding' ('coding' refers to c_CDS and 'non-coding' refers to the other ones).



This command will define ORF sequences if they are at least 50 nucleotides
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -orf_len 50
```


This command will consider an ORF sequence as overlapping with an element (e.g. CDS) if at least 60 % of its sequence overlap with the element or if this element is totally included within the ORF sequence
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -co_ovp 0.6
```





