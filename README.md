# COSPEDTree2 (or COSPEDBTree)
Improved COSPEDTree (Couplet based supertree) with much lower running time, and option for binary supertree generation.

Description
-----------------

COSPEDTree-II / COSPEDTree2 (or COSPEDBTree) is an extended and much improved version of our earlier proposed 
couplet supertree algorithm COSPEDTree (https://github.com/souryacs/COSPEDTree/).

Input
-------

Phylogenetic trees possibly with overlapping taxa subsets, and having conflicting topologies among its taxa subsets. 
Input trees can be either in NEWICK format or in NEXUS format (but all the trees should follow the same format).

Output
------

A supertree covering all the input taxa. Output tree is generated in the NEWICK format.

Improvements compared to the earlier COSPEDTree algorithm
---------------------------------------------------------

1) Novel fractional and dynamic values of couplet based relations.

2) Resolving directed acyclic graph (DAG) to a supertree using the internode count measure (Liu et. al. 2011).

3) Optional feature to generate strict binary supertree. The technique uses a reference package named thSPR/thTBR 
(Lin et. al. 2009). Such package is also included in this archieve.

4) Much lower computation and running time.

Dependencies / Installation Requirements
--------------------------

COSPEDTree-II / COSPEDTree2 (or COSPEDBTree) is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. 
User needs to install following libraries before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: Current version does not support Python 3. A future release would cover it.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: New releases of Dendropy (like 4.1.0 or later) are not fully tested. Future release of this code would 
cover it.

3) Numpy ( available on the link: http://www.numpy.org/ )

Available via pip (python software downloader tool).


UBUNTU version issues
-------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

Execution
------------

python COSPEDBTree.py [options]

-h, --help

	show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME 

	name of the input text file containing source trees. 

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

	1 - input file format is NEWICK (default)
	2 - input file format is NEXUS

-b, --binary

	This is a boolean flag option. 
	Specifying this option toggles the default configuration.
	if TRUE, it produces a strictly binary supertree.
	Otherwise, the tree can be non-binary. Default FALSE.

-u, --underscore

	this is a boolean flag option, used for Dendropy based tree processing. 
	if TRUE, then this option preserves the underscores of the names of taxa. 
	Otherwise, the underscores in the taxon names are not preserved. Default TRUE.

*** NOTE: Last two options signify toggle / complement of their corresponding DEFAULT values. 

Example of a command (followed for the results published in the manuscript)

python COSPEDBTree -I source_tree_input.txt -p1

Depending on the value of -b option (TRUE / FALSE), following folder is created in the same directory containing 
the input treelist file:

	- if -b is FALSE, a folder "COSPEDTree2" is created.
	- if -b is TRUE, a folder "COSPEDTree2_Resolved" is created.

The output supertree is placed in a text file "cospedtree2_newick.tre" within the generated folder. The tree is 
stored in newick format.

Citation
--------

COSPEDTree-II is an implementation of the following paper:

1) COSPEDTree-II: Improved Couplet Based Phylogenetic Supertree, IEEE international conference on bioinformatics and biomedicine (BIBM), 2016, Shenzhen, China, pp. 98-101.

In addition, details of COSPEDTree algorithm is provided in the following articles:

1) Sourya Bhattacharyya, and Jayanta Mukherjee, "COSPEDTree: COuplet Supertree by Equivalence 
Partitioning of taxa set and DAG formation", IEEE/ACM Transactions on Computational Biology and Bioinformatics, 
Vol 12, No 3, pp. 590-603.

2) Sourya Bhattacharyya, and Jayanta Mukhopadhyay, "COuplet Supertree by Equivalence Partitioning of 
taxa set and DAG formation", Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology 
and Health Informatics (ACM-BCB), Newport, California, September 2014, pages 259-268.

For any queries, please contact
-------------------------------

Sourya Bhattacharyya 

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

email: <sourya.bhatta@gmail.com>

