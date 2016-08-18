# COSPEDTree2 (or COSPEDBTree)
Improved COSPEDTree (Couplet based supertree) with much lower running time, and option for binary supertree generation.

Description
-----------------
COSPEDTree-II / COSPEDTree2 (or COSPEDBTree) is an extended version of our earlier proposed 
COSPEDTree algorithm [1,2], which produces supertrees from input phylogenetic trees. 
Input trees may contain overlapping taxa set. These trees often exhibit conflicting topologies among its taxa subsets. 

Input source trees can be either in NEWICK format or in NEXUS format. However, all the source trees should have identical input formats. 
Output tree is generated in the NEWICK format.

COSPEDBTree proposes the following improvements with respect to the earlier COSPEDTree algorithm:
(Details are mentioned in the corresponding paper):

1) Fractional and dynamic freqency measures between individual couplets are proposed, according to the 
coverage of taxa of the corresponding input trees.

2) The MPP problem [1,2] is dealt with a deterministic parent selection strategy. The internode count measure among 
individual couplets is used for such purpose.

3) Option to produce a strict binary supertree is also provided. The method uses a reference package named thSPR/thTBR (Lin et. al. 2009), 
and is also provided with the current archieve.

4) The running time of COSPEDTree-II is considerably lower than COSPEDTree, by resolving the pair of taxa clusters, 
rather than individual couplets.

Dependencies / Installation Requirements
--------------------------

COSPEDTree-II / COSPEDTree2 (or COSPEDBTree) is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following libraries before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: We have not tested the code on Python 3. Any user having Python 3 environment needs to check the 
correct execution of our code, and optionally needs to upgrade it accordingly.

We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.1.0 but we have used 3.12.0 for the implementation. 
We did not upgrade the code for supporting the latest version of DendroPy, so the users 
need to use the old version of the library to execute this code.

Support for the latest version of Dendropy (and corresponding update of the code) 
will be provided in some future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, 
which contains the latest Numpy module in it. We found that Numpy module in the 
traditional Apt-Get repository is of lower version.

UBUNTU version issues
-------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

Execution
------------

Assuming the current working directory contains the python source code of COSPEDBTree, 
following commands are to be provided from a terminal, to execute the code:

chmod +x COSPEDBTree.py (To change its permission to make it an executable file)

./COSPEDBTree.py [options]

NOTE:
-----

Last two options signify toggle / complement of their corresponding DEFAULT values. 
First option (help) displays these command line parameters. 
Details of the options are mentioned below:

-h, --help

	show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME 

	name of the input text file containing source trees. 
	ONE VALID INPUT FILE NEEDS TO BE PROVIDED.

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

	1 - input file format is NEWICK (default)
	2 - input file format is NEXUS       
	USER MUST PROVIDE either p = 1 or 2

-b, --binary

	This is a boolean flag option. 
	Specifying this option toggles the default configuration.
	if TRUE, it produces a strictly binary supertree.
	Otherwise, the tree can be non-binary. Default FALSE.

-u, --underscore

	this is a boolean flag option, used for Dendropy based tree processing. 
	if TRUE, then this option preserves the underscores of the names of taxa. 
	Otherwise, the underscores in the taxon names are not preserved. Default TRUE.

Example of a command (followed for the results published in the manuscript)

./COSPEDBTree -I source_tree_input.txt -p1

command descriptions: 

  1) -I specifies the input filename.

  2) source_tree_input.txt : contains the input trees.

  3) -p option specifies the input tree format (p = 1 means NEWICK format).
 
Depending on the value of -b option (TRUE / FALSE), following folder is created in the same directory containing 
the input treelist file:

	- if -b is FALSE, a folder "COSPEDTree2" is created.
	- if -b is TRUE, a folder "COSPEDTree2_Resolved" is created.

The output supertree is placed in a text file "cospedtree2_newick.tre" within the generated folder. The tree is 
stored in newick format.

Citation
--------

Upon using this package, users need to cite the following articles:

1) Sourya Bhattacharyya, and Jayanta Mukherjee, "COSPEDTree: COuplet Supertree by Equivalence 
Partitioning of taxa set and DAG formation", IEEE/ACM Transactions on Computational Biology and Bioinformatics, 
Vol 12, No 3, pp. 590-603.

2) Sourya Bhattacharyya, and Jayanta Mukhopadhyay, "COuplet Supertree by Equivalence Partitioning of 
taxa set and DAG formation", Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology 
and Health Informatics (ACM-BCB), Newport, California, September 2014, pages 259-268.

For any queries, please contact
-------------------------------

Sourya Bhattacharyya 
Department of Computer Science and Engineering 
Indian Institute of Technology Kharagpur 
email: sourya.bhatta@gmail.com

