# COSPEDBTree
Improved COSPEDTree (Couplet based supertree) with much lower running time, and option for binary supertree generation

Description
-----------------
COSPEDBTree is an extended version of our proposed COSPEDTree algorithm, which produces supertrees from input phylogenetic trees. Input phylogenetic trees may contain overlapping taxa set. These trees often exhibit different topologies among constituent taxa set. The objective is to produce a supertree covering all the input taxa so that individual taxa subsets exhibit consensus relationships as much as possible. This is known as satisfying "Maximum agreement property"

However, given the conflicting nature of input trees, often the consensus (most freqent) relations among individual taxa subsets may not be reflected in the final supertree. This is because, consensus relation among a taxa subset may conflict with the consensus relation of another taxa subset.

In the basic COSPEDTree algorithm, supertree computation is formed using a greedy approach. The method is based on partitioning the input set of taxa based on equivalence relation. The relation is defined between individual taxa pairs (couplet), leading to the proposed couplet based supertree technique.

The relationship among a couplet can be of the following three types: 1) ancestor / descendent 2) sibling 3) inconclusive relation characteristics. The sibling relation is in fact, an equivalence relation. Definition of these relations can be found in our paper (reference provided below). According to the equivalence relation, input taxa set is partitioned. After that, a Directed Acyclic Graph (DAG) is formed initially. From this DAG, the output tree is generated.

Input source trees can be either in NEWICK format or in NEXUS format. However, all the source trees should have identical input formats. Output tree is generated in the NEWICK format.

COSPEDBTree proposes improvement of the earlier COSPEDTree algorithm, by the following ways:

1) It improves the running time performance (lowers the running time) by applying a suitable clustering approach during finalizing the couplet relationships in the final supertree. This saves running time considerably. Along with, input phylogenetic tree processing and derivation of the couplet relationships are performed much quicker, leading to the FASTEST supertree construction (as compared to the reference approaches) !!

2) The algorithm has an option to produce strict binary supertree from the input source trees. Earlier COSPEDTree algorithm produced unresolved non-binary supertrees. COSPEDBTree employs a binary refinement technique by proposing a mapping between the (intermediate) non-binary supertree to individual source trees, and resolves multifurcations using an agglomerative clustering approach.

3) The multiple parent problem in the basic COSPEDTree technique has been resolved with a deterministic depth first strategy. The priority measures of individual couplets have been used to select the parent node of a node, for using in the final supertree.

Dependencies / Installation Requirements
--------------------------

COSPEDBTree is developed in Linux Systems (Ubuntu 12.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 
2) Dendropy ( available on the link: https://pythonhosted.org/DendroPy/ ) 
3) Numpy ( available on the link: http://www.numpy.org/ )

For systems having Ubuntu with lower versions, please notify in case of any errors.

Execution
------------

COSPEDBTree is to be executed with the following command line options, from a terminal: (assuming the present working directory contains the source codes)

chmod +x COSPEDBTree.py (To change its permission to make it an executable file)

./COSPEDBTree.py [options]

NOTE:

All the options except the first three, signify toggle / complement of their corresponding DEFAULT values. First option (help) displays these command line parameters.

It Is Preferable For A Beginner, To Not Use Any Option Other Than The Second And Third Options. Second option is for specifying the input filename (mandatory) Third option is for specifying the corresponding file format.

Details of the options are mentioned below:

-h, --help
show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME 

                name of the input file containing candidate source trees (a text file) USER MUST PROVIDE ONE VALID                 INPUT FILE CONTAINING THE TREE DATA OTHERWISE PROGRAM WILL BREAK FROM EXECUTION

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME
                name of the output file which will contain the target supertree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                    1 - input file format is NEWICK (default)
                    2 - input file format is NEXUS       
                    USER MUST PROVIDE either p = 1 or 2

-q NO_OF_QUEUES, --queues=NO_OF_QUEUES

                    1 - only a single max priority queue is used for
                    storing the score metrics                        
                    2 - two separate queues are used to store the conflicting
                    and non conflicting taxa pairs and corresponding score metrics (default)

-d, --dfsref          

                    if TRUE, Multiple parent problem (C2) is tackled -
                    before appying DFS for arbitrary parenting
                    information, it selects the most probable parent candidate. Default = True

-b, --binary          

                    if TRUE, it produces a strictly binary supertree.
                    Otherwise, the tree can be non-binary. Default FALSE.


Example of a command (followed for the results published in the manuscript)

./COSPEDBTree -I source_tree_input.txt -p1 -b > out.txt

command descriptions: 

  1) -I specifies the input filename 
  2) source_tree_input.txt : contains the input collection of trees 
  3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, as specified by the option (-p1) (1 stands for newick)
  4) -b option enables the binary supertree production.
  
  Here the no of queues = 2 and deterministic depth first strategy for final supertree production is employed, as shown in the default configurations.



The output tree and all the results are printed at console. User can redirect the output results to any standard text file by using standard redirection operation (>). For example, in the above command, all the detailed results (textual descriptions) are redirected to file out.txt.

In addition, one output file "output_supertree_newick.tre" is created in the current directory it contains the derived supertree information (in both newick string format as well as tree plot). The tree can be used subsequently for performance metric computation 


Utilities
-----------

COSPEDTree requires O(N^3 + MN^2) time and O(N^2) space complexity, for N input taxa and M input trees. Both the time and space complexities are equal or lower than existing species tree construction methods.

Citation
--------

Upon using this package, users need to cite the following articles:

1) Sourya Bhattacharyya, and Jayanta Mukherjee, "COSPEDTree: COuplet Supertree by Equivalence Partitioning of taxa set and DAG formation", IEEE/ACM Transactions on Computational Biology and Bioinformatics, pp. 1-1, doi: http://doi.ieeecomputersociety.org/10.1109/TCBB.2014.2366778

2) Sourya Bhattacharyya, and Jayanta Mukhopadhyay, "COuplet Supertree by Equivalence Partitioning of taxa set and DAG formation", Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology and Health Informatics (ACM-BCB), Newport, California, September 2014, pages 259-268.

For any queries, please contact

Sourya Bhattacharyya 
Department of Computer Science and Engineering 
Indian Institute of Technology Kharagpur 
email: sourya.bhatta@gmail.com
