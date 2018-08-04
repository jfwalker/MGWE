# Maximum Gene Wise Edge (MGWE)

This Github contains the MGWE program, a two topology program and the data used in the analysis.

`August 4th 2018: Further edge work has been moved to a new github https://github.com/jfwalker/PHAIL 
This is all left up since it was used in the paper but the constraint method mentioned in the paper
is no implemented in the program EdgeTest.py from PHAIL. Any questions please send an email to
jfwalker@umich.edu`




### Procedure

Details regarding the MGWE procedure may be found in the paper:

SysBio: [Analyzing contentious relationships and outlier genes in phylogenomics](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy043/5034973)

BioRxiv: [Analyzing contentious relationships and outlier genes in phylogenomics](https://www.biorxiv.org/content/early/2018/02/19/115774)


### Requirements

The phylogenetics suite: [phyx](https://github.com/FePhyFoFum/phyx/)

The ML program: [raxml](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

### Running the program

Once you've set up the config file (instructions below) you can run the program by using the command 

`EX: perl MGWE_Calc.pl Config.config`

### The Config file

Running the program should be relatively trivial and can be tested quickly on the vertebrate dataset from the paper. 
The configuration file first needs to be set up with the following options. Each option should have a space after the
colon. An example config file here is called Config.config and it's the one I used to analyze the vertebrate dataset

location of the phyx programs (pxbp and pxrmt), the location of raxml, and a name for the outfile.

You will also need to give the relationship you are curious about, to do this put all species contained within that clade as a comma
separated list. So when asking if turtles are sister to crocodiles in the paper the relationship we used was: alligator,caiman,phrynops,caretta,chelonoidis_nigra,emys_orbicularis.

The number of topologies you are testing can also be specified. This will cause the program to test multiple trees in your treeset by running an SSLL test and a test
where the topology is fixed and the likelihoods are calculated, similar to -M in raxml. If you specify 2 it will perform these analyses on the first 2 trees in your treeset
if you specify three it will perform them on the first three etc...

Next you will want to give the concatenated supermatrix file and then the partition file. Partitions should be set up as they are in the examples

You also need to give a treeset, this can be generated in a number of ways (e.g gene trees) but these trees should have complete taxon sampling.

Threads should then be specified and if you want to run a quick test of two genes that can be done by putting True after test

The program also creates an output folder where files are moved to so they can be double checked.

Their is a secret mode that does not recalculate the the SSLL for the concatenated matrices but this assumes you have already run the program once and are re-running with the same
settings.

### The Output

The output file will give a decent amount of data and each section is divided by a number of # in order to identify the different sections. The first section
prints out all the settings, please check over this because if anything is wrong that's probably a problem. The next section just mentions a few of the output
files.

The section called "your conflicts" has all the conflicts the program identified

The section called "your trees and their conflicts" has each tree in your treeset and which conflict they support. Due to nested conflict from a bipartition analysis
see [Smith et. al 2015](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-015-0423-0) you can get a single tree with multiple conflicts.

The next section will describe the missing data in the genes you use, these are removed from the trees when likelihood calculations are made. This does influence parameters
and is printed to make sure everything is correct.

The next section is the gene counts and this will show how many genes support a conflicting relationship.

Site counts is the number of sites that make up those genes

The parameter info is the number of parameters used in the analysis. Since this program uses GTR+G in raxml each gene tree has (9 + (2n-3)) parameters where the
n is the number of taxa. The 9 comes from 5 parameters for GTR, 1 for GAMMA, 3 for base frequencies. The 2n-3 is the number of branches (n is number of taxa).
In the case of having free branch lengths each gene has (9 + (2n-3)) parameters, if it is a concatenated analysis then the whole supermatrix only has one
set of branch lengths.

The raw likelihoods are listed here, although it is important to not that the concatenated analysis likelihoods (the first ones) are not comparable to those where the branch
lengths are allowed to vary called (Indv BrLen) or the edges.

The AIC scores correcting for the greater number of parameters used during the likelihood calc is listed next.

Finally the Delta AIC is given to help find which relationship proved to be the best.




