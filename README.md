# MGWE

This Github contains the MGWE program, a two topology program and the data used in the analysis.

### Procedure

Details regarding the MGWE procedure may be found in the paper:

[Analyzing contentious relationships and outlier genes in phylogenomics](https://www.biorxiv.org/content/early/2018/02/19/115774)


### Requirements:

The phylogenetics suite: [phyx](https://github.com/FePhyFoFum/phyx/)

The ML program: [raxml](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

### Running the program

Running the program should be relatively trivial and can be tested quickly on the vertebrate dataset from the paper. 
The configuration file first needs to be set up with the following options. Each option should have a space after the
colon.

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


