# MGWE

This Github contains the MGWE program, a two topology program and the data used in the analysis.

###Procedure

Details regarding the MGWE procedure may be found in the paper:

[Analyzing contentious relationships and outlier genes in phylogenomics](https://www.biorxiv.org/content/early/2018/02/19/115774)


###requirements:

The phylogenetics suite: [phyx](https://github.com/FePhyFoFum/phyx/)

The ML program: [raxml](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

###Running the program

Running the program should be relatively trivial and can be tested quickly on the vertebrate dataset from the paper. 
The configuration file first needs to be set up with the following options.

location of the phyx programs (pxbp and pxrmt), the location of raxml, and a name for the outfile.

You will also need to give the relationship you are curious about, to do this put all species contained within that clade as a comma
separated list. So when asking if turtles are sister to crocodiles in the paper the relationship we used was: alligator,caiman,phrynops,caretta,chelonoidis_nigra,emys_orbicularis.
