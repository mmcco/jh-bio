repeatgenome
============

While still under heavy development, this package is in what will likely be its long-term and stable form.

It implements a slightly modified version of Kraken, a metagenomics algorithm. For some constant *k*, Kraken associates all *k*-long subsequences (*"kmers"*) of a set of reference genomes with their least-common ancestor (LCA) in the taxonomy tree of the reference genomes used. For example, a certain kmer may exist only in a given genus. Kraken estimates which taxonomic class (e.g. a family, genus, or species) a sequencing read came from based on the classifications of its kmers. For a full description, see the [journal article](http://genomebiology.com/2014/15/3/R46).

RepeatGenome adapts this technique to estimate where in the genome of a known species a sequencing read came from. This is done using repeat sequences, which are prolific and of low complexity. Like taxonomic classifications, repeat sequence types are organized in a heirarchical tree. For example, the tree of repeat types may contain a subtree with the root node "Satellite", whose child nodes would be different types of [Satellite DNA](https://en.wikipedia.org/wiki/Satellite_DNA). These heirarchies are supplied by projects like [Repbase](http://www.girinst.org/repbase/index.html). This allows us to associate kmers with nodes in a repeat classification tree like Kraken does with a taxonomy tree.

We use [RepeatMasker](http://repeatmasker.org/) to determine what parts of a reference genome are repeat sequences. Therefore, RepeatGenome's input comprises a RepeatMasker output file and the associated reference genome in FASTA format.

This is an admittedly brief description of the program. We may publish a more thorough coverage such as a white paper eventually. However, if you want to fully understand the code, it would be best to read the Kraken article, the apparently relevant parts of the [RepeatMasker documentation](http://repeatmasker.org/webrepeatmaskerhelp.html), and the source code's function comments. The RepeatGenome source code is heavily commented, and the result of the [godoc](http://blog.golang.org/godoc-documenting-go-code) documentation generater will soon be posted publicly.

Integers are usually explicitly declared as having type int64 or uint64, but due to the size of the minimizer maps and of genomes in general, RepeatGenome probably won't run on 32-bit systems. You may be able to make it work for small *m* and *k* values and small reference genomes, but there are no guarantees.
