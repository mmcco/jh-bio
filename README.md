repeatgenome
============

Introduction
------------

While still under heavy development, this package is in what will likely be its long-term and stable form.

It implements a slightly modified version of Kraken, a metagenomics algorithm. For some constant *k*, Kraken associates each *k*-long subsequence (*"k-mer"*) of a set of reference genomes with its least-common ancestor (LCA) in the taxonomy tree of the reference genomes used. For example, a certain k-mer may exist only in a given genus. The resulting library is then used to estimate which subtree of the taxonomy tree each of a set of sequencing reads came from. This is done by querying the reads' k-mers against the library. For a full description, see the [journal article](http://genomebiology.com/2014/15/3/R46).

RepeatGenome adapts this technique to estimate where in a single chosen species' genome a sequencing read originated. This is done using [repeat sequences](http://biol.lf1.cuni.cz/ucebnice/en/repetitive_dna.htm), which are prolific and of low complexity. Like taxonomic classifications, repeat sequence types are organized in a heirarchical tree. For example, the tree of repeat types may contain a subtree with the root node "Satellite", whose child nodes would be different types of [Satellite DNA](https://en.wikipedia.org/wiki/Satellite_DNA). These heirarchies are supplied by projects like [Repbase](http://www.girinst.org/repbase/index.html). This allows us to associate k-mers with subtrees in a repeat classification tree like Kraken does with a taxonomy tree.

We use [RepeatMasker](http://repeatmasker.org/) to determine what parts of a reference genome are repeat sequences. Therefore, RepeatGenome's input comprises a RepeatMasker output file and the associated reference genome in FASTA format.

This is an admittedly brief description of the program logic. We may publish a more thorough coverage such as a white paper eventually. However, if you want to fully understand the code, it would be best to read the Kraken article, the apparently relevant parts of the [RepeatMasker documentation](http://repeatmasker.org/webrepeatmaskerhelp.html), and the source code's function comments. The RepeatGenome source code is heavily commented, and the result of the [godoc](http://blog.golang.org/godoc-documenting-go-code) documentation generater will soon be posted publicly.

Integers are usually explicitly declared as having type int64 or uint64, but due to the size of the minimizer maps and of genomes in general, RepeatGenome probably won't run on 32-bit systems. You may be able to make it work for small *m* and *k* values and small reference genomes, but there are no guarantees.

Input
------------

The requirements for file input are currently somewhat strict. If you find them limiting, let me know and I can quickly add more flexible options.

The main argument to `repeatgenome` is the genome name. This should be the abbreviation of the reference genome being used - for example, most of our tests were run on the most recent *Drosophila melanogaster* (fruitfly) genome, dm3.

In the below file and directory names, `<genome-name>` is used to denote the above-mentioned genome name.

There is expected to be a subdirectory of the current working directory bearing the genome name. This directory should contain:
* all reference sequence files, in fasta format, with the `.fasta` suffix.
* a RepeatMasker output file with the `.fa.out` suffix.

Additionally, if any reads are to be processed, a subdirectory `<genome-name>-reads` of the current working directory is expected to contain all reads in SAM format, with the `.sam` suffix.

Output
------------

`<genome-name>-stats`: Contains any statistics data written in the course of processing.
