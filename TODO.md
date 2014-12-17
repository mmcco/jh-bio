*This is more of a brain-dump for me so that I don't forget potential
additions, simplifications, or issues. Some of it is likely
outdated, and most of it is probably of little use to someone who isn't
very familiar with the code.

Bold items are algorithmic, the rest are minor*

* Premature commenting is the root of all evil, and I have sinned. Please read
comments skeptically - they are being audited.

* Should use iota values for different library and classification types.

* Should make a library for parsing and checking consistency of library
  format.

* Should do memory checking with a thread in its own function, or in a
  thread spawned directly from run-rg.

* Repeat.Instances is populated in an incorrect manner - it seems like it gets
one per ClassNode. (fixed, I believe - should recheck)

* Are we storing everything we need to in the Kraken data file? Could some of the
data that parsed libraries rely on (RepeatGenome.Matches, RepeatGenome.Repeats,
etc.) be non-deterministic?

* Functions for pointer accesses? (Done, with a few odd exceptions)

* Reconcile numRawKmers() and krakenFirstPass(), determine when to filter 'n's.

* Ensure that Match.RepeatStart and Match.RepeatEnd are zero-indexed.

* A lot of explicit variable types could be safely removed.

* Consider renaming Kmer to KmerPair or KrakenPair.

* Could populate RepeatGenome.ClassTree.ClassNodes in a separate loop.

* Non-exported globals or constants for bit offsets (e.g. that of a Kmer's LCA?)

* A ClassifyReadsFile() function that dispatches to different functions based on
the reads filetype

* We really need a map to get min's indexes/offsets in O(1) time. getMin()
currently does a binary search.

* KmerInt.Minimize() logic could be changed now that minimizers are 32 bits

* Should a Seq's first field be a *byte to discard the extra two fields? If not,
we could probably use len() in Seq manipulations.

* Should probably make a file solely for type defs.

* Reads are currently kept in TextSeq form until the bitter end because, with
Go's referenced based slices, there's no compelling reason not to, and because
they're easier (and probably faster) to manipulate than Seqs. This may change
at some point, though.

* If a minimizer is associated with a single repeat type, can we use that
heuristically?

* Error handling should be updated with a custom ParseError type - panics should
be removed, excepting performance-cricial sequence manipulation functions

* Should consider splitting at hyphenated class names like TcMar-Tc1

* The concurrent read-kmer generator could be reintroduced using a select
statement.

* Should probably restrict activity of chans with directionals

* It would make sense to discard kmers associated with ClassNodes greater than a
certain size.

* Kmer counting should be re-added eventually - it's currently excluded for
performance reasons because we aren't using it.

* We should test a version that doesn't cache minimizers, as that seems to be a
needless bottleneck. It could also be conditional on the number of CPUs
available.

* All sequences containing Ns are currently ignored.

* We should review how to deal with m <= len(match) < k.

* For caching efficiency, we should change the minimizer data structure to a
map-indexed 1D slice of Kmers (not *Kmers). (This technique originated in
Kraken.)

* Ensure that the stored Kraken data is referred to as the database, not
the library, to prevent confusion (see also no_write_lib).
