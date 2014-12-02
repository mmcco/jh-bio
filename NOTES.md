* The first round of pointer method calls had no significant effect on
  runtime.
* When moving the Match struct into the bioutils library, I remembered
  that it contained pointers to RepeatGenome structs like Repeat and
  ClassNode. When finding the references to these in the RepeatGenome
  code, there were a few places in which they were written or printed,
  and quite a few places in which they were set or updated, but only a
  single instance (in RepeatGenome.getMatchKmers) where one was actually
  used.
