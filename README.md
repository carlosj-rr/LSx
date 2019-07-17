# LS<sup>X</sup>

LS<sup>X</sup> is a script in R that runs the LS³ and LS⁴ algorithms of data subsampling for multigene phylogenetic inference. Both of these algorithms do a gene-by-gene inspection of the heterogeneity of evolutionary rates among user-defined lineages of interest (LOI). Then, using criteria that differ in both algorithms (see details [here](https://github.com/carlosj-rr/LSx/wiki/Introduction#is-lsx-for-me) or in the [papers](https://github.com/carlosj-rr/LSx/wiki/Citations)), they try to find a subsample of sequences that evolve at a homogeneous rate across all LOIs. If this subset is found, an alignment of the gene is produced with only the sequences that evolve homogeneously. At the same time, a table is also produced showing which sequences were “flagged” (the sequences that were removed), and which sequences were kept. If a subset of sequences that evolve at a homogeneous rate is not found, the gene is flagged entirely.

____ 

## [Manual](https://github.com/carlosj-rr/LSx/wiki)  

----

## Published uses
* Long Branch Attraction (LBA) in Metazoa and Mammalia (LS³-bash)
  * https://doi.org/10.1093/molbev/msw043
* LBA in catfish phylogenetics (LS³ and LS⁴ with LS<sup>X</sup>)
  * https://doi.org/10.1016/j.ympev.2018.06.004
* Benchmarks and presentation of LS⁴ and LS<sup>X</sup>
  * https://doi.org/10.1101/220053
