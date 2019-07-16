# LSx
A script in R to run the LS³ and LS⁴ phylogenetic data subsampling algorithms for reducing lineage rate heterogeneity

For the R script and LS⁴, please cite:
Rivera-Rivera CJ and Montoya-Burgos JI (accepted). LSx: Automated reduction of gene-specific lineage evolutionary rate heterogeneity for multi-gene phylogeny inference. BMC Bioinformatics

For LS³, please also cite:
Rivera-Rivera CJ and Montoya-Burgos JI (2016). LS³: A Method for Improving Phylogenomic Inferences When Evolutionary Rates Are Heterogeneous among Taxa. Molecular Biology and Evolution 33(6):1625-1634

1 Quick Start

1.1 Download and installation


  Download the R script LSx_v1.1.R and the LSX input file LSx_input_file.txt from
https://genev.unige.ch/research/laboratory/Juan-Montoya (under the “MORE” tab).


  DONE! (LSX runs directly through the Rscript or Rscript.exe scripting front-end, see
point 1.3)
  Check dependencies in point 3.
  
  
1.2 Input files

  • Gene sequence alignments, each as an independent file in PHYLIP interleave
format. Taxa names of less than 10 characters are identical across gene datasets,
unique within a gene dataset, and no name is a subset of another name (format
details in section 4.1).

  • The list of gene alignment files to be analyzed, with (optionally) the model of
sequence evolution; in raw text, with comma-delimited columns if needed (format
details in section 4.2).

  • Guide tree in Newick format, rooted, with a polytomy at the base of the lineages of
interest (format details in section 4.3).

  • Lineage-taxon file in raw text format, with comma-delimited columns (format
details in section 4.4).

  • LSX input file with the values for your current run (details in section 4.5).
  
  • PAML control files. These are baseml.ctl or codeml.ctl (aaml.ctl in the LSX
example files).


1.3 Running

  Linux/Mac command:
  
    Rscript LSx_v1.1.R LSx_input_file.txt
    
  Windows command:
  
    Rscript.exe LSx_v1.1.R LSx_input_file.txt
