===============
Genomic context
===============

This section summarizes the genomic context that reads have mapped to. This section
serves as a rough indicator to where reads have aligned to. The counting is naive:

* Counts are in terms of alignments. Thus a single read might
  contribute several counts to the same context or different contexts.

* Some genomic contexts can be overlapping, thus some alignments might
  be counted several times.

* An alignment needs to map to a context over at least 50% of its
  bases.  Thus some alignments spanning several contexts might be
  dropped.

All regions
===========

The following table lists all the genomic contexts that reads map to. 

.. report:: Context.ContextSummary
   :render: table
   :force:

   Number of alignments that align in a certain genomic context

Significance
------------

The plots below show enrichments or depletions for particular
genomic annotations. Those which are signifant (FDR=5%) are shown
in red.

.. report:: Gat.GatLogFoldContext
   :render: interleaved-bar-plot
   :groupby: track
   :colour: colour
   :layout: column-3
   :width: 300

   Log fold change. Those that are significant (qvalue < 0.05) are shown
   in red.

.. report:: Gat.GatSignificantResultsContext
   :render: table
   :force:

   Significant gat results at an FDR of 5%


Ribosomal regions
=================

Ribosomal RNA is one of the most abundant transcripts in a cell and
dominates RNASeq samples until it is removed. The following plots and
tables examine the number of intervals overlapping repetitive
RNA. Repetetive RNA annotation is taken from the UCSC repeatmasker
tracks.

.. report:: Context.ContextSummary
   :render: table
   :slices: total,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA,ribosomal_coding

   Number of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track

.. report:: Context.ContextSummary
   :render: pie-plot
   :pie-first-is-total: notRNA
   :groupby: track
   :slices: total,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA,ribosomal_coding
   :layout: column-3
   :width: 200

   Proportion of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track

Protein coding genes
====================

The following plots list the number of alignments to protein coding and (protein coding) 
pseudogene exons. The annotations are taken from the ENSEMBL gene set.

.. report:: Context.ContextSummary
   :render: table
   :slices: mapped,protein_coding,pseudogene

   Number of alignments that align to protein coding genes or pseudo genes.

.. report:: Context.ContextSummary
   :render: pie-plot
   :pie-first-is-total: genomic
   :groupby: track
   :slices: total,protein_coding,pseudogene
   :layout: column-3
   :width: 200

   Proportion of alignments that align to protein coding genes or pseudo genes.

Repeats
=======

The following plots list the number of alignments to protein coding and (protein coding) 
pseudogene exons. The annotations are taken from the ENSEMBL gene set.

.. report:: Context.ContextSummary
   :render: table
   :slices: mapped,repeats

   Number of alignments that align to repeats

.. report:: Context.ContextSummary
   :render: pie-plot
   :pie-first-is-total: genomic
   :groupby: track
   :slices: total,repeats
   :layout: column-3
   :width: 200

   Proportion of alignments that align to repeats

