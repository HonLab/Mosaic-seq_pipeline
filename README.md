Description
=============
Mosaic-seq is the technique developed in Gary Hon Lab which allows acquisition of sgRNA information and transcriptome simultaneously from the same single cell during single-cell RNA-seq. This is the pipeline used in Gary Hon Lab for Mosaic-seq analysis. 

Authors
==============
Gary Hon (@GaryhonLab), Jialei Duan (@jlduan), Shiqi Xie (@russellxie)

Correspondence: Gary Hon (gary.hon@utsouthwestern.edu)

Notice
===============
This pipeline is written to satify the computational environment in the BioHPC system of UT Southwestern Medical Center, therefore customization is required if you want to run this pipeline in your own system. 

Requirements: 
==============
  * Drop-seq pipeline: http://mccarrolllab.com/dropseq/
  * star: https://github.com/alexdobin/STAR
  * picard: https://broadinstitute.github.io/picard/
  * samtools: http://samtools.sourceforge.net/
  * Matlab
  * Python package requirement: Numpy, Scipy, Matplotlib, pysam
  * R package requirement: ggplot2, reshape2
    
General Workflow:
=================

