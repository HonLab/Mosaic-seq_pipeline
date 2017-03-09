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
  * [Drop-seq pipeline](http://mccarrolllab.com/dropseq/)
  * [Star](https://github.com/alexdobin/STAR)
  * [Picard](https://broadinstitute.github.io/picard/)
  * [Samtools](http://samtools.sourceforge.net/)
  * [Bedtools](http://bedtools.readthedocs.io/en/latest/)
  * Matlab 2016a
  * Python package requirement: [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/), [Pysam](http://pysam.readthedocs.io/en/latest/api.html)
  * R package requirement: [ggplot2](http://ggplot2.org/), [reshape2](https://cran.r-project.org/web/packages/reshape2/reshape2.pdf)
    
General Workflow:
=================
