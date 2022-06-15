# Saha Lab- Michigan State University

This repository contains all code for processing stimulus evoked neural signals from antennal lobe projection neurons.



1. Data preprocessing and saving:
	Adjust the experimental parameters in the file master_preprocess.m. This script invokes helper functions, including read_intan.m and rms_construct.m, which ought to be on the same local or system path.
		read_intan.m: preprocesses RHD files and saves stimulus-specific, position-specific MAT files
		rms_construct.m: saves master file of root-mean square transformed data


2. Data processing and analysis:
	Adjust the experimental parameters in the file master_process.m. This script invokes a number of helper functions, which are dependent upon desired analytical techniques/experimental parameters.


	Single day:
		PSVT.m: Peri-stimulus voltage trace
		Raster.m: Raster plot
		PSTH.m: Peri-stimulus time histogram
	Population-based:
		PCA.m: Principal component analysis
		LDA.m: Linear discriminant analysis
		Confusion.m: Confusion matrix
		
Neural network- ANNs and SNNs