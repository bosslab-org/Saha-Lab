% Code provided by Saha Lab, Michigan State University


FOR TETRODE DATA:

Preprocess RMS data: MasterPreprocessTetrodes.m calls helper functions, read_intan_tetrodes.m and rms_construct_tetrodes.m
	1. read_intan_tetrodes.m- converts RHD files to RMS-transformed, position-specific and master MAT files
		a. Raw voltage traces high pass filtered at 300 Hz cutoff
	2. rms_construct_tetrodes.m
		a. 500 point continuous moving root-mean-square filter
		b. 500 point continuous moving average/smoothing filter
		c. Baseline subtraction
		d. Stimulus-specific baseline as trial- and sample-average for 2 seconds pre-stimulus
		e. Samples binned according to specified size (in msecs)
			aa) Need to construct separate files for different time windows of interest to avoid adjacent bin overlap effects
		f. Average of each bin calculated
		g. Tetrode channels averaged

Preprocess SS data: MasterPreprocessTetrodes.m calls helper functions, read_intan_tetrodes.m and rms_construct_tetrodes.m
	1. read_intan_tetrodes.m- converts RHD files to position-specific MAT files
		a. Raw voltage traces high pass filtered at 300 Hz cutoff
	2. mat_to_Text.m- converts MAT files to txt files
	3. IGOR Pro spike sorting
	4. TextToMat.m- converts spike sorted txt files to stimulus-specific, position-specific and master MAT files

Process data: MasterProcessLocust.m calls helper functions, PCA.m, LDA.m, TrainTestLocust.m and/or LeaveTrialOutLocust.m, Confusion.m
	1. Baseline subtraction (SS ONLY)
	2. Binned and averaged according to user specified bin size
	3. PCA.m- generates principal component analysis-based neural trajectories
	4. LDA.m- generates linear discriminant analysis-specific stimulus clusters
	5. Stimulus prediction
		a. TrainTestLocust.m- Separates out specified trials for training template and classifies remaining test trials
		b. LeaveTrialOutLocust.m- Cycles through trials using one trial data as training template and classifying remaining trials
	6. Confusion.m- generates confusion matrices for visualizing stimulus prediction accuracy

------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------

FOR INDEPENDENT CHANNELS DATA:
