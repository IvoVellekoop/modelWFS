"Model-based wavefront shaping microscopy"

Matlab scripts to generate all the data underlying the article "Model-based wavefront shaping microscopy" [ Correspondence to: a.thendiyammal@utwente.nl]

**********************************************************************************************
Comparison of Model-based wavefront shaping and feedback-based wavefront shaping

Running sequence:

1. Generate refractive index map from TPM intensity images.

   (a) Run /refractive_index_mapping/ShearTPMimage.m for coverting the acquired TPM frames to actual size in micrometer. 
   
       Input file: Image from ScanImage (Here, eg. PDMS_diffuser_surface_1X512.tif) 
       Output file: TPM_3D.m

   (b) Run Refractive_Indexmap_3D.m to generate refractive index distribution. 
   
       Required input: TPM_3D from the previous step.
       Output file: Sample properties (eg. n_sample, PDMS_thickness)

2. Run Model_WFS_PhaseConjugation.m for simulating light propagation through the PDMS diffuser to find correction wavefront.

       Input: SampleProperties.m 
       Output: SLMCorrection.m

3. Run feedback_optimization.m for finding the optimal wavefront.

       Output: ideal_wavefront.m

4. Run WFScomparison.m for generating the Phase pattern on the SLM and compare the model-based wavefront shaping with 
   feedback-based wavefront shaping.
 
   Inputs: SLMCorrection and ideal_wavefront

**********************************************************************************************


Data Analysis codes:

1. ModelWFS_DataStitching.m is used to stitch 3D substacks to form a combined 3D stack. 

2. ModelWFS_DataAnalysis.m - Generate maximum intensity projection. Compute mean intensity of the beads as a function of depth.

3. plotSLMpatterns: Generate SLM patterns at different depths. Four patterns corresponding to four depths are used in the paper.

4. TPM_Image_Calibration: This code is used to find the conversion matrix which converts TPM frames to correct size in micrometer with correct orientation.
   This code also generates SLM coordinates corresponding to a spatial frequency (kNA) of the pupil plane. 
   