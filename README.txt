Scripts used for extendedgastruloid project and manuscript 'Extended culture of 2D gastruloids to model human mesoderm development'. 

Most imaging analyses and figures were generated using MATLAB version 2023b. 

Single cell RNA sequencing and part of the segmentations were done with python packages.

To set up the path of your MATLAB scripts as this repository, run the 'setup.m' script in the main folder.

Template scripts for image preprocessing, segmentation, and quantification of fixed sample data as shown in the manuscript are in the folder 'template_scripts'. We use both ilastik and cellpose (v1) for segmentation. After getting the cell masks from each z slices, run 'combineSegmentations.m' to consolidate the cell masks and final segmentation files will be saved in the designated folder where the data are stored. Then run 'singleCellQuantification.m' to get the quantitative values of the masks as the readout for imaging quantification. All the plots and visualizations are then generated with 'analyze_fixed_integrated_Z.m'. 'external' and 'image_analysis' folders contain the functions that are used for the imaging analysis mentioned above of the fixed data . 

Scripts for live and single cell tracking analysis are in the folder 'single_cell_tracking_analysis'. Scripts for singe cell RNA sequencing analysis are in the folder 'scRNA-seq analysis'.

   