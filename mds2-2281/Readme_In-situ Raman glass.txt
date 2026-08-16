The Readme file 'Readme_In-situ Raman glass.txt' was generated on 20200724 by Y.B. Gerbig and C.A. Michaels
-------------------
GENERAL INFORMATION
-------------------
1. Title of Dataset: 
	In-situ Raman spectra of indented fused silica  

2. Formats of the files: 
	csv and tiff

3. Author Information
  Principal Investigator Contact Information
        Name: Yvonne B. Gerbig
		Institution: Material Measurement Laboratory, NIST 
		Address: 100 Bureau Dr., Gaithersburg, MD 20899
		Email: yvonne.gerbig@nist.gov
  
  Co-investigator Contact Information
        Name: Chris A. Michaels
		Institution: Material Measurement Laboratory, NIST
           	Address: 100 Bureau Dr., Gaithersburg, MD 20899
           	Email: chris.michaels@nist.gov
  
  Alternate Contact Information
           Name:
           	Institution:
           	Address:
           	Email:

4. Date of data collection: 
	2018/11 to 2019/02

5. Geographic location of data collection: 
	NIST Gaithersburg, MD USA

--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 
6. Other than the NIST statements for Copyright, Fair Use, and Licensing found at 
https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software, list the other license or restrictions which are placed on this data:
	None

7. Are there any restrictions or guidelines on how to use the data, (ex. information about access restrictions based on privacy, security, or other policies)?
	No

8. Is there a documentary standard that applies to this dataset? 
	No

9. Links to publications that cite or use the data: 
	https://doi.org/10.1016/j.jnoncrysol.2019.119828

10. Links to other publicly accessible locations of the data:
	None

11. Links/relationships to ancillary data sets:
	None

12. Was data derived from another source? (example, Open FEMA)
	No

13. Recommended citation for the data (See NIST guidance for citation https://inet.nist.gov/nvl/howdoi/cite-nist-datasets): 
	https://doi.org/10.18434/M32281

---------------------
DATA & FILE INFORMATION
---------------------
14. File List
A. Dataset: Z profile_Raman spectra
	Short description:
 		Raman spectra collected in-situ in the center of the indentation imprint on the FS specimen with the objective focused on the surface of the sample and with the focal plane shifted 1 µm to 7 µm away from that position along the z axis into the specimen bulk 
		for different objective/specimen thickness pairings:
		(a) CO/thin FS, 
		(b) MO/thin FS, 
		(c) TO/thick FS, 
		(d) MO/thick FS,

	Format of the files: 
		(a-d) 10 column csv files

	Column assignments: 
		(a-d) C1: Wavenumber, C2: Intensity (counts) of baseline spectrum, C3 through C10: Intensities (counts) at z positions in order from the specimen surface to 7 micrometers into the specimen bulk.

	Raw data files:
		(a) CO_thin_FS_Zprofile_raw_data.csv
      		(b) MO_thin_FS_Zprofile_raw_data.csv
       		(c) TO_thick_FS_Zprofile_raw_data.csv
       		(d) MO_thick_FS_Zprofile_raw_data.csv

B. Dataset: Z profile_indentation curves
	Short description: 
		Force and displacement data recorded as function of time during the indentation of the FS specimens in the z-profile experiments for different objective/specimen thickness pairings:
		(a) CO/thin FS, 
		(b) MO/thin FS, 
		(c) TO/thick FS, 
		(d) MO/thick FS,

	Format of the files: 
		(a-d) 3 column csv files

	Column assignments: 
		(a-d) C1: Time (s), C2: Indentation force (mN), C3: Displacement (nm).
		
	Raw data files:
		(a) CO_thin_FS_Zprofile_indentation.csv
		(b) MO_thin_FS_Zprofile_indentation.csv
		(c) TO_thick_FS_Zprofile_indentation.csv
		(d) MO_thick_FS_Zprofile_indentation.csv

C. Dataset: Load sequence_TO_Raman spectra
	Short description:
		Raman spectra collected in-situ with the TO during the indentation of the thick FS specimen in load-sequence experiment:
		(a) prior to the indentation (0 mN) and at different indentation loads in loading segment,	
		(b) in the residual contact impression (0 mN) and at different indentation loads in unloading segment.

	Format of the files: 
		(a) 8 column csv file
		(b) 7 column csv file

	Column assignments: 
		(a) C1: Wavenumber, C2 through C8: Intensities (counts) at indentation forces in ascending order from 0 mN to 300 mN.
		(b) C1: Wavenumber, C2 through C7: Intensities (counts) at indentation forces in ascending order from 0 mN to 250 mN.

	Raw data files:
		(a) TO_load_curve_raw_data.csv
		(b) TO_unload_curve_raw_data.csv

D. Dataset: Load sequence_TO_indentation curve
	Short description: 
		Force and displacement data recorded as function of time during the indentation of the thick FS specimen in the load-sequence experiment with TO configuration.

	Format of the file: 
		3 column csv file
	
	Column assignments: 
		C1: Time (s), C2: Indentation force (mN), C3: Displacement (nm).
	
	Raw file name: 
		TO_thick_FS_loadsequence_indentation.csv

E. Dataset: Load sequence_TO_WL images
	Short description: 
		White light images of the contact area of diamond probe and thick FS specimen collected in the load-sequence experiment for TO configuration:
		 - during loading at 50 mN, 100mN, 150 mN, 200 mN, 250 mN and 300 mN,
		 - during unloading at 250 mN, 200 mN, 150 mN, 100 mN and 50 mN, 
		 - after complete unloading (residual) of the residual indentation. 
		For calibration purposes, a white light image of a Ronchi grating (100 lines per mm) collected in the TO configuration was added.
	
	Format of the files: 
		tiff files (2000 pixels x 1600 pixels)

	Raw Data Files: 
		1. Loading_50mN_TO.tiff
		2. Loading_100mN_TO.tiff
		3. Loading_150mN_TO.tiff
		4. Loading_200mN_TO.tiff
		5. Loading_250mN_TO.tiff
		6. Loading_300mN_TO.tiff
		7. Unloading_250mN_TO.tiff
		8. Unloading_200mN_TO.tiff
		9. Unloading_150mN_TO.tiff
		10. Unloading_100mN_TO.tiff
		11. Unloading_50mN_TO.tiff
		12. Residual_TO.tiff
		13. Ronchi grating_TO.tiff

F. Dataset: Load sequence_CO_Raman spectra
	Short description: 
		Raman spectra collected in-situ with the CO during the indentation of the thin FS specimen in load-sequence experiment:
		(a) prior to the indentation (0 mN) and at different indentation loads in loading segment,	
		(b) in the residual contact impression (0 mN) and at different indentation loads in unloading segment.

	Format of the files: 
		(a) 8 column csv file
		(b) 7 column csv file

	Column assignments: 
		(a) C1: Wavenumber, C2 through C8: Intensities (counts) at indentation forces in ascending order from 0 mN to 300 mN.
		(b) C1: Wavenumber, C2 through C7: Intensities (counts) at indentation forces in ascending order from 0 mN to 250 mN.

	Raw data files: 
		(a) CO_load_curve_raw_data.csv
		(b) CO_unload_curve_raw_data.csv 

G. Dataset: Load sequence_CO_Indentation curve
	Short description:
 		Force and displacement data recorded as function of time during the indentation of the thin FS specimen during the load-sequence experiment in CO configuration.
	
	Format of the file: 
		3 column csv file

	Column assignment: 
		C1: Time (s), C2: Indentation force (mN), C3: Displacement (nm).

	Raw data file: 
		CO_thin_FS_loadsequence_indentation.csv

H. Dataset: Load sequence_CO_WL images
	Short description: 
		White light images of the contact area of diamond probe and thin FS specimen collected in the load-sequence experiment for CO configuration:
		 - during loading at 50 mN, 100 mN, 150 mN, 200 mN, 250 mN and 300 mN, 
		 - during unloading at 250 mN, 200 mN, 150 mN, 100 mN and 50 mN, 
		 - after complete unloading (residual) of the residual indentation. 
		For calibration purposes, a white light image of a Ronchi grating (100 lines per mm) collected in the CO configuration was added.
 
	Format of the files: 
		tiff files (2000 pixels x 1600 pixels)

	Raw data files: 
		1. Loading_50mN_CO.tiff
		2. Loading_100mN_CO.tiff
		3. Loading_150mN_CO.tiff
		4. Loading_200mN_CO.tiff
		5. Loading_250mN_CO.tiff
		6. Loading_300mN_CO.tiff
		7. Unloading_250mN_CO.tiff
		8. Unloading_200mN_CO.tiff
		9. Unloading_150mN_CO.tiff
		10. Unloading_100mN_CO.tiff
		11. Unloading_50mN_CO.tiff
		12. Residual_CO.tiff
		13. Ronchi grating_CO.tiff

I. Filename: Figure 1 (Raman Intensity vs wavenumber for different probes)
	Short description: 
		Raman spectra in range 150 cm-1 to 1350 cm-1 for 
		(1a) thin FS (FS only) and 
		(1b) through (1d) for thin FS with three different diamond probes (FS+Probe1, FS+Probe2, FS+Probe3) hovering near the sample surface.     
	
	Format of the file: 
		6 column csv file
	
	Column assignments: 
		C1: Wavenumber for (a) and (b), C2: (a) Normalized intensity (FS only), C3: (b) Normalized intensity (FS+Probe1), C4: Wavenumber for (c) and (d), C5: (c) Normalized intensity (FS+Probe2), C6: (d) Normalized intensity (FS+Probe3)

	Specific file name: 
		Dataset_Fig1.csv
	
	Format of raw data file: 
		8 column csv file

	Column assignment of raw data file: 
		C1: Wavenumber (a), C2:Intensity (a), C3: Wavenumber (b), C4:Intensity (b), C5:Wavenumber (c), C6:Intensity (c), C7: Wavenumber (d), C8:Intensity (d).

	Raw data file: 
		Figure_1_raw_data.csv

J. Dataset: Figure 2 (Intensity vs wavenumber for different objectives)
	Short description: 
		Raman spectra in range 150 cm-1 to 1450 cm-1 for FS (FS only) and FS with the same diamond probe hovering near the sample surface (FS+probe) for four objective/specimen thickness pairings: 
		(2a) CO/thin FS, 
		(2b) MO/thin FS, 
		(2c) TO/thick FS, 
		(2d) MO/thick FS,

	Format of the files: 
		3 column csv files

	Column assignments: 
		C1: Wavenumber, C2: Normalized intensity (FS only), C3: Normalized intensity (FS+probe).
	
	Specific file names: 
		(2a) Dataset_Fig2a.csv
		(2b) Dataset_Fig2b.csv
		(2c) Dataset_Fig2c.csv
		(2d) Dataset_Fig2d.csv
	
	Format of raw data file: 
		12 column csv file

	Column assignments of raw data file: 
		C1: Wavenumber (a), C2: Intensity (a), C3: Intensity (a), C4: Wavenumber (b), C5: Intensity (b), C6: Intensity (b), C7: Wavenumber (c), C8: Intensity (c), C9:Intensity (c), C10: Wavenumber (d), C11:Intensity (d), C12: Intensity (d).

	Raw data file:
		(2a-d) Figure_2_raw_data.csv


K. Dataset: Figure 3 (Raman spectra as function of z for MO and CO configurations)
	Short description: 
		(3a) CO Raman spectra in range 150 cm-1 to 1350 cm-1 for thin fused silica recorded at a series of 8 z focal positions ranging from the surface to 7 micrometers below the surface along with a fused silica spectrum recorded with no probe present. 
		(3b) MO Raman spectra in range 150 cm-1 to 1350 cm-1 recorded at the specimen surface for thin fused silica in its unperturbed state (probe hovering above surface) and while being indented at 300 mN (collected in center of contact impression). 
		(3c) Center mass wavenumber of main band plotted against the z axis focal position for the CO and MO configurations for thin fused silica.    

	Format of the files: 
		(3a) 15 column csv file 
		(3b) 3 column csv file
		(3c) 5 column csv file 

	Column assignments: 
		(3a) C1: Wavenumber, C2 through C9: Normalized intensity spectra of indented FS for focal positions, C11 through C15: Normalized intensity for baseline spectrum and its offset copies for focal positions in ascending order. 
		(3b) C1: Wavenumber, C2: Normalized intensity for surface spectrum, C3: Normalized intensity for baseline. 
		(3c) C1: Focal position, C2: Center mass wavenumber in surface spectrum for CO, C3: Center mass wavenumber in surface spectrum for MO, C4: Center mass wavenumbers in the baseline spectrum for CO, C5: Center mass wavenumbers in the baseline spectrum for MO.

	Specific file names: 
		(3a) Dataset_Fig3a.csv
		(3b) Dataset_Fig3b.csv
		(3c) Dataset_Fig3c.csv

	Format of the raw data files: 
		(3a, 3bc) 10 column csv files
	
	Column assignments of raw data files: 
		(3a, 3bc) C1: Wavenumber, C2: Intensity (counts)of baseline spectrum, C3 through C10: Intensities (counts) at z positions in order from the specimen surface to 7 micrometers into the specimen bulk.

	Raw data files:
		(3a) Figure_3a_raw_data.csv
		(3b-c) Figure_3bc_raw_data.csv

L. Dataset: Figure 4 (Raman spectra as function of z for TO and CO configurations)
	Short description: 
		(4a) TO Raman spectra in range 150 cm-1 to 1350 cm-1 recorded at the specimen surface for thick fused silica in its unperturbed state (probe hovering above surface) and while being indented at 300 mN (collected in center of contact impression).  
		(4b) MO Raman spectra in range 150 cm-1 to 1350 cm-1 recorded at the specimen surface for thick fused silica in its unperturbed state (probe hovering above surface) and while being indented at 300 mN (collected in center of contact impression). 
		(4c) Center mass wavenumber of the main band plotted against the z axis focal position for the TO and MO configurations.    

	Format of the files: 
		(4a) 3 column csv file 
		(4b) 3 column csv file 
		(4c) 5 column csv file
	
	Column assignments: 
		(4a) C1: Wavenumber, C2: Normalized intensity (surface spectrum), C3: Normalized intensity (baseline).
		(4b) C1: Wavenumber, C2: Normalized intensity (surface spectrum), C3: Normalized intensity (baseline).
		(4c) C1: Focal position, C2: Center mass wavenumber of TO, C3: Center mass wavenumber of MO, C4: Center mass wavenumbers (TO baseline spectrum), C5: Center mass wavenumbers (MO baseline spectrum).
	
	Specific file names: 
		(4a) Dataset_Fig4a.csv
		(4b) Dataset_Fig4b.csv
		(4c) Dataset_Fig4c.csv
       	 		
	Format of the raw data files: 
		(4ac, 4bc) 10 column csv files

	Column assignments of raw data files: 
		(4ac, 4bc) C1: Wavenumber, C2: Intensity (counts) of baseline spectrum, C3 through C10: Intensities (counts) at z positions in order from the specimen surface to 7 micrometers into the specimen bulk.
	
	Raw data files:
       		(4ac) Figure_4ac_raw_data.csv
		(4bc) Figure_4bc_raw_data.csv

M. Dataset: Figure 5 (Raman spectra as function of z for CO and TO configurations in 300 mN indentation)
	Short description:  
		(5a) CO and TO configurations Raman spectra in range 150 cm-1 to 1350 cm-1 for thin and thick fused silica, respectively, recorded at the specimen surface in the unperturbed state (probe hovering above surface) and while being indented at 300mN (collected in the center of contact impression). 	
		(5b) Relative shift of main band peak wavenumber (MP) and center of mass (CM) plotted against the z axis focal position for the CO and TO configurations including trendlines of the combined (MP+CM) data for the CO and TO.   

	Format of the files:
		(5a) 6 column csv file 
		(5b) 8 column csv file  

	Column assignments: 
		(5a) C1: Wavenumber (for C2 and C3), C2: Normalized intensity of spectrum collected at 300 mN for CO/thin FS, C3: Normalized intensity (baseline for CO/thin FS), C4: Wavenumber (for C5 and C6), C5: Normalized intensity of spectrum collected at 300 mN for TO/thick FS, C6: Normalized intensity (baseline for TO/thick FS). 
		(5b) C1: Focal position, C2: Relative MP shift (CO/thin FS), C3: Relative CM shift (CO/thin FS), C4: Relative MP shift (TO/thick FS), C5: Relative CM shift (TO/thick FS), C6: Horizontal axis data of trendlines, C7: Vertical axis data of trendline (CO/thin FS and TO/thick FS), C8: Vertical axis data of trendline (TO/thick FS).

	Specific file names: 
		(5a): Dataset_Fig5a.csv
		(5b): Dataset_Fig5b.csv

       	Format of the raw data file: 
		(5ab) 20 column csv file

	Column assignments of raw data file: 
		(5ab) C1: Wavenumber CO, C2: Intensity (counts) of baseline CO spectrum, C3 through C10: Intensities (counts) at z positions in order from the specimen surface to 7 micrometers into the specimen bulk for CO objective.
                       C11: Wavenumber TO, C12: Intensity (counts) of baseline TO spectrum, C13 through C20: Intensities (counts) at z positions in order from the specimen surface to 7 micrometers into the specimen bulk for TO objective.
       
	Raw data file:
       		(5ab) Figure_5ab_raw_data.csv		


N. Dataset: Figure 6 (Raman spectra as function of increasing load for CO configuration)
	Short description: 
		(6a) CO Raman spectra in range 150 cm-1 to 1350 cm-1 for thin fused silica recorded in center of contact impression for a series of 5 increasing loads ranging from 0 mN to 300 mN (laoding segment of load-sequence experiment). 
		(6b) Relative shift of main band peak wavenumber (MP) and center of mass (CM) plotted against the indentation force including trendlines of the MP and CM data.
	
	Format of the files: 
		(6a) 6 column csv file.
		(6b) 6 column csv file.
	
	Column assignments: 
		(6a) C1: Wavenumber, C2 through C6: Normalized intensities for indentation forces in ascending order. 
		(6b) C1: Indentation force, C2: Relative MP shift, C3: Relative CM shift, C4: Horizontal axis data of trendlines, C5: Vertical axis data of MP trendline, C6: Vertical axis data of CM trendline.
	
	Specific file names: 
		(6a) Dataset_Fig6a.csv
		(6b) Dataset_Fig6b.csv
       
 	Format of the raw data files: 
		(6ab) 8 column csv files

	Column assignments of raw data files: 
		(6ab) C1: Wavenumber, C2 through C8: Intensities (counts) at indentation forces in ascending order from 0 mN to 300 mN.

	Raw data file:
       		(6ab) Figure_6ab_raw_data.csv

O. Dataset: Figure 7 (Correlations between indentation force and contact pressure and area for loading and unloading sequence)
	Short description:  
		Plot showing changes in contact area and contact pressure with increasing (loading) and decreasing force (unloading) for indentation with a pyramidal diamond probe into fused silica (including corresponding trendlines) determined for load-sequence experiment with CO configuration.

	Format of the file: 
		13 column csv file

	Column assignments: 
		C1: Indentation force, C2: Contact area (loading), C3: Contact area (unloading), C4: Contact pressure (loading), C5: Contact pressure (unloading), C6: Horizontal axis values of force-contact area trendline (loading), C7: Vertical axis values of force-contact trendline (loading), C8: Horizontal axis values of pressure-force trendline (loading), C9: Vertical axis values of pressure-force trendline (loading), C10: Horizontal axis values of force-contact area trendline (unloading), C11: Vertical axis values of force-contact area trendline (unloading), C12: Horizontal axis values of  pressure-force trendline (unloading), C13: Vertical axis values of  pressure-force trendline (unloading).

	Specific file name: 
		Dataset_Fig7.csv
		
P. Dataset: Figure 8 (Raman spectra as function of decreasing load for CO configuration)
	Short description: 
		(8a) CO Raman spectra in range 150 cm-1 to 1350 cm-1 for thin fused silica recorded in center of contact impression for a series of 5 decreasing loads ranging from 300 mN to 0 mN (unloading segment of load-sequence experiment). 
		(8b) Relative shift of main band peak wavenumber and center of mass plotted against the indentation force including trendlines.   

	Format of the files: 
		(8a) 6 column csv file 
		(8b) 5 column csv file 

	Column assignments: 
		(8a) C1: Wavenumber, C2 through C6: Normalized intensities for indentation forces in descending order. 
		(8b) C1: Indentation force, C2: Relative MP shift, C3: Relative CM shift, C4: Vertical axis data of MP trendline, C5: Vertical axis data of CM trendline. (The vertical axis values for the trendlines are contained in column 1).

	Specific file names: 
		(8a) Dataset_Fig8a.csv
		(8b) Dataset_Fig8b.csv
	
        Format of the raw data files: 
		(8ab) 7 column csv files

	Column assignments of raw data files: 
       		(8ab) C1: Wavenumber, C2 through C7: Intensities (counts) at indentation forces in ascending order from 0 mN to 250 mN.
       
       Raw data file:
       		(8ab) Figure_8ab_raw_data.csv


15. Relationship between files:        
	This document is a list of the datasets that serve as a companion to: Y.B. Gerbig and C.A. Michaels, J. Non-Cryst. Solids 530 119828 (2020), https://doi.org/10.1016/j.jnoncrysol.2019.119828.
	Datasets A through H contain the raw data files on which the cited publication is based.
	Datasets I through P contain the specific data depicted in Figures 1 through 8 of the cited publication both in raw form and plotted form.
	Dataset A (Raman spectra) and dataset B (indentation curves) contain the experimental data collected during the z-profile experiments with different objective/sample pairings.
	Dataset C (Raman spectra), dataset D (indentation curves) and dataset E (images) contain the experimental data collected in the load-sequence experiments in the TO configuration on the thick FS specimen.
	Dataset F (Raman spectra), dataset G (indentation curves) and dataset H (images) contain the experimental data collected in the load-sequence experiments in the CO configuration on the thin FS specimen.
	 
	
16. Additional related data collected that was not included in the current data package:
	No

17. Are there multiple versions of the dataset? 
	No
	
--------------------------
METHODOLOGICAL INFORMATION
--------------------------
18. Description of methods used for collection/generation of data: 
	Details about the collection/generation of the datasets are described in Y.B. Gerbig and C.A. Michaels, J. Non-Cryst. Solids 530 119828 (2020), https://doi.org/10.1016/j.jnoncrysol.2019.119828.

19. Methods for processing the data: 
	Details about processing of the spectral datasets are described in Y.B. Gerbig and C.A. Michaels, J. Non-Cryst. Solids 530 119828 (2020), https://doi.org/10.1016/j.jnoncrysol.2019.119828..
	All of the Raman spectra presented here were corrected for instrumental response nonuniformity.  The indentation curves (datasets B, D and G) were corrected for instrument-specific machine compliance and experiment-specific thermal drift.
	
20. Instrument- or software-specific information needed to interpret the data:
	NA

21. Standards and calibration information, if appropriate:
	The standards and calibrations procedures utilized in the handling of this data are described in detail in Y.B. Gerbig and C.A. Michaels, J. Non-Cryst. Solids 530 119828 (2020), https://doi.org/10.1016/j.jnoncrysol.2019.119828. 

22. Environmental/experimental conditions:
	Experiments were conducted in an environment-controlled laboratory (room temperature: 20°C +/- 2°C, relative humidity: 50 % +/- 5%).

23. Describe any quality-assurance procedures performed on the data:
	NA
