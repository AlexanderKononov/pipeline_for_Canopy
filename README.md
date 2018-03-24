# pipeline_for_Canopy
**pipeline_for_Canopy** is set of scripts to prepare WGS data from HERCULES project to Canopy analysis and  to run this analyze.

There is python script *file_adapter_to_Canopy.py* which can used to formating and preprosessing data. To run it, 3 or 4 arguments should be noted in command line during run.
1. First argument is way to SNA file (for example SNA from MuTect2 GATK)
2. Second argument is wey to CNA file, file with data about a segmentation of DNA.
3. Third argument is way to the work directory, way to the output matrices.
4. By fourth argument can be noted a name of normal sample. Normal sample and it's data should be listed in SNA file as other tumor samples of patient. If fourth argument wasn't noted the script take into account all SANs position from input (doesn't filter it to normal polymorphisms). 

`$ python file_adapter_to_Canopy.py ./test_data/SNV_test.csv ./test_data/CNA_test.csv ./test_data/work_dir H002_B1144`
Examples of input file are presented in *test_data* directory. *H002_B1144* is name of normal sample in this example.

*file_adapter_to_Canopy.py* creates a matrices needed to Canopy analyze  (X, R, WM, Wm, Y, C) http://www.pnas.org/content/113/37/E5528 

There is R script *pipeline_Canopy_with_clustering.R* which can used to run Canopy analyze.
If you plan to use Rscript  the two arguments should be noted in command line during run.
1. Work directory, where all matrices place.
2. Project name.

To optimize analyze you can optimized different parameters inside this script (num_cluster, num_ru, K, numchain etc.)

