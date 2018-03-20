import sys
import functions

# As input two file name should be noted: with SNA data and with CNA segmentation data.
arg = sys.argv
print('---Start---')

# Preparing of SNV data
col_for_SNV, sample_name = functions.sampleNameFinderSNV(arg[1])
print('---sampleNameFinderSNV done---')
R_matrix, X_matrix, SNV_list = functions.fullMatrixXR(arg[1], col_for_SNV)
print('---fullMatrixXR done---')
functions.writeFileXRByMatrix(X_matrix, R_matrix, sample_name)
print('----SNV data were done---')

# Preparing of CNA data
target_samples_file = functions.sampleFilerToCNAFile(arg[2], sample_name)
print('----CNA samples were filtrated---')
sampleCNA_list = functions.sampleNameCNAListCreater(target_samples_file)
CNA_dict, CNA_list = functions.fullDictCNA(target_samples_file, sampleCNA_list)
CNA_dict, CNA_list = functions.redesigneCNADict(CNA_dict, CNA_list)
functions.writeFileWsByDictCNA(CNA_dict, CNA_list)
print('----CNA data were done---')

# Preparing of CNA and SNV overlapping
overlap_dict, CNA_regions = functions.finderCNAregioonsOverlaping(CNA_list)
functions.createCByRegionsFile(CNA_list, CNA_regions, overlap_dict)

# Preparing of CNA and SNV overlapping
functions.createYFile(SNV_list, CNA_regions)
print('----createYFile was done---')

