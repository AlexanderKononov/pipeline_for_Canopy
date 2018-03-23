import sys
import functions as f

# As input two file name should be noted: with SNA data and with CNA segmentation data.
# The third argument should be notes a work directory
# The fourth argument should note the name of normal non-tumor sample to filtration SNAs
arg = sys.argv
if arg[3][-1] != '/':
    arg[3] = arg[3] + '/'
print('---Start---')
if len(arg) < 5:
    arg.append('')

# Preparing of SNV data
col_for_SNV, sample_name, col_norm = f.sampleNameFinderSNV(arg[1], IDfilter = True, NormFilter = arg[4])
print('---sampleNameFinderSNV done---')
R_matrix, X_matrix, SNV_list = f.fullMatrixXR(arg[1], col_for_SNV, IDfilter = True, col_norm = col_norm)
print('---fullMatrixXR done---')
f.writeFileXRByMatrix(X_matrix, R_matrix, sample_name, X = arg[3] + 'Xout.txt', R = arg[3] + 'Rout.txt')
print('----SNV data were done---')

# Preparing of CNA data
target_samples_file = f.sampleFilerToCNAFile(arg[2], sample_name, f_wName = arg[3] + 'CNA_for_target_patient.csv')
print('----CNA samples were filtrated---')
sampleCNA_list = f.sampleNameCNAListCreater(target_samples_file)
CNA_dict, CNA_list = f.fullDictCNA(target_samples_file, sampleCNA_list)
CNA_dict, CNA_list = f.redesigneCNADict(CNA_dict, CNA_list)
f.writeFileWsByDictCNA(CNA_dict, CNA_list, WM = arg[3] + 'WMout.txt', Wm = arg[3] + 'Wmout.txt')
print('----CNA data were done---')

# Preparing of CNA and SNV overlapping
overlap_dict, CNA_regions = f.finderCNAregioonsOverlaping(CNA_list)
f.createCByRegionsFile(CNA_list, CNA_regions, overlap_dict, f_wName = arg[3] + 'Cout.txt')

# Preparing of CNA and SNV overlapping
f.createYFile(SNV_list, CNA_regions, f_wName = arg[3] + 'Yout.txt')
print('----createYFile was done---')

