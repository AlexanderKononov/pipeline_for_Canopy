import sys
import functions as f
import argparse as ap

# As input two file name should be noted: with SNA data and with CNA segmentation data.
# The third argument should be notes a work directory
# The fourth argument should note the name of normal non-tumor sample to filtration SNAs
parser = ap.ArgumentParser(description = "Create input file for canopy using list of SNA and CNA observations")
parser.add_argument("sna_file", help = "file with SNA data")
parser.add_argument("cna_file", help = "file with CNA data")
parser.add_argument("-o", "--output_dir", default = './', help = "directory to save output files")
parser.add_argument("--IDfilter", action = "store_true", help = "remove observation without any annatation")
parser.add_argument("-n", "--normal_sample", default = '', help = "name or part of name of normal tissue sample to excluding non-somatic mutations" )
parser.add_argument("--include", nargs = '*',default = [], help = "noted list samples which needeeto include to analiz")

arg = parser.parse_args()

if arg.output_dir[-1] != '/': arg.output_dir = arg.output_dir + '/'


print('---Start---')

# Preparing of SNV data
col_for_SNV, sample_name, col_norm = f.sampleNameFinderSNV(arg.sna_file, IDfilter = arg.IDfilter, NormFilter = arg.normal_sample, includes=arg.include)

print(sample_name)
print(col_for_SNV)
print('---sampleNameFinderSNV done---')
R_matrix, X_matrix, SNV_list = f.fullMatrixXR(arg.sna_file, col_for_SNV, IDfilter = arg.IDfilter, col_norm = col_norm)
print('---fullMatrixXR done---')
f.writeFileXRByMatrix(X_matrix, R_matrix, sample_name, X = arg.output_dir + 'Xout.txt', R = arg.output_dir + 'Rout.txt')
print('----SNV data were done---')

# Preparing of CNA data
target_samples_file = f.sampleFilerToCNAFile(arg.cna_file, sample_name, f_wName = arg.output_dir + 'CNA_for_target_patient.csv')
print('----CNA samples were filtrated---')
sampleCNA_list = f.sampleNameCNAListCreater(target_samples_file)
CNA_dict, CNA_list = f.fullDictCNA(target_samples_file, sampleCNA_list)
CNA_dict, CNA_list = f.redesigneCNADict(CNA_dict, CNA_list)
f.writeFileWsByDictCNA(CNA_dict, CNA_list, WM = arg.output_dir + 'WMout.txt', Wm = arg.output_dir + 'Wmout.txt')
print('----CNA data were done---')

# Preparing of CNA and SNV overlapping
overlap_dict, CNA_regions = f.finderCNAregioonsOverlaping(CNA_list)
f.createCByRegionsFile(CNA_list, CNA_regions, overlap_dict, f_wName = arg.output_dir + 'Cout.txt')

# Preparing of CNA and SNV overlapping
f.createYFile(SNV_list, CNA_regions, f_wName = arg.output_dir + 'Yout.txt')
print('----createYFile was done---')

