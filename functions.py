def relayRaceSerch(ch, pos, cna_data, start):
    current = start
    while ch != cna_data[current][0]:
        if current >= len(cna_data) - 1:
            return [-1, current]
        current += 1
    while pos >= cna_data[current][1]:
        if pos <= cna_data[current][2]:
            return [current, current]
        if current >= len(cna_data)-1:
            return [-1, current]
        if cna_data[current][0] != ch:
            return [-1, current]
        current += 1
    return [-1, current]

# Create list of samples and for each sample
def sampleNameCNAListCreater(f_rName):
    f_r = open(f_rName, 'r')
    line = f_r.readline()
    line = f_r.readline().strip().split()
    file_name = ''
    file_list = []
    while len(line) > 1:
        if file_name != line[0].strip('"'):
            file_name = line[0].strip('"')
            file_list.append(file_name)
        line = f_r.readline().strip().split()
    f_r.close()
    return file_list

# full dictionary by cna_list
def fullDictCNA(f_rName, samples_list):

    # Create dictionary of dictionary
    CNA_dict = {}
    for i in samples_list:
        CNA_dict[i] = {}

    f_r = open(f_rName, 'r')
    line = f_r.readline()
    line = f_r.readline().strip().split()
    sample_name = ''
    cna_list = []
    start = 0
    while len(line) > 1:

        # Check whether it is a new sample
        if sample_name != line[0].strip('"'):
            start = 0
            cna_list.sort(key=lambda x: [x[0], x[1], x[2]])
            sample_name = line[0].strip('"')

        # Adecvati check
        if line[5].strip('"') == '1' and line[6].strip('"') == '1':
            line = f_r.readline().strip().split()
            break
        # Save coordinates of segment
        cna = [int(line[2].strip('"').replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')), int(line[3].strip('"')), int(line[4].strip('"'))]

        # if cna_list is empti and it is ferst iteration there is problem to use "relayRaceSerch" function
        if len(cna_list) == 0:
            cna_list.append(cna)
            for i in CNA_dict.keys():
                CNA_dict[i][str(cna[0]) + ':' + str(cna[1]) + ':' + str(cna[2])] = ['1', '1']
            CNA_dict[sample_name][str(cna[0]) + ':' + str(cna[1]) + ':' + str(cna[2])] = [line[5].strip('"'), line[6].strip('"')]
            line = f_r.readline().strip().split()
            continue
        chek = relayRaceSerch(cna[0], cna[1], cna_list, start)
        start = chek[1]
        if chek[0] != -1:
            if cna_list[chek[1]][2] == cna[2]:
                CNA_dict[sample_name][str(cna[0]) + ':' + str(cna[1]) + ':' + str(cna[2])] = [line[5].strip('"'), line[6].strip('"')]
                line = f_r.readline().strip().split()
                continue
        cna_list.append(cna)
        for i in CNA_dict.keys():
            CNA_dict[i][str(cna[0]) + ':' + str(cna[1]) + ':' + str(cna[2])] = ['1', '1']
        CNA_dict[sample_name][str(cna[0]) + ':' + str(cna[1]) + ':' + str(cna[2])] = [line[5].strip('"'), line[6].strip('"')]
        line = f_r.readline().strip().split()
    f_r.close()
    return CNA_dict, cna_list

# Take dictionary with CNA and write WM and Wm file
def writeFileWsByDictCNA(file_dict, cna_list, WM = 'WMout.txt', Wm = 'Wmout.txt'):
    f_wWM = open(WM, 'w')
    f_wWm = open(Wm, 'w')
    f_wWM.write('CNA_mutation')
    f_wWm.write('CNA_mutation')
    for i in file_dict.keys():
        f_wWM.write('\t' + i)
        f_wWm.write('\t' + i)
    f_wWM.write('\n')
    f_wWm.write('\n')
    for i in cna_list:
        f_wWM.write(str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))
        f_wWm.write(str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))
        for j in file_dict.keys():
            f_wWM.write('\t' + file_dict[j][str(i[0]) + ':' + str(i[1]) + ':' + str(i[2])][0])
            f_wWm.write('\t' + file_dict[j][str(i[0]) + ':' + str(i[1]) + ':' + str(i[2])][1])
        f_wWM.write('\n')
        f_wWm.write('\n')
    f_wWM.close()
    f_wWm.close()

# Take file with SNV and return list of target columns and list of samples name
def sampleNameFinderSNV(f_rName, IDfilterc =  True):
    f_r = open(f_rName, 'r')
    heder = f_r.readline().strip().split()
    if IDfilterc:
        col_map = [0, 1, 4]
    else:
        col_map = [0, 1, 4]
    name_list = []
    for i in range(len(heder)):
        if heder[i].find('.AD') != -1:
            col_map.append(i)
            name_list.append(heder[i].replace('.AD', ''))
        if heder[i].find('.DP') != -1:
            col_map.append(i)
    f_r.close()
    return col_map, name_list

# Take file with SNV, list of column and samples names and return X and R matrixs
def fullMatrixXR(f_rName, col_map, IDfilter = True):
    f_r = open(f_rName, 'r')
    SNV_list = []
    R_matrix = []
    X_matrix = []
    line = f_r.readline()
    line = f_r.readline()
    while line != '':
        l = line.strip().split()
        if (l[col_map[0]]+l[col_map[1]]).find('NA') != -1:
            line = f_r.readline()
            continue
        if IDfilter and len(l[col_map[2]].strip('"')) <= 2:
            line = f_r.readline()
            continue
        R = []
        X = []
        SNV_list.append([int(l[col_map[0]].strip('"').replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')), int(l[col_map[1]].strip('"'))])
        R.append(str(SNV_list[-1][0]) + ':' + str(SNV_list[-1][1]))
        X.append(str(SNV_list[-1][0]) + ':' + str(SNV_list[-1][1]))
        if IDfilter:
            i = 3
        else:
            i = 2
        isNA = False
        while i < len(col_map):
            if len(l[col_map[i]].split(',')) <= 1:
                isNA = True
                break
            R.append(int(l[col_map[i]].split(',')[1].strip('"')))
            X.append(int(l[col_map[i + 1]].strip('"')))
            i += 2
            if R[-1] == X[-1]:
                isNA = True
                break
        if isNA:
            SNV_list.pop(-1)
            line = f_r.readline()
            continue
        R_matrix.append(R)
        X_matrix.append(X)
        line = f_r.readline()
    f_r.close()
    return R_matrix, X_matrix, SNV_list

# Take X and R matrixs SNV and write X and R file
def writeFileXRByMatrix(X_matrix, R_matrix, name_list, X = 'Xout.txt', R = 'Rout.txt'):
    f_wX = open(X, 'w')
    f_wR = open(R, 'w')
    f_wX.write('mut')
    f_wR.write('mut')
    for i in name_list:
        f_wX.write('\t' + i)
        f_wR.write('\t' + i)
    f_wX.write('\n')
    f_wR.write('\n')
    for i in X_matrix:
        line = ''
        for j in i:
            line += str(j)
            line += '\t'
        f_wX.write(line.strip() + '\n')
    f_wX.close()
    for i in R_matrix:
        line = ''
        for j in i:
            line += str(j)
            line += '\t'
        f_wR.write(line.strip() + '\n')
    f_wR.close()

# Create a file with CNA data with samples which presented in sample name list (exaple list of samples from SNV data)
# Return the name of new CNA data file
def sampleFilerToCNAFile(f_rName, sample_names, f_wName = 'CNA_for_target_patient.csv'):
    f_r = open(f_rName, 'r')
    f_w = open(f_wName, 'w')
    line = f_r.readline()
    f_w.write(line)
    for line in f_r:
        l = line.strip().split()
        for i in sample_names:
            if l[1].strip().strip('"').find(i.strip().strip('"')) != -1:
                f_w.write(line)
                break
    f_r.close()
    f_w.close()
    return f_wName

# create Y file with ovelopping of SNV and CNA segments by SNV and CNA lists
def createYFile(SNV_list, CNA_regions, f_wName = 'Yout.txt'):
    f_w = open(f_wName, 'w')
    f_w.write('mut\tnon-cna_region')
    for i in CNA_regions:
        f_w.write('\t' + str(i[0]) + ':' + str(i[1])+ ':' + str(i[2]))
    f_w.write('\n')
    for i in SNV_list:
        f_w.write(str(i[0]) + ':' + str(i[1]))
        line = ''
        non_CNA = True
        for j in CNA_regions:
            if i[0] == j[0] and i[1] >= j[1] and i[1] <= j[2]:
                line += '\t1'
                non_CNA = False
            else:
                line += '\t0'
        if non_CNA:
            f_w.write('\t' + '1')
        else:
            f_w.write('\t' + '0')
        f_w.write(line)
        f_w.write('\n')
    f_w.close()

# create C file with ovelopping of CNA segments by CNA lists
def createCFile(CNA_list, f_wName = 'Cout.txt'):
    f_w = open(f_wName, 'w')
    f_w.write('CNAs')
    for i in CNA_list:
        f_w.write('\t' + str(i[0]) + ':' + str(i[1])+ ':' + str(i[2]))
    f_w.write('\n')
    for i in CNA_list:
        f_w.write(str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))
        for j in CNA_list:
            if i[0] != j[0] or (i[1] > j[2] or i[2] < j[1]):
                f_w.write('\t' + '0')
            else:
                f_w.write('\t' + '1')
        f_w.write('\n')
    f_w.close()

# create C file with ovelopping of CNA segments by CNA regions
def createCByRegionsFile(CNA_list, CNA_regions, overlap_dict, f_wName = 'Cout.txt'):
    f_w = open(f_wName, 'w')
    f_w.write('CNAs')
    for i in CNA_list:
        f_w.write('\t' + str(i[0]) + ':' + str(i[1])+ ':' + str(i[2]))
    f_w.write('\n')
    for i in CNA_regions:
        f_w.write(str(i[1]) + ':' + str(i[2])+ ':' + str(i[3]))
        for j in CNA_list:
            if j in overlap_dict[i[0]]:
                f_w.write('\t' + '1')
            else:
                f_w.write('\t' + '0')
        f_w.write('\n')
    f_w.close()

# Redesigne CNA dictionary and CNA_list to merge similar segments
def redesigneCNADict(CNA_dict, CNA_list, step = 1000):
    association_dict = {}
    merge_CNA = [[-1,0,0,0,0,0]]
    number_code = 0
    for i in CNA_list:
        find_cluster = False
        for j in merge_CNA:
            if i[0] == j[1] and ((j[2] - step) <= i[1] <= (j[3] + step)) and ((j[4] - step) <= i[2] <= (j[5] + step)):
                find_cluster = True
                association_dict[j[0]].append(i)
                j[2] = min(j[2], i[1])
                j[3] = max(j[3], i[1])
                j[4] = min(j[4], i[2])
                j[5] = max(j[5], i[2])
                break
        if find_cluster:
            continue
        merge_CNA.append([number_code, i[0], i[1], i[1], i[2], i[2]])
        association_dict[number_code] = [i]
        number_code += 1
    merge_CNA.pop(0)
    new_CNA_list = []
    for  i in merge_CNA:
        new_CNA_list.append([i[1], i[2], i[5]])
    new_CNA_dict = {}
    for i in CNA_dict.keys():
        new_CNA_dict[i] = {}
        for j in merge_CNA:
            minor = '1'
            major = '1'
            for k in association_dict[j[0]]:
                key_CNA = str(k[0]) + ':' + str(k[1]) + ':' + str(k[2])
                if CNA_dict[i][key_CNA][0] != '1':
                    major = CNA_dict[i][key_CNA][0]
                if CNA_dict[i][key_CNA][1] != '1':
                    minor = CNA_dict[i][key_CNA][1]
            new_CNA_dict[i][str(j[1]) + ':' + str(j[2]) + ':' + str(j[5])] = [major, minor]
    print('CNA data was redesign from ' + str(len(CNA_list)) +' CNAs to ' + str(len(merge_CNA)) + ' CNAs')
    return new_CNA_dict, new_CNA_list

# clastering  overlapped CNAs in CNA regions by CNA list
def finderCNAregioonsOverlaping(CNA_list):
    CNA_regions = [[-1, 0, 0, 0]]
    overlap_dict = {}
    number_code = 0
    for i in CNA_list:
        find_overlap = False
        for j in CNA_regions:
            if i[0] == j[1] and (j[2] <= i[1] <= j[3] or j[2] <= i[2] <= j[3]):
                find_overlap = True
                overlap_dict[j[0]].append(i)
                j[2] = min(j[2], i[1])
                j[3] = max(j[3], i[2])
                break
        if find_overlap:
            continue
        CNA_regions.append([number_code, i[0], i[1], i[2]])
        overlap_dict[number_code] = [i]
        number_code += 1
    CNA_regions.pop(0)
    return overlap_dict, CNA_regions



