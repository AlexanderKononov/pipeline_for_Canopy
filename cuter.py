import sys
arg = sys.argv

'''
coint=1
name_sample = []
f_r = open(arg[2], 'r')
line = f_r.readline()
line = f_r.readline().strip().split()
name_sample.append(line[0].strip('"'))
while coint < 4:
    if line[0].strip('"').find(name_sample[-1]) != -1
        name_sample.append(line[0].strip('"'))
        coint += 1
f_r.close()
'''

f_r = open(arg[1], 'r')
f_w = open(arg[1].replace('.csv', '_cut.csv'), 'w')
line = f_r.readline()
f_w.write(line)
line = f_r.readline()
while line != '':
    l = line.strip('"').split()
    if 4 <= int(l[0].strip('"').replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25').replace('NA', '1')) <= 6:
        f_w.write(line)
    line = f_r.readline()
f_r.close()
f_w.close()

f_r = open(arg[2], 'r')
f_w = open(arg[2].replace('.csv', '_cut.csv'), 'w')
line = f_r.readline()
f_w.write(line)
line = f_r.readline()
while line != '':
    l = line.strip('"').split()
    if 4 <= int(l[2].strip('"').replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25').replace('NA', '1')) <= 6:
        f_w.write(line)
    line = f_r.readline()
f_r.close()
f_w.close()