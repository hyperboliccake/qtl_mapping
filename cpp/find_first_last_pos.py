# use this to find last position of a snp in the data file, or just
# use the sgd file instead
f = open('/net/akey/vol2/aclark4/joseph_data/subset_302_homozygotes.map', 'r')
f_out = open('first_last_pos_all_chr.txt', 'w')
line = f.readline()
prev_chr = '1'
prev_ps = '-1'
while line != '':
    line = line.split()
    current_chr = line[0]
    if current_chr != prev_chr:
        f_out.write(prev_chr + ',1,' + prev_ps + '\n')
        prev_chr = current_chr
    prev_ps = line[3]
    line = f.readline()
f_out.write(prev_chr + ',1,' + prev_ps + '\n')
f.close()
f_out.close()
