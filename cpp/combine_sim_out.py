import sys
import os

ids = [str(x) for x in range(1,11)]

prefix = 'out/sim_output_summary_'
suffix = '.txt'


h2 = ['0.001', '0.005', '0.01', '0.05', '0.1']
fin_size = ['1000', '5000', '10000']#, '20000']

for h in h2:
    for size in fin_size:
        current_prefix = prefix + size + '_' + size + '_1_' + h + '_'
        f_out = open(current_prefix + 'all' + suffix, 'w')
        try:
            f_out.write(open(current_prefix + ids[0] + suffix, 'r').readline())    
            for i in ids:
                f = open(current_prefix + i + suffix, 'r')
                f.readline() # header
                line = f.readline()
                while line != '':
                    f_out.write(line)
                    line = f.readline()
                f.close()
            f_out.close()
            print current_prefix
        except:
            print 'not_found :', current_prefix

# combine LOD files

# ,chr,chr,chr
# ,ps,ps,ps
# sim,score,score,score
# sim,score,score,score

# to:
# chr,ps,0,1,2,3
# chr,ps,LOD,LOD,LOD
# chr,ps,LOD,LOD,LOD

prefix = 'out/sim_output_LOD_'
nsims = 100

for h in h2:
    for size in fin_size:
        print h, size
        current_prefix = prefix + size + '_' + size + '_1_' + h + '_'
        f_out = open(current_prefix + 'all' + suffix, 'w')

        f_out.write('chr\tps')
        for i in range(nsims):
            f_out.write('\tsim' + str(i))
        f_out.write('\n')

        f = open(current_prefix + ids[0] + suffix, 'r')
        chrm = [str(int(x) + 1) for x in f.readline().strip().split(',')[1:]]
        #chrm = f.readline().strip().split(',')[1:]
        ps = f.readline().strip().split(',')[1:]
        f.close()

        lines = [[] for x in range(len(ps))]
             
        for i in ids:
            f = open(current_prefix + i + suffix, 'r')
            f.readline() #chr
            f.readline() #ps
            line = f.readline()
            c = 0
            while line != '':
                scores = line.strip().split(',')[1:]
                for loc in range(len(scores)):
                    lines[loc].append(scores[loc])
                line = f.readline()
            f.close()

        for j in range(len(lines)):
            f_out.write(chrm[j] + '\t' + ps[j])
            for l in lines[j]:
                f_out.write('\t' + l)
            f_out.write('\n')

        f_out.close()



"""
h2 = ['0.001', '0.005', '0.01', '0.05', '0.1', '0.2']
int_size = ['1', '10', '100', '1000', '5000', '10000', '20000']
fin_size = ['1000', '5000', '10000', '20000']

for h in h2:
    for isize in int_size:
        for size in fin_size:
            current_prefix = prefix + size + '_' + isize + '_1_' + h + '_'
            f_out = open(current_prefix + 'all' + suffix, 'w')
            try:
                if h == '0.05' and size == '10000' and isize != '100' and isize != '5000' and isize != '1' and isize != '10':
                    current_prefix = prefix + '10_' + size + '_' + isize + '_1_' + h + '_'
                    print '**********'
                f_out.write(open(current_prefix + ids[0] + suffix, 'r').readline())
                
                for i in ids:
                    f = open(current_prefix + i + suffix, 'r')
                    f.readline() # header
                    line = f.readline()
                    while line != '':
                        f_out.write(line)
                        line = f.readline()
                    f.close()
                
                f_out.close()
                print current_prefix
            except:
                print 'not_found :', current_prefix
                try:
                    os.remove(current_prefix + 'all' + suffix)
                except:
                    os.remove(prefix + size + '_' + isize + '_1_' + h + '_all' + suffix)
"""
