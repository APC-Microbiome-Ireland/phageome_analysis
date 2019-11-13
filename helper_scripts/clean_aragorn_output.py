fin = '../data/megaphage_trna.tab'
fout = '../data/megaphage_trna_clean.tab'

wr = open(fout, 'w')
with open(fin, 'r') as infile:
    for line in infile:
        if line[0] == '>':
            id = line.split('>')[1].split('\n')[0]
        elif line.find('genes found') == -1:
            line = line.split('\n')[0] + '\t' + id
            wr.write(line)
            wr.write('\n')
wr.close()
