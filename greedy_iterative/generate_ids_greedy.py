import numpy as np

import pandas as pd 

from multiprocessing import Pool
import pickle
combined_id_addr = '/sc/arion/projects/kennylab/roohy/hprc/1kg/phase1/combined_ids2.txt'
phase1_id_addr = '/sc/arion/projects/kennylab/roohy/hprc/1kg/phase1/phase1_ids2.txt'
oos_id_addr = '/sc/arion/projects/kennylab/roohy/hprc/1kg/phase1/oos_ids2.txt'
new_ids_addr='/sc/arion/projects/kennylab/roohy/hprc/1kg/greedy/v2/nids.txt'

def read_vcf(chr_num):
    present_snps = set()
    print(f'starting chr{chr_num}')
    with open(f'./temp/sample_chr{chr_num}.bim','r') as bim_file:
        for line in bim_file:
            present_snps.add(int(line.strip().split()[-3]))
    with open(f'./temp/oos_chr{chr_num}.vcf','r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#CHROM'):
                ids = line.strip().split()[9:]
                break
        sample_count = len(ids)
        info = np.zeros((sample_count,2))
        for line in vcf_file:
            data = line.strip().split()
            sample_data = data[9:]
            pos = int(data[1])
            present_flag = pos not in present_snps
            for index,gd in enumerate(sample_data):
                if gd[0] == '1' or gd[2] == '1':
                    info[index,0] += 1
                    if present_flag:
                        info[index,1] += 1
    return info,ids

def read_id_file(addr):
    id_list = []
    with open(addr,'rt') as id_file:
        for line in id_file:
            id_list.append(line.strip())
    print(f'Read {len(id_list)} items from the file:{addr}')
    id_set = set(id_list)
    print(f'It included {len(id_set)} unique ids')
    return id_set
def read_parents(addr,pids):
    parent_ids = []
    with open('./1kGP.3202_samples.pedigree_info.txt') as  pedigree_file:
        pedigree_file.readline()
        for line in pedigree_file:
            data = line.strip().split()
            if data[0] in pids:
                parent_ids.append(data[1])
                parent_ids.append(data[2])
    return set(parent_ids)

def write_id_file(id_set,file_addr):
    with open(file_addr,'w') as output_file:
        for item in id_set:
            print(f'0 {item}',file=output_file)


if __name__ == '__main__':
    #new_ids_addr = 'NA'
    # main_ids = read_id_file(main_id_addr)
    oos_ids = read_id_file(oos_id_addr)
    phase1_ids = read_id_file(phase1_id_addr)
    euro_12_id_addr = '/sc/arion/projects/kennylab/roohy/hprc/1kg/euro_basic/12/ids_1col.txt'
    euro12_ids = read_id_file(euro_12_id_addr)
    nids = read_id_file(new_ids_addr)
    temp_id_set = euro12_ids | phase1_ids

    with Pool(10) as p:
        info_list = p.map(read_vcf,np.arange(1,23))
    counts = np.array([item[0] for item in info_list]).sum(axis=0)
    max_ind = np.argmax(counts[:,1])
    new_id= info_list[0][1][max_ind]
    pickle.dump(info_list,open(f'info_list{len(nids)}.pkl','wb'))
    nids.add(new_id)
    temp_id_set = temp_id_set | nids
    with open(new_ids_addr,'w') as new_id_file:
        for nid in nids:
            print(nid,file=new_id_file)
    write_id_file(temp_id_set,'./temp/ids.txt')
    write_id_file(oos_ids-temp_id_set,'./temp/oos_ids.txt')

