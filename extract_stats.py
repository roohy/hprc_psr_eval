import sys
import numpy as np
import pickle
import os.path

def vcf_single_stat(vcf_addr): #this function is used to calculate the total number of snvs carried by oos samples
    snp_list = set()
    with open(vcf_addr,'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#CHROM'):
                ids = line.strip().split()[9:]
                break
        sample_count = len(ids)
        id_mapper = {}
        for index,id in enumerate(ids):
            id_mapper[id] = index
        
        info = np.zeros((sample_count))
        for line in vcf_file:
            data = line.strip().split()
            sample_data = data[9:]
            # pos = int(data[1])
            # present_flag = pos in present_snps
            for index,gd in enumerate(sample_data):
                if gd[0] == '1' or gd[2] == '1':
                    info[index] += 1

                    
    return info,ids


def vcf_stat_dict( vcf_addr): #this one also creates reference sets to calculate the effects of adding new samples
    # present_snps = set()
    snp_dict = {}
    print(f'starting {vcf_addr}')
    # 
    with open(vcf_addr,'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#CHROM'):
                ids = line.strip().split()[9:]
                break
        sample_count = len(ids)
        id_mapper = {}
        for index,id in enumerate(ids):
            id_mapper[id] = index
        
        info = np.zeros((sample_count,2))
        for line in vcf_file:
            data = line.strip().split()
            sample_data = data[9:]
            pos = data[2]
            snp_dict[pos] = []
            # present_flag = pos in present_snps
            for index,gd in enumerate(sample_data):
                if gd[0] == '1' or gd[2] == '1':
                    snp_dict[pos].append(index)
                    info[index,1] += 1
                    info[index,0] += 1
                    
    return info,ids,snp_dict

def read_snp_list(file_name):
    snp_list = set()
    with open(file_name) as snplist_file:
        for line in snplist_file:
            snp_list.add(line.strip())
    return snp_list
class Strat:
    def __init__(self,label,file_name,chr_num):
        self.label = label
        self.file_name = file_name
        self.ids = []
        self.chr_num = chr_num
    def read_ids(self):
        with open(self.file_name) as id_file:
            for line in id_file:
                self.ids.append(line.strip())
    def set_snp_list(self,snp_list):
        self.snp_list = snp_list
    def check_newones(self):
        self.head = self.snp_list.copy() #make a deep copy
        self.newones = []
        self.iterativenewones = []
        for index,id in enumerate(self.ids):
            
            new_snp_count = 0#self.snp_list-self.
            iterative_new_snp_count = 0
            file_name = f'./inds/allchr/{id}/chr{self.chr_num}.snplist'
            if not os.path.isfile(file_name):
                continue
            with open(file_name,'r') as snplistfile:
                
                for line in snplistfile:
                    snp = line.strip()
                    if snp in self.snp_list:
                        new_snp_count += 1
                    if snp in self.head:
                        iterative_new_snp_count += 1 
                        self.head.remove(snp)
            self.newones.append(new_snp_count)
            self.iterativenewones.append(iterative_new_snp_count)
        self.newones = np.array(self.newones)
        self.iterativenewones = np.array(self.iterativenewones)
        self.argsorted = np.argsort(self.newones)
    def check_overlaps(self,info,snp_dict):
        improvements = []
        removed_snp_set = set()
        titem = info[:,1].copy()
        for index,id in enumerate(self.ids):
            file_name = f'./inds/allchr/{id}/chr{self.chr_num}.snplist'
            if not os.path.isfile(file_name):
                continue
            with open(file_name,'r') as snplistfile:
                for line in snplistfile:
                    snp = line.strip()
                    if snp in snp_dict and snp not in removed_snp_set:
#                         print('found new one')
                        removed_snp_set.add(snp)
                        for item in snp_dict[snp]:
                            titem[item] -= 1
            improvements.append(titem.copy())
        self.improvements = np.array(improvements)
        return self.improvements






if __name__ == '__main__':
    chr_num = sys.argv[1]
    vcf_addr = ''
    print(f'reading whole vcf {chr_num}')
    total_snvs,tsids = vcf_single_stat(vcf_addr=f'./oos_all/chr{chr_num}.vcf') #PLINK VCF file containing all common SNVs in the evaluation set
    print(f'reading small vcf {chr_num}')
    removed_sns,rsnvids,snp_dict = vcf_stat_dict(vcf_addr=f'./oos_phase1_removed/chr{chr_num}.vcf') #Plink VCF files containting all common SNVs not covered by Phase 1 Samples
    
    file_sets = [Strat('AFR 32','./id_sets/a32.txt',chr_num),Strat('EUR12/AFR24','./id_sets/a24e12.txt',chr_num),
             Strat('EUR24/AFR12','./id_sets/a12e24.txt',chr_num),
             Strat('EUR32','./id_sets/e32.txt',chr_num),Strat('Greedy Iterative','./id_sets/gi32.txt',chr_num)] #One strategy instance for every selection strategy
    print(f'making snp list {chr_num}')
    snp_list = set(snp_dict.keys())
    for strat in file_sets:
        print(f'start {strat.label} {chr_num}')
        strat.read_ids()
        strat.set_snp_list(snp_list)

        strat.check_newones()
        strat.check_overlaps(removed_sns,snp_dict)
    res = [file_sets,total_snvs,removed_sns,snp_dict,snp_list]
    pickle.dump(res,open(f'./basic_res/res{chr_num}.pkl','wb'))
