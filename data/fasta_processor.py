import esm
import random
import pandas as pd
from distfit import distfit
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class FastaProcessor():
    def __init__(self, path=''):
        if not path:
            return # empty spliter
        self.data_group = self.read_group(path)
        self.data_seq = [item for group in self.data_group.values() for item in group]
        self.length = len(self.data_seq)
        self.summarize()
        self.filter(min_samples=10)
        self.sample()
    
    def distribution(self, data):
        plt.figure(figsize=(3,2))
        sns.distplot(data)
        plt.legend()
        plt.grid(linestyle='--')
        plt.savefig("/zhouyuyang/Enzyme-multiclassification_lightning/figure/distribution.png")

    def read_group(self, path):
        group = {}
        with open(path, "r") as fp:
            all_content = fp.read()
            # drop the first since the first is empty
            for sequence in all_content.split('>')[1:]:
                # there's no \n in the middle of sentence, indicating it's the start
                if not "\n" in sequence.strip():
                    group[sequence.strip()] = []
                    group_name = sequence.strip()
                else:
                    # save name and sequence
                    name = sequence.strip().split("\n")[0]
                    label = sequence.strip().split("\n")[0].split('(')[-1].split(')')[0]
                    label = label.lower()
                    label = label.split('/')[0]
                    group[group_name].append((name, "".join(sequence.strip().split('\n')[1:]), label))
        self.distribution([len(part) for part in group.values()])
        return group

    def filter(self, min_samples):
        whole_data = [{"label":item[2], "seq":item[1], "tag":item[0], "group_name":name, "group_size":len(group)} for name, group in self.data_group.items() for item in group]
        whole_data = sorted(whole_data, key=lambda x:(x["label"], x["group_size"]))
        whole_data = pd.DataFrame(whole_data)
        self.whole_data = whole_data
        whole_data.to_csv("whole.txt", sep='\t', index=None)
        
        # filter
        filtered_data = self.whole_data[self.whole_data['label'].isin(self.all_label_cnts[self.all_label_cnts>=min_samples].index)]
        self.filtered_data = filtered_data
        filtered_data.to_csv("filtered.txt", sep='\t', index=None)
        
        # rebuild group data
        self.data_group = filtered_data.groupby('group_name').apply(lambda x: list(zip(x['tag'], x['seq'], x['label']))).to_dict()
        
        # rebuild sequence data
        self.data_seq = [item for group in self.data_group.values() for item in group]

    def sample(self, rate=0.25):
        def sample_func(all_data):
            min_data = []
            other_data = []
            names = []
            cur_label = ''
            left_df = []
            for row in all_data:
                if row[0] != cur_label:
                    cur_label=row[0]
                    group_name = row[3]
                    if not group_name in names:
                        names.append(group_name)
                        min_data+=[item for item in self.data_group[group_name]]
                else:
                    group_name = row[3]
                    if not group_name in names:
                        names.append(group_name)
                        other_data+=[item for item in self.data_group[group_name]]
                        left_df+=[data for data in all_data for item in self.data_group[group_name] if item[1]==data[1]]
            return min_data, other_data, left_df
        all_data=[tuple(x) for x in self.filtered_data.values]
        train_data=all_data
        val_data=[]
        while(len(val_data) < len(train_data)*rate):
            added_val, train_data, all_data = sample_func(all_data)
            val_data += added_val
        pd.DataFrame(train_data).to_csv("train.txt", sep='\t', index=None)
        pd.DataFrame(val_data).to_csv("val.txt", sep='\t', index=None)
        
        
    def summarize(self):
        self.all_label_cnts = pd.Series([data[2] for data in self.data_seq]).value_counts()
        pd.DataFrame(self.all_label_cnts).to_csv("summary.txt", sep='\t')

processor = FastaProcessor("/zhouyuyang/Enzyme-multiclassification_lightning/sequence/3644CoreRegion_group.fasta")