import pandas as pd, sys, os

folder_path = sys.argv[1] #path to folder with reports
valid_mut_dict = dict()
valid_ann_dict = dict()
sample_count = 0
for file in os.listdir(folder_path):
    if 'ann.csv' in file:
        sample_count +=1
        if len(valid_mut_dict) == 0:
            df = pd.read_csv(f'{folder_path}/{file}').applymap(str)
            cur_mut_set = list(df.AMINO_ACID_CHANGE)
            cur_nt_list = list(df.MUTATION)
            valid_nt_dict = dict(zip(cur_mut_set, cur_nt_list))
            cur_gene_list = list(df.GENE)
            valid_gene_dict = dict(zip(cur_mut_set, cur_gene_list))
            # cur_ann_list = list(df.ANNOTATION)
            # valid_ann_dict = dict(zip(cur_mut_set, cur_ann_list))
            while 'nan' in cur_mut_set: cur_mut_set.remove('nan')
            valid_mut_dict = dict(zip(cur_mut_set, [1 for _ in range(len(cur_mut_set))]))
        else:
            df = pd.read_csv(f'{folder_path}/{file}').applymap(str)
            cur_mut_set = df.AMINO_ACID_CHANGE
            for mut in cur_mut_set:
                if mut in valid_mut_dict.keys():
                    valid_mut_dict[mut] += 1
                elif mut not in valid_mut_dict.keys() and mut != 'nan':
                    valid_mut_dict[mut] = 1
                # if mut not in valid_ann_dict.keys():
                #     condition = df['AMINO_ACID_CHANGE'] == mut
                #     valid_ann_dict[mut] = df.iloc[df.index[condition], 3].values[0]
                if mut not in valid_gene_dict.keys():
                    condition = df['AMINO_ACID_CHANGE'] == mut
                    valid_gene_dict[mut] = df.iloc[df.index[condition], 1].values[0]
                    # print(dir(df.iloc[df.index[condition], 1]))
                if mut not in valid_nt_dict.keys():
                    condition = df['AMINO_ACID_CHANGE'] == mut
                    valid_nt_dict[mut] = df.iloc[df.index[condition], 0].values[0]

df = pd.DataFrame(data={'Mutation':list(valid_mut_dict.keys()), "Number_of_samples":list(valid_mut_dict.values())})
# df['Annotation'] = df['Mutation'].map(valid_ann_dict)
df['Nt_change'] = df['Mutation'].map(valid_nt_dict)
df['Gene'] = df['Mutation'].map(valid_gene_dict)
df['%_of_samples'] = round(100*(df['Number_of_samples'] / sample_count),2)
df = df[['Nt_change', 'Mutation', 'Number_of_samples', '%_of_samples', 'Gene']]
df.to_csv(f'{folder_path}/mutation_statistics.csv', index=False)