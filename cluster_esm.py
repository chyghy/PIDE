import pandas as pd
import re

def max_average_subarray(nums, MinPNum, PIScore):
    max_avg = sum(nums)/len(nums)
    start_index = 0
    end_index = len(nums)
    if max_avg > PIScore:
        return [start_index,end_index], max_avg
    for window_length in range(len(nums)-1,MinPNum-1,-1):
        for i in range(0,len(nums)-window_length+1):
            avg=sum(nums[i:i+window_length])/window_length
            if avg>max_avg:
                start_index=i
                end_index=i+window_length
                max_avg=avg
        if max_avg > PIScore:
            return [start_index,end_index], max_avg
    print('Warning: The sequence all combo less than '+ str(PIScore))
    return [start_index,end_index], max_avg
	
def add_b(row):
    seq = row['B']
    for key,value in info_dict.items():
        if (value['Start'] >= row['Start']) & (value['End'] <= row['End']):
            seq = seq +'['+key+':'+str(value['Start'])+'-'+str(value['End'])+']'
        else:
            seq=seq + ''
    return seq

def cluster_start(input_data,output_data,MinPNum,Distance,PIScore):	
    df_all = pd.read_table(input_data,sep=',')
    pat = re.compile(r'_\d+$')
    df_all['contig'] = [re.sub(pat,'',i) for i in df_all['ORF']]
    chy=True
    for contig in df_all['contig'].unique():
        df = df_all.copy().loc[df_all['contig'] == contig,]
        df['type'] = 'P'
        df.loc[df['B']>0.5,'type'] = 'B'
        seq = df['type'].str.cat()
        pattern_b = re.compile('B+')
        temp = re.finditer(pattern_b,seq)
        global info_dict
        info_dict = {}
        for i,j in enumerate(temp):
            info_dict['B'+str(i+1)] = {'Start':df['Start'].tolist()[j.start()],'End':df['End'].tolist()[j.end()-1],'length':df['End'].tolist()[j.end()-1] - df['Start'].tolist()[j.start()]  }      
        seq_new = re.sub(pattern_b,'B',seq) 
        Interval = Distance
        ind = 0
        for i in range(len(seq_new)):
            if seq_new[i] =='B' :
                ind+=1
                if info_dict['B'+str(ind)]['length'] < Interval :
                    seq_new=seq_new[0:i]+'p'+seq_new[i+1:]            
        pattern = re.compile('[Pp]+')   
        temp = re.findall(pattern,seq_new)
        start=0
        cluster_name = []
        start_list=[]
        end_list=[]
        score_list=[]
        total_orf_list = [] 
        b_count_list = []
        for i,j in enumerate(temp):
            if j.count('P') >=MinPNum:
                prob = df.loc[df['type']=='P','P'].tolist()[start:(start+j.count('P'))]
                new_idx, avg_max = max_average_subarray(prob,MinPNum,PIScore)
                if avg_max > PIScore:
                    start_idx = df[df['type']=='P'].index[start+new_idx[0]]
                    end_idx = df[df['type']=='P'].index[start+new_idx[1]-1]
                    cluster_df = df.loc[range(start_idx, end_idx + 1)]
                    b_count = (cluster_df['type'] == 'B').sum()
                    total_orf = len(cluster_df)
                    cluster_name.append('Cluster' + str(i+1)) 
                    start_list.append(df.loc[df['type']=='P','Start'].tolist()[start+new_idx[0]])
                    end_list.append(df.loc[df['type']=='P','End'].tolist()[start+new_idx[1]-1])
                    score_list.append(avg_max)
                    total_orf_list.append(total_orf)
                    b_count_list.append(b_count)
            start = j.count('P')+start

        if len(cluster_name) >0:
            Out = pd.DataFrame()
            Out["Cluster"]=cluster_name
            Out["Contig"] = contig
            Out['Start']=start_list
            Out['End'] = end_list
            Out['Score'] = [round(i,2) for i in score_list]
            Out['B'] = ''
            Out['B'] = Out.apply(add_b,axis=1)
            Out['Total_ORFs'] = total_orf_list
            Out['B_count'] = b_count_list
            Out['B_ratio'] = Out['B_count'] / Out['Total_ORFs']

            if chy:
                out_csv = Out
                chy=False
            else:
                out_csv = pd.concat([out_csv,Out])
    if not chy:     
        out_csv.drop('Cluster', axis=1).to_csv(output_data, index=False)
