import os
import pandas as pd
from Bio import SeqIO
import subprocess
import torch
import torch.nn as nn
import time
import argparse
from token_esm import get_token
from binary_esm import binary_test
from cluster_esm import cluster_start
import model

start = time.perf_counter()

#default parameter
parser = argparse.ArgumentParser(description="PIDE")
parser.add_argument('input', type=str, help="Path of the input fasta file")
parser.add_argument('model', type=str, help="Path of the model parameter file")
parser.add_argument('-o', '--output', type=str, default='results', help="Path of the output directory")
parser.add_argument('-g', '--GPU', type=str, default='', help="Determine which GPU(s) to use, see README for more information")
parser.add_argument('-b', '--BatchSize', type=int, help='Define the batch size used in the prediction', default=2)
parser.add_argument('-n', '--MinPNum', type=int, default=5, help="The min prophage ORF number of a PI")
parser.add_argument('-d', '--Distance', type=int, default=3000, help="The clustering distance to use")
parser.add_argument('-s', '--PIScore', type=float, default=0.7, help="The threshold of PI score")
parser.add_argument('-m', '--meta', action='store_true', help="Use meta mode during ORF prediction")
args = parser.parse_args()

output = args.output
if output[-1] == '/':
    output = output[:-1]

if os.path.exists(output) is False:
    subprocess.run(f"mkdir {output}", shell=True)

#ORF ratio threshold
thresholds = {
    500: 1.0,
    1000: 1.0,
    3000: 0.5714285714285714,
    5000: 0.5555555555555556,
    10000: 0.5333333333333333,
    20000: 0.51724
}

# Preprocessing of contigs
print('Process 1: Split different contigs into separate files')
num = 0
name_length = {}
for seq_record in SeqIO.parse(args.input, "fasta"):
    name_length[num] = (seq_record.id, len(seq_record.seq))
    with open(f"{output}/{num}.fa", "w") as output_handle:
        SeqIO.write(seq_record, output_handle, "fasta")
    num += 1
print('Process 1 finish')

print('Process 2: Predict ORFs from contigs')
for i in range(num):
    if args.meta:
        subprocess.run(f"prodigal -i {output}/{i}.fa -a {output}/{i}.faa -o {output}/{i}.gbk -p meta -q", shell=True)
    else:
        subprocess.run(f"prodigal -i {output}/{i}.fa -a {output}/{i}.faa -o {output}/{i}.gbk -q", shell=True)
print('Process 2 finish')   

print('Process 3: Load the AI model')
#load the model and determine which GPU will be used
if args.GPU == '':
    device = torch.device('cpu')
else:
    os.environ['CUDA_VISIBLE_DEVICES'] = args.GPU
    device=torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
esmorf = model.ESMORF_binary(2)
esmorf = esmorf.to(device)
esmorf.load_state_dict(torch.load(args.model, map_location=device)) # 模型的parameter存放的位置
if args.GPU != '':
    if torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        esmorf = nn.DataParallel(esmorf)
print('Process 3 finish')   

print('Process 4: Predict prophage islands and virus contigs')
for i in range(num):
    df = get_token(f"{output}/{i}.faa")
    binary_test(df, esmorf, device, f'{output}/{i}.Predict.csv', f'{output}/{i}.Max.csv',args.BatchSize)
    subprocess.run(f"grep '>' {output}/{i}.faa | awk -F '#' -v OFS=',' '{{print $1, $2, $3}}' | sed 's/>//g' | sed 's/ //g' > {output}/{i}_tmp0",shell=True)
    subprocess.run(f"awk -F ',' -v OFS=',' 'NR > 1 {{print 1-$2,$2 }}' {output}/{i}.Predict.csv > {output}/{i}_tmp1",shell=True)
    subprocess.run(f"paste -d ',' {output}/{i}_tmp0 {output}/{i}_tmp1|sed '1i ORF,Start,End,B,P' > {output}/{i}.bed",shell=True)
    cluster_start(f'{output}/{i}.bed', f'{output}/{i}.cluster', args.MinPNum, args.Distance, args.PIScore)

    #phage contig
    with open(f'{output}/virus_contig.txt', 'w') as file:
        file.write('Contig,length,proportion\n')
    tmp = pd.read_csv(f'{output}/{i}.Max.csv')
    proportion = tmp.iloc[:, 1].value_counts(normalize=True).get(1.0, 0)
    contig_name = name_length[i][0]
    length_value = name_length[i][1]
    if length_value <= 1000 and proportion == 1:
        condition_met = True
    elif 1000 < length_value <= 3000 and proportion > 0.57:
        condition_met = True
    elif 3000 < length_value <= 5000 and proportion > 0.55:
        condition_met = True
    elif 5000 < length_value <= 10000 and proportion > 0.53:
        condition_met = True
    elif length_value > 10000 and proportion > 0.52:
        condition_met = True
    else:
        condition_met = False
    if condition_met:
        with open(f'{output}/virus_contig.csv', 'a') as file:
            file.write(f'{contig_name},{length_value},{proportion}\n')
print('Process 4 finish')

with open(f'{output}/cluster.csv', 'w') as fo:
    fo.write('Contig,Start,End,Score,B\n')
    for i in range(num):
        if os.path.exists(f'{output}/{i}.cluster'):
            with open(f'{output}/{i}.cluster') as f:
                f.readline()
                while 1:
                    line = f.readline()
                    if line == '':
                        break
                    fo.write(line)

subprocess.run(f"awk -F',' -v OFS=' ' '{{print $1, $2-1, $3}}' {output}/cluster.csv > {output}/location.bed",shell=True)
subprocess.run(f"seqtk subseq {args.input} {output}/location.bed > {output}/prophage.fasta",shell=True)
subprocess.run(f"prodigal -i {output}/prophage.fasta -a {output}/prophage.faa -p meta -q",shell=True)

end =  time.perf_counter()
run_time = round(end-start)
hour = run_time // 3600
minute = (run_time - 3600 * hour) // 60
second = run_time - 3600 * hour - 60 * minute

for i in range(num):
    subprocess.run(f'rm -rf {output}/{i}* {output}/location.bed', shell=True)
print(f'All finish, use time: {hour}h{minute}m{second}s')   