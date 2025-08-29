#!/bin/bash

fasta_file="transcript_sequences.fasta"
output_file="mouse_MFE.csv"

echo "transcript_id,MFE" > $output_file

total=$(grep -c "^>" $fasta_file)  # 计算序列总数
count=0
seq=""  # 存储序列内容

# 逐行处理 fasta 文件
while IFS= read -r line; do
    if [[ $line == ">"* ]]; then  # 遇到新的 ID
        if [[ -n "$seq" ]]; then  # 如果之前有序列数据，计算 MFE
            mfe=$(echo -e "$header\n$seq" | RNAfold --noPS | awk '/)/ {print $NF}' | tr -d '>()')
            echo "$transcript_id,$mfe" >> $output_file  # 写入 CSV
            ((count++))
            echo -ne "Processing: $count / $total sequences\r"
        fi
        header="$line"  # 记录新 ID
        transcript_id=${line#>}  # 提取 ID
        seq=""  # 重置序列
    else
        seq+="$line\n"  # 累加序列
    fi
done < "$fasta_file"

# 处理最后一个序列
if [[ -n "$seq" ]]; then
    mfe=$(echo -e "$header\n$seq" | RNAfold --noPS | awk '/)/ {print $NF}' | tr -d '>()')
    echo "$transcript_id,$mfe" >> $output_file
    ((count++))
    echo -ne "Processing: $count / $total sequences\r"
fi

echo -e "\nProcessing complete. Results saved to $output_file."
