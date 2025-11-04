#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
该脚本用于处理 UTR5 和 UTR3 距离测量数据，生成基因、转录本、细胞等不同层次的 UTR 信息矩阵。

输入：
    - 解析 UTR5 和 UTR3 测量结果，计算不同基因、转录本、细胞的 UTR 长度与相对值。
输出：
    - 根据基因、转录本和细胞生成多个结果文件，包含 UTR5 和 UTR3 信息。
"""

import sys
import os
import numpy as np
import argparse


def load_utr_data(file_path, utr_type):
    """
    解析 UTR 文件，加载数据到相应字典中。
    文件中的数据被存储在以下结构中：
        read_dict：存储每个 read 的 UTR 信息
        gene_dic：存储每个基因在不同细胞中的 UTR 信息
        trans_dic：存储每个转录本在不同细胞中的 UTR 信息
        cell_rel_dic：存储每个细胞的 UTR 信息
        sample_cell_set：存储所有样本细胞的集合
    """
    read_dict = {}
    gene_dic = {}
    trans_dic = {}
    cell_rel_dic = {}
    sample_cell_set = set()

    with open(file_path) as f:
        f.readline()  # 跳过表头
        for line in f:
            line = line.rstrip().split('\t')
            if len(line) < 12:
                print(f"Skipping line due to insufficient columns: {line}")
                continue
            gene_id = line[2]
            transcript_id = line[3]
            sample_cell = line[4]
            read_name = line[5]
            utr_rel = 1 - float(line[6]) if utr_type == 'utr5' else float(line[6]) - 2
            utr_length = 0 - int(line[8]) if utr_type == 'utr5' else int(line[11])

            index = f'{read_name}\t{sample_cell}\t{transcript_id}\t{gene_id}'
            read_dict.setdefault(index, {})
            read_dict[index][utr_type] = [utr_length, utr_rel]

            # 更新 gene_dic
            gene_dic.setdefault(gene_id, {})
            gene_dic[gene_id].setdefault(sample_cell, [[], []])
            gene_dic[gene_id][sample_cell][0].append(utr_length)
            gene_dic[gene_id][sample_cell][1].append(utr_rel)

            # 更新 trans_dic
            trans_dic.setdefault(transcript_id, {})
            trans_dic[transcript_id].setdefault(sample_cell, [[], []])
            trans_dic[transcript_id][sample_cell][0].append(utr_length)
            trans_dic[transcript_id][sample_cell][1].append(utr_rel)

            # 更新 cell_rel_dic
            cell_rel_dic.setdefault(sample_cell, [])
            cell_rel_dic[sample_cell].append(utr_rel)

            # 收集样本细胞信息
            sample_cell_set.add(sample_cell)

    return read_dict, gene_dic, trans_dic, cell_rel_dic, sample_cell_set


def write_matrix_file(matrix_file, data_dic, sample_cell_set, data_type):
    """
    写入基因或转录本的 UTR 长度和相对值矩阵文件
    """
    with open(matrix_file, "w") as f:
        outline = [data_type] + list(sample_cell_set)
        f.write('\t'.join(outline) + '\n')
        for data_id in data_dic:
            outline_data = [data_id]
            for sample_cell in sample_cell_set:
                mean_utr_length = np.mean(data_dic[data_id].get(sample_cell, [[], []])[0]) if sample_cell in data_dic[data_id] else '-'
                mean_utr_rel = np.mean(data_dic[data_id].get(sample_cell, [[], []])[1]) if sample_cell in data_dic[data_id] else '-'
                outline_data.append(str(mean_utr_length) + "," + str(mean_utr_rel))
            f.write("\t".join(outline_data) + "\n")


def write_cell_utr_file(cell_rel_dic, output_dir, utr_type):
    """
    写入细胞级别的 UTR 相对值文件
    """
    cell_rel_file = f'{output_dir}/cell_{utr_type}_relative.txt'
    with open(cell_rel_file, "w") as f:
        for cell_id in cell_rel_dic:
            mean_utr_length = np.mean(cell_rel_dic[cell_id])
            counts_utr = len(cell_rel_dic[cell_id])
            f.write(cell_id + "\t" + str(mean_utr_length) + "\t" + str(counts_utr) + "\n")


def write_read_utr_file(read_dict, output_dir):
    """
    写入 read 的 UTR5 和 UTR3 信息
    """
    reads_utr_file = f'{output_dir}/reads_utr5_utr3_length.txt'
    with open(reads_utr_file, 'w') as f:
        for index in read_dict:
            outline = [index]
            outline += read_dict[index].get("utr5", ['-','-'])
            outline += read_dict[index].get("utr3", ['-','-'])
            f.write('\t'.join(map(str, outline)) + '\n')


def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="Process UTR5 and UTR3 distance measures and generate matrices.")
    parser.add_argument('-utr5', '--utr5_file', required=True, help="UTR5 distance measures file")
    parser.add_argument('-utr3', '--utr3_file', required=True, help="UTR3 distance measures file")
    parser.add_argument('-o', '--output_dir', default="utr_results", help="Output directory (default: utr_results)")

    args = parser.parse_args()

    # 创建输出目录
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 处理 UTR5 文件
    print("[INFO] Processing UTR5 data...")
    read_dict, gene_dic, trans_dic, cell_rel_dic, sample_cell_set = load_utr_data(args.utr5_file, 'utr5')

    # 写入文件
    write_matrix_file(f'{args.output_dir}/gene_utr5_matrix.txt', gene_dic, sample_cell_set, "gene_id")
    write_matrix_file(f'{args.output_dir}/transcript_utr5_matrix.txt', trans_dic, sample_cell_set, "transcript_id")
    write_cell_utr_file(cell_rel_dic, args.output_dir, 'utr5')
    write_read_utr_file(read_dict, args.output_dir)

    # 清理字典为下一次处理
    gene_dic = {}
    trans_dic = {}
    cell_rel_dic = {}
    sample_cell_set = set()

    # 处理 UTR3 文件
    print("[INFO] Processing UTR3 data...")
    read_dict, gene_dic, trans_dic, cell_rel_dic, sample_cell_set = load_utr_data(args.utr3_file, 'utr3')

    # 写入文件
    write_matrix_file(f'{args.output_dir}/gene_utr3_matrix.txt', gene_dic, sample_cell_set, "gene_id")
    write_matrix_file(f'{args.output_dir}/transcript_utr3_matrix.txt', trans_dic, sample_cell_set, "transcript_id")
    write_cell_utr_file(cell_rel_dic, args.output_dir, 'utr3')
    write_read_utr_file(read_dict, args.output_dir)

    print("[INFO] Processing complete.")


if __name__ == '__main__':
    main()
