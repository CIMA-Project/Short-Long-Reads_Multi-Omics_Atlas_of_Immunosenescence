#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多进程版注释信息筛选脚本
----------------------------------
功能概述：
  - 从 long_file 载入基因到转录本的映射信息；
  - 读取 annot_file，将每一行中的 reads 与参考转录本匹配；
  - 若匹配成功则输出；
  - 支持多进程并行处理以提升性能。
"""

import sys
import os
import multiprocessing as mp


def load_adic(long_file):
    """
    从 long_file 读取映射关系，生成字典 adic：
        key: gene_name
        value: [transcript_name, transcript_length]
    long_file 格式：
        transcript_name    gene_name    length
    """
    adic = {}
    with open(long_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                # fields[0]: 转录本名称；fields[1]: 基因名称；fields[2]: 转录本长度
                adic[fields[1]] = [fields[0], int(fields[2])]
    return adic


def process_lines(lines, adic):
    """
    处理 annot_file 的部分行（由多进程并行执行）
    核心逻辑：
        - 从每行解析 reads 的转录本信息；
        - 若 reads 转录本名缺失（"."），则尝试用 adic 中的映射补全；
        - 当 reads 转录本与参考转录本一致时，输出该行简化字段。
    """
    results = []
    split_semi = lambda s: s.split(';')
    split_bar = lambda s: s.split('|')

    for line in lines:
        line = line.strip()
        lineL = line.split('\t')
        if len(lineL) < 5:
            continue  # 跳过列数不足的行

        # 第4列信息中包含多个字段，用分号分隔
        reads_info_list = split_semi(lineL[3])
        reads_trans_name = reads_info_list[-3]  # 从末尾第3个字段提取 reads 对应转录本名

        # 若转录本名缺失，则尝试从 gene_name 列表推断
        if reads_trans_name == ".":
            gene_name_L = reads_info_list[1:-3]  # 基因名列表（不含前后附加字段）

            if len(gene_name_L) == 1:
                # 单基因情况，直接查表
                gene_name = gene_name_L[0]
                reads_trans_name = adic.get(gene_name, ["."])[0]
            else:
                # 多基因情况，选择转录本长度最长者
                best = max(
                    (adic[gene] for gene in gene_name_L if gene in adic),
                    key=lambda x: x[1],
                    default=["."]
                )
                reads_trans_name = best[0]

        # 从最后第4列提取参考转录本名称
        ref_trans_name = split_bar(lineL[-4])[0]

        # 比较 reads 转录本名与参考转录本名
        if reads_trans_name == ref_trans_name:
            # 输出部分字段：前6列 + 最后7列
            out_list = lineL[:6] + lineL[-7:]
            results.append("\t".join(out_list))

    return results


def chunkify(file_path, num_chunks):
    """
    将输入文件分割为若干块（用于多进程处理）
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    chunk_size = max(1, len(lines) // num_chunks)
    return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]


def main():
    """
    主函数：解析命令行参数，执行多进程数据处理
    使用方式：
        python script.py long_file annot_file output_file [num_processes]
    参数说明：
        long_file      - 含转录本与基因对应关系的文件
        annot_file     - 注释信息文件（待筛选）
        output_file    - 输出结果文件
        num_processes  - 可选，多进程数量（默认使用 CPU 核心数）
    """
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <long_file> <annot_file> <output_file> [num_processes]")
        sys.exit(1)

    long_file = sys.argv[1]
    annot_file = sys.argv[2]
    output_file = sys.argv[3]
    num_processes = int(sys.argv[4]) if len(sys.argv) > 4 else mp.cpu_count()

    print(f"[INFO] Loading mapping from {long_file} ...")
    adic = load_adic(long_file)
    print(f"[INFO] Loaded {len(adic)} gene-transcript mappings.")

    print(f"[INFO] Splitting {annot_file} into {num_processes} chunks ...")
    chunks = chunkify(annot_file, num_processes)

    print(f"[INFO] Starting {num_processes} parallel processes ...")
    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(process_lines, [(chunk, adic) for chunk in chunks])

    print(f"[INFO] Writing results to {output_file} ...")
    with open(output_file, 'w') as out:
        for sublist in results:
            if sublist:  # 跳过空结果
                out.write('\n'.join(sublist) + '\n')

    print(f"[INFO] Done! Output written to {output_file}")


if __name__ == '__main__':
    main()
