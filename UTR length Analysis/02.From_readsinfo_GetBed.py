#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
from pathlib import Path

def load_read_assignments(read_assignments_file):
    """加载 read_assignments 文件，返回 read_name → transcript_name 的映射字典"""
    adic = {}
    with gzip.open(read_assignments_file, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            lineL = line.strip().split("\t")
            if len(lineL) < 4:
                continue
            read_name = lineL[0]
            trans_name = lineL[3]
            adic[read_name] = trans_name
    print(f"Loaded {len(adic)} read assignments from {read_assignments_file}")
    return adic


def process_info_file(info_file, adic, output_dir):
    """解析单个 info 文件并生成对应的 TSS 和 TES BED 文件"""
    base_name = Path(info_file).name
    prefix = base_name[:-7] if base_name.endswith(".txt.gz") else Path(base_name).stem

    tss_filename = Path(output_dir) / f"{prefix}.TSS.bed"
    tes_filename = Path(output_dir) / f"{prefix}.TES.bed"

    with open(tss_filename, "w") as tss_file, open(tes_filename, "w") as tes_file:
        with gzip.open(info_file, "rt", encoding="utf-8") as f:
            for line in f:
                line_parts = line.strip().split(" ")
                if len(line_parts) < 13:
                    continue

                read_name = line_parts[0]
                chrom = line_parts[1]
                strand = line_parts[2]
                try:
                    tss_site = int(line_parts[3])
                    tes_site = int(line_parts[5])

                except ValueError:
                    continue
                tss_flag = "T" if line_parts[4] != "nan" else "F"
                tes_flag = "T" if line_parts[6] != "nan" else "F"

                sample_cell = line_parts[7]
                gene_name = line_parts[8]

                trans_name = adic.get(read_name, ".")

                tss_file.write(
                    f'{chrom}\t{tss_site}\t{tss_site + 1}\t{sample_cell};{gene_name};{trans_name};{read_name};{tss_flag}\t.\t{strand}\n'
                )
                tes_file.write(
                    f'{chrom}\t{tes_site}\t{tes_site + 1}\t{sample_cell};{gene_name};{trans_name};{read_name};{tes_flag}\t.\t{strand}\n'
                )

    print(f"Processed {info_file} → {tss_filename.name}, {tes_filename.name}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract TSS and TES information from info files and match reads to transcripts."
    )
    parser.add_argument(
        "-r", "--read-assignments",
        required=True,
        help="Path to the gzipped read assignments file (e.g., read_assignments.txt.gz)"
    )
    parser.add_argument(
        "-i", "--info-files",
        required=True,
        nargs="+",
        help="List of gzipped info files to process (e.g., sample1.txt.gz sample2.txt.gz)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=".",
        help="Output directory for generated BED files (default: current directory)"
    )

    args = parser.parse_args()

    # 检查输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    # 加载 read → transcript 映射
    adic = load_read_assignments(args.read_assignments)

    # 逐个处理 info 文件
    for info_file in args.info_files:
        process_info_file(info_file, adic, args.output_dir)

    print("All info files processed successfully.")


if __name__ == "__main__":
    main()
