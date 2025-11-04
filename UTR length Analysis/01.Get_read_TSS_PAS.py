import pysam
from pyfaidx import Fasta
import regex
import multiprocessing
import argparse

def get_sequence(read, ref):
    chrom = read.reference_name
    strand = "-" if read.is_reverse else "+"
    NB_tag = read.get_tag('NB')
    GN_tag = read.get_tag('GN')

    # 确认 Full-Lengh-Reads
    # 计算 poly(A) site 位置
    # 提取 poly(A) site 上游是否为 A-rich 区域
    if strand == "+":
        APA_pos = read.reference_end  # 正链：poly(A) site是比对末端
        
        ref_upstream_start = max(APA_pos - 15, 0)
        ref_upstream_end = APA_pos
        ref_downstream_start = APA_pos
        ref_downstream_end = APA_pos + 15
        
        reads_downstream_start = read.query_alignment_end
        reads_downstream_end = reads_downstream_start + 20
        
        read_start_end = read.query_alignment_start
        read_start_start = max(read_start_end-10, 0)

        read_start_seq = read.query_sequence[read_start_start:read_start_end]
        read_downstream_seq = read.query_sequence[reads_downstream_start:reads_downstream_end]
        
        TSS_site = read.reference_start
        
    else:
        APA_pos = read.reference_start  # 负链：poly(A) site是比对起点

        ref_upstream_start = APA_pos
        ref_upstream_end = APA_pos + 15
        ref_downstream_start =  max(APA_pos - 15, 0)
        ref_downstream_end = APA_pos

        reads_downstream_end  = read.query_alignment_start
        reads_downstream_start = max(reads_downstream_end - 20, 0)
        
        read_start_start = read.query_alignment_end
        read_start_end = read_start_start + 10

        read_start_seq = read.query_sequence[read_start_start:read_start_end]
        read_start_seq = read_start_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
        read_downstream_seq = read.query_sequence[reads_downstream_start:reads_downstream_end]
        read_downstream_seq = read_downstream_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))    
        TSS_site = read.reference_end
        
    try:
        ref_upstream_seq = ref[chrom][ref_upstream_start:ref_upstream_end].seq
        ref_downstream_seq = ref[chrom][ref_downstream_start:ref_downstream_end].seq
    except KeyError:
        return  # 染色体不在参考中，跳过
        
    if strand == "-":
        ref_upstream_seq = ref_upstream_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
        ref_downstream_seq = ref_downstream_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))    

    return (read.query_name, chrom, APA_pos, strand, 
            ref_upstream_seq, ref_downstream_seq, read_downstream_seq, 
            read_start_seq, TSS_site,
            read.query_alignment_start, read.query_alignment_end, read.query_sequence)

def counts_base(seq):
    counts_dic = {}
    seq = seq.upper()
    for i in seq:
        counts_dic[i] = counts_dic.get(i,0) +1
    for base in counts_dic:
        counts = counts_dic[base]
        precentage = counts/len(seq)
        counts_dic[base] = [counts,precentage]
    return counts_dic

def fuzzy_match(seq, pattern, max_mismatch=2):
    """Return True if pattern is found in seq allowing max_mismatch edits"""
    fuzzy_pattern = f"({pattern}){{e<={max_mismatch}}}"
    return regex.search(fuzzy_pattern, seq) is not None
    
def filter_FalsePosivite_APA_reads(ref_upstream_seq, reads_tail, ployA_tail = "CTCTGCGTTG"):
    ref_upstream_seq = ref_upstream_seq.upper()
    reads_tail = reads_tail.upper()
    fl_flag = False
    basecounts_reads_tail_dic = counts_base(reads_tail)
    if fuzzy_match(reads_tail, ployA_tail, max_mismatch=2):
        fl_flag = True

    if fuzzy_match(reads_tail, "AAAAAAAAA", max_mismatch=2):
        fl_flag = True

    ref_flag = True
    if fuzzy_match(ref_upstream_seq, 'AAAAAAAAA', max_mismatch=2):
        ref_flag = False
    if fl_flag and ref_flag:
        return True

def filter_FalsePosivite_TSS_reads(read_start_seq):
    TSO = 'TTATATGGG'
    FalsePosivite_Flag = False
    if len(read_start_seq) < 10:
        FalsePosivite_Flag = True
    if fuzzy_match(read_start_seq, TSO, max_mismatch=2):
        FalsePosivite_Flag = True
    return FalsePosivite_Flag

def bam2ployA_site(bam_file, ref, chr_data):
    chrom, read_list = chr_data
    results = []
    
    for read in read_list:
        if read.is_unmapped:
            continue

        read_name, chrom, APA_site, strand, ref_upstream_seq, ref_downstream_seq, read_tail, read_start_seq, TSS_site, read_alignment_start, read_alignment_end, read_sequence = get_sequence(read, ref)

        TSS_filter = 'nan'
        PAS_filter = 'nan'
        
        # Apply filters and record valid TSS/PAS
        if filter_FalsePosivite_APA_reads(ref_upstream_seq, read_tail):
            PAS_filter = APA_site
        
        if filter_FalsePosivite_TSS_reads(read_start_seq):
            TSS_filter = TSS_site
        
        # Output the results
        results.append((read_name, chrom, strand, TSS_site, TSS_filter, APA_site, PAS_filter))
        
    return results

def parallel_process_bam(bam_file, ref):
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # 按照染色体分组
    chr_data = {}
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        chrom = read.reference_name
        if chrom not in chr_data:
            chr_data[chrom] = []
        chr_data[chrom].append(read)
    
    # 使用多进程并行处理染色体
    with multiprocessing.Pool(processes=8) as pool:
        chr_data_items = list(chr_data.items())
        results = pool.starmap(bam2ployA_site, [(bam_file, ref, chr_data_item) for chr_data_item in chr_data_items])

    # 合并所有结果
    all_results = [item for sublist in results for item in sublist]
    return all_results

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description='TSS-APA位点鉴定')
    parser.add_argument('-i','--input', required=True, help="输入bam文件路径")
    parser.add_argument('-o','--outprefix', required=True, help="输出文件名称")
    parser.add_argument('-r','--reference', required=True, help="参考基因组路径")
    args = parser.parse_args()
    
    bam_file = args.input
    output_prefix = args.outprefix
    ref = Fasta(args.reference)

    # 并行处理 BAM 文件
    results = parallel_process_bam(bam_file, ref)

    # 输出结果到文件
    output_file = f'{output_prefix}_TSS_PAS.txt'
    with open(output_file, "w") as out_f:
        for out_list in results:
            out_f.write(f"{"\t".join(out_list)}\n")

if __name__ == "__main__":
    main()
