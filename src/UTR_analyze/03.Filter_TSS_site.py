import sys
import multiprocessing as mp

def process_chunk(lines):
    """
    处理一组数据行，根据特定规则筛选和去重读段数据
    
    参数:
        lines: 包含制表符分隔数据的字符串列表，每行代表一个读段的注释信息
    
    返回值:
        list: 包含元组的列表，每个元组包含(读段名称, 原始行数据)，按照筛选规则选出的唯一读段
    """
    seen_reads = set()
    selected_lines = []

    # 遍历所有行数据，根据匹配规则筛选读段
    for line in lines:
        line = line.strip()
        lineL = line.split('\t')

        # 解析读段信息字段，提取读段名称、转录本ID和基因名称列表
        read_infolist = lineL[3].split(';')
        read_name = read_infolist[-2]
        if read_name in seen_reads:
            continue
        transcript_id = read_infolist[-3]
        gene_name_list = read_infolist[1:-3]

        # 解析参考信息字段，提取参考转录本/基因ID和名称
        ref_infolist = lineL[9].split(';')
        Flag,ref_gene_name,ref_transgene_id = False,"nan","nan"
        if ref_infolist[0] == '.':
            pass
        elif len(ref_infolist) == 1 :
            Flag = True
        else:
            ref_transgene_id = ref_infolist[3]
            ref_gene_name = ref_infolist[1]

        # 根据匹配规则选择读段：转录本ID匹配或基因名称匹配的情况
        if Flag or transcript_id == ref_transgene_id:
            selected_lines.append((read_name, line))
            seen_reads.add(read_name)
        elif transcript_id =='.' and ref_gene_name in gene_name_list:
            selected_lines.append((read_name, line))
            seen_reads.add(read_name)
        else:
            continue

    return selected_lines

def chunkify_lines(lines, num_chunks):
    chunk_size = len(lines) // num_chunks
    return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

def main():
    TSS_FilterInfo_file = sys.argv[1]
    output_file = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else mp.cpu_count()

    # 读取所有行
    with open(TSS_FilterInfo_file, 'r') as f:
        lines = f.readlines()

    # 分块并并行处理
    chunks = chunkify_lines(lines, num_processes)
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_chunk, chunks)

    # 合并结果并全局去重
    global_seen = set()
    final_lines = []
    for chunk_result in results:
        for read_name, line in chunk_result:
            if read_name not in global_seen:
                final_lines.append(line)
                global_seen.add(read_name)

    # 写入输出文件
    with open(output_file, 'w') as out:
        out.write('\n'.join(final_lines) + '\n')

if __name__ == '__main__':
    main()
