import sys
import multiprocessing as mp

def process_lines(lines):
    """
    处理输入的行数据，过滤并返回有效的行记录
    
    参数:
        lines: 可迭代对象，包含待处理的行数据
        
    返回值:
        list: 包含元组的列表，每个元组包含(read_name, line)，表示有效的读段名称和对应的行数据
    """
    seen_reads = set()
    valid_lines = []

    # 遍历每一行数据，解析并过滤有效记录
    for line in lines:
        line = line.strip()
        lineL = line.split('\t')

        # 解析读段信息列表，提取关键字段
        read_infolist = lineL[3].split(';')

        read_name = read_infolist[-2]
        # 如果该读段已经处理过，则跳过
        if read_name in seen_reads:
            continue

        Flag = read_infolist[-1]
        transcript_id = read_infolist[-3]
        gene_name_list = read_infolist[1:-3]

        # 尝试解析参考信息列表，提取参考转录基因ID和基因名称
        try:
            ref_infolist = lineL[9].split(';')
            ref_transgene_id = ref_infolist[3]
            gene_name = ref_infolist[1]
        except IndexError:
            ref_transgene_id = "nan"
            gene_name = "nan"

        # 根据条件判断是否为有效行记录
        if Flag == 'T' or transcript_id == ref_transgene_id:
            # 当标志为'T'或转录ID匹配时，添加到有效行列表
            if read_name not in seen_reads:
                valid_lines.append((read_name, line))
                seen_reads.add(read_name)
        elif transcript_id =='.' and gene_name in gene_name_list:
            # 当转录ID为'.'且基因名在基因名列表中时，添加到有效行列表
                valid_lines.append((read_name, line))
                seen_reads.add(read_name)
        else:
            continue

    return valid_lines

def chunkify_lines(lines, num_chunks):
    chunk_size = len(lines) // num_chunks
    return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

def main():
    """
    主函数，用于处理TES过滤信息文件，通过多进程去重后输出结果
    
    参数:
        sys.argv[1]: TES_FilterInfo_file - 输入的TES过滤信息文件路径
        sys.argv[2]: output_file - 输出文件路径
        sys.argv[3]: num_processes - 可选参数，指定进程数量，默认为CPU核心数
    
    返回值:
        无返回值，直接将处理结果写入输出文件
    """
    TES_FilterInfo_file = sys.argv[1]
    output_file = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else mp.cpu_count()

    # 读取所有行
    with open(TES_FilterInfo_file, 'r') as f:
        lines = f.readlines()

    # 分块并多进程处理
    chunks = chunkify_lines(lines, num_processes)
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_lines, chunks)

    # 合并结果并全局去重
    global_seen = set()
    final_lines = []
    for res in results:
        for read_name, line in res:
            if read_name not in global_seen:
                final_lines.append(line)
                global_seen.add(read_name)

    # 输出到文件
    with open(output_file, 'w') as out:
        out.write('\n'.join(final_lines) + '\n')

if __name__ == '__main__':
    main()
