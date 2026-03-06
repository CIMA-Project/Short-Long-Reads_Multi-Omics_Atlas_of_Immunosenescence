#!/bin/bash

# ================= 配置区域 =================
# 1. 工作目录
WORK_DIR="/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/1_IgBLAST"
THRESHOLD_FILE="/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/csv/all_threshold_results.csv"
DEFAULT_THRESHOLD=0.12
IMGT_VDJ="/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/IMGT_reference"
MERGE_DIR="/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation"
# ===========================================

# 检查必要文件是否存在
if [ ! -d "$WORK_DIR" ]; then echo "错误: 工作目录不存在"; exit 1; fi
if [ ! -f "$THRESHOLD_FILE" ]; then echo "错误: 阈值文件不存在"; exit 1; fi

cd "$WORK_DIR" || exit

# 遍历样本目录
for sample_dir in */ ; do
    sample_name=$(basename "$sample_dir")
    
    echo "------------------------------------------------------"
    echo "正在处理样本: $sample_name"
    
    # 从 CSV 获取该样本的阈值
    # tr -d '\r' 是为了防止 Windows 格式的换行符导致变量读取错误
    CURRENT_THRESHOLD=$(awk -F, -v s="$sample_name" '$2 == s {print $1}' "$THRESHOLD_FILE" | tr -d '\r' | head -n 1)

    # 检查是否成功获取到了阈值
    if [[ -z "$CURRENT_THRESHOLD" || "$CURRENT_THRESHOLD" == "NA" || ! "$CURRENT_THRESHOLD" =~ ^[0-9.]+$ ]]; then
        echo "  -> ⚠️  警告: 未在CSV中找到有效阈值，将使用默认值: $DEFAULT_THRESHOLD"
        FINAL_DIST="$DEFAULT_THRESHOLD"
    else
        echo "  -> ✅ 成功: 从文件加载阈值: $CURRENT_THRESHOLD"
        FINAL_DIST="$CURRENT_THRESHOLD"
    fi
    echo "------------------------------------------------------"

    # 进入样本目录
    cd "$sample_dir" || continue

    input_file=$(ls filtered_contig_db.tsv 2>/dev/null | head -n 1)
    if [ -z "$input_file" ]; then
        echo "跳过: 未找到输入文件"
        cd "$WORK_DIR" && continue
    fi

    # --- Step 1: 筛选 Productive ---
    echo "[Step 1] 筛选 Productive 序列..."
    ParseDb.py select -d "$input_file" -f productive -u T > /dev/null
    step1_output="${input_file%.tsv}_parse-select.tsv"
    
    if [ ! -f "$step1_output" ]; then echo "Step 1 失败"; cd "$WORK_DIR"; continue; fi

    # --- Step 2: 分离重链和轻链 ---
    echo "[Step 2] 分离重链和轻链..."
    ParseDb.py select -d "$step1_output" -f v_call j_call c_call -u "IGH" --logic all --regex --outname "${sample_name}_heavy" > /dev/null
    ParseDb.py select -d "$step1_output" -f v_call j_call c_call -u "IG[LK]" --logic all --regex --outname "${sample_name}_light" > /dev/null

    heavy_file="${sample_name}_heavy_parse-select.tsv"
    light_file="${sample_name}_light_parse-select.tsv"

    if [[ ! -f "$heavy_file" || ! -f "$light_file" ]]; then echo "Step 2 失败"; cd "$WORK_DIR"; continue; fi

    # --- Step 3: 定义重链克隆 (使用动态阈值) ---
    echo "[Step 3] 定义重链克隆 (Threshold: $FINAL_DIST)..."
    
    DefineClones.py -d "$heavy_file" --act set --model ham \
        --norm len --dist "$FINAL_DIST" --outname "$sample_name" > /dev/null

    heavy_cloned_file="${sample_name}_clone-pass.tsv"

    if [ ! -f "$heavy_cloned_file" ]; then echo "Step 3 失败"; cd "$WORK_DIR"; continue; fi

    # --- Step 4: 轻链聚类校正 ---
    echo "[Step 4] 结合轻链信息..."
    light_cluster="${sample_name}_Final_Heavy_Light_Cloned.tsv"

    /media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/light_cluster.py -d "$heavy_cloned_file" \
        -e "$light_file" -o "$light_cluster"

    # --- Step 5: 重建种系序列 ---
    echo "[Step 5] 添加种系序列..."
    final_output="${sample_name}_Final_Heavy_Light_Cloned_germ-pass.tsv"
    
    /media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/CreateGermlines.py -d "$light_cluster" \
        -g dmask --cloned \
        -r "$IMGT_VDJ/IMGT_Human_IGHV.fasta" "$IMGT_VDJ/IMGT_Human_IGHD.fasta"  "$IMGT_VDJ/IMGT_Human_IGHJ.fasta" \
        -o "$final_output" > /dev/null

    echo "样本 $sample_name 全部处理完成"
    cd "$WORK_DIR" || exit
done

echo "=================================================="
echo "所有样本处理完毕，开始合并..."
# ================= 合并代码 (保持不变) =================
MERGED_FILE="${MERGE_DIR}/All_Samples_Merged_germ-pass.tsv"
find "$WORK_DIR" -name "*_Final_Heavy_Light_Cloned_germ-pass.tsv" > "${WORK_DIR}/file_list.txt"

if [ ! -s "${WORK_DIR}/file_list.txt" ]; then
    echo "错误：未找到文件，无法合并。"
    exit 1
fi

first_file=$(head -n 1 "${WORK_DIR}/file_list.txt")
head -n 1 "$first_file" > "$MERGED_FILE"

count=0
while read -r file; do
    tail -n +2 "$file" >> "$MERGED_FILE"
    ((count++))
done < "${WORK_DIR}/file_list.txt"

rm "${WORK_DIR}/file_list.txt"
echo "合并成功！共处理 $count 个样本。"
echo "最终文件: $MERGED_FILE"

csv_files <- list.files(
  path = "/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/data/", 
  pattern = "\\.tsv$",      
  full.names = TRUE        
)

data <- read.table("/media/AnalysisDisk2/Lifupeng/Lifupeng/ALSM/BCR/Immcantation/heavy_parse-select.tsv", header = TRUE, sep = "\t")