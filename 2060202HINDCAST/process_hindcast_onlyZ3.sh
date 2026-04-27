#!/bin/bash
# =============================================================================
# process_hindcast_parallel.sh
#
# 特性：
# 1. 并行处理：利用多核 CPU 加速。
# 2. 进程控制：通过 MAX_JOBS 控制并发数量，防止 I/O 拥堵。
# 3. 智能识别：沿用之前的智能前缀识别和多段合并逻辑。
# =============================================================================

set -u

# ---------- 配置区域 ----------
INPUT_BASE="/mnt/backup_ETH/lens"
OUTPUT_BASE="/home/weiji/restart_exam/hindcast_data"
TMP_BASE="${OUTPUT_BASE}/_tmp"

# 【关键配置】并发数量
# 建议从 20 开始尝试。如果磁盘读写速度很快，可以加到 50 或 100。
# 不要直接设为 200，除非你的存储是高性能并行文件系统 (Lustre/GPFS)。
MAX_JOBS=2

# 变量和压层
VARS=(Z3)
PLEVS="10,50,100,200,300,500,1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,85000,92500,100000"

# 调试模式
DRY_RUN=0

# 导出变量供函数内部使用
export OUTPUT_BASE TMP_BASE PLEVS DRY_RUN
# 数组无法直接 export，我们在函数内部重新定义或通过参数传递，
# 但由于函数在同一脚本内，子shell会继承这些变量的一份拷贝，所以通常可以直接读。

# ---------- 核心处理函数 (单个 Member) ----------
process_member() {
    # 接收参数
    local prefix_path="$1"
    local dir_name="$2"
    local case_tmp_dir="$3"
    
    # 重新定义数组 (子 Shell 安全性)
    local LOC_VARS=(Z3)
    local var_list
    var_list=$(echo "${LOC_VARS[*]}" | tr ' ' ',')

    local prefix_name
    prefix_name=$(basename "$prefix_path")
    
    # 查找属于该 Member 的所有时间段文件 (seg1, seg2, seg3...)
    # 必须重新 find/ls 确保即时性
    local files=( $(ls "${prefix_path}".*.nc | sort) )
    local num_segs=${#files[@]}

    # 简单日志 (由于并行，echo 可能会乱序，这是正常的)
    # echo "  [Start] $prefix_name ($num_segs segs)"

    if [ "$num_segs" -lt 2 ]; then
        echo "  [警告] $prefix_name 文件少于 2 个，跳过。"
        return
    fi

    # 1. 准备临时文件
    # 注意：在并行中，必须确保文件名绝对唯一。
    # 我们使用 prefix_name (它包含 unique id) 加上进程ID $$ 确保安全
    local tmp_files=()
    local t_file
    
    for f in "${files[@]}"; do
        local f_base
        f_base=$(basename "$f")
        t_file="${case_tmp_dir}/${f_base}.${$}.tmp.nc"
        
        if [ "$DRY_RUN" -eq 0 ]; then
            cp "$f" "$t_file"
            # ncatted 修改
            ncatted -O -a bounds,lev,c,c,ilev "$t_file"
        fi
        tmp_files+=("$t_file")
    done

    local merged_tmp="${case_tmp_dir}/${prefix_name}.merged.${$}.nc"

    # 2. CDO 合并 + 插值
    if [ "$DRY_RUN" -eq 0 ]; then
        cdo -s -O -mergetime \
            -ml2pl,"${PLEVS}" \
            -select,name="${var_list}" \
            "${tmp_files[@]}" "$merged_tmp" 2>/dev/null
        
        if [ $? -ne 0 ]; then
            echo "  [Error] CDO failed: $prefix_name"
            rm -f "${tmp_files[@]}" "$merged_tmp"
            return
        fi
    fi

    # 3. 拆分变量输出
    for var in "${LOC_VARS[@]}"; do
        local dest_dir="${OUTPUT_BASE}/${dir_name}/${var}"
        # 确保目录存在 (mkdir -p 是原子操作，多进程安全)
        mkdir -p "$dest_dir"

        local final_out="${dest_dir}/${prefix_name}.${var}.nc"

        if [ "$DRY_RUN" -eq 0 ]; then
            # -s 静默模式，减少日志冲突
            cdo -s -O -select,name="$var" "$merged_tmp" "$final_out" 2>/dev/null
        fi
    done

    # 4. 清理
    if [ "$DRY_RUN" -eq 0 ]; then
        rm -f "${tmp_files[@]}" "$merged_tmp"
    fi
    
    echo "  [Done] $prefix_name"
}

# 导出函数以便 xargs 或子 shell 调用 (虽然这里直接调用不需要 export -f，但这是好习惯)
export -f process_member


# ---------- 主逻辑 ----------
echo "正在扫描 $INPUT_BASE ..."
echo "并行任务数 (MAX_JOBS): $MAX_JOBS"

dirs=( ${INPUT_BASE}/00[0-9][0-9]-[0-9][0-9]* )

if [ ${#dirs[@]} -eq 0 ]; then
    echo "未找到目录！"
    exit 1
fi

for case_dir in "${dirs[@]}"; do
    [ -d "$case_dir" ] || continue
    dir_name=$(basename "$case_dir")
    
    echo "======================================================="
    echo "Processing Directory: $dir_name"
    
    # 创建临时前缀列表
    prefix_list_file="${TMP_BASE}/prefixes_${dir_name}.txt"
    mkdir -p "$TMP_BASE"
    
    # 提取 Member 前缀
    find "$case_dir" -maxdepth 1 -name "*.cam.h3.*.nc" | \
        sed -E 's/\.[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}\.nc$//' | \
        sort | uniq > "$prefix_list_file"
    
    num_members=$(wc -l < "$prefix_list_file")
    echo "  -> 发现 $num_members 个 Member，开始并行处理..."

    # 创建该 case 的临时目录
    case_tmp_dir="${TMP_BASE}/${dir_name}"
    mkdir -p "$case_tmp_dir"

    # --- 并行循环控制 ---
    job_count=0
    
    while read -r prefix_path; do
        # 在后台启动任务 (&)
        # 将必要的参数传递给函数
        process_member "$prefix_path" "$dir_name" "$case_tmp_dir" &
        
        # 计数器 +1
        job_count=$((job_count + 1))
        
        # 检查后台任务数量
        # `jobs -r -p` 列出正在运行的后台进程 PID
        # 如果数量 >= MAX_JOBS，则等待任意一个完成 (-n)
        if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
            wait -n
        fi
        
    done < "$prefix_list_file"

    # 等待该目录下剩余的所有后台任务完成，再处理下一个目录
    wait
    
    # 清理该目录的临时文件
    rm -rf "$case_tmp_dir"
    rm -f "$prefix_list_file"
    echo "  -> $dir_name 处理完毕。"

done

echo "=============== 所有任务完成 ==============="