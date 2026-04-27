#!/bin/bash
# =============================================================================
# 整合版 Hindcast 并行处理脚本
# 功能：按变量提取、合并并压缩 Hindcast 子实验数据
# =============================================================================

set -euo pipefail

# ---------- 1. 配置区域 ----------
INPUT_BASE="/mnt/backup_ETH/lens"
OUTPUT_BASE="/mnt/soclim0/public_data/weiji/Hindcast"
TMP_BASE="${OUTPUT_BASE}/_tmp"

MAX_JOBS=32  # 根据服务器性能调整

# 变量配置（同步自 BWCN 脚本）
CORE_VARS=("U" "V" "T" "CLOX" "H2O" "Z3" "O3" "PS")
COORD_VARS="P0,hyai,hyam,hybi,hybm,date,time,datesec,time_bnds,lat,lon,lev,ilev,gw"

# 导出变量和函数供子 shell 使用
export OUTPUT_BASE TMP_BASE COORD_VARS
export CORE_VARS_STR=$(IFS=,; echo "${CORE_VARS[*]}")

# ---------- 2. 核心处理函数 ----------
process_member() {
    local prefix_path="$1"   # 完整路径前缀
    local dir_name="$2"      # 子实验名 (例如 0001-01)
    local member_name        # Member 文件名前缀
    member_name=$(basename "$prefix_path")
    
    IFS=',' read -r -a v_array <<< "$CORE_VARS_STR"

    # 查找该 Member 的所有分段文件
    shopt -s nullglob
    local files=( "${prefix_path}".*.nc* )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        return
    fi

    for var in "${v_array[@]}"; do
        local dest_dir="${OUTPUT_BASE}/${dir_name}/${var}"
        mkdir -p "$dest_dir"
        
        local final_out="${dest_dir}/${member_name}.${var}.nc"
        
        # 跳过已存在的文件
        if [[ -s "${final_out}" ]]; then continue; fi

        # 提取、合并与压缩 (参考 BWCN 逻辑)
        # 使用 .tmp.${$} 确保进程间文件名不冲突
        if ncrcat -O -v "${var},${COORD_VARS}" "${files[@]}" "${final_out}.tmp.${$}" > /dev/null 2>&1; then
            ncks -4 -L 1 -O "${final_out}.tmp.${$}" "${final_out}"
            rm -f "${final_out}.tmp.${$}"
        else
            echo "    [Error] Variable ${var} not found in ${member_name}"
        fi
    done
    
    echo "  [Done] ${member_name}"
}
export -f process_member

# ---------- 3. 主逻辑 ----------
mkdir -p "$TMP_BASE"

echo "Scanning $INPUT_BASE ..."
# 匹配类似 0001-01 这种子实验目录
dirs=( ${INPUT_BASE}/[0-9][0-9][0-9][0-9]-[0-9][0-9]* )

if [ ${#dirs[@]} -eq 0 ]; then
    echo "❌ 未找到符合条件的子实验目录！"
    exit 1
fi

for case_dir in "${dirs[@]}"; do
    [ -d "$case_dir" ] || continue
    dir_name=$(basename "$case_dir")
    
    echo "======================================================="
    echo "Processing Case: $dir_name"
    
    # 提取 Member 前缀列表
    # 逻辑：找 h3 文件 -> 去掉结尾的时间戳和后缀 -> 取唯一值
    prefix_list=$(find "$case_dir" -maxdepth 1 -name "*.cam.h3.*.nc*" | \
        sed -E 's/\.[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}\.nc.*$//' | \
        sort | uniq)
    
    if [[ -z "$prefix_list" ]]; then
        echo "  [Skip] No h3 files found in $dir_name"
        continue
    fi

    # 并行控制
    for prefix in $prefix_list; do
        process_member "$prefix" "$dir_name" &
        
        # 限制并发数
        if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
            wait -n
        fi
    done

    wait # 等待当前子实验的所有 Member 处理完再进下一个
    echo "Successfully processed case: $dir_name"
done

# 清理
rm -rf "$TMP_BASE"
echo "--------------- 所有任务完成 ---------------"