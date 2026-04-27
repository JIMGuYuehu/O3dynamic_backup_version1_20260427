#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# B2000WCN 多变量提取脚本 (确保正确性优先)
# 目标：提取 U, V, T, OMEGA, PS, Z3, O3 及其坐标变量，拆分文件夹存储
# ==============================================================================

# --- 配置路径 ---
SRC_RUN001="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.001/atm/hist"
SRC_RUN002="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.002/atm/hist"
OUT_ROOT="/mnt/soclim0/public_data/weiji/B2000WCN001002"

OFFSET_RUN2=104
MAX_JOBS=32

# 变量列表：每个核心变量将拥有独立的文件夹
CORE_VARS=("U" "V" "T" "OMEGA" "PS" "Z3" "O3")
# 必须保留的坐标变量（确保单个文件可用性）
COORD_VARS="P0,hyai,hyam,hybi,hybm,date,time,datesec,time_bnds,lat,lon,lev,ilev,gw"

# ------------------------------------------------------------------------------
# 内部处理函数
# ------------------------------------------------------------------------------
process_year_var() {
    local run_label="$1"
    local src_dir="$2"
    local y_orig_str="$3"
    local var_name="$4"
    local year_offset="$5"

    # 计算平移后的年份
    local y_orig_int=$((10#${y_orig_str}))
    local y_new_int=$((y_orig_int + year_offset))
    local y_new_str
    printf -v y_new_str "%04d" "${y_new_int}"

    local out_dir="${OUT_ROOT}/${var_name}"
    local out_file="${out_dir}/B2000WCN.sample.cam.h3.${y_new_str}.${var_name}.nc"

    # 1. 断点续传：如果文件已存在且大小不为0，跳过
    if [[ -s "${out_file}" ]]; then
        return
    fi

    # 2. 匹配原始文件
    local pattern="${src_dir}/B2000WCN.e122.f19_g16.${run_label}.cam.h3.${y_orig_str}-"*
    shopt -s nullglob
    local files=( ${pattern} )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files found for Year ${y_orig_str} in ${run_label}" >&2
        return
    fi

    # 3. 提取并拼接 (使用临时文件确保原子性，防止程序中断产生破损文件)
    local tmp_file="${out_file}.tmp"
    
    # 核心动作：ncrcat 提取变量并拼接时间维
    if ncrcat -O -v "${var_name},${COORD_VARS}" "${files[@]}" "${tmp_file}" > /dev/null 2>&1; then
        # 4. 压缩处理
        ncks -4 -L 1 -O "${tmp_file}" "${tmp_file}"
        mv "${tmp_file}" "${out_file}"
        echo "[SUCCESS] ${var_name} Year ${y_new_str} (from Run${run_label} Y${y_orig_str})"
    else
        echo "[ERROR] Failed to process ${var_name} in Year ${y_orig_str}" >&2
        rm -f "${tmp_file}"
    fi
}

export -f process_year_var
export OUT_ROOT COORD_VARS

# ------------------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------------------

# 创建文件夹
echo "Creating directories..."
for var in "${CORE_VARS[@]}"; do mkdir -p "${OUT_ROOT}/${var}"; done

# 获取年份列表函数
get_years() {
    ls "$1"/B2000WCN.e122.f19_g16.*.cam.h3.[0-9][0-9][0-9][0-9]-*.nc.extr.nc 2>/dev/null | \
    sed -E 's/.*\.cam\.h3\.([0-9]{4})-.*/\1/' | sort -u
}

echo "Preparing task list..."
TASK_FILE=$(mktemp)

# Run 001
for Y in $(get_years "${SRC_RUN001}"); do
    for V in "${CORE_VARS[@]}"; do
        echo "process_year_var 001 ${SRC_RUN001} ${Y} ${V} 0" >> "${TASK_FILE}"
    done
done

# Run 002
for Y in $(get_years "${SRC_RUN002}"); do
    for V in "${CORE_VARS[@]}"; do
        echo "process_year_var 002 ${SRC_RUN002} ${Y} ${V} ${OFFSET_RUN2}" >> "${TASK_FILE}"
    done
done

echo "Starting parallel execution (Max ${MAX_JOBS} jobs)..."
cat "${TASK_FILE}" | xargs -P "${MAX_JOBS}" -I {} bash -c "{}"

rm -f "${TASK_FILE}"
echo "All done! Data path: ${OUT_ROOT}"