#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# B2000WCN.NOCOUPL 多变量提取脚本 (保持与 B2000WCN001002 完全一致的结构)
# 目标：提取 U, V, T, OMEGA, PS, Z3, O3 及其坐标变量，拆分文件夹存储
# 说明：
#   - run001: 输出年号保持原内部年份
#   - run002: 输出年号 = 原内部年份 + 104
#   - 这样可与现有 B2000WCN001002 的目录结构/命名习惯保持一致
# ==============================================================================

# --- 配置路径 ---
SRC_RUN001="/mnt/backup_ETH/EXTR_2000/EXTR_2000/B2000WCN.NOCOUPL.e122.f19_g16.001/atm/hist"
SRC_RUN002="/mnt/backup_ETH/EXTR_2000/EXTR_2000/B2000WCN.NOCOUPL.e122.f19_g16.002/atm/hist"

# 建议单独输出，避免和 interactive 数据混淆
OUT_ROOT="/mnt/soclim0/public_data/weiji/B2000WCN_NOCOUPL001002"

# 与你原始 B2000WCN 保持一致：run002 文件年号整体 +104
OFFSET_RUN2=104
MAX_JOBS=32

# 核心变量：每个变量单独文件夹
CORE_VARS=("U" "V" "T" "OMEGA" "PS" "Z3" "O3")

# 必须保留的坐标变量（保证单个文件可直接用于后续分析）
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

    # 计算平移后的输出年份
    local y_orig_int=$((10#${y_orig_str}))
    local y_new_int=$((y_orig_int + year_offset))
    local y_new_str
    printf -v y_new_str "%04d" "${y_new_int}"

    local out_dir="${OUT_ROOT}/${var_name}"
    local out_file="${out_dir}/B2000WCN.NOCOUPL.sample.cam.h3.${y_new_str}.${var_name}.nc"

    # 断点续传：若文件已存在且非空，则跳过
    if [[ -s "${out_file}" ]]; then
        return
    fi

    # 匹配该年所有 h3 文件
    # 例如：
    # B2000WCN.NOCOUPL.e122.f19_g16.001.cam.h3.0001-01-01-00000.nc
    local pattern="${src_dir}/B2000WCN.NOCOUPL.e122.f19_g16.${run_label}.cam.h3.${y_orig_str}-"*

    shopt -s nullglob
    local files=( ${pattern} )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No files found for Year ${y_orig_str} in Run ${run_label}" >&2
        return
    fi

    # 使用临时文件，避免中断后留下损坏文件
    local tmp_file="${out_file}.tmp"

    # 提取变量并按时间拼接
    if ncrcat -O -v "${var_name},${COORD_VARS}" "${files[@]}" "${tmp_file}" > /dev/null 2>&1; then
        # 压缩
        ncks -4 -L 1 -O "${tmp_file}" "${tmp_file}"
        mv "${tmp_file}" "${out_file}"
        echo "[SUCCESS] ${var_name} Year ${y_new_str} (from Run${run_label} Y${y_orig_str})"
    else
        echo "[ERROR] Failed to process ${var_name} in Year ${y_orig_str} (Run ${run_label})" >&2
        rm -f "${tmp_file}"
    fi
}

export -f process_year_var
export OUT_ROOT COORD_VARS

# ------------------------------------------------------------------------------
# 获取年份列表
# ------------------------------------------------------------------------------
get_years() {
    # 优先匹配 .nc；若你的文件仍是 .nc.extr.nc，这个 glob 同样可覆盖
    ls "$1"/B2000WCN.NOCOUPL.e122.f19_g16.*.cam.h3.[0-9][0-9][0-9][0-9]-*.nc.extr.nc 2>/dev/null | \
    sed -E 's/.*\.cam\.h3\.([0-9]{4})-.*/\1/' | sort -u
}

# ------------------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------------------
echo "Creating directories..."
for var in "${CORE_VARS[@]}"; do
    mkdir -p "${OUT_ROOT}/${var}"
done

echo "Preparing task list..."
TASK_FILE=$(mktemp)

# Run 001: 年号保持不变
for Y in $(get_years "${SRC_RUN001}"); do
    for V in "${CORE_VARS[@]}"; do
        echo "process_year_var 001 ${SRC_RUN001} ${Y} ${V} 0" >> "${TASK_FILE}"
    done
done

# Run 002: 年号 +104，保持与原 B2000WCN001002 结构一致
for Y in $(get_years "${SRC_RUN002}"); do
    for V in "${CORE_VARS[@]}"; do
        echo "process_year_var 002 ${SRC_RUN002} ${Y} ${V} ${OFFSET_RUN2}" >> "${TASK_FILE}"
    done
done

echo "Starting parallel execution (Max ${MAX_JOBS} jobs)..."
cat "${TASK_FILE}" | xargs -P "${MAX_JOBS}" -I {} bash -c "{}"

rm -f "${TASK_FILE}"
echo "All done! Data path: ${OUT_ROOT}"