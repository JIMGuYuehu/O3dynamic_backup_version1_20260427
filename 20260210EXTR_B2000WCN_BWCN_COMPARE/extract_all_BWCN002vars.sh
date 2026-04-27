#!/usr/bin/env bash
set -euo pipefail

# --- 1. 严格检查输入路径是否存在 ---
SRC_DIR="/mnt/backup_ETH/EXTR_2000/EXTR_2000/BWCN.e122.f19_g16.002/atm/hist"
if [ ! -d "$SRC_DIR" ]; then
    echo "❌ 错误: 输入目录不存在: $SRC_DIR"
    exit 1
fi

OUT_ROOT="/mnt/soclim0/public_data/weiji/BWCN"
MAX_PARALLEL_YEARS=8
CORE_VARS=("U" "V" "T" "OMEGA" "Z3" "O3" "PS")
COORD_VARS="P0,hyai,hyam,hybi,hybm,date,time,datesec,time_bnds,lat,lon,lev,ilev,gw"

echo "Creating directories in ${OUT_ROOT}..."
for var in "${CORE_VARS[@]}"; do mkdir -p "${OUT_ROOT}/${var}"; done

# --- 2. 核心处理函数 ---
process_year() {
    local y_str=$1
    local src=$2
    local out=$3
    local vars_str=$4
    local coords=$5

    IFS=',' read -r -a v_array <<< "$vars_str"
    echo ">>> [BWCN Year ${y_str}] Starting..."
    
    # 更加宽松的通配符匹配
    local pattern="${src}/BWCN.e122.f19_g16.002.cam.h3.${y_str}-"*.nc*
    shopt -s nullglob
    local files=( ${pattern} )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "⚠️ Year ${y_str}: No files found with pattern."
        return
    fi

    for var in "${v_array[@]}"; do
        local target_file="${out}/${var}/BWCN.cam.h3.${y_str}.${var}.nc"
        if [[ -s "${target_file}" ]]; then continue; fi

        # 提取并压缩
        if ncrcat -O -v "${var},${coords}" "${files[@]}" "${target_file}.tmp" > /dev/null 2>&1; then
            ncks -4 -L 1 -O "${target_file}.tmp" "${target_file}"
            rm -f "${target_file}.tmp"
        else
            echo "    Error: Variable ${var} not found in Year ${y_str}"
        fi
    done
}
export -f process_year

# --- 3. 改进的年份探测 ---
echo "Detecting years..."
# 只要文件名里有 .h3. 和 4位数字- 就能抓出来
YEARS=$(ls "${SRC_DIR}"/*.h3.*.nc* 2>/dev/null | grep -oE '\.h3\.[0-9]{4}-' | grep -oE '[0-9]{4}' | sort -u || true)

if [[ -z "$YEARS" ]]; then
    echo "❌ 错误: 未能探测到任何年份！"
    echo "请检查该目录下是否有满足 *.h3.*.nc 格式的文件。"
    ls "${SRC_DIR}" | head -n 5
    exit 1
fi

echo "Found years: $(echo $YEARS | xargs)"
VARS_JOINED=$(IFS=,; echo "${CORE_VARS[*]}")

# --- 4. 执行并行任务 ---
echo "$YEARS" | xargs -n 1 -P "${MAX_PARALLEL_YEARS}" -I {} \
    bash -c "process_year '{}' '${SRC_DIR}' '${OUT_ROOT}' '${VARS_JOINED}' '${COORD_VARS}'"

echo "----------------------------------------------------"
echo "Done."