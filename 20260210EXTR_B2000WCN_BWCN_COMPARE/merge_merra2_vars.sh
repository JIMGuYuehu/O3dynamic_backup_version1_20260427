#!/usr/bin/env bash
set -euo pipefail

# --- 1. 配置路径与变量 ---
SRC_DIR="/mnt/soclim0/public_data/weiji/MERRA2M2I6NPANA"
OUT_ROOT="/mnt/soclim0/public_data/weiji/MERRA2_Processed_FixedTime"

CORE_VARS=("U" "V" "T" "O3")
KEEP_VARS_COMMON="time,time_bnds,lat,lon,lev"

# 128核服务器，保守一些
MAX_PARALLEL_YEARS=16

if [ ! -d "$SRC_DIR" ]; then
    echo "❌ 错误: 输入目录不存在: $SRC_DIR"
    exit 1
fi
mkdir -p "$OUT_ROOT"

process_year() {
    local y=$1
    local src=$2
    local out=$3
    local vars_str=$4

    IFS=',' read -r -a v_array <<< "$vars_str"

    shopt -s nullglob
    local all_files_in_year=( $(ls "${src}"/*.nc | grep -E "\.${y}[0-9]{4}\.SUB\.nc$") )
    local files=( $(printf '%s\n' "${all_files_in_year[@]}" | sort) )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "⚠️  [Year ${y}] 未找到文件，跳过。"
        return
    fi

    echo ">>> [Year ${y}] 正在处理 (${#files[@]} 天)..."

    for var in "${v_array[@]}"; do
        mkdir -p "${out}/${var}"

        local target="${out}/${var}/MERRA2.${var}.${y}.nc"
        local tmp_merge="${target}.mergetime.tmp.nc"
        local tmp_sel="${target}.sel.tmp.nc"

        if [[ -s "${target}" ]]; then
            echo "    [${y}][${var}] 已存在，跳过"
            continue
        fi

        echo "    [${y}][${var}] cdo mergetime ..."
        # 先合并时间
        if ! cdo -L -f nc4c -z zip_1 mergetime "${files[@]}" "${tmp_merge}" ; then
            echo "    ❌ [Year ${y}] ${var} mergetime 失败"
            rm -f "${tmp_merge}" "${tmp_sel}"
            continue
        fi

        echo "    [${y}][${var}] ncks select ..."
        # 再只保留需要的变量和坐标
        if ! ncks -4 -L 1 -O -v "${var},${KEEP_VARS_COMMON}" "${tmp_merge}" "${tmp_sel}" ; then
            echo "    ❌ [Year ${y}] ${var} 变量筛选失败"
            rm -f "${tmp_merge}" "${tmp_sel}"
            continue
        fi

        mv "${tmp_sel}" "${target}"
        rm -f "${tmp_merge}"

        # 简单验证
        echo "    [${y}][${var}] 时间检查:"
        cdo showtimestamp "${target}" | head -2 || true
    done

    echo ">>> [Year ${y}] 完成."
}
export -f process_year

echo "正在探测数据年份..."
YEARS=$(ls "${SRC_DIR}"/*.nc | grep -oE '\.[0-9]{8}\.' | sed 's/\.//g' | cut -c1-4 | sort -u)

if [[ -z "$YEARS" ]]; then
    echo "❌ 无法识别年份，请检查文件名格式。"
    exit 1
fi

echo "识别到年份: $(echo $YEARS | xargs)"
VARS_JOINED=$(IFS=,; echo "${CORE_VARS[*]}")

echo "启动并行处理 (Threads: ${MAX_PARALLEL_YEARS})..."
echo "$YEARS" | xargs -n 1 -P "${MAX_PARALLEL_YEARS}" -I {} \
    bash -c "process_year '{}' '${SRC_DIR}' '${OUT_ROOT}' '${VARS_JOINED}'"

echo "----------------------------------------------------"
echo "✅ 全部变量按年整合完毕！"
echo "结果目录: ${OUT_ROOT}"