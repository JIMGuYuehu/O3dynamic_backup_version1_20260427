#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# MERRA2 daily SUB files -> yearly files (single variable per file)
#
# 输出：
#   /mnt/soclim0/public_data/weiji/MERRA2_Processed/U/MERRA2.U.YYYY.nc
#   /mnt/soclim0/public_data/weiji/MERRA2_Processed/V/MERRA2.V.YYYY.nc
#   /mnt/soclim0/public_data/weiji/MERRA2_Processed/T/MERRA2.T.YYYY.nc
#   /mnt/soclim0/public_data/weiji/MERRA2_Processed/O3/MERRA2.O3.YYYY.nc
#
# 核心流程：
#   1) cdo mergetime 合并全年 daily 文件
#   2) ncks 抽取单变量 + 空间坐标 + time
#   3) cdo settaxis 重建标准 daily 时间轴
#   4) ncatted 修正全局 attrs
#
# 时间轴：
#   daily, 从 YYYY-01-01 开始，步长 1day
# ============================================================

SRC_DIR="/mnt/soclim0/public_data/weiji/MERRA2M2I6NPANA"
OUT_ROOT="/mnt/soclim0/public_data/weiji/MERRA2_Processed"

MAX_PARALLEL_YEARS=8
TMP_ROOT="${OUT_ROOT}/_tmp"

mkdir -p "${OUT_ROOT}" "${TMP_ROOT}"

if [[ ! -d "${SRC_DIR}" ]]; then
    echo "❌ 输入目录不存在: ${SRC_DIR}"
    exit 1
fi

process_year() {
    local y="$1"
    local core_vars=("U" "V" "T" "O3")

    shopt -s nullglob
    local files=( "${SRC_DIR}"/MERRA2_*.inst6_3d_ana_Np."${y}"????.SUB.nc )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "⚠️  [Year ${y}] 未找到文件，跳过。"
        return
    fi

    IFS=$'\n' files=( $(printf '%s\n' "${files[@]}" | sort) )

    local nfiles="${#files[@]}"
    local first_date="${y}-01-01"
    local last_date
    last_date=$(date -d "${first_date} +$((nfiles-1)) day" +%F)

    echo ">>> [Year ${y}] 正在处理 (${nfiles} 天)..."
    echo "    First: ${files[0]}"
    echo "    Last : ${files[-1]}"

    for var in "${core_vars[@]}"; do
        mkdir -p "${OUT_ROOT}/${var}"

        local tmp_merge="${TMP_ROOT}/MERRA2.${var}.${y}.merge.nc"
        local tmp_var="${TMP_ROOT}/MERRA2.${var}.${y}.var.nc"
        local target="${OUT_ROOT}/${var}/MERRA2.${var}.${y}.nc"

        rm -f "${tmp_merge}" "${tmp_var}"

        if [[ -s "${target}" ]]; then
            echo "    [${y}][${var}] 已存在，跳过"
            continue
        fi

        echo "    [${y}][${var}] cdo mergetime ..."
        if ! cdo mergetime "${files[@]}" "${tmp_merge}" ; then
            echo "    ❌ [Year ${y}] ${var} mergetime 失败"
            rm -f "${tmp_merge}" "${tmp_var}"
            continue
        fi

        echo "    [${y}][${var}] ncks extract ..."
        if ! ncks -4 -L 1 -O -v "${var},lev,lat,lon,time" "${tmp_merge}" "${tmp_var}" ; then
            echo "    ❌ [Year ${y}] ${var} ncks 提取失败"
            rm -f "${tmp_merge}" "${tmp_var}"
            continue
        fi

        echo "    [${y}][${var}] cdo settaxis ..."
        if ! cdo settaxis,"${first_date}",00:00:00,1day "${tmp_var}" "${target}" ; then
            echo "    ❌ [Year ${y}] ${var} settaxis 失败"
            rm -f "${tmp_merge}" "${tmp_var}" "${target}"
            continue
        fi

        rm -f "${tmp_merge}" "${tmp_var}"

        echo "    [${y}][${var}] fix global attrs ..."
        # 删除不再适用于年文件的 granule 级 attrs
        ncatted -O \
            -a Filename,global,d,, \
            -a GranuleID,global,d,, \
            -a ProductionDateTime,global,d,, \
            -a RangeBeginningDate,global,o,c,"${first_date}" \
            -a RangeBeginningTime,global,o,c,"00:00:00.000000" \
            -a RangeEndingDate,global,o,c,"${last_date}" \
            -a RangeEndingTime,global,o,c,"00:00:00.000000" \
            -a TemporalRange,global,o,c,"${first_date} -> ${last_date}" \
            -a history,global,a,c," | yearly merged by cdo mergetime; settaxis=${first_date},00:00:00,1day; variable=${var}" \
            "${target}"

        echo "    [${y}][${var}] quick check:"
        cdo showtimestamp "${target}" | head -2 || true
    done

    echo ">>> [Year ${y}] 完成."
}

export -f process_year
export SRC_DIR OUT_ROOT TMP_ROOT

echo "正在探测数据年份..."
YEARS=$(ls "${SRC_DIR}"/MERRA2_*.inst6_3d_ana_Np.*.SUB.nc 2>/dev/null \
    | grep -oE '\.[0-9]{8}\.SUB\.nc$' \
    | sed -E 's/^\.([0-9]{4}).*/\1/' \
    | sort -u)

if [[ -z "${YEARS}" ]]; then
    echo "❌ 无法识别年份，请检查文件名格式。"
    exit 1
fi

echo "识别到年份: $(echo ${YEARS} | xargs)"
echo "启动并行处理 (Threads: ${MAX_PARALLEL_YEARS})..."

echo "${YEARS}" | xargs -n 1 -P "${MAX_PARALLEL_YEARS}" -I {} bash -c 'process_year "$@"' _ {}

echo "----------------------------------------------------"
echo "✅ 全部变量按年整合完毕！"
echo "结果目录: ${OUT_ROOT}"