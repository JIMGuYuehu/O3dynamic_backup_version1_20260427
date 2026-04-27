#!/usr/bin/env bash
set -euo pipefail

# ==========================
# B2000WCN Z3 提取脚本（run001 + run002）
#
# 目标：
#   - 从两个 run:
#       B2000WCN.e122.f19_g16.001
#       B2000WCN.e122.f19_g16.002
#     的 h3 文件中抽取 Z3（以及 hybrid→pressure 需要的坐标变量）
#   - 每个模型年拼成一个年度文件（保留原时间分辨率）
#   - 对 run002 年份整体平移 +104 年，使两个 run 样本年份连续：
#       run001: 0001–0104
#       run002: 0105–0210
#
# 输出目录：
#   /home/weiji/restart_exam/longrun_B2000WCN_withchem_data/Z3
# ==========================

SRC_RUN001="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.001/atm/hist"
SRC_RUN002="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.002/atm/hist"

OUT_ROOT="/home/weiji/restart_exam/longrun_B2000WCN_withchem_data"
OUT_Z3="${OUT_ROOT}/Z3"
mkdir -p "${OUT_Z3}"

OFFSET_RUN2=104

# Z3 + hybrid→pressure 必备 + 坐标/时间信息
VARS_KEEP="Z3,P0,PS,hyai,hyam,hybi,hybm,date,time,datesec,time_bnds,lat,lon,lev,ilev,gw"

detect_year_range () {
  local src_dir="$1"
  # 文件名形如 B2000WCN.e122.f19_g16.00X.cam.h3.YYYY-*.nc.extr.nc
  ls "${src_dir}"/B2000WCN.e122.f19_g16.*.cam.h3.[0-9][0-9][0-9][0-9]-*.nc.extr.nc \
    2>/dev/null | \
    sed -E 's/.*\.cam\.h3\.([0-9]{4})-.*/\1/' | \
    sort -u
}

process_run () {
  local run_label="$1"   # "001" 或 "002"
  local src_dir="$2"
  local year_offset="$3"

  echo "==> Processing run ${run_label} in ${src_dir}"
  local years
  mapfile -t years < <(detect_year_range "${src_dir}")

  if [[ "${#years[@]}" -eq 0 ]]; then
    echo "    ⚠️ No h3 files found in ${src_dir}, skip."
    return
  fi

  echo "    Detected years: ${years[0]} ... ${years[-1]} (N=${#years[@]}); year_offset=${year_offset}"

  for YYYY in "${years[@]}"; do
    local y_orig_str="${YYYY}"
    local y_orig_int=$((10#${YYYY}))

    local y_new_int=$((y_orig_int + year_offset))
    local y_new_str
    printf -v y_new_str "%04d" "${y_new_int}"

    echo "  -> Year ${y_orig_str} (run${run_label})  -->  sample year ${y_new_str}"

    local pattern="${src_dir}/B2000WCN.e122.f19_g16.${run_label}.cam.h3.${y_orig_str}-"*
    shopt -s nullglob
    local files=( ${pattern} )
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
      echo "     ⚠️ No files for year ${y_orig_str}, skip."
      continue
    fi

    # 用“新年份”命名，方便两段拼起来连续
    local out_file="${OUT_Z3}/B2000WCN.sample.cam.h3.${y_new_str}.Z3.hybrid.nc"

    ncrcat -O -v "${VARS_KEEP}" "${files[@]}" "${out_file}"
    ncks -4 -L 1 -O "${out_file}" "${out_file}"

    local nt
    nt=$(ncks -m "${out_file}" | awk '/^time dimension/ {print $NF}')
    echo "     --> ${out_file} (time=${nt})"
  done
}

process_run "001" "${SRC_RUN001}" 0
process_run "002" "${SRC_RUN002}" "${OFFSET_RUN2}"

echo "All done (B2000WCN Z3 yearly hybrid files, both runs)."
