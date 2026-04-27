#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# B2000WCN EP Flux 变量提取脚本 (run001 + run002)
#
# 目标：
#   - 从两个 run (001, 002) 提取计算 EP Flux 所需的全部变量
#   - 变量包括: U, V, T, OMEGA, PS, P0 (以及混合坐标系数和时间/空间坐标)
#   - 两个 Run 的年份处理：
#       run001: 0001–0104 (保持原样)
#       run002: 0001–0106 -> 平移 +104 年 -> 0105–0210
#
# 输出结构：
#   /home/weiji/restart_exam/longrun_B2000WCN_withchem_data/EP_Daily_Inputs/
#     └── B2000WCN.sample.cam.h3.{YYYY}.EP_vars.hybrid.nc
# ==============================================================================

# --- 输入路径 ---
SRC_RUN001="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.001/atm/hist"
SRC_RUN002="/mnt/backup_ETH/extr_2000/extr_2000/B2000WCN.e122.f19_g16.002/atm/hist"

# --- 输出路径 ---
OUT_ROOT="/home/weiji/restart_exam/longrun_B2000WCN_withchem_data"
OUT_DIR="${OUT_ROOT}/EP_Daily_Inputs"  # 专门存放 EP Flux 计算用的多变量文件
mkdir -p "${OUT_DIR}"

# --- 核心设置 ---
OFFSET_RUN2=104

# 关键变量列表 (用于 Python Block C 计算 EP Flux)
# 物理量: U, V, T, OMEGA, PS
# 坐标/辅助: P0, hyai, hyam, hybi, hybm (垂直坐标), lat, lon (水平坐标), time, date...
VARS_KEEP="U,V,T,OMEGA,PS,P0,hyai,hyam,hybi,hybm,date,time,datesec,time_bnds,lat,lon,lev,ilev,gw"

# ------------------------------------------------------------------------------
# 函数定义
# ------------------------------------------------------------------------------

detect_year_range () {
  local src_dir="$1"
  # 查找 h3 文件并提取年份
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

    # 计算新的连续年份
    local y_new_int=$((y_orig_int + year_offset))
    local y_new_str
    printf -v y_new_str "%04d" "${y_new_int}"

    # 输出文件路径
    local out_file="${OUT_DIR}/B2000WCN.sample.cam.h3.${y_new_str}.EP_vars.hybrid.nc"
    
    # 如果文件已存在且大小合理，跳过（断点续传）
    if [[ -f "${out_file}" ]]; then
        echo "    -> Year ${y_new_str} exists, skipping."
        continue
    fi

    echo "    -> Extracting Year ${y_orig_str} (run${run_label}) --> ${y_new_str}"

    # 匹配该年的所有文件
    local pattern="${src_dir}/B2000WCN.e122.f19_g16.${run_label}.cam.h3.${y_orig_str}-"*
    shopt -s nullglob
    local files=( ${pattern} )
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
      echo "     ⚠️ No files for year ${y_orig_str}, skip."
      continue
    fi

    # 核心操作：ncrcat 拼接时间维 + 提取特定变量
    # -O: 覆盖
    # -v: 只保留指定变量
    ncrcat -O -v "${VARS_KEEP}" "${files[@]}" "${out_file}"
    
    # 可选：压缩以节省空间 (Deflate level 1, Shuffle on)
    ncks -4 -L 1 -O "${out_file}" "${out_file}"

    # 简单的完整性检查
    local nt
    nt=$(ncks -m "${out_file}" | awk '/^time dimension/ {print $NF}')
    echo "       [OK] Saved to ${out_file} (days=${nt})"
  done
}

# ------------------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------------------

echo "Start extraction for EP Flux variables..."
echo "Variables: ${VARS_KEEP}"
echo "Target Dir: ${OUT_DIR}"
echo "----------------------------------------------------"

process_run "001" "${SRC_RUN001}" 0
process_run "002" "${SRC_RUN002}" "${OFFSET_RUN2}"

echo "----------------------------------------------------"
echo "All done. Files are ready in ${OUT_DIR}"