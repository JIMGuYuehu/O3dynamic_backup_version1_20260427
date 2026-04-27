#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

# 【死锁修复 1】：关闭底层数学库的隐式多线程，防止与 multiprocessing 打架
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# 【死锁修复 2】：关闭 HDF5 强制文件锁，防止多个进程读取 NetCDF 时互相锁死硬盘
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import re
import glob
import datetime
import multiprocessing  # 新增
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import gc
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# 【死锁修复 3】：将多进程的启动方式从默认的 'fork' 改为 'spawn'
# spawn 会启动一个完全干净的 Python 解释器，不会继承主进程的脏锁状态
try:
    multiprocessing.set_start_method('spawn', force=True)
except RuntimeError:
    pass

DATA_ROOT = "/mnt/soclim0/public_data/weiji"

HINDCAST_ROOT = os.path.join(DATA_ROOT, "Hindcast")
BWCN_ROOT = os.path.join(DATA_ROOT, "BWCN")
B2000_ROOT = os.path.join(DATA_ROOT, "B2000WCN001002")

OUT_ROOT = "/home/weiji/restart_exam/code/20260415egu/plots/hindcast/O3"

USE_CACHE = True
FORCE_REBUILD_O3 = False
FORCE_REBUILD_T = False
FORCE_REBUILD_U = False
FORCE_REBUILD_BWCN = False
FORCE_REBUILD_CLIM = False

# [防爆内存并行配置]
# 气象三维插值极其消耗内存。为了绝对安全，这里严格锁定为 4。
# 千万不要设成 10，否则一定会触发系统的 OOM Killer 导致进程死锁！
MAX_WORKERS = 4  

XLIM_START = 0
XLIM_END = 150

# plot styling
MEAN_LINEWIDTH = 5.0
MEMBER_LINEWIDTH_RATIO = 0.5
MEMBER_ALPHA = 0.28
SIGMA_ALPHA = 0.70   

O3_YLIM = (70, 160)
TMIN50_YLIM = (175, 255)
U60N10_YLIM = (-40, 80)

# O3 only: pressure range
O3_PTOP = 30.0
O3_PBOT = 70.0
O3_TAG = "30_70hPa"

# output subdirs
CACHE_ROOT = os.path.join(OUT_ROOT, "cache")
FIG_ROOT = os.path.join(OUT_ROOT, "figures")

for sub in [
    CACHE_ROOT,
    FIG_ROOT,
    os.path.join(CACHE_ROOT, "O3"),
    os.path.join(CACHE_ROOT, "Tmin50"),
    os.path.join(CACHE_ROOT, "U60N10"),
    os.path.join(CACHE_ROOT, "U60N50"),
    os.path.join(FIG_ROOT, "O3"),
    os.path.join(FIG_ROOT, "Tmin50"),
    os.path.join(FIG_ROOT, "U60N10"),
]:
    os.makedirs(sub, exist_ok=True)


# ============================================================
# FIGURE SPECS
# ============================================================

FIGURE_SPECS = [
    {
        "year": "0008",
        "prefix": "O3_evolution_fig1",
        "experiments": [
            ("0008-01", "0008Jan_small_pert", 0,  "forestgreen"),
            ("0008-02", "0008Feb_small_pert", 31, "royalblue"),
            ("0008-03", "0008Mar_small_pert", 59, "hotpink"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig2",
        "experiments": [
            ("0008-02_v2", "0008Feb_large_pert", 31, "forestgreen"),
            ("0008-03_v2", "0008Mar_large_pert", 59, "royalblue"),
            ("0008-04_v2", "0008Apr_large_pert", 90, "hotpink"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig3",
        "experiments": [
            ("0008-02",    "0008Feb_small_pert", 31, "forestgreen"),
            ("0008-02_v2", "0008Feb_large_pert", 31, "royalblue"),
            ("0008-02_v3", "0008Feb_moist_pert", 31, "hotpink"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig4",
        "experiments": [
            ("0008-03",    "0008Mar_small_pert", 59, "forestgreen"),
            ("0008-03_v2", "0008Mar_large_pert", 59, "royalblue"),
            ("0008-03_v3", "0008Mar_moist_pert", 59, "hotpink"),
        ],
    },
    {
        "year": "0003",
        "prefix": "O3_evolution_fig5",
        "experiments": [
            ("0003-02", "0003Feb_small_pert", 31, "forestgreen"),
            ("0003-03", "0003Mar_small_pert", 59, "royalblue"),
        ],
    },
    {
        "year": "0013",
        "prefix": "O3_evolution_fig6",
        "experiments": [
            ("0013-02", "0013Feb_small_pert", 31, "forestgreen"),
            ("0013-03", "0013Mar_small_pert", 59, "royalblue"),
        ],
    },
    {
        "year": "0014",
        "prefix": "O3_evolution_fig7",
        "experiments": [
            ("0014-02", "0014Feb_small_pert", 31, "forestgreen"),
            ("0014-03", "0014Mar_small_pert", 59, "royalblue"),
        ],
    },
    {
        "year": "0019",
        "prefix": "O3_evolution_fig8",
        "experiments": [
            ("0019-02", "0019Feb_small_pert", 31, "forestgreen"),
            ("0019-03", "0019Mar_small_pert", 59, "royalblue"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig9",
        "experiments": [
            ("0008-02_NOCOUPL", "0008Feb_nocouple", 31, "forestgreen"),
            ("0008-03_NOCOUPL", "0008Mar_nocouple", 59, "royalblue"),
        ],
    },
    {
        "year": "0013",
        "prefix": "O3_evolution_fig10",
        "experiments": [
            ("0013-02_NOCOUPL", "0013Feb_nocouple", 31, "forestgreen"),
            ("0013-03_NOCOUPL", "0013Mar_nocouple", 59, "royalblue"),
        ],
    },
    {
        "year": "0014",
        "prefix": "O3_evolution_fig11",
        "experiments": [
            ("0014-02_NOCOUPL", "0014Feb_nocouple", 31, "forestgreen"),
            ("0014-03_NOCOUPL", "0014Mar_nocouple", 59, "royalblue"),
        ],
    },
    {
        "year": "0019",
        "prefix": "O3_evolution_fig12",
        "experiments": [
            ("0019-02_NOCOUPL", "0019Feb_nocouple", 31, "forestgreen"),
            ("0019-03_NOCOUPL", "0019Mar_nocouple", 59, "royalblue"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig13",
        "experiments": [
            ("0008-02",         "0008Feb_small_pert", 31, "forestgreen"),
            ("0008-02_NOCOUPL", "0008Feb_nocouple",   31, "royalblue"),
        ],
    },
    {
        "year": "0008",
        "prefix": "O3_evolution_fig14",
        "experiments": [
            ("0008-03",         "0008Mar_small_pert", 59, "forestgreen"),
            ("0008-03_NOCOUPL", "0008Mar_nocouple",   59, "royalblue"),
        ],
    },
    {
        "year": "0013",
        "prefix": "O3_evolution_fig15",
        "experiments": [
            ("0013-02",         "0013Feb_small_pert", 31, "forestgreen"),
            ("0013-02_NOCOUPL", "0013Feb_nocouple",   31, "royalblue"),
        ],
    },
    {
        "year": "0013",
        "prefix": "O3_evolution_fig16",
        "experiments": [
            ("0013-03",         "0013Mar_small_pert", 59, "forestgreen"),
            ("0013-03_NOCOUPL", "0013Mar_nocouple",   59, "royalblue"),
        ],
    },
    {
        "year": "0014",
        "prefix": "O3_evolution_fig17",
        "experiments": [
            ("0014-02",         "0014Feb_small_pert", 31, "forestgreen"),
            ("0014-02_NOCOUPL", "0014Feb_nocouple",   31, "royalblue"),
        ],
    },
    {
        "year": "0014",
        "prefix": "O3_evolution_fig18",
        "experiments": [
            ("0014-03",         "0014Mar_small_pert", 59, "forestgreen"),
            ("0014-03_NOCOUPL", "0014Mar_nocouple",   59, "royalblue"),
        ],
    },
    {
        "year": "0019",
        "prefix": "O3_evolution_fig19",
        "experiments": [
            ("0019-02",         "0019Feb_small_pert", 31, "forestgreen"),
            ("0019-02_NOCOUPL", "0019Feb_nocouple",   31, "royalblue"),
        ],
    },
    {
        "year": "0019",
        "prefix": "O3_evolution_fig20",
        "experiments": [
            ("0019-03",         "0019Mar_small_pert", 59, "forestgreen"),
            ("0019-03_NOCOUPL", "0019Mar_nocouple",   59, "royalblue"),
        ],
    },
]


# ============================================================
# HELPERS
# ============================================================

def parse_member_id(path):
    base = os.path.basename(path)
    m = re.search(r"\.(\d{3})\.(?:[^\.]+\.)?cam\.h3", base)
    if m:
        return int(m.group(1))
    return None

def parse_b2000_file_year(path):
    base = os.path.basename(path)
    m = re.search(r"\.cam\.h3\.(\d{4})\.", base)
    if not m:
        raise ValueError(f"Cannot parse file year from: {base}")
    return int(m.group(1))

def get_lat_weights(ds_or_da, lat_slice):
    if "gw" in ds_or_da:
        return ds_or_da["gw"].sel(lat=lat_slice)
    lat = ds_or_da["lat"].sel(lat=lat_slice)
    return xr.DataArray(np.cos(np.deg2rad(lat)), dims=["lat"], coords={"lat": lat})

def find_bwcn_var_files_for_year(year_str, varname):
    root = os.path.join(BWCN_ROOT, varname)
    patterns = [
        os.path.join(root, f"BWCN.cam.h3.{year_str}.{varname}.nc"),
        os.path.join(root, f"BWCN.sample.cam.h3.{year_str}.{varname}.nc"),
        os.path.join(root, f"*{year_str}*.{varname}.nc"),
    ]
    for pat in patterns:
        files = sorted(glob.glob(pat))
        if files:
            return files
    raise FileNotFoundError(
        f"No BWCN {varname} files found for year {year_str} under {root}. "
    )

def get_hindcast_files(exp_name, varname):
    root = os.path.join(HINDCAST_ROOT, exp_name, varname)
    files = sorted(glob.glob(os.path.join(root, f"*.{varname}.nc")))
    if not files:
        raise FileNotFoundError(f"No {varname} files found under {root}")
    return files

def get_b2000_files(varname):
    root = os.path.join(B2000_ROOT, varname)
    files = sorted(glob.glob(os.path.join(root, f"B2000WCN.sample.cam.h3.*.{varname}.nc")))
    if not files:
        raise FileNotFoundError(f"No B2000WCN {varname} files found under {root}")
    return files


# ============================================================
# METRICS HELPERS (FWD & O3 Minimum)
# ============================================================

FW_LAT, FW_PLEV_HPA, FW_THRESHOLD, FW_MAX_CONSEC_WESTERLY = 60.0, 50.0, 7.0, 10

_base_date = datetime.datetime(2001, 1, 1)
MONTH_DAY_365 = [((_base_date + datetime.timedelta(days=i)).month, (_base_date + datetime.timedelta(days=i)).day) for i in range(365)]

def find_fw_for_one_year(vals_daily):
    vals_jj = np.asarray(vals_daily[:181], dtype=np.float64)
    for day_idx in range(len(vals_jj)):
        if not np.isfinite(vals_jj[day_idx]): continue
        if vals_jj[day_idx] < FW_THRESHOLD:
            has_long_westerly = False
            for i in range(day_idx + 1, max(day_idx + 1, len(vals_jj) - FW_MAX_CONSEC_WESTERLY + 1)):
                seg = vals_jj[i:i + FW_MAX_CONSEC_WESTERLY]
                if len(seg) < FW_MAX_CONSEC_WESTERLY: break
                if np.all(seg > FW_THRESHOLD):
                    has_long_westerly = True; break
            if not has_long_westerly:
                doy = day_idx + 1
                return doy, MONTH_DAY_365[doy - 1][0], MONTH_DAY_365[doy - 1][1]
    return np.nan, np.nan, np.nan

def get_o3_ma_min_for_member(o3_da_1d, apply_o3_5d=True):
    if apply_o3_5d:
        o3_da_1d = o3_da_1d.rolling(time=5, center=True, min_periods=5).mean()
    
    mask = o3_da_1d.time.dt.month.isin([3, 4])
    seg = o3_da_1d.where(mask, drop=True)
    
    vals = seg.values
    valid = np.isfinite(vals)
    if not np.any(valid):
        return np.nan, np.nan, np.nan
    
    idx = np.nanargmin(vals)
    min_val = float(vals[idx])
    min_month = int(seg.time.dt.month.values[idx])
    min_day = int(seg.time.dt.day.values[idx])
    
    return min_val, min_month, min_day


# ============================================================
# HYBRID CORE
# ============================================================

def calc_parc_o3_hybrid(ds_sub, p_top_hpa, p_bot_hpa, verbose=False):
    g = 9.80665
    M_air = 28.964 / 1000.0
    Na = 6.02214e23
    DU = 2.687e20
    factor = Na / (g * M_air * DU)

    P0 = ds_sub["P0"].squeeze(drop=True)
    PS = ds_sub["PS"]

    P_interface = ds_sub["hyai"] * P0 + ds_sub["hybi"] * PS
    p_i = P_interface.isel(ilev=slice(0, -1)).rename({"ilev": "lev"})
    p_ip1 = P_interface.isel(ilev=slice(1, None)).rename({"ilev": "lev"})

    if "lev" in ds_sub.coords:
        lev_vals = ds_sub["lev"].values
        p_i = p_i.assign_coords(lev=lev_vals)
        p_ip1 = p_ip1.assign_coords(lev=lev_vals)

    p_layer_top = xr.where(p_i < p_ip1, p_i, p_ip1)
    p_layer_bot = xr.where(p_i > p_ip1, p_i, p_ip1)

    pT = p_top_hpa * 100.0
    pB = p_bot_hpa * 100.0

    upper = xr.where(p_layer_top > pT, p_layer_top, pT)
    lower = xr.where(p_layer_bot < pB, p_layer_bot, pB)
    overlap = xr.where(lower > upper, lower - upper, 0.0)

    if "lev" in ds_sub.coords:
        overlap = overlap.assign_coords(lev=ds_sub["lev"].values)

    O3_col = (ds_sub["O3"] * overlap * factor).sum(dim="lev")

    if verbose:
        dp_eff_hpa = float((overlap.sum("lev") / 100.0).mean().compute())
        pass

    return O3_col

def _interp_profile_logp(var_prof, p_prof, p_target_pa):
    mask = np.isfinite(var_prof) & np.isfinite(p_prof) & (p_prof > 0)
    if mask.sum() < 2:
        return np.nan

    x = np.log(p_prof[mask])
    y = var_prof[mask]

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    xt = np.log(p_target_pa)
    if xt < x[0] or xt > x[-1]:
        return np.nan

    return np.interp(xt, x, y)

def interp_hybrid_to_pressure(ds, varname, p_target_hpa):
    P0 = ds["P0"].squeeze(drop=True)
    p_mid = ds["hyam"] * P0 + ds["hybm"] * ds["PS"]

    out = xr.apply_ufunc(
        _interp_profile_logp,
        ds[varname],
        p_mid,
        input_core_dims=[["lev"], ["lev"]],
        output_core_dims=[[]],
        vectorize=True,
        # 【修改 1】：移除 dask="parallelized" 避免多进程与 Dask 调度冲突
        kwargs={"p_target_pa": p_target_hpa * 100.0},
        output_dtypes=[float],
    )
    return out


# ============================================================
# DIAGNOSTIC FUNCTIONS (优化切片：防爆内存)
# ============================================================

def compute_o3_pc_ts(ds):
    # O3 需要积分，切片 60-90
    lat_slice = slice(60, 90)
    ds_sub = ds.sel(lat=lat_slice)
    O3_pc_3d = calc_parc_o3_hybrid(ds_sub, O3_PTOP, O3_PBOT)
    O3_zm = O3_pc_3d.mean(dim="lon")
    weights = get_lat_weights(ds_sub, lat_slice)
    O3_ts = O3_zm.weighted(weights).mean(dim="lat")
    O3_ts.name = "O3_pc_30_70hPa_60_90N"
    return O3_ts

def compute_polar_min_T50_ts(ds):
    # 【核心优化】：在三维插值前，砍掉全球 5/6 的无用数据
    # 【关键修改】：加上 .load()，将切片后的数据直接读入纯内存，完美避开 Dask 冲突
    ds_polar = ds.sel(lat=slice(60, 90)).load() 
    T50 = interp_hybrid_to_pressure(ds_polar, "T", 50.0)
    T50_zm = T50.mean(dim="lon")
    Tmin50 = T50_zm.min(dim="lat")
    Tmin50.name = "Tmin50_60_90N"
    return Tmin50

def compute_U60N10_ts(ds):
    # 【核心优化】：同上，仅截取 60N 附近的数据并加上 .load()
    ds_60N = ds.sel(lat=slice(55, 65)).load()
    U10 = interp_hybrid_to_pressure(ds_60N, "U", 10.0)
    U10_zm = U10.mean(dim="lon")
    U60N10 = U10_zm.interp(lat=60.0)
    U60N10.name = "U60N_10hPa"
    return U60N10

def compute_U60N50_ts(ds):
    # 【核心优化】：同上，加上 .load()
    ds_60N = ds.sel(lat=slice(55, 65)).load()
    U50 = interp_hybrid_to_pressure(ds_60N, "U", 50.0)
    U50_zm = U50.mean(dim="lon")
    U60N50 = U50_zm.interp(lat=60.0)
    U60N50.name = "U60N_50hPa"
    return U60N50


# ============================================================
# DATA LOADERS / CACHE
# ============================================================

def load_hindcast_diag(exp_name, varname, diag_name, compute_func,
                       use_cache=True, force_rebuild=False):
    cache_dir = os.path.join(CACHE_ROOT, diag_name)
    out_nc = os.path.join(cache_dir, f"{diag_name}_Hindcast_{exp_name}.nc")

    if use_cache and (not force_rebuild) and os.path.exists(out_nc):
        return xr.open_dataarray(out_nc)

    files = get_hindcast_files(exp_name, varname)
    member_ids = [parse_member_id(f) for f in files]
    
    # 【修改 2】：使用 with 上下文管理器确保文件句柄释放，显式指定 engine
    with xr.open_mfdataset(files, concat_dim="member", combine="nested", parallel=False, engine="netcdf4") as ds:
        ds = ds.assign_coords(member=("member", member_ids)).sortby("member")
        ts = compute_func(ds).compute()
        
    ts.to_netcdf(out_nc)
    
    # 【修改 3】：计算完毕后立刻清理内存
    del ds
    gc.collect()
    
    return ts

def load_bwcn_reference_diag(year_str, varname, diag_name, compute_func,
                             use_cache=True, force_rebuild=False):
    cache_dir = os.path.join(CACHE_ROOT, diag_name)
    out_nc = os.path.join(cache_dir, f"{diag_name}_BWCN_{year_str}.nc")

    if use_cache and (not force_rebuild) and os.path.exists(out_nc):
        return xr.open_dataarray(out_nc)

    files = find_bwcn_var_files_for_year(year_str, varname)
    
    with xr.open_mfdataset(files, combine="by_coords", parallel=False, engine="netcdf4") as ds:
        ts = compute_func(ds).compute()
        
    ts.to_netcdf(out_nc)
    
    del ds
    gc.collect()
    
    return ts

def build_b2000_climatology_diag(varname, diag_name, compute_func,
                                 use_cache=True, force_rebuild=False):
    cache_dir = os.path.join(CACHE_ROOT, diag_name)
    ts_nc = os.path.join(cache_dir, f"{diag_name}_B2000WCN_ts.nc")
    clim_nc = os.path.join(cache_dir, f"{diag_name}_B2000WCN_climatology.nc")

    if use_cache and (not force_rebuild) and os.path.exists(ts_nc) and os.path.exists(clim_nc):
        ts = xr.open_dataarray(ts_nc)
        clim = xr.open_dataarray(clim_nc)
        return ts, clim

    print(f"      -> Building {diag_name} climatology year-by-year (Safe RAM mode)...")
    all_files = get_b2000_files(varname)

    # 提取目标年份的文件
    valid_files = []
    for f in all_files:
        yr = parse_b2000_file_year(f)
        if (1 <= yr <= 103) or (105 <= yr <= 210):
            valid_files.append((yr, f))
    valid_files.sort(key=lambda x: x[0])

    ts_list = []

    # 【防爆内存核心】：逐年打开，逐年计算，算完立刻释放清理垃圾
    for yr, f in tqdm(valid_files, desc=f"      Processing B2000 {diag_name}"):
        with xr.open_dataset(f, engine="netcdf4") as ds:  # 使用上下文管理器防止缓存泄漏
            
            # 处理时间偏移
            if yr >= 105:
                offset = 103
                new_times = [t.replace(year=t.year + offset) for t in ds.time.values]
                ds = ds.assign_coords(time=new_times)

            # 在子集上计算时间序列
            ts_yr = compute_func(ds).compute()
            ts_list.append(ts_yr)
        
        # 强制回收内存
        gc.collect()

    # 最后只拼接一维的诊断序列 (非常小)
    ts = xr.concat(ts_list, dim="time").sortby("time")
    ts.to_netcdf(ts_nc)

    clim = ts.groupby("time.dayofyear").mean("time")
    clim.to_netcdf(clim_nc)

    return ts, clim


# ============================================================
# PLOTTING
# ============================================================

def plot_series_figure(ref_data, clim_data, experiments, year,
                       ylabel, ylim, out_subdir, output_filename_prefix,
                       xlim_start=0, xlim_end=150):
    all_months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    all_month_ticks = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

    xticks = [tick for tick in all_month_ticks if xlim_start <= tick <= xlim_end]
    xtick_labels = [all_months[i] for i, tick in enumerate(all_month_ticks) if tick in xticks]

    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)

    # reference
    nref = ref_data.sizes["time"]
    ref_end = min(xlim_end, nref)
    ax.plot(
        np.arange(xlim_start, ref_end),
        ref_data.isel(time=slice(xlim_start, ref_end)).values,
        color="black", linewidth=MEAN_LINEWIDTH, label="Reference"
    )

    # climatology
    clim_dim = list(clim_data.dims)[0]
    nclim = clim_data.sizes[clim_dim]
    clim_end = min(xlim_end, nclim)
    ax.plot(
        np.arange(xlim_start, clim_end),
        clim_data.isel({clim_dim: slice(xlim_start, clim_end)}).values,
        color="black", linestyle=":", linewidth=MEAN_LINEWIDTH, label="Climatology"
    )

    member_lw = MEAN_LINEWIDTH * MEMBER_LINEWIDTH_RATIO

    for exp in experiments:
        data = exp["data"]
        label = exp["label"]
        offset = exp["offset"]
        color = exp["color"]

        if offset >= xlim_end: continue

        total = data.sizes["time"]
        start_idx = max(0, xlim_start - offset)
        end_idx = min(total, xlim_end - offset)
        if end_idx <= start_idx: continue

        x = np.arange(offset + start_idx, offset + end_idx)
        sub = data.isel(time=slice(start_idx, end_idx))

        if "member" in sub.dims:
            plot_data = sub.transpose("member", "time")
            mean_line = plot_data.mean("member")
            std_line = plot_data.std("member")
            lower = mean_line - std_line
            upper = mean_line + std_line

            for i in range(plot_data.sizes["member"]):
                ax.plot(x, plot_data.isel(member=i).values, color=color, linewidth=member_lw, alpha=MEMBER_ALPHA)

            ax.fill_between(x, lower.values, upper.values, color=color, alpha=SIGMA_ALPHA)
            ax.plot(x, mean_line.values, color=color, linewidth=MEAN_LINEWIDTH, label=label)
        else:
            ax.plot(x, sub.values, color=color, linewidth=MEAN_LINEWIDTH, label=label)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize=16)
    ax.set_xlim(xlim_start, xlim_end)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_ylim(*ylim)
    ax.tick_params(axis="y", labelsize=16)

    ax.text(
        0.02, 0.95, f"Year: {year}",
        transform=ax.transAxes, fontsize=18, verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    )

    ax.legend(fontsize=13)

    fig_dir = os.path.join(FIG_ROOT, out_subdir)
    os.makedirs(fig_dir, exist_ok=True)

    fig.savefig(os.path.join(fig_dir, f"{output_filename_prefix}_{year}.pdf"))
    fig.savefig(os.path.join(fig_dir, f"{output_filename_prefix}_{year}.png"), dpi=200)

    return fig, ax


# ============================================================
# PARALLEL WRAPPER
# ============================================================

def process_single_hindcast(args):
    exp_name, varname, diag_name, compute_func_name, force_rebuild = args
    func_map = {
        "compute_o3_pc_ts": compute_o3_pc_ts,
        "compute_polar_min_T50_ts": compute_polar_min_T50_ts,
        "compute_U60N10_ts": compute_U60N10_ts,
        "compute_U60N50_ts": compute_U60N50_ts
    }
    compute_func = func_map[compute_func_name]
    
    try:
        data = load_hindcast_diag(
            exp_name=exp_name, varname=varname, diag_name=diag_name,
            compute_func=compute_func, use_cache=USE_CACHE, force_rebuild=force_rebuild
        )
        return exp_name, data
    except Exception as e:
        return exp_name, f"Error: {str(e)}"

def build_suite_parallel(varname, diag_name, compute_func_name, force_rebuild_flag, ylabel, ylim, out_subdir):
    print(f"\n================ {diag_name} SUITE (Parallel) ================\n")
    func_map = {
        "compute_o3_pc_ts": compute_o3_pc_ts,
        "compute_polar_min_T50_ts": compute_polar_min_T50_ts,
        "compute_U60N10_ts": compute_U60N10_ts
    }
    compute_func = func_map[compute_func_name]
    
    print(f"[{diag_name}] Step 1/3: Loading Climatology and Reference Data...")
    _, clim = build_b2000_climatology_diag(
        varname=varname, diag_name=diag_name, compute_func=compute_func,
        use_cache=USE_CACHE, force_rebuild=FORCE_REBUILD_CLIM or force_rebuild_flag
    )

    years_needed = sorted(set(spec["year"] for spec in FIGURE_SPECS))
    all_exp_names = sorted(set(exp_name for spec in FIGURE_SPECS for exp_name, _, _, _ in spec["experiments"]))
    
    ref_cache = {}
    for year in years_needed:
        ref_cache[year] = load_bwcn_reference_diag(
            year_str=year, varname=varname, diag_name=diag_name, compute_func=compute_func,
            use_cache=USE_CACHE, force_rebuild=FORCE_REBUILD_BWCN or force_rebuild_flag
        )

    print(f"[{diag_name}] Step 2/3: Parallel processing {len(all_exp_names)} Hindcast experiments on {MAX_WORKERS} cores...")
    exp_cache = {}
    tasks = [(exp, varname, diag_name, compute_func_name, force_rebuild_flag) for exp in all_exp_names]
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_single_hindcast, task): task[0] for task in tasks}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Calculating {diag_name}"):
            exp_name, result = future.result()
            if isinstance(result, str) and result.startswith("Error"):
                print(f"FAILED: {exp_name} with {result}")
            else:
                exp_cache[exp_name] = result

    print(f"[{diag_name}] Step 3/3: Generating figures...")
    for spec in tqdm(FIGURE_SPECS, desc=f"Plotting {diag_name}"):
        experiments = []
        for exp_name, label, offset, color in spec["experiments"]:
            if exp_name in exp_cache:
                experiments.append({
                    "data": exp_cache[exp_name],
                    "label": label, "offset": offset, "color": color,
                })

        prefix = spec["prefix"]
        if diag_name != "O3":
            prefix = prefix.replace("O3_evolution", f"{diag_name}_evolution")

        fig, ax = plot_series_figure(
            ref_data=ref_cache[spec["year"]], clim_data=clim, experiments=experiments,
            year=spec["year"], ylabel=ylabel, ylim=ylim, out_subdir=out_subdir,
            output_filename_prefix=prefix, xlim_start=XLIM_START, xlim_end=XLIM_END
        )
        plt.close(fig)


# ============================================================
# METRICS TEXT EXPORTER
# ============================================================

def process_metrics_for_exp(exp_name):
    try:
        o3_data = load_hindcast_diag(exp_name, "O3", "O3", compute_o3_pc_ts, use_cache=True, force_rebuild=False)
        u50_data = load_hindcast_diag(exp_name, "U", "U60N50", compute_U60N50_ts, use_cache=True, force_rebuild=False)
        
        out_file = os.path.join(HINDCAST_ROOT, exp_name, f"{exp_name}_metrics.txt")
        
        with open(out_file, "w") as f:
            f.write("Member\tFWD_DOY\tFWD_Date\tO3_MA_Min_Val\tO3_Min_Date\n")
            
            for m in o3_data.member.values:
                o3_m = o3_data.sel(member=m)
                u50_m = u50_data.sel(member=m)
                
                fwd_doy, fw_month, fw_day = find_fw_for_one_year(u50_m.values)
                fwd_date_str = f"{int(fw_month):02d}-{int(fw_day):02d}" if not np.isnan(fw_month) else "NaN"
                
                o3_val, o3_month, o3_day = get_o3_ma_min_for_member(o3_m)
                o3_date_str = f"{int(o3_month):02d}-{int(o3_day):02d}" if not np.isnan(o3_month) else "NaN"
                
                f.write(f"{m:03d}\t{fwd_doy}\t{fwd_date_str}\t{o3_val:.2f}\t{o3_date_str}\n")
        
        return exp_name, True
    except Exception as e:
        return exp_name, str(e)


# ============================================================
# MAIN
# ============================================================

def main():
    build_suite_parallel(
        varname="O3", diag_name="O3", compute_func_name="compute_o3_pc_ts", 
        force_rebuild_flag=FORCE_REBUILD_O3, ylabel="Partial ozone column, 30–70 hPa (DU)", 
        ylim=O3_YLIM, out_subdir="O3"
    )
    
    build_suite_parallel(
        varname="T", diag_name="Tmin50", compute_func_name="compute_polar_min_T50_ts", 
        force_rebuild_flag=FORCE_REBUILD_T, ylabel="Polar minimum T at 50 hPa (K)", 
        ylim=TMIN50_YLIM, out_subdir="Tmin50"
    )
    
    build_suite_parallel(
        varname="U", diag_name="U60N10", compute_func_name="compute_U60N10_ts", 
        force_rebuild_flag=FORCE_REBUILD_U, ylabel="Zonal-mean U at 60°N, 10 hPa (m s$^{-1}$)", 
        ylim=U60N10_YLIM, out_subdir="U60N10"
    )

    print("\n================ METRICS EXPORT (Parallel) ================\n")
    all_exp_names = sorted(set(exp_name for spec in FIGURE_SPECS for exp_name, _, _, _ in spec["experiments"]))
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_metrics_for_exp, exp): exp for exp in all_exp_names}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Exporting Metrics TXT"):
            exp_name, result = future.result()
            if result is not True:
                print(f"FAILED to export metrics for {exp_name}: {result}")

    print("\nAll tasks completed successfully.")

if __name__ == "__main__":
    main()