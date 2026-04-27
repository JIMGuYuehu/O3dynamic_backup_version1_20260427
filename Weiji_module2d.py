import numpy as np
import xarray as xr

def calc_parc_o3(var, p_top=30, p_bottom=100):
    """
    计算指定压力层范围的臭氧柱
    
    参数:
    var: xarray.DataArray, 臭氧数据
    p_top: float, 顶部压力层 (hPa)
    p_bottom: float, 底部压力层 (hPa)
    """
    m_air = 28.964/(6.022E23)
    g = 980.6

    if 'plev' in var.dims: 
        plev=var.plev
        level='plev'
    if 'lev' in var.dims: 
        plev=var.lev
        level='lev'
    if 'level' in var.dims: 
        plev=var.level
        level='level'
    
    time=var.time
    delta_p = np.zeros((len(time), len(plev)))
    
    if plev[0]<plev[len(plev)-1] and plev[len(plev)-1] <= 1000: 
        factor=100
        factor_2=1 # conversion Pa->hPa
    if plev[0]>plev[len(plev)-1] and plev[0] <=1000: 
        factor=100
        factor_2=1
    if plev[0]<plev[len(plev)-1] and plev[len(plev)-1] >1000: 
        factor=1
        factor_2=100
    if plev[0]>plev[len(plev)-1] and plev[0] >1000: 
        factor=1
        factor_2=100
    
    if plev[0]<plev[len(plev)-1]: # for pressure levels in hPa
        for levelx in range(1,len(plev)): 
            delta_p[:,levelx].fill(plev[levelx] - plev[levelx-1])

        weights_p = xr.DataArray(delta_p*factor, dims=['time',level], coords=[time,plev])
 
        O3 = var * weights_p * 10/ (g * m_air)
    
        if level=='plev': O3=O3.sel(plev=slice(p_top*factor_2, p_bottom*factor_2)) 
        if level=='lev': O3=O3.sel(lev=slice(p_top*factor_2, p_bottom*factor_2))
        if level=='level': O3=O3.sel(level=slice(p_top*factor_2, p_bottom*factor_2))

        O3 = O3.sum(dim=level)
        O3 = O3/2.687E16  #calculate DU
        
    if plev[0]>plev[len(plev)-1]: # for pressure levels in hPa
        for levelx in range(0,len(plev)-1): 
            delta_p[:,levelx].fill(plev[levelx] - plev[levelx+1])

        weights_p = xr.DataArray(delta_p*factor, dims=['time',level], coords=[time,plev])

        O3 = var * weights_p * 10/ (g * m_air)
        
        if level =='plev': O3=O3.sel(plev=slice(p_bottom*factor_2, p_top*factor_2)) 
        if level=='lev': O3=O3.sel(lev=slice(p_bottom*factor_2, p_top*factor_2)) 
        if level=='level': O3=O3.sel(level=slice(p_bottom*factor_2, p_top*factor_2)) 
            
        O3 = O3.sum(dim=level)
        O3 = O3/2.687E16  #calculate DU

    return O3.where(O3 != 0)