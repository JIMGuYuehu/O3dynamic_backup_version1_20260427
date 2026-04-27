begin
  ;--------------------------------------------------
  ; 0. 创建输出目录（如果不存在）
  ;--------------------------------------------------
  outdir = "/home/weiji/restart_exam/plots/epflux/test"
  rc     = systemfunc("mkdir -p " + outdir)

  ;--------------------------------------------------
  ; 1. 打开数据并读取第一个时间步的 U, V, T
  ;--------------------------------------------------
  filename = "/mnt/backup_ETH/Marina/WACCM/CHEM_2000_restart/BWCN.e122.f19_g16.002_0008/Feb/BWCN.e122.f19_g16.002.0008-02.002.nc"
  f = addfile(filename, "r")

  u0 = f->U(0, :, :, :)    ; (plev, lat, lon)
  v0 = f->V(0, :, :, :)
  t0 = f->T(0, :, :, :)
  pl = f->plev             ; 压强 (Pa)
  la = f->lat              ; 纬度 (deg)

  ;--------------------------------------------------
  ; 2. 调用 epflux 计算 EP 通量及散度
  ;--------------------------------------------------
  opt = True
  opt@raw            = False
  opt@scale_sqrt_p   = True
  opt@magf           = False
  opt@print_var_info = False

  ep = epflux(u0, v0, t0, pl, la, opt)
  EPdiv = ep[2]         ; EP 通量散度 (pl, lat)

  ;--------------------------------------------------
  ; 3. 打开绘图设备
  ;--------------------------------------------------
  wks = gsn_open_wks("png", outdir + "/testplot")

  ;--------------------------------------------------
  ; 4. 配置等值线填色资源
  ;--------------------------------------------------
  res = True
  res@gsnMaximize       = True         ; 最大化画布
  res@cnFillOn          = True         ; 打开填色
  res@cnLinesOn         = False        ; 关闭轮廓线
  res@cnFillPalette     = "BlWhRe"     ; 蓝-白-红配色
  res@tiMainString      = "EP Flux Divergence"
  res@tiMainFontHeightF = 0.02
  res@tiXAxisString     = "Latitude (°)"
  res@tiYAxisString     = "Pressure (hPa)"
  res@trYReverse        = True         ; 压强轴从上到下增大

  ; 如果 pl 单位是 Pa，转换到 hPa 以便画轴标签
  EPdiv_hPa = EPdiv / 100.0

  ;--------------------------------------------------
  ; 5. 绘制填色等值图并输出
  ;--------------------------------------------------
  plot = gsn_contour(wks, EPdiv_hPa, res)
  draw(plot)
  frame(wks)
end