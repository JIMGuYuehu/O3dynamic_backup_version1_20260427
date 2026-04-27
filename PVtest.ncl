;--------------------------------------------------------------
; PVtest_new.ncl
; 直接计算全层 PV 并绘制指定气压层数据
; 选取第50个时间点（索引49）的数据，提取第六层数据，
; 并在标题中显示该层对应的实际气压（单位 hPa）
; 对比 nocouple 与 couple 数据
; 输出目录: /home/weiji/restart_exam/20250221/ncl_plot
;--------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ;-------------------------------
  ; 1. 设置输出目录
  ;-------------------------------
  outdir = "/home/weiji/restart_exam/20250221/ncl_plot"
  system("mkdir -p " + outdir)

  ;-------------------------------
  ; 2. 定义数据文件路径（单个文件示例）
  ;-------------------------------
  noc_file  = "/home/weiji/restart_exam/nocouple_data/0008/Feb_NOCOUPL/0008-02/BWCN.e122.f19_g16.002.0008-02.001.nc"
  coup_file = "/mnt/backup_ETH/Marina/WACCM/CHEM_2000_restart/BWCN.e122.f19_g16.002_0008/Feb/BWCN.e122.f19_g16.002.0008-02.001.nc"

  f_noc = addfile(noc_file, "r")
  f_coup = addfile(coup_file, "r")

  ;-------------------------------
  ; 3. 读取变量（数据维度：(time, plev, lat, lon)）
  ;-------------------------------
  T_noc = f_noc->T       ; 温度，单位 K
  U_noc = f_noc->U       ; U 风，单位 m/s
  V_noc = f_noc->V       ; V 风，单位 m/s

  T_coup = f_coup->T
  U_coup = f_coup->U
  V_coup = f_coup->V

  plev  = f_noc->plev    ; 压力层，单位 Pa
  lat   = f_noc->lat     ; 纬度（degrees，须按从南到北排序）
  lon   = f_noc->lon     ; 经度（degrees）

  ;-------------------------------
  ; 4. 选取第50个时间点（索引49）的数据，得到 3D 数据 (plev, lat, lon)
  ;-------------------------------
  time_index = 49
  T_noc_slice = T_noc(time_index,:,:,:)   
  U_noc_slice = U_noc(time_index,:,:,:)
  V_noc_slice = V_noc(time_index,:,:,:)

  T_coup_slice = T_coup(time_index,:,:,:)
  U_coup_slice = U_coup(time_index,:,:,:)
  V_coup_slice = V_coup(time_index,:,:,:)

  ;-------------------------------
  ; 5. 调用 pot_vort_isobaric 计算全层 PV
  ;-------------------------------
  gridType = 1
  opt = 0

  pv_noc = pot_vort_isobaric( plev, U_noc_slice, V_noc_slice, T_noc_slice, lat, (/gridType/), (/opt/) )
  pv_coup = pot_vort_isobaric( plev, U_coup_slice, V_coup_slice, T_coup_slice, lat, (/gridType/), (/opt/) )
  ; pv_noc 和 pv_coup 的维度均为 (plev, lat, lon)

  ;-------------------------------
  ; 6. 提取第六层数据（索引5）用于绘图，并确定该层对应的准确气压层
  ;-------------------------------
  PV_noc_sixth = pv_noc(6,:,:)
  PV_coup_sixth = pv_coup(6,:,:)

  ; 设置为PVU,乘十的六次方，单位也更新为PVU
  PV_noc_sixth = PV_noc_sixth * 1e6
  PV_coup_sixth = PV_coup_sixth * 1e6
  PV_noc_sixth@units = "PVU"  ; 原单位×10⁻⁶，通常称为 PVU（位涡度单位）
  PV_coup_sixth@units = "PVU"


  pressure_level = plev(6)         ; 第六层对应的气压值，单位为 Pa
  pressure_level_hPa = pressure_level / 100.0  ; 转换为 hPa

  ;-------------------------------
  ; 7. 绘图：利用 gsn_csm_contour_map 绘制 PV 场，并面板显示
  ;    使用极地投影（stereographic），中心为北极
  ;-------------------------------
  wks = gsn_open_wks("png", outdir + "/pv_polar_level")

  res = True
  res@cnFillOn      = True
  res@cnLinesOn     = False
  res@cnFillPalette = "BlWhRe"


  ; 设置色标范围
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = 0.0
  res@cnMaxLevelValF       = 1000.0
  res@cnLevelSpacingF      = 50.0  ; 等值线间隔，可根据需要调整

  ; 设置极地投影参数
  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
  res@mpCenterLatF      = 90.               ; choose center lat

  ; 设置标题，显示准确气压层数值（hPa）和时间点（第六层）
  title1 = "NoCouple PV (p = " + sprintf("%4.1f", pressure_level_hPa) + " hPa, Time index 50, Level 6)"
  title2 = "Couple PV (p = " + sprintf("%4.1f", pressure_level_hPa) + " hPa, Time index 50, Level 6)"
  res@tiMainString  = title1
  plot1 = gsn_csm_contour_map(wks, PV_noc_sixth, res)

  res@tiMainString  = title2
  plot2 = gsn_csm_contour_map(wks, PV_coup_sixth, res)

  panelres = True
  panelres@gsnPanelLabelBarOn         = True
  panelres@gsnPanelLabelBarOrientation = "Horizontal"
  panelres@gsnPanelFigureWidth        = 0.9
  panelres@gsnPanelFigureHeight       = 0.5

  plots = (/ plot1, plot2 /)
  gsn_panel(wks, plots, (/1,2/), panelres)

end
