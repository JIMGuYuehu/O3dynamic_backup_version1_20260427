;--------------------------------------------------------------
; PV Plotting for ERA5 Reference Data using NCL
; This program reads PV data from a reference NetCDF file,
; plots the PV field at 50 hPa and 10 hPa pressure levels across all time steps,
; scales the PV to PVU (×1e6), and generates contour maps with a polar projection.
;
; Output images are saved as PNG files in:
;    /home/weiji/restart_exam/20250221/PV_NCL/ERA52020reference
;--------------------------------------------------------------

; Load necessary NCL libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  ;-------------------------------
  ; 1. Set up the output directory
  ;-------------------------------
  parent_outdir = "/home/weiji/restart_exam/20250221/PV_NCL/ERA52020reference"
  system("mkdir -p " + parent_outdir)

  ;-------------------------------
  ; 2. Define the input file and pressure levels
  ;-------------------------------
  ref_file = "/home/weiji/restart_exam/9bf72024082375cc81b7c8ef32698c98.nc"
  pv_levels = (/50, 10/)   ; Pressure levels available in the reference data (hPa)

  ; Open the reference data file
  f_ref = addfile(ref_file, "r")

  ; Read the PV variable (dimensions: valid_time, pressure_level, latitude, longitude)
  pv_ref = f_ref->pv

  ; Read coordinate variables
  valid_time = f_ref->valid_time
  pressure_level = f_ref->pressure_level  ; in hPa (50 and 10)
  lat = f_ref->latitude                   ; 90N to 60N
  lon = f_ref->longitude                  ; -180 to 179.75

  ; Number of time steps
  ntime = dimsizes(valid_time)  ; 121 time steps

  ;-------------------------------
  ; 3. Loop over pressure levels and time steps
  ;-------------------------------
  do k = 0, dimsizes(pv_levels)-1
    target_hPa = pv_levels(k)
    
    ; Set contour parameters based on target pressure level
    if (target_hPa .eq. 10) then
      cn_min = 0.0
      cn_max = 1000.0
      cn_spacing = 50.0
    else if (target_hPa .eq. 50) then
      cn_min = 0.0
      cn_max = 100.0
      cn_spacing = 5.0
    else
      cn_min = 0.0
      cn_max = 1000.0
      cn_spacing = 50.0
    end if
    end if
    

    ; Create a subdirectory for the current pressure level
    subdir = parent_outdir + "/" + target_hPa + "hPa"
    system("mkdir -p " + subdir)

    ; Find the index of the target pressure level
    diff = abs(pressure_level - target_hPa)
    target_index = ind(diff .eq. min(diff))  ; Closest level index

    time_units = "days since 2020-02-01 00:00:00"
    time_indices = tofloat(ispan(0, 121, 1))   ; 转换为浮点型数组
    time_indices@units = time_units   ; 为数组设置 units 属性
    ; 选项 0 返回一个二维数组，每行 [year, month, day, hour, minute, second]
    date_mapping = ut_calendar(time_indices, 0)
    ;print(date_mapping)

    ;-----------------------------------------------
    ; 定义一个月份名称数组，用于将数值月份转换为缩写
    ;-----------------------------------------------
    month_names = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)   



    ; Loop over each time step
    do t = 0, ntime-1
      ; Extract the PV field at the target level and current time step
      PV_target = pv_ref(t, target_index, :, :)
      
      ; Scale PV to PVU (multiply by 1e6, assuming original units are K m**2 kg**-1 s**-1)
      PV_target = PV_target * 1e6
      PV_target@units = "PVU"

      ;-------------------------------
      ; Set up plotting resources
      ;-------------------------------
      res = True
      res@cnFillOn = True
      res@cnLinesOn = False
      res@cnFillPalette = "MPL_BuGn"
      res@cnLevelSelectionMode = "ManualLevels"
      res@cnMinLevelValF = cn_min
      res@cnMaxLevelValF = cn_max
      res@cnLevelSpacingF = cn_spacing

      ; Set polar projection parameters
      res@mpProjection = "Orthographic"
      res@mpPerimOn = False
      res@mpFillOn = False
      res@mpCenterLatF = 90.0

      ; Create title string

      current_year = date_mapping(t, 0)
      current_month_num = toint(date_mapping(t, 1))
      current_day = toint(date_mapping(t, 2))
      current_month_str = month_names(current_month_num - 1)
      
      title = "ERA5 PV (p = " + target_hPa + " hPa, " + current_month_str + "" +current_day+")"
      res@tiMainString = title

      ;-------------------------------
      ; Open a workstation for the current frame
      ;-------------------------------
      frame_prefix = subdir + "/ERA5_" + target_hPa + "hPa_time" + sprinti("%0.3i", t+1)
      wks = gsn_open_wks("png", frame_prefix)

      ; Plot the PV field
      plot = gsn_csm_contour_map(wks, PV_target, res)

      ; Close the workstation to save the image and free memory
      delete(wks)
    end do  ; end time loop
  end do  ; end pressure level loop

  ; Clean up
  delete(f_ref)
  delete(pv_ref)
  delete(valid_time)
  delete(pressure_level)
  delete(lat)
  delete(lon)
end