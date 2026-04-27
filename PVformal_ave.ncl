; Load necessary NCL libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ;-------------------------------
  ; 1. Set up the output parent directory
  ;-------------------------------
  parent_outdir = "/home/weiji/restart_exam/20250307/PV_NCL"
  system("mkdir -p " + parent_outdir)

  ;-------------------------------
  ; 2. Define experiment parameters: years, months, and target pressure levels (in hPa)
  ;-------------------------------
  years = (/"0008", "0013", "0014", "0019"/)
  months = (/"Mar", "Feb"/)
  pv_levels = (/10, 50/)   ; target pressure levels in hPa

  ;-------------------------------
  ; 3. Loop over years, months, and target pressure levels
  ;-------------------------------
  do i = 0, dimsizes(years)-1
    year = years(i)
    do j = 0, dimsizes(months)-1
      month = months(j)
      
      ; Construct file paths for nocouple and couple data based on the month
      if (month .eq. "Feb") then
        noc_file = "/home/weiji/restart_exam/merged_data/" + year + "_Feb_nocoupl_mean.nc"
        coup_file = "/home/weiji/restart_exam/merged_data/" + year + "_Feb_coupl_mean.nc"
      else if (month .eq. "Mar") then
        noc_file = "/home/weiji/restart_exam/merged_data/" + year + "_Mar_nocoupl_mean.nc"
        coup_file = "/home/weiji/restart_exam/merged_data/" + year + "_Mar_coupl_mean.nc"
      else
        print("Invalid month: " + month)
        continue
      end if
      end if

      ; Open the data files
      f_noc = addfile(noc_file, "r")
      f_coup = addfile(coup_file, "r")

      ; Read variables (dimensions: time, plev, lat, lon)
      T_noc = f_noc->T
      U_noc = f_noc->U
      V_noc = f_noc->V

      T_coup = f_coup->T
      U_coup = f_coup->U
      V_coup = f_coup->V

      

      plev  = f_noc->plev    ; Pressure levels in Pa
      lat   = f_noc->lat     ; Latitude (degrees)
      lon   = f_noc->lon     ; Longitude (degrees)

      if (month .eq. "Feb") then
        desiredDays = 120
      else 
        desiredDays = 92
      end if

      ; Determine number of time steps from the first dimension of T_noc
      ;DIM_T = dimsizes(T_noc)
      ;ntime = DIM_T(0)

      ; Loop over target pressure levels
      do k = 0, dimsizes(pv_levels)-1
        target_hPa = pv_levels(k)
        target_Pa  = target_hPa * 100.0  ; Convert hPa to Pa
        
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

        ; Create a subdirectory for the current combination (year, month, target pressure)
        subdir = parent_outdir + "/" + year + "_" + month + "_" + target_hPa + "hPa"
        system("mkdir -p " + subdir)

        ;-----------------------------------------------
        ; 根据 month 指定起始日期（注意添加具体的起始时刻）
        ;-----------------------------------------------
        if (month .eq. "Feb") then
          time_units = "days since 2020-02-01 00:00:00"
        else if (month .eq. "Mar") then
          time_units = "days since 2020-03-01 00:00:00"
        else
          time_units = "days since 2020-01-01 00:00:00"
        end if
        end if

          ;-----------------------------------------------
          ; 预先计算整个时间序列的日期信息
          ;-----------------------------------------------
          ;-----------------------------------------------
          ; 预先计算整个时间序列的日期信息
          ;-----------------------------------------------
          time_indices = tofloat(ispan(0, desiredDays, 1))   ; 转换为浮点型数组
          time_indices@units = time_units   ; 为数组设置 units 属性
          ; 选项 0 返回一个二维数组，每行 [year, month, day, hour, minute, second]
          date_mapping = ut_calendar(time_indices, 0)
          ;print(date_mapping)

          ;-----------------------------------------------
          ; 定义一个月份名称数组，用于将数值月份转换为缩写
          ;-----------------------------------------------
        month_names = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)
  

        ; Loop over each time step ntime cant be used
        do t = 0, desiredDays
          ; Slice the 3D fields at the current time step (dimensions: plev, lat, lon)
          T_noc_slice = T_noc(t,:,:,:)
          U_noc_slice = U_noc(t,:,:,:)
          V_noc_slice = V_noc(t,:,:,:)
          T_coup_slice = T_coup(t,:,:,:)
          U_coup_slice = U_coup(t,:,:,:)
          V_coup_slice = V_coup(t,:,:,:)

          ; Options for PV calculation
          gridType = 1
          opt = 0

          ; Compute full-layer PV for nocouple and couple data using pot_vort_isobaric
          pv_noc = pot_vort_isobaric(plev, U_noc_slice, V_noc_slice, T_noc_slice, lat, gridType, opt)
          pv_coup = pot_vort_isobaric(plev, U_coup_slice, V_coup_slice, T_coup_slice, lat, gridType, opt)
          ; Both pv_noc and pv_coup have dimensions: (plev, lat, lon)

          ; Find the index of the pressure level closest to the target pressure
          diff = abs(plev - target_Pa)
          target_index = ind(diff .eq. min(diff))  ; Use min difference to find closest level
            
          ; Extract the PV field at the target level and convert to PVU (multiply by 1e6)
          PV_noc_target = pv_noc(target_index,:,:)
          PV_coup_target = pv_coup(target_index,:,:)
          PV_noc_target = PV_noc_target * 1e6
          PV_coup_target = PV_coup_target * 1e6
          PV_noc_target@units = "PVU"
          PV_coup_target@units = "PVU"

          ; Get the actual pressure (in hPa) for title display
          pressure_level = plev(target_index)
          pressure_level_hPa = pressure_level / 100.0

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

        ; 根据时间步 t 取出对应的日期信息
          current_year = date_mapping(t, 0)
          current_month_num = toint(date_mapping(t, 1))
          current_day = toint(date_mapping(t, 2))
          ;print(current_day)

             ; 利用 month_names 数组得到月份缩写（注意：月份数值从1开始）
          current_month_str = month_names(current_month_num - 1)

             ; 输出调试信息
          ;print("Debug: t = " + sprintf("%d", t) + " -> " + current_month_str + " " + sprintf("%d", current_day) + ", " + sprintf("%d", current_year))

          ;-----------------------------------------------
          ; 构造图标题（显示动态日期）
          ;-----------------------------------------------
          ps = sprintf("%g", pressure_level_hPa)
  
          title1 = "NoCouple PV (p = " + ps + " hPa, " + current_month_str + "" +current_day+")"
          title2 = "Couple PV (p = " + ps + " hPa, " + current_month_str  + "" +current_day+")"
          res@tiMainString = title1







          ;-------------------------------
          ; Open a workstation for the current frame
          ;-------------------------------
          frame_prefix = subdir + "/" + year + "_" + month + "_" + target_hPa + "hPa_time" +  t
          wks = gsn_open_wks("png", frame_prefix)

          ; Plot the NoCouple PV field
          plot1 = gsn_csm_contour_map(wks, PV_noc_target, res)

          ; Update the title for the Couple plot and create the second plot
          res@tiMainString = title2
          plot2 = gsn_csm_contour_map(wks, PV_coup_target, res)

          ; Panel the two plots side by side
          panelres = True
          panelres@gsnPanelLabelBar = False
          panelres@lbOrientation = "Horizontal"
          panelres@gsnPanelScalePlotIndex = 1
          plots = (/plot1, plot2/)
          gsn_panel(wks, plots, (/1, 2/), panelres)

          ; Save the current frame and delete the workstation to free memory
          delete(wks)
        end do  ; end time loop
        delete(time_indices)
        delete(date_mapping)
      end do  ; end target pressure level loop
      delete(T_noc)
      delete(U_noc)
      delete(V_noc)
      delete(T_coup)
      delete(U_coup)
      delete(V_coup)
      delete(plev)
      delete(lat)
      delete(lon)
    end do  ; end month loop
  end do  ; end year loop
end