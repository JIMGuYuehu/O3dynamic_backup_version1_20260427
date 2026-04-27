;--------------------------------------------------------------
; Extended PV Plotting using NCL (Corrected Version)
; This program loops over 4 years, 2 months, and 2 target pressure levels.
; For each combination, it reads the corresponding "nocouple" and "couple" data,
; computes the full-layer potential vorticity (PV) using pot_vort_isobaric with full data,
; then selects the PV field at the level closest to the target pressure and subsets
; to 60N-90N only for plotting, scales the PV to PVU (×1e6), and plots two contour
; maps side-by-side with a polar stereographic projection.
;
; Each time step is plotted and saved as a PNG image.
; The output images are stored in subdirectories under the parent directory:
;    /home/weiji/restart_exam/20250221/PV_NCL
;--------------------------------------------------------------

; Load necessary NCL libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  ;-------------------------------
  ; 1. Set up the output parent directory
  ;-------------------------------
  parent_outdir = "/home/weiji/restart_exam/20250221/PV_NCL"
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
        noc_file = "/home/weiji/restart_exam/nocouple_data/" + year + "/Feb_NOCOUPL/" + year + "-02/BWCN.e122.f19_g16.002." + year + "-02.001.nc"
        coup_file = "/mnt/backup_ETH/Marina/WACCM/CHEM_2000_restart/BWCN.e122.f19_g16.002_" + year + "/Feb/BWCN.e122.f19_g16.002." + year + "-02.001.nc"
      else if (month .eq. "Mar") then
        noc_file = "/home/weiji/restart_exam/nocouple_data/" + year + "/Mar_NOCOUPL/" + year + "-03/BWCN.e122.f19_g16.002." + year + "-03.001.nc"
        coup_file = "/mnt/backup_ETH/Marina/WACCM/CHEM_2000_restart/BWCN.e122.f19_g16.002_" + year + "/Mar/BWCN.e122.f19_g16.002." + year + "-03.001.nc"
      else
        print("Invalid month: " + month)
        continue
      end if
      end if

      ; Open the data files
      f_noc = addfile(noc_file, "r")
      f_coup = addfile(coup_file, "r")

      ; Read variables (dimensions: time, plev, lat, lon) without immediate subsetting
      T_noc = f_noc->T
      U_noc = f_noc->U
      V_noc = f_noc->V
      T_coup = f_coup->T
      U_coup = f_coup->U
      V_coup = f_coup->V
      plev  = f_noc->plev    ; Pressure levels in Pa
      lat   = f_noc->lat     ; Full latitude range
      lon   = f_noc->lon     ; Full longitude range

      ; Define latitude subset for plotting (60N-90N)
      lat_ind = ind(lat .ge. 60 .and. lat .le. 90)
      if (dimsizes(lat_ind) .lt. 2) then
        print("Error: Not enough latitude points between 60N and 90N.")
        continue
      end if
      lat_subset = lat(lat_ind)

      if (month .eq. "Feb") then
        desiredDays = 120
      else 
        desiredDays = 92
      end if

      ; Debug: Print dimensions of input data
      print("Dimensions of T_noc: " + str_join(tostring(dimsizes(T_noc)), ", "))
      print("Dimensions of plev: " + str_join(tostring(dimsizes(plev)), ", "))
      print("Dimensions of lat: " + str_join(tostring(dimsizes(lat)), ", "))

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

        ; Loop over each time step
        do t = 0, desiredDays-1
          ; Slice the 3D fields at the current time step (dimensions: plev, lat, lon) with full data
          T_noc_slice = T_noc(t, :, :, :)
          U_noc_slice = U_noc(t, :, :, :)
          V_noc_slice = V_noc(t, :, :, :)
          T_coup_slice = T_coup(t, :, :, :)
          U_coup_slice = U_coup(t, :, :, :)
          V_coup_slice = V_coup(t, :, :, :)

          ; Debug: Print dimensions of slices
          print("Time step: " + t + ", Dimensions of T_noc_slice: " + str_join(tostring(dimsizes(T_noc_slice)), ", "))

          ; Options for PV calculation
          gridType = 1
          opt = 0

          ; Compute full-layer PV for nocouple and couple data using pot_vort_isobaric with full data
          pv_noc = pot_vort_isobaric(plev, U_noc_slice, V_noc_slice, T_noc_slice, lat, gridType, opt)
          pv_coup = pot_vort_isobaric(plev, U_coup_slice, V_coup_slice, T_coup_slice, lat, gridType, opt)

          ; Debug: Print dimensions of pv_noc
          print("Dimensions of pv_noc: " + str_join(tostring(dimsizes(pv_noc)), ", "))

          ; Find the index of the pressure level closest to the target pressure
          diff = abs(plev - target_Pa)
          target_index = ind(diff .eq. min(diff))  ; Use min difference to find closest level
          print("Target pressure index: " + target_index)
            
          ; Extract the PV field at the target level and subset to 60N-90N for plotting
          PV_noc_target = pv_noc(target_index, lat_ind, :)
          PV_coup_target = pv_coup(target_index, lat_ind, :)
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
          res@cnFillPalette = "GMT_cool"  ; Blue-Yellow-Red color palette
          res@cnLevelSelectionMode = "ManualLevels"
          res@cnMinLevelValF = cn_min
          res@cnMaxLevelValF = cn_max
          res@cnLevelSpacingF = cn_spacing

          ; Set polar stereographic projection, restrict to 60N-90N
          res@mpProjection = "Stereographic"
          res@mpLimitMode = "LatLon"
          res@mpMinLatF = 60.0
          res@mpMaxLatF = 90.0
          res@mpCenterLatF = 90.0
          res@mpCenterLonF = 0.0  ; Center at 0° longitude
          res@mpPerimOn = True    ; Add perimeter to clearly define the polar region
          res@mpFillOn = False

          ; Create title strings for the two panels
          title1 = "NoCouple PV (p = " + pressure_level_hPa + " hPa, Time " + (t+1) + ", " + target_hPa + " hPa)"
          title2 = "Couple PV (p = " + pressure_level_hPa + " hPa, Time " + (t+1) + ", " + target_hPa + " hPa)"
          res@tiMainString = title1

          ;-------------------------------
          ; Open a workstation for the current frame
          ;-------------------------------
          frame_prefix = subdir + "/" + year + "_" + month + "_" + target_hPa + "hPa_time" + sprinti("%0.3i", t+1)
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

          ; Save the current frame and delete the workstation and temporary variables
          delete(wks)
          delete(pv_noc)
          delete(pv_coup)
          delete(PV_noc_target)
          delete(PV_coup_target)
        end do  ; end time loop
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