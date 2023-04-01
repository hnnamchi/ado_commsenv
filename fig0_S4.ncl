;*************************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************
begin

; ==============================================================
; ==============================================================
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
;latS   = -40.
;latN   =  40.
latS   = -70.
latN   =  70.
lonL   = -180.
lonR   =  180.
;
yrStrt = 1900
yrLast = 2009
nmos = 12
;
;
plev_up = 400
plev_dn = 1000
;
day2sec = 86400.
print(day2sec)
;
scale_factor = 100000.			;e5
;
p95 = 0.05
lag = 1
mxlag = 6
;
;
yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
 ; print(yyyymm) 
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================

 f      = addfile("/home/hnnamchi/work/datasets/e20c/ITCZ_nsss_e20c.oper.fc.sfc.3hr.128_181_nsss.regn80sc.1900_2010.nc", "r")

  TIME   = f->INITIAL_TIME0_HOURS
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
 
   iecmwf   = dble2flt(f->ITCZ(iYYYY,:))
   iecmwf!0 = "time"
   iecmwf!1 = "lon"
   printVarSummary(iecmwf)                              ; variable overview
   iecmwf = lonFlip(iecmwf)
   printVarSummary(iecmwf)                              ; variable overview

  print(iecmwf&lon)
  
   
   amoc = dim_avg_n_Wrap(iecmwf(time|:,{lon|-30:-20}),1)
  ;  amoc = dim_avg_n_Wrap(iecmwf(time|:,{lon|-30:-10}),1)
  ; amoc = dim_avg_n_Wrap(iecmwf(time|:,{lon|-32:-28}),1)
   printVarSummary(amoc)                              ; variable overview
   delete(iecmwf)
   iecmwf = amoc
   delete(amoc)

  ; print(iecmwf&lon)
; print(xxxxxxxx)

  iecmwf = rmAnnCycle1D(iecmwf)
  iecmwf = dtrend_n(iecmwf,False,0)
  printMinMax(iecmwf,True)


 

delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)

; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hnnamchi/work/geomar/adv/sub0/eof_ts.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

 pc   = f->eof_ts(1,iYYYY)
  printVarSummary(pc)                              ; variable overview

  ;pc = dim_standardize_n(n_psmsl,0,0)      ; normalize


  pc = rmAnnCycle1D(pc)


delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)

;=========================== 
;==============================================================
;===============================================================================================================
;===============================================================================================================
 pc_d = new(dimsizes(pc),typeof(pc),pc@_FillValue)
 x   = pc&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = pc(tt::nmos)
    work = dtrend_msg_n(x1,work,True,False,0)
    pc_d(tt::nmos) = work
 end do

 printVarSummary(pc_d)
 pc = pc_d

;===
delete([/x,x1,pc_d,work,tt/])


;************************************************


; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
  ;f      = addfile ("/home/hyacinth/work/datasets/amoc/Reconstructed_SAMOC_HadSST.nc", "r")
;  f      = addfile ("/home/hyacinth/work/datasets/20cr/ITCZ_20cr_vflx.mon.mean.nc.nc", "r")
;  f      = addfile ("/home/hyacinth/work/datasets/20cr/ITCZ2_20cr_vflx.mon.mean.nc.nc", "r")

;f      = addfile ("/home/hyacinth/work/datasets/20cr/v2c/itcz_maskland_vflx.mon.mean.nc", "r")
 f      = addfile("/home/hnnamchi/work/datasets/20cr/v2c/ITCZ_20cr_v2c_vflx.mon.mean.nc", "r")

  TIME   = f->TIME
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

 ; inoaa   = f->Raw(iYYYY,1)
 ; printVarSummary(inoaa)                              ; variable overview
 
   inoaa   = dble2flt(f->ITCZ(iYYYY,:))
   inoaa!0 = "time"
   inoaa!1 = "lon"
   printVarSummary(inoaa)                              ; variable overview
   inoaa = lonFlip(inoaa)
   printVarSummary(inoaa)                              ; variable overview

  print(inoaa&lon)
   ;xx = itcz(TIME|:,{LON|330:335})
   ;printVarSummary(xx)
   
   amoc = dim_avg_n_Wrap(inoaa(time|:,{lon|-30:-20}),1)
  ;  amoc = dim_avg_n_Wrap(inoaa(time|:,{lon|-30:-10}),1)
  ; amoc = dim_avg_n_Wrap(inoaa(time|:,{lon|-32:-28}),1)
   printVarSummary(amoc)                              ; variable overview
   delete(inoaa)
   inoaa = amoc
   delete(amoc)

  ; print(inoaa&lon)
; print(xxxxxxxx)

  inoaa = rmAnnCycle1D(inoaa)
  inoaa = dtrend_n(inoaa,False,0)
  printMinMax(inoaa,True)

 
;************************************************
; set function arguments
;************************************************
  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 3   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
;************************************************
; calculate spectrum
;************************************************
   secmwf = specx_anal(iecmwf,d,sm,pct)
   decmwf = specx_ci(secmwf, 0.05, 0.95)             ; calc confidence interval
   printVarSummary(decmwf)
;  print(xxxxxxxxxxxxx)
   sig = where(secmwf@spcx.ge.decmwf(3,:)-2.,1/secmwf@frq,secmwf@_FillValue)
   print(sig)
;;;;;;;;;;;;;;;

  
   snoaa = specx_anal(inoaa,d,sm,pct)
   dnoaa = specx_ci(snoaa, 0.05, 0.95)             ; calc confidence interval
   printVarSummary(dnoaa)
  sig2= where(snoaa@spcx.ge.dnoaa(3,:),1/snoaa@frq,snoaa@_FillValue)
  print(sig2)
  
  
;===============================================================
;===============================================================
;==================================================
;
;
;
sm = 7
;
ITCZ = month_to_annual_weighted(yyyymm,iecmwf,1)
ITCZ = dim_standardize_n(ITCZ,0,0)      ; normalize 
printVarSummary(ITCZ)
;delete(iecmwf)
;
PC = month_to_annual_weighted(yyyymm,pc,1)
PC = dim_standardize_n(PC,0,0)      ; normalize 
printVarSummary(PC)
;delete(pc)
;
;
;--------------------------------------------------------------
sip = specxy_anal(ITCZ,PC,d,sm,pct)
d_sip = sip@coher       ;transfers metadata/quadrature
d_sip = where(sip@coher.ge.sip@coher_probability(1),sip@coher,sip@_FillValue)           ;tick on the cospectrum
print(d_sip)

 print("============================================")

;--------------------------------------------------------------
ip = new((/2,dimsizes(sip@quspc)/),float,sip@_FillValue)
printVarSummary(ip)
ip(0,:) = sip@quspc
ip(1,:) = sip@cospc
ip!0 = "lev"
printVarSummary(ip)
;
;


delete(ITCZ)
delete(PC)

ITCZ = iecmwf
ITCZ!0 = "year"
PC = pc
PC!0 = "year"
delete(iecmwf)
delete(pc)
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile ("/home/hnnamchi/work/datasets/e20c/q1q2.mon.e20c.YANAI.apparent_heat.1900.2010.nc", "r")
 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  qu   = dble2flt(f->q1(iYYYY,{plev_up},:,:))
  qd   = dble2flt(f->q1(iYYYY,{plev_dn},:,:))

  q1 = qu
 q1 = qu - qd
 delete(qu)
 delete(qd)

 q1 = q1*day2sec				 ;changes the unit to K/s
 q1 = q1/scale_factor				 ;for plotting
 printVarSummary(q1)                              ; variable overview
;
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
; ==============================================================
; ==============================================================
 ;q1    = lonFlip(q1)
  printVarSummary(q1)                              ; note the longitude coord

  ; ==============================================================
; compute climatology and Anomalies
; ==============================================================
  q1Clm   = clmMonTLL(q1)                             ; (12,lat,lon)
  printVarSummary(q1Clm)

  q1  = calcMonAnomTLL(q1,q1Clm)                ; (time, lat,lon)
  printVarSummary(q1)
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
 q1_d = new(dimsizes(q1),typeof(q1),q1@_FillValue)
 x   = q1&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = q1(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    q1_d(tt::nmos,:,:) = work
 end do

 printVarSummary(q1_d)
 q1 = q1_d

;===
delete([/x,x1,q1_d,work,tt/])
;=========================================================================================
;-----------------------------------------------------------------------------------------
;Q1 = month_to_annual_weighted(yyyymm,q1,1)
;printVarSummary(Q1)

Q1 = q1
Q1!0 = "year"
delete(q1)
;
;==========================================================================================
;==========================================================================================
;==========================================================================================
;==========================================================================================
d_ITCZQ1 = regCoef_n(ITCZ,Q1,0,0)
copy_VarCoords(Q1(0,:,:),d_ITCZQ1)
printVarSummary(d_ITCZQ1)
;==============================================================
;--------------test for significance
 d_ITCZQ1d = d_ITCZQ1
 df = d_ITCZQ1
 lag_ts = esacr(ITCZ(year|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 lag_ts_map   = esacr(Q1(lat|:,lon|:,year|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)
 tval = onedtond(d_ITCZQ1@tval,dimsizes(d_ITCZQ1))
           df(:,:) = floattointeger(dimsizes(Q1(:,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag))))
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)
d_ITCZQ1d = where(abs(prob).lt.p95,d_ITCZQ1d,d_ITCZQ1d@_FillValue) ; Sig. mask
printVarSummary(d_ITCZQ1d)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(ITCZ)
printMinMax(d_ITCZQ1,True)
;
;
;
;==========================================================================================
;==========================================================================================
;==========================================================================================
;==========================================================================================
d_PCQ1 = regCoef_n(PC,Q1,0,0)
copy_VarCoords(Q1(0,:,:),d_PCQ1)
printVarSummary(d_PCQ1)
;==============================================================
;--------------test for significance
 d_PCQ1d = d_PCQ1
 df = d_PCQ1
 lag_ts = esacr(PC(year|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 lag_ts_map   = esacr(Q1(lat|:,lon|:,year|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)
 tval = onedtond(d_PCQ1@tval,dimsizes(d_PCQ1))
           df(:,:) = floattointeger(dimsizes(Q1(:,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag))))
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)
d_PCQ1d = where(abs(prob).lt.p95,d_PCQ1d,d_PCQ1d@_FillValue) ; Sig. mask
printVarSummary(d_ITCZQ1d)
;
delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(Q1)
printMinMax(d_PCQ1,True)
;
;
;
;************************************************
; plotting parameters
;************************************************
   wks  = gsn_open_wks("eps","fig0_S4")             ; send graphics to PNG file
   plot = new(4,graphic)                         ; create graphic array 

   r               = True                        ; plot mods desired
   r@gsnDraw       = False                       ; do not draw
   r@gsnFrame      = False                       ; do not advance frame
   r@tiMainString  = ""                       ; title
   r@tiXAxisString = "log[Period(month]"  ; xaxis
   r@tiYAxisString = "Variance (degrees~S~2~N~ month)"                  ; yaxis
   r@xyLineThicknesses   = (/3.,1.,1.,1./)       ; Define line thicknesses 
   r@xyDashPatterns      = (/0,0,16,16/)           ; Dash patterns 
   r@xyLineColors        = (/"deepskyblue","deepskyblue","deepskyblue","deepskyblue"/)
   ;r@tmYLLabelFontColor = "blue"   ; and right axes are same.
   ;r@tiYAxisFontColor  = "blue"
  ; r@trYLog 		 = True 
   r@trXLog 		 = True 
   r@trXReverse 	= False
    ;r@trXMinF 		= sdof@frq
   ;r@trXMaxF 		=  1.
 ;  r@tmXTLabels         = 1./sdof@frq
    r@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    r@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    r@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    r@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    r@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    r@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    r@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    r@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward

    r@tmXBMode   = "Explicit"
    r@tmXBValues = (/2.,5.,10.,20.,50.,100.,300.,1000./)
    r@tmXBLabels  =  r@tmXBValues
    r@trYMaxF 		=  40.
 
    r@gsnYRefLine         = 0.
    r@tmXBLabelFontHeightF = 0.022
    r@tmYLLabelFontHeightF = 0.022
    r@tiXAxisFontHeightF   = 0.022
    r@tiYAxisFontHeightF   = 0.022
 
;------------------------------------------------
;==========================================================================================
;==========================================================================================
   r0               = True                        ; plot mods desideeppink4
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@trYLog         =  False                    ; log scale
   r0@trXLog        =   True
   r0@trYReverse    = False    ; Reverse X axis values
   r0@trYMinF              =  0.                ; manually set lower limit
   r0@trYMaxF              =  1.                ;   "          upper
   r0@tmXBStyle           = "Log"
   r0@vpHeightF          = 0.40                   ; change aspect ratio of plot
   r0@vpWidthF           = 0.40
   r0@tmYLFormat = "0@*+^sg"
   r0@tmXBFormat = "@*+^sg"
   r0@tiYAxisFontThicknessF  = 1.
   r0@tmXBMajorThicknessF    = 1.
   r0@tmXBMinorThicknessF    = 1.
   r0@tmYLMajorThicknessF    = 1.
   r0@tmYLMinorThicknessF    = 1.
   r0@tmYLLabelFontThicknessF= 1.
   r0@tmBorderThicknessF     = 1.
   r0@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@xyLineThicknesses   = (/1,2,2/)       ; Define line thicknesses 
   r0@xyDashPatterns      = (/0,0,0/)           ; Dash patterns 
   r0@xyLineColors        = (/"seagreen3","seagreen3","deeppink4"/)
   r0@gsnYRefLine         = 0.
   r0@tmXBLabelFontHeightF = 0.015
   r0@tmYLLabelFontHeightF = 0.015
   r0@tiXAxisFontHeightF   = 0.015
   r0@tiYAxisFontHeightF   = 0.015
   r0@tiMainString  = ""                       ; title
   r0@tiXAxisString = ""  ; xaxis
   r0@tiYAxisString = ""                  ; yaxis
   r0@tiMainString           = ""                       ; title
   r0@tiXAxisString = "log [Period (year)]"  ; xaxis
   r0@tmYLLabelFontColor = "seagreen3"   ; and right axes are same.
   r0@tiYAxisFontColor   = "seagreen3"
   r0@tmXTLabelsOn  = False    ; have tick mark labels
   r0@tmXTOn        = True      ; have tick marks
   r0@tiYAxisString         = "Coherence, ~F33~g~S~2~N~" 
   ;r0@tiYAxisString         = "(Coherence)~S~2~N~"
   ;r0@tiYAxisString         = "Covariance (~F33~s~S~2~N~)"
   ;r0@tiYAxisString         = "Cospectrum (normalized)"
   r0@tmXBMode   = "Explicit"

   r0@tmXBValues = (/2.02,3.,5.,8.,20.,50.,100./)
   r0@tmXBLabels  = (/2.,3.,5.,8.,20.,50.,100./)
;===============================================
   r1               = True                        ; plot mods desideeppink4
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@trYLog         =  False                    ; log scale
   r1@trXLog        =  True
   r1@trXReverse         = False     ; Reverse X axis values
   r1@trYMinF              = -180                ; manually set lower limit
   r1@trYMaxF              =  180                ;   "          upper
   r1@tmYUseLeft          = False
   r1@tmYRMode   = "Explicit"
   r1@tmYRValues = ispan(-180,180,90)
   r1@tmYRLabels  =  (/"-~F33~p~","-~F33~p~N~/2","0","~F33~p~N~/2","~F33~p~N~"/)
   ;r1@tmYRValues = ispan(-180,180,45)
   ;r1@tmYRLabels  =  (/"-~F33~p~","-3~F33~p~N~/4","-~F33~p~N~/2","-~F33~p~N~/4","0","~F33~p~N~/4","~F33~p~N~/2","3~F33~p~N~/4","~F33~p~N~"/)
   r1@tmYRLabelFontHeightF = 0.015
   r1@tmYLLabelFontHeightF = 0.015
   r1@tiYAxisString = "Phase lag (rad)"
   r1@tmXTOn        = True      ; have tick marks
   r1@tmYRLabelsOn  = True    ; have tick mark labels
   r1@tmYROn        = True      ; have tick marks
   r1@tmYRLabelFontColor = "deeppink4"   ; and right axes are same.
   r1@tiYAxisFontColor   = "deeppink4"
   r1@tmYLBorderOn = False
   r1@tmYRBorderOn = True
   r1@tmXTBorderOn = True
   r1@tmXBBorderOn = False
   r1@tmYLFormat = "0@*+^sg"
   r1@tiYAxisFontThicknessF  = 1.
   r1@tmXBMajorThicknessF    = 1.
   r1@tmXBMinorThicknessF    = 1.
   r1@tmYRMajorThicknessF    = 1.
   r1@tmYRMinorThicknessF    = 1.
   r1@tmYRLabelFontThicknessF= 1.
   r1@tmBorderThicknessF     = 1.
   r1@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r1@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r1@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r1@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r1@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@xyLineThicknesses   = (/0.5,1,1/)       ; Define line thicknesses 
   r1@xyDashPatterns      = (/0,14,14/)           ; Dash patterns 
   r1@xyLineColors        = (/"deeppink4","black","deeppink4"/)
   r1@tiXAxisFontHeightF   = 0.015
   r1@tiYAxisFontHeightF   = 0.015
   r1@tiMainString  = ""                       ; title
   r1@tiXAxisString = ""  ; xaxis
   r1@tmYRFormat             = "@*+^f"          ;Don't force 0's                                        ; and the left Y-axis.  
 ;===================================================================                              
   plot(0) = gsn_csm_xy(wks,1./secmwf@frq,decmwf,r)

;  r0@tiMainString  = "~F21~ITCZ ~F34~.~F21~ PC2"                           ; title
   plot(1) =gsn_csm_xy2(wks,1./sip@frq,sip@coher,sip@phase,r0,r1);
;-----------------------------------------------------------------------------
;============================================================
; Plot maps
;============================================================
  ; 
  res                      = True
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -0.5            ; set min contour level
  res@cnMaxLevelValF       =  0.5             ; set max contour level
  res@cnLevelSpacingF      =  0.025            ; set contour spacing     
  res@mpFillOn             = True        ; turn on map fill
  res@mpGeophysicalLineThicknessF = 0.1
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR
  res@tmXBLabelFontHeightF = 0.010
  res@tmYLLabelFontHeightF = 0.010
  res@pmTickMarkDisplayMode = "Always"


  res@mpFillDrawOrder      = "preDraw"

  res@gsnAddCyclic         = True
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  ;res@gsnContourNegLineDashPattern  = 1
  ;res@gsnContourZeroLineThicknessF  = 0
  res@lbLabelBarOn         = False      ; turn off individual lb's
  ;res@lbOrientation        = "vertical"
  ;res@lbTitleFontHeightF   = 0.015 
  ;res@lbLabelFontHeightF   = 0.015
  ;res@lbLabelFont          = "Helvetica-Bold"     ; label font
  ;res@lbTitleString        = ""
  ;res@lbTitleAngleF        = 90.
  ;res@lbTitlePosition      = "Right"
  ;res@lbTitleDirection     = "Across"
  ;res@lbRightMarginF       = 0.40
  ;res@lbLeftMarginF        = 0.35

  ;res@mpProjection      = "Robinson"       ; choose projection
  ;res@mpPerimOn         = False             ; turn off box around plot
 ; res@mpGridLatSpacingF = 30.               ; spacing for lat lines
 ; res@mpGridLonSpacingF = 30.               ; spacing for lon lines
 ; res@mpOceanFillColor   = "grey"

  res@lbLabelStride        = 20
;  res@lbTitleString         = "Annual rainfall anomaly (s.d.)"


  res@gsnCenterString         = ""
  res@gsnLeftString           = ""
  res@gsnRightString          = ""

  res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward


;======================================================
;==============================================================
;=========================================================;===============
;===============
;-----------------------------------------------------------------------------
  res1                      = True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
 ; res1@cnMinLevelValF       =  -1.0            ; set min contour level
 ; res1@cnMaxLevelValF       =  1.0             ; set max contour level
  res1@cnLevelSpacingF      = 0.10            ; set contour spacing  
  res1@gsnAddCyclic         = False
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString   = ""
  res1@gsnLeftString   = ""
  res1@cnLineLabelDensityF =  0.
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #4
  opt@gsnShadeHigh = 17                     ; use fill pattern #4
  ;opt@gsnShadeFillScaleF = 2.0
  opt@gsnShadeFillDotSizeF = 0.002
  opt@gsnShadeFillScaleF = 0.9
;-----------------------------------------------------------------------------
;*******************************************
; Create the plots
;*******************************************
  res@gsnLeftStringOrthogonalPosF = -0.95
  res@gsnLeftStringParallelPosF =   0.05
  res@gsnLeftStringFontColor = "grey25"
  res@gsnCenterStringFontHeightF  = 0.013
  res@gsnCenterString  = "Diabatic heating regressed on intertropical convergence zone"

   plotA0 = gsn_csm_contour_map(wks,d_ITCZQ1,res)      
   plotB0 = gsn_csm_contour(wks,d_ITCZQ1d,res1)       
   plotBO = gsn_contour_shade(plotB0,d_ITCZQ1d@_FillValue,p95,opt)
   overlay(plotA0,plotB0)
   plot(2) = plotA0
;
   res@gsnCenterString  = "Diabatic heating regressed on second principal component"
   plotA1 = gsn_csm_contour_map(wks,d_PCQ1,res)
   plotB1 = gsn_csm_contour(wks,d_PCQ1d,res1)       
   plotB1 = gsn_contour_shade(plotB1,d_PCQ1d@_FillValue,p95,opt)
   overlay(plotA1,plotB1)
   plot(3) = plotA1
;----------------------------------------------------------------------
;===============================================
xa = (/120.,180.,180.,120.,120./)
xb = (/0.,0.,120.,120.,0./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0 =   gsn_add_polygon(wks,plot(0),xa,xb,pr0)
d1 =   gsn_add_polygon(wks,plot(1),xa,xb,pr0)

;===============================================================
;======================================================================== 
;===============================================================
;==============================================================


xa = (/10.,15.,15.,10.,10./)
xb = (/-120.,-120.,120.,120.,-120./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0m = new(6,graphic)                ; create graphic array

;---Attach the polygon
     do i = 0,0,1
        d0m(i) = gsn_add_polygon(wks,plot(i+1),xa,xb,pr0)	;mark 10-15 yr period
     end do 
 ;delete(i)
;======================================================================== 
;========================================================================
;========================================================================
  x_0 = 5.
  x_1 = 5.
  y_0 = -30.
  y_1 = 30.

;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 1.0            ; line thickness   
  pres1@gsLineColor      = "seagreen3"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "preDraw"
 ;pres1@gsLineColor      = "black"         ; color of lines
;---------- 

  dum0 = new(2,graphic)
  dum1 = new(2,graphic)

;---Attach the VERTICAL polyline
;     do i = 0,nplots-1,1
;      dum0(i) = gsn_add_polyline(wks,plot(i),(/x_0,x_1/),(/y_0,y_1/),pres1)
;     end do
;        delete(i)
;
  x0 = 2.
  x1 = 100.
  y0 = sip@coher_probability(1)
  ;---Attach the ZERO polyline
     do i = 0,0,1
      dum1(i) = gsn_add_polyline(wks,plot(i+1), (/x0,x1/),(/y0,y0/) ,pres1)
     end do
  delete(i)
;
;===============================================================
;===============================================================
;===============================================================  
  txres               = True
  amres                  = True
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.040   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
  amres@amParallelPosF   =  -0.67  ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =   -0.55    ; This is the center of the plot.
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
   txidc = gsn_create_text(wks, "c", txres);
; Move text string to center/right edge of plot.;
    
  amres@amJust           = "TopLeft"
  annoidc = gsn_add_annotation(plot(0),txidc,amres)  ; Attach string to plot
;===========
  txres@txFontHeightF = 0.025   ; Use same font height as left axis
  txidd = gsn_create_text(wks, "d", txres);
  annoidd = gsn_add_annotation(plot(1),txidd,amres)  ;
;
  amres@amParallelPosF   =  -0.57  ; 0.5 is the right edge of the plot, so

  txres@txFontHeightF = 0.025   ; Use same font height as left axis
  txida = gsn_create_text(wks, "a", txres);
  annoida = gsn_add_annotation(plot(2),txida,amres)  ;
;
  txres@txFontHeightF = 0.025   ; Use same font height as left axis
  txidb = gsn_create_text(wks, "b", txres);
  annoidb = gsn_add_annotation(plot(3),txidb,amres)  ;





;=====================================================
;=====================================================
;***********************************************
   resP             			  = True                        ; panel mods desired
   resP@gsnDraw            		 = True       ; don't draw yet
   resP@gsnFrame            		 = False        ; don't advance frame yet
   resP@gsnMaximize 			 = True                        ; blow up plot
   resP@gsnPanelTop                      = 0.45
   resP@gsnPanelBottom                   = 0.0
   resP@gsnPanelLeft                     = 0.10
   resP@gsnPanelRight                    = 0.40
   resP@gsnPaperOrientation = "portrait"   ; force portrait
  ;resP@gsnPanelFigureStrings   = (/"ERA-20C","NOAA-CIRES"/)
   resP@gsnPanelFigureStringsFontHeightF = 0.010
   resP@gsnPanelFigureStringsJust    =     "TopLeft"
   gsn_panel(wks,plot(0:0),(/1,1/),resP)               ; create panel plots
;-----------------------------------------------------
   res_P                                  = True                        ; panel mods desired
   res_P@gsnDraw              		 = True       ; don't draw yet
   res_P@gsnFrame            		 = False        ; don't advance frame yet
   res_P@gsnMaximize                      = True                        ; blow up plot
   res_P@gsnPanelTop                      = 0.45
   res_P@gsnPanelBottom                   = 0.0
   res_P@gsnPanelLeft                     = 0.47
   res_P@gsnPanelRight                    = 0.85
   res_P@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_P@gsnPanelFigureStrings   = (/"ERA-20C","NOAA-CIRES"/)
   res_P@gsnPanelFigureStringsFontHeightF = 0.010
   res_P@gsnPanelFigureStringsJust    =     "TopLeft"
   gsn_panel(wks,plot(1:1),(/1,1/),res_P) 
;=============================================================================
  resM                                  = True         ; modify the panel plot
  resM@gsnDraw                          = True       ; don't draw yet
  resM@gsnFrame                         = False        ; don't advance frame yet
  resM@gsnPanelLabelBar                 = True        ; add common colorbar
  resM@lbLabelAutoStride                = True         ; auto stride on labels
  resM@lbOrientation                    = "vertical"
  resM@gsnPanelYWhiteSpacePercent       = 2.0    ; default is 1.0
  resM@gsnPanelXWhiteSpacePercent       = 0.0    ; default is 1.0
  resM@gsnPanelTop                      = 0.999
  resM@gsnPanelBottom                   = 0.40
  resM@gsnPanelLeft                     = 0.02
  resM@gsnPanelRight                    = 1.0
  resM@lbTitleFontHeightF               = 0.010
  resM@lbLabelFontHeightF               = 0.010
  resM@lbLabelFont                      = "Helvetica-Bold"     ; label font
  resM@lbTitleString                    = "Diabatic heating (10~S~5~N~ K s~S~-1~N~)"
  ;resM@lbTitleString                    = "Diabatic heating (e~S~5~N~ K s~S~-1~N~)"
  resM@lbOrientation                    = "vertical"
  resM@lbTitleFontHeightF               = .012
  resM@lbLabelFontHeightF               = 0.010
  resM@lbLabelFont                      = "Helvetica-Bold"     ; label font
  resM@lbTitleAngleF                    = 90.
  resM@lbTitleDirection                 = "Across"
  resM@lbTitlePosition                  = "Right"
  resM@lbLeftMarginF                    =  0.1
  resM@lbRightMarginF                   =  0.8
  resM@lbTopMarginF                     = -0.5
  resM@lbBottomMarginF                  = -0.5
  resM@lbLabelStride                    = 20
;  resM@lbTitleOffsetF			= 0.2

  ;resM@lbTitleFontHeightF               = .010
  ;resM@lbLabelFontHeightF               = 0.010
  ;resM@lbLabelFont                      = "Helvetica-Bold"     ; label font
  ;resM@lbTitleString                    = "Diabatic heating (e~S~5~N~K s~S~-1~N~)"
 ; resM@lbTitleAngleF                    = 90.
  ;resM@lbLabelStride                    = 20
  ;resM@lbRightMarginF		        = 0.40
  ;resM@lbLeftMarginF       		= 0.35
  ;resM@lbTitleDirection                 = "Across"
 ;resM@lbTitleString                   = " (s.d.~S~2~N~)Q1 deviations (~S~o~N~C)"
  ;resM@lbTitlePosition                  = "Right"
  resM@gsnMaximize                      = True         ; large format
  resM@gsnPaperOrientation              = "portrait"   ; force portrait
; resM@txString                         = title+": "+yrStrt+"-"+yrLast
; resM@gsnPanelMainFontHeightF          = 0.015
; resM@gsnPanelMainString               = "Annual sea surface height leads inter-tropical convergence zone at negative lags"
;========================================================================
 gsn_panel(wks,(/plot(2:)/),(/2,1/),resM)     ; now draw as one plot
;=============================================================================
;=============================================================================
;=============================================================================
 end

