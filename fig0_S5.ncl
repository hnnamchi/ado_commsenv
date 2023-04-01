;*************************************************
begin
; ==============================================================
; ==============================================================
;==============================================================================
;============================================================================== 
;************************************************
; set function arguments
;************************************************
  yrStrt = 1900
 ; yrStrt = 1920
  yrLast = 2007
  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 7   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  nmos = 12
  nplots = 1



; ==============================================================
; ==============================================================
f      = addfile ("/home/hnnamchi/work/datasets/gulev/qnet_qe5drer_qh5drer_8007_ssc_btr.dat.nc", "r")
serg   = f->qnet({year|yrStrt:yrLast},season|:,lat|:,lon|:)
printVarSummary(serg)                              ; variable overview
delete(f)

SERG = wgt_areaave_Wrap(serg(year|:,season|3, {lon|-60:-20}, {lat |40:50}),1.0, 1.0, 0)     ;take season of interest
printVarSummary(SERG)
;print(xxxxxxxxxxxxxxxxx)
;==============================================================================
f    = addfile("/home/hnnamchi/work/datasets/SN/SN_y_tot_V2.0.txt_1700_2021.txt.nc", "r")
SSN   = f->ssn({year|yrStrt:yrLast})          		;sunspot number
printVarSummary(SSN)                                   ; variable overview
delete(f) 

 
;==============================================================================
f    = addfile("/home/hnnamchi/work/datasets/nao/HURRELL_nao_pc_monthly_1899_2021.txt.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
nao   = f->nao(iYYYY)          ; "overturning transport " ;
printVarSummary(nao)                                   ; variable overview
;===============================================================================================================
; nao_d = new(dimsizes(nao),typeof(nao),nao@_FillValue)
;x   = nao&time
;x1  = x(::12)
;do tt = 0,nmos-1
;   work = nao(tt::nmos)
;   work = dtrend_msg_n(x1,work,True,False,0)
;   nao_d(tt::nmos) = work
; end do
; printVarSummary(nao_d)
; nao = nao_d

;===
;delete([/x,x1,nao_d,work,tt/])

;===============================================================================================================
;NAO = month_to_annual_weighted(yyyymm,nao,1)   ; Converts to annual-mean
;NAO = month_to_season (nao, "DJF")
;NAO!0 = "year"
;printVarSummary(NAO)
;NAO = dtrend_msg_n(NAO&year,NAO,True,False,0)
;print(NAO)
;print(NAO&year)
print("----------------------------------")
;print(avg(NAO))
;delete(amoc)
delete(TIME)
delete(iYYYY)
delete(YYYY)
delete(f)
;delete(nao)
list_vars()
;print(xxxxxxxxxxxxxxx)
; ==============================================================

;==============================================================================
;n_psmsl = dim_standardize_n(n_psmsl,0,0)      ; normalize 
;============================================================================== 
; f0      = addfile("/home/hnnamchi/work/datasets/20cr/v2c/ITCZ_20cr_v2c_vflx.mon.mean.nc", "r")
 f0      = addfile("/home/hnnamchi/work/datasets/e20c/ITCZ_nsss_e20c.oper.fc.sfc.3hr.128_181_nsss.regn80sc.1900_2010.nc", "r")
  TIME   = f0->INITIAL_TIME0_HOURS   
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_erai    = dble2flt(f0->ITCZ(iYYYY,:))
  
  i_erai!0 = "TIME"
  i_erai!1 = "LON"
  i_erai = lonFlip(i_erai)
  printVarSummary(i_erai) 
  printMinMax(i_erai,True)
   itcz = dim_avg_n_Wrap(i_erai(TIME|:,{LON|-30:-20}),1)
  printVarSummary(itcz)
itcz = rmAnnCycle1D(itcz)

delete(f0)
delete(TIME)
delete(iYYYY)
delete(YYYY)
;delete(i_erai)
;===============================================================================================================
;===============================================================================================================
; itcz_d = new(dimsizes(itcz),typeof(itcz),itcz@_FillValue)
; x   = itcz&TIME
; x1  = x(::12)
; do tt = 0,nmos-1
;    work = itcz(tt::nmos)
;    work = dtrend_msg_n(x1,work,True,False,0)
;    itcz_d(tt::nmos) = work
; end do

; printVarSummary(itcz_d)
; itcz = itcz_d
;===
;delete([/x,x1,itcz_d,work,tt/])  
;************************************************
;itcz = dim_standardize_n(itcz,0,0)      ; normalize 
;itcz = runave(itcz,13,1)
;************************************************
;===============================================================
; ==============================================================
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
; pc_d = new(dimsizes(pc),typeof(pc),pc@_FillValue)
; x   = pc&time
; x1  = x(::12)
; do tt = 0,nmos-1
;    work = pc(tt::nmos)
;    work = dtrend_msg_n(x1,work,True,False,0)
;    pc_d(tt::nmos) = work
; end do
; printVarSummary(pc_d)
; pc = pc_d
;===
;delete([/x,x1,pc_d,work,tt/])  
;===============================================================
;============================================================
;===============================================================
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hnnamchi/work/datasets/e20c/a_RE_qnet_e20c_1900_2010.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  qnet   = f->lh_sh(iYYYY,0,0)
  printVarSummary(qnet)                              ; variable overview
  qnet = -qnet
  ;pc = dim_standardize_n(n_psmsl,0,0)      ; normalize
;  qnet = rmAnnCycle1D(qnet)
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;=========================== 
;==============================================================
;===============================================================================================================
;===============================================================================================================
; qnet_d = new(dimsizes(qnet),typeof(qnet),qnet@_FillValue)
; x   = qnet&time
; x1  = x(::12)
; do tt = 0,nmos-1
;    work = qnet(tt::nmos)
;    work = dtrend_msg_n(x1,work,True,False,0)
;    qnet_d(tt::nmos) = work
; end do
; printVarSummary(qnet_d)
; qnet = qnet_d
;===
;delete([/x,x1,qnet_d,work,tt/])
;===============================================================
;============================================================
;============================================================
;============================================================
;yrStrt = 1901
;;yrStrt = 1921
;yrLast = 2010

;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_eastern_USA_1856_2017.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
n_psmsl = -n_psmsl
printVarSummary(n_psmsl)                              ; variable overview
;=============================================testing/USA
;psmsl   = f->psmsl(iYYYY,23:23)
;psmsl   = f->psmsl(iYYYY,13:)
;psmsl   = f->psmsl(iYYYY,0:14)
;printVarSummary(psmsl)                              ; variable overview
;n_psmsl = dim_avg_n_Wrap(psmsl,1)
;printVarSummary(n_psmsl)
;delete(psmsl)
;===================================================

;n_psmsl = -rmAnnCycle1D(n_psmsl)
;n_psmsl = -n_psmsl 
;printMinMax(n_psmsl,True)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;===============================================================================================================
;===============================================================================================================
; n_psmsl_d = new(dimsizes(n_psmsl),typeof(n_psmsl),n_psmsl@_FillValue)
; x   = n_psmsl&time
; x1  = x(::12)
; do tt = 0,nmos-1
;    work = n_psmsl(tt::nmos)
;    work = dtrend_msg_n(x1,work,True,False,0)
;    n_psmsl_d(tt::nmos) = work
; end do
; printVarSummary(n_psmsl_d)
; n_psmsl = n_psmsl_d
;===
;delete([/x,x1,n_psmsl_d,work,tt/])  
;================================
;************************************************
yyyymm = yyyymm_time(yrStrt, yrLast, "integer")
;print(yyyymm) 
;
;N_PSMSL = month_to_annual_weighted(yyyymm,n_psmsl,1) 
N_PSMSL = month_to_season(n_psmsl,"DJF")
N_PSMSL!0 = "year"
N_PSMSL = dtrend(N_PSMSL,False)
N_PSMSL = dim_standardize_n(N_PSMSL,0,0)      ; normalize 
printVarSummary(N_PSMSL)
;
ITCZ = month_to_annual_weighted(yyyymm,itcz,1) 
;ITCZ = month_to_season(itcz,"MAM")
ITCZ!0 = "year"
ITCZ = dim_standardize_n(ITCZ,0,0)      ; normalize 
printVarSummary(ITCZ)
delete(itcz)
;
PC = month_to_annual_weighted(yyyymm,pc,1) 
;PC = month_to_season(pc,"DJF")
PC!0 = "year"
PC = dim_standardize_n(PC,0,0)      ; normalize 
printVarSummary(PC)
delete(pc)
;
;NAO = month_to_annual_weighted(yyyymm,nao,1)   ; Converts to annual-mean
NAO = month_to_season (nao, "DJF")
NAO!0 = "year"
printVarSummary(NAO)
NAO = dtrend_msg_n(NAO&year,NAO,True,False,0)
print(NAO)
print(NAO&year)
;
;
;QNET = month_to_annual_weighted(yyyymm,qnet,1)   ; Converts to annual-mean
QNET  = month_to_season(qnet, "DJF")
QNET !0 = "year"
printVarSummary(QNET)
QNET  = dtrend_msg_n(QNET&year,QNET,True,False,0)
;
;
;SSN = month_to_annual_weighted(yyyymm,ssn,1)   ; Converts to annual-mean
;SSN = month_to_season(ssn, "DJF")
;SSN!0 = "year"
;printVarSummary(NAO)
SSN = dtrend_msg_n(SSN&year,SSN,True,False,0)
print(SSN)
print(SSN&year)
;
SERG = dtrend_msg_n(SERG&year,SERG,True,False,0)
print("----------------------------------")

;===============================================================
;************************************************
; calculate cospectrum
;************************************************
print("============================================")
 
;--------------------------------------------------------------
shp = specxy_anal(SERG,N_PSMSL,d,sm,pct)
d_shp = shp@coher	;transfers metadata/quadrature
d_shp = where(shp@coher.ge.shp@coher_probability(1),shp@coher,shp@_FillValue)
print(d_shp)

 print("============================================")
 
;--------------------------------------------------------------
hp = new((/2,dimsizes(shp@quspc)/),float,shp@_FillValue)
printVarSummary(hp) 
hp(0,:) = shp@quspc
hp(1,:) = shp@cospc
hp!0 = "lev"
printVarSummary(hp)
;
;===============================================================
;--------------------------------------------------------------
shi = specxy_anal(SERG,PC,d,sm,pct)
d_shi = shi@coher	;transfers metadata/quadrature
d_shi = where(shi@coher.ge.shi@coher_probability(1),shi@coher,shi@_FillValue)
print(d_shi)

print("============================================")
 
;--------------------------------------------------------------
hi = new((/2,dimsizes(shi@quspc)/),float,shi@_FillValue)
printVarSummary(hi) 
hi(0,:) = shi@quspc
hi(1,:) = shi@cospc
hi!0 = "lev"
printVarSummary(hi)
;
;===============================================================
;===============================================================
;--------------------------------------------------------------
sip = specxy_anal(ITCZ,PC,d,sm,pct)
d_sip = sip@coher	;transfers metadata/quadrature
d_sip = where(sip@coher.ge.sip@coher_probability(1),sip@coher,sip@_FillValue)		;tick on the cospectrum
print(d_sip)

 print("============================================")
 
;--------------------------------------------------------------
ip = new((/2,dimsizes(sip@quspc)/),float,sip@_FillValue)
printVarSummary(hi) 
ip(0,:) = sip@quspc
ip(1,:) = sip@cospc
ip!0 = "lev"
printVarSummary(ip)
;
;
delete(pct)
delete(yyyymm)
delete(i_erai)
delete(n_psmsl)
delete(ITCZ)
delete(d)
delete(yrLast)
delete(N_PSMSL)
delete(PC)
delete(hi)
delete(yrStrt)
delete(sm)
;list_vars()
;print(xxxxxxxxxxxxxxxxxxx)
;===============================================================
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
;============================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("eps","fig0_S5")         ; send graphics to PNG file

  plot = new(nplots,graphic)                ; create graphic array ; only needed if paneling

   r0               = True                        ; plot mods desideeppink4
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@trYLog         =  False 			 ; log scale
   r0@trXLog        =  	True
   r0@trYReverse    = False    ; Reverse X axis values
   r0@trYMinF              =  0.                ; manually set lower limit
   r0@trYMaxF              =  1.                ;   "          upper
   r0@tmXBStyle 	  = "Log"  
   r0@vpHeightF          = 0.52                   ; change aspect ratio of plot
   r0@vpWidthF           = 0.42  
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
   r0@gsnYRefLine	  = 0.
   r0@tmXBLabelFontHeightF = 0.020
   r0@tmYLLabelFontHeightF = 0.020
   r0@tiXAxisFontHeightF   = 0.020
   r0@tiYAxisFontHeightF   = 0.020
   r0@tiMainString  = ""                       ; title
   r0@tiXAxisString = ""  ; xaxis
   r0@tiYAxisString = ""                  ; yaxis
   r0@tiMainString 	     = ""                       ; title
   r0@tiXAxisString = "log [Period (year)]"  ; xaxis
   r0@tmYLLabelFontColor = "seagreen3"   ; and right axes are same.
   r0@tiYAxisFontColor   = "seagreen3"
   r0@tmXTLabelsOn  = False    ; have tick mark labels
   r0@tmXTOn        = True      ; have tick marks
   r0@tiYAxisString         = "Coherence, ~F33~g~S~2~N~"
   ;r0@tiYAxisString         = "(Coherence)~S~2~N~"
   ;r0@tiYAxisString 	     = "Covariance (~F33~s~S~2~N~)"
   ;r0@tiYAxisString 	     = "Cospectrum (normalized)"
   r0@tmXBMode   = "Explicit"
   r0@tmXBValues = (/2.02,5.,10.,20.,50.,100./)
   r0@tmXBLabels  = (/2.,5.,10.,20.,50.,100./)

;=======
   r1               = True                        ; plot mods desideeppink4
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@trYLog         =  False 			 ; log scale
   r1@trXLog        =  True
   r1@trXReverse         = False     ; Reverse X axis values
   r1@trYMinF              = -180                ; manually set lower limit
   r1@trYMaxF              =  180                ;   "          upper
   r1@tmYUseLeft	  = False
   r1@tmYRMode   = "Explicit"
   r1@tmYRValues = ispan(-180,180,90)
   r1@tmYRLabels  =  (/"-~F33~p~","-~F33~p~N~/2","0","~F33~p~N~/2","~F33~p~N~"/)
   r1@tmYRMinorValues   = ispan(-180,180,30)
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
   r1@tiXAxisFontHeightF   = 0.020
   r1@tiYAxisFontHeightF   = 0.020
   r1@tiMainFontHeightF    = 0.018
   r1@tiMainString  = ""                       ; title
   r1@tiXAxisString = ""  ; xaxis
   r1@tmYRFormat             = "@*+^f"		;Don't force 0's 					; and the left Y-axis.  
 
;------------------------------------------------------------------------
;************************************************
;------------------------------------------------------------------------  
 
   ;r0@tiMainString  = "~F21~Sea-level index ~F34~.~F21~ PC2"                       ; title

   ;r0@tiMainString  = "Gulev et al reconstruction"                       ; title
   plot(0) =gsn_csm_xy2(wks,1./shp@frq,shp@coher,shp@phase,r0,r1);

   ;r0@tiMainString  = "~F21~Sea-level index ~F34~.~F21~ ITCZ"                          ; title
   r0@tiMainString  = "ECMWF reanalysis"
   ;plot(1) =gsn_csm_xy2(wks,1./shi@frq,shi@coher,shi@phase,r0,r1);

   ;r0@tiMainString  = "~F21~ITCZ ~F34~.~F21~ PC2"                           ; title
;   plot(2) =gsn_csm_xy2(wks,1./sip@frq,sip@coher,sip@phase,r0,r1);
;-----------------------------------------------------------------------------
;=========================================================================
;=========================================================================
;=========================================================================
;************************************************
;=====================================================================
;=========================================================================
;=========================================================================
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise
  txres@txFontHeightF = 0.040   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "DT anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres);
; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.72 ; 0.5 is the right edge of the plot, so
  amres@amOrthogonalPosF =  -0.57    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
 tt = (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"/)
 tt2 = (/"j","k","l","m","n","o","p","q","r"/)
create_tt = gsn_create_text(wks,tt,txres);
create_tt2 = gsn_create_text(wks,tt2,txres);
tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p2 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
;     do i = 0,nplots-1,1
;      tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres)
;     end do
; delete(i)
;
;
 ;===============================================================
;==============================================================


xa = (/10.,15.,15.,10.,10./)
xb = (/-120.,-120.,120.,120.,-120./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0 = new(6,graphic)                ; create graphic array

;---Attach the polygon
     do i = 0,nplots-1,1
      d0(i) = gsn_add_polygon(wks,plot(i),xa,xb,pr0)	;mark 10-15 yr period
     end do
 delete(i)

;--------------
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
  dum0 = new(9,graphic)
  dum1 = new(9,graphic)

;---Attach the VERTICAL polyline
;     do i = 0,nplots-1,1
;      dum0(i) = gsn_add_polyline(wks,plot(i),(/x_0,x_1/),(/y_0,y_1/),pres1)
;     end do
;        delete(i)
;
  x0 = 2.
  x1 = 100.
  y0 = shp@coher_probability(1)
  ;---Attach the ZERO polyline
;     do i = 0,nplots-1,1
;      dum1(i) = gsn_add_polyline(wks,plot(i), (/x0,x1/),(/y0,y0/) ,pres1)
;     end do
;  delete(i)
;
;===============================================================
;===============================================================
xa = (/10.,15.,15.,10.,10./)
xb = (/-120.,-120.,120.,120.,-120./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0 = new(6,graphic)                ; create graphic array

;---Attach the polygon
     do i = 0,nplots-1,1
      d0(i) = gsn_add_polygon(wks,plot(i),xa,xb,pr0)	;mark 10-15 yr period
     end do
 delete(i)

;--------------
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

  dum0 = new(9,graphic)
  dum1 = new(9,graphic)

;---Attach the VERTICAL polyline
;     do i = 0,nplots-1,1
;      dum0(i) = gsn_add_polyline(wks,plot(i),(/x_0,x_1/),(/y_0,y_1/),pres1)
;     end do
;        delete(i)
;
  x0 = 2.
  x1 = 100.
  y0 = shp@coher_probability(1)
  ;---Attach the ZERO polyline
     do i = 0,nplots-1,1
      dum1(i) = gsn_add_polyline(wks,plot(i), (/x0,x1/),(/y0,y0/) ,pres1)
     end do
  delete(i)
;
;=========================================================================
;=========================================================================
;=======================================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True       ; don't draw yet
  resP@gsnFrame                         = True        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = False        ; add common colorbar
  resP@lbLabelAutoStride              	= True         ; auto stride on labels
  ;resP@lbOrientation                    = "vertical"
  resP@gsnPanelYWhiteSpacePercent      	= 2.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 2.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.95
  resP@gsnPanelBottom                   = 0.35
  resP@gsnPanelLeft                     = 0.005
  resP@gsnPanelRight                    = 0.995
  resP@lbTitleFontHeightF	        = 0.008
  resP@lbLabelFontHeightF               = 0.008
  resP@lbLabelFont   		        = "Helvetica-Bold"     ; label font
  resP@lbTitleString                  	= "Sea surface temperature (K)"
  ;resP@lbTitleAngleF                    = 90.
  resP@lbLabelStride       		= 10
  resP@lbTitleDirection                 = "Across"
 ; resP@lbTitleString                  		= " (s.d.~S~2~N~)SST deviations (~S~o~N~C)"
  ;resP@lbTitlePosition                  = "Right"
  ;resP@lbLeftMarginF                    = 0.2
  ;resP@lbRightMarginF                   = 0.2
  resP@lbTopMarginF			= 0.5
 ; resP@lbBottomMarginF			= 0.0
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  ;resP@gsnPanelTop                      = 0.99
  ;resP@gsnPanelBottom                   = 0.55
  ;resP@gsnPanelLeft                     = 0.33
  ;resP@gsnPanelRight                    = 1.00
  ;resP@txString            = title+": "+yrStrt+"-"+yrLast
  ;resP@txString            =	"PC2 leads at positive lags [SST and wind stress lead at negative lags]"
  resP@gsnPanelFigureStringsFontHeightF = 0.007
  resP@gsnPanelFigureStringsJust        = "BottomLeft"
  ;resP@gsnPanelFigureStrings            = (/"lag -9","lag -6","lag -3","lag 0",\
  ;                                              "lag +3","lag +6","lag +9","lag +4","lag +5","lag +6"/)
  resP@gsnPanelCenter   		= False                 ; =True .... must move labels
  ;resP@gsnPanelRowSpec  		= True                  ; tell panel what order to plot
;========================================================================
;gsn_panel(wks,(/plot(pc)/),(/1,6/),resP)     ; now draw as one plot
;gsn_panel(wks,(/plot(pc)/),(/1,7/),resP)     ; now draw as one plot
;gsn_panel(wks,(/plot(0:npanels-1)/),(/5,4/),resP)     ; now draw as one plot
gsn_panel(wks,plot(0:),(/1,1/),resP)     ; now draw as one plot
;;===============================================================
;;=============================================================================
;=======================================================================
;frame(wks) 
end
