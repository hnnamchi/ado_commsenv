
begin
; ==============================================================
; ==============================================================
    yrStrt = 1920
    yrLast = 2009
;   yrStrt = 1900
;   yrLast = 1985

   nmos = 12
   yyyymm = yyyymm_time(yrStrt, yrLast, "integer")
   print(yyyymm) 
 
   xlabel = ispan(yrStrt,yrLast,1)
   printVarSummary(xlabel)
 ;  yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
  ; print(yyyymm) 
   
   ;xlabel = ispan(yrStrt,yrLast,1) 
    

  ee = exp(log(100000.))
  print(ee)		;for log scaling
  ;print(xxxxxxxx) 
; ==============================================================
; ==============================================================
; Open & read mask data
; ==============================================================
;===============================================================================================================
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_eastern_USA_1856_2017.nc", "r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview

s_psmsl   = f->s_psmsl(iYYYY)
printVarSummary(s_psmsl)                              ; variable overview


delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(f)
;===============================================================================================================

s_psmsl = s_psmsl - n_psmsl
n_psmsl =  -n_psmsl
;=========================================================================================
;=========================================================================================                         
;=========================================================================================
;=========================================================================================                         
NN = month_to_annual_weighted(yyyymm,n_psmsl,1)
;NN = month_to_season(n_psmsl,"DJF")
NN!0 = "year"
NN = dtrend_msg_n(NN&year,NN,False,False,0)
NN = dim_standardize_n(NN, 1, 0)  ; new variable
delete(n_psmsl)
;
SS = month_to_annual_weighted(yyyymm,s_psmsl,1)
;SS = month_to_season(s_psmsl,"DJF")
SS!0 = "year"
SS = dtrend_msg_n(SS&year,SS,False,False,0)
SS = dim_standardize_n(SS, 1, 0)  ; new variable
delete(s_psmsl)
print("----------------------------------")
;======================================================================================================
;===============================================================================================================
rr = escorc(SS,NN)
print(rr)
prob   = rtest(rr,dimsizes(NN),0)
print(prob)
;print(xxxxxxxxx)
;============================================================
;------------------------------------------------------------------
all_PSMSL =  new ((/2,dimsizes(NN)/),float)
all_PSMSL(0,:) = NN		;AMOC index
all_PSMSL(1,:) = SS             ;sea level index
all_PSMSL!0 = "lev"
printVarSummary(all_PSMSL)
delete(NN)
print(SS)
;===============================================================================================================

  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 3   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  
;===========================  
; ==============================================================
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_eastern_USA_1856_2017.nc", "r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview

s_psmsl   = f->s_psmsl(iYYYY)
printVarSummary(s_psmsl)                              ; variable overview
;
ns  = n_psmsl 
ns = s_psmsl - n_psmsl
;
n_psmsl = rmAnnCycle1D(n_psmsl)
;
delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(f)
;===============================================================================================================
;===============================================================================================================
 n_psmsl_d = new(dimsizes(n_psmsl),typeof(n_psmsl),n_psmsl@_FillValue)
 x   = n_psmsl&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = n_psmsl(tt::nmos)
    work = dtrend_msg_n(x1,work,True,False,0)
    n_psmsl_d(tt::nmos) = work
 end do

 printVarSummary(n_psmsl_d)
 n_psmsl = n_psmsl_d

;===
delete([/x,x1,n_psmsl_d,work,tt/])  
n_psmsl = dim_standardize_n(n_psmsl,0, 0)  ; new variable 
;================================
;===============================================================================================================
;===============================================================================================================
 ns_d = new(dimsizes(ns),typeof(ns),ns@_FillValue)
 x   = ns&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = ns(tt::nmos)
    work = dtrend_msg_n(x1,work,True,False,0)
    ns_d(tt::nmos) = work
 end do

 printVarSummary(ns_d)
 ns = ns_d

;===
delete([/x,x1,ns_d,work,tt/])
ns = dim_standardize_n(ns,0, 0)  ; new variable 
;************************************************
; calculate spectrum
;************************************************
  sMSL = specx_anal(n_psmsl,d,sm,pct)
  dMSL = specx_ci(sMSL,0.05, 0.95)             ; calc confidence interval
  printVarSummary(sMSL)
  printVarSummary(dMSL)
 
  sMSLm = specx_anal(ns,d,sm,pct)
  dMSLm = specx_ci(sMSLm,0.05, 0.95)             ; calc confidence interval
  printVarSummary(sMSLm)
  printVarSummary(dMSLm)
;==============================================================
;============================================================
  wks = gsn_open_wks("eps","fig0_S2")         ; send graphics to PNG file

  plot = new(2,graphic)                ; create graphic array ; only needed if paneling
;====================================================
;====================================================
   r               = True                        ; plot mods desired
   r@gsnDraw       = False                       ; do not draw
   r@gsnFrame      = False                       ; do not advance frame
   r@tiMainString  = ""                       ; title
   r@tiXAxisString = "log [Period (month)]"  ; xaxis
   r@tiYAxisString = "Variance [Monthly sea-level (s.d.~S~2~N~~F21~ month)]"                  ; yaxis
   ;r@tiYAxisString = "Variance [Sea-level index (~F33~s~S~2~N~~F21~ month)]"                  ; yaxis
   ;r@tiYAxisString = "Variance [Sea level (e~S~5~N~ mm~S~2~N~ month)]"                  ; yaxis
   r@xyLineThicknesses   = (/2.0,0.5,0.5,0.5/)       ; Define line thicknesses 
   r@xyDashPatterns      = (/0,0,16,16/)           ; Dash patterns 
   r@xyLineColors        = (/"blue","blue","blue","blue"/)
   ;r@tmYLLabelFontColor = "blue"   ; and right axes are same.
   ;r@tiYAxisFontColor  = "blue"
   r@vpHeightF		   = 0.45                   ; change aspect ratio of plot
   r@vpWidthF		   = 0.45         
   ; r@trYLog 		 = True 
   r@trXLog 		 = True 
   r@trXReverse 	= False
    ;r@trXMinF 		= min(1/sMSL@frq)
   ; r@trXMaxF 		= max(1/sMSL@frq)
    r@trYMinF 		= 0.
    r@trYMaxF 		= 40
    ;r@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    ;r@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    ;r@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    ;r@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    ;r@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    ;r@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    ;r@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    ;r@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
    r@tiYAxisFontThicknessF  = 1.
    r@tmXBMajorThicknessF    = 1.
    r@tmXBMinorThicknessF    = 1.
    r@tmYLMajorThicknessF    = 1.
    r@tmYLMinorThicknessF    = 1.
    r@tmYLLabelFontThicknessF= 1.
    r@tmBorderThicknessF     = 1.
    ;r@tmXTBorderOn            =  False
    ;r@tmYRBorderOn            =  False
    ;r@tmXTOn                  =  False
    ;r@tmYROn                  =  False
    r@tmXBMode   = "Explicit"
    r@tmXBValues = (/2.,12.,30.,60.,120.,300.,900./)
    r@tmXBLabels  =  r@tmXBValues
   ; r@tmXTLabelsOn = True
    ;r@tmXUseBottom  = False
    r@tmXTMode   = "Explicit"
    r@tmXTValues = r@tmXBValues
    ;r@tmXTLabels  =  r@tmXTValues/12.
    r@tiYAxisFontHeightF = 0.018
    r@tiXAxisFontHeightF = 0.018
    r@tmXBLabelFontHeightF = 0.018
    r@tmYLLabelFontHeightF = 0.018
;------------------------------------------------
    plotA0 = gsn_csm_xy(wks,1./sMSL@frq,dMSL,r)
    r@xyLineColors        = (/"red","red","red","red"/)
    plotB0 = gsn_csm_xy(wks,1./sMSLm@frq,dMSLm,r)
    overlay(plotA0,plotB0)
    plot(0) = plotA0
;====================================================
  ; l               = True                        ; plot mods desired
  ; l@gsnDraw       = False                       ; do not draw
  ; l@gsnFrame      = False                       ; do not advance frame
  ; l@tiMainString  = ""                       ; title
  ; l@tiXAxisString = "Time (year)"  ; xaxis
  ; l@tiYAxisString = "Sea level anomaly (mm)"                  ; yaxis
  ; l@tiYAxisFontHeightF = 0.008
  ; l@tiXAxisFontHeightF = 0.008
  ; l@xyLineThicknesses   = (/1.,2.,2./)       ; Define line thicknesses 
  ; l@xyDashPatterns      = (/0,0,0/)           ; Dash patterns 
  ; l@xyLineColors        = (/"grey50","red","deepskyblue"/)
  ; l@vpHeightF		 = 0.20                   ; change aspect ratio of plot
  ; l@vpWidthF		 = 0.40           
  ; l@tmXBLabelFontHeightF = 0.008
  ; l@tmYRLabelFontHeightF = 0.008
  ; l@tmYLFormat 	 = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; l@trYMinF 		 = -60.
  ; l@trYMaxF 		 =  60.
  ; l@trXMinF 		 =  1900.
  ; l@trXMaxF 		 =  2020.
  ; l@gsnYRefLine         = (/0.0/)             ; create a reference line
  ; l@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmXBMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmYLMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmXTMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; l@tmYRMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  
;************************************************ 
;------------------------------------------------
 ;plot(1) = gsn_csm_xy(wks,xlabel,all_PSMSL,l)
;-----------------------------------------------------------------------------
   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.020
   r2@tiYAxisFontHeightF = 0.020
   r2@tiXAxisString = "Time (year)"  ; xaxis
   r2@tiYAxisString = "Annual indices (s.d.)"                  ; yaxis
   ;r2@tiYAxisString = "Sea-level index (~F33~s~F21~)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.020
   r2@tmYLLabelFontHeightF = 0.020
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   r2@trXMinF             = yrStrt -5.                 ; manually set lower limit
   r2@trXMaxF             = yrLast +5.                ;   "          upper
   r2@trYMinF             = -3.0                 ; manually set lower limit
   r2@trYMaxF             =  3.0                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/1./)
   r2@xyLineColors        = (/"blue","red","transparent"/)
   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   r2@gsnYRefLine           =  1900. 
   
   r2@vpWidthF                = 0.80
   r2@vpHeightF               = 0.40
   ;r2@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r2@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r2@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r2@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r2@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r2@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r2@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r2@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward

   r2@tiYAxisFontThicknessF  = 1.
   r2@tmXBMajorThicknessF    = 1.
   r2@tmXBMinorThicknessF    = 1.
   r2@tmYLMajorThicknessF    = 1.
   r2@tmYLMinorThicknessF    = 1.
   r2@tmYLLabelFontThicknessF= 1.
   r2@tmBorderThicknessF     = 1.
 
   r2@xyLineThicknesses   = (/1.,1.,1./)       ; Define line thicknesses 
   r2@xyDashPatterns      = (/1.,0.,0./)           ; Dash patterns 
   ;r2@tmYLLabelFontColor = "blue"   ; and right axes are same.
   ;r2@tiYAxisFontColor  = "blue"
   ;r2@tmBorderLineColor = "deepskyblue"
   ;r2@tmYLMajorLineColor = "red" ; YR axis uses YL resources, can
   ;r2@tmYLMinorLineColor = "red" ; change by setting tmYUseLeft = False 
   ;r2@tmBorderLineColor =  "red"

   ;===legends
   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   r2@pmLegendSide           = "Top"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.27                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.40                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.08                ; height of legend.
   r2@lgLabelFontHeightF     = 0.016                ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = ""           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"This study","McCarthy et al. (2015)","ORAS4"/)         ; explicit labels 
;===============================================================
;=============================================================== 
;***********************************************
   r1               = True                        ; plot mods desired
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@tiXAxisString = ""  ; xaxis
   r1@tiYAxisString = "RAPID-MOCHA AMOC (Sv)"      
   r1@tiYAxisAngleF  =  270. 
   r1@xyLineColors   =  (/"red","red","transparent"/)
   r1@xyLineThicknesses  =  (/1.,1.,1./)
   r1@xyDashPatterns     =  (/1.,0.,0./)
   r1@tiXAxisFontHeightF = 0.020
   r1@tiYAxisFontHeightF = 0.020
   r1@tmXBLabelFontHeightF = 0.020
   r1@tmYLLabelFontHeightF = 0.020
       
   r1@tmYUseLeft          =  True
   ;r1@trYMinF             =  13.7                 ; manually set lower limit
   ;r1@trYMaxF             =  20.0                 ;   "          upper
   r1@trYMinF             =  -3.0                 ; manually set lower limit
   r1@trYMaxF             =  3.0                 ;   "          upper


   r1@trXMinF             = 1900.                 ; manually set lower limit
   r1@trXMaxF             = 2020.                ;   "          upper
   
     
   r1@vpWidthF                = 0.80
   r1@vpHeightF               = 0.40
   r1@tmXTBorderOn            =  True
   r1@tmXTOn                  =  True
    
   
  r1@tmYRLabelFontColor = "red"   ; and right axes are same.
  r1@tiYAxisFontColor  = "red"
  ;r1@tiYAxisString   = "SSH anomaly (cm)"
  
  r1@tmXTLabelsOn        = False    ; have tick mark labels
  r1@tmXTOn              = True      ; have tick marks
  ;r1@tmYRLabelFontColor = "deepskyblue"   ; and right axes are same.
  ;r1@tiYAxisFontColor  = "deepskyblue"
 ; r1@tmBorderLineColor =  "deepskyblue"
  r1@tmYLBorderOn 	 = False
  r1@tmYRBorderOn	 = True
  r1@tmXTBorderOn	 = False
  r1@tmXBBorderOn	 = False
  ;r1@tmBorderLineColor = "deepskyblue"
  ;r1@tmYRMajorLineColor = "deepskyblue" ; YR axis uses YL r1ources, can
  ;r1@tmYRMinorLineColor = "deepskyblue" ; change by setting tmYUseLeft = False 
  ;r1@gsnYRefLine           = avg(RAPID)
  r1@gsnYRefLineColor     = "red"
  r1@tmYRFormat   	=  "@*+^sg"          ; trim the uneccecary "0" in labels. 
  
  r1@xyMarker	        = 	16     
  r1@xyMarkLineMode	=	"MarkLines"
  r1@xyMarkerColor	=	"red"
  r1@xyMarkerSizeF	=	0.008
   ;r1@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r1@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r1@tiYAxisFontThicknessF  = 1.
   r1@tmXTMajorThicknessF    = 1.
   r1@tmXTMinorThicknessF    = 1.
   r1@tmYRMajorThicknessF    = 1.
   r1@tmYRMinorThicknessF    = 1.
   r1@tmYRLabelFontThicknessF= 1.
   r1@tmBorderThicknessF     = 1.
   
   
   plot(1)  =  gsn_csm_xy(wks,xlabel,all_PSMSL(lev|:,year|:),r2)

    

   ;yy_data = new(dimsizes(N_PSMSL),typeof(nao),nao@_FillValue)

;   plot(1)  =  gsn_csm_x2y2(wks,xlabel-1,xlabel,all_PSMSL,RAPID,r2,r1)
 
   ;plot(1)  =  gsn_csm_x2y2(wks,xlabel-1,xlabel,N_PSMSL,RAPID,r2,r1)
   ;print(xlabel-1)
   
   
xa = (/120.,180.,180.,120.,120./)
xb = (/0.,0.,120000.,120000.,0./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0 =   gsn_add_polygon(wks,plot(0),xa,xb,pr0)


;pr0@gsFillColor =  "lightpink"
xa = (/1996.,1999.,1999.,1996.,1996./)
xb = (/-1200.,-1200.,1200.,1200.,-1200./)
;d1 = gsn_add_polygon(wks,plot(1),xa,xb,pr0)
;
xa = (/2009.,2010.,2010.,2009.,2009./)
;d2 = gsn_add_polygon(wks,plot(1),xa,xb,pr0)
;===============================================================
;=====================================================
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.035   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
  txid0 = gsn_create_text(wks, "b", txres);
; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.80  ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  annoid0 = gsn_add_annotation(plot(0),txid0,amres)  ; Attach string to plot
;===========
 txres@txFontHeightF = 0.028   ; Use same font height as left axis
 amres@amParallelPosF   =  -0.65  ; 0.5 is the right edge of the plot, so
 txid1 = gsn_create_text(wks, "a", txres);
 annoid1 = gsn_add_annotation(plot(1),txid1,amres)  ;

;----------------------------------------------------------------------
; Code for customizing a legend via simple_legend.
;----------------------------------------------------------------------
;======================================================== 
;========================================================
 var_names1 = (/"This study","McCarthy et al. (2015)"/)
print(var_names1)
;
i1_col = (/"blue","red","black","red","mediumpurple3"/)
i2_col = (/"red","red"/)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.009
  txres@txJust        = "CenterLeft"
  txres@txFont                = "Helvetica"
  txres@txAngleF      =  0.

  labels  = var_names1
  ypos    = fspan(0.93, 0.95,dimsizes(var_names1))

 do i = 0,dimsizes(labels)-1,1
      txres@txFontColor = i1_col(i)
     gsn_text_ndc(wks,labels(i),0.72,ypos(i),txres)
  end do
delete(i)

;---------------------------------------------------------------------- 
;----------------------------------------------------------------------
; Code for customizing a legend via simple_legend.
;----------------------------------------------------------------------
;======================================================== 
;========================================================
var_names1 = (/"~F10~r~F21~  = 0.17"/)
print(var_names1)
;
var_names2 = (/"1921-2010","1941-2010"/)
print(var_names2)
;
i1_col = (/"blue","deepskyblue","black","red","mediumpurple3"/)
i2_col = (/"red","red"/)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.010
  txres@txJust        = "CenterLeft"
  txres@txFont                = "Helvetica"
  txres@txAngleF      =  0.

  labels  = var_names1
 ; ypos    = fspan(0.87, 0.92,dimsizes(var_names1))

 do i = 0,dimsizes(labels)-1,1
      txres@txFontColor =  "black"
     gsn_text_ndc(wks,labels(i),0.50,0.92,txres)
  end do
delete(i)
;labels  = var_names2
; do i = 0,dimsizes(labels)-1,1
;      txres1@txFontColor = i2_col(i)
;     gsn_text_ndc(wks,labels(i),xpos(i),0.36,txres1)
;  end do
;delete(i)
;delete(txres1)
;
;

;===============================================================
;==============================================================
;=====================================================================
  
;===============================================================
;***********************************************
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelTop                      = 0.99
  res_P@gsnPanelBottom                   = 0.65
  res_P@gsnPanelLeft                     = 0.60
  res_P@gsnPanelRight                    = 0.99
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_P@gsnPanelFigureStrings   = (/"~F0~Peak~F34~;~F0~74 yr~S~-1~N~","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  ;res_P@gsnPanelFigureStringsFontHeightF = 0.015
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(0:0),(/1,1/),res_P)     ; now draw as one plot
;======================================================
;======================================================

  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True       ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	           	= False        ; add common colorbar
  resP@lbLabelAutoStride              		= True         ; auto stride on labels
  resP@lbOrientation                    = "vertical"
  resP@gsnPanelYWhiteSpacePercent      	= 0.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 1.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.99
  resP@gsnPanelBottom                   = 0.65
  resP@gsnPanelLeft                     = 0.01
  resP@gsnPanelRight                    = 0.65
  resP@lbTitleFontHeightF	        = .016 
  resP@lbLabelFontHeightF               = 0.014
  resP@lbLabelFont   		        = "Helvetica-Bold"     ; label font
  resP@lbTitleString                  	= "SST anomaly (s.d.)"
  resP@lbTitleAngleF                    = 90.
  resP@lbLabelStride       		= 5
  resP@lbTitleDirection                 = "Across"
 ; resP@lbTitleString                  		= " (s.d.~S~2~N~)SST deviations (~S~o~N~C)"
  resP@lbTitlePosition                  = "Right"
  resP@lbLeftMarginF                    = 0.5
  resP@lbRightMarginF                   = 0.5
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait

   resP@gsnLeftString  = ""
 ; resP@txString            = title+": "+yrStrt+"-"+yrLast
;========================================================================
 gsn_panel(wks,(/plot(1:1)/),(/1,1/),resP)     ; now draw as one plot
;gsn_panel(wks,plotMAP(0:1),(/1,2/),resP)     ; now draw as one plot
;===============================================================
;========================================================================
;===============================================================
frame(wks) 
end
