
begin
; ==============================================================
; ==============================================================
 ;yrStrt = 2010; up to 2010
 yrStrt = 2005	
 yrLast = 2019
 ;yrLast = 2014		;start from 2014

 n1 = (yrLast - yrStrt)+1
 print(n1)
 ;print(xxxxxxxxxxxxxx)
 yyyymm = yyyymm_time(yrStrt, yrLast, "integer")
 print(yyyymm) 
 
nmos = 12
mxlag    = 4

 xlabel2 = ispan(yrStrt,yrLast,1)
 printVarSummary(xlabel2)

; ==============================================================
f    = addfile("/home/hnnamchi/work/datasets/amoc/moc_transports_2020.nc", "r")
;f    = addfile("/home/hyacinth/work/datasets/amoc/moc_transports_dec2019.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
AMOC   = dble2flt(f->moc_mar_hc10(iYYYY))		; "overturning transport " ;
printVarSummary(AMOC)                           	; variable overview
amoc = calculate_monthly_values(AMOC, "avg", 0, True)
printVarSummary(amoc) 
printVarSummary(amoc)                              ; variable overview
print(amoc(0:11))
delete(AMOC)
amoc = rmAnnCycle1D(amoc)
;===============================================================================================================
;==============================================================
yyyymmA = cd_calendar(amoc&time,-1)
yrfracA = yyyymm_to_yyyyfrac(yyyymmA, 0.0); not used here
;=============================================================
;===============================================================================================================
 amoc_d = new(dimsizes(amoc),typeof(amoc),amoc@_FillValue)
 x   = amoc&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = amoc(tt::nmos)
    work = dtrend_msg_n(x1,work,False,False,0)
    amoc_d(tt::nmos) = work
 end do
 printVarSummary(amoc_d)
 amoc = amoc_d
;===
delete([/x,x1,amoc_d,work,tt/])  

amoc = dim_standardize_n(amoc,0, 0)  ; normalise
;
AMOC = month_to_season(amoc,"DJF")
print(AMOC)
;
;delete(amoc)
delete(TIME)
delete(iYYYY)
delete(YYYY)
delete(f)
list_vars()
; ==============================================================
; Open & read mask data
; ==============================================================
;===============================================================================================================
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_eastern_USA_1856_2017.nc", "r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview
n_psmsl = rmAnnCycle1D(n_psmsl)
n_psmsl = -n_psmsl 
;printMinMax(n_psmsl,True)

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
;==================================
    print(amoc)
    print(n_psmsl)
    Damoc = dim_stddev_n(amoc,0)
    Dn_psmsl = dim_stddev_n(n_psmsl,0)
    print(Damoc)
    print(Dn_psmsl)
    ;print(xxxxxxxxxx)
     rr = escorc(amoc,n_psmsl)
     print(rr)     
     ;print(N_PSMSL)
     ;print(RAPID)
 ;    print(xxxxxxxxxxxxxxx)
     mxlag    = 4
     x_Lead_y = esccr(amoc,n_psmsl,mxlag)
     y_Lead_x = esccr(n_psmsl,amoc,mxlag)    ; switch the order of the series
     print("======================AMOC/SEA LEVEL")
     print(x_Lead_y)
     print(y_Lead_x)

delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(yyyymm)

;delete(n_psmsl)


;spc = spcorr_n(RAPID(0:10),N_PSMSL(1:11),0)    
;print(spc)
;print(xxxxxxxxxxxxxx)
;==============================================================
;==============================================================
;--------------test for significance
lag = 1
;p95 = 0.05
;dAMOC = regCoef(RAPID(0:n1-1),N_PSMSL(0:n1-1))		;m/m
 lag_ts = esacr(n_psmsl,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 lag_ts_map   = esacr(amoc,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)
;  tval = onedtond(dAMOC@tval,dimsizes(dAMOC))
       df = floattointeger(dimsizes(amoc)*(1. -(lag_ts(lag)*lag_ts_map(lag)))\
        /(1. + (lag_ts_map(lag)*lag_ts(lag)))) 
;	printVarSummary(tval)
;   printVarSummary(df)
;   b = tval    ; b must be same size as tval (and df)
;   b = 0.5
;   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
;   prob@long_name = "probability"
;   printVarSummary(prob)
;delete(n1)
;print(prob)
;print(df)
;dAMOCd = dAMOC
;dAMOCd = x_Lead_y(1)
;dAMOCd = where(prob.le.p95,dAMOCd,dAMOCd@_FillValue) ; Significant values
;printVarSummary(dAMOCd)
;printMinMax(dAMOCd,True)
;print(dAMOCd)

print(df)
pr     = rtest(x_Lead_y,df+2,0)
print(pr)
;print(xxxxxxxxxxxxxxx)
;;delete([/tval,df,prob,b/])
;delete(lag_ts)
;delete(lag_ts_map)
;delete(tval)
;delete(df)
;delete(prob)
;delete(b)
;delete(N_PSMSL)
;delete(RAPID)
;delete(mxlag)
delete(n_psmsl)
;delete(x_Lead_y)
;delete(y_Lead_x)
;=========================================================================================
;=========================================================================================
;print(xxxxxxxxxx)
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
   yrStrt = 1900
   yrLast = 2019
   nmos = 12

   xlabel = ispan(yrStrt,yrLast,1)
   printVarSummary(xlabel)

; ==============================================================
; Open & read mask data
; ==============================================================
;===============================================================================================================
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_eastern_USA_1856_2017.nc", "r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview
n_psmsl = rmAnnCycle1D(n_psmsl)
n_psmsl = -n_psmsl 
printMinMax(n_psmsl,True)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;==============================================================
yyyymm = cd_calendar(n_psmsl&time,-1)
yrfrac = yyyymm_to_yyyyfrac(yyyymm, 0.0); not used here
;=============================================================
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
;----------
m_psmsl = n_psmsl	;make a temporary copy for plotting
m_psmsl = dim_standardize_n(m_psmsl,0, 0)  ; normalise
;
N_PSMSL = month_to_season(m_psmsl,"DJF")
print(N_PSMSL)
;
delete(n_psmsl)
;=========================================================================================
;=========================================================================================                         
std = stddev(m_psmsl)
print(std)


list_vars()

;===============================================================================================================
;===============================================================================================================
   yrStrt = 1900
   yrLast = 2009

  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 3   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  
;===============================================================================================================  
;===============================================================================================================
;===============================================================================================================
 base_yrs = (/1900,1910,1920,1930,1940/)
 print(base_yrs)
  do ii = 0,0,1

 ;do ii = 0,dimsizes(base_yrs)-1,1
   yrStrt = base_yrs(ii)
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/rlr_ne/ne_psmsl_ordered_eastern_USA_Canada_1895_2021.nc", "r")

f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview
n_psmsl = rmAnnCycle1D(n_psmsl)
n_psmsl = -n_psmsl 
printMinMax(n_psmsl,True)
delete(TIME)
delete(YYYY)
delete(iYYYY)
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
;************************************************
; calculate spectrum
;************************************************
  sMSL = specx_anal(n_psmsl,d,sm,pct)
  dMSL = specx_ci(sMSL,0.05, 0.95)             ; calc confidence interval
  printVarSummary(sMSL)
  printVarSummary(dMSL)
  print(1./sMSL@frq)
  print("===============================================================")
  print(sMSL@spcx)
  ;print(xxxxxxxxx)
;
 system("/bin/rm -f sMSL_"+yrStrt+".nc")
 system("/bin/rm -f dMSL_"+yrStrt+".nc")
 ncdf = addfile("sMSL_"+yrStrt+".nc","c")
 ncdf ->sMSL = sMSL
 delete(ncdf)
 ncdf = addfile("dMSL_"+yrStrt+".nc","c")
 ncdf ->dMSL = dMSL
; 
 delete(n_psmsl)
 delete(sMSL)
 delete(dMSL)
 end do

; delete(d)
; delete(sm)
; delete(pct)
;===============================================================================================================
;===============================================================================================================
f      = addfile ("sMSL_1900.nc", "r")
sMSL_1900   = f->sMSL
printVarSummary(sMSL_1900)    
delete(f)
f      = addfile ("dMSL_1900.nc", "r")
dMSL_1900   = f->dMSL
printVarSummary(dMSL_1900)    
delete(f)
f      = addfile ("sMSL_1910.nc", "r")
sMSL_1910   = f->sMSL
printVarSummary(sMSL_1910)    
delete(f)
f      = addfile ("dMSL_1910.nc", "r")
dMSL_1910   = f->dMSL
printVarSummary(dMSL_1910)    
delete(f)
;
f      = addfile ("sMSL_1920.nc", "r")
sMSL_1920   = f->sMSL
printVarSummary(sMSL_1920)    
delete(f)
f      = addfile ("dMSL_1920.nc", "r")
dMSL_1920   = f->dMSL
printVarSummary(dMSL_1920)    
delete(f)
;
f      = addfile ("sMSL_1930.nc", "r")
sMSL_1930   = f->sMSL
printVarSummary(sMSL_1930)    
delete(f)
f      = addfile ("dMSL_1930.nc", "r")
dMSL_1930   = f->dMSL
printVarSummary(dMSL_1930)    
delete(f)
;
f      = addfile ("sMSL_1940.nc", "r")
sMSL_1940   = f->sMSL
printVarSummary(sMSL_1940)    
delete(f)
f      = addfile ("dMSL_1940.nc", "r")
dMSL_1940   = f->dMSL
printVarSummary(dMSL_1940)    
delete(f)
;print(xxxxxxxxxxxxxxxx)
;
;
;
 do ii = 0,dimsizes(base_yrs)-1,1
   yrStrt = base_yrs(ii)
;===============================================================================================================
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
n_psmsl   = f->n_psmsl(iYYYY)
printVarSummary(n_psmsl)                              ; variable overview
n_psmsl = rmAnnCycle1D(n_psmsl)
n_psmsl = -n_psmsl
printMinMax(n_psmsl,True)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
N_PSMSL0 = month_to_season(n_psmsl,"DJF")
N_PSMSL0 = dtrend_msg_n(N_PSMSL0&time,N_PSMSL0,True,False,0)
N_PSMSL0 = dim_standardize_n(N_PSMSL0, 0, 0)  ; new variable
;================================================
;************************************************
; calculate spectrum
;************************************************
  sMSL0 = specx_anal(N_PSMSL0,d,sm,pct)
  dMSL0 = specx_ci(sMSL0,0.05, 0.95)             ; calc confidence interval
  printVarSummary(sMSL0)
  printVarSummary(dMSL0)
  
;
 system("/bin/rm -f sMSL0_"+yrStrt+".nc")
 system("/bin/rm -f dMSL0_"+yrStrt+".nc")
 ncdf = addfile("sMSL0_"+yrStrt+".nc","c")
 ncdf ->sMSL0 = sMSL0
 delete(ncdf)
 ncdf = addfile("dMSL0_"+yrStrt+".nc","c")
 ncdf ->dMSL0 = dMSL0
;
 delete(n_psmsl)
 delete(sMSL0)
 delete(dMSL0)
 delete(N_PSMSL0)
end do
  delete(d)
  delete(sm)
  delete(pct)
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
f      = addfile ("sMSL0_1900.nc", "r")
sMSL0_1900   = f->sMSL0
printVarSummary(sMSL0_1900)
delete(f)
f      = addfile ("dMSL0_1900.nc", "r")
dMSL0_1900   = f->dMSL0
printVarSummary(dMSL0_1900)
delete(f)
f      = addfile ("sMSL0_1910.nc", "r")
sMSL0_1910   = f->sMSL0
printVarSummary(sMSL0_1910)
delete(f)
f      = addfile ("dMSL0_1910.nc", "r")
dMSL0_1910   = f->dMSL0
printVarSummary(dMSL0_1910)
delete(f)
;
f      = addfile ("sMSL0_1920.nc", "r")
sMSL0_1920   = f->sMSL0
printVarSummary(sMSL0_1920)
delete(f)
f      = addfile ("dMSL0_1920.nc", "r")
dMSL0_1920   = f->dMSL0
printVarSummary(dMSL0_1920)
delete(f)
;
f      = addfile ("sMSL0_1930.nc", "r")
sMSL0_1930   = f->sMSL0
printVarSummary(sMSL0_1930)
delete(f)
f      = addfile ("dMSL0_1930.nc", "r")
dMSL0_1930   = f->dMSL0
printVarSummary(dMSL0_1930)
delete(f)
;
f      = addfile ("sMSL0_1940.nc", "r")
sMSL0_1940   = f->sMSL0
printVarSummary(sMSL0_1940)
delete(f)
f      = addfile ("dMSL0_1940.nc", "r")
dMSL0_1940   = f->dMSL0
printVarSummary(dMSL0_1940)
delete(f)
;print(xxxxxxxxxxxxxxxx)
;===============================================================================================================
;===============================================================================================================
; ==============================================================
;============================================================

; ==============================================================
;============================================================
  wks = gsn_open_wks("eps","fig0_1")         ; send graphics to PNG file

  plot = new(4,graphic)                ; create graphic array ; only needed if paneling
;====================================================
;====================================================
   r               = True                        ; plot mods desired
   r@gsnDraw       = False                       ; do not draw
   r@gsnFrame      = False                       ; do not advance frame
   r@tiMainString  = ""                       ; title
   r@tiXAxisString = "log [Period (month)]"  ; xaxis
   r@tiYAxisString = "Variance [Sea-level (s.d.~S~2~N~~F21~ month)]"                  ; yaxis
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
    r@tmYLFormat = "@*+^sg"     ; use "e" symbol, don't force 0's (see documentation)
    r@trYMinF 		= 0
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
    r@tiYAxisFontHeightF = 0.020
    r@tiXAxisFontHeightF = 0.020
    r@tmXBLabelFontHeightF = 0.020
    r@tmYLLabelFontHeightF = 0.020
;------------------------------------------------
    plotA0 = gsn_csm_xy(wks,1./sMSL_1900@frq,dMSL_1900,r)
    r@xyLineColors        = (/"deepskyblue","deepskyblue","deepskyblue","deepskyblue"/)
    plotB0 = gsn_csm_xy(wks,1./sMSL_1910@frq,dMSL_1910,r)
    r@xyLineColors        = (/"black","black","black","black"/)
    plotC0 = gsn_csm_xy(wks,1./sMSL_1930@frq,dMSL_1930,r)
    r@xyLineColors        = (/"red","red","red","red"/)
    plotD0 = gsn_csm_xy(wks,1./sMSL_1920@frq,dMSL_1920,r)
    r@xyLineColors        = (/"mediumpurple3","mediumpurple3","mediumpurple3","mediumpurple3"/)
    plotE0 = gsn_csm_xy(wks,1./sMSL_1940@frq,dMSL_1940,r)
    overlay(plotA0,plotB0)
    overlay(plotA0,plotC0)
    overlay(plotA0,plotD0)
    overlay(plotA0,plotE0)
    plot(0) = plotA0
;====================================================
;----------------------------------------------------
    r@trYMinF           = 0
    r@trYMaxF           = 8
   ;r@trXMaxF           = 60
    delete(r@tiYAxisString)
    r@tiXAxisString = "log [Period (year)]"  ; xaxis
    r@tiYAxisString = "Variance [Sea-level (s.d.~S~2~N~~F21~year)]"                  ; yaxis

    plotA2 = gsn_csm_xy(wks,1./sMSL0_1900@frq,dMSL0_1900,r)
    r@xyLineColors        = (/"deepskyblue","deepskyblue","deepskyblue","deepskyblue"/)
    plotB2 = gsn_csm_xy(wks,1./sMSL0_1910@frq,dMSL0_1910,r)
    r@xyLineColors        = (/"black","black","black","black"/)
    plotC2 = gsn_csm_xy(wks,1./sMSL0_1930@frq,dMSL0_1930,r)
    r@xyLineColors        = (/"red","red","red","red"/)
    plotD2 = gsn_csm_xy(wks,1./sMSL0_1920@frq,dMSL0_1920,r)
    r@xyLineColors        = (/"mediumpurple3","mediumpurple3","mediumpurple3","mediumpurple3"/)
    plotE2 = gsn_csm_xy(wks,1./sMSL0_1940@frq,dMSL0_1940,r)
    overlay(plotA2,plotB2)
    overlay(plotA2,plotC2)
    overlay(plotA2,plotD2)
    overlay(plotA2,plotE2)
    plot(2) = plotA2

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
   r2@tiXAxisString = "Time (year)"  ; xaxis
   r2@tiXAxisFontHeightF = 0.017
   r2@tiYAxisFontHeightF = 0.017
   r2@tmXBLabelFontHeightF = 0.012
   r2@tmYLLabelFontHeightF = 0.012
   r2@tiYAxisString = "Monthly indices (s.d.)"                  ; yaxis
   ;r2@tiYAxisString = "Sea-level index (~F33~s~F21~)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.020
   r2@tmYLLabelFontHeightF = 0.020
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   r2@trXMinF             = 1900.                 ; manually set lower limit
   r2@trXMaxF             = 2020.                ;   "          upper
   r2@trYMinF             = -4.0                 ; manually set lower limit
   r2@trYMaxF             =  4.0                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/1./)
   r2@xyLineColors        = (/"blue","blue","transparent"/)
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
   r2@xyDashPatterns      = (/0.,0.,0./)           ; Dash patterns 
   r2@tmYLLabelFontColor = "black"   ; and right axes are same.
   r2@tiYAxisFontColor  = "black"
   ;r2@tmBorderLineColor = "deepskyblue"
   ;r2@tmYLMajorLineColor = "red" ; YR axis uses YL resources, can
   ;r2@tmYLMinorLineColor = "red" ; change by setting tmYUseLeft = False 
   ;r2@tmBorderLineColor =  "red"

   ;===legends
   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   r2@pmLegendSide           = "Top"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.20                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.35                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.08                ; height of legend.
   r2@lgLabelFontHeightF     = 0.014                ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = ""           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"Sea level index","Sea level index","ORAS4"/)         ; explicit labels 
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
   r1@tiXAxisFontHeightF = 0.012
   r1@tiYAxisFontHeightF = 0.012
   r1@tmXBLabelFontHeightF = 0.012
   r1@tmYLLabelFontHeightF = 0.012
       
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
   
 
   plotA1 = gsn_csm_xy (wks,yrfrac,m_psmsl,r2)
  
  
   r2@tiYAxisString = "December-February indices (s.d.)"                  ; yaxis 
   r2@trYMinF             =  -3.0                 ; manually set lower limit
   r2@trYMaxF             =  3.0                 ;   "          upper
   plotA3 = gsn_csm_xy (wks,xlabel,N_PSMSL,r2)

   r2@xyLineColors        = (/"red","blue","transparent"/)
   r2@xyExplicitLegendLabels = (/"RAPID-MOCHA AMOC","Sea level index","ORAS4"/)         ; explicit labels
   r2@xyDashPatterns      = (/2.,0.,0./)           ; Dash patterns 
   r2@pmLegendOrthogonalPosF = -0.42                ; move units down
   r2@pmLegendParallelPosF   = 0.24                 ; move units right
   plotB1 = gsn_csm_xy (wks,yrfracA,amoc,r2)
   overlay(plotA1,plotB1)
   plot(1) = plotA1
   
   plotB3 = gsn_csm_xy(wks,xlabel2,AMOC,r2)
   overlay(plotA3,plotB3)
   plot(3) = plotA3


;==============================================================
;==============================================================   
xa = (/120.,180.,180.,120.,120./)
xb = (/0.,0.,120000.,120000.,0./)

pr0  = True
pr0@gsFillColor =  "grey75"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

d0 =   gsn_add_polygon(wks,plot(0),xa,xb,pr0)

xaa = (/10.,15.,15.,10.,10./)
xbb = (/0.,0.,12.,12.,0./)

dd0 =   gsn_add_polygon(wks,plot(2),xaa,xbb,pr0)


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
  txres@txFontHeightF = 0.037   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
  txid0 = gsn_create_text(wks, "b", txres);
  txid0d = gsn_create_text(wks, "d", txres);
; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.75  ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  annoid0 = gsn_add_annotation(plot(0),txid0,amres)  ; Attach string to plot
  annoid0d = gsn_add_annotation(plot(2),txid0d,amres)  ; Attach string to plot
;===========
 txres@txFontHeightF = 0.032   ; Use same font height as left axis
 amres@amParallelPosF   =  -0.65  ; 0.5 is the right edge of the plot, so
 txid1 = gsn_create_text(wks, "a", txres);
 txid1c = gsn_create_text(wks, "c", txres);
 annoid1 = gsn_add_annotation(plot(1),txid1,amres)  ;
 annoid1c = gsn_add_annotation(plot(3),txid1c,amres)  ;
;----------------------------------------------------------------------
; Code for customizing a legend via simple_legend.
;----------------------------------------------------------------------
;======================================================== 
;========================================================
var_names1 = (/"1900-2009","1910-2009","1920-2009","1930-2009","1940-2009"/)
print(var_names1)
;
var_names2 = (/"1921-2010","1941-2010"/)
print(var_names2)
;
i1_col = (/"blue","deepskyblue","black","red","mediumpurple3"/)
i2_col = (/"red","red"/)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.008
  txres@txJust        = "CenterLeft"
  txres@txFont                = "Helvetica"
  txres@txAngleF      =  0.

  labels  = var_names1
  ypos    = fspan(0.87, 0.92,dimsizes(var_names1))

 do i = 0,dimsizes(labels)-1,1
      txres@txFontColor = i1_col(i)
     gsn_text_ndc(wks,labels(i),0.77,ypos(i),txres)
  end do
delete(i)



;----------------------------------------------------------------------
; Code for customizing a legend via simple_legend.
;----------------------------------------------------------------------
;======================================================== 
;========================================================
var_names2 = (/"~F10~r~F21~  = 0.72","~F10~r~F21~  = 0.50"/)
print(var_names2)
;
;var_names2 = (/"1921-2010","1941-2010"/)
;print(var_names2)
;
i1_col = (/"black","deepskyblue","black","red","mediumpurple3"/)
i2_col = (/"red","red"/)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.010
  txres@txJust        = "CenterLeft"
  txres@txFont        = "Helvetica"
  txres@txAngleF      =  0.

  labels2  = var_names2
  ypos    = fspan(0.59, 1.97,dimsizes(var_names1))

 do i = 0,dimsizes(var_names2)-1,1
      txres@txFontColor = "black"
     gsn_text_ndc(wks,labels2(i),0.58,ypos(i),txres)
     print(i)
end do

;===============================================================
;==============================================================
;=====================================================================
;***********************************************
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelTop                      = 0.99
  res_P@gsnPanelBottom                   = 0.65
  res_P@gsnPanelLeft                     = 0.70
  res_P@gsnPanelRight                    = 0.99
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_P@gsnPanelFigureStrings   = (/"~F0~Peak~F34~;~F0~74 yr~S~-1~N~","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  ;res_P@gsnPanelFigureStringsFontHeightF = 0.015
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(0:0),(/1,1/),res_P)     ; now draw as one plot
;======================================================
;===============================================================
;==============================================================
;=====================================================================
;***********************************************
; panel plot only resources
  res_P0                                  = True            ; modify the panel plot
  res_P0@gsnDraw                          = True       ; don't draw yet
  res_P0@gsnFrame                         = False        ; don't advance frame yet
  res_P0@gsnPanelTop                      = 0.65
  res_P0@gsnPanelBottom                   = 0.31
  res_P0@gsnPanelLeft                     = 0.70
  res_P0@gsnPanelRight                    = 0.99
  res_P0@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_P0@gsnPanelFigureStrings   = (/"~F0~Peak~F34~;~F0~74 yr~S~-1~N~","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  ;res_P0@gsnPanelFigureStringsFontHeightF = 0.015
  ;res_P0@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P0@amJust               =   "CenterRight"
  res_P0@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(2:2),(/1,1/),res_P0)     ; now draw as one plot
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
  resP@gsnPanelRight                    = 0.75
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
 ;======================================================
  resP0                                  = True         ; modify the panel plot
  resP0@gsnDraw                          = True       ; don't draw yet
  resP0@gsnFrame                         = False        ; don't advance frame yet
  resP0@gsnPanelLabelBar                         = False        ; add common colorbar
  resP0@lbLabelAutoStride                        = True         ; auto stride on labels
  resP0@lbOrientation                    = "vertical"
  resP0@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
  resP0@gsnPanelXWhiteSpacePercent       = 1.0    ; default is 1.0
  resP0@gsnPanelTop                      = 0.65
  resP0@gsnPanelBottom                   = 0.31
  resP0@gsnPanelLeft                     = 0.01
  resP0@gsnPanelRight                    = 0.75
  resP0@lbTitleFontHeightF               = .016
  resP0@lbLabelFontHeightF               = 0.014
  resP0@lbLabelFont                      = "Helvetica-Bold"     ; label font
  resP0@lbTitleString                    = "SST anomaly (s.d.)"
  resP0@lbTitleAngleF                    = 90.
  resP0@lbLabelStride                    = 5
  resP0@lbTitleDirection                 = "Across"
 ; resP0@lbTitleString                           = " (s.d.~S~2~N~)SST deviations (~S~o~N~C)"
  resP0@lbTitlePosition                  = "Right"
  resP0@lbLeftMarginF                    = 0.5
  resP0@lbRightMarginF                   = 0.5
  resP0@gsnMaximize                      = True         ; large format
  resP0@gsnPaperOrientation              = "portrait"   ; force portrait
  resP0@gsnLeftString  = ""
 ; resP0@txString            = title+": "+yrStrt+"-"+yrLast
;========================================================================
 gsn_panel(wks,(/plot(3:3)/),(/1,1/),resP0)     ; now draw as one plot
;gsn_panel(wks,plotMAP(0:1),(/1,2/),resP)     ; now draw as one plot
;===============================================================
;========================================================================
;===============================================================
frame(wks) 
end
