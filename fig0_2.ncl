begin
; ==============================================================
; ==============================================================
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
  latS   = -50.
  latN   =  70. 
  lonL   = -90.
  lonR   =  30.
   yrStrt = 1900
   ;yrStrt = 1920
   yrLast = 2009
 ;  yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
 ; print(yyyymm) 
   ;xlabel = ispan(yrStrt,yrLast,1) 
  ;season = "MAM"    ; choose Dec-Jan-Feb seasonal mean
  neof   = 2       ; number of EOFs
  optEOF = True       
  optEOF@jopt = 0   ; This is the default; most commonly used; no need to specify.
;;optEOF@jopt = 1   ; **only** if the correlation EOF is desired
  optETS = False
; ==============================================================
  nmos = 12
  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 3   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  mxlag  =  8
  lag   =   1
  p95   =   0.05  
 nmos = 12
 tau_const = 3600.*12.
; ==============================================================
; Open & read mask data
; ==============================================================
  f0      = addfile ("/home/hnnamchi/work/datasets/mask/tmask_atl.nc", "r")
  atl     = dble2flt(f0->TMASK(0,:,:))
  atl!0 = "lat"
  atl!1 = "lon"
  printVarSummary(atl) 
  atl = atl(::-1,:) ; Flip the latitude dimension 
  atl    = lonFlip(atl)
  printVarSummary(atl)   
 ; print(atl)
 ; atl3 = where(atl.eq.0,atl,atl@_FillValue) 
  atl = mask(atl, atl.gt.1,False)     
  ; print(atl)
   printMinMax(atl,True)
;print(xxxxxxxxxxxxxxxxxxxxxxxxxxxxx)    
; ======================================
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile ("/home/hnnamchi/work/datasets/sst/HadISST_sst.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  hadi   = f->sst(iYYYY,:,:)
  printVarSummary(hadi)                              ; variable overview
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
  printMinMax(hadi,True)
  hadi = where(hadi.gt.-30.,hadi,hadi@_FillValue)
  printMinMax(hadi,True)
  hadi!0 = "time"
  hadi!1 = "lat"
  hadi!2 = "lon"
  printVarSummary(hadi)   
; ==============================================================
  atl@_FillValue  =  hadi@_FillValue
  atl@missing_value =  hadi@_FillValue
; ==============================================================  
  printMinMax(hadi,True)
  hadi  = hadi*conform(hadi,atl,(/1,2/))
  printVarSummary(hadi)
  printMinMax(hadi,True) 
; ==============================================================
; ==============================================================
;  hadi    = lonFlip(hadi)
  printVarSummary(hadi)                              ; note the longitude coord
; ==============================================================
; ==============================================================
; compute climatology and Anomalies
; ==============================================================
  hadiClm   = clmMonTLL(hadi)                             ; (12,lat,lon)
  printVarSummary(hadiClm)
  hadi  = calcMonAnomTLL(hadi,hadiClm)                ; (time, lat,lon)
  printVarSummary(hadi)         
  printMinMax(hadi,True)
  ;===============================================================================================================
;===============================================================================================================
 hadi_d = new(dimsizes(hadi),typeof(hadi),hadi@_FillValue)
 x   = hadi&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = hadi(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    hadi_d(tt::nmos,:,:) = work
 end do
 printVarSummary(hadi_d)
 hadi = hadi_d
;===
delete([/x,x1,hadi_d,work,tt/])  
;=========================================================================================
; =================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
; =================================================================
  rad    = 4.*atan(1.)/180.
  clat   = hadi&lat            
  clat   = sqrt(cos(rad*clat))                 ; gw for gaussian grid
  printVarSummary(clat)
; =================================================================
; weight all observations 
; =================================================================
  xw     = hadi*conform(hadi,clat,1)
  copy_VarMeta(hadi,xw)
  xw@long_name = "Wgt: "+hadi@long_name
 printVarSummary(xw)
; =================================================================
; =================================================================
  wx     = xw({lat|latS:latN},{lon|lonL:lonR},time|:)	 ; convenience, cleaner code
  printVarSummary(wx)
  eof    = eofunc_Wrap(wx, neof,optEOF)      
  eof_ts = eofunc_ts_Wrap(wx,eof,optETS)
  printVarSummary( eof )                         ; examine EOF variables
  printVarSummary( eof_ts )
  eof_ts = -dim_standardize_n(eof_ts,0,1)      ; normalize
 delete([/wx,clat,hadiClm/])
;=====================================================================
 
 list_vars() 
; print(xxxxxxxxxxx)
 ;------------------------------------------
    ;   system("/bin/rm -f eof.nc")   ; remove any pre-existing file
     ;  ncdf = addfile("eof.nc" ,"c")  ; open output netCDF file
      ; ncdf->eof  = eof                          ; 4D
;-------------------------------------------
 ;*******************************************
; North significance test: any pcvar, eval_transpose, eval can be used
;*******************************************
  print("---")
  print("--- eofunc_north ---")
  print("---")
  
  dimp   = dimsizes(hadi)
  ntim   = dimp(0)                                            ; max # eigenvalues possible
  
  prinfo = True
  sig    = eofunc_north(eof@pcvar, ntim, prinfo)            
                                                              
  print("---")
  ;delete(dimp)
;*******************************************
;EOF_TSx = month_to_annual_weighted(yyyymm,eof_ts(time|:,evn|:),1) 
;HADI = month_to_annual_weighted(yyyymm,x_hadi,1) 
;EOF_TS = EOF_TSx(evn|:,year|:)
;delete(EOF_TSx)
;EOF_TS!1 = "time"
;HADI!0 = "time"
;dimp   = dimsizes(HADI)


do ne=0,neof-1
   if(ne.eq.1)then
     eof_ts(ne,:) = -eof_ts(ne,:)    
   end if    
end do
; =================================================================
; Regress
; =================================================================
   eof_regres  =  new((/neof,dimp(1),dimp(2)/),typeof(eof))
   x_var = eof_regres

 ; eof_regres = eof                               ; create an array w meta data
  do ne=0,neof-1
  ; if(ne.eq.1)then
  ;   eof_ts(ne,:) = -eof_ts(ne,:)    
  ; end if    
     eof_regres(ne,:,:) = regCoef(eof_ts(ne,:),hadi(lat|:,lon|:,time|:)) 
     x_var(ne,:,:) = 100.*((escorc(eof_ts(ne,:),hadi(lat|:,lon|:,time|:)))^2) 
   end do

;------------------------------------------
    ;   system("/bin/rm -f eof_ts.nc")   ; remove any pre-existing file
     ;  ncdf = addfile("eof_ts.nc" ,"c")  ; open output netCDF file
      ; ncdf->eof_ts  = eof_ts                          ; 4D               
;-------------------------------------------
 
 eof_regres!1 = "lat"
 eof_regres!2 = "lon"
 
 eof_regres&lat = hadi&lat
 eof_regres&lon = hadi&lon
 
  copy_VarAtts(eof,eof_regres)
 
 
  printVarSummary(eof_regres)   
  
   printVarSummary(eof)   

  
   copy_VarCoords(eof_regres,x_var)  
 
   printVarSummary(x_var)
  
   do yt = 0,neof-1  
   printMinMax(x_var(yt,:,:),True)   
   end do
   delete(yt)

;delete(eof)        ;will be used again

;==============================================================
;==============================================================
;--------------test for significance
 eofd = eof_regres
 
 df = eof_regres

 lag_ts_PC = esacr(eof_ts(evn|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_PC)
 
 lag_ts_HADI   = esacr(hadi(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_HADI)

 tval = onedtond(eof_regres@tval,dimsizes(eof_regres))
 printVarSummary(tval)
 ;print(xxxxxxxxxxx)

    do ne = 0,neof-1,1
           df(ne,:,:) = floattointeger(dimsizes(hadi(:,0,0))*(1. -(lag_ts_PC(ne,lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(ne,lag)))) 
    end do
   printMinMax(df,True)
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printMinMax(prob,True)
   printVarSummary(prob)


eofd = where(abs(prob).lt.p95,eofd,eofd@_FillValue) ; Significant TEMP
printVarSummary(eofd)

delete([/tval,df,prob,b/])
delete(dimp)
;=========================================================================================

f= addfile("/home/hnnamchi/work/datasets/e20c/maskland_RE_ewss_e20c.oper.fc.sfc.3hr.128_180_ewss.regn80sc.1900_2010.nc","r")

 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  u    = f->ewss(iYYYY,:,:)
  u = u/tau_const
  printVarSummary(u)       ; variable overview

delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)

u!0 = "time"
u!1 = "lat"
u!2 = "lon"
printVarSummary(u)  

;===================================
;=====================================
;   u   =  -u
    
   u    = lonFlip(u)
   printVarSummary(u)   	; note the longitude coord
;===================================
; ==============================================================  
  printMinMax(u,True)
 
  u  = u*conform(u,atl,(/1,2/))
    
  printVarSummary(u)
  
  printMinMax(u,True) 
     
; ==============================================================
; ==============================================================
;=====================================
 uClm      = clmMonTLL(u)
   printVarSummary(uClm)                   ; (12,nlat,nlon)   
   u         = calcMonAnomTLL (u, uClm)    ; replace with anonamlies
   u@long_name = "ANOMALIES: "+u@long_name
;=====================================

;=====================================
;===================================
 u_d = new(dimsizes(u),typeof(u),u@_FillValue)
 x   = u&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = u(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    u_d(tt::nmos,:,:) = work
 end do

 printVarSummary(u_d)
 u = u_d

;===
delete([/x,x1,u_d,work/])  
;=================================
;=================================
;------------------------------------------
   dimp   = dimsizes(u)
   E_U  =  new((/neof,dimp(1),dimp(2)/),typeof(u))


        do ne=0,neof-1
                E_U(ne,:,:) = regCoef(eof_ts(evn|ne,time|:),u(lat|:,lon|:,time|:))
        end do
 E_U!1 = "lat"
 E_U!2 = "lon"

 E_U&lat = u&lat
 E_U&lon = u&lon
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 E_Ud = E_U
 
 df = E_U

 lag_ts = esacr(eof_ts(evn|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(u(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(E_U@tval,dimsizes(E_U))
    

    do ne = 0,neof-1,1
           df(ne,:,:) = floattointeger(dimsizes(u(:,0,0))*(1. -(lag_ts(ne,lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(ne,lag)))) 
    end do
   printMinMax(df,True)
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printMinMax(prob,True)
   printVarSummary(prob)

delete(E_Ud)
E_Ud = prob

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(u)
delete(dimp)
 
list_vars()

; ==============================================================
f= addfile("/home/hnnamchi/work/datasets/e20c/maskland_RE_nsss_e20c.oper.fc.sfc.3hr.128_181_nsss.regn80sc.1900_2010.nc","r")

 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  v    = f->nsss(iYYYY,:,:)
  v = v/tau_const
  printVarSummary(v)       ; variable overview

delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)

printVarSummary(v)  

v!0 = "time"
v!1 = "lat"
v!2 = "lon"
printVarSummary(v)  

;===================================
;=====================================
;For v
;   v   =  -v
   v    = lonFlip(v)
   printVarSummary(v)   	; note the longitude coord
;=================================== 
 ; ==============================================================  
  printMinMax(v,True)
 
  v  = v*conform(v,atl,(/1,2/))
    
  printVarSummary(v)
  
  printMinMax(v,True) 
     
  
;=====================================

   vClm      = clmMonTLL(v)
   printVarSummary(vClm)                   ; (12,nlat,nlon)   
   v         = calcMonAnomTLL (v, vClm)    ; replace with anonamlies
   v@long_name = "ANOMALIES: "+v@long_name
;=====================================

;=====================================
;===================================
 v_d = new(dimsizes(v),typeof(v),v@_FillValue)
 x   = v&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = v(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    v_d(tt::nmos,:,:) = work
 end do

 printVarSummary(v_d)
 v = v_d

;===
delete([/x,x1,v_d,work/])  
;=============================
;------------------------------------------
 dimp   = dimsizes(v)
 E_V  =  new((/neof,dimp(1),dimp(2)/),typeof(v))

         do ne=0,neof-1
                E_V(ne,:,:) = regCoef(eof_ts(evn|ne,time|:),v(lat|:,lon|:,time|:))
        end do
 E_V!1 = "lat"
 E_V!2 = "lon"

 E_V&lat = v&lat
 E_V&lon = v&lon
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 E_Vd = E_V
 
 df = E_V

 lag_ts = esacr(eof_ts(evn|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(v(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(E_V@tval,dimsizes(E_V))


    do ne = 0,neof-1,1
           df(ne,:,:) = floattointeger(dimsizes(v(:,0,0))*(1. -(lag_ts(ne,lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(ne,lag)))) 
    end do
   printMinMax(df,True)
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printMinMax(prob,True)
   printVarSummary(prob)


delete(E_Vd)
E_Vd = prob
;E_Vd = where(abs(prob).ge.p95,E_Vd,E_Vd@_FillValue) ; Significant TEMP
;printVarSummary(E_Vd)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(v)
delete(dimp)
;print(xxxxxxxxxxxxxxxx)
;==============================================================
;==============================================================
;d_Ud = E_Ud
;d_Vd = E_Vd
;
E_U  = where((E_Ud.lt.p95.or.E_Vd.lt.p95),E_U,E_U@_FillValue) 
E_V  = where((E_Vd.lt.p95.or.E_Ud.lt.p95),E_V,E_V@_FillValue)
;
printMinMax(E_U,True)
printMinMax(E_V,True) 
;===============================================================================================================
;===============================================================================================================

delete(atl)


;==============================================================
;==================================================================
  sPC0 = specx_anal(eof_ts(0,:),d,sm,pct)
  dPC0 = specx_ci(sPC0,0.05,0.95)             ; calc confidence interval
    
;----
  
  sPC1 = specx_anal(eof_ts(1,:),d,sm,pct)
  dPC1 = specx_ci(sPC1,0.05,0.95)             ; calc confidence interval
  
  ;sPC2 = specx_anal(eof_ts(2,:),d,sm,pct)
  ;dPC2 = specx_ci(sPC2,0.05,0.95)             ; calc confidence interval

 ; sigAMV = where(sdof0x@spcx.ge.splt0x(3,:),1/sdof0x@frq,sdof0x@_FillValue)
 ; print(sigAMV)
 
; =================================================================
; Extract the YYYYMM from the time coordinate 
; associated with eof_ts [same as x&time] 
; =================================================================
  yyyymm = cd_calendar(eof_ts&time,-1)  
  yrfrac = yyyymm_to_yyyyfrac(yyyymm, 0.0); not used here

;============================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("eps","fig0_2")         ; send graphics to PNG file

  plot = new(3,graphic)                ; create graphic array ; only needed if paneling
     
  
; EOF patterns

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -0.3            ; set min contour level
  res@cnMaxLevelValF       =  0.3             ; set max contour level
  res@cnLevelSpacingF      =  0.03            ; set contour spacing     
 ; res@gsnPolar             = "SH"
  res@mpFillOn             = True        ; turn off map fill
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN -8
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR-15
   res@pmTickMarkDisplayMode = "Never"
  ;res@tmXBLabelsOn          = False      ;no lon ticks
  res@tmYLMode   = "Explicit"
  res@tmYLValues = ispan(-40,60,20)
  res@tmYLLabels = (/"40~S~o~N~S","20~S~o~N~S", "0~S~o~N","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N"/)
  res@tmYLMinorValues  = ispan(-50,70,10) 
  ;res@tmYLMinorValues = res@tmYLValues
  res@tmYLLabelFontHeightF = 0.012
  res@tmXBMode   = "Explicit"
  res@tmXBValues = ispan(-90,30,30)
  res@tmXBLabels = (/"90~S~o~N~W", "60~S~o~N~W","30~S~o~N~W","0~S~o~N", \
                     "30~S~o~N~E"/)
  res@tmXBMinorValues  = ispan(-90,30,15)
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012
  res@mpGeophysicalLineColor  = "transparent" 
 ; res@mpFillDrawOrder      = "preDraw"
  res@gsnAddCyclic         = False
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False      ; turn off individual lb's
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF   = .020 
  res@lbLabelFontHeightF   = 0.020
  res@lbLabelFont   	   = "Helvetica-Bold"     ; label font
  res@lbTitleString        = ""
  res@lbTitleAngleF        = 90.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  res@lbLabelStride        = 5
  res@tiYAxisFontThicknessF  = 1.
  res@tmXBMajorThicknessF    = 1.
  res@tmXBMinorThicknessF    = 1.
  res@tmYLMajorThicknessF    = 1.
  res@tmYLMinorThicknessF    = 1.
  res@tmYLLabelFontThicknessF= 1.
  res@tmBorderThicknessF     = 1.
  res@vcRefAnnoOrthogonalPosF   = -0.35             ; move ref vector up
  res@vcRefAnnoParallelPosF	= 0.35
  res@vcRefMagnitudeF           = 0.01             ; define vector ref mag
  res@vcRefAnnoString1          = "0.01 N m~S~-2~N~ "
  res@vcRefLengthF            = 0.045            ; define length of vec ref
 ;; res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.035
  res@vcLineArrowColor        = "Foreground"           ; change vector color
  res@vcRefAnnoBackgroundColor = "transparent"
  res@vcLineArrowThicknessF   = 2.0
  res@vcLineArrowHeadMinSizeF = 0.008
  ;res@vcRefAnnoPerimColor     = 1
  res@vcRefAnnoPerimOn		= False
  res@gsnScalarContour	      = True
;===============
;-----------------------------------------------------------------------------
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
 ; res1@cnMinLevelValF       =  -1.0            ; set min contour level
 ; res1@cnMaxLevelValF       =  1.0             ; set max contour level
  res1@cnLevelSpacingF      =  .05            ; set contour spacing  
  res1@gsnAddCyclic         = False   
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString   = ""
  res1@gsnLeftString   = ""
  ;res1@cnLineLabelDensityF =  5
;=============================================================================
;=============================================================================
;=============================================================================
;-----------------------------------------------------------------------------
  res2                      =   True
  res2@gsnDraw              = False        ; don't draw yet
  res2@gsnFrame             = False        ; don't advance frame yet
  res2@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res2@cnMinLevelValF       =  10.0            ; set min contour level
  res2@cnMaxLevelValF       =  100.0             ; set max contour level
  res2@cnLevelSpacingF      =  10.0            ; set contour spacing  
  res2@gsnAddCyclic         = False   
  res2@cnLinesOn            = True       ; True is default
  res2@cnLineLabelsOn       = True        ; True is default
  res2@cnInfoLabelOn        = False
  res2@cnLineColor	    = "white"
  res2@cnLineThicknessF	    =  1.0
  res2@gsnRightString   = ""
  res2@gsnLeftString   = ""
  ;res2@cnLineLabelDensityF =  5
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #4
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;=============================================================================
  ; res@gsnCenterStringFontHeightF   = 0.013
  ; res@tiMainFontHeightF	    = 0.015
 
 
  texts = (/"Atlantic Multidecadal Variability, AMV","Atlantic Decadal Variability, ADV"/)
  print(texts)

    n = 1
        
   res@gsnLeftString  = "Second empirical orthogonal function"
  ;res@gsnLeftString  = "EOF "+(n+1)
   res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
   res@gsnCenterString   = ""                       ; title

   plotA0 = gsn_csm_vector_scalar_map_ce(wks,E_U(n,:,:),E_V(n,:,:),eof_regres(n,:,:),res)     
;  plotA0 = gsn_csm_contour_map_ce(wks,eof_regres(n,:,:),res)
   plotB0 = gsn_csm_contour(wks,x_var(n,:,:),res2)
   plotC0 = gsn_csm_contour(wks,eofd(n,:,:),res1)
   plotC0 = gsn_contour_shade(plotC0,10.,p95,opt)
   
   
  overlay(plotA0,plotB0)
  overlay(plotA0,plotC0)
  plot(0) = plotA0

;*******************************************
   r		   = True
   r@gsnDraw       = False                       ; do not draw
   r@gsnFrame      = False                       ; do not advance frame
   r@tiXAxisString = "log[Period(month)]"  ; xaxis
   ;r@tiYAxisString = "Variance (~F33~s~S~2~N~~F21~ month)"         
   r@tiYAxisString = "Variance (s.d.~S~2~N~~F21~ month)"   
   r@xyLineThicknesses   = (/3.,1.,1.,1./)       ; Define line thicknesses 
   r@xyDashPatterns      = (/0,0,16,16/)           ; Dash patterns 
   r@xyLineColors        = (/"seagreen3","seagreen3","seagreen3","seagreen3"/)
   r@trXLog        = True                 ; log scaling           ; yaxis
   r@trYMaxF       =    240.
   ;r@trYMinF       =    0.
   r@trXMaxF       =        1000.
   r@vpHeightF		   = 0.42                   ; change aspect ratio of plot
   r@vpWidthF		    = 0.40          
   r@trXReverse 	= False 
   r@tmXBLabelFontHeightF = 0.012
   r@tmYLLabelFontHeightF = 0.012
   r@tiXAxisFontHeightF   = 0.012
   r@tiYAxisFontHeightF   = 0.012
   r@tiMainFontHeightF   = 0.013
   r@tiMainString	  =  "Second principal component"
    r@tiYAxisFontThicknessF  = 1.
    r@tmXBMajorThicknessF    = 1.
    r@tmXBMinorThicknessF    = 1.
    r@tmYLMajorThicknessF    = 1.
    r@tmYLMinorThicknessF    = 1.
    r@tmYLLabelFontThicknessF= 1.
    r@tmBorderThicknessF     = 1.
   r@tmXBMode   = "Explicit"
   r@tmXBValues = (/2.,12.,30.,60.,120.,300.,900./)
   r@tmXBLabels  =  r@tmXBValues
   r@pmLegendDisplayMode    = "Never"            ; turn on legend
   r@pmLegendSide           = "Top"               ; Change location of 
   r@pmLegendParallelPosF   = 0.60                  ; move units right
   r@pmLegendOrthogonalPosF = -0.40                ; move units down
   r@pmLegendWidthF         = 0.08                ; Change width and
   r@pmLegendHeightF        = 0.11                ; height of legend.
   r@lgPerimOn              = False               ; turn off box around
   r@lgLabelFontHeightF     = .015                 ; label font height
   r@xyExplicitLegendLabels = (/"95% confidence level","SST indices, unfiltered","PCs, 5-year low-pass"/)         ; create explicit labels

;************************************************
;
;************************************************
   r@gsnCenterStringFontHeightF  =  0.020
 ;r@tiMainString  = "Second principal component"                       ; title
  r@gsnCenterString  = ""     
  ;plotR(0) = gsn_csm_xy(wks,1./sPC0@frq,dPC0,r)  ; create plot 
   
   ;r@tiMainString  = ""                       ; title
   r@trYMaxF       =    100.                ; title
   r@pmLegendDisplayMode    = "Never"            ; turn on legend
   plot(1) = gsn_csm_xy(wks,1./sPC1@frq,dPC1,r)  ; create plot 
   
   
   r@tiMainString  = ""   
;************************************************
;; EOF time series  [bar form]
  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling               
  rts@vpHeightF = 0.36        ; Changes the aspect ratio
  rts@vpWidthF  = 0.40
  ;rts@vpXF      = 0.10        ; change start locations
  ;rts@vpYF      = 0.75        ; the plot
  rts@trYMaxF       =    4.
  rts@trYMinF       =    -4.
  rts@trXMaxF       =    yrLast+10
  rts@trXMinF       =    yrStrt-10
  rts@gsnYRefLine           = 0.              ; reference line
  rts@gsnYRefLineThicknesses       = (/0./)
  rts@tiXAxisString = "Time (year)"  ; xaxis   
  rts@tmXBLabelAngleF	   =  45.
  rts@tmXBLabelFontHeightF = 0.025
  rts@tmYLLabelFontHeightF = 0.025
  rts@tiXAxisFontHeightF   = 0.025
  rts@tiYAxisFontHeightF   = 0.025
  rts@tmYLFormat 	   = "@*+^sg"
  rts@gsnAboveYRefLineColor = "red"              ; above ref line fill red
  rts@gsnBelowYRefLineColor = "deepskyblue"             ; below ref line fill blue
  rts@xyLineColor	    = "transparent"
  rts@tmXTBorderOn          =  False
  rts@tmXTOn                =  False
  rts@tmYRBorderOn          =  False
  rts@tmYROn                =  False
  ;rts@gsnAboveYRefLineColor = "red"              ; above ref line fill red
  ;rts@gsnBelowYRefLineColor = "deepskyblue"             ; below ref line fill blue
  ;rts@xyLineColor	    = "transparent"
  ;rts@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;rts@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ;rts@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;rts@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ;rts@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;rts@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ;rts@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;rts@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  rts@tiYAxisString = "Normalized PCs" 	           ; y-axis label   
  rts@tiYAxisFontThicknessF  = 1.
  rts@tmXBMajorThicknessF    = 1.
  rts@tmXBMinorThicknessF    = 1.
  rts@tmYLMajorThicknessF    = 1.
  rts@tmYLMinorThicknessF    = 1.
  rts@tmYLLabelFontThicknessF= 1.
  rts@tmBorderThicknessF     = 1.  

;------------------------------
; create individual plots
; rts@tiYAxisString = "Principal component (~F33~s~N~)"        ; y-axis label      

  rts@tiYAxisString = "s. d."        ; y-axis label 
  plot(2) = gsn_csm_xy (wks,yrfrac,eof_ts(n,:),rts)
 
;=====================================================

;xa = (/(1./(12.*12.)),(1./(13.*12.)),(1./(13.*12.)),(1./(12.*12.)),(1./(12.*12.))/)

xa = (/150.,150.,150.,150.,150./)
xb = (/0.,0.,1000.,1000.,0./)

pr0  = True
pr0@gsFillColor =  "black"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

;d0 =   gsn_add_polygon(wks,plotR(0),xa,xb,pr0)
;d1 =   gsn_add_polyline(wks,plotR(1),xa,xb,pr0)


;=====================================================
;=====================================================
;=====================================================
;=====================================================
  txres		= True
  txres@txFontColor      =  "black"
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.025   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
  txid0 = gsn_create_text(wks, "a", txres);
; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.68  ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  annoid0 = gsn_add_annotation(plot(0),txid0,amres)  ; Attach string to plot

 txres@txFontHeightF = 0.020   ; Use same font height as left axis
 txid1 = gsn_create_text(wks, "c", txres);
 annoid1 = gsn_add_annotation(plot(1),txid1,amres)  

  txres@txFontHeightF = 0.055   ; Use same font height as left axis
  amres@amParallelPosF   =  -0.80
  txid2 = gsn_create_text(wks, "b", txres);
  annoid2 = gsn_add_annotation(plot(2),txid2,amres)  
;===============================================================
;===============================================================

xa = (/400.,1800.,1800.,400.,400./)	;>60 years
xb = (/0.,0.,120000.,120000.,0./)

pr0  = True
pr0@gsFillColor =  "grey85"
pr0@tfPolyDrawOrder = "PreDraw"
pr0@gsFillIndex  =  0

;d0 =   gsn_add_polygon(wks,plotR(0),xa,xb,pr0)
xa = (/120.,180.,180.,120.,120./)	;10-15 year
d1 =   gsn_add_polygon(wks,plot(1),xa,xb,pr0)
;========================================================================
;========================================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True       ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	       	= True        ; add common colorbar
  resP@lbLabelAutoStride             	= True         ; auto stride on labels
  resP@lbOrientation                    = "horizontal"
  resP@gsnPanelYWhiteSpacePercent      	= 0.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 1.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.95
  resP@gsnPanelBottom                   = 0.35
  resP@gsnPanelLeft                     = 0.01
  resP@gsnPanelRight                    = 0.50
  resP@lbTitleFontHeightF	        = 0.013 
  resP@lbLabelFontHeightF               = 0.013
  resP@lbLabelFont   		        = "Helvetica-Bold"     ; label font
  resP@lbTitleString                  	= "Sea surface temperature (K)"
  ;resP@lbTitleAngleF                    = 90.
  resP@lbTitleDirection                 = "Across"
 ; resP@lbTitleString                  		= " (s.d.~S~2~N~)SST deviations (~S~o~N~C)"
  ;resP@lbTitlePosition                  = "Right"
  ;resP@lbTopMarginF                     = 0.40
  ;resP@lbBottomMarginF                  = 0.35
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  ;resP@lbLeftMarginF                    = 0.8
  ;resP@lbRightMarginF                   = 0.2
  resP@lbTopMarginF                     = 0.30
  resP@lbLabelStride       		= 10

   res@gsnLeftString  = ""
 ; resP@txString            = title+": "+yrStrt+"-"+yrLast

;*******************************************
;  gsn_panel(wks,plot,(/1,neof/),resP)     ; now draw as one plot
;=====================================================
gsn_panel(wks,plot(0:0),(/1,1/),resP)     ; now draw as one plot
;=====================================================
;=====================================================
;=====================================================
;===============================================================
;===============================================================
;===============================================================
;***********************************************
; panel plot only resources
  res_P                     = True            ; modify the panel plot
  res_P@gsnDraw             = True       ; don't draw yet
  res_P@gsnFrame            = False        ; don't advance frame yet
  res_P@gsnPanelTop         = 0.95
  res_P@gsnPanelBottom      = 0.43
  res_P@gsnPanelLeft        = 0.50
  res_P@gsnPanelRight       = 0.99
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_P@gsnPanelFigureStrings   = (/"~F0~Peak~F34~;~F0~74 yr~S~-1~N~","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  ;res_P@gsnPanelFigureStringsFontHeightF = 0.015
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust              =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(1:1),(/1,1/),res_P)     ; now draw as one plot
;***********************************************
; panel plot only resources
  res_L                     = True            ; modify the panel plot
  res_L@gsnDraw             = True       ; don't draw yet
  res_L@gsnFrame            = False        ; don't advance frame yet
  res_L@gsnPanelTop         = 0.90
  res_L@gsnPanelBottom      = 0.70
  res_L@gsnPanelLeft        = 0.60
  res_L@gsnPanelRight       = 0.80
  res_L@gsnPaperOrientation = "portrait"   ; force portrait
  ;res_L@gsnPanelFigureStrings   = (/"~F0~Peak~F34~;~F0~74 yr~S~-1~N~","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  ;res_L@gsnPanelFigureStringsFontHeightF = 0.015
  ;res_L@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_L@amJust               =   "CenterRight"
  res_L@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(2:2),(/1,1/),res_L)     ; now draw as one plot
;======================================================
;======================================================
;======================================================
frame(wks) 
end
