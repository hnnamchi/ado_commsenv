;*************************************************
begin
; ==============================================================
; ==============================================================
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
  latS   = -50.
  latN   =  70. 
  lonL   = -90.
  lonR   =  25.
 
  nmos = 12
  nlags = 7
  mxlag = 7
  lag   =   1
  tau_const = 3600.*12.
  p95   =   0.05
  nmos = 12
  yrStrt = 1900
  yrLast = 2009	;
  ;yrLast = 2007         ;

  season_w = "DJF"
  
  yyyymm = yyyymm_time(yrStrt, yrLast, "integer")
;=========================================================================================
;=========================================================================================
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
  atl = where(atl.eq.1,atl,atl@_FillValue) 
 ; atl = where(atl.eq.2,atl@_FillValue,atl) 
  ;atl2 = atl
  ;atl = mask(atl, atl.gt.1,True)     
  ; print(atl)
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
  ;hadi = hadi/100.
  hadi = where(hadi.gt.-10.and.hadi.lt.40.,hadi,hadi@_FillValue)
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
 ;hadi = center_finite_diff_n(hadi,1,False,0,0)
  hadi!0 = "time"
  hadi!1 = "lat"
  hadi!2 = "lon"
  printVarSummary(hadi)   
; ==============================================================
 ; ==============================================================
  atl@_FillValue  =  hadi@_FillValue
  atl@missing_value =  hadi@_FillValue
; hadi    = lonFlip(hadi)
; ==============================================================  
  printMinMax(hadi,True)
  hadi  = hadi*conform(hadi,atl,(/1,2/))
  printVarSummary(hadi)
  printMinMax(hadi,True) 
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
;=========================================================================================
;------------------------------------------------------------------
hadi_1 = hadi
delete(hadi)
hadi = month_to_annual_weighted(yyyymm,hadi_1,1) 
;hadi = month_to_season(hadi_1,season_w)
hadi!0 = "time"
;hadi = dtrend_n(hadi,False,0)
hadi = dtrend_msg_n(hadi&time,hadi,True,False,0)
printVarSummary(hadi)
;hadi = dim_standardize_n(hadi,0,0)      ; normalize
delete(hadi_1)

list_vars()
;----------------------------------------------------------------------------------------------
;=========================================================================================
f= addfile("/home/hnnamchi/work/datasets/e20c/maskland_RE_ewss_e20c.oper.fc.sfc.3hr.128_180_ewss.regn80sc.1900_2010.nc","r")

 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  u    = f->ewss(iYYYY,:,:)
  u = u/tau_const	;converts to N/m2
  ;u   = u*1000.		;scales to e**-3 
  printVarSummary(u)       ; variable overview


;copy_VarCoords(hadi(0,:,:),u(0,:,:)) 
;printVarSummary(u)       ; variable overview

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
;For u
 ;  u   =  -u
    
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

;For OI ERSST
   uClm      = clmMonTLL(u)
   printVarSummary(uClm)                   ; (12,nlat,nlon)   
   u         = calcMonAnomTLL (u, uClm)    ; replace with anonamlies
   u@long_name = "ANOMALIES: "+u@long_name
;=====================================
;=========================================================================================
;=========================================================================================
u_1 = u
delete(u)
u = month_to_annual_weighted(yyyymm,u_1,1) 
;u = month_to_season(u_1,season_w)
;u = dtrend_n(u,False,0)
u!0 = "time"
u = dtrend_msg_n(u&time,u,True,False,0)
printVarSummary(u)
delete(u_1)
;==============================================================
f= addfile("/home/hnnamchi/work/datasets/e20c/maskland_RE_nsss_e20c.oper.fc.sfc.3hr.128_181_nsss.regn80sc.1900_2010.nc","r")

 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  v    = f->nsss(iYYYY,:,:)
  v = v/tau_const	;converts to N/m2
  ;v = v*1000.		;scales to e**-3
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
;=====================================
   vClm      = clmMonTLL(v)
   printVarSummary(vClm)                   ; (12,nlat,nlon)   
   v         = calcMonAnomTLL (v, vClm)    ; replace with anonamlies
   v@long_name = "ANOMALIES: "+v@long_name
;=====================================
;===================================
;---------------------------------------------------------------------------
v_1 = v
delete(v)
v = month_to_annual_weighted(yyyymm,v_1,1) 
;v = month_to_season(v_1,season_w)
;v = dtrend_n(v,False,0)
v!0 = "time"
v = dtrend_msg_n(v&time,v,True,False,0)
printVarSummary(v)
delete(v_1)
print("==============================================================")

;----------------------------------------------------------------------------------------------
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
; f      = addfile ("/home/hnnamchi/work/datasets/slp/HadSLP2_slp.mnmean.real.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/20cr/v2c/prmsl.mon.mean.nc", "r")
f      = addfile("/home/hnnamchi/work/datasets/e20c/msl_e20c.oper.an.sfc.3hr.128_151_msl.regn80sc.1900_2010.nc", "r")
; TIME   = f->time
 TIME   = f->initial_time0_hours
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  slp   = f->MSL_GDS4_SFC(iYYYY,:,:)
 ; slp   = f->prmsl(iYYYY,:,:)
  ;slp   = short2flt(f->slp(iYYYY,:,:))
  printVarSummary(slp)                              ; variable overview
  slp = slp/100.       ;converts from Pa to hPa


delete([/f,TIME,YYYY,iYYYY/])
slp!0 = "time"
slp!1 = "lat"
slp!2 = "lon"
printVarSummary(slp)
;===================================
 ;  slp   =  -slp
   slp    = lonFlip(slp)
   printVarSummary(slp)        ; note the longitude coord
;===================================
   ;===================================
m_slp = dim_avg_n_Wrap(slp,0)
printVarSummary(m_slp)
printMinMax(m_slp,True)
;===================================

; ==============================================================  
  printMinMax(slp,True)
 ; slp  = slp*conform(slp,atl,(/1,2/))
  printVarSummary(slp)
  printMinMax(slp,True)
;===================================
;=========================================================================================
slp_1 = slp
delete(slp)
;slp = month_to_annual_weighted(yyyymm,slp_1,1)
slp = month_to_season(slp_1,season_w)
slp!0 = "time"
;===================================
;===================================
slp = dtrend_msg_n(slp&time,slp,True,False,0)
printVarSummary(slp)
delete(slp_1)
print("==============================================================")
; ==============================================================
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
; ==============================================================
; ==============================================================
;f      = addfile ("/home/hnnamchi/work/datasets/gulev/qnet_qe5drer_qh5drer_8007_ssc_btr.dat.nc", "r")
;qnet   = f->qnet({year|yrStrt:yrLast},season|0,lat|:,lon|:)
;printVarSummary(qnet)                              ; variable overview
;delete(f)
;============================================================== 
;==============================================================
;==============================================================
f= addfile("/home/hnnamchi/work/datasets/e20c/RE_lh_sh_e20c_1900_2010.nc","r")
;f= addfile("/home/hnnamchi/work/datasets/e20c/RE_qnet_e20c_1900_2010.nc","r")
 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  qnet    = f->lh_sh(iYYYY,:,:)
;  qnet    = f->qnet(iYYYY,:,:)
  printVarSummary(qnet)       ; variable overview
delete([/f,TIME,YYYY,iYYYY/])
qnet!0 = "time"
qnet!1 = "lat"
qnet!2 = "lon"
printVarSummary(qnet)
;===================================
;; qnet   =  -qnet
   qnet    = lonFlip(qnet)
   printVarSummary(qnet)        ; note the longitude coord
;===================================
; ==============================================================  
  printMinMax(qnet,True)
  qnet  = qnet*conform(qnet,atl,(/1,2/))
  printVarSummary(qnet)
  printMinMax(qnet,True)
;=====================================
   qnetClm      = clmMonTLL(qnet)
   printVarSummary(qnetClm)                   ; (12,nlat,nlon)   
   qnet         = calcMonAnomTLL(qnet, qnetClm)    ; replace with anonamlies
   qnet@long_name = "ANOMALIES: "+qnet@long_name
;===================================
;=========================================================================================
qnet_1 = qnet
delete(qnet)
qnet = month_to_annual_weighted(yyyymm,qnet_1,1)
;qnet = month_to_season(qnet_1,season_w)
;;qnet = dtrend_n(qnet,False,0)
qnet!0 = "time"
qnet = dtrend_msg_n(qnet&time,qnet,True,False,0)
printVarSummary(qnet)
delete(qnet_1)
print("==============================================================")
;======================================================================================
;======================================================================================
;======================================================================================
;======================================================================================
; Open the file: Read only the user specified period
;======================================================================================
f      = addfile ("/home/hnnamchi/work/datasets/ora20c/zos_opa0_9_ora20c_1m_1900_2009_grid_1x1.nc", "r")
 TIME   = f->time
  ;TIME   = f->initial_time0_hours
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  ssh   = f->zos(iYYYY,:,:)
  printVarSummary(ssh)                              ; variable overview
  ssh!0 = "time"
  ssh!1 = "lat"

  ssh  =  ssh*100.	;converts to cm

delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;==============================================================
  printVarSummary(ssh)
  ssh = ssh(:,::-1,:) 		; Flip the latitude dimension 
;
  ssh    = lonFlip(ssh)
  printVarSummary(ssh)                              ; note the longitude coord
;==============================================================  
  printMinMax(ssh,True)
  ssh  = ssh*conform(ssh,atl,(/1,2/))
  printVarSummary(ssh)
  printMinMax(ssh,True)
;=================================================================
;=================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
  rad    = 4.*atan(1.)/180.
  clat   = ssh&lat
  clat   = sqrt(cos(rad*clat))                 ; gw for gaussian grid
; weight all observations 
  wssh   = ssh                                 ; copy meta data
  wssh   = ssh*conform(ssh,clat,1)
  wssh@long_name = "Wgt: "+wssh@long_name
  printVarSummary(wssh)
;==============================================================
 g_wssh = wgt_areaave_Wrap(wssh,clat,1.0,0)
 print(g_wssh)
 printVarSummary(g_wssh)
 ;
 delete(rad)
 delete(wssh)
;
 printMinMax(ssh,True)
; ssh = ssh - conform(ssh,g_wssh,0)
 delete(g_wssh)
 printMinMax(ssh,True)
; ==============================================================
; ==============================================================
; compute climatology and Anomalies
; ==============================================================
  sshClm   = clmMonTLL(ssh)                             ; (12,lat,lon)
  printVarSummary(sshClm)
  ssh  = calcMonAnomTLL(ssh,sshClm)                ; (time, lat,lon)
  printVarSummary(ssh)
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;=========================================================================================
ssh_1 = ssh
delete(ssh)
ssh = month_to_annual_weighted(yyyymm,ssh_1,1)
;ssh = month_to_season(ssh_1,season_w)
;ssh = dtrend_n(ssh,False,0)
ssh!0 = "time"
ssh = dtrend_msg_n(ssh&time,ssh,True,False,0)
printVarSummary(ssh)
delete(ssh_1)
print("==============================================================")

;======================================================================================
;======================================================================================
;======================================================================================
;======================================================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hnnamchi/work/geomar/adv/sub0/eof_ts.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  pc   = f->eof_ts(1,iYYYY)
  printVarSummary(pc)                              ; variable overview
  pc = dim_standardize_n(pc,0,0)      ; normalize
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;
PC = month_to_annual_weighted(yyyymm,pc,1) 
;PC = month_to_season(pc,season_w)
;delete(pc)
PC!0 = "time"
hadi!0 = "time"
;=====================================================================
;
df = hadi(0,:,:)
lag_ts_PC = esacr(PC(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(hadi(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(hadi(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag)))) 
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================
;
dimp   = dimsizes(hadi)
;
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PC(0:dimsizes(PC)-ii-1),hadi(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_var(ii,:,:) = 100.*((escorc_n(PC(0:dimsizes(PC)-ii-1),hadi(ii:,:,:),0,0))^2)     ;
print(ii)
;printMinMax(x_Lead_y@tval,True)
;tval = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
;print("==================")
;printMinMax(tval,True)
;printVarSummary(tval)
end do
printVarSummary(x_Lead_y_tval)
delete(ii)
copy_VarCoords(hadi(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(hadi(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)
;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
printVarSummary(y_Lead_x)
y_Lead_x_tval = y_Lead_x

do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PC(ii:),hadi(0:dimsizes(PC)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_var(ii,:,:) = 100.*((escorc_n(PC(ii:),hadi(0:dimsizes(PC)-ii-1,:,:),0,0))^2)     ;
print(ii)
;printVarSummary(y_Lead_x(ii,:,:))
;printMinMax(y_Lead_x@tval,True)
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))	;takes out the t-values
end do
printVarSummary(y_Lead_x)
delete(ii)
copy_VarCoords(hadi(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;
copy_VarCoords(hadi(0:nlags-1,:,:),y_Lead_x_tval)
printVarSummary(y_Lead_x_tval)
printMinMax(y_Lead_x_tval,True)
;---------------------------------------------
copy_VarCoords(hadi(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
;
d_hadi = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(hadi),hadi@_FillValue)
printVarSummary(d_hadi)
d_hadi(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_hadi(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_hadi!0 = "lag"
;d_hadi&lag = ispan(-nlags,nlags-1,1)
;print(d_hadi&lag)
printVarSummary(d_hadi)
;
d_hadi_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(hadi),hadi@_FillValue)		;re-arrange the t-values
printVarSummary(d_hadi)
d_hadi_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_hadi_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_hadi_tval!0 = "lag"
;d_hadi_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_hadi_tval&lag)
printVarSummary(d_hadi_tval)
;----------------------------------
x = conform_dims(dimsizes(d_hadi_tval),df,(/1,2/)) 
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_hadi_tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+d_hadi_tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_hadi_d = d_hadi	;transfers metadata
d_hadi_d = where(abs(prob).lt.p95,d_hadi,d_hadi_d@_FillValue) ; Significant TEMP
printVarSummary(d_hadi_d)
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_hadi_tval)
delete(b)
delete(prob)
delete(hadi)
delete(x_Lead_y_var)
delete(y_Lead_x_var)

























;=====================================================================
;
df = qnet(0,:,:)
lag_ts_PC = esacr(PC(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(qnet(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(qnet(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag))))
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================

;=================================================================
;
dimp   = dimsizes(qnet)
;
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PC(0:dimsizes(PC)-ii-1),qnet(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_var(ii,:,:) = 100.*((escorc_n(PC(0:dimsizes(PC)-ii-1),qnet(ii:,:,:),0,0))^2)     ;
print(ii)
;printMinMax(x_Lead_y@tval,True)
;tval = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
;print("==================")
;printMinMax(tval,True)
;printVarSummary(tval)
end do
printVarSummary(x_Lead_y_tval)
delete(ii)
copy_VarCoords(qnet(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(qnet(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)
;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
printVarSummary(y_Lead_x)
y_Lead_x_tval = y_Lead_x

do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PC(ii:),qnet(0:dimsizes(PC)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_var(ii,:,:) = 100.*((escorc_n(PC(ii:),qnet(0:dimsizes(PC)-ii-1,:,:),0,0))^2)     ;
print(ii)
;printVarSummary(y_Lead_x(ii,:,:))
;printMinMax(y_Lead_x@tval,True)
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))	;takes out the t-values
end do
printVarSummary(y_Lead_x)
delete(ii)
copy_VarCoords(qnet(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;
copy_VarCoords(qnet(0:nlags-1,:,:),y_Lead_x_tval)
printVarSummary(y_Lead_x_tval)
printMinMax(y_Lead_x_tval,True)
;---------------------------------------------
copy_VarCoords(qnet(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
;
d_qnet = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(qnet),qnet@_FillValue)
printVarSummary(d_qnet)
d_qnet(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_qnet(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_qnet!0 = "lag"
;d_qnet&lag = ispan(-nlags,nlags-1,1)
;print(d_qnet&lag)
printVarSummary(d_qnet)
;
d_qnet_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(qnet),qnet@_FillValue)		;re-arrange the t-values
printVarSummary(d_qnet)
d_qnet_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_qnet_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_qnet_tval!0 = "lag"
;d_qnet_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_qnet_tval&lag)
printVarSummary(d_qnet_tval)
;----------------------------------
x_qnet = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(qnet),qnet@_FillValue)             ;re-arrange the t-values
printVarSummary(x_qnet)
x_qnet(0:nlags-1,:,:) = y_Lead_x_var(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
x_qnet(nlags:,:,:)    = x_Lead_y_var(1:nlags-1,:,:)     ; "positive lag"
x_qnet!0 = "lag"
copy_VarCoords(d_qnet,x_qnet) 
printVarSummary(x_qnet)
;----------------------------------
x = conform_dims(dimsizes(d_qnet_tval),df,(/1,2/)) 
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_qnet_tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+d_qnet_tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_qnet_d = d_qnet	;transfers metadata
d_qnet_d = where(abs(prob).lt.p95,d_qnet,d_qnet_d@_FillValue) ; Significant TEMP
printVarSummary(d_qnet_d)
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_qnet_tval)
delete(b)
delete(prob)
delete(qnet)
delete(x_Lead_y_var)
delete(y_Lead_x_var)
;==============================================================
;=================================================================
;=====================================================================


















;PC = month_to_annual_weighted(yyyymm,pc,1)
PCs = month_to_season(pc,season_w)		;uses DJF PC for slp
delete(pc)
PCs!0 = "time"


;=====================================================================
;
df = slp(0,:,:)
lag_ts_PC = esacr(PCs(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(slp(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(slp(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag))))
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================

;=================================================================
;
dimp   = dimsizes(slp)
;
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PCs(0:dimsizes(PC)-ii-1),slp(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_var(ii,:,:) = 100.*((escorc_n(PCs(0:dimsizes(PC)-ii-1),slp(ii:,:,:),0,0))^2)     ;
print(ii)
;printMinMax(x_Lead_y@tval,True)
;tval = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
;print("==================")
;printMinMax(tval,True)
;printVarSummary(tval)
end do
printVarSummary(x_Lead_y_tval)
delete(ii)
copy_VarCoords(slp(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(slp(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)
;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
printVarSummary(y_Lead_x)
y_Lead_x_tval = y_Lead_x

do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PCs(ii:),slp(0:dimsizes(PCs)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_var(ii,:,:) = 100.*((escorc_n(PCs(ii:),slp(0:dimsizes(PCs)-ii-1,:,:),0,0))^2)     ;
print(ii)
;printVarSummary(y_Lead_x(ii,:,:))
;printMinMax(y_Lead_x@tval,True)
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))	;takes out the t-values
end do
printVarSummary(y_Lead_x)
delete(ii)
copy_VarCoords(slp(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;
copy_VarCoords(slp(0:nlags-1,:,:),y_Lead_x_tval)
printVarSummary(y_Lead_x_tval)
printMinMax(y_Lead_x_tval,True)
;---------------------------------------------
copy_VarCoords(slp(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
;
d_slp = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(slp),slp@_FillValue)
printVarSummary(d_slp)
d_slp(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_slp(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_slp!0 = "lag"
;d_slp&lag = ispan(-nlags,nlags-1,1)
;print(d_slp&lag)
printVarSummary(d_slp)
;
d_slp_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(slp),slp@_FillValue)		;re-arrange the t-values
printVarSummary(d_slp)
d_slp_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_slp_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_slp_tval!0 = "lag"
;d_slp_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_slp_tval&lag)
printVarSummary(d_slp_tval)
;----------------------------------
x = conform_dims(dimsizes(d_slp_tval),df,(/1,2/)) 
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_slp_tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+d_slp_tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_slp_d = d_slp	;transfers metadata
d_slp_d = where(abs(prob).lt.p95,d_slp,d_slp_d@_FillValue) ; Significant TEMP
printVarSummary(d_slp_d)
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_slp_tval)
delete(b)
delete(prob)
delete(slp)
delete(x_Lead_y_var)
delete(y_Lead_x_var)

delete(PCs)
;==============================================================















;=====================================================================
;
df = ssh(0,:,:)
lag_ts_PC = esacr(PC(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(ssh(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(ssh(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag))))
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================

;=================================================================
;
dimp   = dimsizes(ssh)
;
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PC(0:dimsizes(PC)-ii-1),ssh(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_var(ii,:,:) = 100.*((escorc_n(PC(0:dimsizes(PC)-ii-1),ssh(ii:,:,:),0,0))^2)     ;
print(ii)
;printMinMax(x_Lead_y@tval,True)
;tval = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
;print("==================")
;printMinMax(tval,True)
;printVarSummary(tval)
end do
printVarSummary(x_Lead_y_tval)
delete(ii)
copy_VarCoords(ssh(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(ssh(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)
;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
printVarSummary(y_Lead_x)
y_Lead_x_tval = y_Lead_x

do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PC(ii:),ssh(0:dimsizes(PC)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_var(ii,:,:) = 100.*((escorc_n(PC(ii:),ssh(0:dimsizes(PC)-ii-1,:,:),0,0))^2)     ;
print(ii)
;printVarSummary(y_Lead_x(ii,:,:))
;printMinMax(y_Lead_x@tval,True)
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))	;takes out the t-values
end do
printVarSummary(y_Lead_x)
delete(ii)
copy_VarCoords(ssh(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;
copy_VarCoords(ssh(0:nlags-1,:,:),y_Lead_x_tval)
printVarSummary(y_Lead_x_tval)
printMinMax(y_Lead_x_tval,True)
;---------------------------------------------
copy_VarCoords(ssh(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
;
d_ssh = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(ssh),ssh@_FillValue)
printVarSummary(d_ssh)
d_ssh(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_ssh(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_ssh!0 = "lag"
;d_ssh&lag = ispan(-nlags,nlags-1,1)
;print(d_ssh&lag)
printVarSummary(d_ssh)
;
d_ssh_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(ssh),ssh@_FillValue)		;re-arrange the t-values
printVarSummary(d_ssh)
d_ssh_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_ssh_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_ssh_tval!0 = "lag"
;
x_ssh = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(ssh),ssh@_FillValue)             ;re-arrange the t-values
printVarSummary(x_ssh)
x_ssh(0:nlags-1,:,:) = y_Lead_x_var(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
x_ssh(nlags:,:,:)    = x_Lead_y_var(1:nlags-1,:,:)     ; "positive lag"
x_ssh!0 = "lag"
copy_VarCoords(d_ssh,x_ssh) 
;d_ssh_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_ssh_tval&lag)
printVarSummary(x_ssh)
;----------------------------------
x = conform_dims(dimsizes(d_ssh_tval),df,(/1,2/)) 
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_ssh_tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+d_ssh_tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_ssh_d = d_ssh	;transfers metadata
d_ssh_d = where(abs(prob).lt.p95,d_ssh,d_ssh_d@_FillValue) ; Significant TEMP
printVarSummary(d_ssh_d)
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_ssh_tval)
delete(b)
delete(prob)
delete(ssh)
delete(x_Lead_y_var)
delete(y_Lead_x_var)
;==============================================================



















;=====================================================================
;
df = u(0,:,:)
lag_ts_PC = esacr(PC(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(u(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(u(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag))))
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================
;

;
dimp   = dimsizes(u)			;U-wind stress
;=================================================================
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PC(0:dimsizes(PC)-ii-1),u(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
print(ii)
end do
delete(ii)
copy_VarCoords(u(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(u(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)

;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
y_Lead_x_tval = y_Lead_x
printVarSummary(y_Lead_x)
do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PC(ii:),u(0:dimsizes(PC)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))       ;takes out the t-values
print(ii)
end do
delete(ii)
copy_VarCoords(u(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;---------------------------------------------
copy_VarCoords(u(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
d_u = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(u),u@_FillValue)
printVarSummary(d_u)
d_u(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_u(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_u!0 = "lag"
;d_u&lag = ispan(-nlags,nlags-1,1)
;print(d_u&lag)
printVarSummary(d_u)
;
;
d_u_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(u),u@_FillValue)          ;re-arrange the t-values
printVarSummary(d_u)
d_u_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_u_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_u_tval!0 = "lag"
;d_hadi_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_hadi_tval&lag)
printVarSummary(d_u_tval)
;----------------------------------
;----------------------------------
x = conform_dims(dimsizes(d_u_tval),df,(/1,2/))
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_u_tval    ; b must be same size as tval (and df)
b = 0.5
probU = betainc(df/(df+d_u_tval^2),df/2.0,b)       ; prob(nlat,nlon)
probU@long_name   = "regression coefficient"
probU@long_name = "probability"
printMinMax(probU,True)
printVarSummary(probU)
;
;----------------------------------
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_u_tval)
delete(b)
;delete(prob)
delete(u)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;
df = v(0,:,:)
lag_ts_PC = esacr(PC(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_PC)
lag_ts_HADI   = esacr(v(lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
printVarSummary(lag_ts_HADI)
df = floattointeger(dimsizes(v(:,0,0))*(1. -(lag_ts_PC(lag)*lag_ts_HADI(:,:,lag)))\
        /(1. + (lag_ts_HADI(:,:,lag)*lag_ts_PC(lag))))
printVarSummary(df)
printMinMax(df,True)
;
delete(lag_ts_PC)
delete(lag_ts_HADI)
;=================================================================
;
;
dimp   = dimsizes(v)                    ;V-wind stress
;
x_Lead_y  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
x_Lead_y_var = x_Lead_y
x_Lead_y_tval = x_Lead_y
printVarSummary(x_Lead_y)
do ii=0,nlags-1,1
x_Lead_y(ii,:,:) = regCoef_n(PC(0:dimsizes(PC)-ii-1),v(ii:,:,:),0,0)          ; PC leads map
x_Lead_y_tval(ii,:,:) = onedtond(x_Lead_y@tval,dimsizes(x_Lead_y(0,:,:)))
print(ii)
end do
delete(ii)
copy_VarCoords(v(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
printMinMax(x_Lead_y,True)
;
copy_VarCoords(v(0:nlags-1,:,:),x_Lead_y_tval)
printVarSummary(x_Lead_y_tval)
printMinMax(x_Lead_y_tval,True)
;---------------------------------------------
y_Lead_x  =  new((/nlags,dimp(1),dimp(2)/),typeof(PC))
y_Lead_x_var = y_Lead_x
y_Lead_x_tval = y_Lead_x
printVarSummary(y_Lead_x)
do ii=0,nlags-1,1
y_Lead_x(ii,:,:) = regCoef_n(PC(ii:),v(0:dimsizes(PC)-ii-1,:,:),0,0)          ; Map leads PC
y_Lead_x_tval(ii,:,:) = onedtond(y_Lead_x@tval,dimsizes(y_Lead_x(0,:,:)))       ;takes out the t-values
print(ii)
end do
delete(ii)
copy_VarCoords(v(0:nlags-1,:,:),y_Lead_x)
printVarSummary(y_Lead_x)
printMinMax(y_Lead_x,True)
;---------------------------------------------
copy_VarCoords(v(0:nlags-1,:,:),x_Lead_y)
printVarSummary(x_Lead_y)
d_v = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(v),v@_FillValue)
printVarSummary(d_v)
d_v(0:nlags-1,:,:) = y_Lead_x(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order 
d_v(nlags:,:,:)    = x_Lead_y(1:nlags-1,:,:)     ; "positive lag"
d_v!0 = "lag"
;d_v&lag = ispan(-nlags,nlags-1,1)
;print(d_v&lag)
printVarSummary(d_v)
;
;
;
d_v_tval = new((/2*nlags-1,dimp(1),dimp(2)/), typeof(v),v@_FillValue)          ;re-arrange the t-values
printVarSummary(d_v)
d_v_tval(0:nlags-1,:,:) = y_Lead_x_tval(0:nlags-1:-1,:,:)  ; "negative lag", -1 reverses order
d_v_tval(nlags:,:,:)    = x_Lead_y_tval(1:nlags-1,:,:)     ; "positive lag"
d_v_tval!0 = "lag"
;d_hadi_tval&lag = ispan(-nlags,nlags-1,1)
;print(d_hadi_tval&lag)
printVarSummary(d_v_tval)
;----------------------------------
;----------------------------------
x = conform_dims(dimsizes(d_v_tval),df,(/1,2/))
printVarSummary(x)
delete(df)
df = x
delete(x)
printVarSummary(df)
b = d_v_tval    ; b must be same size as tval (and df)
b = 0.5
probV = betainc(df/(df+d_v_tval^2),df/2.0,b)       ; prob(nlat,nlon)
probV@long_name   = "regression coefficient"
probV@long_name = "probability"
printMinMax(probV,True)
printVarSummary(probV)
;
;----------------------------------
;----------------------------------
delete(y_Lead_x)
delete(x_Lead_y)
delete(y_Lead_x_tval)
delete(x_Lead_y_tval)
delete(dimp)
delete(df)
delete(d_v_tval)
delete(b)
;delete(prob)
delete(v)
;==============================================================
;==============================================================
;==============================================================
;
;
;
print("==============================================")
printVarSummary(probV)
printVarSummary(probU)
printVarSummary(d_u)
printVarSummary(d_v)
print("==============================================")
d_u  = where((probU.lt.p95.or.probV.lt.p95),d_u,d_u@_FillValue) 
d_v  = where((probV.lt.p95.or.probU.lt.p95),d_v,d_v@_FillValue)
;
printMinMax(d_u,True)
printMinMax(d_v,True) 
;==============================================================
;============================================================
; PLOTS
;============================================================
  npanels = 2*mxlag-1    ;number of panels

  wks = gsn_open_wks("eps","fig0_4")         ; send graphics to PNG file

  plotM = new(npanels,graphic)                ;
  plotX2 = new(npanels,graphic)
  plotX3 = new(npanels,graphic)
  plot = new(npanels,graphic)                ; create graphic array ; only needed if paneling
  plotS = new(npanels,graphic)                ;

  plot1 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
  plotS1 = new(npanels,graphic)                ;

  plot2 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
  plotS2 = new(npanels,graphic)                ;

  
  plot3 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
  plotS3 = new(npanels,graphic)                ;


; EOF patterns

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -.3            ; set min contour level
  res@cnMaxLevelValF       =  .3             ; set max contour level
  res@cnLevelSpacingF      =  0.03            ; set contour spacing     
 ; res@gsnPolar             = "SH"
 ; res@mpProjection      = "LambertConformal"       ; choose projection
 
  res@mpFillOn             = True        ; turn off map fill
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR
  ;res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.030
  res@tmXBLabelFontHeightF = 0.030
  res@pmTickMarkDisplayMode = "Always"
  res@tmXBMajorThicknessF   =  1.0
  res@tmYLMajorThicknessF   =  1.0
  res@tmBorderThicknessF    =  1.0
  ;res@tmXBLabelsOn          = False      ;no lon ticks
  ;res@tmYLMode   = "Explicit"
  ;res@tmYLValues = ispan(-40,60,20)
  ;res@tmYLLabels = (/"40~S~o~N~S","20~S~o~N~S", "0~S~o~N","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N"/)
  ;res@tmYLMinorValues  = ispan(-50,70,10) 
  ;res@tmYLMinorValues = res@tmYLValues
  ;res@tmYLLabelFontHeightF = 0.015

  ;res@tmXBMode   = "Explicit"
  ;res@tmXBValues = ispan(-90,30,30)
  ;res@tmXBLabels = (/"90~S~o~N~W", "60~S~o~N~W","30~S~o~N~W","0~S~o~N", \
  ;                   "30~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-90,30,15)
  ;res@tmXBMinorValues  = res@tmXBValues
  res@mpGeophysicalLineColor  = "transparent" 
  res@mpFillDrawOrder      = "preDraw"

  res@gsnAddCyclic         = True
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  
  res@lbLabelBarOn         = False     ; turn off individual lb's
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF   = 0.020 
  res@lbLabelFontHeightF   = 0.020
  res@lbLabelFont          = "Helvetica-Bold"     ; label font
  res@lbTitleString        = ""
  ;res@lbTitleAngleF       = 90.
  ;res@lbTitlePosition     = "Right"
  ;res@lbTitleDirection    = "Across"
  res@lbTopMarginF         = 0.40
  res@lbBottomMarginF      = 0.35

  ;res@vpWidthF            = 0.40
  ;res@vpHeightF           = 0.80
  
  res@lbLabelStride        = 10
  res@lbTitleString         = "Ocean temperature (K)"
  
  res@gsnCenterString  	      = ""   
  res@gsnLeftString           = ""  
  res@gsnRightString          = ""

  ;res@tmXTBorderOn            =  False
  ;res@tmXTOn                  =  False
  ;res@tmYRBorderOn            =  False
  ;res@tmYROn                  =  False

  ;res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  
  ;res@tmXBMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmYLMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmXTMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmYRMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  

  ;res@gsnMaskLambertConformal = True            ; turn on lc masking

  ;res@mpGridAndLimbOn         = True  
  ;res@mpGridLatSpacingF       = 10
  ;res@mpGridLonSpacingF       =  100

  res@vcRefAnnoOrthogonalPosF   = -0.35             ; move ref vector up
  res@vcRefAnnoParallelPosF	= 0.45
  res@vcRefMagnitudeF           = 0.005            ; define vector ref ma
  ;res@vcRefMagnitudeF           = 0.005            ; define vector ref mag
  res@vcRefAnnoString1          = "0.005 N m~S~-2~N~ "
  ;res@vcRefAnnoString1          = "0.005 N m~S~-2~N~ "
  ;res@vcRefAnnoString1		= "~F10~e~S~~F21~-4~N~ N m~S~-2~N~"
  res@vcRefLengthF            = 0.045            ; define length of vec ref
 ;; res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.035
  res@vcLineArrowColor        = "grey25"           ; change vector color
  res@vcRefAnnoBackgroundColor = "transparent"
  res@vcLineArrowThicknessF   = 1.0
  res@vcLineArrowHeadMinSizeF = 0.015
  ;res@vcRefAnnoPerimColor     = 1
  res@vcRefAnnoPerimOn		= False
  res@gsnScalarContour	      = True
  
;======================================================
;==============================================================
;-----------------------------------------------------------------------------
  res0                      =   True
  res0@gsnDraw              = False        ; don't draw yet
  res0@gsnFrame             = False        ; don't advance frame yet
  res0@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res0@cnMinLevelValF       =  -0.10            ; set min contour level
  res0@cnMaxLevelValF       =  0.10             ; set max contour level
  res0@cnLevelSpacingF      =  0.01            ; set contour spacing  
  res0@gsnAddCyclic         = True   
  res0@cnLinesOn            = True       ; True is default
  res0@cnLineLabelsOn       = True        ; True is default
  res0@cnInfoLabelOn        = False
  res0@gsnRightString   = ""
  res0@gsnLeftString   = ""
  res0@cnLineLabelDensityF =  1
  res0@gsnContourZeroLineThicknessF = 0.
  res0@cnLineLabelBackgroundColor   = "transparent"
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res2                      =   True
  res2@gsnDraw              = False        ; don't draw yet
  res2@gsnFrame             = False        ; don't advance frame yet
  ;res2@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  ;res2@cnMinLevelValF       = 900.            ; set min contour level
  ;res2@cnMaxLevelValF       =  1500.             ; set max contour level
  ;res2@cnLevelSpacingF      =  3.            ; set contour spacing  
  res2@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res2@cnLevels    = (/1009.,1019./)   ; set levels
  res2@cnMonoLineDashPattern  = False
  res2@cnLineDashPatterns     = (/16,0/)
  ;res2@gsnContourNegLineDashPattern  = 1
  res2@cnLineDashPattern  = 2           ;dashed low-pressure cells    
  res2@gsnAddCyclic         = True
  res2@cnLinesOn            = True       ; True is default
  res2@cnLineColor          = "white"       ; True is default
  res2@cnLineLabelsOn       = False        ; True is default
  res2@cnInfoLabelOn        = False
  res2@cnLineLabelPerimColor = "transparent"
  res2@cnLineLabelFontHeightF = 0.010
  res2@cnLineThicknessF       = 1.5
  res2@gsnRightString         = ""
  res2@gsnLeftString          = ""
  res2@cnLabelMasking         = True
  res2@cnLineLabelDensityF    =  2
  ;res2@gsnContourNegLineDashPattern  = 1
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  ;res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
 ; res1@cnMinLevelValF       =  -1.0            ; set min contour level
 ; res1@cnMaxLevelValF       =  1.0             ; set max contour level
 ; res1@cnLevelSpacingF      =  .05            ; set contour spacing  
  res1@gsnAddCyclic         = True   
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString   = ""
  res1@gsnLeftString   = ""
  res1@cnLineLabelDensityF =  5
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #17
  opt@gsnShadeHigh = 17                     ; use fill pattern #17
  opt@gsnShadeFillScaleF = 0.6
  opt@gsnShadeFillDotSizeF = 0.001
  ;

;*******************************************
; first plot
;*******************************************
    res@lbTitleString         = "SST anomaly (K)"
    res@gsnStringFont		= "Helvetica-Bold"     ; label font
    res@gsnCenterStringFontHeightF   = 0.015
 
    do n=0,npanels-1,1
     ;do n=0,11,3
    ;res@gsnLeftString  = "EOF "+(n+1)
   ;res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
    ;res@gsnLeftString  = "lag" +  ccr&lag(n)
 ;   res@gsnLeftString  = DATA_NAMES(n)
    ;if(n.gt.npanels/2)
   ;   res@gsnCenterString  = "lag +"+labels(n)
   ;else
   ;  res@gsnCenterString  = "lag "+labels(n)
   ;end if
;if (n.ge.npanels/2) then
;   res@tmXBLabelsOn          =  True
;else
;   res@tmXBLabelsOn          =  False
;end if
;---
;       if (n.eq.0)then
;   res@tmYLLabelsOn          =  True
;       else if(n.eq.6)
;  res@tmYLLabelsOn          =  True
;       else
;   res@tmYLLabelsOn          =  False
;       end if
;       end if


  ;     if (n.eq.0)then
  ; res@tmYLLabelsOn          =  True
  ;     else if(n.gt.0)
  ;res@tmYLLabelsOn          =  False
  ;     end if
  ;     end if

  ; plot(n) = gsn_csm_vector_scalar_map_ce(wks,d_u(n+2,:,:),d_v(n+2,:,:),d_hadi(n+2,:,:),res)  ;               ; gsn_csm_contour_map(wks,r_ssh(lat|:,lon|:,lag|n),res)      ;ssh/sea-level index
  ; plotS(n) = gsn_csm_contour(wks,d_hadi_d(lat|:,lon|:,lag|n+2),res1)       ;ssh/sea-level index
  ; plotS(n) = gsn_contour_shade(plotS,p95,d_hadi_d@_FillValue,opt)
  ; plotQ(n) = gsn_csm_contour(wks,d_qnet(lat|:,lon|:,lag|n+2),res2)

	   
    print(n)
   plot(n) = gsn_csm_vector_scalar_map_ce(wks,d_u(n,:,:),d_v(n,:,:),d_hadi(n,:,:),res)  ; 		; gsn_csm_contour_map(wks,r_ssh(lat|:,lon|:,lag|n),res)      ;ssh/sea-level index
   plotS(n) = gsn_csm_contour(wks,d_hadi_d(lat|:,lon|:,lag|n),res1)       ;ssh/sea-level index
   plotS(n) = gsn_contour_shade(plotS(n),p95,d_hadi_d@_FillValue,opt)
   ;plotQ(n) = gsn_csm_contour(wks,d_qnet(lat|:,lon|:,lag|n),res2)
   overlay(plot(n),plotS(n))
   ;overlay(plot(n),plotQ(n))
   plot(n) = plot(n)
end do

;delete(n)
;=========================================================================
;-----------------------------------------------------------------------------
  delete(res@vcRefAnnoOrthogonalPosF)
  delete(res@vcRefAnnoParallelPosF)
  delete(res@vcRefMagnitudeF)
  delete(res@vcRefAnnoString1)
  delete(res@vcRefLengthF)
  delete(res@vcMinDistanceF)
  delete(res@vcLineArrowColor)
  delete(res@vcRefAnnoBackgroundColor)
  delete(res@vcLineArrowThicknessF)
  delete(res@vcLineArrowHeadMinSizeF)
  delete(res@vcRefAnnoPerimOn)
  delete(res@gsnScalarContour)
;*******************************************
;*******************************************
  delete(res@cnFillPalette)
  res@cnFillPalette        = "GMT_panoply"                        ;different colour bar  
  res@cnMinLevelValF       = -0.5            ; set min contour level
  res@cnMaxLevelValF       =  0.5             ; set max contour level
  res@cnLevelSpacingF      =  .05
  res@mpGeophysicalLineColor  = "foreground"


  

  ;d_slp_d(lat|:,lon|:,lag|12) =  d_slp_d(lat|:,lon|:,lag|15)
  ;d_slp(lat|:,lon|:,lag|11) =  d_slp(lat|:,lon|:,lag|14)

do n=0,npanels-1,1

   print(n)

   ;       if (n.eq.11)then
  ;         d_slp_d(lat|:,lon|:,lag|n) = 0.
;   res@tmYLLabelsOn          =  True
;       else if(n.eq.6)
;  res@tmYLLabelsOn          =  True
;       else
;   res@tmYLLabelsOn          =  False
    ;   end if
;       end if
   plot1(n) = gsn_csm_contour_map_ce(wks,d_slp(n,:,:),res)  ;           ; gsn_csm_contour_map(wks,r_ssh(lat|:,lon|:,lag|n),res)      ;ssh/sea-level index
   plotS1(n) = gsn_csm_contour(wks,d_slp_d(lat|:,lon|:,lag|n),res1)       ;ssh/sea-level index
   plotS1(n) = gsn_contour_shade(plotS1(n),d_slp_d@_FillValue,p95,opt)
   plotM(n) = gsn_csm_contour(wks,m_slp,res2)  ;
  
   overlay(plot1(n),plotS1(n))
   overlay(plot1(n),plotM(n))
   plot1(n) = plot1(n)
end do

delete(n)

print(d_slp_d(lat|:,lon|:,lag|11))
;=======================================================================
;======================================================================= 
 delete(res@cnFillPalette)
 delete(res2@cnLevels)
 delete(res2@cnLineDashPatterns)
  ;res2@cnLevels    = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100./)   ; set levels 
  res2@cnLevels    = (/3.,6.,9.,12.,15.,18.,21.,24.,27.,30./)   ; set levels 
  ;res2@cnLevels    = (/5.,10.,15.,20.,25.,30.,35.,40.,45.,50./)   ; set levels 
  res2@cnLineDashPattern  = 0           ;dashed low-pressure cells   
  res2@cnLineDashPatterns     =    (/0,0,0,0,0,0,0,0,0,0/)
  res2@cnLineThicknessF       =    0.5

 ;resQ@cnFillPalette        = "temp_19lev"       ;different colour bar 
  res@cnFillPalette        = "amwg_blueyellowred"       ;different colour bar
  res@cnMinLevelValF       = -10.            ; set min contour level
  res@cnMaxLevelValF       =  10.             ; set max contour level
  res@cnLevelSpacingF      =  1.
  res@mpGeophysicalLineColor  = "transparent"
  res@mpFillDrawOrder      = "postDraw"



do n=0,npanels-1,1

   print(n)
   plot2(n) = gsn_csm_contour_map_ce(wks,d_qnet(n,:,:),res)  ;           ; 
   plotS2(n) = gsn_csm_contour(wks,d_qnet_d(lat|:,lon|:,lag|n),res1)       ;
   plotS2(n) = gsn_contour_shade(plotS2(n),d_qnet_d@_FillValue,p95,opt)
   ;plotX2(n) = gsn_csm_contour(wks,x_qnet(n,:,:),res2)  ;
   overlay(plot2(n),plotS2(n))
   ;overlay(plot2(n),plotX2(n))
   plot2(n) = plot2(n)
end do

delete(n)
;=======================================================================
;;===============================================================
  delete(res@cnFillPalette)
  ;res@cnFillPalette        = "ncl_default"       ;different colour bar 
  res@cnFillPalette        = "amwg_blueyellowred"       ;different colour bar  
  res@cnMinLevelValF       = -2.            ; set min contour level
  res@cnMaxLevelValF       =  2.             ; set max contour level
  res@cnLevelSpacingF      =  .2
  res@mpGeophysicalLineColor  = "transparent"
  res@mpFillDrawOrder      = "postDraw"

do n=0,npanels-1,1

   print(n)
   plot3(n) = gsn_csm_contour_map_ce(wks,d_ssh(n,:,:),res)  ;           ; gsn_csm_contour_map(wks,r_ssh(lat|:,lon|:,lag|n),res)      ;ssh/sea-level index
   plotS3(n) = gsn_csm_contour(wks,d_ssh_d(lat|:,lon|:,lag|n),res1)       ;ssh/sea-level index
   plotS3(n) = gsn_contour_shade(plotS3(n),d_ssh_d@_FillValue,p95,opt)
   plotX3(n) = gsn_csm_contour(wks,x_ssh(n,:,:),res2)  ;
   overlay(plot3(n),plotS3(n))
   ;overlay(plot3(n),plotX3(n))
   plot3(n) = plot3(n)
end do

delete(n)
;-----------------------------------------------------------------------------
;=============================================================================
;=========================================================================
 pc = ispan(0,npanels-1,2)
 print(pc)




 ;=========================================================================
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise
  txres@txFontHeightF = 0.070   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "DT anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres);
; Move text string to center/right edge of plot.;
  txres                  = True
  txres@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so
  txres@amOrthogonalPosF =  -0.65    ; This is the center of the plot.
  txres@amJust           = "TopLeft"

 tt0 = (/"a","","b","","c","","d","","e","","f","","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z" /)
 create_tt0 = gsn_create_text(wks,tt0,txres);
 tt_p0 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
;
 tt1 = (/"h","","i","","j","","k","","l","","m","","n","","o","","p","q","r","s","t","u","v","w","x","y","z" /)
 create_tt1 = gsn_create_text(wks,tt1,txres);
 tt_p1 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
;
;
 ;
 tt2 = (/"o","","p","","q","","r","","s","","t","","u","","v","","w","","x","y","z" /)
 create_tt2 = gsn_create_text(wks,tt2,txres);
 tt_p2 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
;

 ;
 tt3 = (/"s","","","t","","","u","","","v","","","w","","","x","","","y","z" /)
 create_tt3 = gsn_create_text(wks,tt3,txres);
 tt_p3 = new(npanels,graphic)                ; create graphic array ; only needed if paneling
;

;---Attach the polyline
     do i = 0,npanels-1,2
	     print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt0(i),txres)
     end do
 delete(i)
;
;
;---Attach the polyline
     do i = 0,npanels-1,2
             print(i)
        tt_p1(i) = gsn_add_annotation(plot3(i),create_tt1(i),txres)
     end do
 delete(i)
;
;---Attach the polyline
     do i = 0,npanels-1,2
             print(i)
        tt_p2(i) = gsn_add_annotation(plot1(i),create_tt2(i),txres)
     end do
 delete(i)
;
;
;---Attach the polyline
;     do i = 0,npanels-1,3
;             print(i)
;        tt_p3(i) = gsn_add_annotation(plot1(i),create_tt3(i),txres)
;     end do
; delete(i)
;
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;---------------------------------------------------------------------
var_names = (/"Atmospheric and oceanic fields regressed on the second principal component during December-February"/)
printVarSummary(var_names)
print(var_names)

  txres1               = True         ; Text resources
  txres1@txFontHeightF = 0.012
  txres1@txJust        = "CenterLeft"
  txres1@txFont                = "Helvetica"
  txres1@txAngleF      =  0.

  labels0  = var_names
  ;xpos0    = fspan(0.150, 0.60,dimsizes(var_names))

;     gsn_text_ndc(wks,var_names,0.15,0.98,txres1)
  
delete(txres1)
;=====================================================================
;=====================================================================


;=========================================================================
;=========================================================================
;=======================================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True       ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = True        ; add common colorbar
  resP@lbLabelAutoStride              	= True         ; auto stride on labels
  ;resP@lbOrientation                    = "vertical"
  resP@gsnPanelYWhiteSpacePercent      	= 2.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 2.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.95
  resP@gsnPanelBottom                   = 0.75
  resP@gsnPanelLeft                     = 0.005
  resP@gsnPanelRight                    = 0.995
  resP@lbTitleFontHeightF	        = 0.009
  resP@lbLabelFontHeightF               = 0.009
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
  resP@gsnPanelFigureStringsFontHeightF = 0.008
  resP@gsnPanelFigureStringsJust        = "BottomLeft"
  resP@gsnPanelFigureStrings            = (/"lag -6","lag -4","lag -2","lag 0",\
                                                "lag +2","lag +4","lag +4","lag +4","lag +5","lag +6"/)

  resP@gsnPanelCenter   		= False                 ; =True .... must move labels
  ;resP@gsnPanelRowSpec  		= True                  ; tell panel what order to plot
;========================================================================
gsn_panel(wks,(/plot(pc)/),(/1,7/),resP)     ; now draw as one plot
;gsn_panel(wks,(/plot(pc)/),(/1,7/),resP)     ; now draw as one plot
;gsn_panel(wks,(/plot(0:npanels-1)/),(/5,4/),resP)     ; now draw as one plot
;gsn_panel(wks,plotMAP(0:1),(/1,2/),resP)     ; now draw as one plot
;;===============================================================
;;=============================================================================
;=======================================================================
  res_P                                  = True         ; modify the panel plot
  res_P@gsnDraw                          = True       ; don't draw yet
  res_P@gsnFrame                         = False        ; don't advance frame yet
  res_P@gsnPanelLabelBar                 = True        ; add common colorbar
  res_P@lbLabelAutoStride                = True         ; auto stride on labels
  res_P@gsnPanelYWhiteSpacePercent       = 2.0    ; default is 1.0
  res_P@gsnPanelXWhiteSpacePercent       = 2.0    ; default is 1.0
  res_P@gsnPanelTop                      = 0.55
  res_P@gsnPanelBottom                   = 0.35
 ; res_P@gsnPanelTop                      = 0.35
 ; res_P@gsnPanelBottom                   = 0.15
  res_P@gsnPanelLeft                     = 0.005
  res_P@gsnPanelRight                    = 0.995
  res_P@lbTitleFontHeightF               = 0.009
  res_P@lbLabelFontHeightF               = 0.009
  res_P@lbLabelFont                      = "Helvetica-Bold"     ; label font
  res_P@lbTitleString                    = "Sea level pressure (hPa)"
  ;res_P@lbTitleAngleF                    = 90.
  res_P@lbLabelStride                    = 10
  res_P@lbTitleDirection                 = "Across"
  res_P@lbTopMarginF                     = 0.5
  res_P@gsnMaximize                      = True         ; large format
  res_P@gsnPaperOrientation              = "portrait"   ; force portrait
  res_P@gsnPanelFigureStringsFontHeightF = 0.008
  res_P@gsnPanelFigureStringsJust        = "BottomLeft"
  res_P@gsnPanelFigureStrings            = (/"lag -6","lag -4","lag -2","lag 0",\
                                                "lag +2","lag +4","lag +6","lag +4","lag +5","lag +6"/)
;========================================================================
 gsn_panel(wks,(/plot1(pc)/),(/1,7/),res_P)     ; now draw as one plot
;=========================================================================
;*******************************************
  res_T                                  = True         ; modify the panel plot
  res_T@gsnDraw                          = True       ; don't draw yet
  res_T@gsnFrame                         = False        ; don't advance frame yet
  res_T@gsnPanelLabelBar                 = True        ; add common colorbar
  res_T@lbLabelAutoStride                = True         ; auto stride on labels
  res_T@gsnPanelYWhiteSpacePercent       = 2.0    ; default is 1.0
  res_T@gsnPanelXWhiteSpacePercent       = 2.0    ; default is 1.0
  res_T@gsnPanelTop                      = 0.55
  res_T@gsnPanelBottom                   = 0.35
  res_T@gsnPanelLeft                     = 0.005
  res_T@gsnPanelRight                    = 0.995
  res_T@lbTitleFontHeightF               = 0.008
  res_T@lbLabelFontHeightF               = 0.008
  res_T@lbLabelFont                      = "Helvetica-Bold"     ; label font
  res_T@lbTitleString                    = "Surface net heat flux (W m ~S~-2~N~)"
  ;res_T@lbTitleAngleF                    = 90.
  res_T@lbLabelStride                    = 10
  res_T@lbTitleDirection                 = "Across"
  res_T@lbTopMarginF                     = 0.5
  res_T@gsnMaximize                      = True         ; large format
  res_T@gsnPaperOrientation              = "portrait"   ; force portrait
  res_T@gsnPanelFigureStringsFontHeightF = 0.007
  res_T@gsnPanelFigureStringsJust        = "BottomLeft"
  res_T@gsnPanelFigureStrings            = (/"lag -6","lag -4","lag -2","lag 0",\
                                                "lag +2","lag +4","lag +6","lag +4","lag +5","lag +6"/)
  
;========================================================================
; gsn_panel(wks,(/plot2(pc)/),(/1,6/),res_T)     ; now draw as one plot
;=========================================================================
;=========================================================================
;*******************************************
  res_H                                  = True         ; modify the panel plot
  res_H@gsnDraw                          = True       ; don't draw yet
  res_H@gsnFrame                         = False        ; don't advance frame yet
  res_H@gsnPanelLabelBar                 = True        ; add common colorbar
  res_H@lbLabelAutoStride                = True         ; auto stride on labels
  res_H@gsnPanelYWhiteSpacePercent       = 2.0    ; default is 1.0
  res_H@gsnPanelXWhiteSpacePercent       = 2.0    ; default is 1.0
  res_H@gsnPanelTop                      = 0.75
  res_H@gsnPanelBottom                   = 0.55
  res_H@gsnPanelLeft                     = 0.005
  res_H@gsnPanelRight                    = 0.995
  res_H@lbTitleFontHeightF               = 0.009
  res_H@lbLabelFontHeightF               = 0.009
  res_H@lbLabelFont                      = "Helvetica-Bold"     ; label font
  res_H@lbTitleString                    = "Sea surface height (10~S~-2~N~m)"
  ;res_H@lbTitleString                    = "Sea surface height (~F10~e~S~~F21~-2~N~m)"
  ;r0@tiYAxisString = "Zonal wind stress (~F10~e~S~~F21~-3~N~ N m~S~-2~N~)"                  ; yaxis
  ;res_H@lbTitleAngleF                    = 90.
  res_H@lbLabelStride                    = 10
  res_H@lbTitleDirection                 = "Across"
  res_H@lbTopMarginF                     = 0.5
  res_H@gsnMaximize                      = True         ; large format
  res_H@gsnPaperOrientation              = "portrait"   ; force portrait
  res_H@gsnPanelFigureStringsFontHeightF = 0.008
  res_H@gsnPanelFigureStringsJust        = "BottomLeft"
  res_H@gsnPanelFigureStrings            = (/"lag -6","lag -4","lag -2","lag 0",\
 						"lag +2","lag +4","lag +6","lag +4","lag +5","lag +6"/)
 ;========================================================================
 gsn_panel(wks,(/plot3(pc)/),(/1,7/),res_H)     ; now draw as one plot
;=========================================================================

frame(wks)
print(xxxxxxxxxxxxxxx)
;;========================================================================
;; gsn_panel(wks,(/plot(0:)/),(/1,1/),resP)     ; now draw as one plot
;;gsn_panel(wks,plotMAP(0:1),(/1,2/),resP)     ; now draw as one plot
;;===============================================================
;;=====================================================
;;===============================================================
;;===============================================================
;;===============================================================
;;***********************************************
    res_P             = True                        ; panel mods desideeppink4
;   res_P@gsnDraw              = True       ; don't draw yet
;   res_P@gsnFrame             = False        ; don't advance frame yet
;  ; res_P@gsnPanelTop                      = 0.95
;   ;res_P@gsnPanelBottom                   = 0.60
;   res_P@amJust   			= "TopRight"
;   ;res_P@gsnPanelTop        	        = 0.52    ; only panel on lower half of page 
;   ;res_P@gsnPanelBottom     	        = 0.03    ; space for label bar
   ;res_P@gsn_PanelRight                   = 0.99
   ;res_P@gsnPanelLeft                     = 0.01
   res_P@gsnPanelFigureStringsPerimOn = False 
   res_P@gsnPaperOrientation = "portrait"
   res_P@gsnPanelDebug       = True
   res_P@gsnPanelTop                      = 0.99
   res_P@gsnPanelBottom                   = 0.55
   res_P@gsnPanelLeft                     = 0.01
   res_P@gsnPanelRight                    = 0.32
   res_P@amJust   			= "TopLeft"
   res_P@gsnPanelFigureStringsPerimOn   = False
   res_P@gsnPanelFigureStringsFontHeightF = 0.020
   res_P@gsnPanelFigureStringsBackgroundFillColor = "deeppink4"
   res_P@gsnPanelFigureStrings = (/"","","","","","","","h","i","j","k","l"/) ; add strings to panel
   ;res_P@gsnPanelMainString    = "~F12~AMOC indices at 30~S~o~N~N and 30~S~o~N~S"
   res_P@gsnMaximize = True                        ; blow up plot
  ; gsn_panel(wks,plot(0:0),(/1,1/),res_P)               ; create panel plots
;==============================================================      ;
;==============================================================
;==============================================================
;======================================================
;======================================================

;end do

frame(wks) 
end
