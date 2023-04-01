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
  yrLast = 2009
  d   = 1    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 7   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  nmos = 12
  nplots = 3
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
 itcz_d = new(dimsizes(itcz),typeof(itcz),itcz@_FillValue)
 x   = itcz&TIME
 x1  = x(::12)
 do tt = 0,nmos-1
    work = itcz(tt::nmos)
    work = dtrend_msg_n(x1,work,True,False,0)
    itcz_d(tt::nmos) = work
 end do

 printVarSummary(itcz_d)
 itcz = itcz_d
;===
delete([/x,x1,itcz_d,work,tt/])  
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
;===============================================================
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
printVarSummary(n_psmsl)                              ; variable overview
;=============================================testing/USA
;===================================================
n_psmsl = -rmAnnCycle1D(n_psmsl)
;n_psmsl = -n_psmsl 
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
;================================
;************************************************
yyyymm = yyyymm_time(yrStrt, yrLast, "integer")
;print(yyyymm) 
;
N_PSMSL = month_to_annual_weighted(yyyymm,n_psmsl,1) 
;N_PSMSL = month_to_season(n_psmsl,"DJF")
N_PSMSL!0 = "year"
N_PSMSL = dtrend(N_PSMSL,False)
N_PSMSL = dim_standardize_n(N_PSMSL,0,0)      ; normalize 
printVarSummary(N_PSMSL)
;
ITCZ = month_to_annual_weighted(yyyymm,itcz,1) 
;ITCZ = month_to_season(itcz,"DJF")
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
;
;===============================================================
;************************************************
; calculate cospectrum
;************************************************
print("============================================")
 
;--------------------------------------------------------------
shp = specxy_anal(N_PSMSL,PC,d,sm,pct)
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
shi = specxy_anal(N_PSMSL,ITCZ,d,sm,pct)
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
list_vars()
;print(xxxxxxxxxxxxxxxxxxx)
;===============================================================
; ==============================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
  latS   = -50.
  latN   =  70. 
  lonL   = -90.
  lonR   =  25.
  lag   =   1
  tau_const = 3600.*12.
  p95   =   0.05
  nmos = 12
  mxlag = 5
 
   yrStrt = 1900
  ;yrStrt = 1920
  yrLast = 2009	;
  yyyymm = yyyymm_time(yrStrt, yrLast, "integer")

  nplots = 3  
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
;------------------------------------------------------------------
;hadi_1 = hadi
;delete(hadi)
;hadi = month_to_annual_weighted(yyyymm,hadi_1,1) 
;printVarSummary(hadi)
;;hadi = dim_standardize_n(hadi,0,0)      ; normalize
;delete(hadi_1)

list_vars()
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
;
;
m_slp = dim_avg_n_Wrap(slp,0)
printVarSummary(m_slp)
printMinMax(m_slp,True)
;===================================
 ;  slp   =  -slp
   slp    = lonFlip(slp)
   printVarSummary(slp)        ; note the longitude coord
;===================================
; ==============================================================  
  printMinMax(slp,True)
 ; slp  = slp*conform(slp,atl,(/1,2/))
  printVarSummary(slp)
  printMinMax(slp,True)
;=====================================
   slpClm      = clmMonTLL(slp)
   printVarSummary(slpClm)                   ; (12,nlat,nlon)   
   slp         = calcMonAnomTLL(slp, slpClm)    ; replace with anonamlies
   slp@long_name = "ANOMALIES: "+slp@long_name
;===================================
 slp_d = new(dimsizes(slp),typeof(slp),slp@_FillValue)
 x   = slp&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = slp(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    slp_d(tt::nmos,:,:) = work
 end do
 printVarSummary(slp_d)
 slp = slp_d
;===
delete([/tt,x,x1,slp_d,work/])
;=========================================================================================
;slp_1 = slp
;delete(slp)
;slp = month_to_annual_weighted(yyyymm,slp_1,1)
;slp!0 = "time"
;printVarSummary(slp)
;delete(slp_1)
print("==============================================================")
; ==============================================================
;=========================================================================================
f= addfile("/home/hnnamchi/work/datasets/e20c/RE_qnet_e20c_1900_2010.nc","r")
 TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
; qnet    = f->lh(iYYYY,:,:)
  qnet    = f->qnet(iYYYY,:,:)
  printVarSummary(qnet)       ; variable overview
delete([/f,TIME,YYYY,iYYYY/])
qnet!0 = "time"
qnet!1 = "lat"
qnet!2 = "lon"
printVarSummary(qnet)
;===================================
 ;  qnet   =  -qnet
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
 qnet_d = new(dimsizes(qnet),typeof(qnet),qnet@_FillValue)
 x   = qnet&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = qnet(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    qnet_d(tt::nmos,:,:) = work
 end do
 printVarSummary(qnet_d)
 qnet = qnet_d
;===
delete([/tt,x,x1,qnet_d,work/])
;=========================================================================================
;qnet_1 = qnet
;delete(qnet)
;qnet = month_to_annual_weighted(yyyymm,qnet_1,1)
;qnet!0 = "time"
;printVarSummary(qnet)
;delete(qnet_1)
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
 ssh_d = new(dimsizes(ssh),typeof(ssh),ssh@_FillValue)
 x   = ssh&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = ssh(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    ssh_d(tt::nmos,:,:) = work
 end do

 printVarSummary(ssh_d)
 ssh = ssh_d

;===
delete([/x,x1,ssh_d,work,tt/])
;=========================================================================================
;ssh_1 = ssh
;delete(ssh)
;ssh = month_to_annual_weighted(yyyymm,ssh_1,1)
;ssh!0 = "time"
;printVarSummary(ssh)
;delete(ssh_1)
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
printMinMax(pc,True)
pc = dim_standardize_n(pc, 1, 0)  ; new variable
printMinMax(pc,True)
;print(xxxxxxxxxx)
;
;PC = month_to_annual_weighted(yyyymm,pc,1) 
;delete(pc)
;PC!0 = "time"
;hadi!0 = "time"
;=====================================================================
d_hadi = regCoef_n(pc,hadi,0,0)
x_var  = 100.*((escorc_n(pc,hadi,0,0))^2)     ;
;-----------------
df = hadi(0,:,:)
lag_ts_PC = esacr(pc(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
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
;
tval = onedtond(d_hadi@tval,dimsizes(d_hadi))
b = tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_hadid = d_hadi                ;transfers metadata
d_hadid = where(abs(prob).lt.p95,d_hadi,d_hadi@_FillValue) ; Significance mask 
printVarSummary(d_hadid)

copy_VarCoords(df,d_hadi)
copy_VarCoords(df,d_hadid)
;
delete(hadi)
delete([/tval,df,prob,b/])
;=====================================================================
;=====================================================================
;=====================================================================
d_qnet = regCoef_n(pc,qnet,0,0)      
x_qnet  = 100.*((escorc_n(pc,qnet,0,0))^2)     ;
;-----------------
df = qnet(0,:,:)
lag_ts_PC = esacr(pc(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
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
;
tval = onedtond(d_qnet@tval,dimsizes(d_qnet))
b = tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_qnetd = d_qnet		;transfers metadata
d_qnetd = where(abs(prob).lt.p95,d_qnet,d_qnet@_FillValue) ; Significance mask 
printVarSummary(d_qnetd)

copy_VarCoords(df,d_qnet)
copy_VarCoords(df,d_qnetd)
copy_VarCoords(df,x_qnet)
;
delete(qnet)
delete([/tval,df,prob,b/])
;=====================================================================
;=====================================================================
d_ssh = regCoef_n(pc,ssh,0,0)      
x_ssh  = 100.*((escorc_n(pc,ssh,0,0))^2)     ;
;-----------------
df = ssh(0,:,:)
lag_ts_PC = esacr(pc(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
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
;
tval = onedtond(d_ssh@tval,dimsizes(d_ssh))
b = tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_sshd = d_ssh		;transfers metadata
d_sshd = where(abs(prob).lt.p95,d_ssh,d_ssh@_FillValue) ; Significance mask 
printVarSummary(d_sshd)

copy_VarCoords(df,d_ssh)
copy_VarCoords(df,d_sshd)
copy_VarCoords(df,x_ssh)
;
delete(ssh)
delete([/tval,df,prob,b/])
;=====================================================================
delete(atl)
delete(clat)
;delete(hadi)
delete(qnetClm)
delete(hadiClm)
delete(slpClm)
delete(sshClm)
;=====================================================================
list_vars()
d_slp = regCoef_n(pc,slp,0,0)      
x_slp  = 100.*((escorc_n(pc,slp,0,0))^2)     ;
;-----------------
df = slp(0,:,:)
lag_ts_PC = esacr(pc(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
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
;
tval = onedtond(d_slp@tval,dimsizes(d_slp))
b = tval    ; b must be same size as tval (and df)
b = 0.5
prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
prob@long_name   = "regression coefficient"
prob@long_name = "probability"
printMinMax(prob,True)
printVarSummary(prob)
;
d_slpd = d_slp		;transfers metadata
d_slpd = where(abs(prob).lt.p95,d_slp,d_slp@_FillValue) ; Significance mask 
printVarSummary(d_slpd)

copy_VarCoords(df,d_slp)
copy_VarCoords(df,d_slpd)
copy_VarCoords(df,x_slp)
;
delete(slp)
delete([/tval,df,prob,b/])
;============================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("eps","fig0_3")         ; send graphics to PNG file

  plot = new(nplots,graphic)                ; create graphic array ; only needed if paneling

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
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015

  res@mpGeophysicalLineColor  = "transparent"
  res@mpFillDrawOrder      = "preDraw"

  
  
  res@gsnAddCyclic         = True
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = True     ; turn off individual lb's
  ;res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF   = 0.020 
  res@lbLabelFontHeightF   = 0.020
  res@lbLabelFont          = "Helvetica-Bold"     ; label font
  res@lbTitleString        = ""
 ; res@lbTitleAngleF       = 90.
  ;res@lbTitlePosition     = "Right"
  ;res@lbTitleDirection    = "Across"
  res@lbTopMarginF         = 0.30
  ;res@lbBottomMarginF      = 0.35
  ;res@vpWidthF            = 0.40
  ;res@vpHeightF           = 0.80
  res@lbLabelStride        = 10
  res@lbTitleString         = "Ocean temperature (K)"
  res@gsnCenterString  	      = ""   
  res@gsnLeftString           = ""  
  res@gsnRightString          = ""
  res@tiYAxisFontThicknessF  = 1.
  res@tmXBMajorThicknessF    = 1.
  res@tmXBMinorThicknessF    = 1.
  res@tmYLMajorThicknessF    = 1.
  res@tmYLMinorThicknessF    = 1.
  res@tmYLLabelFontThicknessF= 1.
  res@tmBorderThicknessF     = 1.
; 
  res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward

  res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
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
  res2@cnLineDashPattern  = 16           ;dashed low-pressure cells    
  res2@gsnAddCyclic         = True
  res2@cnLinesOn            = True       ; True is default
  res2@cnLineColor          = "grey95"       ; True is default
  res2@cnLineLabelsOn       = False        ; True is default
  res2@cnInfoLabelOn        = False
  res2@cnLineLabelPerimColor = "transparent"
  res2@cnLineLabelFontHeightF = 0.010
  res2@cnLineThicknessF       = 2.5
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
;-----------------------------------------------------------------------------

;*******************************************
; first plot
;*******************************************
    res@lbTitleString         = "Sea level pressure (hPa)"
    res@gsnStringFont		= "Helvetica-Bold"     ; label font
    res@gsnCenterStringFontHeightF   = 0.015
;*******************************************
;*******************************************
  delete(res@cnFillPalette)
  res@cnFillPalette        = "GMT_panoply"                        ;different colour bar  
  res@cnMinLevelValF       = -0.5            ; set min contour level
  res@cnMaxLevelValF       =  0.5             ; set max contour level
  res@cnLevelSpacingF      =  .05
  res@mpGeophysicalLineColor  = "foreground"
;
  ;plotA2 = gsn_csm_contour_map_ce(wks,d_slp,res)  ;           ; 
  ;plotB2 = gsn_csm_contour(wks,d_slpd,res1)       
  ;plotB2 = gsn_contour_shade(plotB2,d_slpd@_FillValue,p95,opt)
  ;plotC2 = gsn_csm_contour(wks,m_slp,res2)  ;
  ;overlay(plotA2,plotB2)
  ;overlay(plotA2,plotC2)
  ;plot(2) = plotA2

;=======================================================================
;=======================================================================
 delete(res@cnFillPalette)
 delete(res2@cnLevels)
 delete(res2@cnLineDashPatterns)
  ;res2@cnLevels    = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100./)   ; set levels 
  res2@cnLevels    = (/3.,6.,9.,12.,15.,18.,21.,24.,27.,30./)   ; set levels 
  res2@cnLineDashPattern  = 0           ;dashed low-pressure cells   
  res2@cnLineDashPatterns     =    (/0,0,0,0,0,0,0,0,0,0/) 
  res2@cnLineThicknessF       =    1.0

  ;res@cnFillPalette        = "temp_19lev"       ;different colour bar 
  res@cnFillPalette        = "amwg_blueyellowred"       ;different colour bar  
  res@cnMinLevelValF       = -10.            ; set min contour level
  res@cnMaxLevelValF       =  10.             ; set max contour level
  res@cnLevelSpacingF      =  1.
  res@mpGeophysicalLineColor  = "transparent"
  res@mpFillDrawOrder      = "postDraw"
  res@lbTitleString        = "Surface net heat flux (W m ~S~-2~N~)"
  ;
  ;plotA1 = gsn_csm_contour_map_ce(wks,d_qnet,res)  ;           ;
  ;plotB1 = gsn_csm_contour(wks,d_qnetd,res1)
  ;plotB1 = gsn_contour_shade(plotB1,d_qnetd@_FillValue,p95,opt)
  ;plotC1 = gsn_csm_contour(wks,x_qnet,res2)  ;
  ;overlay(plotA1,plotB1)
  ;overlay(plotA1,plotC1)
  ;plot(1) = plotA1

;=======================================================================
;=======================================================================
  delete(res@cnFillPalette)
  ;res@cnFillPalette        = "ncl_default"       ;different colour bar 
  res@cnFillPalette        = "amwg_blueyellowred"       ;different colour bar  
  res@cnMinLevelValF       = -1.            ; set min contour level
  res@cnMaxLevelValF       =  1.             ; set max contour level
  res@cnLevelSpacingF      =  .1
  res@mpGeophysicalLineColor  = "transparent"
  res@mpFillDrawOrder      = "preDraw"
;  res@lbTitleString                    = "Sea surface height (~F10~e~S~~F21~-2~N~m)"
  res@lbTitleString                    = "Sea surface height (10~S~-2~N~m)"
  ; 
  plotA0 = gsn_csm_contour_map_ce(wks,d_ssh,res)  ;           ;
  plotB0 = gsn_csm_contour(wks,d_sshd,res1)
  plotB0 = gsn_contour_shade(plotB0,d_sshd@_FillValue,p95,opt)
  plotC0 = gsn_csm_contour(wks,x_ssh,res2)  ;
  overlay(plotA0,plotB0)
  overlay(plotA0,plotC0)
  plot(0) = plotA0

   print("==========================ssh expl var++++++++++++++++")
   printMinMax(x_ssh,True)
   print("==========================qnet expl var++++++++++++++++")
   printMinMax(x_qnet,True)
   print("======================================++++++++++++++++")
 ;==========================================================================================
;==========================================================================================
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
   ;r1@tmYRValues = ispan(-180,180,45)
   ;r1@tmYRLabels  =  (/"-~F33~p~","-3~F33~p~N~/4","-~F33~p~N~/2","-~F33~p~N~/4","0","~F33~p~N~/4","~F33~p~N~/2","3~F33~p~N~/4","~F33~p~N~"/)

   r1@tmYRMinorValues	= ispan(-180,180,30)
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
   r1@tiMainString  = ""                       ; title
   r1@tiXAxisString = ""  ; xaxis
   r1@tmYRFormat             = "@*+^f"		;Don't force 0's 					; and the left Y-axis.  
 
;------------------------------------------------------------------------
;************************************************
;------------------------------------------------------------------------  
 
   ;r0@tiMainString  = "~F21~Sea-level index ~F34~.~F21~ PC2"                       ; title
   plot(1) =gsn_csm_xy2(wks,1./shp@frq,shp@coher,shp@phase,r0,r1);

   ;r0@tiMainString  = "~F21~Sea-level index ~F34~.~F21~ ITCZ"                          ; title
   plot(2) =gsn_csm_xy2(wks,1./shi@frq,shi@coher,shi@phase,r0,r1);

   ;r0@tiMainString  = "~F21~ITCZ ~F34~.~F21~ PC2"                           ; title
;   plot(2) =gsn_csm_xy2(wks,1./sip@frq,sip@coher,sip@phase,r0,r1);
;-----------------------------------------------------------------------------
;=========================================================================
;=========================================================================
;=========================================================================
;************************************************
; add station labels
   y = (/36.822,36.947,37.247,38.317,38.873,38.983,39.267, 38.573, 37.165, 38.782, 39.933, 38.968, 39.355, 40.467, 40.700, 40.793, 41.173,\
        41.048, 41.360, 41.505, 41.807, 41.523, 41.285, 42.353/)

   x = (/-76.293,-76.330,-76.500,-76.452,-77.022,-76.480,-76.578,-76.068,-75.988,-75.120,-75.142,-74.960,-74.418,-74.008,-74.013,-73.782,-73.182,\
        -71.960,-72.090,-71.327,-71.400,-70.672,-70.097,-71.053/)
print(x)
print(y)
;
  npts = dimsizes(y)
  print(npts)
;--------------------
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 16                          ; index for circle
  polyres@gsMarkerColor     = "black"
  polyres@gsMarkerSizeF     = 0.005                        ; size
  polyres@gsMarkerOpacityF  = 0.15

  ;txres@txFontHeightF = 0.005   ;reduce fontsize

  centers_n   = new(npts,graphic)
  print(centers_n)

  do i = 0,npts-1,1
     ;txres@txFontColor = icol(i,:)
     ;polyres@gsMarkerColor = icol(i,:)
   ;  centers_n(i)   = gsn_add_text(wks,plot(0),snums(i),x(i)+0.10,y(i),txres)
     centers_n(i)   = gsn_add_polymarker(wks,plot(0),x(i),y(i),polyres)
  end do
delete(i)

;=====================================================================
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
     do i = 0,nplots-1,1
      tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres)
     end do
 delete(i)
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
     do i = 1,nplots-1,1
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
     do i = 1,nplots-1,1
      dum1(i) = gsn_add_polyline(wks,plot(i), (/x0,x1/),(/y0,y0/) ,pres1)
     end do
  delete(i)
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
     do i = 1,nplots-1,1
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
     do i = 1,nplots-1,1
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
  resP@gsnPanelBottom                   = 0.55
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
gsn_panel(wks,plot(0:),(/1,3/),resP)     ; now draw as one plot
;;===============================================================
;;=============================================================================
;=======================================================================
;frame(wks) 
end
