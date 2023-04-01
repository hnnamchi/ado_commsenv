
begin
; ==============================================================
; ==============================================================
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
  latS   = -45.
  latN   =  60. 
  lonL   = -100.
  lonR   =  25.
  
  

npts = 24	; # of tide gauges
; ==============================================================
; Open & read mask data
; ==============================================================

snames = (/"01. Portsmouth","02. Sewells Point","03. Gloucester Point","04. Solomon's Island","05. Washington DC","06. Annapolis",\
	"07. Baltimore","08. Cambridge II","09. Kiptopeke Beach","10. Lewes","11. Philadelphia","12. Cape May","13. Atlantic City",\
	"14. Sandy Hook","15. New York","16. Willets Point","17. Bridgeport","18. Montauk","19. New London","20. Newport",\
	"21. Providence","22. Woods Hole","23. Nantucket Island","24. Boston"/)

snums = (/"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"/)

y = (/36.822,36.947,37.247,38.317,38.873,38.983,39.267, 38.573, 37.165, 38.782, 39.933, 38.968, 39.355, 40.467, 40.700, 40.793, 41.173,\
	41.048, 41.360, 41.505, 41.807, 41.523, 41.285, 42.353/)

x = (/-76.293,-76.330,-76.500,-76.452,-77.022,-76.480,-76.578,-76.068,-75.988,-75.120,-75.142,-74.960,-74.418,-74.008,-74.013,-73.782,-73.182,\
	-71.960,-72.090,-71.327,-71.400,-70.672,-70.097,-71.053/)
print(snames)
print(snums)
print(x)
print(y)
;print(xxxxxxxxxxxx)
;===============================================================================================================
;===============================================================================================================

yrStrt = 1900
yrLast = 2010
;===============================================================================================================
;===============================================================================================================
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/n_psmsl_ordered_eastern_USA_1856_2020.nc", "r")
;f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_all_eastern_USA_1856_2017.nc", "r")
f      = addfile ("/home/hnnamchi/work/datasets/amoc/psmsl/psmsl_ordered_eastern_USA_1856_2017.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
 
psmsl   = f->psmsl(iYYYY,:)
printVarSummary(psmsl)   
;
psmsl = where(psmsl.gt.0.,psmsl,psmsl@_FillValue)
;------------------------------------------------------------------
do ii = 0,dimsizes(psmsl(0,:))-1,1
	print(psmsl(:,ii))
end do

;delete(psmsl)

;psmsl@_FillValue = psmsl@missing_value
;print(psmsl@_FillValue)
;print(xxxxxxxxxxxxxxx)
;-----------------------------------------------------------------
psmsl = psmsl - psmsl


;printMinMax(PSMSL,True)

;
do ii = 0,dimsizes(psmsl(0,:))-1,1
psmsl(:,ii) = psmsl(:,ii)+ii
end do 
printMinMax(psmsl,True)
;
delete(ii)
delete(f)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;===============================================================================================================
yyyymm = cd_calendar(psmsl&time,-1)  
yrfrac = yyyymm_to_yyyyfrac(yyyymm, 0.0)

;===============================================================================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("eps","fig0_S1")         ; send graphics to PNG file

  plot = new(2,graphic)                ; create graphic array ; only needed if paneling
;====================================================
levels = ispan(1,200,1)
  icol = span_color_rgba("cmp_flux",24)
 ; print(icol)
;====================================================
;====================================================
   l               = True                        ; plot mods desired
   l@gsnDraw       = False                       ; do not draw
   l@gsnFrame      = False                       ; do not advance frame
   l@tiMainString  = ""                       ; title
   l@tiXAxisString = "Time (year)"  ; xaxis
   l@tiYAxisString = ""                  ; yaxis
   l@tiYAxisFontHeightF = 0.008
   l@tiXAxisFontHeightF = 0.008
   l@xyLineThicknessF   = 5.
   ;l@xyLineThicknesses   = 5.       ; Define line thicknesses 
   l@xyDashPatterns      = 0           ; Dash patterns 
   l@xyLineColors        = icol
   l@vpHeightF		 = 0.29                   ; change aspect ratio of plot
   l@vpWidthF		 = 0.35           
   l@tmXBLabelFontHeightF = 0.008
   l@tmYRLabelFontHeightF = 0.008
   l@tmYLFormat 	 = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   l@trYMinF 		 =  4.
   l@trYMaxF 		 =  31.
   l@trXMinF 		 =  1900.
   l@trXMaxF 		 =  2010.
   ;l@gsnXRefLine         = 1905.0 
   l@gsnXRefLineDashPattern    =  16           ; create a reference line
  ;l@gsnYRefLineThicknesses    = (/0.0/)  

   l@tmYRLabelsOn	 = False
   l@tmYLLabelsOn	 = False
   l@tmYROn		 = False
   l@tmYLOn		 = False

  ; l@tmYLMode 		= "Explicit"
  ; l@tmYLValues		= fspan(0,30,30)
  ; print(l@tmYLValues)
      
   l@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmXBMinorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmYLMinorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmXTMinorOutwardLengthF = 0.0               ; draw tickmarks inward
   l@tmYRMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  
;************************************************ 
;------------------------------------------------

 print(icol)

; plot(1) = gsn_csm_xy(wks,yrfrac,psmsl(stations|0:,time|:),l)

  plot(1) = gsn_csm_xy(wks,yrfrac,psmsl(stations|6:,time|:),l)

;-----------------------------------------------------------------------------
  res                             = True
  res@gsnDraw     		  = False                       ; do not draw
  res@gsnFrame  		  = False                       ; do not advance frame
  res@mpFillOn                    = False        ; turn off gray fill
  res@mpOutlineBoundarySets       = "National"   ; turn on country boundaries
  res@mpGeophysicalLineColor      = "Navy"       ; color of cont. outlines
  res@mpGeophysicalLineThicknessF = 0.5          ; thickness of outlines
  res@vpHeightF			 = 0.90                  ; change aspect ratio of plot
  res@vpWidthF			 = 0.20         
  res@mpMaxLatF                   = 43.           ; choose subregion           
  res@mpMinLatF                   = 36.
  res@mpMaxLonF                   = -69.5
  res@mpMinLonF                   = -77.5
  res@tmYLMode   = "Explicit"
  res@tmYLValues = (/36,38,40,42/)
  res@tmYLLabels = (/"36~S~o~N~N","38~S~o~N~N","40~S~o~N~N","42~S~o~N~N"/)
  res@tmYLLabelFontHeightF = 0.0045
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/-78,-76,-74,-72,-70/)
  res@tmXBLabels = (/"78~S~o~N~W","76~S~o~N~W","74~S~o~N~W","72~S~o~N~W","70~S~o~N~W"/)
  res@tmXBLabelFontHeightF = 0.0045
  res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXBMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYLMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMinorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@gsnLeftString  		= ""
  res@gsnRightString 		= ""
  res@gsnCenterString   	= ""                       ; title
  res@cnLineLabelsOn       = False        ; True is default
  res@cnInfoLabelOn        = False
  res@cnLineDashPattern       = 1
  res@cnLevelSelectionMode   = "ExplicitLevels" 
  res@cnLevels		     =  (/0./)  
  res@gsnContourNegLineDashPattern = 1
  res@gsnContourZeroLineThicknessF = -1

  plot(0) = gsn_csm_map(wks,res)                 ; draw map

 ;plot(0) = gsn_csm_contour_map(wks,m_adt,res)                 ; draw map
;=====================================================================
;************************************************
; add station labels
   yy = fspan(0.535,0.827,npts)   ; 
   print(yy)
   xx = 0.89

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.008
  txres@txJust        = "CenterLeft"

  do i = 0,npts-1,1
     gsn_text_ndc(wks,snames(i),xx,yy(i),txres)
  end do

delete(i)
;--------------------
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 16                          ; index for circle
  ;polyres@gsMarkerColor	    = "cyan"
  polyres@gsMarkerSizeF     = .006                        ; size
  polyres@gsMarkerOpacityF  = 1.0

  txres@txFontHeightF = 0.005	;reduce fontsize

  centers_n   = new(npts,graphic)
  print(centers_n)

  do i = 0,npts-1,1
     txres@txFontColor = icol(i,:)
     polyres@gsMarkerColor = icol(i,:)
     centers_n(i)   = gsn_add_text(wks,plot(0),snums(i),x(i)+0.10,y(i),txres)
     centers_n(i)   = gsn_add_polymarker(wks,plot(0),x(i),y(i),polyres)
  end do
delete(i) 

;=====================================================================
;=====================================================================
nplots = 2
;-----------------------------------------------------------------------------
  txres1              = True
  txres1@txAngleF      = 0.          ; Rotate string clockwise 
  txres1@txFontHeightF = 0.010   ; Use same font height as left axis
  txres1@txFont                = "Helvetica-Bold"
  amres                  = True
  amres@amParallelPosF   =  -0.42 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.42    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  
  txid0 = gsn_create_text(wks, "a", txres1);
  annoid0 = gsn_add_annotation(plot(0),txid0,amres)  ; Attach string to plot

  txres1@txFontHeightF = 0.015   ; Use same font height as left axis
  txid1 = gsn_create_text(wks, "b", txres1);
  annoid1 = gsn_add_annotation(plot(1),txid1,amres)  ; Attach string to plot
;=====================================================================
;=====================================================================
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelTop                      = 0.99
  res_P@gsnPanelBottom                   = 0.35
  res_P@gsnPanelLeft                     = 0.01
  res_P@gsnPanelRight                    = 0.45
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
  res_P@gsnMaximize                      = True         ; large format
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
  resP@gsnPanelBottom                   = 0.36
  resP@gsnPanelLeft                     = 0.44
  resP@gsnPanelRight                    = 0.89
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
