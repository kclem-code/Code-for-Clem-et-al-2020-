begin

  clev = (/"500","300","200","150","100","50","30"/)
  clev1 = (/"surface","500hPa","300hPa","200hPa","150hPa","100hPa","50hPa","30hPa"/)

  cseas = (/"Annual","MAM","JJA","SON","DJF"/)

  
  t_data = new((/8,5,63/),"float")
  z_data = new((/7,5,63/),"float")



;  t_data = -999.9
;  z_data = -999.9

  sigt=new((/40,33,33/),float)
  slopet=new((/40,33,33/),float)    
  sigz=new((/35,33,33/),float)
  slopez=new((/35,33,33/),float)    


  slopet@_FillValue=-999.9
  sigt@_FillValue=-999.9
  slopez@_FillValue=-999.9
  sigz@_FillValue=-999.9


  cnt = 0
  time1=ispan(1957,2018,1)
  printVarSummary(time1)

  in1 = asciiread("../radiosonde/plot_data/surft.txt",(/63,6/),"float")

  t_data(0,0,:) = in1(:,1)
  t_data(0,1,:) = in1(:,2)
  t_data(0,2,:) = in1(:,3)
  t_data(0,3,:) = in1(:,4)
  t_data(0,4,:) = in1(:,5)
  ann_mean_t = in1(:61,1)
  year = in1(:61,0)
  delete(in1)

;  printVarSummary(ann_mean_t)
;  print(year)
;  exit
  ann_mean_t@_FillValue = -999.9
  run_trend_ann = ann_mean_t
  run_trend_ann = -999.9
  run_trend_ann@_FillValue = -999.9
  cis = new((/2,dimsizes(run_trend_ann)/),"float")
  cis =-999.9
  cis@_FillValue = -999.9

  do ilev = 0 , 6

  in1 = asciiread("../radiosonde/plot_data/"+clev(ilev)+"t.txt",(/59,6/),"float")  
  in2 = asciiread("../radiosonde/plot_data/"+clev(ilev)+"z.txt",(/59,6/),"float")  
  
  do iseas = 0 , 4

  t_data(ilev+1,iseas,4:) = in1(:,iseas+1)
  z_data(ilev,iseas,4:) = in2(:,iseas+1)

  end do
  end do

  t_data@_FillValue = -999.9
  z_data@_FillValue = -999.9

;data are read in calculate trends!

  cnt=0

  do i = 0,32
  sval = i
  eval = i+29
  ttime = ispan(sval,eval,1)
  tmpt1 = ann_mean_t(sval:eval)
  if (num(.not.ismissing(tmpt1)).eq.30) then
  bt = regline(ttime,tmpt1)
  run_trend_ann(i) = bt*10.
  cis(1,i) = (bt@rstd*2.05 + bt)*10.
  cis(0,i) = (bt - bt@rstd*2.05)*10.
  end if
;  print((i+1957)+" "+run_trend_ann(i)+" "+uci(i)+" "+lci(i))
  end do
;  exit

  ann_mean_t = dim_standardize(ann_mean_t,0)

  do ilev = 0,6

  do iseas = 0,4
 
  do i=0,32
   do j=0,32
  sval=i
  eval=29+j
  diff = eval - sval + 1
  if (diff.ge.30) then

  tmp1 = t_data(ilev+1,iseas,sval:eval)
  tmp2 = z_data(ilev,iseas,sval:eval)
  time = time1(sval:eval)
  
;  print(i+" "+j+" "+sval+" "+eval)
;  print(time)

;  print(ilev+" "+iseas+" "+time1(sval)+" "+time1(eval)+" "+num(ismissing(tmp1)))
  if(num(.not.ismissing(tmp1)).ge.30) then
  
  b1=regline(time,tmp1)
  df1 = b1@nptxy-2
  prob1 = (betainc(df1/(df1+b1@tval^2),df1/2.0,0.5))

  slopet(cnt+5,i,j) = b1*10
  sigt(cnt+5,i,j) = 1-prob1
  end if

  if(num(.not.ismissing(tmp2)).ge.30) then
  
  b2=regline(time,tmp2)
  df2 = b2@nptxy-2
  prob2 = (betainc(df2/(df2+b2@tval^2),df2/2.0,0.5))

  slopez(cnt,i,j) = b2*10
  sigz(cnt,i,j) = 1-prob2
  end if

  delete([/tmp1,tmp2,time/])
  end if ;check if 30 years



  end do ;all times
  end do ;all times

;  exit

  cnt = cnt + 1
  end do ;seasons
  end do ;levels


;repeat for surface temperature
  cnt = 0
  do ilev = 0,0

  do iseas = 0,4
 
  do i=0,32
   do j=0,32
  sval=i
  eval=29+j
  diff = eval - sval + 1
  if (diff.ge.30) then

  tmp1 = t_data(ilev,iseas,sval:eval)

  if(num(.not.ismissing(tmp1)).ge.30) then
  
  b1=regline(time1(sval:eval),tmp1)
  df1 = b1@nptxy-2
  prob1 = (betainc(df1/(df1+b1@tval^2),df1/2.0,0.5))

  slopet(cnt,i,j) = b1*10
  sigt(cnt,i,j) = 1-prob1
  end if

  delete(tmp1)
  end if ;check if 30 years



  end do ;all times
  end do ;all times


  cnt = cnt + 1
  end do ;seasons
  end do ;levels


  y1_start = 1957
  y1_end   = 1989
  nyrs1    = y1_end-y1_start+1
  y2_start = 1986
  y2_end   = 2018
  nyrs2    = y2_end-y2_start+1

; ---------- Plot parameters ---------- ;

  wks = gsn_open_wks("ps","NCC_figure_revised_final")
  gsn_define_colormap(wks,"cmocean_balance")    ; choose colormap
  colors = gsn_retrieve_colormap(wks)  ; retrieve color map for editing.
  gsn_define_colormap(wks,colors)

  plot = new(9,graphic)



  res3                    = True              ; plot mods desired
  res3@gsnDraw            = False             ; don't draw yet
  res3@gsnFrame           = False             ; don't advance frame yet


  res3@vpHeightF 	  = 0.11               ; change aspect ratio of plot
  res3@vpWidthF 	  = 0.8
  res3@vpXF = 0.12
  res3@vpYF = 0.75
  res3@trXMinF	          = 1957              ; set x-axis minimum
  res3@trXMaxF            = 1989		     ; set x-axis max
  res3@trYMinF	          = -1               ; set y-axis minimum
  res3@trYMaxF            = 1.		     ; set y-axis max
  res3@tmXBMaxTicks	  = 10
  res3@tmXBMinorPerMajor  = 4
  res3@tmYMajorGrid	  = True		     
  res3@tmYMajorGridThicknessF      = 1.0	     ; properties for y-grid
  res3@tmYMajorGridLineDashPattern = 2.0
  res3@tmXMajorGridThicknessF      = 1.0	     ; properties for y-grid
  res3@tmXMajorGridLineDashPattern = 2.0
  res3@tmYLLabelFontHeightF 	   = .01
  res3@tmYRLabelFontHeightF 	   = .01
  res3@tmXBLabelFontHeightF 	   = .01
  res3@tmXTOn = False  
  res3@tmXMajorGrid      = True
  res3@tmXBMode			  = "Manual"
  res3@tmXBTickStartF		  = "1955"
  res3@tmXBTickSpacingF		  = "5"
  res3@tmXBTickEndF		  = "2020"
  res3@gsnYRefLine = 0.0
  res3@tmYLMode			  = "Manual"
  res3@tmYLTickStartF		  = "-3"
  res3@tmYLTickSpacingF		  = ".25"
  res3@tmYLTickEndF		  = "3"
  res3@tmYRMode			  = "Manual"
  res3@tmYRTickStartF		  = "-3"
  res3@tmYRTickSpacingF		  = ".25"
  res3@tmYRTickEndF		  = "3"
  res3@tmYRLabelsOn = True
;  res3@tiYAxisString	 = "30-yr Running Trend (~F34~0~F~C / decade)"		 ; add a axis title    
  res3@tiXAxisString	 = "Starting Year"     	     	 ; add a axis title    
  res3@tiXAxisFontHeightF = 0.012           		 ; change axis title font heights
  res3@tiYAxisFontHeightF = 0.012           		 ; change axis title font heights

  res3@gsnXYFillColors = "grey"
  res3@xyLineColor = "transparent"
plot_bot  = gsn_csm_xy(wks,year,cis,res3)       ; create plot 
;  res3@gsnxyFillColors = "transparent"
  delete(res3@gsnXYFillColors)

  res3@xyLineThicknessF  = 3
  res3@xyLineColor   = "black"
  res3@xyDashPattern	  = 0      ; make all lines solid
plot_mid  = gsn_csm_xy(wks,year,run_trend_ann,res3)
overlay(plot_bot,plot_mid)

  res3@tmYLMode			  = "Manual"
  res3@tmYLTickStartF		  = "-3"
  res3@tmYLTickSpacingF		  = "1.0"
  res3@tmYLTickEndF		  = "3"
  res3@tmYRMode			  = "Manual"
  res3@tmYRTickStartF		  = "-3"
  res3@tmYRTickSpacingF		  = "1.0"
  res3@tmYRTickEndF		  = "3"
  res3@tmYRLabelsOn = True
  res3@vpYF = 0.955
  res3@xyLineThicknessF  = 4
  res3@xyLineColor   = "red"
  res3@trXMinF	          = 1957              ; set x-axis minimum
  res3@trXMaxF            = 2018		     ; set x-axis max
  res3@trYMinF	          = -3               ; set y-axis minimum
  res3@trYMaxF            = 3.		     ; set y-axis max
;  res3@tiYAxisString	 = "Standardized Temperature (~F34~0~F~C)"		 ; add a axis title    
  res3@tiXAxisString	 = "Year"     	     	 ; add a axis title    


plot_top = gsn_csm_xy(wks,year,ann_mean_t,res3)       ; create plot 

;overlay(plot_bot,plot_top)
  poly2 = new(2,graphic)

  pxres1 = True
  pxres1@txFontHeightF = 0.01
  pxres1@txAngleF      = 0.
  pxres1@txJust  = "CenterLeft"
  pxres1@txPerimOn = True
  pxres1@txPerimColor = "white"
  pxres1@txBackgroundFillColor = "white"
  pres1 = True
  pres1@gsLineColor      = "black"                     ; color of lines
  pres1@gsLineThicknessF = 2.0   
;  dum = gsn_add_text(wks,plot_top,"30-yr Running Trend + 95% CI",1966,2.6,pxres1)
;  dum = gsn_add_text(wks,plot_top,"Surface Temperature",1966,2.1,pxres1)
;   poly2(0) = gsn_add_polyline(wks,plot_top,(/1960,1965/),(/2.6,2.6/),pres1)  
;  pres1@gsLineColor      = "red"                     ; color of lines
;   poly2(1) = gsn_add_polyline(wks,plot_top,(/1960,1965/),(/2.1,2.1/),pres1)  

  draw(wks)

  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@cnFillOn            = True
  res@cnLinesOn           = False
  res@cnLineLabelsOn      = False
  res@cnInfoLabelOn       = False
 ; res@cnFillColors	  = (/2,7,11,16,21,26,30,35,40,45,50,54,59,64,68,73,78,81,85,87,92,0,146,151,156,161,166,171,177,181,186,191,196,201,206,211,216,221,226,231,236,241/)
  res@gsnSpreadColors     = True
  res@gsnSpreadColorStart   = 2
  res@gsnSpreadColorEnd     = 255
  res@cnMissingValFillColor = "white"
  res@cnLevelSelectionMode = "ManualLevels"

  res@cnFillMode          = "RasterFill"
  res@tmXBPrecision       = 4
  res@tmYLPrecision       = 4
  res@gsnStringFontHeightF = 0.035
  res@tmXTOn              = False
  res@tmYROn              = False
  res@tiXAxisString       = "End year"
  res@tiYAxisString       = "Start year"
  res@tiYAxisFontHeightF  = 0.04;0.035
  res@tiXAxisFontHeightF  = 0.04;0.035
  res@tiXAxisOffsetYF  = 0.007
  res@lbLabelBarOn        = False
  res@tmXBLabelFontHeightF = 0.035;0.025
  res@tmYLLabelFontHeightF = 0.035;0.025
  res@tmXBLabelAngleF     = 45

  res@tmXBMode            = "Explicit"
  res@tmXBValues          = ispan(y2_start,y2_end,4)
  res@tmXBLabels          = ispan(y2_start,y2_end,4)
  res@tmXBMinorValues     = ispan(y2_start+5,y2_end-5,2)

  res@tmYLMode            = "Explicit"
  res@tmYLValues          = ispan(y1_start,y1_end,4)
  res@tmYLLabels          = ispan(y1_start,y1_end,4)
  res@tmYLMinorValues     = ispan(y1_start+5,y1_end-5,2)

  res@sfXArray = ispan(y2_start,y2_end+1,1)-0.5
  res@sfYArray = ispan(y1_start,y1_end+1,1)-0.5


  NSEAS = 5
  NLEV1  = 8
  NLEV2  = 7
; ------------------------------------------------------

  season_list   = (/"Annual","DJF","MAM","JJA","SON"/)

  plot11 = new(8,graphic)
;  plot21 = new(35,graphic)



  res2 = True
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frame
  res2@gsnRightString       = ""             ; turn off special titles
  res2@gsnLeftString        = ""
  res2@tiMainString         = ""
  res2@gsnCenterString      = ""

  res2@cnLinesOn           = False
  res2@cnLineLabelsOn      = False

  res2@cnFillOn 	   = True       ; turn on color/pattern fill
  res2@cnMonoFillPattern   = False      ; allow different patterns
  res2@cnMonoFillColor 	   = True       ; only use one color (black)
  res2@lbLabelBarOn        = False
  res2@cnLevelSelectionMode = "ExplicitLevels" ; set explicit cnlev
  res2@cnLevels   	    = (/.90,0.95/)    ; only have 1 contour level
  res2@cnFillPatterns	    = (/-1,3,17/) ; stipple <0.05, don't fill >=0.05
  res2@cnInfoLabelOn       = False       ; turn off info label

  res2@sfXArray = ispan(y2_start,y2_end+1,1)-0.5
  res2@sfYArray = ispan(y1_start,y1_end+1,1)-0.5



  plot12 = new(40,graphic)
  plot22 = new(35,graphic)

  pres = True
  pres@sfXArray = ispan(y2_start,y2_end+1,1)-0.5
  pres@sfYArray = ispan(y1_start,y1_end+1,1)-0.5
  pres@gsLineColor      = "black"                     ; color of lines
  pres@gsLineThicknessF = 2.0   
 
  poly3 = new(3,graphic)
  poly41 = new(40,graphic)
  poly42 = new(40,graphic)

  xvals = (/0.12,0.32,0.52,0.72,0.12,0.32,0.52,0.72/)
  yvals = (/0.65,0.65,0.65,0.65,0.2,0.2,0.2,0.2/)
  iplot1 = 0

;  res3@vpHeightF 	  = 0.001               ; change aspect ratio of plot
;  res3@vpWidthF 	  = 0.001

    do pval = 0,3

  res@cnMinLevelValF      = -1.
  res@cnMaxLevelValF      =  1.
  res@cnLevelSpacingF     = 0.05
;  res@vpXF = xvals(pval)
;  res@vpYF = yvals(pval)
  
      plot(pval) = gsn_csm_contour(wks,slopet(pval+1,:,:),res)
       plot11(pval) = gsn_csm_contour(wks,sigt(pval+1,:,:),res2)
       overlay(plot(pval),plot11(pval))
       poly41(pval) = gsn_add_polyline(wks,plot(pval),(/y2_start,y2_end/),(/y1_start,y1_end/),pres)
    end do

  res@cnMinLevelValF      = -40
  res@cnMaxLevelValF      =  40
  res@cnLevelSpacingF     = 2

      plot(4) = gsn_csm_contour(wks,slopez(0,:,:),res)
       plot11(4) = gsn_csm_contour(wks,sigz(0,:,:),res2)
       overlay(plot(4),plot11(4))
       poly41(4) = gsn_add_polyline(wks,plot(4),(/y2_start,y2_end/),(/y1_start,y1_end/),pres)

  res@cnMinLevelValF      = -1.
  res@cnMaxLevelValF      =  1.
  res@cnLevelSpacingF     = 0.05

      plot(5) = gsn_csm_contour(wks,slopet(5,:,:),res)
       plot11(5) = gsn_csm_contour(wks,sigt(5,:,:),res2)
       overlay(plot(5),plot11(5))
       poly41(5) = gsn_add_polyline(wks,plot(5),(/y2_start,y2_end/),(/y1_start,y1_end/),pres)

  res@cnMinLevelValF      = -40
  res@cnMaxLevelValF      =  40
  res@cnLevelSpacingF     = 2

      plot(6) = gsn_csm_contour(wks,slopez(15,:,:),res)
       plot11(6) = gsn_csm_contour(wks,sigz(15,:,:),res2)
       overlay(plot(6),plot11(6))
       poly41(6) = gsn_add_polyline(wks,plot(6),(/y2_start,y2_end/),(/y1_start,y1_end/),pres)

  res@cnMinLevelValF      = -1
  res@cnMaxLevelValF      =  1
  res@cnLevelSpacingF     = 0.05

      plot(7) = gsn_csm_contour(wks,slopet(20,:,:),res)
       plot11(7) = gsn_csm_contour(wks,sigt(20,:,:),res2)
       overlay(plot(7),plot11(7))
       poly41(7) = gsn_add_polyline(wks,plot(7),(/y2_start,y2_end/),(/y1_start,y1_end/),pres)


   poly3(0) = gsn_add_polyline(wks,plot(0),(/1990,1998/),(/1980,1972/),pres)
   poly3(1) = gsn_add_polyline(wks,plot(0),(/1998,1995/),(/1972,1972/),pres)
   poly3(2) = gsn_add_polyline(wks,plot(0),(/1998,1998/),(/1975,1972/),pres)


  pxres = True
  pxres@txFontHeightF = 0.035
  pxres@txAngleF      = 0.
  pxres@txJust  = "CenterLeft"


  pxres@txPerimOn = False
  pxres@txJust        = "CenterLeft"
  dum = gsn_add_text(wks,plot(0),"Shorter Periods",1988,1982,pxres)
;  pxres@txJust        = "CenterRight"
  dum1 = gsn_add_text(wks,plot(0),"Longer ~C~Periods",1986,1969,pxres)

  pxres@txAngleF      = 45.
  pxres@txJust        = "CenterCenter"
  dum2 = gsn_add_text(wks,plot(1),"30-year Trends",1998,1972,pxres)



  resPanel                     = True
  resPanel@gsnMaximize         = True
  resPanel@gsnPanelLabelBar    = True
  resPanel@gsnPaperOrientation = "portrait"
  resPanel@lbOrientation       = "horizontal"
  resPanel@lbLabelAutoStride   = True
  resPanel@lbLabelStride       = 5
  resPanel@pmLabelBarHeightF   = 0.07
  resPanel@pmLabelBarWidthF    = 0.43
  resPanel@lbLabelFontHeightF  = 0.012
  resPanel@lbBoxLinesOn        = False
  resPanel@gsnPanelBottom      = 0.05
  resPanel@gsnPanelTop      = 0.55
  resPanel@gsnFrame = False  
;  resPanel@gsnPanelLeft        = 0.01
  resPanel@lbTitleFontHeightF  = .014                 ; make title smaller
  resPanel@pmLabelBarOrthogonalPosF = -.0           ; move whole thing down
  resPanel@gsnPanelXWhiteSpacePercent = 3
  resPanel@gsnPanelYWhiteSpacePercent = 5
  resPanel@lbTitleString       = "Trend (decade~S~-1~N~)"
;  resPanel@gsnPanelFigureStrings    = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)
  resPanel@amJust                   = "TopLeft"
  resPanel@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  resPanel@gsnPanelFigureStringsFontHeightF = 0.012
  resPanel@gsnPanelFigureStringsPerimOn = False

  gsn_panel(wks,plot(:7),(/2,4/),resPanel)

;=========================================================

  txres               = True
  txres@txJust        = "CenterCenter"
  txres@txFontColor   = "black"
  txres@txFontHeightF = 0.015
;  txres@txAngleF      = 90.
;  txres@tfPolyDrawOrder = "postdraw"
;  drawNDCGrid(wks)

   gsn_text_ndc(wks,"MAM",0.13,0.525,txres)
   gsn_text_ndc(wks,"JJA",0.37,0.525,txres)
   gsn_text_ndc(wks,"SON",0.605,0.525,txres)
   gsn_text_ndc(wks,"DJF",0.84,0.525,txres)

   gsn_text_ndc(wks,"500 hPa ~C~Height",0.15,0.29,txres)
   gsn_text_ndc(wks,"500 hPa ~C~Temp",0.39,0.29,txres)
   gsn_text_ndc(wks,"150 hPa ~C~Height",0.61,0.29,txres)
   gsn_text_ndc(wks,"150 hPa ~C~Temp",0.85,0.29,txres)

   txres@txJust = "CenterLeft"
   gsn_text_ndc(wks,"a) Amundsen-Scott Annual Mean Surface Temperature",0.1,0.97,txres)
   gsn_text_ndc(wks,"b) 30-yr Running Temperature Trends + 95% CI",0.1,0.765,txres)
   gsn_text_ndc(wks,"c) Seasonal Time-Varying Surface Temperature Trends",0.1,0.555,txres)
   gsn_text_ndc(wks,"d) Annual Mean Height and Temperature Trends",0.1,0.33,txres)   
  
  txres@txFontHeightF = 0.011
  txres@txAngleF = 90
  txres@txJust = "CenterCenter"
  gsn_text_ndc(wks,"Standardized",0.03,0.9,txres)
  gsn_text_ndc(wks,"Temperature",0.045,0.9,txres)
  gsn_text_ndc(wks,"30-yr Trend",0.03,0.7,txres)
  gsn_text_ndc(wks,"(~F34~0~F~C decade~S~-1~N~)",0.045,0.7,txres)
;  res3@tiYAxisString	 = "30-yr Running Trend (~F34~0~F~C / decade)"		 ; add a axis title    
  txres@txFontHeightF = resPanel@lbLabelFontHeightF
  txres@txAngleF = 0
  txres@txJust = "CenterRight"
  gsn_text_ndc(wks,"Temperature (~F34~0~F~C)",0.28,0.06,txres)
  gsn_text_ndc(wks,"Height (m)",0.28,0.03,txres)

  txres@txJust = "CenterCenter"
  gsn_text_ndc(wks,"-40",0.3,0.03,txres)
  gsn_text_ndc(wks,"-30",0.35,0.03,txres)
  gsn_text_ndc(wks,"-20",0.4,0.03,txres)
  gsn_text_ndc(wks,"-10",0.45,0.03,txres)
  gsn_text_ndc(wks,"0",0.5,0.03,txres)
  gsn_text_ndc(wks,"40",0.7,0.03,txres)
  gsn_text_ndc(wks,"30",0.65,0.03,txres)
  gsn_text_ndc(wks,"20",0.6,0.03,txres)
  gsn_text_ndc(wks,"10",0.55,0.03,txres)

end
