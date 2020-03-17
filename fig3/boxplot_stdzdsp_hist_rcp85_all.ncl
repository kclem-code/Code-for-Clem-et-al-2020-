;----------------------------------------------------------------------
; box_6.ncl
;
; Concepts illustrated:
;   - Drawing box plots using statistics calculated with stat_dispersion
;   - Using text function codes to generate accented characters
;----------------------------------------------------------------------
; This script was contributed by Guilherme Martins of the
; National Institute for Space Research.
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

atxt     = "data_stdzdsp_hist_rcp85_all.txt"  ; Nome do arquivo no formato txt.
nrows  = numAsciiRow(atxt) ; NÃºmero de rows do arquivo.
ncols = numAsciiCol(atxt) ; NÃºmero de cols do arquivo.

f = asciiread(atxt,(/nrows,ncols/),"float") ; Abertura do arquivo texto. Esse arquivo
                                                 ; possui 78 rows e 4 cols.

vars = (/"All","All","Max","Max","SAM","Nino3.4","IPO","Weddell"/) ; Nome dos rÃ³tulos.

opt           = True ; Habilita personalizaÃ§Ã£o da estatÃ­stica.
opt@PrintStat = True ; Mostra na tela o resultado da estÃ¡tistica da funÃ§Ã£o stat_dispersion. 

; A linha abaixo imprime na tela o resultado da estatÃ­stica para cada variÃ¡vel.

y = new((/dimsizes(vars),nrows,ncols+1/),float)
x = ispan(1,ncols,1) ; Valores para o eixo x.

do j = 0,dimsizes(vars)-1
   g = stat_dispersion(f(:,j),opt)
   y(j,:,0) = g(2)  ; minimum value
   y(j,:,1) = g(6)  ; lower quartile
   y(j,:,2) = g(8)  ; median
   y(j,:,3) = g(10) ; upper quartile
   y(j,:,4) = g(14) ; max value
end do

wks = gsn_open_wks("pdf","box_stdzdsp_hist_rcp85_all")

res                 = True         ; PersonalizaÃ§Ã£o do boxplot.
res@vpHeightF= 0.5                    ; change aspect ratio of plot
res@vpWidthF = 0.85                  
res@tmXBLabels      = vars      ; RÃ³tulos do eixo x.
;res@tiMainString    = "CMIP5 PI-Control 30-Year South Pole Warming"
res@tiMainFontHeightF = 0.025
res@trYMinF         = -4.5          ; MÃ­nimo valor do eixo y.      
res@trYMaxF         = 4.5         ; MÃ¡ximo valor do eixo y.
res@tmYLMode        = "Manual"	   ; Eixo y esquerdo (YL) de forma manual.
res@tmYLTickStartF  = -4  ; Onde comeÃ§am os valores do eixo y.
res@tmYLTickEndF    = 4  ; Onde terminam os valores do eixo y.
res@tmYLTickSpacingF= 1.0            ; EspaÃ§amento entre os rÃ³tulos.
res@tmYLMinorOn     = False        ; Desliga o minortick do eixo y esquerdo.
res@tmYRBorderOn    = True         ; Habilita a borda do eixo y direito (YR). 
res@tmXTBorderOn    = True         ; Habilita a borda do eixo x superior (XB). 
res@tmXBLabelFontHeightF = 0.019   ; Tamanho da fonte do eixo x inferior (XB).
res@tmYLLabelFontHeightF = 0.019   ; Tamanho da fonte do eixo x inferior (XB).
res@tiYAxisFontHeightF = 0.022
res@tiYAxisString   = "Trend (30-Year~S~-1~N~)" ; TÃ­tulo do eixo y.
res@vpXF     = 0.12                   ; start plot at x ndc coord 


llres                   = True     ; PersonalizaÃ§Ã£o das rows do boxplot.
llres@gsLineThicknessF  = 5.0      ; Espessura da linha.

opti           = True ; Habilita personalizaÃ§Ã£o da caixa.
opti@boxWidth  = 0.5 ; Largura da caixa. Quanto maior o valor, mais larga serÃ¡ a caixa.
opti@boxColors = (/"red","red","red","red","black","blue","sea green","orange"/) ; Cores de cada caixa.
res@gsnMaximize = True

plot = boxplot(wks,x,y(:,dimsizes(vars),:),opti,res,llres) ; GeraÃ§Ã£o do boxplot.

; draw observed values as filled circles
  mres               = True                     ; marker mods desired
  mres@gsMarkerColor = "navy blue"              ; change color
  mres@gsMarkerIndex = 1                        ; change style
  mres@gsMarkerSizeF = 55.                      ; change size

  dum1 = gsn_add_polymarker(wks,plot,1.,2.39,mres) 
  dum2 = gsn_add_polymarker(wks,plot,2.,2.39,mres) 
  dum3 = gsn_add_polymarker(wks,plot,3.,2.39,mres) 
  dum4 = gsn_add_polymarker(wks,plot,4.,2.39,mres) 
  dum5 = gsn_add_polymarker(wks,plot,5.,2.14,mres)
  dum6 = gsn_add_polymarker(wks,plot,6.,-0.036,mres) 
  dum7 = gsn_add_polymarker(wks,plot,7.,-0.481,mres) 
  dum8 = gsn_add_polymarker(wks,plot,8.,-1.585,mres)


; draw multi-model mean values as asterisks
  mres@gsMarkerColor = "navy blue"              ; change color
  mres@gsMarkerIndex = 2                        ; change style
  mres@gsMarkerSizeF = 15.                      ; change size
  mres@gsMarkerThicknessF = 5

  dum9 = gsn_add_polymarker(wks,plot,1.,0.24,mres) 
  dum10 = gsn_add_polymarker(wks,plot,2.,0.01,mres) 
  dum11 = gsn_add_polymarker(wks,plot,3.,1.72,mres) 
  dum12 = gsn_add_polymarker(wks,plot,4.,1.90,mres) 
  dum13 = gsn_add_polymarker(wks,plot,5.,-1.15,mres) 
  dum14 = gsn_add_polymarker(wks,plot,6.,0.11,mres) 
  dum15 = gsn_add_polymarker(wks,plot,7.,-0.01,mres) 
  dum16 = gsn_add_polymarker(wks,plot,8.,0.36,mres)

; draw black line at zero axis (y=0)
  lres               = True                     ; marker mods desired
  lres@gsLineColor = "black"              ; change color
  lres@gsLineThicknessF = 3
  lx = (/ 0., 14./)
  ly = (/ 0., 0./)

  dum17 = gsn_add_polyline(wks,plot,lx,ly,lres) 

  lx = (/ 4.5, 4.5/)
  ly = (/ -4.5, 4.5/)

  dum18 = gsn_add_polyline(wks,plot,lx,ly,lres) 

; create legend
  mres@gsMarkerColor = "navy blue"              ; change color
  mres@gsMarkerIndex = 2                        ; change style
  mres@gsMarkerSizeF = .015                      ; change size
  mres@gsMarkerThicknessF = 5
  gsn_polymarker(wks,plot,2.6,-3.01,mres)

  mres@gsMarkerIndex = 1                        ; change style
  mres@gsMarkerSizeF = .04                      ; change size
  gsn_polymarker(wks,plot,2.6,-3.55,mres)

  tres               = True                   ; text mods desired
  tres@txFontHeightF = 0.016                   ; change font height
  gsn_text_ndc(wks,"ensemble mean",0.455,0.385, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"observed",0.423,0.355, tres)  ; add text in ndc coordinates

  tres@txFontHeightF = 0.018   
  gsn_text_ndc(wks,"S. Pole SAT Trends",0.355,0.775, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"Climate Index Trends",0.73,0.775, tres)  ; add text in ndc coordinates                ; change font height
  gsn_text_ndc(wks,"(Hist)",0.215,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.31,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(Hist)",0.403,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.497,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.592,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.687,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.782,0.238, tres)  ; add text in ndc coordinates
  gsn_text_ndc(wks,"(PI)",0.877,0.238, tres)  ; add text in ndc coordinates
;  gsn_text_ndc(wks,"(PI)",0.463,0.238, tres)  ; add text in ndc coordinates
;  gsn_text_ndc(wks,"(PI)",0.578,0.238, tres)  ; add text in ndc coordinates
;  gsn_text_ndc(wks,"(PI)",0.691,0.238, tres)  ; add text in ndc coordinates
;  gsn_text_ndc(wks,"(PI)",0.804,0.238, tres)  ; add text in ndc coordinates


draw(wks)
frame(wks)

end
