pro final_color
  scale=0.05
  rad=findgen(60)+1
  minrad=findgen(60)
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  fits_read,'./galfit_output/m59co_475_galfit_model.fits',gimg,exten_no=1
  fits_read,'./galfit_output/m59co_850_galfit_model.fits',zimg,exten_no=1
  fits_read,'./galfit_output/m59co_475_galfit_model.fits',gfree,exten_no=2
  fits_read,'./galfit_output/m59co_475_galfit_modelfixed.fits',gfixed,exten_no=2
  fits_read,'./galfit_output/m59co_850_galfit_model.fits',zfree,exten_no=2
  fits_read,'./galfit_output/m59co_850_galfit_modelfixed.fits',zfixed,exten_no=2
  xc=59 & yc=48
  gdata=fltarr(n_elements(rad))
  zdata=fltarr(n_elements(rad))
  gmagfixed=fltarr(n_elements(rad))
  zmagfixed=fltarr(n_elements(rad))
  gmagfree=fltarr(n_elements(rad))
  zmagfree=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     area=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))
     maxgdata=djs_phot(xc,yc,rad[i],0.,gimg,skyval=skyval)
     mingdata=djs_phot(xc,yc,minrad[i],0.,gimg,skyval=skyval)
     gflux=maxgdata-mingdata
     gdata[i]=zeropoint_g + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(gflux)
     maxzdata=djs_phot(xc,yc,rad[i],0.,zimg,skyval=skyval)
     minzdata=djs_phot(xc,yc,minrad[i],0.,zimg,skyval=skyval)
     zflux=maxzdata-minzdata
     zdata[i]=zeropoint_z + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(zflux)
     
     maxgfree=djs_phot(xc,yc,rad[i],0.,gfree,skyval=skyval)
     mingfree=djs_phot(xc,yc,minrad[i],0.,gfree,skyval=skyval)
     gfluxfree=maxgfree-mingfree
     gmagfree[i]=zeropoint_g + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(gfluxfree)
     
     maxzfree=djs_phot(xc,yc,rad[i],0.,zfree,skyval=skyval)
     minzfree=djs_phot(xc,yc,minrad[i],0.,zfree,skyval=skyval)
     zfluxfree=maxzfree-minzfree
     zmagfree[i]=zeropoint_z + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(zfluxfree)
     
     maxgfixed=djs_phot(xc,yc,rad[i],0.,gfixed,skyval=skyval)
     mingfixed=djs_phot(xc,yc,minrad[i],0.,gfixed,skyval=skyval)
     gfluxfixed=maxgfixed-mingfixed
     gmagfixed[i]=zeropoint_g + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(gfluxfixed)
     
     maxzfixed=djs_phot(xc,yc,rad[i],0.,zfixed,skyval=skyval)
     minzfixed=djs_phot(xc,yc,minrad[i],0.,zfixed,skyval=skyval)
     zfluxfixed=maxzfixed-minzfixed
     zmagfixed[i]=zeropoint_z + 5*alog10(scale) + 2.5*alog10(area) - 2.5*alog10(zfluxfixed)
     
     
     
  endfor
  colordata=gdata-zdata
  colorfree=gmagfree-zmagfree
  gcolorfixed=gmagfixed-zmagfree
  zcolorfixed=gmagfree-zmagfixed
  stop
                                ;GFREE
 ; sersic    : (  307.47,   301.85)   19.40      3.13    1.06    0.97   -65.24
 ; sersic    : ( {307.47}, {301.85})  18.38     12.85    1.09    0.98    88.42
  bn1=1.999*1.06-0.327
  bn2=1.999*1.09-0.327
  mue1=19.4+5*alog10(3.13*scale)+2.5*alog10(2*!PI*1.06*(exp(bn1)/((bn1)^(2*1.06)))*GAMMA(2*1.06))
  mue2=18.38+5*alog10(12.85*scale)+2.5*alog10(2*!PI*1.09*(exp(bn2)/((bn2)^(2*1.09)))*GAMMA(2*1.09))

  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(3.13))^(1./1.06))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(12.85))^(1./1.09))-1)
  int1=(10^(-0.4*(mu1-zeropoint_g)))
  int2=(10^(-0.4*(mu2-zeropoint_g)))
  inttot=int1+int2
  mutot_gfree=zeropoint_g-2.5*alog10(inttot)

                                ;ZFREE
; sersic    : (  308.57,   300.97)   18.17      2.96    1.02    0.99    34.06
; sersic    : ( {308.57}, {300.97})  16.72     12.25    1.21    0.98    17.69
  bn1850=1.999*1.02-0.327
  bn2850=1.999*1.21-0.327
  mue1850=18.17+5*alog10(2.96*scale)+2.5*alog10(2*!PI*1.02*(exp(bn1850)/((bn1850)^(2*1.02)))*GAMMA(2*1.02))
  mue2850=16.72+5*alog10(12.25*scale)+2.5*alog10(2*!PI*1.21*(exp(bn2850)/((bn2850)^(2*1.21)))*GAMMA(2*1.21))

  mu1850=mue1850+((2.5*bn1850)/(alog(10)))*((((rad)/(2.96))^(1./1.02))-1)
  mu2850=mue2850+((2.5*bn2850)/(alog(10)))*((((rad)/(12.25))^(1./1.21))-1)
  int1850=(10^(-0.4*(mu1850-zeropoint_z)))
  int2850=(10^(-0.4*(mu2850-zeropoint_z)))
  inttot850=int1850+int2850
  mutot_zfree=zeropoint_z-2.5*alog10(inttot850)
  colormu_free=mutot_gfree-mutot_zfree

  set_plot,'ps'
  device,filename='color_all.ps',/color
  djs_plot,rad*scale,colormu_free,ytitle='(\mu_{F475W} - \mu_{F850LP}) [mag/arcsec^2]',xtitle='Radius ["]',yran=[1.3,2.2],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,xstyle=8,ymargin=[4,4],linestyle=2
  djs_oplot,rad*scale,colorfree,thick=3
  djs_oplot,rad*scale,colordata,thick=3,psym=4

                                ;GFIXED
;  sersic    : (  307.47,   301.85)   19.56     [2.96]  [1.02]  [0.99]  [34.06]
;  sersic    : ( {307.47}, {301.85})  18.32    [12.25]  [1.21]  [0.98]  [17.69]
  bn1=1.999*1.02-0.327
  bn2=1.999*1.21-0.327
  mue1=19.56+5*alog10(2.96*scale)+2.5*alog10(2*!PI*1.02*(exp(bn1)/((bn1)^(2*1.02)))*GAMMA(2*1.02))
  mue2=18.32+5*alog10(12.25*scale)+2.5*alog10(2*!PI*1.21*(exp(bn2)/((bn2)^(2*1.21)))*GAMMA(2*1.21))

  mu1=mue1+((2.5*bn1)/(alog(10)))*((((rad)/(2.96))^(1./1.02))-1)
  mu2=mue2+((2.5*bn2)/(alog(10)))*((((rad)/(12.25))^(1./1.21))-1)
  int1=(10^(-0.4*(mu1-zeropoint_g)))
  int2=(10^(-0.4*(mu2-zeropoint_g)))
  inttot=int1+int2
  mutot_gfixed=zeropoint_g-2.5*alog10(inttot)

  colormu_gfixed=mutot_gfixed-mutot_zfree
  djs_oplot,rad*scale,colormu_gfixed,color='blue',thick=3,linestyle=2
  djs_oplot,rad*scale,gcolorfixed,thick=3,color='blue'
  
                                ;ZFIXED
; sersic    : (  308.57,   300.98)   17.98     [3.13]  [1.06]  [0.97] [-65.24]
; sersic    : ( {308.57}, {300.98})  16.78    [12.85]  [1.09]  [0.98]  [88.42]
  bn1850=1.999*1.06-0.327
  bn2850=1.999*1.09-0.327
  mue1850=17.98+5*alog10(3.13*scale)+2.5*alog10(2*!PI*1.06*(exp(bn1850)/((bn1850)^(2*1.06)))*GAMMA(2*1.06))
  mue2850=16.78+5*alog10(12.85*scale)+2.5*alog10(2*!PI*1.09*(exp(bn2850)/((bn2850)^(2*1.09)))*GAMMA(2*1.09))

  mu1850=mue1850+((2.5*bn1850)/(alog(10)))*((((rad)/(3.13))^(1./1.06))-1)
  mu2850=mue2850+((2.5*bn2850)/(alog(10)))*((((rad)/(12.85))^(1./1.09))-1)
  int1850=(10^(-0.4*(mu1850-zeropoint_z)))
  int2850=(10^(-0.4*(mu2850-zeropoint_z)))
  inttot850=int1850+int2850
  mutot_zfixed=zeropoint_z-2.5*alog10(inttot850)

  colormu_zfixed=mutot_gfree-mutot_zfixed
  djs_oplot,rad*scale,colormu_zfixed,color='red',thick=3,linestyle=2
  djs_oplot,rad*scale,zcolorfixed,thick=3,color='red'

  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=3,xtitle='Radius [Pixels]',/xsty,xran=[1,60]
  items=['Data','Convolved Model','Unconvolved Model']
  lines=[0,0,2]
  sym=[4,0,0]
  al_legend,items,linestyle=lines,psym=sym,/window,background_color='white',charthick=4,thick=3,/bottom,/right;,charsize=1.5

  device,/close
  set_plot,'x'
  stop

  
END

pro surface_profile
  scale=0.05
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  fits_read,'./galfit_output/m59co_475_galfit_model.fits',gimg,exten_no=1
  fits_read,'./galfit_output/m59co_850_galfit_model.fits',zimg,exten_no=1
  fits_read,'./galfit_output/m59co_475_galfit_model.fits',gfree,exten_no=2
  fits_read,'./galfit_output/m59co_475_galfit_modelfixed.fits',gfixed,exten_no=2
  fits_read,'./galfit_output/m59co_850_galfit_model.fits',zfree,exten_no=2
  fits_read,'./galfit_output/m59co_850_galfit_modelfixed.fits',zfixed,exten_no=2

  imgsize=size(gimg,/dim)

  gmag=fltarr(imgsize[0],imgsize[1])
  zmag=fltarr(imgsize[0],imgsize[1])
  color=fltarr(imgsize[0],imgsize[1])
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if (gimg[i,j] gt 0.) then begin
           gmag[i,j]=zeropoint_g+5*alog10(scale)-2.5*alog10(gimg[i,j])
        endif
        if (zimg[i,j] gt 0.) then begin
           zmag[i,j]=zeropoint_z+5*alog10(scale)-2.5*alog10(zimg[i,j])
        endif
        color[i,j]=gmag[i,j]-zmag[i,j]
     endfor
  endfor

  
  x=(findgen(101)) # (fltarr(101)+1)
  y=(fltarr(101)+1) # (fltarr(101))
  klevels=[0.016,0.04,0.1,0.25,0.65]*max(color)
  set_plot,'ps'
  device,get_decomposed=currentmode
  device,decomposed=0
  loadct,33,ncolors=20,bottom=1
  device,filename='2dcolor.ps',/color
  contour,color,x,y,levels=klevel,/fill,/xs,/yc,charsize=2.5,/nodata,xtitle='X offset ["]',ytitle='Y offset ["]',xmargin=[6,8],ymargin=[3.5,3.5]
  colorrange=max(color)-min(color)
  outcol=FIX((color-min(color))/colorrange*254)+1
  loadct,34
  tvlct,r,g,b,/get
  r[0]=0 & g[0]=0 & b[0]=0
  tvlct,r,g,b
  imgunder,outcol
  contour,color,x,y,/over,levels=klevels,c_color=0,c_thick=4
  device,/close
  set_plot,'x'
  stop
END

