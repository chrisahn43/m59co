pro colormap_orig
  rad=40.
  fits_read, 'data/HST_9401_09_ACS_WFC_F475W_drz.fits',img_g,h_g
  img_g=img_g[4260:4372,3942:4054]
  fit_g=gauss2dfit(img_g,gparam)
  xc_g=gparam[4]
  yc_g=gparam[5]
  xmin=fix(xc_g-rad) & xmax=fix(xc_g+rad+0.5)
  ymin=fix(yc_g-rad) & ymax=fix(yc_g+rad+0.5)
  img_g=img_g[xmin:xmax,ymin:ymax]
  img_g=img_g/750.
  sizeg=size(img_g,/dim)

  fits_read, 'data/HST_9401_09_ACS_WFC_F850LP_drz.fits',img_z,h_z
  img_z=img_z[4260:4372,3942:4054]
  fit_z=gauss2dfit(img_z,zparam)
  xc_z=zparam[4]
  yc_z=zparam[5]
  xmin=fix(xc_z-rad) & xmax=fix(xc_z+rad+0.5)
  ymin=fix(yc_z-rad) & ymax=fix(yc_z+rad+0.5)  
  img_z=img_z[xmin:xmax,ymin:ymax]
  img_z=img_z/1210.
  sizez=size(img_z,/dim)
  flratio=img_g/img_z
  writefits,'fluxratio_orig.fits',flratio
  scale=0.05
  zeropoint_g = 26.05923
  zeropoint_z = 24.84245
  mapg = fltarr(sizeg[0],sizeg[1])
  mapz = fltarr(sizez[0],sizez[1])
  for i=0,sizeg[0]-1 do begin
     for j=0,sizeg[1]-1 do begin
        mag_g = zeropoint_g + 5*alog10(scale) - 2.5*alog10(img_g[i,j])
        mapg[i,j] = mag_g
        mag_z = zeropoint_z + 5*alog10(scale) - 2.5*alog10(img_z[i,j])
        mapz[i,j] = mag_z
     endfor
  endfor
  color=mapg - mapz
  size=size(color,/dim)
  xcoord=(findgen(size[0])-rad)*0.05
  ycoord=(findgen(size[1])-rad)*0.05
  klevels=[0.016,0.04,0.1,0.25,0.65]*max(color)
;  device, get_decomposed=currentmode
;  device, decomposed=0
;  loadct, 33, Ncolors=20, bottom=1
  contour,color,xcoord,ycoord,levels=klevels,/fill,/xs,/ys,charsize=2.5,xminor=2,/nodata,xtitle='X offset ["]',ytitle='Y offset ["]',xmargin=[6,8],ymargin=[3.5,3.5]
  colrange=max(color)-min(color)
  outcol=FIX((color-min(color))/colrange*254)+1
  loadct,34
  tvlct,r,g,b,/GET
  r[0]=0 & g[0]=0 & b[0]=0
  tvlct,r,g,b
  imgunder,outcol
  contour,color,xcoord,ycoord,/over,levels=klevels,c_color=0,c_thick=4


  writefits,'colormap_orig.fits',color
  stop


  
END
