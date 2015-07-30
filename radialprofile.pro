pro radialprofile
  rad=40.
  fits_read, 'data/m59_g_lucy.fits', img_gd, h_gd
  fit_g=gauss2dfit(img_gd,gparam)
  xc_gd=gparam[4]
  yc_gd=gparam[5]
  undefine,gparam
  xmin=fix(xc_gd-rad) & xmax=fix(xc_gd+rad+0.5)
  ymin=fix(yc_gd-rad) & ymax=fix(yc_gd+rad+0.5)
  img_gd=img_gd[xmin:xmax,ymin:ymax]
  fit_g=gauss2dfit(img_gd,gparam)
  xc_gd=gparam[4]
  yc_gd=gparam[5]
  img_gd=img_gd/375.

  fits_read, 'data/m59_z_lucy.fits',img_zd,h_zd
  fit_z=gauss2dfit(img_zd,zparam)
  xc_zd=zparam[4]
  yc_zd=zparam[5]
  undefine,zparam
  xmin=fix(xc_zd-rad) & xmax=fix(xc_zd+rad+0.5)
  ymin=fix(yc_zd-rad) & ymax=fix(yc_zd+rad+0.5)  
  img_zd=img_zd[xmin:xmax,ymin:ymax]
  fit_z=gauss2dfit(img_zd,zparam)
  xc_zd=zparam[4]
  yc_zd=zparam[5]
  img_zd=img_zd/560.

  fits_read, 'data/HST_9401_09_ACS_WFC_F475W_drz.fits',img_go,h_go
  img_go=img_go[4260:4372,3942:4054]
  fit_go=gauss2dfit(img_go,gparamo)
  xc_go=gparamo[4]
  yc_go=gparamo[5]
  undefine,gparamo
  xmin=fix(xc_go-rad) & xmax=fix(xc_go+rad+0.5)
  ymin=fix(yc_go-rad) & ymax=fix(yc_go+rad+0.5)
  img_go=img_go[xmin:xmax,ymin:ymax]
  fit_go=gauss2dfit(img_go,gparamo)
  xc_go=gparamo[4]
  yc_go=gparamo[5]
  img_go=img_go/750.

  fits_read, 'data/HST_9401_09_ACS_WFC_F850LP_drz.fits',img_zo,h_zo
  img_zo=img_zo[4260:4372,3942:4054]
  fit_zo=gauss2dfit(img_zo,zparamo)
  xc_zo=zparamo[4]
  yc_zo=zparamo[5]
  undefine,zparamo
  xmin=fix(xc_zo-rad) & xmax=fix(xc_zo+rad+0.5)
  ymin=fix(yc_zo-rad) & ymax=fix(yc_zo+rad+0.5)  
  img_zo=img_zo[xmin:xmax,ymin:ymax]
  fit_zo=gauss2dfit(img_zo,zparamo)
  xc_zo=zparamo[4]
  yc_zo=zparamo[5]
  img_zo=img_zo/1210.

  scale=0.05
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  rad=findgen(25)+1
  minrad=findgen(25)
  mag_gd=fltarr(n_elements(rad))
  mag_zd=fltarr(n_elements(rad))
  mag_go=fltarr(n_elements(rad))
  mag_zo=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     maxflux_gd=djs_phot(xc_gd,yc_gd,rad[i],0.,img_gd,skyval=skyval)
     minflux_gd=djs_phot(xc_gd,yc_gd,minrad[i],0.,img_gd,skyval=skyval)
     flux_gd=maxflux_gd-minflux_gd
     maxflux_zd=djs_phot(xc_zd,yc_zd,rad[i],0.,img_zd,skyval=skyval)
     minflux_zd=djs_phot(xc_zd,yc_zd,minrad[i],0.,img_zd,skyval=skyval)
     flux_zd=maxflux_zd-minflux_zd
     maxflux_go=djs_phot(xc_go,yc_go,rad[i],0.,img_go,skyval=skyval)
     minflux_go=djs_phot(xc_go,yc_go,minrad[i],0.,img_go,skyval=skyval)
     flux_go=maxflux_go-minflux_go
     maxflux_zo=djs_phot(xc_zo,yc_zo,rad[i],0.,img_zo,skyval=skyval)
     minflux_zo=djs_phot(xc_zo,yc_zo,minrad[i],0.,img_zo,skyval=skyval)
     flux_zo=maxflux_zo-minflux_zo
     mag_gd[i]=zeropoint_g + 5*alog10(scale) - 2.5*alog10(flux_gd)
     mag_zd[i]=zeropoint_z + 5*alog10(scale) - 2.5*alog10(flux_zd)
     mag_go[i]=zeropoint_g + 5*alog10(scale) - 2.5*alog10(flux_go)
     mag_zo[i]=zeropoint_z + 5*alog10(scale) - 2.5*alog10(flux_zo)


  endfor
  
  colormagd=mag_gd-mag_zd
  colormago=mag_go-mag_zo

  djs_plot,rad,mag_gd,psym=2,ytitle='Magnitude',xtitle='radius',title='Deconvolved Images'
  djs_oplot,rad,mag_zd,psym=2,color='blue'
  stop
  djs_plot,rad,mag_go,psym=2,ytitle='Magnitude',xtitle='radius',title='Original Images',color='green'
  djs_oplot,rad,mag_zo,psym=2,color='red'
  stop
  djs_plot,rad,mag_gd,psym=2,ytitle='Magnitude',xtitle='radius',title='All images',yran=[5,30]
  djs_oplot,rad,mag_zd,psym=2,color='blue'
  djs_oplot,rad,mag_go,psym=2,color='green'
  djs_oplot,rad,mag_zo,psym=2,color='red'
  stop
  djs_plot,rad,colormagd,psym=2,ytitle='Magnitdue',xtitle='radius',title='Color of Deconvolved and Original'
  djs_oplot,rad,colormago,psym=2,color='blue'
  stop
  

END
