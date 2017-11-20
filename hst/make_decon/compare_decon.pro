PRO COMPARE_DECON
  scale=0.05
  zeropt=26.05923
  extinct=0.107
  fits_read,'mgefixmodelg.fits',fixmodelg
  fits_read,'serfixmodelg.fits',fixmodelgser
  fits_read,'mgefreemodelg.fits',freemodelg
  fits_read,'serfreemodelg.fits',freemodelgser
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xc,yc
  deconimg=img/375.

  find_galaxy,fixmodelg,m,e,a,xc_m,yc_m

    radius=(10^(findgen(40)*0.05+0.05))
  minradius=[0,(10^(findgen(39)*0.05+0.05))]
  deconmag=fltarr(n_elements(radius))
  modelmag=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  frmodelmag=fltarr(n_elements(radius))
  frsersicmag=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,deconimg,xc,yc,maxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,deconimg,xc,yc,minflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelg,xc_m,yc_m,maxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelg,xc_m,yc_m,minmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelgser,xc_m,yc_m,maxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelgser,xc_m,yc_m,minmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelg,xc_m,yc_m,fmaxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelg,xc_m,yc_m,fminmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelgser,xc_m,yc_m,fmaxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelgser,xc_m,yc_m,fminmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] lt 0.025) then begin
        flux=maxflux
        model=maxmodel
        sermodel=maxmodelser
        fmodel=fmaxmodel
        fsermodel=fmaxmodelser
     endif else begin
        flux=maxflux-minflux
        model=maxmodel-minmodel
        sermodel=maxmodelser-minmodelser
        fmodel=fmaxmodel-fminmodel
        fsermodel=fmaxmodelser-fminmodelser
     endelse

     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     deconmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(flux)-extinct
     modelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(model);-extinct
     sersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(sermodel) -extinct
     frmodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fmodel) ;-extinct
     frsersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fsermodel) -extinct

  endfor

  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='deconvolved_vs_sersic_mge_475.ps',/color
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. MGE fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frmodelmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,frsersicmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-sersicmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frsersicmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Fixed Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_plot,radius*scale,sersicmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
;  stop
  djs_plot,radius*scale,frsersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Free Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,frsersicmag-frmodelmag,psym=2,color='red',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  device,/close
  set_plot,'x'
  stop
  zeropt=24.84245
  extinct=0.041
  fits_read,'mgefixmodelz.fits',fixmodelz
  fits_read,'serfixmodelz.fits',fixmodelzser
  fits_read,'mgefreemodelz.fits',freemodelz
  fits_read,'serfreemodelz.fits',freemodelzser
  hrotate,oldimg,head,img,newhead,1
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_z_lucy.fits',oldimgz,headz
  hrotate,oldimgz,headz,imgz,newheadz,1
  find_galaxy,imgz,m,e,a,xc,yc
  deconimg=imgz/560.

  find_galaxy,fixmodelz,m,e,a,xc_m,yc_m
  radius=(10^(findgen(40)*0.05+0.05))
  minradius=[0,(10^(findgen(39)*0.05+0.05))]
  deconmag=fltarr(n_elements(radius))
  modelmag=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  frmodelmag=fltarr(n_elements(radius))
  frsersicmag=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,deconimg,xc,yc,maxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,deconimg,xc,yc,minflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelz,xc_m,yc_m,maxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelz,xc_m,yc_m,minmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,fixmodelzser,xc_m,yc_m,maxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,fixmodelzser,xc_m,yc_m,minmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelz,xc_m,yc_m,fmaxmodel,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelz,xc_m,yc_m,fminmodel,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,freemodelzser,xc_m,yc_m,fmaxmodelser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,freemodelzser,xc_m,yc_m,fminmodelser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] lt 0.025) then begin
        flux=maxflux
        model=maxmodel
        sermodel=maxmodelser
        fmodel=fmaxmodel
        fsermodel=fmaxmodelser
     endif else begin
        flux=maxflux-minflux
        model=maxmodel-minmodel
        sermodel=maxmodelser-minmodelser
        fmodel=fmaxmodel-fminmodel
        fsermodel=fmaxmodelser-fminmodelser
     endelse

     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     deconmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(flux)-extinct
     modelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(model);-extinct
     sersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(sermodel) -extinct
     frmodelmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fmodel) ;-extinct
     frsersicmag[i]=zeropt+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(fsermodel) -extinct

  endfor
  !P.MULTI=[0,1,2]
  set_plot,'ps'
  device,filename='deconvolved_vs_sersic_mge_850.ps',/color
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. MGE fits',charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frmodelmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,deconmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Deconvolved Image vs. Sersic fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,sersicmag,color='blue',thick=3
  djs_oplot,radius*scale,frsersicmag,color='red',thick=3
  djs_plot,radius*scale,deconmag-sersicmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.35,0.6] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  djs_oplot,radius*scale,deconmag-frsersicmag,psym=2,color='red'
;  stop
  djs_plot,radius*scale,sersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Fixed Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,modelmag,color='blue',thick=3
  djs_plot,radius*scale,sersicmag-modelmag,psym=2,color='blue',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
;  stop
  djs_plot,radius*scale,frsersicmag,yran=[25,13],/ysty,xtitle='Radius ["]',ytitle='\mu [Mag/square arcsecond]',ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],xran=[0,1.5],/xsty,psym=2,title='Free Sersic Fits vs. MGE fits' ,charthick=4,xthick=3,ythick=3,thick=3
  djs_oplot,radius*scale,frmodelmag,color='red',thick=3
  djs_plot,radius*scale,frsersicmag-frmodelmag,psym=2,color='red',ymargin=[4,0],position=[0.1,0.05,0.95,0.25],xtitle='Radius ["]',ytitle='Residual',ycharsize=0.6,xran=[0,1.5],/xsty,yran=[-0.05,0.05] ,charthick=4,xthick=3,ythick=3,charsize=1.5
  device,/close
  set_plot,'x'
  !P.MULTI=[0,1,1]
  stop


  
END
