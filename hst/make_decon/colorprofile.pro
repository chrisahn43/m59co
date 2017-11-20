PRO COLORPROFILE
  scale=0.05
  radius=findgen(42)+1         ;(10^(findgen(40)*0.05+0.05))
  minradius=findgen(42);[0,(10^(findgen(39)*0.05+0.05))]
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  eg=0.107
  ez=0.041
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_475_galfit_model.fits',gimg,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_850_galfit_model.fits',zimg,exten_no=1
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_475_galfit_model.fits',gfree,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_475_galfit_modelfixed.fits',gfixed,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_850_galfit_model.fits',zfree,exten_no=2
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_850_galfit_modelfixed.fits',zfixed,exten_no=2
  find_galaxy,gimg,m,e,a,xcg,ycg
  find_galaxy,zimg,m,e,a,xcz,ycz

  fits_read,'serfixmodelg.fits',serfixg
  fits_read,'serfreemodelg.fits',serfreeg
  fits_read,'serfixmodelz.fits',serfixz
  fits_read,'serfreemodelz.fits',serfreez
  find_galaxy,serfreeg,m,e,a,xcg_m,ycg_m
  find_galaxy,serfixz,m,e,a,xcz_m,ycz_m
  fits_read,'mgedirectfit_475.fits',mgeg
  fits_read,'mgedirectfit_850.fits',mgez
  find_galaxy,mgeg,m,e,a,xcg_mge,ycg_mge
  find_galaxy,mgez,m,e,a,xcz_mge,ycz_mge
  gdata=fltarr(n_elements(radius))
  zdata=fltarr(n_elements(radius))
  gmagfixed=fltarr(n_elements(radius))
  zmagfixed=fltarr(n_elements(radius))
  gmagfree=fltarr(n_elements(radius))
  zmagfree=fltarr(n_elements(radius))
  gmagserfixed=fltarr(n_elements(radius))
  gmagserfree=fltarr(n_elements(radius))
  zmagserfixed=fltarr(n_elements(radius))
  zmagserfree=fltarr(n_elements(radius))
  gmagmge=fltarr(n_elements(radius))
  zmagmge=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))

     aper,gimg,xcg,ycg,maxgdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gimg,xcg,ycg,mingdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,zimg,xcz,ycz,maxzdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zimg,xcz,ycz,minzdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,gfree,xcg,ycg,maxgfree,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gfree,xcg,ycg,mingfree,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,gfixed,xcg,ycg,maxgfixed,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gfixed,xcg,ycg,mingfixed,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,zfree,xcz,ycz,maxzfree,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zfree,xcz,ycz,minzfree,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,zfixed,xcz,ycz,maxzfixed,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zfixed,xcz,ycz,minzfixed,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfixg,xcg_m,ycg_m,maxgfixser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixg,xcg_m,ycg_m,mingfixser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfreeg,xcg_m,ycg_m,maxgfreeser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreeg,xcg_m,ycg_m,mingfreeser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfixz,xcz_m,ycz_m,maxzfixser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixz,xcz_m,ycz_m,minzfixser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,serfreez,xcz_m,ycz_m,maxzfreeser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreez,xcz_m,ycz_m,minzfreeser,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,mgeg,xcg_mge,ycg_mge,maxgmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgeg,xcg_mge,ycg_mge,mingmgeflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,mgez,xcz_mge,ycz_mge,maxzmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgez,xcz_mge,ycz_mge,minzmgeflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        gflux=maxgdata
        zflux=maxzdata
        gfreeflux=maxgfree
        gfixedflux=maxgfixed
        zfreeflux=maxzfree
        zfixedflux=maxzfixed
        gfixserflux=maxgfixser
        gfreeserflux=maxgfreeser
        zfixserflux=maxzfixser
        zfreeserflux=maxzfreeser
        gmgeflux=maxgmgeflux
        zmgeflux=maxzmgeflux
     endif else begin
        gflux=maxgdata-mingdata
        zflux=maxzdata-minzdata
        gfreeflux=maxgfree-mingfree
        gfixedflux=maxgfixed-mingfixed
        zfreeflux=maxzfree-minzfree
        zfixedflux=maxzfixed-minzfixed
        gfixserflux=maxgfixser-mingfixser
        gfreeserflux=maxgfreeser-mingfreeser
        zfixserflux=maxzfixser-minzfixser
        zfreeserflux=maxzfreeser-minzfreeser
        gmgeflux=maxgmgeflux-mingmgeflux
        zmgeflux=maxzmgeflux-minzmgeflux
     endelse
     gdata[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gflux)-eg
     zdata[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zflux)-ez
     gmagfixed[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gfixedflux)-eg
     zmagfixed[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zfixedflux)-ez
     gmagfree[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gfreeflux)-eg
     zmagfree[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zfreeflux)-ez
     gmagserfixed[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gfixserflux)-eg
     gmagserfree[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gfreeserflux)-eg
     zmagserfixed[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zfixserflux)-ez
     zmagserfree[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zfreeserflux)-ez
     gmagmge[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gmgeflux)
     zmagmge[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zmgeflux)
  endfor

  colordata=gdata-zdata
  colorfree=gmagfree-zmagfree
  gcolorfixed=gmagfixed-zmagfree
  zcolorfixed=gmagfree-zmagfixed
  sercolorfree=gmagserfree-zmagserfree
  gsercolorfixed=gmagserfixed-zmagserfree
  zsercolorfixed=gmagserfree-zmagserfixed
  mgecolor=gmagmge-zmagmge
  set_plot,'ps'
  device,filename='color_all.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;  myplot,filename='color_all.ps'
  djs_plot,radius*scale,colordata,ytitle='(\mu_{F475W} - \mu_{F850LP}) [Mag/asec^2]',xtitle='Radius ["]',yran=[1.3,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,xstyle=8,ymargin=[4,4],psym=4,xran=[0,1.5]
  djs_oplot,radius*scale,colorfree,thick=8
  djs_oplot,radius*scale,sercolorfree,thick=8,linestyle=2
  djs_oplot,radius*scale,gcolorfixed,thick=8,color='blue'
  djs_oplot,radius*scale,gsercolorfixed,thick=8,color='blue',linestyle=2
  djs_oplot,radius*scale,zcolorfixed,thick=8,color='red'
  djs_oplot,radius*scale,zsercolorfixed,thick=8,color='red',linestyle=2
;  djs_oplot,radius*scale,mgecolor,thick=4,color='green',linestyle=2
  rad=findgen(42)+1
  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=4,xtitle='Radius [Pixels]',/xsty,xran=[1,30]
  items=['Data','Convolved Model','Unconvolved Model'];,'Direct MGE Fit']
  lines=[0,0,2];,2]
  sym=[4,0,0];,0]
  colors=['black','black','black'];,'green']
  al_legend,items,linestyle=lines,psym=sym,color=colors,/window,background_color='white',charthick=4,thick=3,/bottom,/right ;,charsize=1.5
  xyouts,[0.02],[1.57],['M59cO'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  stop

END

PRO CROSS_CON

  scale=0.05
  radius=findgen(42)+1
  minradius=findgen(42)
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  eg=0.107
  ez=0.041
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_475_skysubtract.fits',gimg
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_850_skysubtract.fits',zimg

  xcg=315 & ycg=298
  xcz=315 & ycz=298
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/psf/temp_psf.fits',gpsf
  gpsf=gpsf/TOTAL(gpsf)
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/psf/temp_psf_850lp.fits',zpsf
  zpsf=zpsf/TOTAL(zpsf)
  gconv=convolve(gimg,zpsf)
  zconv=convolve(zimg,gpsf)
  xcg_c=324 & ycg_c=296
  xcz_c=325 & ycz_c=295
  stop

  gdata=fltarr(n_elements(radius))
  zdata=fltarr(n_elements(radius))
  gcross=fltarr(n_elements(radius))
  zcross=fltarr(n_elements(radius))

  for i=0,n_elements(radius)-1 do begin
     area=((!PI*(radius[i])^2)-(!PI*(minradius[i])^2))
     
     aper,gimg,xcg,ycg,maxgdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gimg,xcg,ycg,mingdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     
     aper,zimg,xcz,ycz,maxzdata,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zimg,xcz,ycz,minzdata,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,gconv,xcg_c,ycg_c,maxgconv,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gconv,xcg_c,ycg_c,mingconv,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     
     aper,zconv,xcz_c,ycz_c,maxzconv,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zconv,xcz_c,ycz_c,minzconv,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minradius[i] lt 0.025) then begin
        gflux=maxgdata
        zflux=maxzdata
        gconvflux=maxgconv
        zconvflux=maxzconv
     endif else begin
        gflux=maxgdata-mingdata
        zflux=maxzdata-minzdata
        gconvflux=maxgconv-mingconv
        zconvflux=maxzconv-minzconv
     endelse
     gdata[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gflux)-eg
     zdata[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zflux)-ez
     gcross[i]=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gconvflux)-eg
     zcross[i]=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zconvflux)-ez
     
  endfor
  
  colordata=gdata-zdata
  colorconv=gcross-zcross

  set_plot,'ps'
  device,filename='crosscolor.ps',/color
  djs_plot,radius*scale,colordata,ytitle='(\mu_{F475W} - \mu_{F850LP}) [mag/arcsec^2]',xtitle='Radius ["]',yran=[1.3,1.6],/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,thick=3,xstyle=8,ymargin=[4,4],psym=4,xran=[0,2.]
  djs_oplot,radius*scale,colorconv,thick=3,psym=4,color='red',linestyle=2
  rad=findgen(42)+1
  axis,xaxis=1,xtickv=rad,xcharsize=1.5,charthick=4,xthick=4,xtitle='Radius [Pixels]',/xsty,xran=[1,42]
  items=['Data','Cross Convolved']
  lines=[0,0]
  sym=[4,4]
  color=['black','red']
  al_legend,items,linestyle=lines,psym=sym,colors=color,/window,background_color='white',charthick=4,thick=3,/bottom,/right;,charsize=1.5
  device,/close
  set_plot,'x'
  stop

END



  
