pro hst_m59co_850
  fits_read,'../data/HST_9401_09_ACS_WFC_F850LP_drz.fits',img,h
  fits_read,'sky_mask850lp.fits',sky
  fits_read,'../data/HST_9401_09_ACS_WFC_F850LP_drz.fits',wht,exten_no=2
  fits_read,'../data/HST_9401_09_ACS_WFC_F850LP_drz.fits',mask,exten_no=3
  mdrizz=[9.92686453499,50.4037169999,49.6092825856]
  expt=[90.,560.,560.]
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img=img+avgdrizz
  img=img[4000:4630,3700:4296]
  totcounts=img*(mean(expt))
  imgsize=size(img,/dim)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if (totcounts[i,j] lt 0.) then totcounts[i,j]=mean(totcounts)
     endfor
  endfor
  err=sqrt(totcounts)
  
  writefits,'m59co_weightmap_850.fits',err
;writefits,'m59co_mask.fits',mask
  img=(img-sky)
  writefits,'m59co_850_skysubtract.fits',img
  temp=img[266:366,249:349]
  ngauss=20
  minlevel=0.
  find_galaxy,temp,majoraxis,eps,ang,xc,yc,fraction=0.8
  
  sectors_photometry,temp,eps,ang,xc,yc,radius,angle,counts,minlevel=minlevel

  readcol,'tinytim_fit_850.dat',normpsf,sigmapsf,format='F,F'
  MGE_fit_sectors,radius,angle,counts,eps,sol=sol,ngauss=ngauss,scale=scale,normpsf=normpsf,sigmapsf=sigmapsf

  zp=24.84245
  scale=0.05
  extinct=0.041
  msun=4.54
  peak=sol[0,*]/(2*!PI*sol[1,*]^2*sol[2,*])
  mu=zp+5*alog10(scale)-2.5*alog10(peak)-extinct
  const=(64800/!PI)^2
  intensity=const*(10^(0.4*(msun-mu)))
  sigmaarc=sol[1,*]*scale
  forprint,intensity,sigmaarc,sol[2,*],format='F,F,F',textout='m59co_mge_output_850.dat'

  fits_read,'./galfit_output/m59co_850_galfit_model.fits',img,exten_no=1
  fits_read,'./galfit_output/m59co_850_galfit_model.fits',modelimg,exten_no=2
  find_galaxy,img,majoraxis,eps,ang,xc,yc
  find_galaxy,modelimg,majoraxis,eps,ang,xc_mod,yc_mod

  fits_read,'./make_decon/serfreemodelz.fits',serfreez
  find_galaxy,serfreez,m,e,a,xc_ser,yc_ser
  radius=(10^(findgen(77)*0.025))*scale
  minrad=[0,(10^(findgen(76)*0.025))*scale]

  radius=radius/scale
  minrad=minrad/scale
  photmag2=fltarr(n_elements(radius))
  modelmag2=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  area=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,img,xc,yc,maxflux_phot,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,img,xc,yc,minflux_phot,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,maxflux_model,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,minflux_model,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreez,xc_ser,yc_ser,maxflux_ser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfreez,xc_ser,yc_ser,minflux_ser,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     if (minrad[i] lt 0.025) then begin
        area[i]=((!PI*(radius[i])^2)) ;*scale
        
        photflux=maxflux_phot
        modelflux=maxflux_model
        sersicflux=maxflux_ser
     endif else begin
        area[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
        photflux=maxflux_phot-minflux_phot
        modelflux=maxflux_model-minflux_model
        sersicflux=maxflux_ser-minflux_ser
     endelse
     
     photmag2[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux)-extinct
     modelmag2[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modelflux)-extinct
     sersicmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(sersicflux)-extinct
  endfor
  readcol,'m59co_mge_outputsersic850_free.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_z_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xcz,ycz
  msun=4.54d
  a=where(sersicpa gt 20.)
  inintensity=sersiclum[a]
  insigma=sersicsig[a]
  inq=sersicq[a]
  inpa=sersicpa[a]
  mge2image,img,xcz,ycz,inintensity,insigma,inq,inpa,inmodel,zeropoint=zp,scale=scale,msun=msun

  b=where(sersicpa lt 20.)
  outintensity=sersiclum[b]
  outsigma=sersicsig[b]
  outq=sersicq[b]
  outpa=sersicpa[b]
  mge2image,img,xcz,ycz,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zp,scale=scale,msun=msun

  readcol,'m59co_mge_output_850.dat',mgelum,mgesig,mgeq,format='F,F,F'
  pa=fltarr(n_elements(mgeq))
  pa[*]=0.
  mge2image,img,xcz,ycz,mgelum,mgesig,mgeq,pa,mgemodel,zeropoint=zp,scale=scale,msun=msun
  writefits,'./make_decon/mgedirectfit_850.fits',mgemodel

  mgemag=fltarr(n_elements(radius))
  sersicmag1=fltarr(n_elements(radius))
  sersicmag2=fltarr(n_elements(radius))
  area2=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,inmodel,xcz,ycz,maxinflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,inmodel,xcz,ycz,mininflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcz,ycz,maxoutflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcz,ycz,minoutflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgemodel,xcz,ycz,maxmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,mgemodel,xcz,ycz,minmgeflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minrad[i] lt 0.025) then begin
        area2[i]=((!PI*(radius[i])^2)) ;*scale
        sersicflux1=maxinflux
        sersicflux2=maxoutflux
        mgeflux=maxmgeflux

     endif else begin
        area2[i]=((!PI*(radius[i])^2)-(!PI*(minrad[i])^2)) ;*scale
        sersicflux1=maxinflux-mininflux
        sersicflux2=maxoutflux-minoutflux
        mgeflux=maxmgeflux-minmgeflux
     endelse
     sersicmag1[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux1)
     sersicmag2[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(sersicflux2)
     mgemag[i]=zp+5*alog10(scale)+2.5*alog10(area2[i])-2.5*alog10(mgeflux)
     
  endfor

  set_plot,'ps'
  !P.Multi=[0,1,2]
  radius=radius*scale
  device,filename='surfbright_hstindivsersic_center850.ps',/color
  djs_plot,radius,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu_{F850LP} [Mag/sqare arcsecond]',yran=[21,14],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty
  djs_oplot,radius,sersicmag1,color='green',thick=3,linestyle=2 ;,psym=2
  djs_oplot,radius,sersicmag2,color='blue',thick=3,linestyle=2
  djs_oplot,radius,sersicmag,color='red',thick=3
  djs_oplot,radius,modelmag2,color='purple',thick=4
;  djs_oplot,radius,mgemag,color='cyan',thick=4
  items=['n=1.02','n=1.21']
  lines=[2,2]
  color=['green','blue']
  al_legend,items,linestyle=lines,colors=color,/window,background_color='white',charthick=4,thick=3,/top,/right
  djs_plot,radius,photmag2-modelmag2,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='\Delta \mu',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.5,xran=[0,1.],yran=[-0.05,0.05],/ysty
  device,/close
  !P.Multi=[0,1,1]
  set_plot,'x'
  stop


END
