pro hst_m59co_475
  fits_read,'../data/HST_9401_09_ACS_WFC_F475W_drz.fits',img,h
  fits_read,'sky_mask.fits',sky
  fits_read,'../data/HST_9401_09_ACS_WFC_F475W_drz.fits',wht,exten_no=2
  fits_read,'../data/HST_9401_09_ACS_WFC_F475W_drz.fits',mask,exten_no=3
  mdrizz=[39.2494828065,39.0788202159]
  expt=375.
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img=img+avgdrizz
  img=img[4000:4630,3700:4296]
  totcounts=img*expt
  imgsize=size(img,/dim)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if (totcounts[i,j] lt 0.) then totcounts[i,j]=mean(totcounts)
     endfor
  endfor
  err=sqrt(totcounts)
  
  writefits,'m59co_weightmap_475.fits',err
;writefits,'m59co_mask.fits',mask
  img=(img-sky)
  writefits,'m59co_475_skysubtract.fits',img
  temp=img[266:366,249:349]
  ngauss=20
  minlevel=0.
  find_galaxy,temp,majoraxis,eps,ang,xc,yc,fraction=0.8
  
  sectors_photometry,temp,eps,ang,xc,yc,radius,angle,counts,minlevel=minlevel

  readcol,'tinytim_fit_475.dat',normpsf,sigmapsf,format='F,F'
  MGE_fit_sectors,radius,angle,counts,eps,sol=sol,ngauss=ngauss,scale=scale,normpsf=normpsf,sigmapsf=sigmapsf

;  stop
  zp=26.05923
  scale=0.05
  extinct=0.107
  msun=5.11
  peak=sol[0,*]/(2*!PI*sol[1,*]^2*sol[2,*])
  mu=zp+5*alog10(scale)-2.5*alog10(peak)-extinct
  const=(64800/!PI)^2
  intensity=const*(10^(0.4*(msun-mu)))
  sigmaarc=sol[1,*]*scale
  forprint,intensity,sigmaarc,sol[2,*],format='F,F,F',textout='m59co_mge_output.dat'
  fits_read,'./galfit_output/m59co_475_galfit_modelfixed.fits',img,exten_no=1
  fits_read,'./galfit_output/m59co_475_galfit_modelfixed.fits',modelimg,exten_no=2
  find_galaxy,img,majoraxis,eps,ang,xc,yc
  find_galaxy,modelimg,majoraxis,eps,ang,xc_mod,yc_mod

  fits_read,'./make_decon/serfixmodelg.fits',serfixg
  find_galaxy,serfixg,m,e,a,xc_ser,yc_ser
  radius=findgen(140)*0.7;(10^(findgen(77)*0.025))*scale
  minrad=[0,findgen(139)*0.7];[0,(10^(findgen(76)*0.025))*scale]

;  radius=radius/scale
;  minrad=minrad/scale
  photmag2=fltarr(n_elements(radius))
  modelmag2=fltarr(n_elements(radius))
  sersicmag=fltarr(n_elements(radius))
  area=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,img,xc,yc,maxflux_phot,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,img,xc,yc,minflux_phot,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,maxflux_model,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,modelimg,xc_mod,yc_mod,minflux_model,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixg,xc_ser,yc_ser,maxflux_ser,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,serfixg,xc_ser,yc_ser,minflux_ser,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
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
     
     photmag2[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux)-0.107
     modelmag2[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modelflux)-0.107
     sersicmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(sersicflux)-0.107
  endfor
  readcol,'m59co_mge_outputsersic.dat',sersiclum,sersicsig,sersicq,sersicpa,format='D,D,D,D'
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xcg,ycg
  msun=5.11d
  a=where(sersicpa gt 20.)
  inintensity=sersiclum[a]
  insigma=sersicsig[a]
  inq=sersicq[a]
  inpa=sersicpa[a]
  mge2image,img,xcg,ycg,inintensity,insigma,inq,inpa,inmodel,zeropoint=zp,scale=scale,msun=msun

  b=where(sersicpa lt 20.)
  outintensity=sersiclum[b]
  outsigma=sersicsig[b]
  outq=sersicq[b]
  outpa=sersicpa[b]
  mge2image,img,xcg,ycg,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zp,scale=scale,msun=msun

  readcol,'m59co_mge_output.dat',mgelum,mgesig,mgeq,format='F,F,F'
  pa=fltarr(n_elements(mgeq))
  pa[*]=0.
  mge2image,img,xcg,ycg,mgelum,mgesig,mgeq,pa,mgemodel,zeropoint=zp,scale=scale,msun=msun
  writefits,'./make_decon/mgedirectfit_475.fits',mgemodel
  mgemag=fltarr(n_elements(radius))
  sersicmag1=fltarr(n_elements(radius))
  sersicmag2=fltarr(n_elements(radius))
  area2=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,inmodel,xcg,ycg,maxinflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,inmodel,xcg,ycg,mininflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcg,ycg,maxoutflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,outmodel,xcg,ycg,minoutflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,mgemodel,xcg,ycg,maxmgeflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
   aper,mgemodel,xcg,ycg,minmgeflux,fluxerr,0.,skyerr,1,minrad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

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
  device,filename='surfbright_hstindivsersic_center.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
;  myplot,filename='surfbright_hstindivsersic_center.ps'
  djs_plot,radius,photmag2,psym=2,xtitle='Radius ["]',ytitle='\mu_{F475W} [Mag/asec^2]',yran=[20.99,14.5],xran=[0,1.],charsize=1.5,charthick=4,xthick=3,ythick=3,ymargin=[0,3],xcharsize=0.000001,position=[0.1,0.25,0.95,0.95],/ysty,ycharsize=0.95
  djs_oplot,radius,sersicmag1,color='green',thick=3,linestyle=2 ;,psym=2
  djs_oplot,radius,sersicmag2,color='blue',thick=3,linestyle=2
  djs_oplot,radius,sersicmag,color='red',thick=3
  djs_oplot,radius,modelmag2,color='cyan',thick=4
;  djs_oplot,radius,mgemag,color='cyan',thick=4
  items=['n=1.02','n=1.21']
  lines=[2,2]
  color=['green','blue']
  al_legend,items,linestyle=lines,colors=color,/window,background_color='white',charthick=4,thick=3,/top,/right
  xyouts,[0.02],[20.5],['M59cO'],charthick=3,charsize=1.5,/data
  djs_plot,radius,photmag2-modelmag2,psym=2,ymargin=[4,0],position=[0.1,0.1,0.95,0.25],xtitle='Radius ["]',ytitle='\Delta \mu',charthick=4,xthick=3,ythick=2,charsize=1.5,ycharsize=0.65,xran=[0,1.],yran=[-0.04,0.04],/ysty
  djs_oplot,radius,fltarr(n_elements(radius)),linestyle=2,thick=3
  device,/close
  !P.Multi=[0,1,1]
  set_plot,'x'
  stop

END
