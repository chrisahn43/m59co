PRO COMPARE_JAM_OUTPUT
  infile='../../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav[0:n_elements(rav)-3]*0.05 ;major axis*0.05 arcsec/pixel
  xbinerror=((rav[0:n_elements(rav)-3] - rin[0:n_elements(rav)-3])/2.)*0.05
  xbinout=rav[7:8]*0.05
  xbinouterror=((rav[7:8]-rin[7:8])/2.)*0.05
  dispout=disp[7:8]
  disp=disp[0:n_elements(rav)-3]
  disperrout=disperr[7:8]
  disperr=disperr[0:n_elements(rav)-3]
  ybin=fltarr(n_elements(xbin)) ;minor axis
  distance=16.5
  mbhs=[0.,1.e5,5.e5,1.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  ml=findgen(35)*0.1+0.1
  inclinations=[60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  fixgauss_mass=mrdfits('fixgausspsf_mass.fits',1)
  fixmoffat_mass=mrdfits('fixmoffatpsf_mass.fits',1)
  fixgauss=mrdfits('fixgausspsf.fits',1)
  freemoffat_mass=mrdfits('freemoffatpsf_mass.fits',1)
  freegauss_mass=mrdfits('freegausspsf_mass.fits',1)
  freegauss=mrdfits('freegausspsf.fits',1)

                                ;COMPUTE LIKELIHOOD FUNCTION FOR EACH
                                ;START WITH FIX MODEL WITH MASS AND NEW PSF
  likelihood=dblarr(nmbhs)
  temp=0.D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(fixmoffat_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixmoffat_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixmoffat_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixmoffat_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixmoffat_mass=likelihood/max(likelihood)
  mbhsfixmoffat_mass=interpol(mbhs,likefixmoffat_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  sigma=[-3,-2,-1,0,1,2,3]
;  stop
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(fixgauss.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss=likelihood/max(likelihood)
  mbhsfixgauss=interpol(mbhs,likefixgauss,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])


  ;stop
  
  
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(fixgauss_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass=likelihood/max(likelihood)
  mbhsfixgauss_mass=interpol(mbhs,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])


  ;stop
 
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(freemoffat_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freemoffat_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freemoffat_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freemoffat_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreemoffat_mass=likelihood/max(likelihood)
  mbhsfreemoffat_mass=interpol(mbhs,likefreemoffat_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])


                                ;stop

  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(freegauss_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freegauss_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss_mass=likelihood/max(likelihood)
  mbhsfreegauss_mass=interpol(mbhs,likefreegauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])


  ;stop
 
  temp=0D
  for i=0,nmbhs-1 do begin
 ;    temp=0D
     mbhind=where(freegauss.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freegauss[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss=likelihood/max(likelihood)
  mbhsfreegauss=interpol(mbhs,likefreegauss,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  set_plot,'ps'
  device,filename='./cummulike.ps',/color
  djs_plot,mbhs,likefixgauss_mass,ytitle='Cummulative Likelihood',xtitle='M_{BH}',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12
  djs_oplot,mbhs,likefixmoffat_mass,thick=4,color='red'
  djs_oplot,mbhs,likefreegauss_mass,thick=4,color='blue'
  djs_oplot,mbhs,likefixgauss,thick=4,color='blue',linestyle=1
  djs_oplot,mbhs,likefreegauss,thick=4,color='blue',linestyle=2

  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
 ;  djs_plot,mbhs,likefixmoffat_mass,ytitle='Cummulative Likelihood',xtitle='M_{BH}',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=4,color='green'
;  djs_oplot,mbhs,likefixgauss,color='purple',thick=4
;  djs_oplot,mbhs,likefixgauss_mass,thick=12
;  djs_oplot,mbhs,likefreemoffat_mass,color='red',thick=4
;  djs_oplot,mbhs,likefreegauss_mass,color='blue',thick=4
;  djs_oplot,mbhs,likefreegauss,color='cyan',thick=4
;  items=['Fix Mass Gauss PSF','Free Mass Gauss PSF','Fix Mass Moffat PSF','Free Mass Moffat PSF','Fix Gauss PSF','Free Gauss PSF']
;  lines=[0,0,0,0,0,0]
;  color=['black','blue','green','red','purple','cyan']
;  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
;  device,/close
  device,filename='./sigmambhs.ps',/color
  djs_plot,sigma,mbhsfixgauss_mass,ytitle='M_{BH}',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[-3.5,3.5],yran=[2.e6,9.e6],/xsty,/ysty
  djs_oplot,sigma,mbhsfixmoffat_mass,color='red',thick=4
  djs_oplot,sigma,mbhsfreegauss_mass,color='blue',thick=4
  djs_oplot,sigma,mbhsfixgauss,color='blue',thick=4,linestyle=1
  djs_oplot,sigma,mbhsfreegauss,color='blue',thick=4,linestyle=2
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close


;  djs_oplot,sigma,mbhsfixmoffat_mass,color='green',thick=4
;  djs_oplot,sigma,mbhsfreemoffat_mass,color='red',thick=4
;  djs_oplot,sigma,mbhsfixgauss,color='purple',thick=4
;  djs_oplot,sigma,mbhsfreegauss,color='cyan',thick=4
;  sym=[0,0,0,0,0,0]
;  al_legend,items,psym=sym,colors=color,background_color='white',charthick=4,thick=3,/bottom,/right
;  device,/close
  
;  set_plot,'x'
;  device,filename='./onedbestfit.ps',/color
;  plotsym,0,1/2.,/fill
;  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.9],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/xlog
;  plotsym,0,1/2.,/fill,color=100
;  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
;  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
;  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  
;  a=where(freegauss_mass.outmbh eq 0.)
;  c=where(freegauss_mass[a].chi2 eq min(freegauss_mass[a].chi2))
;  djs_oplot,xbin,freegauss_mass[a[c]].rms,color='red',thick=4
;  b=where(freegauss_mass.chi2 eq min(freegauss_mass.chi2))
;  djs_oplot,xbin,freegauss_mass[b].rms,color='blue',thick=4
;  device,/close
  
;  set_plot,'x'
  
;  stop
  stop
  mbhs=[0.,6318284.4]
  mbhs=[0.,5861087.5];fix model no mass
  inclinations=90.
  beta=0.
  nmbhs=n_elements(mbhs)
  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
  readcol,'~/research/code/gemini15/m59co/kinematics/kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../m59co_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../m59co_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
  surf_pot = intensity ;no mass
  sigma_pot = sigmaarc
  qobs_pot = q
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        fitml=ml[j]
        mbh=mbhs[i]/fitml
        jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inclinations, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
        final[i,j].chi2=chi2*FLOAT(n_elements(xbin))
        final[i,j].inmbh=mbh*fitml
        final[i,j].ml=fitml
        final[i,j].rms=rmsmodel
     endfor
  endfor
;  stop
  scale=0.05
  zeropoint=26.05923
  msun=5.11
  mlgin=2.50556d
  mlgout=5.45d
  rad=50.;6.337

  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xc,yc
  inind=where(pa gt 20.)
  inintensity=intensity[inind]
  insigma=sigmaarc[inind]
  inq=q[inind]
  inpa=pa[inind]
  mge2image,img,xc,yc,inintensity,insigma,inq,inpa,inmodel,zeropoint=zeropoint,scale=scale,msun=msun

  outind=where(pa lt 20.)
  outintensity=intensity[outind]
  outsigma=sigmaarc[outind]
  outq=q[outind]
  outpa=pa[outind]
  mge2image,img,xc,yc,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zeropoint,scale=scale,msun=msun

  aper,inmodel,xc,yc,influx,influxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,outmodel,xc,yc,outflux,outfluxerr,0.,skyerr,1,rad,-1,[1,1],/silent,setskyval=0.,/flux,/exact

  bhmlind=where(final.chi2 eq min(final.chi2))
  bhml=final[bhmlind].ml
  inmag=zeropoint-2.5*alog10(influx)
  inabsmag=inmag-5*(alog10(16.5e6)-1)
  inlum=10^((inabsmag-msun)/(-2.5))
  outmag=zeropoint-2.5*alog10(outflux)
  outabsmag=outmag-5*(alog10(16.5e6)-1)
  outlum=10^((outabsmag-msun)/(-2.5))
  bhmass=bhml*((mlgin*inlum)+(mlgout*outlum))
  temp=where(final.inmbh eq 0.)
  nobhmlind=where(final[temp].chi2 eq min(final[temp].chi2))
  nobhml=final[temp[nobhmlind]].ml
  nobhmass=nobhml*((mlgin*inlum)+(mlgout*outlum))
  aveml=((influx*mlgin)+(outflux*mlgout))/(influx+outflux)
  
  print,'Total mass with black hole = ', final[bhmlind].inmbh, ' = ',bhmass
  print,'Dynamical M/L with black hole = ', bhml;*aveml
  print,'Total mass w/o black hole = ', nobhmass
  print,'Dynamical M/L w/o black hole = ',nobhml;*aveml
  
  intdisp=32900.
  r=((16.5e6)*(3.086e16))*((rad*0.05)*(4.848e-6))
  g=6.67e-11
  dynmass=(((intdisp^2)*r)/(g))/(1.989e30)
  print,'Total mass from integrated dispersion = ', dynmass
  stop
;  set_plot,'ps'
  device,filename='./onedbestfit_BEST.ps',/color
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/xlog
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  a=where(final.inmbh eq 0.)
  c=where(final[a].chi2 eq min(final[a].chi2))
  djs_oplot,xbin,final[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  device,/close
  device,filename='./onedbestfit_BEST_linear.ps',/color
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  a=where(final.inmbh eq 0.)
  c=where(final[a].chi2 eq min(final[a].chi2))
  djs_oplot,xbin,final[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  device,/close

  set_plot,'x'

  stop
END
