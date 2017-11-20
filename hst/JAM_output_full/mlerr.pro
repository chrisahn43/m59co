PRO MLERR

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
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  ml=findgen(35)*0.1+0.1
  inclinations=[14.5,20.,30.,40.,50.,60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  ninclinations=n_elements(inclinations)
  fixgauss_mass=mrdfits('./temp/fixgausspsf_mass.fits',1)
  fixmoffat_mass=mrdfits('./temp/fixmoffatpsf_mass.fits',1)
  fixgauss=mrdfits('./temp/fixgausspsf.fits',1)
  freemoffat_mass=mrdfits('./temp/freemoffatpsf_mass.fits',1)
  freegauss_mass=mrdfits('./temp/freegausspsf_mass.fits',1)
  freegauss=mrdfits('./temp/freegausspsf.fits',1)
  ind=where(fixgauss_mass.inbeta eq 0.)
  fixgauss_mass=fixgauss_mass[ind]
  ind=where(fixmoffat_mass.inbeta eq 0.)
  fixmoffat_mass=fixmoffat_mass[ind]
  ind=where(fixgauss.inbeta eq 0.)
  fixgauss=fixgauss[ind]
  ind=where(freemoffat_mass.inbeta eq 0.)
  freemoffat_mass=freemoffat_mass[ind]
  ind=where(freegauss_mass.inbeta eq 0.)
  freegauss_mass=freegauss_mass[ind]
  ind=where(freegauss.inbeta eq 0.)
  freegauss=freegauss[ind]
  mlpop=4.2872776
  likelihood=dblarr(nmls)
  temp=0.D
  for i=0,nmls-1 do begin
     mlind=where((fixmoffat_mass.ml) eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixmoffat_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixmoffat_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixmoffat_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixmoffat_mass=likelihood/max(likelihood)
  mlfixmoffat_mass=(interpol(ml,likefixmoffat_mass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop

  sigma=[-3,-2,-1,0,1,2,3]
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixgauss.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixgauss[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss=likelihood/max(likelihood)
  mlfixgauss=interpol(ml,likefixgauss,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(fixgauss_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(fixgauss_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass=likelihood/max(likelihood)
  mlfixgauss_mass=(interpol(ml,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(freemoffat_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(freemoffat_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(freemoffat_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freemoffat_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreemoffat_mass=likelihood/max(likelihood)
  mlfreemoffat_mass=(interpol(ml,likefreemoffat_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(freegauss_mass.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(freegauss_mass[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss_mass[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss_mass[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss_mass=likelihood/max(likelihood)
  mlfreegauss_mass=(interpol(ml,likefreegauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  temp=0D
  for i=0,nmls-1 do begin
     mlind=where(freegauss.ml eq ml[i])
     for j=0,nmbhs-1 do begin
        mbhind=where(freegauss[mlind].outmbh eq mbhs[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss[mlind[mbhind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss[mlind[mbhind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss=likelihood/max(likelihood)
  mlfreegauss=interpol(ml,likefreegauss,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  set_plot,'ps'
  device,filename='./cummulike_ml.ps',/color
  
  djs_plot,ml*mlpop,likefixgauss_mass,ytitle='Cummulative Likelihood',xtitle='M/L',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[0.,4.],xstyle=8,/nodata
  mlratio=(ml/mlpop)
  xt='M/L_{dyn} / M/L_{pop}'
  xt=textoidl(xt)

  axis,xaxis=1,xtickv=mlratio,xcharsize=1.5,charthick=4,xthick=4,xtitle=xt,/xsty,xran=[0.0,0.932993]

  arrow,1.44,0.5,1.44,1.,/data,thick=8,color=100,hsize=0.01
  arrow,1.44,0.5,1.44,0.,/data,thick=8.,color=100,hsize=0.01
  arrow,0.9677,0.16,0.9677,1.,/data,thick=4,color=100,hsize=0.01
  arrow,0.9677,0.16,0.9677,0.,/data,thick=4.,color=100,hsize=0.01
  arrow,1.846,1.,1.846,0.,/data,thick=4.,color=100,hsize=0.01;0.84
  arrow,0.431,1.,0.431,0.,/data,thick=2.,color=100,hsize=0.01;0.00135
  arrow,2.704,1.,2.704,0.,/data,thick=2.,color=100,hsize=0.01 ;0.99865
  djs_oplot,ml*mlpop,likefixgauss_mass,thick=12
  djs_oplot,ml*mlpop,likefixgauss_mass,thick=12
  djs_oplot,ml*mlpop,likefixmoffat_mass,thick=4,color='red'
  djs_oplot,ml*mlpop,likefreegauss_mass,thick=4,color='blue'
  djs_oplot,ml,likefixgauss,thick=4,color='blue',linestyle=1
  djs_oplot,ml,likefreegauss,thick=4,color='blue',linestyle=2
;  mlratio=(ml/mlpop)
;  xt='M/L_{dyn} / M/L_{pop}'
;  xt=textoidl(xt)
;  axis,xaxis=1,xtickv=mlratio,xcharsize=1.5,charthick=4,xthick=4,xtitle=xt,/xsty,xran=[0.0,0.932993]
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/bottom,/right
  device,/close
  device,filename='./sigmambhs_ml.ps',/color
  djs_plot,sigma,mlfixgauss_mass,ytitle='M/L',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[-3.5,3.5],yran=[0.1,3.],/xsty,/ysty
  djs_oplot,sigma,mlfixmoffat_mass,color='red',thick=4
  djs_oplot,sigma,mlfreegauss_mass,color='blue',thick=4
  djs_oplot,sigma,mlfixgauss,color='blue',thick=4,linestyle=1
  djs_oplot,sigma,mlfreegauss,color='blue',thick=4,linestyle=2
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close



  stop

END
