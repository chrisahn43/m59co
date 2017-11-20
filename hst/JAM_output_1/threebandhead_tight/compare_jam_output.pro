PRO COMPARE_JAM_OUTPUT
  infile='../../../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav*0.05                ;major axis*0.05 arcsec/pixel
  ybin=fltarr(n_elements(rav)) ;minor axis
  distance=16.5
  mbhs=[0.,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6]
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
  stop
  set_plot,'ps'
  device,filename='./cummulike.ps',/color
                                ;FIX MODEL WITH MASS OLD PSF
  djs_plot,mbhs,likefixmoffat_mass,ytitle='Cummulative Likelihood',xtitle='M_{BH}',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=4,color='green'
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

  djs_oplot,mbhs,likefixgauss,color='purple',thick=4

  ;stop
                                ;FIX MODEL WITH MASS GAUSS PSF
  
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

  djs_oplot,mbhs,likefixgauss_mass,thick=4

  ;stop
                                ;FREE MODEL WITH MASS WITH NEW PSF
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

  djs_oplot,mbhs,likefreemoffat_mass,color='red',thick=4

                                ;stop
;  ml=findgen(26)*0.1+0.5
;  nmls=n_elements(ml)
                                ;FREE MODEL WITH NEW PSF NO MASS
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

  djs_oplot,mbhs,likefreegauss_mass,color='blue',thick=4

  ;stop
                                ;FIX MODEL WITH NEW PSF NO MASS
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

  djs_oplot,mbhs,likefreegauss,color='cyan',thick=4
  items=['Fix Mass Gauss PSF','Free Mass Gauss PSF','Fix Mass Moffat PSF','Free Mass Moffat PSF','Fix Gauss PSF','Free Gauss PSF']
  lines=[0,0,0,0,0,0]
  color=['black','blue','green','red','purple','cyan']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
  device,filename='./sigmambhs.ps',/color
  djs_plot,sigma,mbhsfixgauss_mass,ytitle='M_{BH}',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=4,xran=[-3.5,3.5],yran=[-150000,7.e6],/xsty,/ysty
  djs_oplot,sigma,mbhsfreegauss_mass,color='blue',thick=4
  djs_oplot,sigma,mbhsfixmoffat_mass,color='green',thick=4
  djs_oplot,sigma,mbhsfreemoffat_mass,color='red',thick=4
  djs_oplot,sigma,mbhsfixgauss,color='purple',thick=4
  djs_oplot,sigma,mbhsfreegauss,color='cyan',thick=4
  sym=[0,0,0,0,0,0]
  al_legend,items,psym=sym,colors=color,background_color='white',charthick=4,thick=3,/bottom,/right
  device,/close
  

  device,filename='./onedbestfit.ps',/color
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.5],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  a=where(freegauss_mass.outmbh eq 0.)
  c=where(freegauss_mass[a].chi2 eq min(freegauss_mass[a].chi2))
  djs_oplot,xbin,freegauss_mass[a[c]].rms,color='red',thick=4
  b=where(freegauss_mass.chi2 eq min(freegauss_mass.chi2))
  djs_oplot,xbin,freegauss_mass[b].rms,color='blue',thick=4
  device,/close

  set_plot,'x'
  
  stop
END
