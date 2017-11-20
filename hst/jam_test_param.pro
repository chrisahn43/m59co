PRO JAM_TEST_PARAM
  infile='../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,dispersion,dispmc,dispersionerr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav[0:n_elements(rav)-3]*0.05
  ybin=fltarr(n_elements(rav)-2)
  disp=dispersion[0:n_elements(rav)-3]
  disperr=dispersionerr[0:n_elements(rav)-3]
  distance=16.5
;  ml=findgen(40)*0.1+0.1
  ml=findgen(35)*0.1+0.1
;  mbhs=[0.,10^(findgen(6)*0.2+6.)]
  mbhs=[5.8e6];,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  betas=[0.];findgen(11)*0.1-0.2;[0.]
  inclinations=[14.5];,20.,30.,40.,50.,60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
                                ;START WITH DOUBLE GAUSSIAN PSF
  readcol,'../kinematics/kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
                                ;LUM ONLY FREE
  fixgausspsf=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)

  readcol,'./make_decon/m59co_mge_outputsersic_k.dat',intensity,sigmaarc,q,format='F,F,F'
  readcol,'m59co_mge_outputsersic_nmass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass/10.
  sigma_pot = sigmaarc
  qobs_pot = q
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for k=0,ninclinations-1 do begin
           for l=0,nmls-1 do begin
              fitml=ml[l]
              mbh=mbhs[i]/fitml
              beta=betas[j]
              inc=inclinations[k]
              jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
              fixgausspsf[i,j,k,l].chi2=chi2*FLOAT(n_elements(xbin))
              fixgausspsf[i,j,k,l].outmbh=mbh*fitml
              fixgausspsf[i,j,k,l].ml=fitml
              fixgausspsf[i,j,k,l].inmbh=mbhs[i]
              fixgausspsf[i,j,k,l].inbeta=betas[j]
              fixgausspsf[i,j,k,l].ininc=inclinations[k]
              fixgausspsf[i,j,k,l].rms=rmsmodel
           endfor
        endfor
     endfor
  endfor
  stop
  mwrfits,fixgausspsf,'./JAM_output_fixedml/fixgausspsf_mass.fits'
  a=where(fixgausspsf.outmbh eq 0.)
  c=where(fixgausspsf[a].chi2 eq min(fixgausspsf[a].chi2))
  rms_nobh_fixgauss=fixgausspsf[a[c]].rms
  b=where(fixgausspsf.chi2 eq min(fixgausspsf.chi2))
  rms_fixgauss=fixgausspsf[b].rms
  djs_plot,rav*0.05,dispersion,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.9],yran=[25,45],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,rav*0.05,dispersion,dispersionerr
  djs_oplot,rav*0.05,rms_nobh_fixgauss,color='red',thick=4
  djs_oplot,rav*0.05,rms_fixgauss,color='blue',thick=4
  help,fixgausspsf[a[c]]
  help,fixgausspsf[b]
  stop

  

END
