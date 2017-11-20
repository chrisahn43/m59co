PRO JAM_ANI
;THIS CODE RUNS JAM CODE AND WITH VARYING ANISOTROPY TO SEE AT WHAT
;POINT WE WILL HAVE A ZERO MASS BLACK HOLE.

  infile='../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,dispersion,dispmc,dispersionerr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  xbin=rav[0:n_elements(rav)-3]*0.05
  xbinerror=((rav[0:n_elements(rav)-3] - rin[0:n_elements(rav)-3])/2.)*0.05
  xbinout=rav[7:8]*0.05
  xbinouterror=((rav[7:8]-rin[7:8])/2.)*0.05
  ybin=fltarr(n_elements(rav)-2)
  disp=dispersion[0:n_elements(rav)-3]
  dispout=dispersion[7:8]
  disperr=dispersionerr[0:n_elements(rav)-3]
  disperrout=dispersionerr[7:8]
  distance=16.5
  ml=findgen(35)*0.1+0.1
  mbhs=[0.]
  betas=findgen(20)*0.1-1.
  inclinations=[60.,70.,80.,90.]
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  readcol,'../kinematics/kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'

  fixgausspsf_mass=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)

  readcol,'m59co_mge_outputsersic.dat',intensity,sigmaarc,q,format='F,F,F'
  readcol,'m59co_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
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
              fixgausspsf_mass[i,j,k,l].chi2=chi2*FLOAT(n_elements(xbin))
              fixgausspsf_mass[i,j,k,l].outmbh=mbh*fitml
              fixgausspsf_mass[i,j,k,l].ml=fitml
              fixgausspsf_mass[i,j,k,l].inmbh=mbhs[i]
              fixgausspsf_mass[i,j,k,l].inbeta=betas[j]
              fixgausspsf_mass[i,j,k,l].ininc=inclinations[k]
              fixgausspsf_mass[i,j,k,l].rms=rmsmodel
           endfor
        endfor
     endfor
  endfor
  mwrfits,fixgausspsf_mass,'./JAM_output_anisotropy/fixgausspsf_mass.fits'

  a=where(fixgausspsf_mass.chi2 eq min(fixgausspsf_mass.chi2))
  rms_anis=fixgausspsf_mass[a].rms
  b=where(fixgausspsf_mass.inbeta eq 0.)
  c=where(fixgausspsf_mass[b].chi2 eq min(fixgausspsf_mass[b].chi2))
  rms_noanis=fixgausspsf_mass[b[c]].rms
  
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.9],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/xlog
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  djs_oplot,xbin,rms_noanis,color='red',thick=4
  djs_oplot,xbin,rms_anis,color='blue',thick=4
  help,fixgausspsf_mass[a]
  help,fixgausspsf_mass[b[c]]


END
