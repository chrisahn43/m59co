pro m59co_jam

  infile='../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'

  xbin=rav[0:7]*0.05
  ybin=fltarr(n_elements(rav)-1)
  disp=disp[0:7]
  disperr=disperr[0:7]
  readcol,'m59co_mge_outputsersic.dat',intensity,sigmaarc,q,format='F,F,F'
  readcol,'m59co_mge_outputsersic_mass.dat',mass,format='D'
                                ;readcol,'m59co_mge_outputsersic_mass.dat',mass,format='F'

  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
  sigma_pot = sigmaarc
  qobs_pot = q
  distance = 16.5              ; Assume Virgo distance in Mpc (Mei et al. 2007)

  ml=findgen(20)*0.1+0.1
;  ml=findgen(40)*0.1+0.1
  mbhs=[6.e6]
;  mbhs=[5.5e6]
  betas=[0.]
  inclinations=[11.5]

  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  out=REPLICATE({inmbh:0.0,inbeta:0.0,ininc:0.0,chi2:0.0,ml:0.0,outmbh:0.0,rms:fltarr(n_elements(disp))},nmbhs,nbetas,ninclinations,nmls)
  readcol,'./kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
  peak=normpsf/sqrt(2.*!PI*sigmapsf^2)
  gauss=fltarr(n_elements(xbin))
  for i=0,n_elements(xbin)-1 do begin
     temp=0.
     for j=0,n_elements(sigmapsf)-1 do begin
        temp+=peak[j]*exp(-0.5*(xbin[i]^2/sigmapsf[j]^2))
     endfor
     gauss[i]=temp
  endfor
  djs_plot,xbin,gauss,psym=2
  stop
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for k=0,ninclinations-1 do begin
           for l=0,nmls-1 do begin
              fitml=ml[l]
              mbh=mbhs[i]/fitml
              beta=betas[j]
              inc=inclinations[k]
              jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inc, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,step=0.005
              out[i,j,k,l].chi2=chi2*FLOAT(n_elements(xbin))
              out[i,j,k,l].outmbh=mbh*fitml
              out[i,j,k,l].ml=fitml
              out[i,j,k,l].inmbh=mbhs[i]
              out[i,j,k,l].inbeta=betas[j]
              out[i,j,k,l].ininc=inclinations[k]
              out[i,j,k,l].rms=rmsmodel
              djs_plot,xbin,disp,psym=4
              djs_oplot,xbin,rmsmodel,color='blue'
           endfor
        endfor
     endfor
  endfor
  a=where(out.outmbh eq 0.)
  c=where(out[a].chi2 eq min(out[a].chi2))
  rms_nobh=out[a[c]].rms
  b=where(out.chi2 eq min(out.chi2))

;  set_plot,'ps'
;  device,filename='oned_bestfit_rms.ps',/color
  djs_plot,xbin,disp,psym=4,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.,0.9],yran=[25,45],/xsty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  oploterr,xbin,disp,disperr
  djs_oplot,xbin,rms_nobh,color='red',thick=4
  djs_oplot,xbin,out[b].rms,color='blue',thick=4
  help,out[a[c]]
  help,out[b]
;  device,/close
;  set_plot,'x'
  stop
  
END
