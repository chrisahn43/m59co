PRO PLOTDISP
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

  data=mrdfits('fixgausspsf_mass.fits',1)
  set_plot,'ps'
  device,filename='disp_anisotropy.ps',/color
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  a=where(data.chi2 eq min(data.chi2))
  rms_anis=data[a].rms
  b=where(data.inbeta eq 0.)
  c=where(data[b].chi2 eq min(data[b].chi2))
  rms_noanis=data[b[c]].rms

  djs_oplot,xbin,rms_noanis,color='red',thick=4
  djs_oplot,xbin,rms_anis,color='blue',thick=4
  device,/close
  set_plot,'x'
END
