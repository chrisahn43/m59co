PRO M59CO_JAMVEL
  infile='../kinematics/vor_out/intspec_wallace_goodmay.dat'
  READCOL,infile,filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'

  xbin=rav*0.05
  ybin=fltarr(n_elements(rav))
;  vel=vel-719.840
  nvel=fltarr(n_elements(rav))
  nvelerr=fltarr(n_elements(rav))
  for i=0,n_elements(rav)-1,2 do begin
     xbin[i]=-xbin[i]
     nvel[i]=(vel[i+1]-vel[i])/2.
     nvelerr[i]=sqrt(velerr[i+1]^2+velerr[i]^2)
  endfor
  sorted=sort(xbin)
  xbin=xbin[sorted]
  vel=vel[sorted]
  velerr=velerr[sorted]
  ind=where(xbin gt 0.)
  xbin=xbin[ind]
  ind=where(nvel gt 0.)
  vel=nvel[ind]
  velerr=nvelerr[ind]
  ybin=fltarr(n_elements(xbin))
  stop
  readcol,'m59co_mge_outputsersic.dat',intensity,sigmaarc,q,format='F,F,F'
  readcol,'m59co_mge_outputsersic_mass.dat',mass,format='D'

  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass*0.45;intensity*1.3;
  sigma_pot = sigmaarc
  qobs_pot = q
  distance = 16.5              ; Assume Virgo distance in Mpc (Mei et al. 2007)
  inclination=11.5
  mbhs=[5.9e6]
  betas=[-0.2,0.,0.4,0.8];findgen(11)*0.1-0.2
  ml=findgen(35)*0.1+0.1
  kappas=findgen(5)*0.5
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  nkappa=n_elements(kappas)
  out=REPLICATE({inmbh:0.0,inbeta:0.0,chi2:0.0,outmbh:0.0,inkappa:0.0,vel:fltarr(n_elements(vel))},nmbhs,nbetas,nkappa)
  readcol,'./kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        for l=0,nkappa-1 do begin
           mbh=mbhs[i]          ;/fitml
           beta=betas[j]
           kappa=kappas[l]
           jam_axisymmetric_vel,surf_lum,sigma_lum,qobs_lum,surf_pot,sigma_pot,qobs_pot,inclination,mbh,distance,xbin,ybin,velmodel,beta=beta,vel=vel,evel=velerr,normpsf=normpsf,sigmapsf=sigmapsf,pixsize=0.05,chi2=chi2,step=0.005,kappa=kappa
           out[i,j,l].chi2=chi2*FLOAT(N_ELEMENTS(XBIN))
           out[i,j,l].outmbh=mbhs[i] ;*fitml
           out[i,j,l].inmbh=mbhs[i]
           out[i,j,l].inbeta=betas[j]
           out[i,j,l].inkappa=kappa;s[l]
           out[i,j,l].vel=velmodel
           if (chi2 lt 10.) then begin
              djs_plot,xbin,vel,psym=4,yran=[min(vel)-5,max(vel)+5]
              djs_oplot,xbin,velmodel,color='blue'
;              stop
           endif
           
        endfor
           
     endfor
  endfor
  set_plot,'ps'
  device,filename='jamvel.ps',/color
;  a=where(out.chi2 eq min(out.chi2))
  djs_plot,xbin,vel,psym=4,yran=[min(vel)-5,max(vel)+5],ytitle='Velocity [km/s]',xtitle='Radius ["]',charsize=1.5,charthick=4,xthick=3,ythick=3
  oploterr,xbin,vel,velerr
  color=['blue','red','purple','cyan']
  lines=[0,0,0,0]
  items=['$\beta$ = -0.2, $\kappa$ = 0.59','$\beta$ = 0.0, $\kappa$ = 0.63','$\beta$ = 0.4, $\kappa$ = 0.89','$\beta$ = 0.8, $\kappa$ = 2.23']
  for i=0,nbetas-1 do begin &$
     ind=where(out.inbeta eq betas[i]) &$
     cind=where(out[ind].chi2 eq min(out[ind].chi2)) &$
     djs_oplot,xbin,out[ind[cind]].vel,color=color[i],thick=4 &$
     help,out[ind[cind]] &$
  endfor
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  device,/close
  set_plot,'x'
  stop
END
