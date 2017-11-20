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
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  ml=findgen(35)*0.1+0.1
  inclinations=[14.5,20.,30.,40.,50.,60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  nbeta=n_elements(betas)
  fixgauss_mass=mrdfits('./temp/fixgausspsf_mass.fits',1)
  fixmoffat_mass=mrdfits('./temp/fixmoffatpsf_mass.fits',1)
  fixgauss=mrdfits('./temp/fixgausspsf.fits',1)
  freemoffat_mass=mrdfits('./temp/freemoffatpsf_mass.fits',1)
  freegauss_mass=mrdfits('./temp/freegausspsf_mass.fits',1)
  freegauss=mrdfits('./temp/freegausspsf.fits',1)
  fixoldgauss_mass=mrdfits('./temp/fixoldgausspsf_mass.fits',1)
  

  minchi=dblarr(nmbhs,nbeta)
  data=reform(fixgauss_mass,nmbhs,nbetas,ninclinations,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        minchi[i,j]=MIN(data[i,j,*,*].chi2)
     endfor
  endfor
  
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betambh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass.outmbh,fixgauss_mass.inbeta,YTITLE='!7b!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(betas)-0.1,max(betas)+0.1],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] ;,/XLOG
;  minlev=MIN(minchi) & maxlev=MAX(minchi) & nlev=2d5
;  LOADCT,0,/SILENT
;  chi2im=(minchi-minlev)/(maxlev-minlev)
;  chi2im=chi2im/TOTAL(chi2im)*nlev
;  imgunder,chi2im
  
;  PLOTSYM,0,/FILL
  OPLOT,fixgauss_mass.outmbh,fixgauss_mass.inbeta,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  xyouts,[1.e5],[-0.25],['M59cO'],charthick=3,charsize=1.5,/data
  min=MIN(fixgauss_mass.chi2,minpos)
  CONTOUR,SMOOTH(minchi,[2,4]),mbhs,betas,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+2.3,fixgauss_mass[minpos].chi2+6.18,fixgauss_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]
;  CONTOUR,SMOOTH(minchi,4),mbhs,betas,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+3.2],C_LINESTYLE=[0],C_THICK=[10]

  xx=REPLICATE(fixgauss_mass[minpos].outmbh,N_ELEMENTS(fixgauss_mass.outmbh))
  yy=REPLICATE(fixgauss_mass[minpos].inbeta,N_ELEMENTS(fixgauss_mass.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  set_plot,'x'
;  stop

  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[i,*,*,j].chi2)
     endfor
  endfor

  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='mlmbh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass.outmbh,fixgauss_mass.ml,YTITLE='!7C!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,1.2],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
;  PLOTSYM,0,/FILL
  OPLOT,fixgauss_mass.outmbh,fixgauss_mass.ml,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  xyouts,[1.e5],[0.05],['M59cO'],charthick=3,charsize=1.5,/data
  min=MIN(fixgauss_mass.chi2,minpos)
  CONTOUR,SMOOTH(minchi,[1,1],/EDGE_TRUNCATE,/NAN),mbhs,ml,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+2.3,fixgauss_mass[minpos].chi2+6.18,fixgauss_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass[minpos].outmbh,N_ELEMENTS(fixgauss_mass.outmbh))
  yy=REPLICATE(fixgauss_mass[minpos].ml,N_ELEMENTS(fixgauss_mass.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  stop

  minchi=dblarr(nbeta,nmls)

  for i=0,nbetas-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[*,i,*,j].chi2)
     endfor
  endfor

  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betaml.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass.inbeta,fixgauss_mass.ml,YTITLE='M/L',XTITLE='Beta',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,1.2],CHARSIZE=1,XSTY = 1,XRANGE=[min(betas)-0.1,max(betas)+0.1],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
;  PLOTSYM,0,/FILL
  OPLOT,fixgauss_mass.inbeta,fixgauss_mass.ml,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass.chi2,minpos)
  CONTOUR,minchi,betas,ml,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+2.3,fixgauss_mass[minpos].chi2+6.18,fixgauss_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass[minpos].inbeta,N_ELEMENTS(fixgauss_mass.inbeta))
  yy=REPLICATE(fixgauss_mass[minpos].ml,N_ELEMENTS(fixgauss_mass.inbeta))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150

  stop
  ind=where(fixgauss_mass.inbeta eq 0.)
  fixgauss_mass=fixgauss_mass[ind]
  minchi=dblarr(nmbhs,ninclinations)
  for i=0,nmbhs-1 do begin
     for j=0,ninclinations-1 do begin
        mbh=mbhs[i] & inc=inclinations[j]
        ind=where(fixgauss_mass.ininc eq inc and fixgauss_mass.outmbh eq mbh)
        minchi[i,j]=MIN(fixgauss_mass[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incmbh_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass.outmbh,fixgauss_mass.ininc,YTITLE='Inclination',XTITLE='MBH',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
;  PLOTSYM,0,/FILL
  OPLOT,fixgauss_mass.outmbh,fixgauss_mass.ininc,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass.chi2,minpos)
  CONTOUR,minchi,mbhs,inclinations,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+2.3,fixgauss_mass[minpos].chi2+6.18,fixgauss_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass[minpos].outmbh,N_ELEMENTS(fixgauss_mass.outmbh))
  yy=REPLICATE(fixgauss_mass[minpos].ininc,N_ELEMENTS(fixgauss_mass.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150

  minchi=dblarr(nmls,ninclinations)
  for i=0,nmls-1 do begin
     for j=0,ninclinations-1 do begin
        mls=ml[i] & inc=inclinations[j]
        ind=where(fixgauss_mass.ininc eq inc and fixgauss_mass.ml eq mls)
        minchi[i,j]=MIN(fixgauss_mass[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incml_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass.ml,fixgauss_mass.ininc,YTITLE='Inclination',XTITLE='M/L',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(ml)-0.1,1.2],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

  
;  PLOTSYM,0,/FILL
  OPLOT,fixgauss_mass.ml,fixgauss_mass.ininc,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass.chi2,minpos)
  CONTOUR,minchi,ml,inclinations,/OVERPLOT,LEVELS=[fixgauss_mass[minpos].chi2+2.3,fixgauss_mass[minpos].chi2+6.18,fixgauss_mass[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass[minpos].ml,N_ELEMENTS(fixgauss_mass.ml))
  yy=REPLICATE(fixgauss_mass[minpos].ininc,N_ELEMENTS(fixgauss_mass.ml))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150


  set_plot,'x'
  stop
  ;DOING MARGINALIZATION ASSUMING ISOTROPY
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
  ind=where(fixoldgauss_mass.inbeta eq 0.)
  fixoldgauss_mass=fixoldgauss_mass[ind]
  likelihood=dblarr(nmbhs)
  temp=0.D

  for i=0,nmbhs-1 do begin
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
  sigma=[-3,-2,-1,0,1,2,3]
  temp=0D
  for i=0,nmbhs-1 do begin
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
  temp=0D
  for i=0,nmbhs-1 do begin
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
  temp=0D
  for i=0,nmbhs-1 do begin
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
  temp=0D
  for i=0,nmbhs-1 do begin
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
  temp=0D
  for i=0,nmbhs-1 do begin
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
  temp=0.D

  for i=0,nmbhs-1 do begin
     mbhind=where(fixoldgauss_mass.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixoldgauss_mass[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixoldgauss_mass[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixoldgauss_mass[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixoldgauss_mass=likelihood/max(likelihood)
  mbhsfixoldgauss_mass=interpol(mbhs,likefixoldgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])

;  set_plot,'ps' & LOADCT,0,/SILENT
;  device,filename='./cummulike.ps',/color
;  device,filename='./sigmambhs.ps',/color
;  djs_plot,sigma,mbhsfixgauss_mass,ytitle='M_{BH}',xtitle='\sigma',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,xran=[-3.5,3.5],yran=[2.e6,9.e6],/xsty,/ysty
;  djs_oplot,sigma,mbhsfixmoffat_mass,color='red',thick=4
;  djs_oplot,sigma,mbhsfreegauss_mass,color='blue',thick=4
;  djs_oplot,sigma,mbhsfixgauss,color='blue',thick=4,linestyle=1
;  djs_oplot,sigma,mbhsfreegauss,color='blue',thick=4,linestyle=2
;  items=['Best Fit Model','PSF Variations','Model Variations']
;  lines=[0,0,0]
;  color=['black','red','blue']
;  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
;  device,/close
;  set_plot,'x'
  

  stop

  mbhs=[findgen(17)*1.e5+5.2e6]
  inclinations=14.5
  beta=0.
  ml=findgen(35)*0.1+0.1
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
  readcol,'~/research/code/gemini15/m59co/kinematics/kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../m59co_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../m59co_mge_outputsersic_mass.dat',mass,format='D'
  surf_lum = intensity
  sigma_lum = sigmaarc
  qobs_lum = q
  surf_pot = mass
;  surf_pot = intensity ;no mass
  sigma_pot = sigmaarc
  qobs_pot = q
;  for i=0,nmbhs-1 do begin
;     for j=0,nmls-1 do begin
;        fitml=ml[j]
;        mbh=mbhs[i]/fitml
;        jam_axisymmetric_rms,surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, inclinations, mbh, distance, xbin, ybin, rmsModel, BETA=beta,RMS=disp,ERMS=disperr,SIGMAPSF=sigmapsf,NORMPSF=normpsf,PIXSIZE=0.05,ml=fitml,chi2=chi2,STEP=0.02
;        final[i,j].chi2=chi2*FLOAT(n_elements(xbin))
;        final[i,j].inmbh=mbh*fitml
;        final[i,j].ml=fitml
;        final[i,j].rms=rmsmodel
;     endfor
;  endfor
;  stop
  final=mrdfits('finalgrid.fits',1)
  final=reform(final,17,35)
  stop
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6,findgen(3)*1.e5+5.2e6,findgen(4)*1.e5+5.6e6,findgen(4)*1.e5+6.1e6,findgen(3)*1.e5+6.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)
  ml=[findgen(35)*0.1+0.1];,findgen(35)*0.1+0.05]
  nmls=n_elements(ml)
  ind=where(fixgauss_mass.ininc eq 14.5)
  fixgauss_mass=fixgauss_mass[ind]
  fixgauss_mass=reform(fixgauss_mass,19,35)
  finalbh=[fixgauss_mass.outmbh,final.inmbh]
  sorted=sort(finalbh)
  finalbh=finalbh[sorted]
  finalml=[fixgauss_mass.ml,final.ml]
  finalml=finalml[sorted]
  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        if (mbhs[i] lt 5.2e6 OR mbhs[i] gt 6.8e6) then begin
           mbhind=where(fixgauss_mass.inmbh eq mbhs[i])
           mlind=where(fixgauss_mass[mbhind].ml eq ml[j])
           minchi[i,j]=fixgauss_mass[mbhind[mlind]].chi2
        endif else begin
           mbhind=where(final.inmbh gt mbhs[i]-1 and final.inmbh lt mbhs[i]+1)
           mlind=where((final[mbhind].ml gt ml[j]-0.09 and final[mbhind].ml lt ml[j]-0.09) OR final[mbhind].ml eq ml[j])
           minchi[i,j]=final[mbhind[mlind]].chi2
        endelse
     endfor
  endfor
            
  finalchi2=[fixgauss_mass.chi2,final.chi2]
  finalchi2=finalchi2[sorted]
  set_plot,'ps' & LOADCT,40,/SILENT
;  device,filename='finalmlmbh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7
;  plot,finalbh,finalml,ytitle='M/L!Idyn!N / M/L!Ipop',xtitle='M!IBH',XCHARS=1.2,/NODATA,YCHARS=1.5,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(finalml)-0.01,1.0],CHARSIZE=1.5,XSTY = 1,XRANGE=[min(finalbh)-2d4,MAX(finalbh)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 

;  plotsym,0,/fill
;  OPLOT,finalbh,finalml,psym=8,thick=5,symsize=1,color=100
;  min=MIN(minchi,minpos)
;  contour,smooth(minchi,[1,1]),mbhs,ml,/OVERPLOT,LEVELS=[minchi[minpos]+2.3,minchi[minpos]+6.18,minchi[minpos]+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]
;  device,/close
;  stop
  likelihood=dblarr(nmbhs)
  
  temp=0D
  for i=0,nmbhs-1 do begin
     if (mbhs[i] lt 5.2e6 OR mbhs[i] gt 6.8e6) then begin
        mbhind=where(fixgauss_mass.inmbh eq mbhs[i])
        for j=0,nmls-1 do begin
           mlind=where(fixgauss_mass[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(fixgauss_mass[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endif else begin
        mbhind=where(final.inmbh gt mbhs[i]-1 and final.inmbh lt mbhs[i]+1)
        for j=0,nmls-1 do begin
           mlind=where((final[mbhind].ml gt ml[j]-0.09 and final[mbhind].ml lt ml[j]-0.09) OR final[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(final[mbhind[mlind]].chi2))
        endfor
        likelihood[i]=temp
     endelse
  endfor
  minchi=dblarr(nmbhs)
  temp=0D

  for i=0,nmbhs-1 do begin
     if (mbhs[i] lt 5.2e6 OR mbhs[i] gt 6.8e6) then begin
        ind=where(fixgauss_mass.outmbh eq mbhs[i])
        temp+=exp(-0.5*(min(fixgauss_mass[ind].chi2)))
     endif else begin
        ind=where(final.inmbh eq mbhs[i])
        temp+=exp(-0.5*(min(final[ind].chi2)))
     endelse
     minchi[i]=temp
  endfor
  likeminchi=minchi/max(minchi)
  testmbh=interpol(mbhs,likeminchi,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  likefixgauss_mass=likelihood/max(likelihood)
  mbhsfixgauss_mass=interpol(mbhs,likefixgauss_mass,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  print,'Best fit= ',mbhsfixgauss_mass[3]
  print,'1-sigma errors = ',mbhsfixgauss_mass[4]-mbhsfixgauss_mass[3],mbhsfixgauss_mass[3]-mbhsfixgauss_mass[2]
  print,'3-sigma errors = ',mbhsfixgauss_mass[6]-mbhsfixgauss_mass[3],mbhsfixgauss_mass[3]-mbhsfixgauss_mass[0]
;  myplot,filename='./cummulike.ps'
  set_plot,'ps'
  device,filename='./cummulike.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  LOADCT,0,/SILENT
  djs_plot,mbhs,likefixgauss_mass,ytitle='Cumulative Likelihood',xtitle='M_{BH} [M_{\odot}]',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,/nodata
  arrow,mbhsfixgauss_mass[3],0.,mbhsfixgauss_mass[3],1.,/data,thick=8.,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass[2],0.,mbhsfixgauss_mass[2],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass[4],0.,mbhsfixgauss_mass[4],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass[0],0.,mbhsfixgauss_mass[0],1.,/data,thick=2,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass[6],0.,mbhsfixgauss_mass[6],1.,/data,thick=2,color=100,hsize=0.01

  djs_oplot,mbhs,likefixgauss_mass,thick=12
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]

  djs_oplot,mbhs,likefixmoffat_mass,thick=4,color='red'
  djs_oplot,mbhs,likefixoldgauss_mass,thick=4,color='red',linestyle=1
  djs_oplot,mbhs,likefreegauss_mass,thick=4,color='blue'
  djs_oplot,mbhs,likefixgauss,thick=4,color='blue',linestyle=1
  djs_oplot,mbhs,likefreegauss,thick=4,color='blue',linestyle=2
  
  items=['Best Fit Model','PSF Variations','Model Variations']
  lines=[0,0,0]
  color=['black','red','blue']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  xyouts,[8.5e6],[0.05],['M59cO'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6,findgen(3)*1.e5+5.2e6,findgen(4)*1.e5+5.6e6,findgen(4)*1.e5+6.1e6,findgen(3)*1.e5+6.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)

  mlpop=4.2872776
  stop
  likelihood=dblarr(nmls)
  temp=0D
  for i=0,nmls-1 do begin
     for j=0,nmbhs-1 do begin
        if (mbhs[j] lt 5.2e6 OR mbhs[j] gt 6.8e6) then begin
           mlind=where(fixgauss_mass.ml eq ml[i])
           mbhind=where(fixgauss_mass[mlind].inmbh eq mbhs[j])
           temp+=exp(-0.5*(fixgauss_mass[mlind[mbhind]].chi2))
        endif else begin
           mlind=where((final.ml gt ml[i]-0.09 and final.ml lt ml[i]-0.09) OR final.ml eq ml[i])
           mbhind=where(final[mlind].inmbh gt mbhs[j]-1. and final[mlind].inmbh lt mbhs[j]+1.,c)
           if (c gt 1 or c lt 1) then stop
           temp+=exp(-0.5*(final[mlind[mbhind]].chi2))
        endelse
     endfor
     likelihood[i]=temp
  endfor
  
  minchi=dblarr(nmls)
  temp=0D

  for i=0,nmls-1 do begin
     ind=where(fixgauss_mass.ml eq ml[i])
     fixgaussmin=MIN(fixgauss_mass[ind].chi2)
     tind=where((final.ml gt ml[i]-0.09 and final.ml lt ml[i]-0.09) or final.ml eq ml[i])
     finalmin=MIN(final[tind].chi2)
     if (fixgaussmin lt finalmin) then begin
        temp+= exp(-0.5*fixgaussmin)
     endif else begin
        temp+= exp(-0.5*finalmin)
     endelse
     minchi[i]=temp
  endfor
  likeminchi=minchi/max(minchi)
  testml=interpol(ml,likeminchi,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])*mlpop

  likefixgauss_mass=likelihood/max(likelihood)
  mlfixgauss_mass=(interpol(ml,likefixgauss_mass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  print,'Best fit M/L = ',mlfixgauss_mass[3]
  print,'1-sigma error M/L = ',mlfixgauss_mass[4]-mlfixgauss_mass[3],mlfixgauss_mass[3]-mlfixgauss_mass[2]
  print,'3-sigma error M/L = ',mlfixgauss_mass[6]-mlfixgauss_mass[3],mlfixgauss_mass[3]-mlfixgauss_mass[0]
  likelihood=dblarr(nmls)
  temp=0D
  a=where(fixgauss_mass.inmbh eq 0.)
  for i=0,nmls-1 do begin
     ind=where((fixgauss_mass[a].ml gt ml[i]-0.01 and fixgauss_mass[a].ml lt ml[i]+0.01) OR fixgauss_mass[a].ml eq ml[i],c)
     if (c gt 1 or c lt 1) then stop
     temp+=exp(-0.5*(fixgauss_mass[a[ind]].chi2))
     likelihood[i]=temp
  endfor
  likefixnorm_nomass=likelihood/max(likelihood)
  mlfixnorm_nomass=(interpol(ml,likefixnorm_nomass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])) *mlpop
  print,'Best fit M/L ratio= ', mlfixnorm_nomass[3]
  print,'3-sigma error M/L ratio= ', mlfixnorm_nomass[6]-mlfixnorm_nomass[3],mlfixnorm_nomass[3]-mlfixnorm_nomass[0]

  stop
  scale=0.05
  zeropoint=26.05923
  msun=5.11
  mlgin=2.26893d
  mlgout=4.94259d
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
  temp=where(fixgauss_mass.outmbh eq 0.)
  nobhmlind=where(fixgauss_mass[temp].chi2 eq min(fixgauss_mass[temp].chi2))
  nobhml=fixgauss_mass[temp[nobhmlind]].ml
  nobhmass=nobhml*((mlgin*inlum)+(mlgout*outlum))
  aveml=((influx*mlgin)+(outflux*mlgout))/(influx+outflux)
  
  print,'Total mass with black hole = ', final[bhmlind].inmbh, ' = ',bhmass
  print,'Dynamical M/L with black hole = ', bhml*aveml
  print,'Total mass w/o black hole = ', nobhmass
  print,'Dynamical M/L w/o black hole = ',nobhml*aveml
  print,'Total Luminosity = ', inlum + outlum
  intdisp=32900.
  r=((16.5e6)*(3.086e16))*((rad*0.05)*(4.848e-6))
  g=6.67e-11
  dynmass=(((intdisp^2)*r)/(g))/(1.989e30)
  print,'Total mass from integrated dispersion = ', dynmass
  stop
;  set_plot,'ps'
;  device,filename='./onedbestfit_BEST.ps',/color
;  myplot,filename='./onedbestfit_BEST.ps'
;  plotsym,0,1/2.,/fill
;  djs_plot,xbin,disp,psym=8,ytitle='Dispersion (km/s)',xtitle='Radius (arcseconds)',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2,/xlog
;  plotsym,0,1/2.,/fill,color=100
;  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
;  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
;  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
;  a=where(final.inmbh eq 0.)
;  c=where(final[a].chi2 eq min(final[a].chi2))
;  djs_oplot,xbin,final[a[c]].rms,color='red',thick=4
;  b=where(final.chi2 eq min(final.chi2))
;  djs_oplot,xbin,final[b].rms,color='blue',thick=4
;  device,/close
;  set_plot,'x'
  set_plot,'ps'
  device,filename='./onedbestfit_BEST_linear.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  
;  myplot,filename='./onedbestfit_BEST_linear.ps'
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  a=where(fixgauss_mass.outmbh eq 0.)
  c=where(fixgauss_mass[a].chi2 eq min(fixgauss_mass[a].chi2))
  djs_oplot,xbin,fixgauss_mass[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  fixgauss_mass=mrdfits('./temp/fixgausspsf_mass.fits',1)
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  ml=findgen(35)*0.1+0.1
  inclinations=[14.5,20.,30.,40.,50.,60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  nbeta=n_elements(betas)
;  a=where(fixgauss_mass.ininc eq 14.5)
;  fixgauss_mass=fixgauss_mass[a]
  fixgauss_mass=reform(fixgauss_mass,nmbhs,nbetas,ninclinations,nmls)
  a=where(fixgauss_mass.outmbh eq 0.)
  b=where(fixgauss_mass[a].chi2 eq min(fixgauss_mass[a].chi2))
  help,fixgauss_mass[a[b]]
  djs_oplot,xbin,fixgauss_mass[a[b]].rms,color='grey',thick=4
  xyouts,[0.02],[26],['M59cO'],charthick=3,charsize=1.5,/data
  device,/close

  set_plot,'x'

  stop

END
