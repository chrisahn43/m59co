PRO COMPARE_JAM_OUTPUT
                                ;THIS CODE CREATES FIGURES 5, 6, 7,
                                ;and 8 for M59cO in VUCD3 & M59cO BH
                                ;PAPER. STARTS BY READING IN
                                ;KINEMATICS AND ALL OF THE DYNAMICAL
                                ;MODELS. CREATES CONTOUR PLOTS FOR THE
                                ;DEFAULT DYNAMICAL MODEL. COMPUTES
                                ;CUMULATIVE LIKELIHOOD FOR ALL MODELS
                                ;ASSUMING ISOTROPY. THEN CREATES
                                ;DISPERSION PROFILE FITS FOR THE
                                ;DEFAULT MODEL.
  ;Read in Kinematics file for later use...
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
  ;Read in dynamical models.
  fixgauss_mass_k=mrdfits('./fixgausspsf_mass_k.fits',1)
  fixmoffat_mass_k=mrdfits('./fixmoffatpsf_mass_k.fits',1)
  fixgauss_k=mrdfits('./fixgausspsf_k.fits',1)
  fixgauss_mass_g=mrdfits('./fixgausspsf_mass_g.fits',1)
  freegauss_mass_k=mrdfits('./freegausspsf_mass_k.fits',1)
  freegauss_k=mrdfits('./freegausspsf_k.fits',1)
  fixold_mass_k=mrdfits('./fixoldpsf_mass_k.fits',1)
  ;CREATE CONTOUR PLOTS
  minchi=dblarr(nmbhs,nbeta)
  data=reform(fixgauss_mass_k,nmbhs,nbetas,ninclinations,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nbetas-1 do begin
        minchi[i,j]=MIN(data[i,j,*,*].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betambh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass_k.outmbh,fixgauss_mass_k.inbeta,YTITLE='!7b!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(betas)-0.1,max(betas)+0.1],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] ;,/XLOG
  OPLOT,fixgauss_mass_k.outmbh,fixgauss_mass_k.inbeta,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  min=MIN(fixgauss_mass_k.chi2,minpos)
  CONTOUR,SMOOTH(minchi,[2,8]),mbhs,betas,/OVERPLOT,LEVELS=[fixgauss_mass_k[minpos].chi2+2.3,fixgauss_mass_k[minpos].chi2+6.18,fixgauss_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]
  xx=REPLICATE(fixgauss_mass_k[minpos].outmbh,N_ELEMENTS(fixgauss_mass_k.outmbh))
  yy=REPLICATE(fixgauss_mass_k[minpos].inbeta,N_ELEMENTS(fixgauss_mass_k.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  a=where(fixgauss_mass_k.ininc eq 40.)
  b=where(fixgauss_mass_k[a].inmbh eq 5.e5)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].inbeta,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=145
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,-.3,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=145
  onepercent=fixgauss_mass_k[a[b[c]]].rms
;  stop
  b=where(fixgauss_mass_k[a].inmbh eq 1.5e6)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].inbeta,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=220
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,-.3,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=220
  fivepercent=fixgauss_mass_k[a[b[c]]].rms
;  stop
  b=where(fixgauss_mass_k[a].inmbh eq 3.e6)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].inbeta,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=200
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,-.3,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].inbeta,/data,thick=10,hsize=0.01,color=200
  tenpercent=fixgauss_mass_k[a[b[c]]].rms
  xyouts,[7.3e6],[0.83],['M59cO'],charthick=8,charsize=1.5,/data

;  stop
  set_plot,'x'
  set_plot,'ps'
  device,filename='onedanisotropy.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  loadct,0,/silent

  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  loadct,40,/silent
  djs_oplot,xbin,onepercent,thick=4,color=145
  djs_oplot,xbin,fivepercent,thick=4,color=220
  djs_oplot,xbin,tenpercent,thick=4,color=200
  items=['1% BH Mass','5% BH Mass', '10% BH Mass']
  lines=[0,0,0]
  color=[145,220,200]
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/right
  xyouts,[0.02],[26],['M59cO'],charthick=3,charsize=1.5,/data

  set_plot,'x'
  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[i,*,*,j].chi2)
     endfor
  endfor

  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='mlmbh.ps',/color,/HELVETICA,bits=4,/cmyk,/encaps,/inches,ysize=7;BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass_k.outmbh,fixgauss_mass_k.ml,YTITLE='!7C!3',XTITLE='M!IBH!N [M!D!9n!3!N]',XCHARS=1.2,/NODATA,YCHARS=1.3,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,1.2],CHARSIZE=1.25,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 
  OPLOT,fixgauss_mass_k.outmbh,fixgauss_mass_k.ml,PSYM=6,THICK=10,SYMSIZE=0.2,color=100
  min=MIN(fixgauss_mass_k.chi2,minpos)
  CONTOUR,SMOOTH(minchi,[1,1],/EDGE_TRUNCATE,/NAN),mbhs,ml,/OVERPLOT,LEVELS=[fixgauss_mass_k[minpos].chi2+2.3,fixgauss_mass_k[minpos].chi2+6.18,fixgauss_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass_k[minpos].outmbh,N_ELEMENTS(fixgauss_mass_k.outmbh))
  yy=REPLICATE(fixgauss_mass_k[minpos].ml,N_ELEMENTS(fixgauss_mass_k.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=0.7,color=150
  a=where(fixgauss_mass_k.ininc eq 40.)
  b=where(fixgauss_mass_k[a].inmbh eq 5.e5)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].ml,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=145
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,0.0,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=145
  onepercent=fixgauss_mass_k[a[b[c]]].rms
  b=where(fixgauss_mass_k[a].inmbh eq 1.5e6)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].ml,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=220
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,0.0,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=220
  fivepercent=fixgauss_mass_k[a[b[c]]].rms
  b=where(fixgauss_mass_k[a].inmbh eq 3.e6)
  c=where(fixgauss_mass_k[a[b]].chi2 eq min(fixgauss_mass_k[a[b]].chi2))
  arrow,0.,fixgauss_mass_k[a[b[c]]].ml,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=200
  arrow,fixgauss_mass_k[a[b[c]]].inmbh,0.0,fixgauss_mass_k[a[b[c]]].inmbh,fixgauss_mass_k[a[b[c]]].ml,/data,thick=10,hsize=0.01,color=200
  tenpercent=fixgauss_mass_k[a[b[c]]].rms
  xyouts,[7.3e6],[1.13],['M59cO'],charthick=8,charsize=1.5,/data

  set_plot,'x'
  minchi=dblarr(nbeta,nmls)

  for i=0,nbetas-1 do begin
     for j=0,nmls-1 do begin
        minchi[i,j]=MIN(data[*,i,*,j].chi2)
     endfor
  endfor

  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='betaml.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass_k.inbeta,fixgauss_mass_k.ml,YTITLE='!7C!3',XTITLE='Beta',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(ml)-0.1,1.2],CHARSIZE=1,XSTY = 1,XRANGE=[min(betas)-0.1,max(betas)+0.1],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 
  OPLOT,fixgauss_mass_k.inbeta,fixgauss_mass_k.ml,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass_k.chi2,minpos)
  CONTOUR,minchi,betas,ml,/OVERPLOT,LEVELS=[fixgauss_mass_k[minpos].chi2+2.3,fixgauss_mass_k[minpos].chi2+6.18,fixgauss_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass_k[minpos].inbeta,N_ELEMENTS(fixgauss_mass_k.inbeta))
  yy=REPLICATE(fixgauss_mass_k[minpos].ml,N_ELEMENTS(fixgauss_mass_k.inbeta))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150
  ind=where(fixgauss_mass_k.inbeta eq 0.)
  fixgauss_mass_k=fixgauss_mass_k[ind]
  minchi=dblarr(nmbhs,ninclinations)
  for i=0,nmbhs-1 do begin
     for j=0,ninclinations-1 do begin
        mbh=mbhs[i] & inc=inclinations[j]
        ind=where(fixgauss_mass_k.ininc eq inc and fixgauss_mass_k.outmbh eq mbh)
        minchi[i,j]=MIN(fixgauss_mass_k[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incmbh_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass_k.outmbh,fixgauss_mass_k.ininc,YTITLE='Inclination',XTITLE='MBH',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(mbhs)-2d4,MAX(mbhs)+2d4],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 
  OPLOT,fixgauss_mass_k.outmbh,fixgauss_mass_k.ininc,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass_k.chi2,minpos)
  CONTOUR,minchi,mbhs,inclinations,/OVERPLOT,LEVELS=[fixgauss_mass_k[minpos].chi2+2.3,fixgauss_mass_k[minpos].chi2+6.18,fixgauss_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass_k[minpos].outmbh,N_ELEMENTS(fixgauss_mass_k.outmbh))
  yy=REPLICATE(fixgauss_mass_k[minpos].ininc,N_ELEMENTS(fixgauss_mass_k.outmbh))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150

  minchi=dblarr(nmls,ninclinations)
  for i=0,nmls-1 do begin
     for j=0,ninclinations-1 do begin
        mls=ml[i] & inc=inclinations[j]
        ind=where(fixgauss_mass_k.ininc eq inc and fixgauss_mass_k.ml eq mls)
        minchi[i,j]=MIN(fixgauss_mass_k[ind].chi2)
     endfor
  endfor
  set_plot,'ps' & LOADCT,40,/Silent
  device,filename='incml_nobeta.ps',BIT=8,XSIZE=7,YSIZE=6,/INCHES,SET_FONT='Times',_EXTRA=extra,/TT,/COLOR
  plot,fixgauss_mass_k.ml,fixgauss_mass_k.ininc,YTITLE='Inclination',XTITLE='M/L',XCHARS=1,/NODATA,YCHARS=1,CHARTHICK=3,XTHICK=3,PSYM=3,Yrange=[min(inclinations)-5,max(inclinations)+5],CHARSIZE=1,XSTY = 1,XRANGE=[min(ml)-0.1,1.2],YTHICK=3,YSTY=1,POSITION=[0.15,0.22,0.98,0.98] 
  OPLOT,fixgauss_mass_k.ml,fixgauss_mass_k.ininc,PSYM=6,THICK=10,SYMSIZE=1,color=100

  min=MIN(fixgauss_mass_k.chi2,minpos)
  CONTOUR,minchi,ml,inclinations,/OVERPLOT,LEVELS=[fixgauss_mass_k[minpos].chi2+2.3,fixgauss_mass_k[minpos].chi2+6.18,fixgauss_mass_k[minpos].chi2+11.83],C_LINESTYLE=[0,0,0],C_THICK=[15,5,5],C_COLORS=[0,50,250]


  xx=REPLICATE(fixgauss_mass_k[minpos].ml,N_ELEMENTS(fixgauss_mass_k.ml))
  yy=REPLICATE(fixgauss_mass_k[minpos].ininc,N_ELEMENTS(fixgauss_mass_k.ml))
  OPLOT,xx,yy,psym=6,thick=15,symsize=2,color=150


  set_plot,'x'
  stop
  ind=where(fixgauss_mass_k.inbeta eq 0.)
  fixgauss_mass_k=fixgauss_mass_k[ind]
  ind=where(fixmoffat_mass_k.inbeta eq 0.)
  fixmoffat_mass_k=fixmoffat_mass_k[ind]
  ind=where(fixgauss_k.inbeta eq 0.)
  fixgauss_k=fixgauss_k[ind]
  ind=where(fixgauss_mass_g.inbeta eq 0.)
  fixgauss_mass_g=fixgauss_mass_g[ind]
  ind=where(freegauss_mass_k.inbeta eq 0.)
  freegauss_mass_k=freegauss_mass_k[ind]
  ind=where(freegauss_k.inbeta eq 0.)
  freegauss_k=freegauss_k[ind]
  ind=where(fixold_mass_k.inbeta eq 0.)
  fixold_mass_k=fixold_mass_k[ind]
  likelihood=dblarr(nmbhs)
  temp=0.D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixgauss_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass_k=likelihood/max(likelihood)
  mbhsfixgauss_mass_k=interpol(mbhs,likefixgauss_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  sigma=[-3,-2,-1,0,1,2,3]
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixgauss_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_k=likelihood/max(likelihood)
  mbhsfixgauss_k=interpol(mbhs,likefixgauss_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixgauss_mass_g.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixgauss_mass_g[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixgauss_mass_g[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixgauss_mass_g[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixgauss_mass_g=likelihood/max(likelihood)
  mbhsfixgauss_mass_g=interpol(mbhs,likefixgauss_mass_g,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freegauss_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freegauss_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss_mass_k=likelihood/max(likelihood)
  mbhsfreegauss_mass_k=interpol(mbhs,likefreegauss_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(freegauss_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(freegauss_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(freegauss_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(freegauss_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefreegauss_k=likelihood/max(likelihood)
  mbhsfreegauss_k=interpol(mbhs,likefreegauss_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0D
  for i=0,nmbhs-1 do begin
     mbhind=where(fixmoffat_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixmoffat_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixmoffat_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixmoffat_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixmoffat_mass_k=likelihood/max(likelihood)
  mbhsfixmoffat_mass_k=interpol(mbhs,likefixmoffat_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  temp=0.D

  for i=0,nmbhs-1 do begin
     mbhind=where(fixold_mass_k.outmbh eq mbhs[i])
     for j=0,nmls-1 do begin
        mlind=where(fixold_mass_k[mbhind].ml eq ml[j])
        for k=0,ninclinations-1 do begin
           incind=where(fixold_mass_k[mbhind[mlind]].ininc eq inclinations[k])
           temp+=exp(-0.5*(fixold_mass_k[mbhind[mlind[incind]]].chi2))
        endfor
     endfor
     likelihood[i]=temp
  endfor
  likefixold_mass_k=likelihood/max(likelihood)
  mbhsfixold_mass_k=interpol(mbhs,likefixold_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  stop

;  mbhs=[findgen(17)*1.e5+5.2e6]
;  inclinations=14.5
;  beta=0.
;  ml=findgen(35)*0.1+0.1
;  nmbhs=n_elements(mbhs)
;  nmls=n_elements(ml)
;  final=REPLICATE({inmbh:0.0,chi2:0.0,ml:0.0,rms:fltarr(n_elements(disp))},nmbhs,nmls)
;  readcol,'~/research/code/gemini15/m59co/kinematics/kinematic_psf_gauss.dat',normpsf,sigmapsf,format='F,F'
  readcol,'../make_decon/m59co_mge_outputsersic_k.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  readcol,'../m59co_mge_outputsersic_nmass.dat',mass,format='D'
;  surf_lum = intensity
;  sigma_lum = sigmaarc
;  qobs_lum = q
;  surf_pot = mass
;  sigma_pot = sigmaarc
;  qobs_pot = q
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
;  mwrfits,final,'finalgrid.fits'
;  stop
  final=mrdfits('finalgrid.fits',1)
  final=reform(final,17,35)
  stop
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6,findgen(3)*1.e5+5.2e6,findgen(4)*1.e5+5.6e6,findgen(4)*1.e5+6.1e6,findgen(3)*1.e5+6.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)
  ml=[findgen(35)*0.1+0.1];,findgen(35)*0.1+0.05]
  nmls=n_elements(ml)
  ind=where(fixgauss_mass_k.ininc eq 14.5)
  fixgauss_mass_k=fixgauss_mass_k[ind]
  fixgauss_mass_k=reform(fixgauss_mass_k,19,35)
  finalbh=[fixgauss_mass_k.outmbh,final.inmbh]
  sorted=sort(finalbh)
  finalbh=finalbh[sorted]
  finalml=[fixgauss_mass_k.ml,final.ml]
  finalml=finalml[sorted]
  minchi=dblarr(nmbhs,nmls)
  for i=0,nmbhs-1 do begin
     for j=0,nmls-1 do begin
        if (mbhs[i] lt 5.2e6 OR mbhs[i] gt 6.8e6) then begin
           mbhind=where(fixgauss_mass_k.inmbh eq mbhs[i])
           mlind=where(fixgauss_mass_k[mbhind].ml eq ml[j])
           minchi[i,j]=fixgauss_mass_k[mbhind[mlind]].chi2
        endif else begin
           mbhind=where(final.inmbh gt mbhs[i]-1 and final.inmbh lt mbhs[i]+1)
           mlind=where((final[mbhind].ml gt ml[j]-0.09 and final[mbhind].ml lt ml[j]-0.09) OR final[mbhind].ml eq ml[j])
           minchi[i,j]=final[mbhind[mlind]].chi2
        endelse
     endfor
  endfor
            
  finalchi2=[fixgauss_mass_k.chi2,final.chi2]
  finalchi2=finalchi2[sorted]
   likelihood=dblarr(nmbhs)
  
  temp=0D
  for i=0,nmbhs-1 do begin
     if (mbhs[i] lt 5.2e6 OR mbhs[i] gt 6.8e6) then begin
        mbhind=where(fixgauss_mass_k.inmbh eq mbhs[i])
        for j=0,nmls-1 do begin
           mlind=where(fixgauss_mass_k[mbhind].ml eq ml[j])
           temp+=exp(-0.5*(fixgauss_mass_k[mbhind[mlind]].chi2))
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
        ind=where(fixgauss_mass_k.outmbh eq mbhs[i])
        temp+=exp(-0.5*(min(fixgauss_mass_k[ind].chi2)))
     endif else begin
        ind=where(final.inmbh eq mbhs[i])
        temp+=exp(-0.5*(min(final[ind].chi2)))
     endelse
     minchi[i]=temp
  endfor
  likeminchi=minchi/max(minchi)
  testmbh=interpol(mbhs,likeminchi,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  likefixgauss_mass_k=likelihood/max(likelihood)
  mbhsfixgauss_mass_k=interpol(mbhs,likefixgauss_mass_k,[.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])
  print,'Best fit= ',mbhsfixgauss_mass_k[3]
  print,'1-sigma errors = ',mbhsfixgauss_mass_k[4]-mbhsfixgauss_mass_k[3],mbhsfixgauss_mass_k[3]-mbhsfixgauss_mass_k[2]
  print,'3-sigma errors = ',mbhsfixgauss_mass_k[6]-mbhsfixgauss_mass_k[3],mbhsfixgauss_mass_k[3]-mbhsfixgauss_mass_k[0]
;  myplot,filename='./cummulike.ps'
  set_plot,'ps'
  device,filename='./cummulike.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  LOADCT,0,/SILENT
  djs_plot,mbhs,likefixgauss_mass_k,ytitle='Cumulative Likelihood',xtitle='M_{BH} [M_{\odot}]',charsize=1.5,charthick=4,xthick=3,ythick=3,thick=12,/nodata
  arrow,mbhsfixgauss_mass_k[3],0.,mbhsfixgauss_mass_k[3],1.,/data,thick=8.,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass_k[2],0.,mbhsfixgauss_mass_k[2],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass_k[4],0.,mbhsfixgauss_mass_k[4],1.,/data,thick=4,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass_k[0],0.,mbhsfixgauss_mass_k[0],1.,/data,thick=2,color=100,hsize=0.01
  arrow,mbhsfixgauss_mass_k[6],0.,mbhsfixgauss_mass_k[6],1.,/data,thick=2,color=100,hsize=0.01

  djs_oplot,mbhs,likefixgauss_mass_k,thick=12
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]

  djs_oplot,mbhs,likefixmoffat_mass_k,thick=4,color='red'
  djs_oplot,mbhs,likefixold_mass_k,thick=4,color='red',linestyle=1
  djs_oplot,mbhs,likefreegauss_mass_k,thick=4,color='blue'
  djs_oplot,mbhs,likefixgauss_k,thick=4,color='blue',linestyle=1
  djs_oplot,mbhs,likefreegauss_k,thick=4,color='blue',linestyle=2
  djs_oplot,mbhs,likefixgauss_mass_g,thick=4,color='cyan'
  items=['Best Fit Model','PSF Variations','Mass Model Variations', 'Luminosity Model Variations']
  lines=[0,0,0,0]
  color=['black','red','blue','cyan']
  al_legend,items,linestyle=lines,colors=color,background_color='white',charthick=4,thick=3,/top,/left
  xyouts,[8.5e6],[0.05],['M59cO'],charthick=3,charsize=1.5,/data
  device,/close
  set_plot,'x'
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6,findgen(3)*1.e5+5.2e6,findgen(4)*1.e5+5.6e6,findgen(4)*1.e5+6.1e6,findgen(3)*1.e5+6.6e6]
  mbhs=mbhs[sort(mbhs)]
  nmbhs=n_elements(mbhs)
  mlpop=4.8340184
    likelihood=dblarr(nmls)
  temp=0D
  for i=0,nmls-1 do begin
     for j=0,nmbhs-1 do begin
        if (mbhs[j] lt 5.2e6 OR mbhs[j] gt 6.8e6) then begin
           mlind=where(fixgauss_mass_k.ml eq ml[i])
           mbhind=where(fixgauss_mass_k[mlind].inmbh eq mbhs[j])
           temp+=exp(-0.5*(fixgauss_mass_k[mlind[mbhind]].chi2))
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
     ind=where(fixgauss_mass_k.ml eq ml[i])
     fixgaussmin=MIN(fixgauss_mass_k[ind].chi2)
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

  likefixgauss_mass_k=likelihood/max(likelihood)
  mlfixgauss_mass_k=(interpol(ml,likefixgauss_mass_k,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865]))*mlpop
  print,'Best fit M/L = ',mlfixgauss_mass_k[3]
  print,'1-sigma error M/L = ',mlfixgauss_mass_k[4]-mlfixgauss_mass_k[3],mlfixgauss_mass_k[3]-mlfixgauss_mass_k[2]
  print,'3-sigma error M/L = ',mlfixgauss_mass_k[6]-mlfixgauss_mass_k[3],mlfixgauss_mass_k[3]-mlfixgauss_mass_k[0]
  likelihood=dblarr(nmls)
  temp=0D
  a=where(fixgauss_mass_k.inmbh eq 0.)
  for i=0,nmls-1 do begin
     ind=where((fixgauss_mass_k[a].ml gt ml[i]-0.01 and fixgauss_mass_k[a].ml lt ml[i]+0.01) OR fixgauss_mass_k[a].ml eq ml[i],c)
     if (c gt 1 or c lt 1) then stop
     temp+=exp(-0.5*(fixgauss_mass_k[a[ind]].chi2))
     likelihood[i]=temp
  endfor
  likefixnorm_nomass=likelihood/max(likelihood)
  mlfixnorm_nomass=(interpol(ml,likefixnorm_nomass,[0.00135,0.02275,0.16,0.5,0.84,0.97725,0.99865])) *mlpop
  print,'Best fit M/L ratio= ', mlfixnorm_nomass[3]
  print,'3-sigma error M/L ratio= ', mlfixnorm_nomass[6]-mlfixnorm_nomass[3],mlfixnorm_nomass[3]-mlfixnorm_nomass[0]

  stop
  scale=0.05
  zeropoint=26.05923
  msun=5.10
  mlgin=2.849d
  mlgout=5.478d
  rad=50.;6.337
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  readcol,'../m59co_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
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
  temp=where(fixgauss_mass_k.outmbh eq 0.)
  nobhmlind=where(fixgauss_mass_k[temp].chi2 eq min(fixgauss_mass_k[temp].chi2))
  nobhml=fixgauss_mass_k[temp[nobhmlind]].ml
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
  set_plot,'ps'
  device,filename='./onedbestfit_BEST_linear.ps',/color,/HELVETICA,bits=8,/cmyk,/encaps,/inches,ysize=7
  
;  myplot,filename='./onedbestfit_BEST_linear.ps'
  plotsym,0,1/2.,/fill
  djs_plot,xbin,disp,psym=8,ytitle='Dispersion [km s^{-1}]',xtitle='Radius ["]',xran=[0.015,0.55],yran=[25,45],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=2
  plotsym,0,1/2.,/fill,color=100
  oploterror,xbin,disp,xbinerror,disperr,/NOHAT,psym=3
  oploterror,xbinout,dispout,xbinouterror,disperrout,/NOHAT,ERRCOLOR='grey',psym=3
  djs_oplot,xbinout,dispout,psym=8,color='grey',charsize=1.5,charthick=4,symsize=2
  a=where(fixgauss_mass_k.outmbh eq 0.)
  c=where(fixgauss_mass_k[a].chi2 eq min(fixgauss_mass_k[a].chi2))
  djs_oplot,xbin,fixgauss_mass_k[a[c]].rms,color='red',thick=4
  b=where(final.chi2 eq min(final.chi2))
  djs_oplot,xbin,final[b].rms,color='blue',thick=4
  fixgauss_mass_k=mrdfits('./fixgausspsf_mass_k.fits',1)
  mbhs=[0.,1.e5,5.e5,1.e6,1.5e6,2.e6,2.5e6,3.e6,3.5e6,4.e6,4.5e6,5.e6,5.5e6,6.e6,6.5e6,7.e6,7.5e6,8.e6,8.5e6]
  ml=findgen(35)*0.1+0.1
  inclinations=[14.5,20.,30.,40.,50.,60.,70.,80.,90.]
  betas=findgen(11)*0.1-0.2
  nmbhs=n_elements(mbhs)
  nmls=n_elements(ml)
  nbetas=n_elements(betas)
  ninclinations=n_elements(inclinations)
  nbeta=n_elements(betas)
  fixgauss_mass_k=reform(fixgauss_mass_k,nmbhs,nbetas,ninclinations,nmls)
  a=where(fixgauss_mass_k.outmbh eq 0.)
  b=where(fixgauss_mass_k[a].chi2 eq min(fixgauss_mass_k[a].chi2))
  help,fixgauss_mass_k[a[b]]
  djs_oplot,xbin,fixgauss_mass_k[a[b]].rms,color='grey',thick=4
  xyouts,[0.02],[26],['M59cO'],charthick=3,charsize=1.5,/data
  device,/close

  set_plot,'x'

  stop

END
