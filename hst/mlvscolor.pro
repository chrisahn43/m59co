pro mlvscolor
  scale=0.05
  zeropoint_g_ab=26.05923
  zeropoint_z_ab=24.84245
  zeropoint_g_ve=26.16252
  zeropoint_z_ve=24.32305
  hstcolorin=1.42
  hstcolorout=1.6
  gsun=5.10 ;MAY NEED TO BE VEGA MAG
  zsun=4.52
  readcol,'sdss_ssp.dat',sdssz,sdssage,sdssmbol,u,g,r,i,z,format='F,F,F,F,F,F,F,F'
  readcol,'hstwfc_ssp.dat',hstz,hstage,hstmbol,f435w,f475w,f502n,f550m,f555w,f606w,f625w,f658n,f660n,f775w,f814w,f850lp,f892n,format='F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  hstcolor=(f475w+zeropoint_g_ab-zeropoint_g_ve)-(f850lp+zeropoint_z_ab-zeropoint_z_ve)
  sdsscolor=g-z
  djs_plot,hstcolor,sdsscolor,psym=4
  fit=linfit(hstcolor,sdsscolor)
  y=hstcolor*fit[1]+fit[0]
  djs_oplot,hstcolor,y,color='blue'
  stop
  colorin=(hstcolorin*fit[1]+fit[0])-(0.107-0.041)
  colorout=(hstcolorout*fit[1]+fit[0])-(0.107-0.041)

  readcol,'~/bc03/bc03/models/Padova1994/chabrier/bc2003_hr_m62_chab_ssp.1ABmag',logage1,mbol,gab,ug,gr,gi,gz,format='F,F,F,F,F,F,F'
  readcol,'~/bc03/bc03/models/Padova1994/chabrier/bc2003_hr_m62_chab_ssp.4color',logage4,mbol,bmag,vmag,mlb,mlv,mstar,mgas,mgal,sfr,format='F,F,F,F,F,F,F,F,F,F'

  inind=where(gz gt colorin-0.003 and gz lt colorin+0.003)
  outind=where(gz gt colorout-0.003 and gz lt colorout+0.003)

  ginlum=(10^(-0.4*(gab[inind]-gsun)))
  goutlum=(10^(-0.4*(gab[outind]-gsun)))

  vsun=4.80
  vinlum=(10^(-0.4*(vmag[inind]-vsun)))
  voutlum=(10^(-0.4*(vmag[outind]-vsun)))
;  min=mstar[inind]
;  mout=mstar[outind]

  mlin=mlv[inind]*(vinlum/ginlum)
  mlin=mlin[0]
  mlout=mlv[outind]*(voutlum/goutlum)
  mlout=mlout[0]
  print,mlin,mlout

  readcol,'m59co_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(pa gt 20.)
  intensity[a]=intensity[a]*mlin
  b=where(pa lt 20.)
  intensity[b]=intensity[b]*mlout
  forprint,intensity,sigma,q,pa,textout='m59co_mge_outputsersic_mass.dat',format='D,D,D,D'
  
  stop
                                ;NOW COMPUTE M/L Directly from HST filters in new BC03 Models
  colorin=1.42-(0.107-0.041)
  colorout=1.6-(0.107-0.041)
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f475w-f850lp
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'
  inind=where(hstcolor gt colorin-0.001 and hstcolor lt colorin+0.001)
  outind=where(hstcolor gt colorout-0.02 and hstcolor lt colorout+0.02)
  ginlum=(10^(-0.4*(f475w[inind]-gsun)))
  goutlum=(10^(-0.4*(f475w[outind]-gsun)))

  vsun=4.80
  vinlum=(10^(-0.4*(vmag[inind]-vsun)))
  voutlum=(10^(-0.4*(vmag[outind]-vsun)))
;  min=mstar[inind]
;  mout=mstar[outind]

  mlin=mlvtot[inind]*(vinlum/ginlum)
  mlin=mlin[0]
  mlout=mlvtot[outind]*(voutlum/goutlum)
  mlout=mlout[0]
  print,mlin,mlout
  readcol,'m59co_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(pa gt 20.)
  intensity[a]=intensity[a]*mlin
  b=where(pa lt 20.)
  intensity[b]=intensity[b]*mlout
  forprint,intensity,sigma,q,pa,textout='m59co_mge_outputsersic_nmass.dat',format='D,D,D,D'

  stop
  readcol,'m59co_mge_outputsersic_free.dat',glum,gsig,gq,gpa,format='D,D,D,D'
  radius=gsig/scale
  minradius=[0,radius[0:5],0,radius[7:12]]
  fits_read,'./make_decon/mgefreemodelg.fits',gimg
  fits_read,'./make_decon/mgefreemodelz.fits',zimg
  find_galaxy,gimg,m,e,a,xcg,ycg
  find_galaxy,zimg,m,e,a,xcz,ycz
  gmag=fltarr(n_elements(radius))
  zmag=fltarr(n_elements(radius))

  for i=0,n_elements(radius)-1 do begin
     area=((!PI*radius[i]^2)-(!PI*minradius[i]^2))

     aper,gimg,xcg,ycg,gmaxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,gimg,xcg,ycg,gminflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     aper,zimg,xcz,ycz,zmaxflux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     aper,zimg,xcz,ycz,zminflux,fluxerr,0.,skyerr,1,minradius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact

     if (minradius[i] eq 0.) then begin
        gflux=gmaxflux
        zflux=zmaxflux
     endif else begin
        gflux=gmaxflux-gminflux
        zflux=zmaxflux-zminflux
     endelse
     if (zflux lt 0.) then stop
     gmag[i]=zeropoint_g_ab+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gflux)-0.107
     zmag[i]=zeropoint_z_ab+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zflux)-0.041
  endfor
  color=gmag-zmag

  color=color*fit[1]+fit[0]

  masstolight=fltarr(n_elements(color))
  for i=0,n_elements(color)-1 do begin
     cind=where(gz lt color[i]+0.001 and gz gt color[i]-0.001,count)
     if (count eq 1.) then begin
        mass=mean(mstar[cind])
        light=mean(10^((-0.4*(gab[cind]-gsun))))
        masstolight[i]=mass/light
     endif else begin
        cind=where(gz lt color[i]+0.002 and gz gt color[i]-0.002,count)
        if (count eq 1.) then begin
           mass=mean(mstar[cind])
           light=mean(10^((-0.4*(gab[cind]-gsun))))
           masstolight[i]=mass/light
        endif else begin
           cind=where(gz lt color[i]+0.004 and gz gt color[i]-0.004,count)
           if (count ge 1.) then begin
              mass=mean(mstar[cind])
              light=mean(10^((-0.4*(gab[cind]-gsun))))
              masstolight[i]=mass/light
           endif else stop
        endelse
     endelse
  endfor



  readcol,'m59co_mge_outputsersic_free.dat',intensity,sigma,q,pa,format='D,D,D,D'
  intensity=intensity*masstolight
  forprint,intensity,sigma,q,pa,textout='m59co_mge_outputsersic_mass_free.dat',format='D,D,D,D'
  
  stop
                                ;FREE MODEL M/L FROM HST BC03
  color=gmag-zmag
  masstolight=fltarr(n_elements(color))
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f475w-f850lp
  gsun=5.10
  vsun=4.80
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'
  for i=0,n_elements(color)-1 do begin
     cind=where(hstcolor lt color[i]+0.001 and hstcolor gt color[i]-0.001,count)
     if (count eq 1.) then begin
        glight=mean(10^((-0.4*(f475w[cind]-gsun))))
        vlight=mean(10^((-0.4*(vmag[cind]-vsun))))
        masstolight[i]=mlvtot[cind]*(vlight/glight)
     endif else begin
        cind=where(hstcolor lt color[i]+0.002 and hstcolor gt color[i]-0.002,count)
        if (count eq 1.) then begin
           vlight=mean(10^((-0.4*(vmag[cind]-vsun))))
           glight=mean(10^((-0.4*(f475w[cind]-gsun))))
           masstolight[i]=mlvtot[cind]*(vlight/glight)
        endif else begin
           cind=where(hstcolor lt color[i]+0.004 and hstcolor gt color[i]-0.004,count)
           if (count ge 1.) then begin
              vlight=mean(10^((-0.4*(vmag[cind]-vsun))))
              glight=mean(10^((-0.4*(f475w[cind]-gsun))))
              masstolight[i]=mean(mlvtot[cind]*(vlight/glight))
           endif else begin
              cind=where(hstcolor lt color[i]+0.016 and hstcolor gt color[i]-0.016,count)
              if (count ge 1.) then begin
                 vlight=mean(10^((-0.4*(vmag[cind]-vsun))))
                 glight=mean(10^((-0.4*(f475w[cind]-gsun))))
                 masstolight[i]=mean(mlvtot[cind]*(vlight/glight))
              endif else stop
           endelse
        endelse
     endelse
  endfor
  readcol,'m59co_mge_outputsersic_free.dat',intensity,sigma,q,pa,format='D,D,D,D'
  intensity=intensity*masstolight
  forprint,intensity,sigma,q,pa,textout='m59co_mge_outputsersic_nmass_free.dat',format='D,D,D,D'
  stop
END
