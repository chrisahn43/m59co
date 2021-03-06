PRO TESTROTATION

  fwhmfile='./../sky/m59co_disp_med.fits'
  contfile='./lum_model/m59co_combine_goodmay_cont.fits'
  fwhm=READFITS(fwhmfile)
  cont=READFITS(contfile)
  infile='./m59co_combine_goodmay.fits'
  cube=DOUBLE(READFITS(infile,head,ext=1,/SILENT))
  varcube=DOUBLE(READFITS(infile,head,ext=2,/SILENT))
  cubesize=size(cube,/dim)
  lambda0=SXPAR(head,'CRVAL3')
  dlambda=SXPAR(head,'CD3_3')
  lambda=FINDGEN(cubesize[2])*dlambda+lambda0
  uoutspec=DBLARR(cubesize[2])
  uoutrawspec=DBLARR(cubesize[2])
  uoutvar=DBLARR(cubesize[2])
  uoutsky=DBLARR(cubesize[2])
  loutspec=DBLARR(cubesize[2])
  loutrawspec=DBLARR(cubesize[2])
  loutvar=DBLARR(cubesize[2])
  loutsky=DBLARR(cubesize[2])

  radii=[1.75,2.5,3.3,4.3,5.5,7.2,12.];,22.] took out first point of 1
  minrad=[0.,1.75,2.5,3.3,4.3,5.5,7.2];,12.] took out second point of 1
  angle=findgen(73)*5.+2.;[252]
  angle=angle*!DTOR
  xcen=38.06
  ycen=45.44
  makex,cont,xarr,yarr,/zero
  rarr=SQRT((xarr-xcen)^2+(yarr-ycen)^2)
;  OPENW,1,'./intspec/spectra_fwhm_goodmay.dat'
  for j=0,n_elements(radii)-1 do begin
     OPENW,1,'./intspec/spectra_fwhm_goodmay.dat'
     for i=0,n_elements(angle)-1 do begin
        ind=where(rarr gt minrad[j] and rarr le radii[j])
        m=tan(angle[i])
        equation=(m*(xarr[ind]-xcen))+ycen
        uind=where((yarr[ind]-equation) ge 0.)
        uxpix=xarr[ind[uind]]
        uypix=yarr[ind[uind]]
        lind=where((yarr[ind]-equation) lt 0.)
        lxpix=xarr[ind[lind]]
        lypix=yarr[ind[lind]]
        djs_plot,uxpix,uypix,psym=4,xran=[25,50],yran=[30,60]
        djs_oplot,lxpix,lypix,psym=4,color='red'
        wait,1
;        stop
        avfwhm=TOTAL(fwhm[ind]*cont[ind])/TOTAL(cont[ind])
        avrad=TOTAL(rarr[ind]*cont[ind])/TOTAL(cont[ind])
        outfile='m59co_up'+STRTRIM(FIX(angle[i]*!RADEG),2)+'-'+STRTRIM(FIX(minrad[j]),2)+'-'+STRTRIM(FIX(radii[j]),2)+'.fits'
        printf,1,minrad[j],radii[j],avrad,avfwhm,' ',outfile,FORMAT='(2F5.1,2F7.3,2A)'
        FOR k=0,cubesize[2]-1 DO BEGIN
           skyval=MEDIAN(cube[5:20,67:81,k])
           uoutsky[k]=skyval*!PI*(radii[j]^2-minrad[j]^2)
           skyval=0
           uflux=total(cube[uxpix,uypix,k])
           uvarflux=total(varcube[uxpix,uypix,k])
           uoutspec[k]=uflux
           uoutvar[k]=uvarflux
           uoutrawspec[k]=uflux
        endfor
        ind=WHERE(lambda LT 2.01e4 or lambda GT 2.43e4)
        uoutspec[ind]=0.0
        uoutvar[ind]=0.0
        uoutsky[ind]=0.0

        uoutarr=[[uoutspec],[uoutvar],[uoutsky]]
        
        mkhdr,headout,uoutarr
        sxaddpar,headout,'CRVAL1',lambda0
        sxaddpar,headout,'CRPIX1',1
        sxaddpar,headout,'CD1_1',dlambda
        SXADDPAR,head,'BUNIT','erg/cm2/s/A'
        WRITEFITS,'./intspec/'+outfile,uoutarr,headout
     endfor
     for i=0,n_elements(angle)-1 do begin
        ind=where(rarr gt minrad[j] and rarr le radii[j])
        m=tan(angle[i])
        equation=(m*(xarr[ind]-xcen))+ycen
        lind=where((yarr[ind]-equation) lt 0.)
        lxpix=xarr[ind[lind]]
        lypix=yarr[ind[lind]]
        avfwhm=TOTAL(fwhm[ind]*cont[ind])/TOTAL(cont[ind])
        avrad=TOTAL(rarr[ind]*cont[ind])/TOTAL(cont[ind])
        outfile='m59co_down'+STRTRIM(FIX(angle[i]*!RADEG),2)+'-'+STRTRIM(FIX(minrad[j]),2)+'-'+STRTRIM(FIX(radii[j]),2)+'.fits'

        printf,1,minrad[j],radii[j],avrad,avfwhm,' ',outfile,FORMAT='(2F5.1,2F7.3,2A)'
        FOR k=0,cubesize[2]-1 DO BEGIN
           skyval=MEDIAN(cube[5:20,67:81,k])
           loutsky[k]=skyval*!PI*(radii[j]^2-minrad[j]^2)
           skyval=0
           lflux=total(cube[lxpix,lypix,k])
           lvarflux=total(varcube[lxpix,lypix,k])
           loutspec[k]=lflux
           loutvar[k]=lvarflux
           loutrawspec[k]=lflux
        endfor
        ind=WHERE(lambda LT 2.01e4 or lambda GT 2.43e4)
        loutspec[ind]=0.0
        loutvar[ind]=0.0
        loutsky[ind]=0.0
        
        loutarr=[[loutspec],[loutvar],[loutsky]]
        
        mkhdr,headout,loutarr
        sxaddpar,headout,'CRVAL1',lambda0
        sxaddpar,headout,'CRPIX1',1
        sxaddpar,headout,'CD1_1',dlambda
        SXADDPAR,head,'BUNIT','erg/cm2/s/A'
        WRITEFITS,'./intspec/'+outfile,loutarr,headout
     endfor
  
     close,1
;     stop
;  endfor
  
;  close,1
  intspec_wallace
  stop
  readcol,'./vor_out/intspec_wallace_goodmay.dat',filename,rin,rout,rav,sn,chi,vel,velmc,velerr,disp,dispmc,disperr,h3,h3mc,h3err,h4,h4mc,h4err,spclass,spclassmc,spclasserr,lumclass,lumclassmc,lumclasserr,nodispchi,FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  a=where(strmatch(filename, '*up*', /FOLD_CASE) eq 1)
  upvel=vel[a]
  b=where(strmatch(filename, '*down*', /FOLD_CASE) eq 1)
  downvel=vel[b]
  l=where((angle*!RADEG) lt 90.)
  deltavel=fltarr(N_ELEMENTS(angle))
  deltavel[l]=upvel[l]-downvel[l]
  m=where((angle*!RADEG) gt 90 AND (angle*!RADEG) lt 270)
  deltavel[m]=downvel[m]-upvel[m]
  u=where((angle*!RADEG) gt 270)
  deltavel[u]=upvel[u]-downvel[u]
  set_plot,'ps'
  device,filename='rotation.ps',/color
  djs_plot,angle*!RADEG,deltavel,psym=2,xran=[0.,362.],yran=[min(deltavel)-0.2,max(deltavel)+0.2],/xsty,/ysty,charsize=1.5,charthick=4,xthick=3,ythick=3,symsize=1.5,xtitle='Angle [ ^\circ ]';,ytitle='\Delta V [km/s] at '+STRTRIM(FIX(radii[j]),2)+' pixel(s)'
  velerr=sqrt(velerr[a]^2+velerr[b]^2)
  oploterr,angle*!RADEG,deltavel,velerr
  A=[2.,1.,1.]
  weights=fltarr(n_elements(angle))+1.
  result=curvefit(angle,deltavel,weights,A,sigma,FUNCTION_NAME='myfunct')
  djs_oplot,angle*!RADEG,result,color='blue'
  b=where(result eq min(result))
  print,result[b]
  print,angle[b]*!RADEG
  device,/close
  set_plot,'x'
  
  spawn, 'rm -rf ./intspec/*up*'
  spawn, 'rm -rf ./intspec/*down*'
  stop
  endfor
  

  
     
END
