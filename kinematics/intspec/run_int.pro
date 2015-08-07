@int_spectra_annulus.pro
PRO RUN_INT


fwhmfile='../../sky/m59co/fullcube_disp_med.fits'
contfile='../lum_model/m59co_combine_goodmay_cont.fits'
fwhm=READFITS(fwhmfile)
cont=READFITS(contfile)

radii=[0,1.5,2.5,4,6,9,12,18]
nradii=N_ELEMENTS(radii)

xcen=38.06
ycen=45.44
makex,cont,xarr,yarr,/zero
rarr=SQRT((xarr-xcen)^2+(yarr-ycen)^2)
OPENW,1,'spectra_fwhm_goodmay.dat'
FOR i=0,nradii-2 DO BEGIN
   ind=WHERE(rarr GT radii[i] AND rarr LE radii[i+1])   
   avfwhm=TOTAL(fwhm[ind]*cont[ind])/TOTAL(cont[ind])
   avrad=TOTAL(rarr[ind]*cont[ind])/TOTAL(cont[ind])
   print,avrad
   outfile='m59co_goodmay_int_fc'+STRTRIM(FIX(radii[i]),2)+'-'+STRTRIM(FIX(radii[i+1]),2)+'.fits'

   printf,1,radii[i],radii[i+1],avrad,avfwhm,' ',outfile,FORMAT='(2F5.1,2F7.3,2A)'
ENDFOR
close,1



FOR i=0,nradii-2 DO int_spectra_annulus,radii[i],radii[i+1],xcen,ycen

STOP
END
