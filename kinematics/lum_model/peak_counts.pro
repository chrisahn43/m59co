PRO PEAK_COUNTS
READCOL,'Catlist',filenames,format='A'
FOR i=0,N_ELEMENTS(filenames)-1 DO BEGIN
   im=READFITS(filenames[i],/SILENT)
   print,filenames[i],MAX(im);,TOTAL(im)
ENDFOR
STOP
END
