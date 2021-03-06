;######################################################################
;
; Copyright (C) 1999-2005, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/software
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to
; `use of the MGE fitting software developed by Cappellari (2002)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;       TEST_MGE_FIT
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Exercise all the routines in the MGE_FIT_SECTORS package.
;       This procedure is intended to be used as a template to be
;       customized for each particular MGE fitting problem.
;
; EXPLANATION:
;       Further information is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       TEST_MGE_FIT
;
; EXAMPLE:
;       The following command will take various minutes to complete,
;       while plotting intermediate results on the screen
;
;           test_mge_fit
;
; PROCEDURES USED:
;       The following procedures are contained in the TEST_MGE_FIT program.
;           FIT_M32     -- Perform an MGE fit of the galaxy M32
;           FIT_NGC4342 -- Perform an MGE fit of the galaxy NGC 4342
;           FIT_NGC4473 -- Perform an MGE fit of the galaxy NGC 4473
;           FIT_DOUBLE_POWERLAW_1D -- Fit a 1D MGE model
;           MGE_FIT_1D_HERNQUIST_MODEL -- Compute the circular velocity of an MGE 1D fit
;           FIT_NGC5831_TWIST -- Perform a twisted MGE fit of the galaxy NGC 5831
;
;       Other routines needed from the MGE_FIT_SECTORS package by
;       Michele Cappellari (http://purl.org/cappellari/idl)
;           SECTORS_PHOTOMETRY       -- perform photometry along sectors
;           MGE_FIT_SECTORS          -- do the actual MGE fit
;           MGE_PRINT_CONTOURS       -- plot the results
;           SECTORS_PHOTOMETRY_TWIST -- perform photometry along sectors with point symmetry
;           MGE_FIT_SECTORS_TWIST    -- do the actual MGE fit with possible isophote twist
;           MGE_PRINT_CONTOURS_TWIST -- plot the results with possible isophote twist
;           WFPC2_MGE_FIT            -- do an MGE fit of a WFPC2 image (PC1 + MOSAIC)
;           FIND_GALAXY              -- find galaxy center and position angle
;           MGE_FIT_1D               -- perform a 1D MGE fit
;           CAP_RANGE                -- returns a sequence of values
;
;       Other IDL routines needed:
;           BVLS  -- Michele Cappellari porting of Lawson & Hanson generalized NNLS
;                    http://purl.org/cappellari/idl
;           MPFIT -- Craig Markwardt porting of Levenberg-Marquardt MINPACK-1
;                    http://purl.com/net/mpfit
;           JAM package -- The Jeans Anisotropic MGE (JAM) package is required only
;                    if one wants to run the optional example MGE_FIT_1D_HERNQUIST_MODEL
;                    http://purl.org/cappellari/idl
;
;       Astronomy User's Library routines needed (http://idlastro.gsfc.nasa.gov):
;           FITS_READ
;           DIST_CIRCLE
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Leiden, January 2000
;       V2.0: Updated all examples, MC, Leiden July 2001
;       V2.1: Added 1D MGE fit example, MC, Leiden, 10 May 2003
;       V2.11: Removed the un-necessary parameter EPS from the routine
;           SECTORS_PHOTOMETRY_TWIST. MC, Leiden, 1 September 2004
;       V2.12: Replaced LOGRANGE keyword with the new MAGRANGE.
;           MC, Leiden, 1 May 2005
;       V2.2: Included a new example MGE_FIT_1D_HERNQUIST_MODEL.
;           MC, Oxford, 30 November 2008
;       V2.21: Uses CAP_RANGE. MC, Paranal, 8 November 2013
;-
;----------------------------------------------------------------------------
PRO fit_m59co
;
fits_read, '../data/m59_g_lucy.fits', img, h
;estimate sky level
;mdrizz=[7.46042280032,7.9962758789,6.9261413824]
expt=375.0
;mdrizzcount=mdrizz/expt
;avgdrizz=mean(mdrizzcount)
img=img/expt

;bl=img[265:365,23:123]
;br=img[1182:1282,125:225]
;tl=img[25:125,1070:1170]
;tr=img[936:1036,1176:1276]
;avgsky=[mean(bl),mean(br),mean(tl),mean(tr)]
;skylev = mean(avgsky) ; counts/pixel
;img = img - skylev    ; subtract sky
imgsize=size(img,/dim)
scale = 0.05
ngauss = 15
minlevel = 0.000;01 ; counts/pixel

readcol,'tinytim_fits.dat',normPSF,sigmaPSF,q,format='F,F,F'

; Here we use FIND_GALAXY directly inside the procedure. Usually you may want
; to experiment with different values of the FRACTION keyword, before adopting
; given values of Eps, Ang, Xc, Yc.

find_galaxy, img, majorAxis, eps, ang, xc, yc,FRACTION=1. ,/PLOT

; Perform galaxy photometry

sectors_photometry, img, eps, ang, xc, yc, radius, angle, counts, MINLEVEL=minlevel

; Do the actual MGE fit 
;set_plot,'ps'
;device,filename='m59co_photfits.ps'
MGE_fit_sectors, radius, angle, counts, eps, $
    NGAUSS=ngauss, SIGMAPSF=sigmaPSF, NORMPSF=normPSF, SOL=sol, SCALE=scale
;device,/close
;set_plot,'x'
stop

peakbright=(sol[0,*])/(2*!PI*((sol[1,*])^2)*sol[2,*])
extinct=0.109
surfbrightI=26.16252+5.*alog10(scale)-2.5*alog10(peakbright)-extinct
;http://www.stsci.edu/hst/acs/analysis/zeropoints/old_page/localZeropoints
Intensity=((64800./!PI)^2)*10^(0.4*(4.08-surfbrightI))
print,intensity
sigmaarc=(sol[1,*])*scale
q=sol[2,*]
forprint,intensity,sigmaarc,q,format='F,F,F',textout='m59co_mge_output.dat'

modelweight=peakbright
modelsig=sol[1,*]
rad=10^((findgen(36)*0.073)+0.073);[findgen(19)*3+3,findgen(30)*10+60];,250,370]
minrad=10^(findgen(36)*0.073);[findgen(20)*3,findgen(29)*10+60];,200,250]
N=n_elements(rad)
photflux=fltarr(N)
modelflux=fltarr(N)
area=fltarr(N)
photmag=fltarr(N)
modelmag=fltarr(N)

zp=26.05923
modelimg=fltarr(imgsize[0],imgsize[1])
x=[reverse(findgen(xc))+1,findgen(imgsize[0]-xc)]
y=[reverse(findgen(yc))+1,findgen(imgsize[1]-yc)]
for i=0,n_elements(x)-1 do begin
   for j=0,n_elements(y)-1 do begin
      temp=0.
      for k=0,n_elements(modelweight)-1 do begin
         temp+= (modelweight[k]*exp(-(1./(2*modelsig[k]^2))*(x[i]^2+(y[j]^2/q[k]^2))))
      endfor
      modelimg[i,j]=temp
    endfor
endfor
writefits,'hst_mge_m59co.fits',modelimg

for i=0,n_elements(rad)-1 do begin
   maxflux_phot=djs_phot(xc,yc,rad[i],0.,img,skyval=skyval)
   minflux_phot=djs_phot(xc,yc,minrad[i],0.,img,skyval=skyval)
   maxflux_model=djs_phot(xc,yc,rad[i],0.,modelimg,skyval=skyval)
   minflux_model=djs_phot(xc,yc,minrad[i],0.,modelimg,skyval=skyval)
   area[i]=((!PI*(rad[i])^2)-(!PI*(minrad[i])^2))*scale
   photflux[i]=maxflux_phot-minflux_phot
   modelflux[i]=maxflux_model-minflux_model
   photmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(photflux[i])
   modelmag[i]=zp+5*alog10(scale)+2.5*alog10(area[i])-2.5*alog10(modelflux[i])
endfor
;set_plot,'ps'
;device,filename='surfbright_hst_fits.ps',/color
djs_plot,rad*scale,modelmag,psym=2,xtitle='Radius ["]',ytitle='\mu [Mag/sqare arcsecond]',yran=[23,10],xran=[0,13],charsize=1.5,charthick=4,xthick=3,ythick=3
djs_oplot,rad*scale,photmag,color='blue',thick=3;,psym=2
djs_oplot,rad*scale,photmag,psym=4
;device,/close
;set_plot,'x'

stop
; Print the data-model contours comparison of the whole image

MGE_print_contours, img>minlevel, ang, xc, yc, sol, $
    FILE='hst_m59co.ps', SCALE=scale, MAGRANGE=9, $
    SIGMAPSF=sigmaPSF, NORMPSF=normPSF, BINNING=7

; Print the data-model contours comparison of the central regions

s = SIZE(img)
img = img[xc-s[1]/9:xc+s[1]/9,yc-s[2]/9:yc+s[2]/9]
MGE_print_contours, img, ang, s[1]/9, s[2]/9, sol, $
    FILE='hst_m59co_nuclear.ps', SCALE=scale, MAGRANGE=9, $
                    SIGMAPSF=sigmaPSF, NORMPSF=normPSF


END
;----------------------------------------------------------------------------
;----------------------------------------------------------------------------
PRO hst_m59co_fit
;
; This is the main routine to call in succession the MGE fits to
; M32, NGC4342, NGC4473, power-law and NGC5831, and measure the execution time.
; A run of this program takes: 
; - 691s on a Pentium III, 1.0GHz PC, with IDL 5.4.
; - 270s on a Pentium M, 1.6GHz PC, with IDL 6.1.
; - 150s on a Core2 Duo, 2.2GHz PC, with IDL 7.0.
; It was tested with IDL 5.6-8.1 under both Windows and Linux.
;
t = SYSTIME(1)
fit_m59co
PRINT, 'Total computation time:', SYSTIME(1) - t, ' Seconds'

END
;----------------------------------------------------------------------------
