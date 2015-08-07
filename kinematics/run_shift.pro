@shift_cube.pro
PRO RUN_SHIFT

reffile='CatfbrgnN20140518S0045.fits'

cube = mrdfits(reffile,1,h1)

imsize=SIZE(cube,/dim)
nlambda=imsize[2]
lambda0=SXPAR(h1,'CRVAL3')
dlambda=SXPAR(h1,'CD3_3')
lambda=FINDGEN(imsize[2])*dlambda+lambda0
;Feb 20 16.783
;May 18 -22.006
;correct everything to May 18 data
corfiles=['CatfbrgnN20140220S0279','CatfbrgnN20140220S0281','CatfbrgnN20140220S0282','CatfbrgnN20140220S0284']
;veldiff=vbarycor_Object-vbarycor_Reference 
veldiff=16.783+22.006
FOR i=0,N_ELEMENTS(corfiles)-1 DO SHIFT_CUBE,corfiles[i],lambda,VELDIFF=veldiff
END
