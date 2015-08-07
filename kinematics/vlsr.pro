PRO VLSR
;http://www.dark-cosmology.dk/~justyn/iraf/obsdb.dat
;Gemini north
longitude=155.+28.142805d0/60.
latitude=19.+49.42809d0/60.
altitude = 4213.4
ra=190.48054167/15.
dec=11.66769444
;May20 2014
jd=56797.330592
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd
print,'May 20',corr
;May3 2015
jd=57144.402524
helcorr,longitude,latitude,altitude,ra,dec,jd,corr,hjd
print,'May 3',corr
;May 20      -22.627695
;May 3      -16.340266

STOP
END
