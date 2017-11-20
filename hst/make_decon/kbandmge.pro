PRO KBANDMGE
  scale=0.05
  zeropoint_g_ab=26.05923
  zeropoint_z_ab=24.84245
  zeropoint_g_ve=26.16252
  zeropoint_z_ve=24.32305
  colorin=1.42-(0.107-0.041)
  colorout=1.6-(0.107-0.041)
  gsun=5.10 ;MAY NEED TO BE VEGA MAG
  zsun=4.52
  readcol,'~/bc03/bc03/hst_add.1ABmag',logage,u,g,r,i,z,f606w,f814w,f475w,f850lp,format='F,F,F,F,F,F,F,F,F,F'
  hstcolor=f475w-f850lp
  readcol,'~/bc03/bc03/hst_add.4color',logage4,mbol,bmag,vmag,kmag,mliv,mrem,mret,mgal,sfr,mtot,mlbtot,mlvtot,format='F,F,F,F,F,F,F,F,F,F,F,F,F'
  inind=where(hstcolor gt colorin-0.001 and hstcolor lt colorin+0.001)
  outind=where(hstcolor gt colorout-0.02 and hstcolor lt colorout+0.02)
  ksun=3.28
  gk=f475w-kmag
  inscale=(10^(-0.4*(gk[inind]-(gsun-ksun))))
  inscale=inscale[0]
  outscale=(10^(-0.4*(gk[outind]-(gsun-ksun))))
  outscale=outscale[0]
  readcol,'../m59co_mge_outputsersic.dat',intensity,sigma,q,pa,format='D,D,D,D'
  a=where(pa gt 20.)
  intensity[a]=intensity[a]*inscale
  b=where(pa lt 20.)
  intensity[b]=intensity[b]*outscale
  forprint,intensity,sigma,q,pa,textout='./m59co_mge_outputsersic_k.dat',format='D,D,D,D'
  stop

END
