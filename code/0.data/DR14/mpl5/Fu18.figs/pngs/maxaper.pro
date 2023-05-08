
; load all DR14 catalog
ncat = mrdfits('../7.nucsample.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]
s_pri = where((ncat.mngtarg1 and 2L^10) ne 0) ; primary
s_sec = where((ncat.mngtarg1 and 2L^11) ne 0) ; secondary
s_coe = where((ncat.mngtarg1 and 2L^12) ne 0) ; color enhanced

; load galaxy pair catalog
pcat = mrdfits('../7.pairsample.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 168
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 128
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 122
	pcat.sep_arcsec le pcat.IFURin)  ; 114
pcat = pcat[s]

;;;;;;;;;;;;;;;
; Distribution of maximum allowed aperture radius in kpc
;;;;;;;;;;;;;;;
; angular scale
scale = dangular(ncat.NSA_Z,/kpc)/(180.*3600./!pi) ; kpc/arcsec
print,minmax(scale) ; 0.22 - 2.62 kpc/arcsec
; offset between IFU center and Target center
offset = sphdist(ncat.OBJRA,ncat.OBJDEC,ncat.IFURA,ncat.IFUDEC,/deg)*3600
help,where(offset gt 1.0) ; 33
; maximum aperture radius in kpc  
maxkpc = scale * (ncat.IFURin - offset)
plothist,maxkpc,xtit='Maximum Aperture Radius Allowed by IFU Size (kpc)'
save_screen,'maxaper1.png'
help,where(maxkpc lt max(scale)/2) ; 1 = 8239-3701 where offset = 11.7" > IFU Rin = 8.1"
help,where(maxkpc lt max(scale)) ; 58 out of 2618 

; pair separation distribution
help,pcat,where(pcat.sep_kpc lt max(scale)) ; 9 out of 114
print,minmax(pcat.sep_kpc) ; 1.2 - 28.6 kpc

;;;;;;;;;;;;;;;;
plothist,scale*ncat.IFURin,xtit='Maximum Pair Separation Allowed by IFU Size (kpc)'
plothist,scale*ncat[s_pri].IFURin,/overplot,color=cgcolor('red')
plothist,scale*ncat[s_sec].IFURin,/overplot,color=cgcolor('blue')
plothist,scale*ncat[s_coe].IFURin,/overplot,color=cgcolor('green')
save_screen,'maxaper2.png'

;;;;;;;;;;;;;;;;
bin = 0.005
plothist,ncat.NSA_z,xtit='Redshift',bin=bin
plothist,ncat[s_pri].NSA_z,bin=bin,/overplot,color=cgcolor('red')
plothist,ncat[s_sec].NSA_z,bin=bin,/overplot,color=cgcolor('blue')
plothist,ncat[s_coe].NSA_z,bin=bin,/overplot,color=cgcolor('green')
save_screen,'maxaper3.png'

; effective radius vs. stellar mass
plot,ncat.NSA_ELPETRO_ABSMAG[5], ncat.NSA_ELPETRO_TH50_R * scale,$
	psym=6,yr=[0,50],xr=[-17,-24],$
	xtit='Absolute Magnitude',ytit='Effective Radius (kpc)'
plot,alog10(ncat.NSA_ELPETRO_MASS), ncat.NSA_ELPETRO_TH50_R * scale,$
	psym=6,yr=[0,50],xr=[8,12],$
	xtit='Stellar Mass',ytit='Effective Radius (kpc)'
save_screen,'maxaper4.png'

; Blandon 07 Fig 18
; i-band mass-to-light ratio
; M_i - M_bol(sun) = -2.5 log(L/Lsun)
Mbsun = 4.74
MtoL = alog10(ncat.NSA_ELPETRO_MASS)+0.4*(ncat.NSA_ELPETRO_ABSMAG[4]-Mbsun)
color = ncat.NSA_ELPETRO_ABSMAG[3] - ncat.NSA_ELPETRO_ABSMAG[4] 
plot,color,MtoL,psym=3,xr=[0.1,1.1],yr=[-0.4,0.5],/xs,/ys,xtit='g-r',ytit='log(M/L_r [Msun/Lsun])'
x = range(0.1,1.1,50)
oplot,x,1.4*x - 0.73,lines=2 ; Bell & de Jong (2001) curve
save_screen,'maxaper4.png'

end
