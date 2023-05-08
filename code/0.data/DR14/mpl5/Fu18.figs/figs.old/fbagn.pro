;+
; Compute f_BAGN among mergers based on comoving volume weights
; Generate a figure as well
;-

; setup plot
loadct,0
psfile = 'figs/fbagn.eps'
setps,psfile,15*2.2,15,font='helvetica'
titsize = 1.5
labelsize = 1.2
syms = 1.5
pos1=[0.09,0.11,0.49+0.01/2,0.96]
pos2=[0.58+0.01/2,0.11,0.99,0.96]

Mi_cut = -21.5
; load pair catalog
pcat = mrdfits('../8.class.fits',1)
Rin = (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 166
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut and $ ; 135
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 107
	pcat.sep_arcsec le Rin) ; 100
pcat = pcat[s]
; Note that big offsets between IFU center and primary object is fine 
; as long as Pair Sep is less than Rin of the IFU, because these pairs
; would remain in the pair sample even if we recenter the IFU on the object

; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targ = target[sa]

; match target catalog to DR14 catalog
drpall = mrdfits('../../sample/matched/drpall.fits',1)
s = where((drpall.mngtarg1 and 2L^10+2L^11+2L^12) ne 0)
drpall = drpall[s]
match2,drpall.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]
Rin = (drpall.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius
maxsepkpc = Rin * dangular(drpall.nsa_z,/kpc)/206264.81d

; Primary / Secondary / Color-Enhanced in pair sample
help,where((pcat.mngtarg1 and 2L^10) ne 0) ; 48 (39.7%)
help,where((pcat.mngtarg1 and 2L^11) ne 0) ; 52 (43.0%)
help,where((pcat.mngtarg1 and 2L^12) ne 0) ; 21 (17.3%)
;print,minmax(pcat.NSA_ELPETRO_ABSMAG[5])   ; -23.256400      -17.512800
; in comparison to DR14, the pair sample contains more Secondarys
; 1245 primary (47.6%)
; 935 secondary (35.7%)
; 438 color enhanced (16.7%)
; Primary / Secondary / Color-Enhanced in BAGN sample
s = where(min(pcat.bclass,dim=1) ge 2)
help,where((pcat[s].mngtarg1 and 2L^10) ne 0) ;  2 (13.3%)
help,where((pcat[s].mngtarg1 and 2L^11) ne 0) ; 10 (66.7%)
help,where((pcat[s].mngtarg1 and 2L^12) ne 0) ;  3 (20.0%)
;print,minmax(pcat[s].NSA_ELPETRO_ABSMAG[5]) ; -22.831 -20.9194

; Comoving Volume Density of BAGN
s = where(strtrim(pcat.class) eq 'BAGN',n_bagn)
cvd_bagn = total(targ[s].ESWEIGHT)/1e6 ; #/Mpc^3 comoving
; Comoving Volume Density of Mergers
n_pair = n_elements(pcat)
cvd_pair = total(targ.ESWEIGHT)/1e6 ; #/Mpc^3 
; BAGN fraction - 5.3% vs. 10% from simple division of numbers 
;print,cvd_bagn,cvd_pair,cvd_bagn/cvd_pair,1.0*n_bagn/n_pair

yr = [0.05,100]

;; plot ESWEIGHT vs. projected angular separation
;plot,[0,1],[0,1],/nodata,xr=[-2,17],yr=yr,/yl,/xs,/ys,$
;	xtit='Pair Angular Seperation (arcsec)',$
;	ytit=textoidl('Density (# per 10^6 Mpc^{-3})'),$
;	xminor=5,xtickinterval=5,$
;	pos=pos1,charsize=titsize
;; full DR14 sample
;plots,Rin,targDR14.esweight,psym=3,color=cgcolor('black')
;; pair sample
;xx = pcat.sep_arcsec
;yy = targ.esweight
;s = where((pcat.mngtarg1 and 2L^10) ne 0) ; primary - black filled
;plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('black')
;s = where((pcat.mngtarg1 and 2L^12) ne 0) ; CE - black filled
;plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('grey')
;s = where((pcat.mngtarg1 and 2L^11) ne 0) ; secondary - white filled
;plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('white')
;; exterior circles
;plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('black')
;color = cgcolor('black')
;bin = 5. ; kpc bins
;for xbin=2.5,12.5,bin do begin
;	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
;	bs = bootstrap_mean(yy[s],nboot=500)
;	print,xbin,total(yy[s])
;	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
;	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
;endfor
;; BAGN subsample
;s = where(min(pcat.bclass,dim=1) ge 2)
;xx = pcat[s].sep_arcsec
;yy = targ[s].esweight
;color = cgcolor('red')
;plots,xx,yy,psym=cgsymcat(9),syms=syms,noclip=0,color=color
;for xbin=2.5,12.5,bin do begin
;	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
;	bs = bootstrap_mean(yy[s],nboot=500)
;	print,xbin,total(yy[s])
;	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
;	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
;endfor

; plot ESWEIGHT vs. Mi 
xr=[-17.5,-25]
plot,[0,1],[0,1],/nodata,xr=xr,yr=yr,/yl,/xs,/ys,$
	xtit='i-band Absolute Mag',$
	ytit=textoidl('Density (# per 10^6 Mpc^{-3})'),$
	xminor=4,xtickinterval=2,$
	pos=pos1,charsize=titsize
; full DR14 sample
s = where((targDR14.manga_target1 and 2L^10+2L^12) ne 0)
plots,targDR14[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7),targDR14[s].esweight,psym=cgSymCat(14),color=cgcolor('black'),noclip=0,syms=0.5
s = where((targDR14.manga_target1 and 2L^11) ne 0)
plots,targDR14[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7),targDR14[s].esweight,psym=cgSymCat(14),color=cgcolor('grey'),noclip=0,syms=0.5
; pair sample
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = targ.esweight
s = where((pcat.mngtarg1 and 2L^10) ne 0) ; primary - black filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('black')
s = where((pcat.mngtarg1 and 2L^12) ne 0) ; CE - black filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('grey')
s = where((pcat.mngtarg1 and 2L^11) ne 0) ; secondary - white filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('white')
; exterior circles
plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('black')
; calculate total weights
bin = 1.0 ; mag bins
color = cgcolor('black')
for xbin=-23.0,-22.0,bin do begin
	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
	bs = bootstrap_mean(yy[s],nboot=500)
	print,xbin,total(yy[s])
	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
endfor
; BAGN subsample
s = where(min(pcat.bclass,dim=1) ge 2)
xx = pcat[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = targ[s].esweight
color = cgcolor('red')
print,minmax(xx)
plots,xx,yy,psym=cgsymcat(9),syms=syms,noclip=0,color=color
for xbin=-23.0,-22.0,bin do begin
	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
	bs = bootstrap_mean(yy[s],nboot=500)
	print,xbin,total(yy[s]),ct
	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
endfor

;;;;;;;;;;;;;;;
; plot ESWEIGHT vs. projected proper separation
plot,[0,1],[0,1],/nodata,/noerase,xr=[-2,17]*2,yr=yr,/yl,/xs,ys=1,$
	xtit='Pair Seperation or IFU Radius (kpc)',$
	;ytit=textoidl('Comoving Volume Density (# per 10^6 Mpc^{-3})'),$
	xminor=5,xtickinterval=5,$
	pos=pos2,charsize=titsize
; full DR14 sample
s = where((targDR14.manga_target1 and 2L^10+2L^12) ne 0)
plots,maxsepkpc[s],targDR14[s].esweight,psym=cgSymCat(14),color=cgcolor('black'),noclip=0,syms=0.5
s = where((targDR14.manga_target1 and 2L^11) ne 0)
plots,maxsepkpc[s],targDR14[s].esweight,psym=cgSymCat(14),color=cgcolor('grey'),noclip=0,syms=0.5
; pair sample
xx = pcat.sep_kpc
yy = targ.esweight
s = where((pcat.mngtarg1 and 2L^10) ne 0) ; primary - black filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('black')
s = where((pcat.mngtarg1 and 2L^12) ne 0) ; CE - grey filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('grey')
s = where((pcat.mngtarg1 and 2L^11) ne 0) ; secondary - white filled
plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('white')
plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('black')
bin = 10. ; kpc bins
color = cgcolor('black')
for xbin=7.5,17.5,bin do begin
	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
	bs = bootstrap_mean(yy[s],nboot=500)
	print,xbin,total(yy[s])
	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
endfor
; BAGN subsample
s = where(min(pcat.bclass,dim=1) ge 2)
xx = pcat[s].sep_kpc
yy = targ[s].esweight
color = cgcolor('red')
plots,xx,yy,psym=cgsymcat(9),syms=syms,noclip=0,color=color
for xbin=7.5,17.5,bin do begin
	s = where(xx ge xbin-bin/2 and xx lt xbin+bin/2,ct)
	bs = bootstrap_mean(yy[s],nboot=500)
	print,xbin,total(yy[s]),ct
	oploterror,xbin,total(yy[s]),bin/2,(bs[2]-bs[1])*ct,/hibar,color=color
	oploterror,xbin,total(yy[s]),bin/2,(bs[1]-bs[0])*ct,/lobar,color=color
endfor



device,/close
spawn,'gv '+psfile+' &'

end
