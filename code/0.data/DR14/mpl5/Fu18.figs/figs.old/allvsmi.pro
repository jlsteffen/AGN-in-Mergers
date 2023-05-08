;+
; Plot redshift, color, and LOIII vs. Mi for full DR14 main galaxy
; sample. Highlight pairs, SAGN, and BAGN subsamples.
;-

PRO mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,_REF_EXTRA=extra
	syms=1.5
	plot,[0,1],[0,1],/nodata,/xs,/ys,_extra=extra
	; full DR14 sample
	s = where((t1 and 2L^10+2L^12) ne 0)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('dark gray'),noclip=0,syms=0.5
	s = where((t1 and 2L^11) ne 0)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('gray'),noclip=0,syms=0.5
	; pair sample
	s = where((tt1 and 2L^10+2L^12) ne 0) ; primary+ - black filled
	plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('black')
	s = where((tt1 and 2L^11) ne 0) ; secondary - grey filled
	plots,xx[s],yy[s],psym=cgsymcat(16),syms=syms,color=cgcolor('gray')
	plots,xx,yy,psym=cgsymcat(9),syms=syms,color=cgcolor('royal blue')
	; SAGN subsample
	plots,xx[idx1],yy[idx1],psym=cgsymcat(9),syms=syms,noclip=0,color=cgcolor('orange')
	; BAGN subsample
	plots,xx[idx2],yy[idx2],psym=cgsymcat(9),syms=syms,noclip=0,color=cgcolor('red')
END

; setup plot
loadct,0
psfile = 'figs/allvsmi.eps'
setps,psfile,15*1.1*3,15,font='helvetica'
titsize = 1.5
labelsize = 1.5
multiplot,/default
positon = [0.01,0.01,0.99,0.99]
multiplot,[3,1],xgap=0.025 ;,position=position

;;;;;;;;;;;;;;;;;;
; Load catalogs
;;;;;;;;;;;;;;;;;;
Mi_cut = -18.0

; three matched DR14 pair catalog - pcat,targ,ncatpair
pcat = mrdfits('../7.pairsample.fits',1)
Rin = (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 166
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 130
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 124
	pcat.sep_arcsec le Rin and $ ; 116
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 116 -> 95 (Mi_cut = -21.5)
pcat = pcat[s]
; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targ = target[sa]
; match DR14 drpall/nucsample to pair catalog
ncat = mrdfits('../7.nucsample.fits',1)
match2,pcat.NSA_NSAID,ncat.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ncatpair = ncat[sa]

; two matched DR14 catalog - ncat, targDR14
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0)
ncat = ncat[s]
; match target catalog to DR14 catalog
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]

;;;;;;;;;;;;;;;;;;
; Making Plots
;;;;;;;;;;;;;;;;;;

; plot z vs. Mi 
xr = [-17.5,-25]
yr = [0,0.17] 
xtit=textoidl('Absolute Magnitude (M_i)')
ytit='Redshift'
; DR14 sample
x = targDR14.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = targDR14.NSA_z
t1 = targDR14.manga_target1 
; pair sample
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = targ.NSA_z
tt1= targ.manga_target1
; single AGN sample
idx1 = where(max(pcat.bclass,dim=1) ge 2)
; BAGN sample
idx2 = where(min(pcat.bclass,dim=1) ge 2)
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,xminor=4,xtickinterval=2,charsize=titsize
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
al_legend,['Primary+ (DR14)','Secondary (DR14)','Binary AGNs','Single AGNs','Other Pairs'],$
	psym=[16,16,9,9,9],$
	color=cgcolor(['dark gray','gray','red','orange','royal blue']),$
	/top,/left,background=cgcolor('white'),charsize=labelsize
multiplot,/dox,/doy

; plot NUV-r vs. Mi 
yr = [1, 7] 
ytit=textoidl('NUV-r')
; DR14 sample
x = targDR14.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = targDR14.NSA_ELPETRO_ABSMAG[1]-targDR14.NSA_ELPETRO_ABSMAG[4]
t1 = targDR14.manga_target1 
; pair sample
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = pcat.NSA_ELPETRO_ABSMAG[1]-pcat.NSA_ELPETRO_ABSMAG[4]
tt1= targ.manga_target1
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,$
	xminor=4,xtickinterval=2,yminor=5,charsize=titsize
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; plot LOIII vs. Mi
yr = [36.5,42.5]
ytit = textoidl('log(L_{[O~III]} / erg s^{-1})')
; DR14 sample
x = targDR14.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = alog10(ncat.Lo3) 
t1 = targDR14.manga_target1 
; pair sample
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = alog10(ncatpair.Lo3)
tt1= targ.manga_target1
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,$
	xminor=4,xtickinterval=2,yminor=5,charsize=titsize
; legend
al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
;multiplot,/dox,/doy


device,/close
spawn,'gv '+psfile+' &'

end
