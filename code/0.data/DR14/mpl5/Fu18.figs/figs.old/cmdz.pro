;+
; Generate BPT and WHAN diagnostic diagrams 
;-

; load DRPALL
drpall = mrdfits('$MANGA_DIR/hfdap/$MANGADRP_VER/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)
; select main galxy sample - Primary, Color Enhanced, and Secondary
; maskbits MANGA_TARGET1    10  PRIMARY_v1_2_0 " "
; maskbits MANGA_TARGET1    11  SECONDARY_v1_2_0 " "
; maskbits MANGA_TARGET1    12  COLOR_ENHANCED_v1_2_0 " "
s = where((drpall.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
drpall = drpall[s]

Mi_cut = -18.0 
; pair catalog
pcat = mrdfits('../8.class.fits',1) ; 174
Rin = (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 166
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 130
	pcat.sep_arcsec le Rin and $ ; 121
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 100
pcat = pcat[s]

; indices for BAGN
s_bagns = where(strtrim(pcat.class) eq 'BAGN') ; 17
s_other = where(strtrim(pcat.class) ne 'BAGN') ; 149

; setup plot
loadct,0
psfile = 'figs/cmdz.eps'
setps,psfile,15*2.2,15,font='helvetica'
titsize = 1.5
labelsize = 1.2
syms = 1.2
pos1=[0.08,0.11,0.49-0.01/2,0.96]
pos2=[0.58+0.01/2,0.11,0.99,0.96]

; CMD diagram
xr=[-17.5,-25]
yr=[1,7]
xtit='Absolute Magnitude in Rest-Frame i-Band'
plot,[0,0],/nodata,/noerase,psym=1,xr=xr,yr=yr,xs=1,ys=1,$
	xtit=xtit,ytit='Rest-Frame NUV-r',$
	xminor=4,xtickinterval=2,$
	yminor=4,ytickinterval=2,$
	pos=pos1,charsize=titsize
; plot all sources
; Mcat = ABSMAG - 5 log(h) = mag - 5 log(cz/(100 km/s/Mpc * 10pc))
; where h = H0/100 = 1.0 in the catalog
; convert to h = 0.7
; ABSMAG = Mcat + 5 log(0.7) = Mcat - 0.775
; [FNugriz]
xall = drpall.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yall = drpall.NSA_ELPETRO_ABSMAG[1]-drpall.NSA_ELPETRO_ABSMAG[4]
plots,xall,yall,psym=cgSymCat(14),color=cgcolor('grey'),noclip=0,syms=0.5
ver,-21.5,lines=2
; plot individual sources
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = pcat.NSA_ELPETRO_ABSMAG[1]-pcat.NSA_ELPETRO_ABSMAG[4]
plots,xx[s_other],yy[s_other],psym=cgSymCat(9),noclip=0,syms=syms
plots,xx[s_bagns],yy[s_bagns],psym=cgSymCat(16),noclip=0,syms=syms,color=cgcolor('red')
plots,xx[s_bagns],yy[s_bagns],psym=cgSymCat(9),noclip=0,syms=syms
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
al_legend,['DR14 MaNGA','Merging Pairs','Binary AGNs'],$
	psym=[14,9,16],color=cgcolor(['gray','black','red']),$
	/top,/left,background=cgcolor('white')

; z vs. M_i diagram
yr = [0.0,0.17]
plot,[0,0],/nodata,/noerase,psym=1,xr=xr,yr=yr,/xs,/ys,$
	xtit=xtit,$
	ytit='Redshift',$
	xminor=4,xtickinterval=2,$
	pos=pos2,charsize=titsize
; plot all sources
xall = drpall.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yall = drpall.NSA_Z
plots,xall,yall,psym=cgSymCat(14),color=cgcolor('grey'),noclip=0,syms=0.5
; plot individual sources
xx = pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = pcat.NSA_Z
plots,xx[s_other],yy[s_other],psym=cgSymCat(9),noclip=0,syms=syms
plots,xx[s_bagns],yy[s_bagns],psym=cgSymCat(16),noclip=0,syms=syms,color=cgcolor('red')
plots,xx[s_bagns],yy[s_bagns],psym=cgSymCat(9),noclip=0,syms=syms
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')

device,/close
spawn,'gv '+psfile+' &'

end
