;+
; Generate diagnostic diagrams for BAGN candidates.
; This version makes images of different IFUs the same size
; It also plots the polygon-extracted spectra' measurements
;-

; load pair sample catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu)
; choose a subclass
Mi_cut = -18.0
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin and $ ; 105
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 105
pcat = pcat[s]

; sort in Mi
;idx = sort(abs(pcat.NSA_ELPETRO_ABSMAG[5]))
; sort in Z
;idx = sort(pcat.NSA_Z)
; sort in Sep
idx = sort(pcat.sep_kpc)
pcat = pcat[idx]

for jj=0,1 do begin

; set up plotting window
nx = 10
ny = 11
psfile = 'figs/pairimg'+strc(jj)+'.eps'
setps,psfile,8.5*nx,8.5*ny,font='helvetica'
multiplot,/default
multiplot,[nx,ny],/nomargin

blank = replicate(' ',10)
for i=0,n_elements(pcat)-1 do begin
	drp = pcat[i]		
	; SDSS finder chart
	; matches the values used by sdss.pro
	scale = 0.25 
	npix = fix(40./scale)
	xr = [-0.5,0.5]*scale*npix
	; SDSS image w/o spec label
	jpgfile = '~/work/manga/cats/sdss/'+drp.plateifu+'b.png'
	read_jpeg,jpgfile,sdss_img
	; make WCS header
	mkhdr,h_sdss,fltarr(npix,npix)
	putast,h_sdss,[[-1,0],[0,1]]*scale/3600d, [npix/2+0.5,npix/2-1.5], $
		[drp.ifura,drp.ifudec], ['RA---TAN','DEC--TAN']	
	; compute position
	position=!p.position
	r_ypos = position[3]-position[1]
	pos=position
	; display image
	tvimage,sdss_img,true=1,position=pos
	; invisible axes
	plot,[1,1],[1,1],/nodata,xr=[0,npix],yr=[0,npix],xs=1+4,ys=1+4,pos=pos 
	;; show NSA ELPETRO Ellipse
	;adxy,h_sdss,pcat[i].objra,pcat[i].objdec,xc,yc
	;;plots,xc,yc,psym=1,color=cgcolor('red'),syms=5,thick=2
	;rmax = pcat[i].NSA_ELPETRO_TH50_R/scale ; Re - arcsec -> pixel
	;rmin = rmax*pcat[i].NSA_ELPETRO_BA ; b/a ratio
	;pos_ang = pcat[i].NSA_ELPETRO_PHI+90 ; PA (east of N)
	;tvellipse,rmax,rmin,xc,yc,pos_ang,/data,color=cgcolor('red'),lines=2,noclip=0
	; show SDSS components
	if jj eq 1 then begin
	adxy,h_sdss,pcat[i].ras,pcat[i].decs,xc,yc
	plots,xc[0]+0.5,yc[0]+0.5,psym=cgsymcat(6),syms=3,color=cgcolor('green')
	plots,xc[1]+0.5,yc[1]+0.5,psym=cgsymcat(9),syms=3,color=cgcolor('green')
	endif
	; plot arcsec-interval tickmarks
	plot,[1,1],[1,1],/nodata,xr=xr,yr=xr,/xs,/ys,/noerase,$
		xtickint=10,xminor=5,ytickint=10,yminor=5,$		
		position=pos,xtickname=blank,color=cgcolor('white')
	; draw IFU around target
	; 3.5" increased due to 3-point dithering
	a = findgen(7)*(!pi*2/6.)
	oplot,(drp.IFUdiam+3.5)*0.5*cos(a),(drp.IFUdiam+3.5)*0.5*sin(a),$
		color=cgcolor('sky blue'),lines=2
	; label the plot
	xyouts,0,xr[1]*0.82,drp.plateifu,$ ; +' z='+string(drp.nsa_z,f='(f5.3)'),$
		chars=1.7,color=cgcolor('white'),align=0.5	
	;xyouts,xr[0]*0.9,xr[0]*0.93,string(pcat[i].sep_arcsec,f='(f4.1)')+$
	;	'"='+string(pcat[i].sep_kpc,f='(f4.1,"kpc")'),$
	;	chars=1.5,color=cgcolor('white'),align=0
	;; plot black outline
	;plot,[0,1],[0,1],/nodata,/xs,/ys,/noerase,$
	;	xtickint=10,xminor=-1,ytickint=10,yminor=-1,$		
	;	position=pos,xtickname=blank,color=cgcolor('black')
	multiplot

endfor
device,/close
;spawn,'gv '+psfile+' &'

endfor

end
