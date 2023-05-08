;+
; Plot redshift, color, and LOIII vs. Mi for full DR14 main galaxy
; sample. Highlight pairs, SAGN, and BAGN subsamples.
;-

PRO mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,_REF_EXTRA=extra
	syms=1.0
	plot,[0,1],[0,1],/nodata,/xs,/ys,_extra=extra
	; full DR14 sample
	s = where((t1 and 2L^10+2L^12) ne 0)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('dark gray'),noclip=0,syms=0.5
	s = where((t1 and 2L^11) ne 0)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('gray'),noclip=0,syms=0.5
	; AGN sample
	;s = where((tt1 and 2L^10+2L^12) ne 0) ; primary+ - black filled
	;plots,xx[s],yy[s],psym=cgsymcat(16),noclip=0,syms=syms,color=cgcolor('dark gray')
	;s = where((tt1 and 2L^11) ne 0) ; secondary - grey filled
	;plots,xx[s],yy[s],psym=cgsymcat(16),noclip=0,syms=syms,color=cgcolor('gray')
	plots,xx,yy,noclip=0,psym=cgsymcat(9),syms=syms,color=cgcolor('royal blue')
	; AGN subsample (Kewley01 selection)
	plots,xx[idx1],yy[idx1],noclip=0,psym=cgsymcat(9),syms=syms,color=cgcolor('red')
END

; setup plot
loadct,0
psfile = 'figs/allvsmi_dr14.eps'
setps,psfile,15*1.1*4,15,font='helvetica'
titsize = 1.5
labelsize = 1.5
multiplot,/default
multiplot,[4,1],xgap=0.018,mxtitle=textoidl('Absolute Magnitude (M_i)'),$
	mxTitSize=titsize,mxTitOffset=1.3 

;;;;;;;;;;;;;;;;;;
; Load catalogs
;;;;;;;;;;;;;;;;;;

; two matched DR14 catalog - ncat, targDR14
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0)
ncat = ncat[s]
; match target catalog to DR14 catalog
; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]

;;;;;;;;;;;;;;;;;;
; Making Plots
;;;;;;;;;;;;;;;;;;
xr = [-17,-25]
xtit='' ;textoidl('Absolute Magnitude (M_i)')

; plot z vs. Mi 
yr = [0,0.17] 
ytit='Redshift'
; DR14 sample
x = targDR14.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = targDR14.NSA_z
t1 = targDR14.manga_target1 
; AGN sample
s = where(ncat.bclass ge 2)
xx = targDR14[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
yy = targDR14[s].NSA_z
tt1= targDR14[s].manga_target1
; Kew01 AGN sample
idx1 = where(ncat[s].bclass ge 3)
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,xminor=4,xtickinterval=2,charsize=titsize
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
al_legend,['Primary+ (DR14)','Secondary (DR14)','AGNs (Kauffmann03)','AGNs (Kewley01)'],$
	psym=[16,16,9,9],$
	color=cgcolor(['dark gray','gray','royal blue','red']),$
	/top,/left,background=cgcolor('white'),charsize=labelsize
multiplot,/dox,/doy

; plot NUV-r vs. Mi 
yr = [1, 7] 
ytit=textoidl('NUV-r')
; DR14 sample
y = targDR14.NSA_ELPETRO_ABSMAG[1]-targDR14.NSA_ELPETRO_ABSMAG[4]
; AGN sample
yy = targDR14[s].NSA_ELPETRO_ABSMAG[1]-targDR14[s].NSA_ELPETRO_ABSMAG[4]
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,$
	xminor=4,xtickinterval=2,yminor=5,charsize=titsize
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; plot LOIII vs. Mi
yr = [37.5,42.5]
ytit = textoidl('log(L_{[O III]} / erg s^{-1})')
; DR14 sample
y = alog10(ncat.Lo3) 
; AGN sample
yy = alog10(ncat[s].Lo3)
; plot symbols
mkplot,x,y,t1,xx,yy,tt1,idx1,idx2,$
	xr=xr,yr=yr,xtit=xtit,ytit=ytit,$
	xminor=4,xtickinterval=2,yminor=5,charsize=titsize
; legend
al_legend,'(c)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

; plot AGN fraction vs. Mi
x0 = targDR14.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y0 = targDR14.ESWEIGHT
s = where(ncat.bclass ge 2)
x1 = targDR14[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y1 = targDR14[s].ESWEIGHT
s = where(ncat.bclass ge 3)
x2 = targDR14[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y2 = targDR14[s].ESWEIGHT
; binning
bin = 1.0
np = 6
mi = range(-18.5,-23.5,np)
phi0 = mi*0
phi1 = mi*0
phi2 = mi*0
;nboot = 500
;bs = bootstrap_mean(y0[s],nboot=nboot)
for i = 0,n_elements(mi)-1 do begin
	s = between(x0,mi[i]-bin/2,mi[i]+bin/2)
	phi0[i] = total(y0[s])
	s = between(x1,mi[i]-bin/2,mi[i]+bin/2)
	phi1[i] = total(y1[s])
	s = between(x2,mi[i]-bin/2,mi[i]+bin/2)
	phi2[i] = total(y2[s])
endfor
;; add edge datapoints
;mi = [mi[0]+bin/2,mi,mi[np-1]-bin/2]
;phi0 = [phi0[0],phi0,phi0[np-1]]
;phi1 = [phi1[0],phi1,phi1[np-1]]
;phi2 = [phi2[0],phi2,phi2[np-1]]
; make plot
yr = [0,1]*0.25
plot,[0,1],[0,1],/nodata,/xs,/ys,xr=xr,yr=yr,ytickinterval=0.1,yminor=5,$
	xtit=xtit,ytit='AGN Fraction',charsize=titsize
xh = fltarr(np*2)
xh[0:*:2] = mi+bin/2
xh[1:*:2] = mi-bin/2
yh = fltarr(np*2)
yh[0:*:2] = phi1/phi0
yh[1:*:2] = phi1/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('royal blue')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/line_fill,orient=45,/data,color=cgcolor('royal blue')
yh = fltarr(np*2)
yh[0:*:2] = phi2/phi0
yh[1:*:2] = phi2/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('red')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/line_fill,orient=-45,/data,color=cgcolor('red')
; legend
al_legend,'(d)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')

device,/close
;spawn,'gv '+psfile+' &'

end
