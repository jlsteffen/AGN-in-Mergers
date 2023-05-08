; plot AGN fraction vs. Mi
xr = [-17,-25]
xtit = textoidl('Absolute Magnitude (M_i)')
x0 = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y0 = ntarg.ESWEIGHT
s = where(ncat.bclass eq 2)
x1 = ntarg[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y1 = ntarg[s].ESWEIGHT
s = where(ncat.bclass eq 3)
x2 = ntarg[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y2 = ntarg[s].ESWEIGHT
s = where(ncat.bclass eq 4)
x3 = ntarg[s].NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y3 = ntarg[s].ESWEIGHT
; binning
bin = 1.0
np = 6
mi = range(-18.5,-23.5,np)
phi0 = mi*0
phi1 = mi*0
phi2 = mi*0
phi3 = mi*0
for i = 0,n_elements(mi)-1 do begin
	s = between(x0,mi[i]-bin/2,mi[i]+bin/2)
	phi0[i] = total(y0[s])
	s = between(x1,mi[i]-bin/2,mi[i]+bin/2)
	phi1[i] = total(y1[s])
	s = between(x2,mi[i]-bin/2,mi[i]+bin/2)
	phi2[i] = total(y2[s])
	s = between(x3,mi[i]-bin/2,mi[i]+bin/2)
	phi3[i] = total(y3[s])
endfor
; make plot
yr = [0,1]*0.25
plot,[0,1],[0,1],/nodata,/xs,/ys,xr=xr,yr=yr,ytickinterval=0.1,yminor=5,$
	xtit=xtit,ytit='AGN Fraction',charsize=titsize
; plot histograms
xh = fltarr(np*2)
xh[0:*:2] = mi+bin/2
xh[1:*:2] = mi-bin/2
yh = fltarr(np*2)
yh[0:*:2] = (phi1+phi2+phi3)/phi0
yh[1:*:2] = (phi1+phi2+phi3)/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('black')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/fill,/data,color=cgcolor('light gray')
yh = fltarr(np*2)
yh[0:*:2] = phi1/phi0
yh[1:*:2] = phi1/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('royal blue')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/line_fill,/data,orient=0,color=cgcolor('royal blue')
yh = fltarr(np*2)
yh[0:*:2] = phi2/phi0
yh[1:*:2] = phi2/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('forest green')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/line_fill,orient=-45,/data,color=cgcolor('forest green')
yh = fltarr(np*2)
yh[0:*:2] = phi3/phi0
yh[1:*:2] = phi3/phi0
oplot,[xh[0],xh,xh[np*2-1]],[0,yh,0],lines=0,color=cgcolor('red')
polyfill,[xh[0],xh,xh[np*2-1]],[0,yh,0],/line_fill,orient=45,/data,color=cgcolor('red')
; legend
;ver,-21,lines=2
al_legend,'(f)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')

end
