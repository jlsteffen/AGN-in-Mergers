
; load all DR14 catalog
ncat = mrdfits('../7.nucsample.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]

; load target catalog w/ new weights based on DR14
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to DR14 catalog
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]
help,where(targDR14.NSA_NSAID ne ncat.NSA_NSAID) ; -1
help,where(targDR14.MANGA_TARGET1 ne ncat.mngtarg1) ; 12 mismatches

; load pair catalog
pcat = mrdfits('../7.pairsample.fits',1) ; 174
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 166
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 130
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600) ; 124
pcat = pcat[s]

wind,0,xsize=1000,ysize=700
multiplot,/default
multiplot,[3,2],xgap=0.04,ygap=0.04 
sauron_colormap

for bc = 2,3 do begin
; ESWEIGHT vs. Mi
x = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = targDR14.ESWEIGHT/1e6
xtit = 'M_i'
ytit = 'ESWeight'
; AGN subsample
iAGN = where(ncat.bclass ge bc) ; 389
xx = x[iAGN]
yy = y[iAGN]
; plotting
plot,x,alog10(y),xr=[-17,-25],yr=[-7.4,-4.3],psym=3,/xs,/ys,xtit=xtit,ytit=ytit
for i=0,n_elements(iAGN)-1 do begin
	color = interpol([25,230],alog10(minmax(ncat[iAGN].lo3)),$
		alog10(ncat[iAGN[i]].lo3))
	plots,x[iAGN[i]],alog10(y[iAGN[i]]),psym=cgsymcat(9),color=color
endfor
multiplot,/dox,/doy

; L[OIII] vs. Mi
x = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = ncat.Lo3
xtit = 'M_i'
ytit = '[OIII] Luminosity'
; AGN subsample
iAGN = where(ncat.bclass ge bc) 
xx = x[iAGN]
yy = y[iAGN]
; plotting
plot,x,alog10(y),xr=[-17,-25],yr=[37,43],psym=3,/xs,/ys,xtit=xtit,ytit=ytit
hor,40,lines=2
ver,-21,lines=2
for i=0,n_elements(iAGN)-1 do begin
	color = interpol([25,230],alog10(minmax(ncat[iAGN].lo3)),$
		alog10(ncat[iAGN[i]].lo3))
	plots,x[iAGN[i]],alog10(y[iAGN[i]]),psym=cgsymcat(9),color=color
endfor
multiplot,/dox,/doy

; L[OIII] histogram
iAGN = where(ncat.bclass ge bc) 
x = ncat[iAGN].Lo3
idx1 = where(pcat.bclass ge bc) ; ALL AGNs in pair sample
xx = (pcat.Lo3)[idx1]
idx2 = where(min(pcat.bclass,dim=1) ge bc) ; BAGN
xxx = pcat[idx2].Lo3
xtit = '[OIII] Luminosity'
ytit = 'Number per bin'
xr=[37,43]
yr=[-0.1,1]*250
nboot = 1000
plot,[0,1],[0,1],/nodata,xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit
plothist,alog10(x),bin=1.0,xr=xr,yr=yr,/overplot,color=254
tmp = bootstrap_median(alog10(x),nboot=nboot)
ver,tmp,color=254,lines=2
print,n_elements(x),tmp
plothist,alog10(xx),bin=1.0,xr=xr,/overplot,color=100
tmp = bootstrap_median(alog10(xx),nboot=nboot)
ver,tmp,color=100,lines=2
print,n_elements(xx),tmp
plothist,alog10(xxx),bin=1.0,xr=xr,/overplot,color=200
tmp = bootstrap_median(alog10(xxx),nboot=nboot)
ver,tmp,color=200,lines=2
print,n_elements(xxx),tmp
multiplot,/dox,/doy

; Bootstrap results, 68% conf intervals
;         389       39.409889       39.454449       39.490997
;          54       39.688398       39.890058       40.033798
;          32       39.856088       40.040789       40.359086
;         114       39.905385       39.954185       40.024246
;          23       40.217978       40.289199       40.359086
;          10       40.448016       40.560171       40.578446

endfor

save_screen,'pngs/lo3_mi.png'

stop

bin = 1.0
for mi = -18.5,-23.5,-1 do begin
;for mi = -22.5,-22.5,-1 do begin
	s1 = between(x,mi-bin/2,mi+bin/2)
	s2 = between(xx,mi-bin/2,mi+bin/2)
	phi_all = total(y[s1])/1e6
	phi_agn = total(yy[s2])/1e6
	print,mi,n_elements(s1),n_elements(s2),alog10(phi_all),alog10(phi_agn),phi_agn/phi_all
endfor

;     -18.5000         438          24     -2.36741     -3.72630    0.0437631
;     -19.5000         519          54     -2.45721     -3.49567    0.0915248
;     -20.5000         569         124     -2.57365     -3.25464     0.208456
;     -21.5000         512         111     -3.01679     -3.67157     0.221423
;     -22.5000         412          68     -3.87747     -4.67644     0.158867
;     -23.5000          48           3     -5.21560     -6.58080    0.0431318

;     -18.5000         438           3     -2.36741     -4.64965   0.00522110
;     -19.5000         519           8     -2.45721     -4.29126    0.0146540
;     -20.5000         569          26     -2.57365     -3.94707    0.0423231
;     -21.5000         512          43     -3.01679     -4.05938    0.0906585
;     -22.5000         412          31     -3.87747     -5.09261    0.0609339
;     -23.5000          48           3     -5.21560     -6.58080    0.0431318

end
