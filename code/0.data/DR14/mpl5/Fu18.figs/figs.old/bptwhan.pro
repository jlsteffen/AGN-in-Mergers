PRO mkplot, x,y,bc, px,py,dpx,dpy,pbc, _REF_EXTRA=extra
	
	; plot frame
	plot,[0,1],[0,1],/nodata,/xs,/ys,_extra=extra

	; all DR14 RGs & SFGs
	syms = 0.5
	s = where(bc eq 0) ; retired galaxies
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('gray'),noclip=0,syms=syms
	s = where(bc eq 1) ; SFGs
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('dark gray'),noclip=0,syms=syms
	; all DR14 AGNs
	syms = 0.8
	s = where(bc eq 2)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('royal blue'),noclip=0,syms=syms
	s = where(bc eq 3)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('forest green'),noclip=0,syms=syms
	s = where(bc eq 4)
	plots,x[s],y[s],psym=cgSymCat(16),color=cgcolor('red'),noclip=0,syms=syms
	
	;; Pair Sample - BAGNs
	;syms = 1.5
	;s = where(min(pbc,dim=1) ge 2)
	;px2 = px[*,s]
	;py2 = py[*,s]
	;pbc2 = pbc[*,s]
	;dpx2 = dpx[*,s]
	;dpy2 = dpy[*,s]
	;i = where(pbc2 eq 2)
	;plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('royal blue')
	;plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	;i = where(pbc2 eq 3)
	;plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('green')
	;plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	;i = where(pbc2 eq 4)
	;plots,px2[i],py2[i],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('red')
	;plots,px2[i],py2[i],psym=cgSymCat(6),noclip=0,syms=syms
	;; add error bars
	;oploterror,px[*,s],py[*,s],dpx[*,s],dpy[*,s],psym=3,/nohat
end

; load all DR14 catalog
ncat = mrdfits('../7.nucsample.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]
; load galaxy pair catalog
pcat = mrdfits('../7.pairsample.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 168
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 128
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 122
	pcat.sep_arcsec le pcat.IFURin)  ; 114
pcat = pcat[s]
; match target catalog to DR14 catalog
; load full target catalog w/ new weights
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ntarg = target[sa]
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
ptarg = target[sa]

;;;;;;;;;;;;;;;;;;;;
; NII/Ha BPT diagram
; SII/Ha BPT
; WHAN diagram
;;;;;;;;;;;;;;;;;;;;;
; setup plot
loadct,0
psfile = 'figs/bptwhan.eps'
setps,psfile,15*1.1*3,15*2,font='helvetica'
titsize = 1.5
labelsize = 1.2
multiplot,/default
multiplot,[3,2],xgap=0.025,ygap=0.02,/dox,/doy,/row

; [N II] BTP diagram
xr=[-1.3,0.6]
yr=[-1.2,1.4]
xtit='' ;'log([N II]/H'+cgSymbol('alpha')+')'
ytit='log([O III]/H'+cgSymbol('beta')+')'
x = ncat.n2ha
y = ncat.o3hb
bc = ncat.bclass
px = pcat.n2ha
py = pcat.o3hb
dpx = pcat.dn2ha
dpy = pcat.do3hb
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; Kewley06 extreme SB line
x = range(-2.5,1.5,40)
oplot,x,1.19+0.61/((x<0.47)-0.47),lines=2 
; Kauffmann03 pure SF line
oplot,x,1.30+0.61/((x<0.05)-0.05),lines=0 
; Schawinski07 Seyfert/Liner line
x = range(-0.18,1.5,40)
oplot,x,1.05*x+0.45,lines=3
; labels
xyouts,-0.02,-1.0,'Composite',align=0.5,/data,charsize=labelsize
xyouts,-0.35,0.95,'Seyfert',align=1.0,/data,charsize=labelsize
xyouts,0.53,0.0,'LINER',align=1.0,/data,charsize=labelsize
xyouts,-1.2,-0.5,'Starforming',align=0,/data,charsize=labelsize
; Reddening arrow
; for A(V)=5 mag => E(B-V)=A(V)/R, R=3.1
lambda = [6585.28,6564.61,5008.240,4862.683] 	
ccm_UNRED,lambda,[1.0,1.0,1.0,1.0],5/3.1,funred,R_V=3.1
px = -1.1 ; xr[1]-0.70*(xr[1]-xr[0])-0.05 
py = -1.0 ; yr[1]-0.9*(yr[1]-yr[0])-0.1
ax = alog10(funred[0]/funred[1])  ; arrow x 
ay = alog10(funred[2]/funred[3])  ; arrow y
l = sqrt(ax^2+ay^2)               ; length of the arrow
x1 = px+ax & y1 = py+ay
rsin = 0.5*sin(22.5*!dtor)
rcos = 0.5*cos(22.5*!dtor)
; 0.5 : Length of arrowhead legs relative to vectorlength.
; 22.5 deg = 1/2 times the angle between the arrowhead legs.
plots,[px,x1,x1-(ax*rcos+ay*rsin), x1,x1-(ax*rcos-ay*rsin)], $
      [py,y1,y1-(ay*rcos-ax*rsin), y1,y1-(ay*rcos+ax*rsin)], THICK=4.0
xyouts,noclip=0,px-0.00,py+0.02,'A!DV!N = 5',align=0.5,chars=1.5,chart=5
; legend
;al_legend,['Kauffmann03','Kewley01','Schawinski07'],lines=[0,2,3],/top,/left,background=cgcolor('white')
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0
al_legend,['Retired','Starforming','Composite','LINER','Seyfert'],$
	psym=[16,16,16,16,16],$
	color=cgcolor(['gray','dark gray','royal blue','forest green','red']),$
	/top,/left,background=cgcolor('white')
multiplot,/dox,/doy

; WHAN diagram
xr = [-1.3,0.6]
xtit = 'log([N II]/H'+cgSymbol('alpha')+')'
yr = [-0.5,2.5]
ytit='log(H'+cgSymbol('alpha')+' EW / '+cgSymbol('angstrom')+')'
xtit='log([N II]/H'+cgSymbol('alpha')+')'
x = ncat.n2ha
y = alog10(ncat.WHa)
bc = ncat.bclass
px = pcat.n2ha
py = alog10(pcat.WHa)
dpx = pcat.dn2ha
dpy = (pcat.dWHa/pcat.WHa)/alog(10) ; d log(y) = 1/ln(10) dy/y
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; dividing lines
hor,alog10(3.0),lines=0
oplot,[-0.4,xr[1]],alog10([6.,6.]),lines=3 ;,color=cgcolor('sky blue')
oplot,[-0.4,-0.4],[alog10(3),yr[1]],lines=2 ;,color=cgcolor('red')
; labels
xyouts,0.53,alog10(2),'Retired',align=1.0,/data,charsize=labelsize
xyouts,0.53,alog10(30),'sAGN',align=1.0,/data,charsize=labelsize
xyouts,0.53,alog10(4),'wAGN',align=1.0,/data,charsize=labelsize
xyouts,-1.2,alog10(10),'Starforming',align=0,/data,charsize=labelsize
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0
multiplot,/dox,/doy

; [S II] BTP diagram
xr=[-1.1,0.6]
yr=[-1.2,1.4]
xtit='' ;'log([S II]/H'+cgSymbol('alpha')+')'
ytit='log([O III]/H'+cgSymbol('beta')+')'
x = ncat.s2ha
y = ncat.o3hb
bc = ncat.bclass
px = pcat.s2ha
py = pcat.o3hb
dpx = pcat.ds2ha
dpy = pcat.do3hb
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,$
	xtit=xtit,ytit=ytit,charsize=titsize
; Kewley06 [S II] classification lines
x = range(xr[0],xr[1],40)
oplot,x,1.30+0.72/((x<0.27)-0.32),lines=2 
oplot,x>(-0.4),1.89*(x>(-0.4))+0.76,lines=3 
; labels
xyouts,-0.6,0.95,'Seyfert',align=1.0,/data,charsize=labelsize
xyouts,0.53,0.6,'LINER',align=1.0,/data,charsize=labelsize
xyouts,-0.5,-1.00+0.05,'Starforming',align=0,/data,charsize=labelsize
xyouts,-0.5,-1.13+0.05,'& Composite',align=0,/data,charsize=labelsize
; Reddening arrow
; for A(V)=5 mag => E(B-V)=A(V)/R, R=3.1
lambda = [6726.0,6564.61,5008.240,4862.683] 	
ccm_UNRED,lambda,[1.0,1.0,1.0,1.0],5/3.1,funred,R_V=3.1
px = -0.9 ; xr[1]-0.70*(xr[1]-xr[0])-0.05 
py = -1.0 ; yr[1]-0.9*(yr[1]-yr[0])-0.1
ax = alog10(funred[0]/funred[1])  ; arrow x 
ay = alog10(funred[2]/funred[3])  ; arrow y
l = sqrt(ax^2+ay^2)               ; length of the arrow
x1 = px+ax & y1 = py+ay
rsin = 0.5*sin(22.5*!dtor)
rcos = 0.5*cos(22.5*!dtor)
; 0.5 : Length of arrowhead legs relative to vectorlength.
; 22.5 deg = 1/2 times the angle between the arrowhead legs.
plots,[px,x1,x1-(ax*rcos+ay*rsin), x1,x1-(ax*rcos-ay*rsin)], $
      [py,y1,y1-(ay*rcos-ax*rsin), y1,y1-(ay*rcos+ax*rsin)], THICK=4.0
xyouts,noclip=0,px-0.00,py+0.02,'A!DV!N = 5',align=0.5,chars=1.5,chart=5
; legend
;al_legend,['Kauffmann03','Kewley01','Schawinski07'],lines=[0,2,3],/top,/left,background=cgcolor('white')
al_legend,'(c)',/top,/right,box=0,chars=2,margin=0
multiplot,/dox,/doy

; S(gas) vs. [S II]/Ha
xr=[-1.1,0.6]
yr=[0,350]
xtit='log([S II]/H'+cgSymbol('alpha')+')'
ytit=textoidl('H\alpha Velocity Dispersion (km s^{-1})')
x = ncat.s2ha
y = ncat.sgas
bc = ncat.bclass
px = pcat.s2ha
py = pcat.Sgas
dpx = pcat.ds2ha
dpy = pcat.dSgas
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,xtickint=0.5,ytickint=100,$
	xtit=xtit,ytit=ytit,charsize=titsize
; legend
al_legend,'(d)',/top,/right,box=0,chars=2,margin=0
multiplot,/dox,/doy

; plot NUV-r vs. Mi 
xr = [-17,-25]
xtit = '' ; textoidl('Absolute Magnitude (M_i)')
yr = [1, 7] 
ytit=textoidl('NUV-r')
x = ncat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7)
y = ncat.NSA_ELPETRO_ABSMAG[1]-ncat.NSA_ELPETRO_ABSMAG[4]
bc = ncat.bclass
px = [1,1] # (pcat.NSA_ELPETRO_ABSMAG[5]+5.0*alog10(0.7))
py = [1,1] # (pcat.NSA_ELPETRO_ABSMAG[1]-pcat.NSA_ELPETRO_ABSMAG[4])
dpx = px*0
dpy = py*0
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,ytickint=2,$
	xtit=xtit,ytit=ytit,charsize=titsize
; legend
al_legend,'(e)',/top,/right,box=0,chars=2,margin=0,background=cgcolor('white')
multiplot,/dox,/doy

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

device,/close
spawn,'gv '+psfile+' &'

end