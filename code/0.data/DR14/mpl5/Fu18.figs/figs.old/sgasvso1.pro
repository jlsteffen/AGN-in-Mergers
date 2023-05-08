;+
; Sigma(gas) vs. Sigma(star)
; Sigma(gas) vs. [O I]/Ha 
;-

PRO mkplot, x,y,bc, px,py,dpx,dpy,pbc, _REF_EXTRA=extra
	
	; plot frame
	plot,[0,1],[0,1],/nodata,/xs,/ys,_extra=extra

	; plot all DR14
	s = where(bc eq 0) ; retired galaxies
	plots,x[s],y[s],psym=cgSymCat(16),$
		color=cgcolor('gray'),noclip=0,syms=0.5
	s = where(bc ge 1)
	plots,x[s],y[s],psym=cgSymCat(16),$
		color=cgcolor('dark gray'),noclip=0,syms=0.5
	
	; plot pair sample
	; non-BAGNs
	; s_other = where(min(pcat.bclass,dim=1) lt 2)
	;plots,xx[*,s_other],yy[*,s_other],psym=cgSymCat(16),color=cgcolor('gray'),noclip=0,syms=syms
	;plots,xx[*,s_other],yy[*,s_other],psym=cgSymCat(9),color=cgcolor('black'),noclip=0,syms=syms
	
	; BAGNs
	syms = 1.5	
	s = where(min(pbc,dim=1) ge 2)
	;for i=0,n_elements(s)-1 do $ ; connecting lines
	;	oplot,px[*,s[i]],py[*,s[i]],lines=0,color=cgcolor('green')
	plots,px[0,s],py[0,s],psym=cgSymCat(15),noclip=0,syms=syms,color=cgcolor('red')
	plots,px[0,s],py[0,s],psym=cgSymCat(6),noclip=0,syms=syms
	plots,px[1,s],py[1,s],psym=cgSymCat(16),noclip=0,syms=syms,color=cgcolor('sky blue')
	plots,px[1,s],py[1,s],psym=cgSymCat(9),noclip=0,syms=syms
	oploterror,px[*,s],py[*,s],dpx[*,s],dpy[*,s],psym=3,/nohat
end

; load galaxy pair catalog
pcat = mrdfits('../7.pairsample.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 168
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 128
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 122
	pcat.sep_arcsec le pcat.IFURin)  ; 114
pcat = pcat[s]
; load all DR14 catalog
ncat = mrdfits('../7.nucsample.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]

; setup plot
loadct,0
psfile = 'figs/sgasvso1.eps'
setps,psfile,15*2.2,15,font='helvetica'
titsize = 1.5
labelsize = 1.2
pos1=[0.08,0.11,0.49-0.01/2,0.96]
pos2=[0.58+0.01/2,0.11,0.99,0.96]

; S(gas) vs. [O I]/Ha
xr=[-2.5,0.0]
yr=[0,350]
xtit='log([O I]/H'+cgSymbol('alpha')+')'
ytit=textoidl('H\alpha Velocity Dispersion (km s^{-1})')
x = ncat.o1ha
y = ncat.sgas
bc = ncat.bclass
px = pcat.o1ha
py = pcat.Sgas
dpx = pcat.do1ha
dpy = pcat.dSgas
pbc = pcat.bclass
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,xtickint=0.5,ytickint=100,$
	xtit=xtit,ytit=ytit,pos=pos1,charsize=titsize
; legend
al_legend,'(a)',/top,/right,box=0,chars=2,margin=0
al_legend,['Retired (DR14)','Others (DR14)','BAGN Primaries','BAGN Secondaries'],$
	psym=[16,16,15,16],color=cgcolor(['gray','dark gray','red','sky blue']),$
	/top,/left,background=cgcolor('white')

; Sig(gas) vs Sig(star)
xr = yr
xtit=textoidl('Stellar Velocity Dispersion (km s^{-1})')
ytit=textoidl('H\alpha Velocity Dispersion (km s^{-1})')
x = ncat.sstar
px = pcat.sstar
dpx = pcat.dSstar
mkplot,x,y,bc, px,py,dpx,dpy,pbc,xr=xr,yr=yr,xtickint=100,ytickint=100,$
	xtit=xtit,ytit=ytit,pos=pos2,charsize=titsize
oplot,xr,xr,lines=2
; legend
al_legend,'(b)',/top,/right,box=0,chars=2,margin=0
;al_legend,['Retired (DR14)','Others (DR14)','BAGN Primaries','BAGN Secondaries'],$
;	psym=[16,16,15,16],color=cgcolor(['gray','dark gray','red','sky blue']),$
;	/top,/left,background=cgcolor('white')

device,/close
spawn,'gv '+psfile+' &'

end
