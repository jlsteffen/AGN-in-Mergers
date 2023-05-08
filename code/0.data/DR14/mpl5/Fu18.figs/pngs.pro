; load DRPALL
drpall = mrdfits('../../sample/matched/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)
; load target catalog
target = mrdfits('../../sample/matched/target.fits',1) ; 2718
; only keep main galaxy sample
s = where((drpall.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 166
drpall = drpall[s]
target = target[s]
; IFU offsets
offset = sphdist(drpall.ifura,drpall.ifudec,drpall.objra,drpall.objdec,/deg)*3600

; load pair catalog
pcat = mrdfits('../8.class.fits',1)
pcat.plateifu = strtrim(pcat.plateifu)
; only keep main galaxy sample
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 166 
pcat = pcat[s]
offs = sphdist(pcat.ifura,pcat.ifudec,pcat.objra,pcat.objdec,/deg)*3600
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targ = target[sa]
; verify matching results
help,where(targ.NSA_NSAID ne pcat.NSA_NSAID)
dis = sphdist(pcat.objra,pcat.objdec,targ.object_ra,targ.object_dec,/deg)*3600
print,minmax(dis)

; identify IFUs showing large offsets
idx = where(offset gt 2) 
manga_regroup,drpall[idx].plateifu,'bigoffset/all/',indir='../../../1.sdsspng/'
idx = where(offs gt 2) 
manga_regroup,pcat[idx].plateifu,'bigoffset/pairs/',indir='../../../1.sdsspng/'

; remove pairs w/ separations greater than inner circle radius
;help,where(pcat.sep_arcsec gt pcat.IFUdiam/2*sqrt(3)/2)
;<Expression>    LONG      = Array[19]
;help,where(pcat.sep_arcsec gt pcat.IFUdiam/2*sqrt(3)/2 and strtrim(pcat.class) eq 'BAGN')
;<Expression>    LONG      = Array[1]

; make some plots
window,0,xsize=500,ysize=500

; Pair separation vs. Effective Radius (Re)
; NSA_ELPETRO_TH50_R	float32[1]	arcsec	Elliptical Petrosian 50% light radius in SDSS r-band (NSA)
plot,targ.NSA_ELPETRO_TH50_R,pcat.sep_arcsec,psym=6,xr=[0,30],yr=[0,30],$
	xtit='ELPETRO effective radius (arcsec)',ytit='Pair Separation (arcsec)'
oplot,[0,30],[0,30],lines=2,color=cgcolor('red')
save_screen,'pngs/Sep_vs_Re.png'

; Re vs. IFU radius
plot,targ.NSA_ELPETRO_TH50_R,pcat.IFUrad,psym=6,xr=[0,30],yr=[0,30],$
	xtit='ELPETRO effective radius (arcsec)',ytit='IFU radius (arcsec)'
oplot,[0,30],[0,30],lines=2,color=cgcolor('red')
help,where(targ.NSA_ELPETRO_TH50_R lt pcat.IFUrad) ; 158
help,where(targ.NSA_ELPETRO_TH50_R*1.5 lt pcat.IFUrad) ; 128
help,where(targ.NSA_ELPETRO_TH50_R*2.5 lt pcat.IFUrad) ; 41
save_screen,'pngs/Rifu_vs_Re.png'

; IFU radius vs. Pair Sep
plot,pcat.sep_arcsec,pcat.IFUrad,psym=6,xr=[0,30],yr=[0,30],$
	xtit='Pair Sep (arcsec)',ytit='IFU radius (arcsec)'
oplot,[0,30],[0,30],lines=2,color=cgcolor('yellow')
s = where(strtrim(pcat.class) eq 'BAGN')
plots,pcat[s].sep_arcsec,pcat[s].IFUrad,psym=8,syms=3,color=cgcolor('red')
s = where(offs gt 2)
plots,pcat[s].sep_arcsec,pcat[s].IFUrad,psym=7,syms=2,color=cgcolor('green')
s = where(pcat.sep_arcsec gt (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2)
plots,pcat[s].sep_arcsec,pcat[s].IFUrad,psym=1,syms=2,color=cgcolor('red')
save_screen,'pngs/Rifu_vs_Sep.png'
;idx = where(pcat.sep_arcsec gt pcat.IFUrad) ; 12
idx = where(pcat.sep_arcsec gt (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2) ; 12
manga_regroup,pcat[idx].plateifu,'bigoffset/sep_gt_Rin/',indir='../../../1.sdsspng/'

; IFU design size vs. IFU diameter
plot,drpall.ifudesignsize,drpall.ifudiam,psym=6,$
	xtit='# of Fibers',ytit='IFU diameter (arcsec)'
x = range(0,140)
oplot,x,sqrt(x)*2.88 ; best-fit relation
save_screen,'pngs/Difu_vs_Nfiber.png'

; Pair Sep distribution vs. IFU size
plothist,pcat.sep_arcsec,bin=1
ver,[12.5,17.5,22.5,27.5,32.5]/2,lines=2,color=cgcolor('red')
save_screen,'pngs/sep_dist.png'

; compute binary fraction vs. kpc separation 
; by counting the # of each directly -- this is not the correct method
bin = 5.0*2 ; kpc
plothist,pcat.sep_kpc,xh,yh,bin=bin,/noplot  
plothist,pcat[s_bagn].sep_kpc,xh2,yh2,bin=bin,/noplot
nbin = n_elements(yh2)
f_bagn = fltarr(nbin,5)
for i=0,nbin-1 do f_bagn[i,*] = binormial_ci(yh2[i],yh[i])
; plot results
setx
wind,xsize=600,ysize=500
plot,xh2,f_bagn[*,0],psym=6,syms=2,xr=[-2,30],yr=[0,0.2],/xs,/ys,$
	xtit='Sep (kpc)',ytit='Fraction of BAGN among Mergers'
oploterror,xh2,f_bagn[*,0],xh2*0+bin/2,f_bagn[*,1],/lobar
oploterror,xh2,f_bagn[*,0],xh2*0+bin/2,f_bagn[*,2],/hibar
plots,pcat[s_bagn].sep_kpc,s_bagn*0+0.02,psym=4,color=cgcolor('red')
plots,pcat.sep_kpc,pcat.sep_kpc*0+0.01,psym=4,color=cgcolor('royal blue')
save_screen,'pngs/BAGN_frac.png'

end
