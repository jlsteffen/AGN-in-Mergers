;+
; Generate tabs/fbagn1-2.tex
; Compute f_BAGN among mergers based on comoving volume weights
; Note: Ellison11 pair selection - r_Petro < 17.77 @ z = 0.0916, Mr_Petro = -20.35
;-

minbclass = 2 ; Kauffmann03
minbclass = 3 ; Kewley01

files = 'tabs/fbagn_quan.tex' ; list of columns

FOR kk=0,1 do begin

CASE kk OF
	0: Mi_cut = -21.0
	1: Mi_cut = -18.0
ENDCASE
texbase = 'tabs/fbagn'+strc(minbclass*10 + kk)

; load pair catalog
pcat = mrdfits('../7.pairsample.fits',1) ; 174
Rin = (pcat.IFUdiam+3.5-2.31)/2*sqrt(3)/2 ; inner circle radius
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 166
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 130
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 124
	pcat.sep_arcsec le Rin and $ ; 116
	pcat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut) ; 95 for -21.5
pcat = pcat[s]
; Note that big offsets between IFU center and primary object is fine 
; as long as Pair Sep is less than Rin of the IFU, because these pairs
; would remain in the pair sample even if we recenter the IFU on the object

; load DR14 observed catalog
drpall = mrdfits('../../sample/matched/drpall.fits',1)
s = where((drpall.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $
	drpall.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut)
drpall = drpall[s] ; 1135

; load target catalog w/ new weights based on DR14
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targpair = target[sa]
; match target catalog to DR14 catalog
match2,drpall.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]
help,where(targDR14.MANGA_TARGET1 ne drpall.mngtarg1) ; only 2 mismatches

; load target catalog w/ new weights based on full tile
fulltile = mrdfits('../../sample/target/DR14_newcosm_fulltile.fits',1)
; keep targets in tiled region, allocated an IFU, and brighter than Mi_cut
s = where(max(fulltile.manga_tileids,dim=1) ge 0 and $
	fulltile.IFUDESIGNSIZEMAIN gt 0 and $
	fulltile.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7) lt Mi_cut)
tiled = fulltile[s]
; targets in both Secondary and Color-Enhanced 
; - remove them from Secondary
s = where((tiled.manga_target1 and 2L^11) ne 0 and (tiled.manga_target1 and 2L^12) ne 0)
tiled[s].manga_target1 = 2L^12

; bit numbers for different subsamples
sambit = [2L^10,2L^10+2L^12,2L^11,2L^10+2L^11+2L^12]
snames = ['Primary','Primary+','Secondary','Combined']
wnames = ['pweight','eweight','sweight','esweight']
nboot = 500 ; # of bootstrap

for i=1,n_elements(sambit)-1 do begin
	texfile = texbase+strmid(wnames[i],0,2)+'.tex'
	files = [files,texfile]
	openw,5,texfile
	; count # of galaxies and get indices
	Itile = where((tiled.MANGA_TARGET1 and sambit[i]) ne 0,Ntile)
	Iobsv = where((targDR14.MANGA_TARGET1 and sambit[i]) ne 0,Nobsv) 
	Ipair = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0,Npair)
	Isagn = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0 and $
		max(pcat.bclass,dim=1) ge minbclass,Nsagn)
	Ibagn = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0 and $
		min(pcat.bclass,dim=1) ge minbclass,Nbagn)
	; weight := Phi * 1e6 Mpc^3
	xx = execute('wtile = tiled[Itile].'+wnames[i])
	xx = execute('wobsv = targDR14[Iobsv].'+wnames[i])
	xx = execute('wpair = targpair[Ipair].'+wnames[i])
	xx = execute('wsagn = targpair[Isagn].'+wnames[i])
	xx = execute('wbagn = targpair[Ibagn].'+wnames[i])
	; bootstrap - lower limit, mean, and upper limit.
	wobsv_bs = bootstrap_mean(wobsv,nboot=nboot)
	wpair_bs = bootstrap_mean(wpair,nboot=nboot)
	wsagn_bs = bootstrap_mean(wsagn,nboot=nboot)
	wbagn_bs = bootstrap_mean(wbagn,nboot=nboot)
	; compute comoving densities
	ncobsv = [total(wobsv),wobsv_bs * Nobsv]/1e6 ; #/Mpc^3
	ncpair = [total(wpair),wpair_bs * Npair]/1e6
	ncsagn = [total(wsagn),wsagn_bs * Nsagn]/1e6
	ncbagn = [total(wbagn),wbagn_bs * Nbagn]/1e6
	; log errors
	; d log(y) = 1/ln(10) dy/y
	dlogncobsv = (ncobsv[3]-ncobsv[1])/2/ncobsv[0]/alog(10)
        dlogncpair = (ncpair[3]-ncpair[1])/2/ncpair[0]/alog(10)
        dlogncsagn = (ncsagn[3]-ncsagn[1])/2/ncsagn[0]/alog(10)
        dlogncbagn = (ncbagn[3]-ncbagn[1])/2/ncbagn[0]/alog(10)
	; log Pair fractions
	logfpair = alog10(ncpair[0])-alog10(ncobsv[0])
	dlogfpair = sqrt(dlogncpair^2+dlogncobsv^2)
	; log SAGN fractions
	flogsagn = alog10(ncsagn[0])-alog10(ncpair[0])
	dlogfsagn = sqrt(dlogncsagn^2+dlogncpair^2)
	; log BAGN fractions
	flogbagn = alog10(ncbagn[0])-alog10(ncpair[0])
	dlogfbagn = sqrt(dlogncbagn^2+dlogncpair^2)
	; Pair fractions
	fpair = eratio(ncpair[0],(ncpair[3]-ncpair[1])/2,ncobsv[0],(ncobsv[3]-ncobsv[1])/2)
	; SAGN fractions
	fsagn = eratio(ncsagn[0],(ncsagn[3]-ncsagn[1])/2,ncpair[0],(ncpair[3]-ncpair[1])/2)
	; BAGN fractions
	fbagn = eratio(ncbagn[0],(ncbagn[3]-ncbagn[1])/2,ncpair[0],(ncpair[3]-ncpair[1])/2)
	; excess from expected random binary AGN fraction 
	; dx^2 = 2x dx
	ex1 = eratio(fbagn[0],fbagn[1],fsagn[0]^2,2*fsagn[0]*fsagn[1])

	; print table
	fmt0 ='("&",i5)'
	fmt1 ='("&",f5.2)'
	fmt2 ='("&",f5.2,"$\pm$",f4.2)'
	fmt3 ='("&",f4.1,"$\pm$",f4.1,"\%")'
	printf,5,f=fmt0,Ntile
	printf,5,f=fmt0,Nobsv
	printf,5,f=fmt0,Npair
	printf,5,f=fmt0,Nsagn
	printf,5,f=fmt0,Nbagn	
	printf,5,f=fmt1,alog10(total(wtile)/1e6)
	printf,5,f=fmt2,alog10(ncobsv[0]),dlogncobsv
	printf,5,f=fmt2,alog10(ncpair[0]),dlogncpair
	printf,5,f=fmt2,alog10(ncsagn[0]),dlogncsagn
	printf,5,f=fmt2,alog10(ncbagn[0]),dlogncbagn
	printf,5,f=fmt3,fpair[0]*1e2,fpair[1]*1e2
	printf,5,f=fmt3,fsagn[0]*1e2,fsagn[1]*1e2
	printf,5,f=fmt3,fbagn[0]*1e2,fbagn[1]*1e2
	close,5
	spawn,'sed -i '''' ''s/+/$+$/g'' '+texfile
	spawn,'sed -i '''' ''s/-/$-$/g'' '+texfile
	spawn,'sed -i '''' ''s/$ /$/g'' '+texfile
endfor

ENDFOR

; combine all the columns
files = [files,'tabs/fbagn_unit.tex']
spawn,'paste '+strjoin(files,' ')+' > '+texbase+'.tex'
spawn,'del '+strjoin(files[1:6],' ')

end
