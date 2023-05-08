;+
; Generate tabs/fbagn.tex
; Note: Ellison11 pair selection - r_Petro < 17.77 @ z = 0.0916, Mr_Petro = -20.35
;-

h = 0.7 ; Hubble constant
files = 'tabs/fbagn_quan.tex' ; list of columns
minbclass = 2 ; Kauffmann03
texbase = 'tabs/fbagn'+strc(minbclass)

; load pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1) ; 176 - -18.28 < Mi < -25.20
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin) ; 105
pcat = pcat[s]
; Note that big offsets between IFU center and primary object is fine 
; as long as Pair Sep is less than Rin of the IFU, because these pairs
; would remain in the pair sample even if we recenter the IFU on the object

; load DR14 observed catalog
ncat = mrdfits('../7.sample/nuc1kpc.fits',1)
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0)
ncat = ncat[s] 
; flag_pair: 1 - pair sample, 0 - control sample
flag_pair = intarr(n_elements(ncat))
match2,ncat.NSA_NSAID,pcat.NSA_NSAID,sa,sb
flag_pair[sb] = 1

; load target catalog w/ new weights based on DR14
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to pair catalog
match2,pcat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targpair = target[sa]
; match target catalog to DR14 catalog
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]
help,where(targDR14.MANGA_TARGET1 ne ncat.mngtarg1) ; only 2 mismatches

; load target catalog w/ new weights based on full tile
fulltile = mrdfits('../../sample/target/DR14_newcosm_fulltile.fits',1)
; keep targets in tiled region 
s = where(max(fulltile.manga_tileids,dim=1) ge 0) 
tiled = fulltile[s]
;; targets in both Secondary and Color-Enhanced 
;; - remove them from Secondary
;s = where((tiled.manga_target1 and 2L^11) ne 0 and (tiled.manga_target1 and 2L^12) ne 0)
;tiled[s].manga_target1 = 2L^12

; AGN probability for ncat 
restore,'pd2agn.idlsav'
mass = alog10(ncat.NSA_ELPETRO_MASS/h^2)
z = ncat.NSA_Z
fagnDR14 =  p[0] * exp(-0.5d*(mass-p[1])^2/p[2]^2) * (1+z)^4
; AGN probability for pcat
z = pcat.NSA_Z
mass_t = pcat.NSA_ELPETRO_MASS/h^2 ; total mass
mass_p = alog10(mass_t * 10.^pcat.mstar[0]/$
	(10.^pcat.mstar[0]+10.^pcat.mstar[1]))
mass_s = alog10(mass_t * 10.^pcat.mstar[1]/$
	(10.^pcat.mstar[0]+10.^pcat.mstar[1]))
;mass_p = alog10(pcat.NSA_ELPETRO_MASS/h^2)
;mass_s = mass_p + (pcat.mstar[1] - pcat.mstar[0])
fagn_p =  p[0] * exp(-0.5d*(mass_p-p[1])^2/p[2]^2) * (1+z)^4
fagn_s =  p[0] * exp(-0.5d*(mass_s-p[1])^2/p[2]^2) * (1+z)^4

; bit numbers for different subsamples
sambit = [2L^10, 2L^10+2L^12,2L^11,2L^10+2L^11+2L^12]
snames = ['Primary','Primary+','Secondary','Combined']
wnames = ['pweight','eweight','sweight','esweight']
nboot = 100 ; # of bootstrap

for i=1,n_elements(sambit)-1 do begin
	texfile = texbase+strmid(wnames[i],0,2)+'.tex'
	files = [files,texfile]
	openw,5,texfile
	
	; count # of galaxies and get indices
	; all galaxies within tiled area, whether or not allocated an IFU
	Itile = where((tiled.MANGA_TARGET1 and sambit[i]) ne 0,Ntile)
	; DR14 sample
	Iobsv = where((targDR14.MANGA_TARGET1 and sambit[i]) ne 0,Nobsv) 
	; Control sample
	Icntr = where((targDR14.MANGA_TARGET1 and sambit[i]) ne 0 and $
		flag_pair eq 0, Ncntr) 
	Iagn  = where((targDR14.MANGA_TARGET1 and sambit[i]) ne 0 and $
		flag_pair eq 0 and ncat.bclass ge minbclass, Nagn)
	; Pair sample
	Ipair = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0,Npair)
	; either as AGN
	Ieagn = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0 and $
		(pcat.bclass[0] ge minbclass or pcat.bclass[1] ge minbclass),Neagn) 
	; single AGNs
	Isagn = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0 and $
		(pcat.bclass[0] ge minbclass or pcat.bclass[1] ge minbclass) and $
		min(pcat.bclass,dim=1) lt minbclass,Nsagn)
	; binary AGNs
	Ibagn = where((targpair.MANGA_TARGET1 and sambit[i]) ne 0 and $
		min(pcat.bclass,dim=1) ge minbclass,Nbagn)
	
	; weight := Phi * 1e6 Mpc^3
	xx = execute('wobsv = targDR14[Iobsv].'+wnames[i])
	xx = execute('wcntr = targDR14[Icntr].'+wnames[i])
	xx = execute('wagn  = targDR14[Iagn].'+wnames[i])
	xx = execute('wpair = targpair[Ipair].'+wnames[i])
	xx = execute('weagn = targpair[Ieagn].'+wnames[i])
	xx = execute('wsagn = targpair[Isagn].'+wnames[i])
	xx = execute('wbagn = targpair[Ibagn].'+wnames[i])
	; bootstrap - lower limit, mean, and upper limit.
	wobsv_bs = bootstrap_mean(wobsv,nboot=nboot)
	wcntr_bs = bootstrap_mean(wcntr,nboot=nboot)
	wagn_bs  = bootstrap_mean(wagn ,nboot=nboot)
	wpair_bs = bootstrap_mean(wpair,nboot=nboot)
	weagn_bs = bootstrap_mean(weagn,nboot=nboot)
	wsagn_bs = bootstrap_mean(wsagn,nboot=nboot)
	wbagn_bs = bootstrap_mean(wbagn,nboot=nboot)
	; observed comoving densities
	ncobsv = [total(wobsv),wobsv_bs * Nobsv]/1e6 ; #/Mpc^3
	nccntr = [total(wcntr),wcntr_bs * Ncntr]/1e6
	ncagn  = [total(wagn) ,wagn_bs  * Nagn ]/1e6 
	ncpair = [total(wpair),wpair_bs * Npair]/1e6
	nceagn = [total(weagn),weagn_bs * Neagn]/1e6
	ncsagn = [total(wsagn),wsagn_bs * Nsagn]/1e6
	ncbagn = [total(wbagn),wbagn_bs * Nbagn]/1e6
	; log errors of comoving densities
	; d log(y) = 1/ln(10) dy/y
	dlogncobsv = (ncobsv[3]-ncobsv[1])/2/ncobsv[0]/alog(10)
	dlognccntr = (nccntr[3]-nccntr[1])/2/nccntr[0]/alog(10)
	dlogncagn  = ( ncagn[3]- ncagn[1])/2/ ncagn[0]/alog(10)
        dlogncpair = (ncpair[3]-ncpair[1])/2/ncpair[0]/alog(10)
        dlognceagn = (nceagn[3]-nceagn[1])/2/nceagn[0]/alog(10)
        dlogncsagn = (ncsagn[3]-ncsagn[1])/2/ncsagn[0]/alog(10)
        dlogncbagn = (ncbagn[3]-ncbagn[1])/2/ncbagn[0]/alog(10)

	; expected AGN densities 
	; either AGNs
	weagn = wpair*fagn_p[Ipair]+wpair*fagn_s[Ipair]-$
		wpair*fagn_p[Ipair]*fagn_s[Ipair]
	tmp = bootstrap_mean(weagn,nboot=nboot)*Npair
	logexpeagn = alog10(tmp[1]/1e6)
	dlogexpeagn = (tmp[2]-tmp[0])/2/tmp[1]/alog(10)
	; single AGNs
	wsagn = wpair*fagn_s[Ipair]*(1-fagn_p[Ipair]) +$
	      	wpair*fagn_p[Ipair]*(1-fagn_s[Ipair])
	tmp = bootstrap_mean(wsagn,nboot=nboot)*Npair
	logexpsagn = alog10(tmp[1]/1e6)
	dlogexpsagn = (tmp[2]-tmp[0])/2/tmp[1]/alog(10)
	; binary AGNs
	wbagn = wpair*fagn_p[Ipair]*fagn_s[Ipair]
	tmp = bootstrap_mean(wbagn,nboot=nboot)*Npair
	logexpbagn = alog10(tmp[1]/1e6)
	dlogexpbagn = (tmp[2]-tmp[0])/2/tmp[1]/alog(10)
	; add systematic uncertainties
	dlogexpeagn = sqrt(dlogexpeagn^2+0.04^2)
	dlogexpsagn = sqrt(dlogexpsagn^2+0.04^2)
	dlogexpbagn = sqrt(dlogexpbagn^2+0.04^2)

	; AGN excesses
	logexe = alog10(nceagn[0])-logexpeagn
	dlogexe = sqrt(dlognceagn^2+dlogexpeagn^2)
	logexs = alog10(ncsagn[0])-logexpsagn
	dlogexs = sqrt(dlogncsagn^2+dlogexpsagn^2)
	logexb = alog10(ncbagn[0])-logexpbagn
	dlogexb = sqrt(dlogncbagn^2+dlogexpbagn^2)

	; print table
	fmt0 ='("&",i5)'
	fmt2 ='("&",f5.2,"$\pm$",f4.2)'
	; among all galaxies
	printf,5,f=fmt0,Ntile
	printf,5,f=fmt0,Nobsv
	printf,5,f=fmt2,alog10(ncobsv[0]),dlogncobsv
	; among Control Sample
	printf,5,' '
	printf,5,' '
	printf,5,' '
	printf,5,f=fmt0,Ncntr
	printf,5,f=fmt0,Nagn
	printf,5,f=fmt2,alog10(nccntr[0]),dlognccntr
	printf,5,f=fmt2,alog10(ncagn[0]),dlogncagn
	; among Mergers
	printf,5,' '
	printf,5,' '
	printf,5,' '
	printf,5,f=fmt0,Npair
	printf,5,f=fmt0,Neagn
	printf,5,f=fmt0,Nsagn
	printf,5,f=fmt0,Nbagn	
	printf,5,f=fmt2,alog10(ncpair[0]),dlogncpair
	printf,5,f=fmt2,alog10(nceagn[0]),dlognceagn
	printf,5,f=fmt2,alog10(ncsagn[0]),dlogncsagn
	printf,5,f=fmt2,alog10(ncbagn[0]),dlogncbagn
	printf,5,f=fmt2,logexpeagn,dlogexpeagn
	printf,5,f=fmt2,logexpsagn,dlogexpsagn
	printf,5,f=fmt2,logexpbagn,dlogexpbagn	
	; BAGN Excesses
	printf,5,' '
	printf,5,' '
	printf,5,' '
	printf,5,f=fmt2,logexe,dlogexe
	printf,5,f=fmt2,logexs,dlogexs
	printf,5,f=fmt2,logexb,dlogexb

	close,5
	spawn,'sed -i '''' ''s/+/$+$/g'' '+texfile
	spawn,'sed -i '''' ''s/-/$-$/g'' '+texfile
	spawn,'sed -i '''' ''s/$ /$/g'' '+texfile
endfor

files = [files,'tabs/fbagn_unit.tex']
; combine all the columns
spawn,'paste '+strjoin(files,' ')+' > tabs/fbagn.tex'
; remove intermediate files
n = n_elements(files)
spawn,'del '+strjoin(files[1:n-2],' ')

end
