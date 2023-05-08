; compute pair fraction in a volume-limited sample
; using the best-fit f_pair(M_s,Rifu) relation

h = 0.7

; load DR14 observed catalog
ncat = mrdfits('../7.sample/nuc1kpc.fits',1)
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $
	ncat.NSA_ELPETRO_MASS/h^2 gt 10.^(10.66-0.5) and $)
	ncat.NSA_ELPETRO_MASS/h^2 lt 10.^(10.66+0.5))
ncat = ncat[s] 
; load target catalog w/ new weights based on DR14
target = mrdfits('../../sample/target/DR14_newcosm.fits',1) 
; match target catalog to DR14 catalog
match2,ncat.NSA_NSAID,target.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targDR14 = target[sa]

; pair fraction
rpmax = 20/h ; kpc
Ms = ncat.NSA_ELPETRO_MASS/h^2 ; Msun
fpair = 0.093 * (rpmax/30.) * Ms/1e11
W = targDR14.esweight
npair = total(W*fpair)
ntot = total(W)
print,'F_pair = ',npair,ntot,npair/ntot,1.3/9.3*npair/ntot

; range in Ms
print,'logM =',minmax(alog10(Ms))
; range in Mi
Mi = ncat.NSA_ELPETRO_ABSMAG[5]+5*alog10(0.7)
print,'M_i = ',minmax(Mi)
; range in z
print,'z = ',minmax(ncat.nsa_z),mean(ncat.nsa_z),median(ncat.nsa_z)

; load pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1) ; 176 - -18.28 < Mi < -25.20
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin) ; 105
pcat = pcat[s]
dMs = abs(pcat.mstar[0]-pcat.mstar[1])
help,where(dMs lt alog10(3)),where(dMs gt alog10(3))

end
