; load all DR14 catalog
ncat = mrdfits('nuc1kpc.fits',1) ; 2718 
s = where((ncat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0) ; 2618
ncat = ncat[s]

; load galaxy pair catalog
pcat = mrdfits('pair1kpc.fits',1)
pcat.plateifu = strtrim(pcat.plateifu) ; 176
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin)  ; 105
pcat = pcat[s]
; Pair flag
match2,pcat.NSA_NSAID,ncat.NSA_NSAID,sa,sb
flg_pair = sb ge 0

; report control sample statistics
s = where(~flg_pair) ; 2513
str = ncat
forprint,['Main Sample','Unclass','RG', 'SF', 'Comp', 'LINER', 'Sey2', 'Sey1'],$
	[n_elements(where(str[s].bclass gt -9)),$
	n_elements(where(str[s].bclass lt 0)),$
	n_elements(where(str[s].bclass eq 0)),$
	n_elements(where(str[s].bclass eq 1)),$
	n_elements(where(str[s].bclass eq 2)),$
	n_elements(where(str[s].bclass eq 3)),$
	n_elements(where(str[s].bclass eq 4)),$
	n_elements(where(str[s].bclass eq 5))],$
	f='(a10,i6)',text='ctrl_sample.txt'

end
