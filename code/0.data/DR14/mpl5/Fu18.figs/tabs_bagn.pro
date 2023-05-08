; generate a TEX table summarizing the BAGN sample properties

; load pair catalog
pcat = mrdfits('../7.sample/pair1kpc.fits',1) ; 156
pcat.plateifu = strtrim(pcat.plateifu)
; restrict to main MaNGA sample
s = where((pcat.mngtarg1 and 2L^10+2L^11+2L^12) ne 0 and $ ; 149
	abs(pcat.mstar[0]-pcat.mstar[1]) le 1.0 and $ ; 119
	abs(pcat.vstar[0]-pcat.vstar[1]) le 600 and $ ; 113
	pcat.sep_arcsec le pcat.IFUrin) ; 105
pcat = pcat[s]

; array for different subsamples
; select main galxy sample - Primary, Color Enhanced, and Secondary
; maskbits MANGA_TARGET1    10  PRIMARY_v1_2_0 " "
; maskbits MANGA_TARGET1    11  SECONDARY_v1_2_0 " "
; maskbits MANGA_TARGET1    12  COLOR_ENHANCED_v1_2_0 " "
samp = intarr(n_elements(pcat))
samp[where((pcat.mngtarg1 and 2L^10) ne 0)] = 1
samp[where((pcat.mngtarg1 and 2L^11) ne 0)] = 2
samp[where((pcat.mngtarg1 and 2L^12) ne 0)] = 3

; match w/ target catalog
targ = mrdfits('../../sample/matched/target.fits',1)
match2,pcat.NSA_NSAID,targ.NSA_NSAID,sa,sb
help,where(sa lt 0) ; should be -1
targ = targ[sa]
; verify matching results
help,where(targ.NSA_NSAID ne pcat.NSA_NSAID)
dis = sphdist(pcat.objra,pcat.objdec,targ.object_ra,targ.object_dec,/deg)*3600
print,minmax(dis)

; replace Lo3 measurements for RGs
ind = where(pcat.bclass[0] eq 0)
pcat[ind].lo3[0] = 0
ind = where(pcat.bclass[1] eq 0)
pcat[ind].lo3[1] = 0

for ii=0,1 do begin
	if ii eq 0 then begin
		; BAGN table
		s = where(min(pcat.bclass,dim=1) ge 2,ct)
		print,'minmax Lo3:',minmax(pcat[s].lo3)
		texfile = 'tabs/bagn.tex'
	endif else begin
		; all pair table
		s = where(min(pcat.bclass,dim=1) lt 2,ct)
		;ct = n_elements(pcat)
		;s = indgen(ct) 
		print,'minmax Lo3:',minmax(pcat[s].lo3)
		texfile = 'tabs/merger.tex'
	endelse
	
	openw,5,texfile
	fmt1 ='(i4,"&",i5,"&",i1,3("&",f5.1),'+$
		'"&",f9.5,"&",f+9.5,"&",f7.5,3("&",f6.1),"&",i1,"\\")' 
	fmt2 ='(4x,"&",5x,"&",1x,3("&",5x  ),'+$
		'"&",f9.5,"&",f+9.5,"&",f7.5,3("&",f6.1),"&",i1,"\\")' 
	for i=0,ct-1 do begin
		si = s[i] ; simplify index
		; primary = manga target
		k = 0
		printf,5,f=fmt1,$
			pcat[si].plate,pcat[si].ifudsgn,samp[si],$
			pcat[si].sep_arcsec,pcat[si].sep_kpc,$
			abs(pcat[si].vstar[0]-pcat[si].vstar[1]),$
			pcat[si].RAs[k],pcat[si].Decs[k],$
			pcat[si].NSA_Z+pcat[si].vstar[k]/c_kms,$
			pcat[si].sstar[k],pcat[si].mstar[k],$
			alog10(pcat[si].lo3[k]),$
			pcat[si].bclass[k]
		; companion
		k = 1
		printf,5,f=fmt2,$
			pcat[si].RAs[k],pcat[si].Decs[k],$
			pcat[si].NSA_Z+pcat[si].vstar[k]/c_kms,$		
			pcat[si].sstar[k],pcat[si].mstar[k],$
			alog10(pcat[si].lo3[k]),$
			pcat[si].bclass[k]
	endfor
	close,5
	spawn,'sed -i '''' ''s/+/$+$/g'' '+texfile
	spawn,'sed -i '''' ''s/-/$-$/g'' '+texfile
	spawn,'sed -i '''' ''s/\*\*\*\*\*\*/  \\nod/g'' '+texfile
endfor

end
