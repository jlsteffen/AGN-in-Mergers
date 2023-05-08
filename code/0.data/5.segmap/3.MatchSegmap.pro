; pair sample
ncomps = mrdfits('../4.pairs/ncomp.fits',1)
pair = ncomps[where(ncomps.ngood ge 2,npair)]
pair.plateifu = strtrim(pair.plateifu)

; specObj catalog that incl. RA/Dec of objects inside IFU
proc = mrdfits('../0.manga_specobj.fits',1)
proc.plateifu = strtrim(proc.plateifu)

for i=0,npair-1 do begin
;for i=526,526 do begin
	; specObj catalog
	sobj = proc[where(proc.plateifu eq pair[i].plateifu,nobj)]
	; sextractor catalog
	scat = rsex('rsex_cleanoff/'+pair[i].plateifu+'.cat')
	; match the two
	CC,sobj.ra,sobj.dec,scat.ALPHA_J2000,scat.DELTA_J2000,0.5,1.5,dis=dis,mid=mid
	scat_matched = scat[mid]
	; deal w/ missing objects in sextractor catalog
	; by emptying the apparent match
	ind = where(dis lt 0 or dis gt 1,nnon)
	if nnon gt 0 then begin
		print,i,' ',pair[i].plateifu
		zero_struct,scat_matched,ind
		scat_matched[ind].number = -1
	endif
	; save result
	mwrfits,scat_matched,'rseg/'+pair[i].plateifu+'_cat.fits',/create

	; generate matched segmentation map
	smap = mrdfits('rsex_cleanoff/'+pair[i].plateifu+'_seg.fits',0,h)
	smap_matched = smap*0 ; to save matched segmap
	for j=0,nobj-1 do begin
		; sanity check
		if j ne sobj[j].index then stop
		
		; select pixels belonging to the object
		; and avoid overwriting other segments
		idx = where(smap eq scat_matched[j].number and $
			smap_matched eq 0,npix)
		if npix gt 0 then begin
			smap_matched[idx] = sobj[j].index+1
		endif else begin
			; deal w/ missing objects in sextractor catalog
			; by adding a circlular region
			adxy,h,sobj[j].ra,sobj[j].dec,xc,yc
			dim = (size(smap))[1:2]
			dist_ellipse,dismap,dim[0:1],xc,yc,1.0,0.0
			s = where(dismap lt 2.5)
			smap_matched[s] = sobj[j].index+1
		endelse
	end
	mwrfits,smap_matched,'rseg/'+pair[i].plateifu+'_seg.fits',h,/create

endfor

end
