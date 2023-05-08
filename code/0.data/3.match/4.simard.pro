; MaNGA SpecObj catalog
spobj = mrdfits('../0.manga_specobj.fits',1)

;;;;;;;;;;;;;;;;
; match w/ Table 3 of Simard+11: Pure Sersic profile fitting
; Sextractor + GIM2D photometry
;---------------
si11 = mrdfits('/Volumes/archive/sdss/Simard11/table3.fit',1)
pobj  = mrdfits('/Volumes/archive/sdss/Simard11/photObj.fit',1)
CC,spobj.ra,spobj.dec,si11._ra,si11._de,0.5,1,dis=dis,mid=mid,/verb
si11_matched = struct_combine(si11[mid],pobj[mid])
; deal w/ non-matches
ind = where(dis lt 0 or dis gt 1)
help,ind 
zero_struct,si11_matched,ind
; save
mwrfits,si11_matched,'simard.fits',/create


end
