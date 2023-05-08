; load DRPALL
drpall = mrdfits('../../sample/matched/drpall.fits',1) ; 2718
drpall.plateifu = strtrim(drpall.plateifu,2)
; load targ catalog
target = mrdfits('../../sample/matched/target.fits',1) ; 2718

; coordinates in DRPALL
; OBJRA	float64	degrees	Right ascension of the science object in J2000
; IFURA	float64	degrees	Right ascension of this IFU in J2000
; CENRA	float64	degrees	Plate center right ascension in J2000

; coordinates in Target catalog
; CATALOG_RA	float64[1]	deg	Right Ascension of measured object center (J2000) as given in the input catalog (NSA for main samples and most ancillaries)
; IFU_RA	float64[1]	deg	The Right Ascension (J2000) of the IFU center
; OBJECT_RA	float64[1]	deg	The best estimate of the Right Ascension (J2000) of the center of the object. Normally the same as CATALOG_RA but can be modified particularly as a result of visual inspection
; TILERA	float64[1]	deg	The Right Ascension (J2000) of the tile to which this object has been allocated

s = where(target.NSA_NSAID ne 0) ; 2657, remove non-matches

; drpall.OBJRA = target.Object_RA
dis = sphdist(drpall.objra,drpall.objdec,target.object_ra,target.object_dec,/deg)
print,minmax(dis[s])*3600 ; < 0.213"

; drpall.OBJRA != target.catalog_RA
dis = sphdist(drpall.objra,drpall.objdec,target.catalog_ra,target.catalog_dec,/deg)
print,minmax(dis[s])*3600 ; < 11.5"

; drpall.IFURA != target.IFU_RA
dis = sphdist(drpall.ifura,drpall.ifudec,target.ifu_ra,target.ifu_dec,/deg)
print,minmax(dis[s])*3600 ; < 11.7

; offset between object and IFU centers
offset = sphdist(drpall.ifura,drpall.ifudec,drpall.objra,drpall.objdec,/deg)*3600
plothist,offset[s],/ylog

end
