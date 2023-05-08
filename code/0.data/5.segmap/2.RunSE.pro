; run sextractor and generate segmentation maps for all images

inpdir = 'rband/'
outdir = 'rsex_cleanon/'

infiles = file_basename(file_search(inpdir+'7815-6102.fits',count=ct))

for i=0,ct-1 do begin

cmd = 'sex '+inpdir+infiles[i]+' -c config/sdss_segmap.sex'+$
	' -CATALOG_NAME '+outdir+repstr(infiles[i],'.fits','.cat')+$
	' -CLEAN YES'+$
	' -CHECKIMAGE_TYPE SEGMENTATION'+$
	' -CHECKIMAGE_NAME '+outdir+repstr(infiles[i],'.fits','_seg.fits')

spawn,cmd

;cs = rsex(outdir+repstr(infiles[i],'.fits','.cat'))
;help,cs,/st

endfor

end
