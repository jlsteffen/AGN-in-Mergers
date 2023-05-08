; load Liu11 pair catalog - 2618 objects
readcol,'apj398354t1_mrt.txt',desig,plate,fiber,mjd,z,Delth,rp,DelV,rmag,FAGN,FTidal,f='a,i,i,i,f,f,f,i,f,i,i'
coord = strmid(desig,1,2)+' '+strmid(desig,3,2)+' '+strmid(desig,5,5)+' '+$
	strmid(desig,10,3)+' '+strmid(desig,13,2)+' '+strmid(desig,15)
stringad,coord,ra,dec

; select only BAGNs to cross-match 
pclass = mrdfits('$MANGA_DIR/merger/8.class.fits',1)
print,' BAGNs with strong W(Ha):'
ind = where(strtrim(pclass.class,2) eq 'BAGN')
CC,pclass[ind].ifura,pclass[ind].ifudec,ra,dec,lindgen(n_elements(ra)),10,20,dis=dis,mid=mid
s = where(dis ge 0)
forprint,'; '+pclass[ind[s]].plateifu,' '+desig[mid[s]],ra[mid[s]],dec[mid[s]]
; 4 MaNGA BAGNs are in Liu11 Catalog 
; 8083-9101  J032033.21-002024.0       50.138375     -0.34000000
; 8133-12704 J073906.17+442410.0       114.77571       44.402778
; 8612-12705 J170024.36+382106.3       255.10150       38.351750
; 9049-12701 J162628.71+240132.6       246.61962       24.025722

print,' BAGNs regardless of W(Ha):'
ind = where(pclass.bpt[0] ge 2 and pclass.bpt[1] ge 2) ; 61 pairs
CC,pclass[ind].ifura,pclass[ind].ifudec,ra,dec,lindgen(n_elements(ra)),10,20,dis=dis,mid=mid
s = where(dis ge 0)
forprint,'; '+pclass[ind[s]].plateifu,' '+desig[mid[s]],ra[mid[s]],dec[mid[s]]
; 5 MaNGA BAGNs are in Liu11 Catalog 
; 8083-9101  J032033.21-002024.0       50.138375     -0.34000000
; 8133-12704 J073906.17+442410.0       114.77571       44.402778
; 8549-12705 J160737.74+450355.2       241.90725       45.065333
; 8612-12705 J170024.36+382106.3       255.10150       38.351750
; 9049-12701 J162628.71+240132.6       246.61962       24.025722

; 38 AGNs in Liu11 pair catalog have MaNGA data
drpall = mrdfits('$MANGA_DIR/hfdap/$MANGADRP_VER/drpall.fits',1) ; 2718
CC,ra,dec,drpall.ifura,drpall.ifudec,lindgen(n_elements(drpall)),10,20,dis=dis,mid=mid
s = where(dis ge 0) ; 38 objects
forprint,'; '+drpall[mid[s]].plateifu,' '+desig[s],dis[s]
;-------------------------------
; both AGN in the same IFU - completely overlaps with the above
; 8083-9101  J032033.21-002024.0 BAGN 
; 8083-9101  J032033.65-002021.2 BAGN
; 8133-12704 J073905.83+442410.4 BAGN
; 8133-12704 J073906.17+442410.0 BAGN
; 8549-12705 J160737.24+450355.2 AGN-RG
; 8549-12705 J160737.74+450355.2 AGN-RG
; 9049-12701 J162628.06+240137.3 BAGN
; 9049-12701 J162628.71+240132.6 BAGN
; 8612-12705 J170024.36+382106.3 BAGN
; 8612-12705 J170024.76+382115.5 BAGN
; 7960-3702  J171309.43+313455.3 SF-RG
; 7960-3702  J171309.92+313452.3 SF-RG
;-------------------------------
; individual matches - mostly only one galaxy in IFU
; 8078-12703 J025016.86+000531.1
; 8083-9101  J032033.21-002024.0
; 8083-9101  J032033.65-002021.2
; 8082-6104  J032305.28+002115.4
; 8083-6103  J032345.63+000026.8
; 8133-12704 J073905.83+442410.4
; 8133-12704 J073906.17+442410.0
; 8140-3702  J074939.00+425801.4
; 8140-12703 J075135.63+425248.4
; 8717-9101  J075322.76+361627.8
; 8718-12701 J075643.72+445124.2
; 8143-3704  J080132.31+420008.4
; 8456-1901  J100452.02+464124.1
; 8568-9101  J102219.39+363458.9
; 8258-6102  J110824.92+430046.7
; 8950-12703 J125929.96+275723.2
; 8318-1901  J130701.21+461833.9
; 8320-6101  J134630.60+224221.7
; 9002-3702  J144845.74+314113.3
; 9042-1901  J153244.74+280347.3
; 8313-9102  J155957.11+412840.0
; 8549-12702 J160505.16+452634.8                               
; 8549-12705 J160737.24+450355.2
; 8549-12705 J160737.74+450355.2
; 9049-12701 J162628.06+240137.3
; 9049-12701 J162628.71+240132.6
; 9029-12704 J162852.06+424843.2
; 8604-12703 J163103.41+395018.5
; 8550-6104  J163339.16+391527.5
; 8588-6101  J163349.62+391547.5
; 8612-12705 J170024.36+382106.3
; 8612-12705 J170024.76+382115.5
; 7960-3702  J171309.43+313455.3
; 7960-3702  J171309.92+313452.3
; 7962-9102  J172043.95+284229.5
; 8611-6101  J172752.52+600550.1
; 8611-3704  J173159.22+595818.0
; 8615-6104  J211936.44+004215.2

end
