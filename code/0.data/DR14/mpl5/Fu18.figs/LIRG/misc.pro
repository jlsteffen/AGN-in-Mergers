;; Arp catalog
;info=QueryVizier('VII/192/arplist','NONE')
;stringad,info.raj2000+' '+info.dej2000,ra,dec
;; MaNGA catalog
;drpall = mrdfits('$MANGA_SPECTRO_REDUX/$MANGADRP_VER/drpall*.fits',1)
;CC,ra,dec,drpall.objra,drpall.objdec,lindgen(n_elements(drpall)),5,10,dis=dis,mid=mid,/verbose
;s = where(mid gt 0)
;forprint,info[s].arp,info[s].name,drpall[mid[s]].plateifu,f='i4,a10,a15'

; Arp#, Other Name, MaNGA ID
; 118  NGC 1143     8154-12702
; 243  NGC 2623     9507-12704
; 243  NGC 2623     9507-12704
;  55 UGC 04881     8250-12704
; 151  MRK 0040      9000-1901
; 183 UGC 08560     8985-12701
;  73   IC 1222     8484-12703
; 125 UGC 10491      8601-1902
; 310   IC 1259      8624-6101
; 323 NGC 7783B      8655-6104

; Sanders2003 catalog
info=QueryVizier('J/AJ/126/1607/table1','NONE',/all)
s = where(info.logir2 gt 11) ; limit to LIRG
info = info[s]
stringad,info.raj2000+' '+info.dej2000,ra,dec
; MaNGA catalog
drpall = mrdfits('$MANGA_SPECTRO_REDUX/$MANGADRP_VER/drpall*.fits',1)
CC,ra,dec,drpall.objra,drpall.objdec,lindgen(n_elements(drpall)),10,20,dis=dis,mid=mid,/verbose
s = where(mid gt 0)
forprint,info[s].IRAS,info[s].raj2000,info[s].dej2000,info[s].name,info[s].logir2,$
	dis[s],drpall[mid[s]].plateifu,f='a,1x,2a10,1x,a15,1x,f5.2,f5.1,a15'

; V2_2_0
; 08 38 23.8 +25 45 17        NGC 2623 11.54  2.9     9507-12704 = Arp 243 (coalesced)

; MPL-5
; F09126+4432 09 15 54.8 +44 19 58       UGC 04881 11.69  7.8 8250-12704 = Arp 55 (9.1 kpc)
; F15163+4255 15 18 06.7 +42 44 41          VV 705 11.89  7.4 7443-12703 = Mrk 848 (5.5 kpc)

end
