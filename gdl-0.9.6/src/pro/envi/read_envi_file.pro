;+
;
; NAME: 
;     READ_ENVI_FILE
; PURPOSE: 
;     Reads an ENVI style formatted image file (both compressed and 
;     uncompressed). 
;
; CATEGORY:
;     File I/O
;
; CALLING SEQUENCE:
;     READ_ENVI_FILE, filename, image=image, info=image_info
;
; KEYWORD PARAMETERS:
;     FILENAME: A string containing the full filepath name of the image.
;     IMAGE: Variable name to contain the image/array
;     INFO: Variable name to contain the image information as a structure
;     HELP: showing how to use and exit
;
; OUTPUTS:
;     image: An array containing the image 
;     image_info: A structure containing information about the image
;                 eg. number of samples, number of lines, number of bands, 
;                 data type, interleave, file description, band names. 
;
; MODIFICATION HISTORY:
; 02-Jul-2012: Written by Josh Sixsmith
; 16-Feb-2014: Better managment of input filenames
;
; LICENCE:
; Copyright (C) 2012, Josh Sixsmith
; This program is free software; you can redistribute it and/or modify  
; it under the terms of the GNU General Public License as published by  
; the Free Software Foundation; either version 3 of the License, or     
; (at your option) any later version.                                   
;
;-
;
function NUM_SAMPLES, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fnsw = WHERE(STRPOS(header, "samples") NE -1, count)
IF (count NE 0) THEN BEGIN
   fns  = STRTRIM(header[fnsw], 2)
   fns1 = STRPOS(fns, '=')
   ns   = STRMID(fns, fns1 + 2)
   ;;
   RETURN, (LONG(ns))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'Number of Samples Not Found.'
ENDELSE
;
END
;
;-----------------------------------------------------------------------
;
function NUM_LINES, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fnlw = WHERE(STRPOS(header, "lines") NE -1, count)
IF (count NE 0) THEN BEGIN
   fnl  = STRTRIM(header[fnlw], 2)
   fnl1 = STRPOS(fnl, '=')
   nl   = STRMID(fnl, fnl1 + 2)
   ;;
   RETURN, (LONG(nl))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'Number of Lines Not Found.'
ENDELSE
;
END
;
;-----------------------------------------------------------------------
;
function NUM_BANDS, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fnbw = WHERE(STRPOS(header, "bands") NE -1, count)
IF (count NE 0) THEN BEGIN
   fnb = STRTRIM(header[fnbw], 2)
   fnb1 = STRPOS(fnb, '=')
   nb = STRMID(fnb, fnb1 + 2)
   ;;
   RETURN, (LONG(nb))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'Number of Bands Not Found.'
ENDELSE
;
END
;
;-----------------------------------------------------------------------
;
function MAP_INFO, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fmiw = WHERE(STRPOS(header, "map") NE -1, count)
IF (count NE 0) THEN BEGIN
   fmi = STRTRIM(header[fmiw], 2)
   b1 = STRPOS(fmi, '{')
   b2  = STRPOS(fmi, '}', /reverse_search)
   b2b = STRSPLIT(STRMID(fmi, b1+1, b2-b1-1), ',', /extract)
   mi  = {MapInfo, $
          ProjectionName : b2b[0], $
          ULRefX : b2b[1], $
          ULRefY : b2b[2],$
          ULXCoord : b2b[3],$
          ULYCoord : b2b[4], $
          XS : b2b[5], $
          YS : b2b[6], $
          Zone : b2b[7], $
          Units : b2b[8]}
ENDIF ELSE BEGIN
   mi = {map_info, empty:''}
ENDELSE
;
RETURN, mi
;
END
;
;-----------------------------------------------------------------------
;
FUNCTION DATATYPE, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
;; Maybe need a check for correct data type, using filesize and array
;; dimensions before trying to read the data.
;
fdtw = WHERE(STRPOS(header, "data type") NE -1, count)
IF (count NE 0) THEN BEGIN
   fdt =  STRTRIM(header[fdtw], 2)
   fdt1 = STRPOS(fdt, '=')
   dtype   = STRMID(fdt, fdt1 + 2)
   ;;
   RETURN, (FIX(dtype))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'No Data Type Found.'
ENDELSE
;
END
;
;-----------------------------------------------------------------------
;
FUNCTION interleave, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
filw = WHERE(STRPOS(header, "interleave") NE -1, count)
IF (count NE 0) THEN BEGIN
   fil = STRTRIM(header[filw], 2)
   fil1 = STRPOS(fil, '=')
   ileave   = STRMID(fil, fil1 + 2)
   ;;
   IF (STRCMP(ileave, 'BSQ', /fold_case)) EQ 1 THEN BEGIN
      intleave = 0
   ENDIF ELSE BEGIN
      IF (strcmp(ileave, 'BIL', /fold_case)) EQ 1 THEN BEGIN
         intleave = 1
      ENDIF ELSE BEGIN
         IF (STRCMP(ileave, 'BIP', /fold_case)) EQ 1 THEN BEGIN
            intleave = 2
         ENDIF ELSE BEGIN
            MESSAGE, 'Unknown Interleaving; Need either BSQ/BIL/BIP.'
         ENDELSE
      ENDELSE
   ENDELSE
   ;;
ENDIF ELSE BEGIN
   MESSAGE, 'No Interleave Found, Assuming BSQ.', /continue
   intleave = 0
ENDELSE
;
RETURN, intleave
;
END
;
;-----------------------------------------------------------------------
;
function READ_HEADER, filename
;
COMPILE_OPT hidden
ON_ERROR, 2
;
OPENR, lun, filename, /get_lun
array = ''
line = ''
;
WHILE NOT EOF(lun) DO BEGIN
   READF, lun, line
   array = [array, line]
ENDWHILE
;
FREE_LUN, lun
RETURN, array[1:*]
;
END
;
;-----------------------------------------------------------------------
;
function SENSOR_TYPE, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fstw = WHERE(STRPOS(header, "sensor type") NE -1, count)
IF (count NE 0) THEN BEGIN
   fst = STRTRIM(header[fstw], 2)
   fst1 = STRPOS(fst, '=')
   stype   = (STRMID(fst, fst1 + 2))[0]
   ;;
ENDIF ELSE BEGIN
   stype = "Unknown"
ENDELSE
;
RETURN, stype
END
;
;-----------------------------------------------------------------------
;
function WAVELENGTH_UNITS, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fwuw = WHERE(STRPOS(header, "wavelength units") NE -1, count)
IF (count NE 0) THEN BEGIN
   fwu = STRTRIM(header[fwuw], 2)
   fwu1 = STRPOS(fwu, '=')
   wvunits   = (STRMID(fwu, fwu1 + 2))[0]
   ;;
ENDIF ELSE BEGIN
   wvunits = "Unknown"
ENDELSE
;
RETURN, wvunits
END
;
;-----------------------------------------------------------------------
;
function BAND_NAMES, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fbnw = WHERE(STRPOS(header, "band names") NE -1, count)
IF (count NE 0) THEN BEGIN
   fbn = STRTRIM(header[fbnw], 2)
   fbn1 = STRPOS(fbn, '{')
   ;;
   eb_array = STRPOS(header[fbnw+1:*], '}')
   names = ''
   FOR i = 1, N_ELEMENTS(eb_array) DO BEGIN
      names = names + header[fbnw+i]
   ENDFOR
   eb = STRPOS(names, '}')
   names  = STRTRIM(STRMID(names, 0, eb), 2)
   b_names = STRSPLIT(names, ',', /extract)
   ;;
ENDIF ELSE BEGIN
   nb   = num_bands(header)
   band = STRING(LONARR(nb))
   number=STRING(LONARR(nb))
   b_names=STRING(LONARR(nb))
   ;;
   ;;create the array with value 'Band' placed in each element
   FOR i=0L, nb[0]-1 DO BEGIN
      band[i]= 'Band '
   ENDFOR
   ;;
   ;;create the array with values of 1 to the total number of files
   FOR i=0L, nb[0]-1 DO BEGIN
      number[i]= STRTRIM(i+1,1)
   ENDFOR
   ;;
   ;;concatenate (join) the band and number arrays into one singular array
   FOR i=0L, nb[0]-1 DO BEGIN
      b_names[i]= band[i] + number[i]
   ENDFOR
ENDELSE
;
RETURN, b_names
END
;
;-----------------------------------------------------------------------
;
function BYTE_ORDER, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fbow = WHERE(STRPOS(header, "byte order") NE -1, count)
IF (count NE 0) THEN BEGIN
   fbo = STRTRIM(header[fbow], 2)
   fbo1 = STRPOS(fbo, '=')
   byt_order   = (FIX(STRMID(fbo, fbo1 + 2)))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'Byte order not found, assuming byte order of current machine.', /continue
   byt_order = (BYTE(1,0,1))[0] ? 0 : 1
ENDELSE
;
RETURN, byt_order
;
END
;
;-----------------------------------------------------------------------
;
function HEADER_OFFSET, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fhow = WHERE(STRPOS(header, "header offset") NE -1, count)
IF (count NE 0) THEN BEGIN
   fho = STRTRIM(header[fhow], 2)
   fho1 = STRPOS(fho, '=')
   offset = (LONG(STRMID(fho, fho1 + 2)))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'No offset found, assuming zero.', /continue
   offset = 0L
ENDELSE
;
RETURN, offset
;
END
;
;-----------------------------------------------------------------------
;
function DESCRIPTION, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fdsw = WHERE(STRPOS(header, "description") NE -1, count)
IF (count NE 0) THEN BEGIN
   fds = STRTRIM(header[fdsw], 2)
   fds1 = STRPOS(fds, '{')
   ;;
   ;; Using the same method as for band names. It seems to work fine.
   eb_array = STRPOS(header[fdsw+1:*], '}')
   desc = ''
   FOR i = 1, N_ELEMENTS(eb_array) DO BEGIN
      desc = desc + header[fdsw+i]
   ENDFOR
   eb = STRPOS(desc, '}')
   descrip  = (STRTRIM(STRMID(desc, 0, eb), 2))[0]
   ;;
ENDIF ELSE BEGIN
   descrip = 'None'
ENDELSE
;
RETURN, descrip
;
END
;
;-----------------------------------------------------------------------
;
FUNCTION FILETYPE, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
fftw = WHERE(STRPOS(header, "file type") NE -1, count)
IF (count NE 0) THEN BEGIN
   ff_t = STRTRIM(header[fftw], 2)
   fft1 = STRPOS(ff_t, '=')
   ftype = (STRMID(ff_t, fft1 + 2))[0]
ENDIF ELSE BEGIN
   MESSAGE, 'File type not found, assumed ENVI Standard.', /continue
   ftype = 'ENVI Standard'
ENDELSE
;
RETURN, ftype
;
END
;
;-----------------------------------------------------------------------
;
function F_COMPRESSION, header
;
COMPILE_OPT hidden
ON_ERROR, 2
;
ffcw = WHERE(STRPOS(header, "file compression") NE -1, count)
IF (count NE 0) THEN BEGIN
   ffc  = STRTRIM(header[ffcw], 2)
   ffc1 = STRPOS(ffc, '=')
   fc   = STRMID(ffc, ffc1 + 2)
   rfc   = FIX(fc[0])
ENDIF ELSE BEGIN
   rfc   = 0
ENDELSE
;
RETURN, rfc
;
END
;
; -----------------------------------------------------
; better managing finding the files (.IMG + .HDR or .IMG.HDR + .IMG)
;
function ENVI_SELECT_FILENAME, filename, hname, fname, test=test, debug=debug
;
suffixe=STRMID(filename, 3, /reverse_offset)
files_found=0
;
; testing first case: HDR file as input.
;
IF (STRCMP(suffixe, '.hdr', /fold_case)) EQ 1 THEN BEGIN
    ;; path
    path=FILE_DIRNAME(filename, /mark_directory)
    if path EQ './' then path=''
    ;; 
    bodyname=FILE_BASENAME(filename, '.hdr')
    fname=FILE_SEARCH(path+bodyname+'.img',/fold)
    if (STRLEN(fname) EQ 0) then begin
        mess='No IMG data file found corresponding to HDR header file : '
        MESSAGE, mess+filename
    endif
    if (N_ELEMENTS(fname) GT 1) then begin
        MESSAGE,/continue, 'More than one IMG data file found, first used !'
        fname=fname[0]
    endif
    hname=filename
    files_found=1
endif
;
; testing first case: IMG file as input.
;
IF (STRCMP(suffixe, '.img', /fold_case)) EQ 1 THEN BEGIN
    ;; path
    path=FILE_DIRNAME(filename, /mark_directory)
    if path EQ './' then path=''
    ;; 
    bodyname=FILE_BASENAME(filename, '.img')
    ;; we use '*.hdr' because we may have '.img.hdr' or '.hdr'
    hname=FILE_SEARCH(path+bodyname+'*.hdr',/fold)
    if (STRLEN(hname) EQ 0) then begin
        mess='No HDR header file found corresponding to IMG data file : '
        MESSAGE, mess+filename
    endif
    if (N_ELEMENTS(hname) GT 1) then begin
        MESSAGE,/continue, 'More than one HDR header file found, first used !'
        hname=hname[0]
    endif
    fname=filename
    files_found=1
endif
;
; if just the "file_basename" is provide, do we have the files ?
;
if (files_found EQ 0) then begin
    ;; path
    path=FILE_DIRNAME(filename, /mark_directory)
    if path EQ './' then path=''
    ;; 
    hname=FILE_SEARCH(filename+'*.hdr',/fold)
    fname=FILE_SEARCH(filename+'*.img',/fold)
    ;
    if (N_ELEMENTS(hname) EQ 0) then $
      MESSAGE, /cont, 'no HDR header file found with given basename pattern'
    ;;
    if (N_ELEMENTS(fname) EQ 0) then $
      MESSAGE, /cont, 'no IMG data file found with given basename pattern'
    ;;
    if ((N_ELEMENTS(hname) EQ 1) AND (N_ELEMENTS(fname) EQ 1)) then begin
        b1=FILE_BASENAME(hname,'.hdr',/fold)
        b1bis=FILE_BASENAME(hname,'.img.hdr',/fold)
        b2=FILE_BASENAME(fname,'.img',/fold)
        if STRCMP(b1,b2,/fold_case) EQ 1 then files_found=1
        if STRCMP(b1bis,b2,/fold_case) EQ 1 then files_found=1
        if files_found EQ 0 then begin
            MESSAGE, /cont, 'bad files names, please check the input filter'
        endif
    endif
    ;;
    txt='found, please check input filename'
    if (N_ELEMENTS(hname) GT 1) then begin
        MESSAGE, /cont, 'More than one HDR header file '+txt
    endif
    if (N_ELEMENTS(fname) GT 1) then begin
        MESSAGE, /cont, 'More than one IMG data file '+txt
    endif
endif
;
; Are really the files around ? Are the files void ?
;
if (files_found EQ 1) then begin
    count=0
    if ~FILE_TEST(fname) then begin
        count++
        MESSAGE, /continue, 'IMG data file not found (no file/bad name ?)'
    endif else begin
        if FILE_TEST(fname,/zero_length) then begin
            count++
            MESSAGE, /continue, 'IMG data file does not contain data !'
        endif
    endelse
    if ~FILE_TEST(hname) then begin
        count++
        MESSAGE, /continue, 'HDR header file not found (no file/bad name ?)'
    endif else begin        
        if FILE_TEST(hname,/zero_length) then begin
            count++
            MESSAGE, /continue, 'HDR header file does not contain data !'
        endif
    endelse
    if (count NE 0) then files_found=0
endif
;
if (files_found GT 0) then begin
    PRINT, 'IMG data filename   = ', fname
    PRINT, 'HDR header filename = ', hname
    PRINT, 'Directory name      = ', FILE_DIRNAME(hname, /mark_directory)
endif
;
if KEYWORD_SET(test) then STOP
;
return, files_found
;
end
;
;-----------------------------------------------------------------------
;
pro READ_ENVI_FILE, filename, image=image, info=info, $
                    help=help, test=test, debug=debug
;
if ~KEYWORD_SET(debug) then ON_ERROR, 2
;
if KEYWORD_SET(help) THEN BEGIN
   PRINT, 'pro TEST_READ_ENVI, filename, image=image, info=info, $'
   PRINT, '                    help=help, test=test, debug=debug'
   PRINT, ''
   PRINT, 'Reads an ENVI style image format. Set info to a variable'
   PRINT, 'that will contain a structure containing the image'
   PRINT, 'information such as samples, lines, bands, data type etc.'
   PRINT, ''
   PRINT, 'Warning : currently, suffixes can be only .IMG for data files and'
   PRINT, '(.IMG.HRD or .HDR) for Header files [with any Up/Low case combi]'
   return
ENDIF
;
txt='FILENAME (.HDR, .IMG or basename)'
;
if N_PARAMS() EQ 0 then filename=DIALOG_PICKFILE(filter='*.hdr')
if N_ELEMENTS(filename) GT 1 then $
  MESSAGE, 'You can provide only one '+txt+' at once'
;
; from a given file name (w/o extension .hrd/.img + w/o path)
; returning ONE and only one pair of Hname+Fname 
;
files_found=ENVI_SELECT_FILENAME(filename, hname, fname, test=test)
if ~files_found then MESSAGE, 'Please check carrefully the input '+txt
;
; after this point, we have one IMG file and one HDR file 
; (existances checked)
;
if KEYWORD_SET(debug) then STOP
;
hdr = READ_HEADER(hname)
;
;; Get the description info
descript = DESCRIPTION(hdr)
;
;; Get samples, lines, bands, datatype, interleave
ns = NUM_SAMPLES(hdr)
nl = NUM_LINES(hdr)
nb = NUM_BANDS(hdr)
;
dtype    = DATATYPE(hdr)
intleave = INTERLEAVE(hdr)
;
CASE intleave OF
   ;; BSQ Interleaving
   0: BEGIN
      CASE dtype OF
         0: MESSAGE, 'Undefined Data Type.'
         1: image = BYTARR(ns, nl, nb, /nozero)
         2: image = INTARR(ns, nl, nb, /nozero)
         3: image = LONARR(ns, nl, nb, /nozero)
         4: image = FLTARR(ns, nl, nb, /nozero)
         5: image = DBLARR(ns, nl, nb, /nozero)
         6: image = COMPLEXARR(ns, nl, nb, /nozero)
         7: image = STRARR(ns, nl, nb)
         8: MESSAGE, 'Structured Arrays Not Supported.'
         9: image = DCOMPLEXARR(ns, nl, nb, /nozero)
         10: image = PTRARR(ns, nl, nb, /nozero)
         11: image = OBJARR(ns, nl, nb, /nozero)
         12: image = UINTARR(ns, nl, nb, /nozero)
         13: image = ULONARR(ns, nl, nb, /nozero)
         14: image = LON64ARR(ns, nl, nb, /nozero)
         15: image = ULON64ARR(ns, nl, nb, /nozero)
         ELSE: MESSAGE, 'Unknown Data Type.'
      ENDCASE
      BREAK
   END
                                ; BIL Interleaving
   1: BEGIN
      CASE dtype OF
         0: MESSAGE, 'Undefined Data Type.'
         1: image = BYTARR(ns, nb, nl, /nozero)
         2: image = INTARR(ns, nb, nl, /nozero)
         3: image = LONARR(ns, nb, nl, /nozero)
         4: image = FLTARR(ns, nb, nl, /nozero)
         5: image = DBLARR(ns, nb, nl, /nozero)
         6: image = COMPLEXARR(ns, nb, nl, /nozero)
         7: image = STRARR(ns, nb, nl)
         8: MESSAGE, 'Structured Arrays Not Supported.'
         9: image = DCOMPLEXARR(ns, nb, nl, /nozero)
         10: image = PTRARR(ns, nb, nl, /nozero)
         11: image = OBJARR(ns, nb, nl, /nozero)
         12: image = UINTARR(ns, nb, nl, /nozero)
         13: image = ULONARR(ns, nb, nl, /nozero)
         14: image = LON64ARR(ns, nb, nl, /nozero)
         15: image = ULON64ARR(ns, nb, nl, /nozero)
         ELSE: MESSAGE, 'Unknown Data Type.'
      ENDCASE
      BREAK
   END
                                ; BIP Interleaving
   2: BEGIN
      CASE dtype OF
         0: MESSAGE, 'Undefined Data Type.'
         1: image = BYTARR(nb, ns, nl, /nozero)
         2: image = INTARR(nb, ns, nl, /nozero)
         3: image = LONARR(nb, ns, nl, /nozero)
         4: image = FLTARR(nb, ns, nl, /nozero)
         5: image = DBLARR(nb, ns, nl, /nozero)
         6: image = COMPLEXARR(nb, ns, nl, /nozero)
         7: image = STRARR(nb, ns, nl)
         8: MESSAGE, 'Structured Arrays Not Supported.'
         9: image = DCOMPLEXARR(nb, ns, nl, /nozero)
         10: image = PTRARR(nb, ns, nl, /nozero)
         11: image = OBJARR(nb, ns, nl, /nozero)
         12: image = UINTARR(nb, ns, nl, /nozero)
         13: image = ULONARR(nb, ns, nl, /nozero)
         14: image = LON64ARR(nb, ns, nl, /nozero)
         15: image = ULON64ARR(nb, ns, nl, /nozero)
         ELSE: MESSAGE, 'Unknown Data Type.'
      ENDCASE
      BREAK
   END
   ELSE: MESSAGE, 'Unknown Interleave.'
ENDCASE
;
;; byte order of file
f_byte_order = BYTE_ORDER(hdr)
;
;; file type
ftype = FILETYPE(hdr)
;
;; byte order of current machine
m_byte_order = (BYTE(1,0,1))[0] ? 0 : 1
;
;; header offset
h_offset = HEADER_OFFSET(hdr)
;
;; file compression
compressed = F_COMPRESSION(hdr)
;
IF (f_byte_order NE m_byte_order) THEN BEGIN
   IF compressed NE 0 THEN BEGIN
      OPENR, lun, fname, /get_lun, /swap_endian, /compress
   ENDIF ELSE BEGIN
      OPENR, lun, fname, /get_lun, /swap_endian
   ENDELSE
   POINT_LUN, lun, h_offset
   READU, lun, image
   FREE_LUN, lun
ENDIF ELSE BEGIN
   ;; Open and read the image file
   IF compressed NE 0 THEN BEGIN
      OPENR, lun, fname, /get_lun, /compress
   ENDIF ELSE BEGIN
      OPENR, lun, fname, /get_lun
   ENDELSE
   POINT_LUN, lun, h_offset
   READU, lun, image
   FREE_LUN, lun
ENDELSE
;
stype = SENSOR_TYPE(hdr)
wl_units = WAVELENGTH_UNITS(hdr)
bnames = BAND_NAMES(hdr)

info = {Samples:ns, Lines:nl, Bands:nb, Data_Type:dtype,$
        Interleave:intleave, Filename:fname, Sensor_Type:stype, $
        Wavelength_Units:wl_units, Band_Names:bnames, $
        Byte_Order:f_byte_order, Header_offset:h_offset, $
        Description:descript, File_Type:ftype}
;
if KEYWORD_SET(test) then STOP
;
END    
