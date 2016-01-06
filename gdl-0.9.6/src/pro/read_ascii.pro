;+
; NAME: READ_ASCII
;
;
; PURPOSE: Reads an ASCII file. The output is a structure whose tags contains
;          columns (1D array) from the file. Its use is flexible: the user can 
;          specify the types, and tag names of the fields and he/she also can 
;          group columns into a single tag (2D array).
;
;
; CATEGORY: IO
;
;
; CALLING SEQUENCE: 
;   structure=read_ascii(filename, count=, data_start=, delimiter=,
;                        missing_value=, comment_symbol=, record_start=
;                        num_records=, template=, header=, verbose=)
;
;
; INPUT:
;   filename       Name of the ASCII file to be read
;
;
; KEYED INPUTS:
;   data_start     Specify the number of lines that constitute the header
;                  These lines will be discarded as records, but are available
;                  with the header keyword
;
;   delimiter      If set (or non equal to ''), the records will be split 
;                  according to the supplied delimiter locations and the length
;                  of the fields is not necessarily the same for all records. 
;                  Otherwise, template.fieldlocations will be used to identify 
;                  the fields (columns). If the template is not provided, the 
;                  delimiter default is ' '
;
;   missing_value  Specify the value that will be used for missing or 
;                  non-numeric values
;
;   comment_symbol Comment symbol. The part of the line that begins with
;                  the comment symbol and ends to the end of line is discarded
;                  Default is ';'
;
;   record_start   Record number of the first record to be read (starts from 0)
;
;   num_record     Number of records to be read
;
;   template       structure that defines how the file will be processed
;                  the tags datastart, delimiter, missingvalue and commentsymbol
;                  can be overridden by the keywords data_start, delimiter,
;                  missing_value and comment_symbol.
;
;   template.VERSION           template version number (not used)
;   template.FIELDCOUNT=n      number of fields
;   template.FIELDNAMES[n]     field names
;   template.FIELDTYPES[n]     field integer types
;   template.FIELDGROUPS[n]    group ID. fields can be grouped into a single tag
;   template.FIELDLOCATIONS[n] start positions of the fields
;   template.DATASTART         see data_start keyed input
;   template.DELIMITER         see delimiter keyed input
;   template.MISSINGVALUE      see missing_value keyed input
;   template.COMMENTSYMBOL     see comment_symbol keyed input
;
;   if field of different types are grouped together, the following priority
;   if assumed to determine the tag type:
;     BYTE<UINT<INT<ULONG<LONG<ULONG64<LONG64<FLOAT<DOUBLE<COMPLEX<DCOMPLEX
;   Exception: grouping double and COMPLEX will result in a DCOMPLEX
;   No boundary check is made, so information could be lost.
;
;
; KEYWORDS:
;   verbose        Not used.
;
;
; OUTPUT:
;   structure      Structure containing the file columns.
;
;
; KEYED OUTPUTS:
;   count          Number of records that have beeen read
;   header         First data_start lines of the file
;
;
; EXAMPLES:
; t = {version:1.0, fieldnames : STRSPLIT('fa,fb,fc,fd,fe,ff',',', /extr), $
;      fieldtypes : [7, 4, 7, 2, 1, 5], fieldgroups : [0, 1, 2, 3, 4, 5], $
;      fieldcount: 6, fieldlocations:[0, 5, 9, 11, 14, 16], datastart:0, $
;      delimiter:'', missingvalue:-999, commentsymbol:';'}
; a = read_ascii('test-read_ascii.txt', template=t, header=header, 
;                data_start=2, count=count)
; 
;
; RESTRICTIONS:
;   no boundary check is performed if the missing value can not be represented
;   by the field type.
;
;
; MODIFICATION HISTORY:
;   13-Jan-2006 : written by Pierre Chanial
;   06-APr-2008 : m_schellens: made data_start independent of header
;   15-Nov-2011 : A. Coulais : better management of dir/file and
;                 missing file
;   05-Feb-2014 : G. Duvert : avoid unlawful tag names 
;
;-
; LICENCE:
; Copyright (C) 2006, P. Chanial; 2011 A. Coulais
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
;
;-
;
pro READ_ASCII_HELPER, tags, tag, structure, variable, default
;
COMPILE_OPT hidden

if N_ELEMENTS(variable) ne 0 then return
if N_ELEMENTS(tags) eq 0 then begin
   variable = default
   return
endif
index = (WHERE(tags eq tag, count))[0]
if count eq 0 then begin
   variable = default
endif else begin
   variable = structure.(index)
endelse;
end
;
; -----------------------------------
;
function READ_ASCII_READ, filename
;
COMPILE_OPT hidden
ON_ERROR, 2
;
if FILE_TEST(filename, /directory) then $
   MESSAGE, 'this is not a file but a directory !'
if not FILE_TEST(filename) then MESSAGE, 'Invalid filename'
;
nline = FILE_LINES(filename) 
text = STRARR(nline)
OPENR, unit, filename, /GET_LUN
READF, unit, text
FREE_LUN, unit
return, text
; 
end
;
; -----------------------------------
;
function READ_ASCII_GETTYPE, types

compile_opt hidden
on_error, 2

priority = [ 0,  $              ; undefined
             1,  $              ; byte
             3,  $              ; int
             5,  $              ; long
             8,  $              ; float
             9,  $              ; double
             10, $              ; COMPLEX
             99, $              ; STRING
             0,  $              ; struct
             11, $              ; DCOMPLEX
             0,  $              ; pointer
             0,  $              ; object
             2,  $              ; uint
             4,  $              ; ulong
             7,  $              ; long64
             6   $              ; ulong64
           ]
 
p = MAX(priority[types], itype, min=status)
if (status eq 0) then MESSAGE, 'Invalid field type.'
type = types[itype]
;
; special case: double+COMPLEX -> DCOMPLEX
if MAX(types eq 5) and MAX(types eq 6) and (type ne 7) then type = 11
;
return, type
;
end
;
; ----------------------------------
;

function READ_ASCII, filename, count=linecount, $
                     data_start=data_start, delimiter=delimiter, $
                     missing_value=missing_value, $
                     comment_symbol=comment_symbol, $
                     num_records=num_records, record_start=record_start,   $
                     template=template, header=header, $
                     help=help, test=test, verbose=verbose
;
ON_ERROR, 2
;
if KEYWORD_SET(help) then begin
   print, 'function READ_ASCII, filename, count=linecount, $'
   print, '                     data_start=data_start, delimiter=delimiter, $'
   print, '                     missing_value=missing_value, $'
   print, '                     comment_symbol=comment_symbol, $
   print, '                     num_records=num_records, record_start=record_start,   $'
   print, '                     template=template, header=header, $'
   print, '                     help=help, test=test, verbose=verbose'
   return, -1
endif
;
if (N_PARAMS() EQ 0) then begin
; ~(FILE_INFO(filename)).exits then begin
   filename=DIALOG_PICKFILE()
endif
;
if SIZE(template, /tname) eq 'STRUCT' then begin
   tags = tag_names(template)
   deldefault = ''
endif else begin
   deldefault = ' 	'
endelse
READ_ASCII_HELPER, tags, 'DATASTART',     template, data_start,     0
READ_ASCII_HELPER, tags, 'DELIMITER',     template, delimiter,      deldefault
READ_ASCII_HELPER, tags, 'MISSINGVALUE',  template, missing_value,  'NaN'
READ_ASCII_HELPER, tags, 'COMMENTSYMBOL', template, comment_symbol, ';'
;
text=READ_ASCII_READ(filename)
;
;----------------
; get the header
;----------------
;
if ARG_PRESENT(header) then begin
   if data_start eq 0 then begin
      if N_ELEMENTS(header) ne 0 then junk = TEMPORARY(header)
   endif else begin
      header = text[0:data_start-1]
      ;; text = text[data_start:*]
   endelse
endif
;
if data_start ne 0 then begin
   if data_start gt N_ELEMENTS(text) - 1 then begin
      MESSAGE, 'DATA_START value >= data length (' $
               + STRTRIM(STRING(data_start), 2) + ' >= ' $
               + STRTRIM(STRING(N_ELEMENTS(text)),2) + ')' 
   endif
   text = text[data_start:*]
endif
;
;-----------------
; remove comments
;-----------------
;
if KEYWORD_SET(comment_symbol) then begin
   pos = STRPOS(text, comment_symbol)
   index = WHERE(pos ne -1, count)
   for i=0, count-1 do begin
      j = index[i]
      text[j] = STRMID(text[j], 0, pos[j])
   endfor
endif

;--------------------
; remove blank lines
;--------------------

text = STRTRIM(text, 2)
index = WHERE(STRLEN(text) ne 0, linecount)
if linecount eq 0 then return, 0
text = text[index]

;---------------------------
; get the requested records
;---------------------------

if KEYWORD_SET(record_start) then begin
   if record_start ge linecount then begin
      MESSAGE, 'Invalid record start ('+STRTRIM(record_start,1)+$
               '): the file only has '+STRTRIM(linecount,1)+' lines.'
   endif
   text = text[record_start:*]
   linecount = N_ELEMENTS(text)
endif

if KEYWORD_SET(num_records) then begin
   if num_records gt linecount then begin
      MESSAGE, 'Excessive number of requested records.'
   endif
   text = text[0:num_records-1]
   linecount = num_records
endif

rnumber = '^[+-]?([0-9]*\.?[0-9]*[ed]?[+-]?[0-9]+\.?[0-9]*|NaN|Inf|Infinity)$'

;------------------
; no-template case
;------------------

if N_ELEMENTS(template) eq 0 then begin
   ncolumns = LONARR(linecount)
   for line=0l, linecount-1 do begin
      ncolumns[line] = N_ELEMENTS(STRSPLIT(text[line], delimiter))
   endfor
   ncolumn = MAX(ncolumns)
   result = MAKE_ARRAY(ncolumn, linecount, /float, value=FLOAT(missing_value))
   for line=0l, linecount-1 do begin
      row = STRSPLIT(text[line], delimiter, /extract)
      index = WHERE(STREGEX(row, rnumber, /fold_case, /boolean), count)
      if count gt 0 then result[index, line] = float(row[index])
   endfor
   return, {field1:TEMPORARY(result)}
endif
;
; should take into account the field keyword, when RSI implements it.
;
fieldcount  = template.fieldcount
fieldtypes  = template.fieldtypes
fieldnames  = template.fieldnames
;fieldnames must not begin with a Number; blanks will be replaced by '_'. fieldnames must not be reserved words.
RWORDS=['AND','BEGIN','BREAK','CASE','COMMON','COMPILE_OPT','CONTINUE','DO',$
'ELSE','END','ENDCASE','ENDELSE','ENDFOR','ENDFOREACH','ENDIF','ENDREP',$
'ENDSWITCH','ENDWHILE','EQ','FOR','FOREACH','FORWARD_FUNCTION','FUNCTION',$
'GE','GOTO','GT','IF','INHERITS','LE','LT','MOD','NE','NOT','OF','ON_IOERROR',$
'OR','PRO','REPEAT','SWITCH','THEN','UNTIL','WHILE','XOR'] 
for ifield=0L,n_elements(fieldnames)-1 do begin
 if total(strmatch(RWORDS,fieldnames[ifield],/FOLD_CASE)) ne 0 then Message,'Illegal field name: '+fieldnames[ifield]
 ;unblank blanks:
  b=byte(fieldnames[ifield])
  for ibyte=0L,n_elements(b)-1 do if b[ibyte] eq 32 then b[ibyte]=95
  fieldnames[ifield]=string(b)
 ; test unconsistencies:
 if stregex(fieldnames[ifield],'[0123456789]') eq 0 then Message,'Illegal field name: '+fieldnames[ifield]
 if stregex(fieldnames[ifield],'[^ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_]',/FOLD_CASE) ne -1 then Message,'Illegal field name: '+fieldnames[ifield]
end
fieldlocs   = template.fieldlocations
fieldgroups = template.fieldgroups

strresult = STRARR(fieldcount, linecount)


;-------------------------------------
; slice the file content into columns
;-------------------------------------

if KEYWORD_SET(delimiter) then begin
   for line=0l, linecount-1 do begin
      row = STRSPLIT(text[line], STRING(delimiter), /extract)
      strresult[0, line] = row
   endfor
endif else begin
   for i=0l, fieldcount-2 do begin
      strresult[i,*] = STRMID(text, fieldlocs[i], fieldlocs[i+1]-fieldlocs[i])
   endfor
   strresult[i,*] = STRMID(text, fieldlocs[i])
endelse
strresult = STRTRIM(strresult,2)


;---------------------------
; get output structure info
;---------------------------

tagcount  = N_ELEMENTS(UNIQ(fieldgroups, SORT(fieldgroups)))
tagncols  = LONARR(tagcount)
tagnames  = STRARR(tagcount)
tagtypes  = INTARR(tagcount)
taggroups = LONARR(tagcount)

; get the group IDs of the tags, which are the UNIQue elements of |fieldgroups|
; with preserved order (we can not use UNIQ)
taggroups[0] = fieldgroups[0]
itag = 1l
for i=1l, fieldcount-1 do begin
   if max(fieldgroups[i] eq taggroups[0:itag-1]) then continue
   taggroups[itag++] = fieldgroups[i]
endfor

for i=0l, tagcount-1 do begin
   index = where(fieldgroups eq taggroups[i], count)
   tagncols[i] = count
   tagnames[i] = fieldnames[index[0]]
   tagtypes[i] = read_ascii_gettype(fieldtypes[index])
endfor


;-----------------------------
; create the output structure
;-----------------------------

; deal with columns that will be grouped as a 2D array into a single tag
dims = replicate(STRTRIM(linecount,1), tagcount)
index = where(tagncols gt 1, count)
if count gt 0 then dims[index] = STRTRIM(tagncols[index],1)+','+dims[index]

; deal with the missing value. If it is not finite, use 0 for integers
values = STRARR(tagcount)
index = where(tagtypes ne 7, count)
if count gt 0 then values[index] = ', value=missing_value'

; construct the statement
arrays = 'MAKE_ARRAY(dim=['+dims+'], type='+STRTRIM(tagtypes,1)+values+')'
str = 'result={'+STRJOIN(tagnames+':'+arrays, ',')+'}'
ok = EXECUTE(str)

;---------------------------
; fill the output structure
;---------------------------

; loop over the output structure tags
for i=0l, tagcount-1 do begin
   icol = WHERE(fieldgroups eq taggroups[i], ncol)
   ;; loop over the columns in the same group
   for j=0l, ncol-1 do begin
      row = REFORM(strresult[icol[j],*])
      if fieldtypes[i] eq 7 then begin
         if ncol eq 1 then begin
            result.(i) = row
         endif else begin
            result.(i)[j,*] = row
         endelse
      endif else begin
         index = WHERE(STREGEX(row, rnumber, /fold_case, /boolean), count)
         if count eq 0 then continue
         if ncol eq 1 then begin
            result.(i)[index] = row[index]
         endif else begin
            result.(i)[j,index] = row[index]
         endelse
      endelse
   endfor
endfor
;
if KEYWORD_SET(test) then STOP
;
return, result
;
end
