pro read_spice, file, cl, mult=l, theta = theta



junk = headfits(file, errmsg=errmsg,/silent)

if (errmsg ne '') then begin
; ascii file

    head=' ' & fline = ' '
    openr,lun,file,/get_lun

    ; read file parameters
    readf,lun,head
    b=strpos(head,'=')
    head = strmid(head,b+1)
    args = long(strsplit(head,' ',/extract))
    nlmax = args[0] & ncor = args[1] & nside = args[2]

    ; find out number of columns
    readf,lun,fline
    fields = strsplit(fline,' ',/extract)
    nf = n_elements(fields)
    free_lun,lun

    case (nf - ncor) of
        1: readcol, file, l,         clt, cle, clb, clx, cly, clz ; ang. power spectrum
        2: readcol, file, theta, ct, clt, cle, clb, clx, cly, clz ; ang. correlation
        else: begin
            message,'wrong number of fields in file '+file
        end
    endcase

    case ncor of
        1: cl = clt ; T only
        4: cl = [[clt],[cle],[clb],[clx]]
        6: cl = [[clt],[cle],[clb],[clx],[cly],[clz]]
        else: begin
            message,'unable to read file '+file
        end
    endcase

endif else begin
; FITS file
    st=mrdfits(file,1,/silent)

    junk = where(tag_names(st) eq 'ANGLE', nj)
    if (nj gt 0) then begin
                                ; ang. correlation
        offset = 1
        na = n_elements(st.(0))
        theta = st.angle

    endif else begin
                                ; ang. power spectrum
        offset = 0
        na = n_elements(st.(0))
        l = dindgen(na)

    endelse

    nf = n_tags(st)-offset
    cl = dblarr(na, nf)
    for i=0, nf-1 do begin
        cl[0,i] = st.(i+offset)
    endfor

endelse

return
end

