function LineTraceShortMulti, Nlines, lon_arr, lat_arr, r_arr, Bph, Bth, Br, $
                              lon0_m, lat0_m, r0_m, ds, $
                              lineclosed_m, L_m, Bmax_m, Bavg_m, B1_m, B2_m, fs=fs
;Computes the parameters of a number of the magnetic field lines in the 
;spherical coordinate system. The difference from the function LineTraceMulti
;is that LineTraceShortMulti does not return the line coordinates themselves -
;only some general information about the lines.
;The magnetic field is specified at the nodes of 3D grid (vs. the longitude,
;latitude, and radial distance). The longitude, latitude and radial distance
;grids are independent on each other. The longitude varies from 0 to 360 
;degrees; the latitude is relative to equator and varies from -90 to 90
;degrees; the radial distance is >0.
;
;Input parameters:
; Nlines - number of the field lines, long;
; lon_arr - 1D array of longitude grid points, in degrees, double;
; lat_arr - 1D array of latitude grid points, in degrees, double;
; r_arr - 1D array of radial distance grid points, in R_Sun units, double;
; Bph, Bth, Br - 3D arrays (double) specifying respectively the longitudinal
;  (Bph), latitudinal (Bth) and radial (Br) components of the magnetic field 
;  vector at the grid nodes;
; lon0_m, lat0_m, r0_m - 1D arrays (with Nlines elements, double) specifying
;  the coordinates of the initial points of the field lines: for the j-th 
;  line, the coordinates  are assumed to be (lon0_m[j], lat0_m[j], r0_m[j]);
;  lon0_m and lat0_m in degrees, r0_m in R_Sun units;
; ds - integration step (length of one element of the field line), in R_Sun
;  units, double.
;
;Return value:
; N_m - 1D array with Nlines elements specifying the number of points in 
;  each field line, long (some elements can be zero if the initial point is 
;  located outside the data cube).
;
;Output parameters:
; lineclosed_m - 1D array with Nlines elements (long) specifying whether the
;  respective field line is closed (lineclosed_m[j]=1) or open
;  (lineclosed_m[j]=0);
; L_m - 1D array with Nlines elements (double) specifying the lengths of the 
;  field lines, in R_Sun units;
; Bmax_m - 1D array with Nlines elements (double) specifying the maximum 
;  magnetic field strengths along the field lines;
; Bavg_m - 1D array with Nlines elements (double) specifying the average 
;  magnetic field strengths along the field lines;
; B1_m, B2_m - 1D arrays with Nlines elements (double) specifying the magnetic 
;  field strengths at the endpoints of the respective field lines.
;
;Note:
; The line-tracing algorithm stops if the number of points in any half of a
; field line exceeds 1000000.

 common ltsm_parms, libname
                                           
 if size(libname, /type) eq 0 then begin                     
  dirpath=file_dirname((routine_info('LineTraceShortMulti', /source, /functions)).path, /mark)
  if !version.os_family eq 'Windows' then $
   libname=dirpath+'RenderSpherical'+$
           ((!version.arch eq 'x86_64') ? '64' : '32')+'.dll' $
   else libname=dirpath+'RenderSpherical.so' 
 endif

 fs_on=long(keyword_set(fs))

 lineclosed_m=lonarr(Nlines)
 l_m=dblarr(Nlines)
 Bmax_m=dblarr(Nlines)
 Bavg_m=dblarr(Nlines)
 B1_m=dblarr(Nlines)
 B2_m=dblarr(Nlines)
 N_m=lonarr(Nlines)
 
 Nlon=n_elements(lon_arr)
 Nlat=n_elements(lat_arr)
 Nr=n_elements(r_arr)
 
 NptsMax=call_external(libname, 'TRACE_SHORT_MULTI', $
                       Nlines, Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, $
                       Bph, Bth, Br, $
                       lon0_m, lat0_m, r0_m, ds, fs_on, $
                       N_m, lineclosed_m, l_m, Bmax_m, Bavg_m, B1_m, B2_m)
 
 return, N_m             
end      