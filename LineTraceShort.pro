function LineTraceShort, lon_arr, lat_arr, r_arr, Bph, Bth, Br, $
                         lon0, lat0, r0, ds, $
                         lineclosed, L, Bmax, Bavg, B1, B2, fs=fs
;Computes the parameters of the magnetic field line in the spherical 
;coordinate system. The difference from the function LineTrace is that
;LineTraceShort does not return the line coordinates themselves - only
;some general information about the line.
;The magnetic field is specified at the nodes of 3D grid (vs. the longitude,
;latitude, and radial distance). The longitude, latitude and radial distance
;grids are independent on each other. The longitude varies from 0 to 360 
;degrees; the latitude is relative to equator and varies from -90 to 90
;degrees; the radial distance is >0.
;
;Input parameters:
; lon_arr - 1D array of longitude grid points, in degrees, double;
; lat_arr - 1D array of latitude grid points, in degrees, double;
; r_arr - 1D array of radial distance grid points, in R_Sun units, double;
; Bph, Bth, Br - 3D arrays (double) specifying respectively the longitudinal
;  (Bph), latitudinal (Bth) and radial (Br) components of the magnetic field 
;  vector at the grid nodes;
; lon0, lat0, r0 - coordinates of the initial point of the field line
;  (lon0 and lat0 in degrees, r0 in R_Sun units), double;
; ds - integration step (length of one element of the field line), in R_Sun
;  units, double.
;
;Return value:
; N - the number of points in the field line, long;
;  can be zero, if the initial point is located outside the data cube.
;
;Output parameters:
; lineclosed - 1 if the feld line is closed and 0 otherwise, long;
; L - length of the field line, in R_Sun units, double;
; Bmax - maximum magnetic field strength along the field line, double;
; Bavg - average magnetic field strength along the field line, double;
; B1, B2 - magnetic field strengths at the field line endpoints, double.
;
;Note:
; The line-tracing algorithm stops if the number of points in any half of the
; field line exceeds 1000000.

 common lts_parms, libname
                                           
 if size(libname, /type) eq 0 then begin                     
  dirpath=file_dirname((routine_info('LineTraceShort', /source, /functions)).path, /mark)
  if !version.os_family eq 'Windows' then $
   libname=dirpath+'RenderSpherical'+$
           ((!version.arch eq 'x86_64') ? '64' : '32')+'.dll' $
   else libname=dirpath+'RenderSpherical.so' 
 endif
 
 fs_on=long(keyword_set(fs))
 
 lineclosed=0L
 l=0d0
 Bmax=0d0
 Bavg=0d0
 B1=0d0
 B2=0d0
 
 Nlon=n_elements(lon_arr)
 Nlat=n_elements(lat_arr)
 Nr=n_elements(r_arr)
 
 N=call_external(libname, 'TRACE_SHORT', $
                 Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, $
                 Bph, Bth, Br, lon0, lat0, r0, ds, fs_on, $
                 lineclosed, l, Bmax, Bavg, B1, B2)
                  
 return, N             
end         