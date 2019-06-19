function RenderSpherical, lon_arr, lat_arr, r_arr, $
                          p1, p2, rmin, $
                          index, dl, lon_ind, lat_ind, r_ind, $
                          closed=closed
;Finds the voxels intersected by a straight line of sight in the spherical
;coordinate system. The voxels are specified by their boundaries with
;respect to the longitude, latitude and radial distance (the longitude,
;latitude and radial distance grids are independent on each other). The
;longitude varies from 0 to 360 degrees; the latitude is relative to equator
;and varies from -90 to 90 degrees; the radial distance is >0.
;
;Input parameters:
; lon_arr - 1D array of longitude grid points, in degrees, double;
; lat_arr - 1D array of latitude grid points, in degrees, double;
; r_arr - 1D array of radial distance grid points, in R_Sun units, double;
; p1=[r1, lon1, lat1] - 
;  3-element array specifying the first point of the line of sight, double;
; p2=[r2, lon2, lat2] - 
;  3-element array specifying the second point of the line of sight, double;
; rmin - effective photosphere radius (where the emission is absorbed),
;  in R_Sun units, double; cannot exceed min(r_arr).
;
;Keywords:
; /closed - if set, the model is assumed to cover the full rotation with
;  respect to longitude (including the region with longitude>max(lon_arr)
;  and longitude<min(lon_arr));
;  otherwise, only the region of min(lon_arr)<longitude<max(lon_arr)
;  is considered.
;
;Return value:
; Nvox - the number of intersected voxels 
;  (can be zero if no intersections found), long.
;
;Output parameters:
; index - 1D indices of the intersected voxels 
;  (1D array with Nvox elements, long);
; dl - projected depths of the intersected voxels, i.e., dl[i] is the length
;  of the fragment of the line of sight falling within the i-th voxel, 
;  in R_Sun units (1D array with Nvox elements, double);
; lon_ind, lat_ind, r_ind - fractional indices of the line-of-sight midpoints
;  within each of the intersected voxels (1D arrays with Nvox elements, double)
;  with respect to the longitude, latitude and radial distance grids;
;  the integer part of an index is the integer index of the corresponding
;  voxel and the fractional part is relative to the size of the voxel in the
;  corresponding dimension.
; The output parameters are defined if Nvox>0. All output arrays are sorted
; according to the emission propagation direction (from p1 to p2).

 common rs_parms, libname, Nmax                      
 
 if size(libname, /type) eq 0 then begin                     
  dirpath=file_dirname((routine_info('RenderSpherical', /source, /functions)).path, /mark)
  if !version.os_family eq 'Windows' then $
   libname=dirpath+'RenderSpherical'+$
           ((!version.arch eq 'x86_64') ? '64' : '32')+'.dll' $
   else libname=dirpath+'RenderSpherical.so' 
 endif
 
 if size(Nmax, /type) eq 0 then Nmax=128L
 
 closed_on=long(keyword_set(closed))

 done=0
 while ~done do begin
  index=lonarr(Nmax)
  dl=dblarr(Nmax)
  r_ind=dblarr(Nmax)
  lon_ind=dblarr(Nmax)
  lat_ind=dblarr(Nmax)
 
  Nr=n_elements(r_arr)
  Nlon=n_elements(lon_arr)
  Nlat=n_elements(lat_arr)
  
  Nvox=call_external(libname, 'RENDER', $
                     Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, $
                     p1, p2, rmin, closed_on, Nmax, $
                     index, dl, lon_ind, lat_ind, r_ind)
                    
  if Nvox le Nmax then done=1 else Nmax=(Nvox/128+1)*128
 endwhile     
 
 if Nvox gt 0 then begin
  index=index[0 : Nvox-1]
  dl=dl[0 : Nvox-1]
  r_ind=r_ind[0 : Nvox-1]
  lon_ind=lon_ind[0 : Nvox-1]
  lat_ind=lat_ind[0 : Nvox-1]
 endif   
 
 return, Nvox   
end                      