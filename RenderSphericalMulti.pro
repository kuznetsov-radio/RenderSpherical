function RenderSphericalMulti, NLOS, lon_arr, lat_arr, r_arr, $
                               p1_m, p2_m, rmin, $
                               index_m, dl_m, lon_ind_m, lat_ind_m, r_ind_m, $
                               closed=closed
;Finds the voxels intersected by a number of straight lines of sight in the 
;spherical coordinate system. The voxels are specified by their boundaries with
;respect to the longitude, latitude and radial distance (the longitude,
;latitude and radial distance grids are independent on each other). The
;longitude varies from 0 to 360 degrees; the latitude is relative to equator
;and varies from -90 to 90 degrees; the radial distance is >0.
;
;Input parameters:
; NLOS - number of the lines of sight
; lon_arr - 1D array of longitude grid points, in degrees, double;
; lat_arr - 1D array of latitude grid points, in degrees, double;
; r_arr - 1D array of radial distance grid points, in R_Sun units, double;
; p1_m, p2_m - 2D arrays (3, NLOS) specifying the lines of sight by their
;  endpoints, double: j-th line of sight is assumed to be from 
;  [r1, lon1, lat1]=p1_m[*, j] to [r2, lon2, lat2]=p2_m[*, j];
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
; Nvox_m - 1D array with NLOS elements specifying the numbers of intersected 
; voxels for each line of sight, long (some elements can be zero if no 
; intersections found).
;
;Output parameters:
; index_m - 2D array (Nmax, NLOS) specifying 1D indices of the intersected 
;  voxels, long;
; dl_m - 2D array (Nmax, NLOS) specifying projected depths of the intersected 
;  voxels (i.e., lengths of the fragments of the respective lines of sight
;  falling within the intersected voxels), in R_Sun units, double;
; lon_ind_m, lat_ind_m, r_ind_m - 2D arrays (Nmax, NLOS) specifying fractional 
;  indices (double) of the line-of-sight midpoints within the intersected 
;  voxels with respect to the longitude, latitude and radial distance grids;
;  the integer part of an index is the integer index of the corresponding
;  voxel and the fractional part is relative to the size of the voxel in the
;  corresponding dimension.
; In all output parameters, Nmax is chosen as an array size sufficient to 
; accommodate all results, i.e., Nmax equals or exceeds max(Nvox_m).
; The value index_m[i, j] represents the 1D index of the i-th intersected
; voxel along j-th line of sight (defined only for i<Nvox_m[j]); the same
; convention is used for all other output arrays.
; Each j-th column of the output arrays is sorted according to the emission
; propagation direction (from p1_m[*, j] to p2_m[*, j]).

 common rs_parmsM, libname, Nmax                      
 
 if size(libname, /type) eq 0 then begin                     
  dirpath=file_dirname((routine_info('RenderSphericalMulti', /source, /functions)).path, /mark)
  if !version.os_family eq 'Windows' then $
   libname=dirpath+'RenderSpherical'+$
           ((!version.arch eq 'x86_64') ? '64' : '32')+'.dll' $
   else libname=dirpath+'RenderSpherical.so' 
 endif
 
 if size(Nmax, /type) eq 0 then Nmax=128L
 
 closed_on=long(keyword_set(closed))

 Nvox_m=lonarr(NLOS)

 done=0
 while ~done do begin
  index_m=lonarr(Nmax, NLOS)
  dl_m=dblarr(Nmax, NLOS)
  r_ind_m=dblarr(Nmax, NLOS)
  lon_ind_m=dblarr(Nmax, NLOS)
  lat_ind_m=dblarr(Nmax, NLOS)
 
  Nr=n_elements(r_arr)
  Nlon=n_elements(lon_arr)
  Nlat=n_elements(lat_arr)
  
  NvoxMax=call_external(libname, 'RENDER_MULTI', $
                        NLOS, Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, $
                        p1_m, p2_m, rmin, closed_on, Nmax, $
                        Nvox_m, index_m, dl_m, lon_ind_m, lat_ind_m, r_ind_m)
                    
  if NvoxMax le Nmax then done=1 else Nmax=(NvoxMax/128+1)*128
 endwhile     
 
 return, Nvox_m
end                      