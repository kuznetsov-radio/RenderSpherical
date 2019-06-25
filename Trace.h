int TraceBline(int Nlon, int Nlat, int Nr, double *lon_arr, double *lat_arr, double *r_arr, 
	           double *Bph, double *Bth, double *Br, 
	           double lon0, double lat0, double r0, double ds, int fs_on, int Nmax,
	           double *l_lon, double *l_lat, double *l_r, 
	           int *isclosed, double *length, double *Bmax, double *Bavg, double *B1, double *B2, int returnLine);