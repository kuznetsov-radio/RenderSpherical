#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif
#include "Render.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER(int argc, void **argv)
#else
extern "C" int RENDER(int argc, void **argv)
#endif
{
 int Nlon=*((int*)argv[0]);
 int Nlat=*((int*)argv[1]);
 int Nr=*((int*)argv[2]);
 double *lon_arr=(double*)argv[3];
 double *lat_arr=(double*)argv[4];
 double *r_arr=(double*)argv[5];
 double *p1=(double*)argv[6];
 double *p2=(double*)argv[7];
 double rmin=*((double*)argv[8]);
 int closed=*((int*)argv[9]);
 int Nmax=*((int*)argv[10]);
 int *index=(int*)argv[11];
 double *dl=(double*)argv[12];
 double *lon_ind=(double*)argv[13];
 double *lat_ind=(double*)argv[14];
 double *r_ind=(double*)argv[15];

 int N=RenderSpherical(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr,
	                   p1, p2, rmin, closed, Nmax, 
	                   index, dl, lon_ind, lat_ind, r_ind);

 return N;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER_MULTI(int argc, void **argv)
#else
extern "C" int RENDER_MULTI(int argc, void **argv)
#endif
{
 int NLOS=*((int*)argv[0]);
 int Nlon=*((int*)argv[1]);
 int Nlat=*((int*)argv[2]);
 int Nr=*((int*)argv[3]);
 double *lon_arr=(double*)argv[4];
 double *lat_arr=(double*)argv[5];
 double *r_arr=(double*)argv[6];
 double *p1_m=(double*)argv[7];
 double *p2_m=(double*)argv[8];
 double rmin=*((double*)argv[9]);
 int closed=*((int*)argv[10]);
 int Nmax=*((int*)argv[11]);
 int *Nvox_m=(int*)argv[12];
 int *index_m=(int*)argv[13];
 double *dl_m=(double*)argv[14];
 double *lon_ind_m=(double*)argv[15];
 double *lat_ind_m=(double*)argv[16];
 double *r_ind_m=(double*)argv[17];

 #ifndef LINUX
 concurrency::parallel_for(0, NLOS, [&](int i)
 {
  Nvox_m[i]=RenderSpherical(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr,
	                        p1_m+i*3, p2_m+i*3, rmin, closed, Nmax, 
	                        index_m+i*Nmax, dl_m+i*Nmax, lon_ind_m+i*Nmax, lat_ind_m+i*Nmax, r_ind_m+i*Nmax);
 });
 #else
 #pragma omp parallel for
 for (int i=0; i<NLOS; i++)
 {
  Nvox_m[i]=RenderSpherical(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr,
	                        p1_m+i*3, p2_m+i*3, rmin, closed, Nmax, 
	                        index_m+i*Nmax, dl_m+i*Nmax, lon_ind_m+i*Nmax, lat_ind_m+i*Nmax, r_ind_m+i*Nmax);
 }
 #endif
 
 int N=0;
 for (int i=0; i<NLOS; i++) if (Nvox_m[i]>N) N=Nvox_m[i];

 return N;
}