#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif
#include "Render.h"
#include "Trace.h"

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
extern "C" __declspec(dllexport) int TRACE(int argc, void **argv)
#else
extern "C" int TRACE(int argc, void **argv)
#endif
{
 int Nlon=*((int*)argv[0]);
 int Nlat=*((int*)argv[1]);
 int Nr=*((int*)argv[2]);
 double *lon_arr=(double*)argv[3];
 double *lat_arr=(double*)argv[4];
 double *r_arr=(double*)argv[5];
 double *Bph=(double*)argv[6];
 double *Bth=(double*)argv[7];
 double *Br=(double*)argv[8];
 double lon0=*((double*)argv[9]);
 double lat0=*((double*)argv[10]);
 double r0=*((double*)argv[11]);
 double ds=*((double*)argv[12]);
 int fs_on=*((int*)argv[13]);
 int Nmax=*((int*)argv[14]);
 double *l_lon=(double*)argv[15];
 double *l_lat=(double*)argv[16];
 double *l_r=(double*)argv[17];
 int *isclosed=(int*)argv[18];
 double *length=(double*)argv[19];
 double *Bmax=(double*)argv[20];
 double *Bavg=(double*)argv[21];
 double *B1=(double*)argv[22];
 double *B2=(double*)argv[23];

 int N=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0, lat0, r0, ds, fs_on, Nmax,
	              l_lon, l_lat, l_r, isclosed, length, Bmax, Bavg, B1, B2, 1);

 return N;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int TRACE_SHORT(int argc, void **argv)
#else
extern "C" int TRACE_SHORT(int argc, void **argv)
#endif
{
 int Nlon=*((int*)argv[0]);
 int Nlat=*((int*)argv[1]);
 int Nr=*((int*)argv[2]);
 double *lon_arr=(double*)argv[3];
 double *lat_arr=(double*)argv[4];
 double *r_arr=(double*)argv[5];
 double *Bph=(double*)argv[6];
 double *Bth=(double*)argv[7];
 double *Br=(double*)argv[8];
 double lon0=*((double*)argv[9]);
 double lat0=*((double*)argv[10]);
 double r0=*((double*)argv[11]);
 double ds=*((double*)argv[12]);
 int fs_on=*((int*)argv[13]);
 int *isclosed=(int*)argv[14];
 double *length=(double*)argv[15];
 double *Bmax=(double*)argv[16];
 double *Bavg=(double*)argv[17];
 double *B1=(double*)argv[18];
 double *B2=(double*)argv[19];

 int N=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0, lat0, r0, ds, fs_on, 0,
	              0, 0, 0, isclosed, length, Bmax, Bavg, B1, B2, 0);

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

#ifndef LINUX
extern "C" __declspec(dllexport) int TRACE_MULTI(int argc, void **argv)
#else
extern "C" int TRACE_MULTI(int argc, void **argv)
#endif
{
 int Nlines=*((int*)argv[0]);
 int Nlon=*((int*)argv[1]);
 int Nlat=*((int*)argv[2]);
 int Nr=*((int*)argv[3]);
 double *lon_arr=(double*)argv[4];
 double *lat_arr=(double*)argv[5];
 double *r_arr=(double*)argv[6];
 double *Bph=(double*)argv[7];
 double *Bth=(double*)argv[8];
 double *Br=(double*)argv[9];
 double *lon0_m=(double*)argv[10];
 double *lat0_m=(double*)argv[11];
 double *r0_m=(double*)argv[12];
 double ds=*((double*)argv[13]);
 int fs_on=*((int*)argv[14]);
 int Nmax=*((int*)argv[15]);
 int *N_m=(int*)argv[16];
 double *l_lon_m=(double*)argv[17];
 double *l_lat_m=(double*)argv[18];
 double *l_r_m=(double*)argv[19];
 int *isclosed_m=(int*)argv[20];
 double *length_m=(double*)argv[21];
 double *Bmax_m=(double*)argv[22];
 double *Bavg_m=(double*)argv[23];
 double *B1_m=(double*)argv[24];
 double *B2_m=(double*)argv[25];

 #ifndef LINUX
 concurrency :: parallel_for(0, Nlines, [&](int i)
 {
  N_m[i]=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0_m[i], lat0_m[i], r0_m[i], ds, fs_on, Nmax,
	                l_lon_m+i*Nmax, l_lat_m+i*Nmax, l_r_m+i*Nmax, isclosed_m+i, length_m+i, Bmax_m+i, Bavg_m+i, B1_m+i, B2_m+i, 1);
 });
 #else
 #pragma omp parallel for
 for (int i=0; i<Nlines; i++)
 {
  N_m[i]=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0_m[i], lat0_m[i], r0_m[i], ds, fs_on, Nmax,
	                l_lon_m+i*Nmax, l_lat_m+i*Nmax, l_r_m+i*Nmax, isclosed_m+i, length_m+i, Bmax_m+i, Bavg_m+i, B1_m+i, B2_m+i, 1);
 }
 #endif

 int NptsMax=0;
 for (int i=0; i<Nlines; i++) if (N_m[i]>NptsMax) NptsMax=N_m[i];

 return NptsMax;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int TRACE_SHORT_MULTI(int argc, void **argv)
#else
extern "C" int TRACE_SHORT_MULTI(int argc, void **argv)
#endif
{
 int Nlines=*((int*)argv[0]);
 int Nlon=*((int*)argv[1]);
 int Nlat=*((int*)argv[2]);
 int Nr=*((int*)argv[3]);
 double *lon_arr=(double*)argv[4];
 double *lat_arr=(double*)argv[5];
 double *r_arr=(double*)argv[6];
 double *Bph=(double*)argv[7];
 double *Bth=(double*)argv[8];
 double *Br=(double*)argv[9];
 double *lon0_m=(double*)argv[10];
 double *lat0_m=(double*)argv[11];
 double *r0_m=(double*)argv[12];
 double ds=*((double*)argv[13]);
 int fs_on=*((int*)argv[14]);
 int *N_m=(int*)argv[15];
 int *isclosed_m=(int*)argv[16];
 double *length_m=(double*)argv[17];
 double *Bmax_m=(double*)argv[18];
 double *Bavg_m=(double*)argv[19];
 double *B1_m=(double*)argv[20];
 double *B2_m=(double*)argv[21];

 #ifndef LINUX
 concurrency :: parallel_for(0, Nlines, [&](int i)
 {
  N_m[i]=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0_m[i], lat0_m[i], r0_m[i], ds, fs_on, 0,
	                0, 0, 0, isclosed_m+i, length_m+i, Bmax_m+i, Bavg_m+i, B1_m+i, B2_m+i, 0);
 });
 #else
 #pragma omp parallel for
 for (int i=0; i<Nlines; i++)
 {
  N_m[i]=TraceBline(Nlon, Nlat, Nr, lon_arr, lat_arr, r_arr, Bph, Bth, Br, lon0_m[i], lat0_m[i], r0_m[i], ds, fs_on, 0,
	                0, 0, 0, isclosed_m+i, length_m+i, Bmax_m+i, Bavg_m+i, B1_m+i, B2_m+i, 0);
 }
 #endif

 int NptsMax=0;
 for (int i=0; i<Nlines; i++) if (N_m[i]>NptsMax) NptsMax=N_m[i];

 return NptsMax;
}