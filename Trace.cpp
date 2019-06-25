#define _USE_MATH_DEFINES
#include <malloc.h>
#include <math.h>
#include <float.h>

inline double sqr(double x)
{
 return x*x;
}

double index1D(int N, double *x_arr, double x)
{
 int i1=0;
 int i2=N-1;
 while (i2-i1>1)
 {
  int k=(i1+i2)>>1;
  if (x_arr[k]>x) i2=k;
  else i1=k;
 }

 double dx=(x-x_arr[i1])/(x_arr[i2]-x_arr[i1]);

 return dx+i1;
}

#define D3(arr, s1, s2, i1, i2, i3) arr[(i1)+((i2)+(i3)*(s2))*(s1)]

double TriLinear(int Nx, int Ny, int Nz, double *a, double x_ind, double y_ind, double z_ind)
{
 int i1=int(x_ind);
 int i2=i1+1;
 int j1=int(y_ind);
 int j2=j1+1;
 int k1=int(z_ind);
 int k2=k1+1;

 double dx=x_ind-i1;
 double dy=y_ind-j1;
 double dz=z_ind-k1;

 double res=D3(a, Nx, Ny, i1, j1, k1)*(1.0-dx)*(1.0-dy)*(1.0-dz)+
            D3(a, Nx, Ny, i2, j1, k1)*dx*(1.0-dy)*(1.0-dz)+
            D3(a, Nx, Ny, i1, j2, k1)*(1.0-dx)*dy*(1.0-dz)+
            D3(a, Nx, Ny, i1, j1, k2)*(1.0-dx)*(1.0-dy)*dz+
            D3(a, Nx, Ny, i2, j1, k2)*dx*(1.0-dy)*dz+
            D3(a, Nx, Ny, i1, j2, k2)*(1.0-dx)*dy*dz+
            D3(a, Nx, Ny, i2, j2, k1)*dx*dy*(1.0-dz)+
            D3(a, Nx, Ny, i2, j2, k2)*dx*dy*dz;

 return res;
}

#define NMAX 1000000

class PointStack
{
 double *lon_arr, *lat_arr, *r_arr;
 int N, Nmax, keepAll;
 double Bmax, Btotal;
 public:
 PointStack(int returnLine);
 ~PointStack();
 int push(double lon, double lat, double r, double Bprev);
 void Update(double Bprev);
 void GetInfo(int *NH, double *BmaxH, double *BtotalH);
 void ExportData(double *l_lon, double *l_lat, double *l_r, int i0, int reverse);
};

PointStack :: PointStack(int returnLine)
{
 lon_arr=lat_arr=r_arr=0;
 N=Nmax=0;
 Bmax=Btotal=0;
 keepAll=returnLine;
}

PointStack :: ~PointStack()
{
 if (lon_arr) free(lon_arr);
 if (lat_arr) free(lat_arr);
 if (r_arr) free(r_arr);
}

int PointStack :: push(double lon, double lat, double r, double Bprev)
{
 if (keepAll)
 {
  if (N>=Nmax)
  {
   Nmax+=1024;
   lon_arr=lon_arr ? (double*)realloc(lon_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
   lat_arr=lat_arr ? (double*)realloc(lat_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
   r_arr=  r_arr   ? (double*)realloc(r_arr, sizeof(double)*Nmax)   : (double*)malloc(sizeof(double)*Nmax);
  }

  lon_arr[N]=lon;
  lat_arr[N]=lat;
  r_arr[N]=r;
 }

 if (Bprev>Bmax) Bmax=Bprev;
 Btotal+=(N>0) ? Bprev : Bprev/2;

 N++;

 return N<NMAX;
}

void PointStack :: Update(double Bprev)
{
 if (Bprev>Bmax) Bmax=Bprev;
 Btotal+=(N>0) ? Bprev : Bprev/2;
}

void PointStack :: GetInfo(int *NH, double *BmaxH, double *BtotalH)
{
 *NH=N;
 *BmaxH=Bmax;
 *BtotalH=Btotal;
}

void PointStack :: ExportData(double *l_lon, double *l_lat, double *l_r, int i0, int reverse)
{
 for (int i=0; i<N; i++)
 {
  int j=reverse ? i0+N-1-i : i0+i;
  l_lon[j]=lon_arr[i];
  l_lat[j]=lat_arr[i];
  l_r[j]=r_arr[i];
 }
}

int TraceBline(int Nlon, int Nlat, int Nr, double *lon_arrX, double *lat_arrX, double *r_arr, 
	           double *Bph, double *Bth, double *Br, 
	           double lon0X, double lat0X, double r0, double dsX, int fs_on, int Nmax,
	           double *l_lon, double *l_lat, double *l_r, 
	           int *isclosed, double *length, double *Bmax, double *Bavg, double *B1, double *B2, int returnLine)
{
 double lonC=(lon_arrX[0]+lon_arrX[Nlon-1])/2;

 double *lon_arr=(double*)malloc(Nlon*sizeof(double));
 for (int i=0; i<Nlon; i++) lon_arr[i]=(lon_arrX[i]-lonC)*M_PI/180;

 double *lat_arr=(double*)malloc(Nlat*sizeof(double));
 for (int j=0; j<Nlat; j++) lat_arr[j]=lat_arrX[j]*M_PI/180;

 double lon0=(lon0X-lonC)*M_PI/180;
 double lat0=lat0X*M_PI/180;

 int N=0;

 if (lon0>lon_arr[0] && lon0<lon_arr[Nlon-1] &&
	 lat0>lat_arr[0] && lat0<lat_arr[Nlat-1] &&
	 r0>r_arr[0] && r0<r_arr[Nr-1])
 {
  PointStack* ps[2];
  for (int d=0; d<2; d++) ps[d]=new PointStack(returnLine);

  double Bedge[2];

  int Nlb=0;

  for (int d=0; d<2; d++)
  {
   double ds=dsX*(d*2-1);

   double p_sph[3], dr[3];

   p_sph[0]=lon0;
   p_sph[1]=lat0;
   p_sph[2]=r0;

   int inbox;
   do
   {
    double lon_ind=index1D(Nlon, lon_arr, p_sph[0]);
	double lat_ind=index1D(Nlat, lat_arr, p_sph[1]);
	double r_ind=index1D(Nr, r_arr, p_sph[2]);
	
	double Bph_loc=TriLinear(Nlon, Nlat, Nr, Bph, lon_ind, lat_ind, r_ind);
	double Bth_loc=TriLinear(Nlon, Nlat, Nr, Bth, lon_ind, lat_ind, r_ind);
	double Br_loc=TriLinear(Nlon, Nlat, Nr, Br, lon_ind, lat_ind, r_ind);
	double B=sqrt(Bph_loc*Bph_loc+Bth_loc*Bth_loc+Br_loc*Br_loc);
	
	if (B>0)
	{
	 double ct=cos(p_sph[1]);
	 dr[0]=Bph_loc/p_sph[2]/ct;
	 dr[1]=Bth_loc/p_sph[2];
	 dr[2]=Br_loc;
	 double drAbs=sqrt(sqr(dr[2])+sqr(p_sph[2])*(sqr(dr[1])+sqr(ct*dr[0]))); //dr^2+r^2*[dtheta^2+(cos(theta)*dphi)^2]

	 for (int k=0; k<3; k++) p_sph[k]+=dr[k]/drAbs*ds;
	}

	inbox=p_sph[0]>lon_arr[0] && p_sph[0]<lon_arr[Nlon-1] &&
		  p_sph[1]>lat_arr[0] && p_sph[1]<lat_arr[Nlat-1] &&
		  p_sph[2]>r_arr[0] && p_sph[2]<r_arr[Nr-1] && B>0;
	if (inbox) inbox=ps[d]->push(p_sph[0], p_sph[1], p_sph[2], B);
	else 
	{
	 ps[d]->Update(B);
	 Bedge[d]=B;
	 if (p_sph[2]<=r_arr[0]) Nlb++;
	}
   }
   while (inbox);
  }

  int NH[2];
  double BmaxH[2], BtotalH[2];
  for (int d=0; d<2; d++) ps[d]->GetInfo(NH+d, BmaxH+d, BtotalH+d);

  N=NH[0]+NH[1]+1;
  *isclosed=(Nlb==2);
  *length=dsX*(N-1);
  *Bmax=(BmaxH[0]>BmaxH[1]) ? BmaxH[0] : BmaxH[1];
  *Bavg=(BtotalH[0]+BtotalH[1])/N;
  *B1=Bedge[0];
  *B2=Bedge[1];

  if (returnLine) if (N<=Nmax)
  {
   ps[0]->ExportData(l_lon, l_lat, l_r, 0, 1);
   l_lon[NH[0]]=lon0;
   l_lat[NH[0]]=lat0;
   l_r[NH[0]]=r0;
   ps[1]->ExportData(l_lon, l_lat, l_r, NH[0]+1, 0);

   for (int i=0; i<N; i++)
   {
	l_lon[i]*=(180.0/M_PI);
	l_lon[i]+=lonC;
	l_lat[i]*=(180.0/M_PI);
   }
  }

  for (int d=0; d<2; d++) delete ps[d];
 }

 free(lon_arr);
 free(lat_arr);

 return N;
}