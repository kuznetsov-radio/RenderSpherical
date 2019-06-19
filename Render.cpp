#define _USE_MATH_DEFINES
#include <malloc.h>
#include <math.h>
#include <float.h>
#include <memory.h>

inline double sqr(double x)
{
 return x*x;
}

#ifndef LINUX
#define finite _finite
#else
#define finite isfinite
#endif

void sort(double *a, int N, int *idx)
{
 if (idx) for (int i=0; i<N; i++) idx[i]=i;

 for (int i=0; i<(N-1); i++) for (int j=i+1; j<N; j++) if (a[i]>a[j])
 {
  double tmp=a[i];
  a[i]=a[j];
  a[j]=tmp;

  if (idx)
  {
   int itmp=idx[i];
   idx[i]=idx[j];
   idx[j]=itmp;
  }
 }
}

class SphGrid
{
 int Nr, Nlon, Nlat;
 double *r_arr, *lon_arr, *lat_arr;
 int *r_cN, *lon_cN, *lat_cN;
 double *r_cL, *lon_cL, *lat_cL;
 double x1, y1, z1, ax, ay, az;
 int FindCrossingR(int i, double *t_c, double *r_c, double *lon_c, double *lat_c);
 int FindCrossingLon(int j, double *t_c, double *r_c, double *lon_c, double *lat_c);
 int FindCrossingLat(int k, double *t_c, double *r_c, double *lon_c, double *lat_c);
 public:
 double L;
 SphGrid(int _Nr, int _Nlon, int _Nlat, double *_r_arr, double *_lon_arr, double *_lat_arr, 
	     double *p1, double *p2);
 ~SphGrid();
 int FindCrossings(int i1, int i2, int j1, int j2, int k1, int k2, double *t_arr);
 int FindCrossingRmin(double rmin, double *t_c);
 void FindIndices(double tc, int i1, int i2, int j1, int j2, int k1, int k2,
	              double *r_idx, double *lon_idx, double *lat_idx);
};

SphGrid :: SphGrid(int _Nr, int _Nlon, int _Nlat, double *_r_arr, double *_lon_arr, double *_lat_arr, 
	               double *p1, double *p2)
{
 Nr=_Nr;
 Nlon=_Nlon;
 Nlat=_Nlat;

 r_arr=(double*)malloc(sizeof(double)*Nr);
 for (int i=0; i<Nr; i++) r_arr[i]=_r_arr[i];

 lon_arr=(double*)malloc(sizeof(double)*Nlon);
 double phi_max=-1e100;
 double phi_min=1e100;
 for (int j=0; j<Nlon; j++)
 {
  lon_arr[j]=_lon_arr[j]*M_PI/180;
  if (lon_arr[j]>phi_max) phi_max=lon_arr[j];
  if (lon_arr[j]<phi_min) phi_min=lon_arr[j];
 }
 double phi0=(phi_min+phi_max)/2;
 for (int j=0; j<Nlon; j++) lon_arr[j]-=phi0;

 lat_arr=(double*)malloc(sizeof(double)*Nlat);
 for (int k=0; k<Nlat; k++) lat_arr[k]=_lat_arr[k]*M_PI/180;

 r_cN=(int*)malloc(sizeof(int)*Nr);
 for (int i=0; i<Nr; i++) r_cN[i]=-1;

 lon_cN=(int*)malloc(sizeof(int)*Nlon);
 for (int j=0; j<Nlon; j++) lon_cN[j]=-1;

 lat_cN=(int*)malloc(sizeof(int)*Nlat);
 for (int k=0; k<Nlat; k++) lat_cN[k]=-1;

 r_cL=(double*)malloc(sizeof(double)*Nr*8);
 lon_cL=(double*)malloc(sizeof(double)*Nlon*4);
 lat_cL=(double*)malloc(sizeof(double)*Nlat*8);

 double phi1=p1[1]*M_PI/180-phi0;
 double theta1=p1[2]*M_PI/180;
 double phi2=p2[1]*M_PI/180-phi0;
 double theta2=p2[2]*M_PI/180;
 x1=p1[0]*sin(phi1)*cos(theta1);
 y1=p1[0]*sin(theta1);
 z1=p1[0]*cos(phi1)*cos(theta1);
 ax=p2[0]*sin(phi2)*cos(theta2)-x1;
 ay=p2[0]*sin(theta2)-y1;
 az=p2[0]*cos(phi2)*cos(theta2)-z1;

 L=sqrt(sqr(ax)+sqr(ay)+sqr(az));
}

SphGrid :: ~SphGrid()
{
 free(r_arr);
 free(lon_arr);
 free(lat_arr);

 free(r_cN);
 free(lon_cN);
 free(lat_cN);

 free(r_cL);
 free(lon_cL);
 free(lat_cL);
}

int SphGrid :: FindCrossingR(int i, double *t_c, double *r_c, double *lon_c, double *lat_c)
{
 int N;

 if (r_cN[i]<0)
 {
  N=0;

  double A=sqr(ax)+sqr(ay)+sqr(az);
  double B=2.0*(ax*x1+ay*y1+az*z1);
  double C=sqr(x1)+sqr(y1)+sqr(z1)-sqr(r_arr[i]);

  double D=sqr(B)-4.0*A*C;
  if (D>0)
  {
   D=sqrt(D);

   for (int s=-1; s<=1; s+=2)
   {
	double t=(-B+D*s)/(2.0*A);
	if (finite(t)) if (t>=0 && t<=1)
	{
	 double x=x1+ax*t;
	 double y=y1+ay*t;
	 double z=z1+az*t;
	 double phi=atan2(x, z);
	 double theta=asin(y/r_arr[i]);

	 double *q=r_cL+i*8+N*4;
	 q[0]=t;
	 q[1]=r_arr[i];
	 q[2]=phi;
	 q[3]=theta;

	 N++;
	}
   }
  }

  r_cN[i]=N;
 }

 N=r_cN[i];
 if (N>0) for (int m=0; m<N; m++)
 {
  double *q=r_cL+i*8+m*4;
  t_c[m]=q[0];
  r_c[m]=q[1];
  lon_c[m]=q[2];
  lat_c[m]=q[3];
 }

 return N;
}

int SphGrid :: FindCrossingRmin(double rmin, double *t_c)
{
 int N=0;

 double A=sqr(ax)+sqr(ay)+sqr(az);
 double B=2.0*(ax*x1+ay*y1+az*z1);
 double C=sqr(x1)+sqr(y1)+sqr(z1)-sqr(rmin);

 double D=sqr(B)-4.0*A*C;
 if (D>0)
 {
  D=sqrt(D);

  for (int s=-1; s<=1; s+=2)
  {
   double t=(-B+D*s)/(2.0*A);

   if (finite(t)) if (t>=0 && t<=1)
   {
	t_c[N]=t;
	N++;
   }
  }
 }

 return N;
}

int SphGrid :: FindCrossingLon(int j, double *t_c, double *r_c, double *lon_c, double *lat_c)
{
 if (lon_cN[j]<0)
 {
  lon_cN[j]=0;

  double cp=cos(lon_arr[j]);
  double sp=sin(lon_arr[j]);
  double t=-(x1*cp-z1*sp)/(ax*cp-az*sp);
  if (finite(t)) if (t>=0 && t<=1)
  {
   double x=x1+ax*t;
   double z=z1+az*t;

   if ((x*sp)>=0 && (z*cp)>=0)
   {
	double y=y1+ay*t;
	double r=sqrt(sqr(x)+sqr(y)+sqr(z));
	double theta=asin(y/r);

	lon_cN[j]=1;
	double *q=lon_cL+j*4;
	q[0]=t;
	q[1]=r;
	q[2]=lon_arr[j];
	q[3]=theta;
   }
  }
 }

 int N=lon_cN[j];
 if (N>0)
 {
  double *q=lon_cL+j*4;
  t_c[0]=q[0];
  r_c[0]=q[1];
  lon_c[0]=q[2];
  lat_c[0]=q[3];
 }

 return N;
}

int SphGrid :: FindCrossingLat(int k, double *t_c, double *r_c, double *lon_c, double *lat_c)
{
 int N;

 if (lat_cN[k]<0)
 {
  N=0;

  double st=sin(lat_arr[k]);
  double ct=cos(lat_arr[k]);
  double st2=sqr(st);
  double ct2=sqr(ct);

  if (st==0)
  {
   double t=-y1/ay;

   if (finite(t)) if (t>=0 && t<=1)
   {
	double x=x1+ax*t;
    double y=y1+ay*t;
    double z=z1+az*t;
    double r=sqrt(sqr(x)+sqr(y)+sqr(z));
    double phi=atan2(x, z);

	double *q=lat_cL+k*8+N*4;
	q[0]=t;
	q[1]=r;
	q[2]=phi;
	q[3]=lat_arr[k];

	N++;
   }
  }
  else
  {
   double A=sqr(ay)*ct2-(sqr(ax)+sqr(az))*st2;
   double B=2.0*(ay*y1*ct2-(ax*x1+az*z1)*st2);
   double C=sqr(y1)*ct2-(sqr(x1)+sqr(z1))*st2;

   double D=sqr(B)-4.0*A*C;
   if (D>0)
   {
	D=sqrt(D);

	for (int s=-1; s<=1; s+=2)
	{
	 double t=(-B+D*s)/(2.0*A);

	 if (finite(t)) if (t>=0 && t<=1)
	 {
	  double y=y1+ay*t;

	  if ((y*st)>0)
	  {
	   double x=x1+ax*t;
       double z=z1+az*t;
       double r=sqrt(sqr(x)+sqr(y)+sqr(z));
       double phi=atan2(x, z);

	   double *q=lat_cL+k*8+N*4;
	   q[0]=t;
	   q[1]=r;
	   q[2]=phi;
	   q[3]=lat_arr[k];

	   N++;
	  }
	 }
	}
   }
  }

  lat_cN[k]=N;
 }

 N=lat_cN[k];
 if (N>0) for (int m=0; m<N; m++)
 {
  double *q=lat_cL+k*8+m*4;
  t_c[m]=q[0];
  r_c[m]=q[1];
  lon_c[m]=q[2];
  lat_c[m]=q[3];
 }

 return N;
}

int SphGrid :: FindCrossings(int i1, int i2, int j1, int j2, int k1, int k2, double *t_arr)
{
 int N=0;
 
 double t_c[2], r_c[2], lon_c[2], lat_c[2];
 int Nl;

 Nl=FindCrossingR(i1, t_c, r_c, lon_c, lat_c);
 if (Nl>0) for (int m=0; m<Nl; m++) 
  if ((j2>j1 ? (lon_c[m]>=lon_arr[j1] && lon_c[m]<lon_arr[j2]) : 
	           (lon_c[m]>=lon_arr[j1] || lon_c[m]<lon_arr[j2])) && 
	  (lat_c[m]>=lat_arr[k1] && lat_c[m]<lat_arr[k2]))
  {
   t_arr[N]=t_c[m];
   N++;
  }

 Nl=FindCrossingR(i2, t_c, r_c, lon_c, lat_c);
 if (Nl>0) for (int m=0; m<Nl; m++) 
  if ((j2>j1 ? (lon_c[m]>=lon_arr[j1] && lon_c[m]<lon_arr[j2]) : 
	           (lon_c[m]>=lon_arr[j1] || lon_c[m]<lon_arr[j2])) && 
	  (lat_c[m]>=lat_arr[k1] && lat_c[m]<lat_arr[k2]))
  {
   t_arr[N]=t_c[m];
   N++;
  }

 Nl=FindCrossingLon(j1, t_c, r_c, lon_c, lat_c);
 if (Nl>0) 
  if ((r_c[0]>=r_arr[i1] && r_c[0]<r_arr[i2]) && (lat_c[0]>=lat_arr[k1] && lat_c[0]<lat_arr[k2]))
  {
   t_arr[N]=t_c[0];
   N++;
  }

 Nl=FindCrossingLon(j2, t_c, r_c, lon_c, lat_c);
 if (Nl>0) 
  if ((r_c[0]>=r_arr[i1] && r_c[0]<r_arr[i2]) && (lat_c[0]>=lat_arr[k1] && lat_c[0]<lat_arr[k2]))
  {
   t_arr[N]=t_c[0];
   N++;
  }

 Nl=FindCrossingLat(k1, t_c, r_c, lon_c, lat_c);
 if (Nl>0) for (int m=0; m<Nl; m++) 
  if ((r_c[m]>=r_arr[i1] && r_c[m]<r_arr[i2]) && 
	  (j2>j1 ? (lon_c[m]>=lon_arr[j1] && lon_c[m]<lon_arr[j2]) : 
	           (lon_c[m]>=lon_arr[j1] || lon_c[m]<lon_arr[j2])))
  {
   t_arr[N]=t_c[m];
   N++;
  }  

 Nl=FindCrossingLat(k2, t_c, r_c, lon_c, lat_c);
 if (Nl>0) for (int m=0; m<Nl; m++) 
  if ((r_c[m]>=r_arr[i1] && r_c[m]<r_arr[i2]) && 
	  (j2>j1 ? (lon_c[m]>=lon_arr[j1] && lon_c[m]<lon_arr[j2]) : 
	           (lon_c[m]>=lon_arr[j1] || lon_c[m]<lon_arr[j2])))
  {
   t_arr[N]=t_c[m];
   N++;
  }  

 return N;
}

void SphGrid :: FindIndices(double tc, int i1, int i2, int j1, int j2, int k1, int k2,
	                        double *r_idx, double *lon_idx, double *lat_idx)
{
 double xc=x1+ax*tc;
 double yc=y1+ay*tc;
 double zc=z1+az*tc;
 double r_c=sqrt(sqr(xc)+sqr(yc)+sqr(zc));
 double lon_c=atan2(xc, zc);
 if (j2<j1 && lon_c<0) lon_c+=(M_PI*2);
 double lat_c=asin(yc/r_c);
      
 *r_idx=i1+(r_c-r_arr[i1])/(r_arr[i2]-r_arr[i1]);
 *lon_idx=j2>j1 ? j1+(lon_c-lon_arr[j1])/(lon_arr[j2]-lon_arr[j1]) :
	              j1+(lon_c-lon_arr[j1])/(lon_arr[j2]-lon_arr[j1]+M_PI*2);
 *lat_idx=k1+(lat_c-lat_arr[k1])/(lat_arr[k2]-lat_arr[k1]);
}

class VoxelStack
{
 int *s;
 int N, Nmax;
 public:
 VoxelStack();
 ~VoxelStack();
 void push(int i1, int i2, int j1, int j2, int k1, int k2);
 void pop(int *i1, int *i2, int *j1, int *j2, int *k1, int *k2);
 int empty();
};

VoxelStack :: VoxelStack()
{
 s=0;
 N=Nmax=0;
}

VoxelStack :: ~VoxelStack()
{
 if (s) free(s);
}

void VoxelStack :: push(int i1, int i2, int j1, int j2, int k1, int k2)
{
 if (N>=Nmax)
 {
  Nmax+=64;
  s=s ? (int*)realloc(s, sizeof(int)*Nmax*6) : (int*)malloc(sizeof(int)*Nmax*6);
 }

 int *q=s+N*6;
 q[0]=i1;
 q[1]=i2;
 q[2]=j1;
 q[3]=j2;
 q[4]=k1;
 q[5]=k2;

 N++;
}

int VoxelStack :: empty()
{
 return N==0;
}

void VoxelStack :: pop(int *i1, int *i2, int *j1, int *j2, int *k1, int *k2)
{
 N--;

 int *q=s+N*6;
 *i1=q[0];
 *i2=q[1];
 *j1=q[2];
 *j2=q[3];
 *k1=q[4];
 *k2=q[5];
}

class PixelStack
{
 int *idx_arr;
 double *tc_arr, *dl_arr, *r_idx_arr, *lon_idx_arr, *lat_idx_arr;
 int N, Nmax;
 public:
 PixelStack();
 ~PixelStack();
 void push(int idx, double tc, double dl, double r_idx, double lon_idx, double lat_idx);
 int exportS(int Bsize, int *index, double *dl, double *r_ind, double *lon_ind, double *lat_ind);
 void Sort();
};

PixelStack :: PixelStack()
{
 idx_arr=0;
 tc_arr=dl_arr=r_idx_arr=lon_idx_arr=lat_idx_arr=0;
 N=Nmax=0;
}

PixelStack :: ~PixelStack()
{
 if (idx_arr) free(idx_arr);
 if (tc_arr) free(tc_arr);
 if (dl_arr) free(dl_arr);
 if (r_idx_arr) free(r_idx_arr);
 if (lon_idx_arr) free(lon_idx_arr);
 if (lat_idx_arr) free(lat_idx_arr);
}

void PixelStack :: push(int idx, double tc, double dl, double r_idx, double lon_idx, double lat_idx)
{
 if (N>=Nmax)
 {
  Nmax+=64;
  idx_arr=idx_arr ? (int*)realloc(idx_arr, sizeof(int)*Nmax) : (int*)malloc(sizeof(int)*Nmax);
  tc_arr=tc_arr ? (double*)realloc(tc_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
  dl_arr=dl_arr ? (double*)realloc(dl_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
  r_idx_arr=r_idx_arr ? (double*)realloc(r_idx_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
  lon_idx_arr=lon_idx_arr ? (double*)realloc(lon_idx_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
  lat_idx_arr=lat_idx_arr ? (double*)realloc(lat_idx_arr, sizeof(double)*Nmax) : (double*)malloc(sizeof(double)*Nmax);
 }

 idx_arr[N]=idx;
 tc_arr[N]=tc;
 dl_arr[N]=dl;
 r_idx_arr[N]=r_idx;
 lon_idx_arr[N]=lon_idx;
 lat_idx_arr[N]=lat_idx;

 N++;
}

void PixelStack :: Sort()
{
 if (N)
 {
  int *s=(int*)malloc(sizeof(int)*N);
  sort(tc_arr, N, s);

  int *_idx_arr=(int*)malloc(sizeof(int)*Nmax);
  double *_dl_arr=(double*)malloc(sizeof(double)*Nmax);
  double *_r_idx_arr=(double*)malloc(sizeof(double)*Nmax);
  double *_lon_idx_arr=(double*)malloc(sizeof(double)*Nmax);
  double *_lat_idx_arr=(double*)malloc(sizeof(double)*Nmax);

  for (int i=0; i<N; i++)
  {
   int j=s[i];
   _idx_arr[i]=idx_arr[j];
   _dl_arr[i]=dl_arr[j];
   _r_idx_arr[i]=r_idx_arr[j];
   _lon_idx_arr[i]=lon_idx_arr[j];
   _lat_idx_arr[i]=lat_idx_arr[j];
  }

  free(idx_arr);
  free(dl_arr);
  free(r_idx_arr);
  free(lon_idx_arr);
  free(lat_idx_arr);

  idx_arr=_idx_arr;
  dl_arr=_dl_arr;
  r_idx_arr=_r_idx_arr;
  lon_idx_arr=_lon_idx_arr;
  lat_idx_arr=_lat_idx_arr;

  free(s);
 }
}

int PixelStack :: exportS(int Bsize, int *index, double *dl, double *r_ind, double *lon_ind, double *lat_ind)
{
 int Nvox=N;
 int N0=0;

 if (Nvox) for (int k=Nvox-1; k>=0; k--) if (idx_arr[k]==-1)
 {
  N0=k+1;
  break;
 }
 Nvox-=N0;

 if (Nvox) if (Nvox<=Bsize)
 {
  memcpy(index, idx_arr+N0, sizeof(int)*Nvox);
  memcpy(dl, dl_arr+N0, sizeof(double)*Nvox);
  memcpy(r_ind, r_idx_arr+N0, sizeof(double)*Nvox);
  memcpy(lon_ind, lon_idx_arr+N0, sizeof(double)*Nvox);
  memcpy(lat_ind, lat_idx_arr+N0, sizeof(double)*Nvox);
 }

 return Nvox;
}

int RenderSpherical(int Nlon, int Nlat, int Nr, double *lon_arr, double *lat_arr, double *r_arr,
	                double *p1, double *p2, double rmin, int closed, int Nmax, 
	                int *index, double *dl, double *lon_ind, double *lat_ind, double *r_ind)
{
 SphGrid *G=new SphGrid(Nr, Nlon, Nlat, r_arr, lon_arr, lat_arr, p1, p2);
 VoxelStack *vS=new VoxelStack();
 PixelStack *pS=new PixelStack();

 vS->push(0, Nr-1, 0, Nlon-1, 0, Nlat-1);
 if (closed) vS->push(0, Nr-1, Nlon-1, 0, 0, Nlat-1);

 while (!vS->empty())
 {
  int i1, i2, j1, j2, k1, k2;
  double t_arr[6];

  vS->pop(&i1, &i2, &j1, &j2, &k1, &k2);
  int N=G->FindCrossings(i1, i2, j1, j2, k1, k2, t_arr);

  if (N>0) if ((i2-i1)==1 && ((j2-j1)==1 || j2<j1) && (k2-k1)==1)
  {
   if (N%2==0)
   {
	sort(t_arr, N, 0);

	int idx=j1+k1*Nlon+i1*Nlon*Nlat;

	for (int m=0; m<N/2; m++)
	{
     double tc=(t_arr[m*2]+t_arr[m*2+1])/2;
     double dl=(t_arr[m*2+1]-t_arr[m*2])*G->L;

	 double r_idx, lon_idx, lat_idx;
	 G->FindIndices(tc, i1, i2, j1, j2, k1, k2, &r_idx, &lon_idx, &lat_idx);

	 pS->push(idx, tc, dl, r_idx, lon_idx, lat_idx);
	}
   }
  }
  else
  {
   int ic=(i1+i2)/2;
   int kc=(k1+k2)/2;

   if (j2>j1)
   {
    int jc=(j1+j2)/2;

    if (i1!=ic && j1!=jc && k1!=kc) vS->push(i1, ic, j1, jc, k1, kc);
    if (i1!=ic && j1!=jc && kc!=k2) vS->push(i1, ic, j1, jc, kc, k2);
    if (i1!=ic && jc!=j2 && k1!=kc) vS->push(i1, ic, jc, j2, k1, kc);
    if (i1!=ic && jc!=j2 && kc!=k2) vS->push(i1, ic, jc, j2, kc, k2);
    if (ic!=i2 && j1!=jc && k1!=kc) vS->push(ic, i2, j1, jc, k1, kc);
    if (ic!=i2 && j1!=jc && kc!=k2) vS->push(ic, i2, j1, jc, kc, k2);
    if (ic!=i2 && jc!=j2 && k1!=kc) vS->push(ic, i2, jc, j2, k1, kc);
    if (ic!=i2 && jc!=j2 && kc!=k2) vS->push(ic, i2, jc, j2, kc, k2);
   }
   else
   {
	if (i1!=ic && k1!=kc) vS->push(i1, ic, j1, j2, k1, kc);
	if (i1!=ic && kc!=k2) vS->push(i1, ic, j1, j2, kc, k2);
	if (ic!=i2 && k1!=kc) vS->push(ic, i2, j1, j2, k1, kc);
	if (ic!=i2 && kc!=k2) vS->push(ic, i2, j1, j2, kc, k2);
   }
  }
 }

 if (rmin>0)
 {
  double t_arr[2];
  int N=G->FindCrossingRmin((rmin<r_arr[0]) ? rmin : r_arr[0], t_arr);

  if (N==2)
  {
   sort(t_arr, 2, 0);
   double tc=(t_arr[0]+t_arr[1])/2;
   double dl=(t_arr[1]-t_arr[0])*G->L;

   pS->push(-1, tc, dl, -1, -1, -1);
  }
 }

 pS->Sort();
 int Nvox=pS->exportS(Nmax, index, dl, r_ind, lon_ind, lat_ind);

 delete G;
 delete vS; 
 delete pS;
 
 return Nvox; 
}