pro Example_MultiThread
 ;Load the PFSS extrapolation
 restore, 'PFSSmodel1.sav' ; => grid points (lon, lat, r) +
                           ; magnetic field components (Bph, Bth, Br) at the nodes
 
 ;Rotate the model by 50 degrees
 lon+=50

 ;--------------------------------------------------------------------------------- 
 ;Draw the Sun
 InitPlot, G, 2.6
 PlotSphere, G, [0d0, 0d0, 0d0], 1d0, 'yellow', 'gray'
 
 ;Draw the boundaries of the model volume
 PlotSphericalBox, G, min(r), max(r), min(lon), max(lon), min(lat), max(lat), 'black'
 
 ;---------------------------------------------------------------------------------
 NLOS=4L ;number of lines of sight
 
 ;Endpoints of the lines of sight (in Cartesian system):
 r1_m=[[0.65, -0.3, -2.5], $
       [0.70, -0.3, -2.5], $
       [0.75, -0.3, -2.5], $
       [0.80, -0.3, -2.5]]
 r2_m=[[0.65, -0.3,  2.5], $
       [0.70, -0.3,  2.5], $
       [0.75, -0.3,  2.5], $
       [0.80, -0.3,  2.5]]
 
 b0=5.205*!dpi/180 ;solar b0 angle
 
 ;Correct the lines of sight to account for the b0 angle
 ;(it is more convenient than to rotate the model around x axis)
 for j=0, NLOS-1 do begin
  r1_m[*, j]=[r1_m[0, j], r1_m[1, j]*cos(b0)+r1_m[2, j]*sin(b0), -r1_m[1, j]*sin(b0)+r1_m[2, j]*cos(b0)]
  r2_m[*, j]=[r2_m[0, j], r2_m[1, j]*cos(b0)+r2_m[2, j]*sin(b0), -r2_m[1, j]*sin(b0)+r2_m[2, j]*cos(b0)]
 endfor 
 
 a_m=r2_m-r1_m 
 for j=0, NLOS-1 do $
  a_m[*, j]=a_m[*, j]/sqrt(a_m[0, j]^2+a_m[1, j]^2+a_m[2, j]^2) ;unit vectors along the lines of sight
 
 ;Draw the lines of sight
 for j=0, NLOS-1 do $
  G=polyline([r1_m[0, j], r2_m[0, j]], [r1_m[1, j], r2_m[1, j]], [r1_m[2, j], r2_m[2, j]], $
             /data, overplot=G, thick=2, color='cyan')
 
 ;Convert the lines-of-sight endpoints to spherical system
 p1_m=dblarr(3, NLOS)
 p2_m=dblarr(3, NLOS)
 for j=0, NLOS-1 do begin
  p1_m[*, j]=[sqrt(r1_m[0, j]^2+r1_m[1, j]^2+r1_m[2, j]^2), $
              atan(r1_m[0, j], r1_m[2, j])/!dpi*180, $
              asin(r1_m[1, j]/sqrt(r1_m[0, j]^2+r1_m[1, j]^2+r1_m[2, j]^2))/!dpi*180]
  p2_m[*, j]=[sqrt(r2_m[0, j]^2+r2_m[1, j]^2+r2_m[2, j]^2), $
              atan(r2_m[0, j], r2_m[2, j])/!dpi*180, $
              asin(r2_m[1, j]/sqrt(r2_m[0, j]^2+r2_m[1, j]^2+r2_m[2, j]^2))/!dpi*180]
 endfor             
 
 forward_function RenderSphericalMulti
 
 ;Call the rendering code
 Nvox_m=RenderSphericalMulti(NLOS, lon, lat, r, p1_m, p2_m, min(r), $
                             index_m, dl_m, lon_ind_m, lat_ind_m, r_ind_m)  
 
 Nmax=max(Nvox_m) ;maximum number of voxels - to create arrays for B, theta, psi
 
 if Nmax gt 0 then begin
  B=dblarr(Nmax, NLOS)
  theta=dblarr(Nmax, NLOS)
  psi=dblarr(Nmax, NLOS)
  l=dblarr(Nmax, NLOS)
 
  for j=0, NLOS-1 do if Nvox_m[j] gt 0 then begin
   ;Compute the interpolated magnetic field values along the line of sight
   r_loc=interpolate(r, r_ind_m[0 : Nvox_m[j]-1, j]) ;radial distance
   lon_loc=interpolate(lon, lon_ind_m[0 : Nvox_m[j]-1, j]) ;longitude
   lat_loc=interpolate(lat, lat_ind_m[0 : Nvox_m[j]-1, j]) ;latitude
   br_loc=interpolate(br, lon_ind_m[0 : Nvox_m[j]-1, j], lat_ind_m[0 : Nvox_m[j]-1, j], $
                      r_ind_m[0 : Nvox_m[j]-1, j]) ;B_r
   bph_loc=interpolate(bph, lon_ind_m[0 : Nvox_m[j]-1, j], lat_ind_m[0 : Nvox_m[j]-1, j], $
                       r_ind_m[0 : Nvox_m[j]-1, j]) ;B_ph=B_lon
   bth_loc=interpolate(bth, lon_ind_m[0 : Nvox_m[j]-1, j], lat_ind_m[0 : Nvox_m[j]-1, j], $
                       r_ind_m[0 : Nvox_m[j]-1, j]) ;B_th=B_lat
  
   ;Convert the field components to Cartesian system
   sp=sin(lon_loc*!dpi/180)
   cp=cos(lon_loc*!dpi/180)
   st=sin(lat_loc*!dpi/180)
   ct=cos(lat_loc*!dpi/180)
   Bx=br_loc*sp*ct+bph_loc*cp-bth_loc*sp*st
   By=br_loc*st+bth_loc*ct
   Bz=br_loc*cp*ct-bph_loc*sp-bth_loc*cp*st
  
   B[0 : Nvox_m[j]-1, j]=sqrt(Bx^2+By^2+Bz^2) ;absolute value of the magnetic field
 
   a=reform(a_m[*, j]) ;unit vector along the current line of sight
   ;viewing angle
   theta[0 : Nvox_m[j]-1, j]=acos((a[0]*Bx+a[1]*By+a[2]*Bz)/B[0 : Nvox_m[j]-1, j])/!dpi*180 
  
   n1=[0d0, a[2], -a[1]]
   n1=n1/sqrt(n1[0]^2+n1[1]^2+n1[2]^2) ;unit vector perpendicular to the line of sight
   n2=[-a[1]^2-a[2]^2, a[0]*a[1], a[0]*a[2]]
   n2=n2/sqrt(n2[0]^2+n2[1]^2+n2[2]^2) ;unit vector perpendicular to both the line of sight and n1
   ;azimuthal angle
   psi[0 : Nvox_m[j]-1, j]=atan(n1[0]*Bx+n1[1]*By+n1[2]*Bz, n2[0]*Bx+n2[1]*By+n2[2]*Bz)/!dpi*180 
   
   l[0 : Nvox_m[j]-1, j]=total(dl_m[0 : Nvox_m[j]-1, j], /cumulative) ;distance along the line of sight
  endif
  
  ;Plot the magnetic field strength vs. coordinate along the line of sight
  window, 1, title='Magnetic field strength'
  wset, 1
  for j=0, NLOS-1 do if Nvox_m[j] gt 0 then begin
   if j eq 0 then plot, l[0 : Nvox_m[j]-1, j], B[0 : Nvox_m[j]-1, j], linestyle=j, $
                        xrange=[min(l), max(l)], yrange=[min(B), max(B)], $
                        xtitle='Distance, R_Sun', ytitle='!8B!3, G' $
   else oplot, l[0 : Nvox_m[j]-1, j], B[0 : Nvox_m[j]-1, j], linestyle=j
  endif      
        
  ;Plot the viewing angle vs. coordinate along the line of sight
  window, 2, title='Viewing angle'
  wset, 2
  for j=0, NLOS-1 do if Nvox_m[j] gt 0 then begin
   if j eq 0 then plot, l[0 : Nvox_m[j]-1, j], theta[0 : Nvox_m[j]-1, j], linestyle=j, $
                        xrange=[min(l), max(l)], yrange=[min(theta), max(theta)], $
                        xtitle='Distance, R_Sun', ytitle='!7h!3, degrees' $
   else oplot, l[0 : Nvox_m[j]-1, j], theta[0 : Nvox_m[j]-1, j], linestyle=j
  endif   
        
  ;Plot the azimuthal angle vs. coordinate along the line of sight
  window, 3, title='Azimuthal angle'
  wset, 3
  for j=0, NLOS-1 do if Nvox_m[j] gt 0 then begin
   if j eq 0 then plot, l[0 : Nvox_m[j]-1, j], psi[0 : Nvox_m[j]-1, j], linestyle=j, $
                        xrange=[min(l), max(l)], yrange=[min(psi), max(psi)], $
                        xtitle='Distance, R_Sun', ytitle='!7w!3, degrees' $
   else oplot, l[0 : Nvox_m[j]-1, j], psi[0 : Nvox_m[j]-1, j], linestyle=j
  endif
 endif
end

pro InitPlot, G, L
 G=surface(dblarr(10, 10), dblarr(10, 10), dblarr(10, 10), xrange=[-L, L], yrange=[-L, L], zrange=[-L, L], $
           axis_style=0, aspect_z=1, aspect_ratio=1, nodata=1, dimensions=[800, 800])     
 if float(!version.release le 8.2) then begin ;for some reason rotation works differently in different IDL versions
  G.Rotate, -30, /xaxis           
  G.Rotate, -120, /yaxis
  G.Rotate, 90, /xaxis
  G.Rotate, 90, /zaxis
 endif else begin
  G.Rotate, 90, /xaxis
  G.Rotate, -30, /yaxis
  G.Rotate, -30, /xaxis
 endelse
 ;basic geometry with z axis to us, x axis to right and y axis to top
end

pro PlotSphere, G, r0, R, color1, color2
 x=dblarr(100, 100)
 y=dblarr(100, 100)
 z=dblarr(100, 100)
 z[*, *]=!values.d_nan
 for i=0, 99 do for j=0, 99 do begin
  theta=!dpi/2*i/99
  q=sin(theta)
  phi=2d0*!dpi*j/99
  x[i, j]=q*cos(phi)
  y[i, j]=q*sin(phi)
  z[i, j]=sqrt((1d0-x[i, j]^2-y[i, j]^2)>0)
 endfor
 x*=R
 y*=R
 z*=R
 G=surface( z+r0[2], x+r0[0], y+r0[1], overplot=G, color=color1)
 G=surface(-z+r0[2], x+r0[0], y+r0[1], overplot=G, color=color1)
 
 phi=!dpi*2*dindgen(360)/359
 G=polyline(cos(phi), sin(phi), replicate(0, 360), /data, color=color2, linestyle=2, overplot=G)
 G=polyline(cos(phi), replicate(0, 360), sin(phi), /data, color=color2, linestyle=2, overplot=G)
 G=polyline(replicate(0, 360), cos(phi), sin(phi), /data, color=color2, linestyle=2, overplot=G)           
end

pro PlotSphericalBox, G, r1, r2, lon1, lon2, lat1, lat2, color
 lon1*=(!dpi/180)
 lat1*=(!dpi/180)
 lon2*=(!dpi/180)
 lat2*=(!dpi/180)
  
 G=polyline(cos(lat1)*sin(lon1)*[r1, r2], $
            sin(lat1)*[r1, r2], $
            cos(lat1)*cos(lon1)*[r1, r2], $
            overplot=G, /data, color=color)
 G=polyline(cos(lat1)*sin(lon2)*[r1, r2], $
            sin(lat1)*[r1, r2], $
            cos(lat1)*cos(lon2)*[r1, r2], $
            overplot=G, /data, color=color)     
 G=polyline(cos(lat2)*sin(lon1)*[r1, r2], $
            sin(lat2)*[r1, r2], $
            cos(lat2)*cos(lon1)*[r1, r2], $
            overplot=G, /data, color=color)
 G=polyline(cos(lat2)*sin(lon2)*[r1, r2], $
            sin(lat2)*[r1, r2], $
            cos(lat2)*cos(lon2)*[r1, r2], $
            overplot=G, /data, color=color)    
            
 phi=lon1+(lon2-lon1)*dindgen(100)/99
 G=polyline(cos(lat1)*sin(phi)*r1, $
            replicate(sin(lat1)*r1, 100), $
            cos(lat1)*cos(phi)*r1, $
            overplot=G, /data, color=color)
 G=polyline(cos(lat1)*sin(phi)*r2, $
            replicate(sin(lat1)*r2, 100), $
            cos(lat1)*cos(phi)*r2, $
            overplot=G, /data, color=color)
 G=polyline(cos(lat2)*sin(phi)*r1, $
            replicate(sin(lat2)*r1, 100), $
            cos(lat2)*cos(phi)*r1, $
            overplot=G, /data, color=color)
 G=polyline(cos(lat2)*sin(phi)*r2, $
            replicate(sin(lat2)*r2, 100), $
            cos(lat2)*cos(phi)*r2, $
            overplot=G, /data, color=color)       
            
 theta=lat1+(lat2-lat1)*dindgen(100)/99
 G=polyline(cos(theta)*sin(lon1)*r1, $
            sin(theta)*r1, $
            cos(theta)*cos(lon1)*r1, $
            overplot=G, /data, color=color)      
 G=polyline(cos(theta)*sin(lon1)*r2, $
            sin(theta)*r2, $
            cos(theta)*cos(lon1)*r2, $
            overplot=G, /data, color=color)
 G=polyline(cos(theta)*sin(lon2)*r1, $
            sin(theta)*r1, $
            cos(theta)*cos(lon2)*r1, $
            overplot=G, /data, color=color)      
 G=polyline(cos(theta)*sin(lon2)*r2, $
            sin(theta)*r2, $
            cos(theta)*cos(lon2)*r2, $
            overplot=G, /data, color=color)
end