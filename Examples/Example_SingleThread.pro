pro Example_SingleThread
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
 ;Endpoints of the line of sight (in Cartesian system):
 r1=[0.75, -0.3, -2.5]
 r2=[0.75, -0.3,  2.5]
 
 b0=5.205*!dpi/180 ;solar b0 angle
 
 ;Correct the line of sight to account for the b0 angle
 ;(it is more convenient than to rotate the model around x axis)
 r1=[r1[0], r1[1]*cos(b0)+r1[2]*sin(b0), -r1[1]*sin(b0)+r1[2]*cos(b0)]
 r2=[r2[0], r2[1]*cos(b0)+r2[2]*sin(b0), -r2[1]*sin(b0)+r2[2]*cos(b0)]
 
 a=r2-r1 
 a=a/sqrt(a[0]^2+a[1]^2+a[2]^2) ;unit vector along the line of sight
 
 ;Draw the line of sight
 G=polyline([r1[0], r2[0]], [r1[1], r2[1]], [r1[2], r2[2]], /data, overplot=G, thick=2, color='cyan')
 
 ;Convert the line-of-sight endpoints to spherical system
 p1=[sqrt(r1[0]^2+r1[1]^2+r1[2]^2), $
     atan(r1[0], r1[2])/!dpi*180, $
     asin(r1[1]/sqrt(r1[0]^2+r1[1]^2+r1[2]^2))/!dpi*180]
 p2=[sqrt(r2[0]^2+r2[1]^2+r2[2]^2), $
     atan(r2[0], r2[2])/!dpi*180, $
     asin(r2[1]/sqrt(r2[0]^2+r2[1]^2+r2[2]^2))/!dpi*180]
 
 forward_function RenderSpherical
 
 ;Call the rendering code
 Nvox=RenderSpherical(lon, lat, r, p1, p2, min(r), index, dl, lon_ind, lat_ind, r_ind)  
 
 if Nvox gt 0 then begin
  ;Compute the interpolated magnetic field values along the line of sight
  r_loc=interpolate(r, r_ind) ;radial distance
  lon_loc=interpolate(lon, lon_ind) ;longitude
  lat_loc=interpolate(lat, lat_ind) ;latitude
  br_loc=interpolate(br, lon_ind, lat_ind, r_ind) ;B_r
  bph_loc=interpolate(bph, lon_ind, lat_ind, r_ind) ;B_ph=B_lon
  bth_loc=interpolate(bth, lon_ind, lat_ind, r_ind) ;B_th=B_lat
  
  ;Convert the field components to Cartesian system
  sp=sin(lon_loc*!dpi/180)
  cp=cos(lon_loc*!dpi/180)
  st=sin(lat_loc*!dpi/180)
  ct=cos(lat_loc*!dpi/180)
  Bx=br_loc*sp*ct+bph_loc*cp-bth_loc*sp*st
  By=br_loc*st+bth_loc*ct
  Bz=br_loc*cp*ct-bph_loc*sp-bth_loc*cp*st
  
  B=sqrt(Bx^2+By^2+Bz^2) ;absolute value of the magnetic field
  
  theta=acos((a[0]*Bx+a[1]*By+a[2]*Bz)/B)/!dpi*180 ;viewing angle
  
  n1=[0d0, a[2], -a[1]]
  n1=n1/sqrt(n1[0]^2+n1[1]^2+n1[2]^2) ;unit vector perpendicular to the line of sight
  n2=[-a[1]^2-a[2]^2, a[0]*a[1], a[0]*a[2]]
  n2=n2/sqrt(n2[0]^2+n2[1]^2+n2[2]^2) ;unit vector perpendicular to both the line of sight and n1
  psi=atan(n1[0]*Bx+n1[1]*By+n1[2]*Bz, n2[0]*Bx+n2[1]*By+n2[2]*Bz)/!dpi*180 ;azimuthal angle
  
  ;Plot the magnetic field strength vs. coordinate along the line of sight
  window, 1, title='Magnetic field strength'
  wset, 1
  plot, total(dl, /cumulative), B, $
        xtitle='Distance, R_Sun', ytitle='!8B!3, G'
        
  ;Plot the viewing angle vs. coordinate along the line of sight
  window, 2, title='Viewing angle'
  wset, 2
  plot, total(dl, /cumulative), theta, $
        xtitle='Distance, R_Sun', ytitle='!7h!3, degrees'        
        
  ;Plot the azimuthal angle vs. coordinate along the line of sight
  window, 3, title='Azimuthal angle'
  wset, 3
  plot, total(dl, /cumulative), psi, $
        xtitle='Distance, R_Sun', ytitle='!7w!3, degrees'          
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