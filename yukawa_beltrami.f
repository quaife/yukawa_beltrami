      program yukawa_beltrami 
c     ------------------------------------------------------------------
c
c
c  Solves integral equation for Yukawa-Beltrami on the Sphere 
c
c  this solves the Dirichlet bvp in multiply connected domains
c  All layer potentials are evaluated directly which means that it
c  requires quadratic CPU time
c
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
c Maximum number of holes and points per hole
      parameter (kmax = 50, npmax = 512, nmax = kmax*npmax)
      parameter (nvortmax = 10)
c
c Geometry of holes
      dimension ak(kmax), bk(kmax), th_k(kmax), phi_k(kmax), cx(kmax),
     1          cy(kmax), cz(kmax)
      dimension xs(nmax), ys(nmax), zs(nmax), dx(nmax), dy(nmax),
     1          dz(nmax), xn(nmax), yn(nmax), zn(nmax)
      dimension dsda(nmax), arcLength(kmax), diag(nmax)
      dimension d2x(nmax), d2y(nmax), d2z(nmax)
c
c  Vortex variables
      dimension x1Vort(nvortmax),x2Vort(nvortmax),x3Vort(nvortmax)
c     vortex location
      dimension vortK(nvortmax)  
c     vortex strength

c
c  Grid variables
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      dimension igrid(ng_max), near_gr(ng_max*kmax), 
     1          th_gr(ng_max), phi_gr(ng_max),
     2          u_gr(ng_max), x_gr(ng_max), y_gr(ng_max), z_gr(ng_max),
     3          uExact_gr(ng_max)
c
c target points are used to check accuracy
      dimension xtar(ng_max), ytar(ng_max), ztar(ng_max)
c
c Density function
      dimension density(nmax)
c
c  Matrix equation variables for GMRES
c  maxl is the maximum nubmer of GMRES iterations performed
c       before restarting.
c  lrwork is the dimension of a real workspace needed by DGMRES.
c  liwork is the dimension of an integer workspace needed by DGMRES.
c  gmwork and iwork are work arrays used by DGMRES
      parameter (maxl = 200,liwork=30,  
     1           lrwork=10+(nmax+kmax)*(maxl+6)+maxl*(maxl+3))
      dimension gmwork(lrwork), igwork(liwork)
      dimension rhs(nmax)
c
c Variables to time different components
      real*4 timep(2), etime
c
c common blocks
      common /sys_size/ k, nd, nbk
      common /sphere_int/ xs, ys, zs, xn, yn, zn, dsda, diag, freq
c
c Initial Hole Geometry is given by reading in data
      call PRINI(6,13)
      call READ_DATA (freq, k, nd, nbk, nth, nphi, ak, bk, 
     1                      cx, cy, cz, th_k, phi_k,
     2                      nvort, x1Vort, x2Vort, x3Vort, vortK)
c
c Construct hole geometry and grid on surface of sphere
      call MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1               dx, dy, dz, d2x, d2y, d2z, xn, yn, zn, dsda, 
     2               arcLength, diag)

c Construct grid on surface of sphere
      call SURFACE_GRID (k, nd, nbk, freq, nth, nphi, ak, bk,
     1                   cx, cy, cz, th_k, phi_k, th_gr, phi_gr, 
     2                   arcLength,
     3                   nvort, x1Vort, x2Vort, x3Vort, vortK, 
     4                   x_gr, y_gr, z_gr, igrid, near_gr, uExact_gr)

c Construct the RHS and solve
      call getRHS (k, nd, nbk, nvort, freq, xs, ys, zs, 
     1    x1Vort, x2Vort, x3Vort, vortK, rhs)

c Find the density function defined on the boundaries of the
c geometry
      call solveBIE (nd, k, nbk, rhs, density, gmwork, 
     1            lrwork, igwork, liwork, dsda, maxl)

c Construct solution on surface grid
       call SOL_GRID (freq, nd, k, nbk, nth, nphi, density, 
     1        xs, ys, zs, xn, yn, zn, dsda, arcLength, 
     2        x_gr, y_gr, z_gr, igrid, near_gr, u_gr)

c Check the infinity norm of the error in the far, intermediate
c and near regions
       call CHECK_ERROR (k, nth, nphi, near_gr, u_gr, uExact_gr)

c Create a matlab file that plots the solution on the surface
c of the sphere
      open (unit = 35, file = 'surfPlot.m')
      call dump_movie (nth, nphi, x_gr, y_gr, z_gr, 
     1          u_gr, uExact_gr, 35)
      close(35)



      stop
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine READ_DATA (freq, k, nd, nbk, nth, nphi, ak, bk, 
     1                      cx, cy, cz, th_k, phi_k,
     2                      nvort, x1Vort, x2Vort, x3Vort, vortK)
c
c  *** DESCRIPTION :
c
c   Load in parameters for the geometry that will be defined on the
c   surface of the unit circle
c
c   *** INPUT PARAMETERS :
c
c   NONE
c
c   *** OUTPUT PARAMETERS :
c
c   freq    = constant in the pde (\Delta - freq^2)u = 0
c   k       = number of connected componenets in geometry
c   nd      = number of points per connected component
c   nbk     = total number of points k*nd
c   nth     = number of points in theta direction for meshgrid
c   nphi    = number of points in phi direction for meshgrid
c   ak      = parameter for ellipses on sphere
c   bk      = parameter for ellipses on sphere
c   cx      = x-coordinate of center ellipses on sphere
c   cy      = y-coordinate of center ellipses on sphere
c   cz      = z-coordinate of center ellipses on sphere
c   th_k    = parameter for ellipses on sphere
c   phi_k   = parameter for ellipses on sphere
c   nvort   = number of vorticies to form exact solution
c   x1Vort  = x-coordinate of the vorticies
c   x2Vort  = y-coordinate of the vorticies
c   x3Vort  = z-coordinate of the vorticies
c   vortK   = strength of the vorticies
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c 
      dimension ak(*), bk(*), cx(*), cy(*), cz(*), th_k(*), phi_k(*)
      dimension x1Vort(*),x2Vort(*),x3Vort(*),vortK(*)
c
      complex *16 eye
c
      eye = dcmplx(0.d0,1.d0)
c
      open (unit = 12, file = 'input.data')
c
c     read in number of holes and number of points per hole
      read (12,*) k, nd, nvort
c     total number of points and number of vorticies for forming
c     exact solution
      nbk = k*nd
      call PRINF (' nbk = *', nbk, 1)
c     read in the one parameter in the PDE
      read (12,*) freq
c     read in the number of points for discretizing the meshgrid
      read(12,*) nth, nphi
      do kbod = 1, k
        read(12,*) ak(kbod), bk(kbod), th_k(kbod), phi_k(kbod)
c       compute center of each componenet curve
        call sph2cart (th_k(kbod),phi_k(kbod), 1.d0, cx(kbod),
     1                 cy(kbod), cz(kbod))
      end do

      do ivort = 1,nvort
        read(12,*) theta,phi,vortK(ivort)
        call sph2cart (theta,phi,1.d0,x1Vort(ivort),
     1      x2Vort(ivort),x3Vort(ivort))
      enddo


c     stop if the parameter in the PDE is less than 1/2
c     Our fundamental solution changes to the other associated
c     Legendre function when this is the case ... future work
      if (freq .le. 5.d-1) then
        call PRIN2 (' freq must be greater than one half *', freq, 1)
        stop
      endif

      close(12)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine sph2cart (theta, phi, r, x, y, z)
c
c  *** DESCRIPTION :
c
c   Convert a spherical coordinate to a cartesian coordinate
c
c   *** INPUT PARAMETERS :
c
c   theta    = azimuthal coordinate
c   phi      = z-uthal coordinate
c   r        = radius (which should always be one for us)
c
c   *** OUTPUT PARAMETERS :
c
c   x        = x-coordinate
c   y        = y-coordinate
c   z        = z-coordinate
c
c***********************************************************************
c
      real *8 theta,phi,r
      real *8 x,y,z
c
      x = r*dcos(phi)*dcos(theta)
      y = r*dcos(phi)*dsin(theta)
      z = r*dsin(phi)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine cross (u, v, u_cross_v)
c
c  *** DESCRIPTION :
c
c   Compute cross product of two three-dimensional vectors
c
c   *** INPUT PARAMETERS :
c
c   u         = first vector
c   v         = second vector
c
c   *** OUTPUT PARAMETERS :
c
c   u_cross_v = cross product of u and v
c
c***********************************************************************
c
      real *8 u(3), v(3), u_cross_v(3)
c
      u_cross_v(1) = u(2)*v(3) - v(2)*u(3)
      u_cross_v(2) = v(1)*u(3) - u(1)*v(3)
      u_cross_v(3) = u(1)*v(2) - v(1)*u(2)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine dot (u, v, u_dot_v)
c
c  *** DESCRIPTION :
c      
c   Compute the dot product of two three-dimensional vectors
c
c   *** INPUT PARAMETERS :
c
c   u         = first vector
c   v         = second vector
c
c   *** OUTPUT PARAMETERS :
c
c   u_dot_v = dot product of u and v
c
c***********************************************************************
c
      real *8 u_dot_v, u(3), v(3)
c
      u_dot_v = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine r_func (alpha, a, b, th, phi, x, y, z)
c
c  *** DESCRIPTION :
c
c   Compute points on the unit sphere using a generic parameterization
c   of ellipses on the sphere
c
c   *** INPUT PARAMETERS :
c
c   alpha     = argument used to parameterize the boundaries
c   a         = parameter for ellipses on sphere
c   b         = parameter for ellipses on sphere
c   th        = parameter for ellipses on sphere
c   phi       = parameter for ellipses on sphere
c
c   *** OUTPUT PARAMETERS :
c
c   x         = x-coordinate of point on sphere
c   y         = y-coordinate of point on sphere
c   z         = z-coordinate of point on sphere
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension zaxis(3), xaxis(3), yaxis(3)
c
      pi = 4.d0*datan(1.d0)
c
      call sph2cart (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
      call sph2cart (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1               xaxis(3))
      call cross (zaxis, xaxis, yaxis)
      xp = a*dcos(alpha)
      yp = b*dsin(alpha)
      zp = dsqrt(abs(1.d0 - xp**2 - yp**2))
      x = xp*xaxis(1) + yp*yaxis(1) + zp*zaxis(1)
      y = xp*xaxis(2) + yp*yaxis(2) + zp*zaxis(2)
      z = xp*xaxis(3) + yp*yaxis(3) + zp*zaxis(3)
c
      return
      end ! d_func
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine dr_func (alpha, a, b, th, phi, dx, dy, dz)
c
c  *** DESCRIPTION :
c
c   Compute first derivative of ellipses on the unit sphere
c
c   *** INPUT PARAMETERS :
c
c   alpha     = argument used to parameterize the boundaries
c   a         = parameter for ellipses on sphere
c   b         = parameter for ellipses on sphere
c   th        = parameter for ellipses on sphere
c   phi       = parameter for ellipses on sphere
c
c   *** OUTPUT PARAMETERS :
c
c   dx        = derivative of x-coordinate of point on sphere
c   dy        = derivative of y-coordinate of point on sphere
c   dz        = derivative of z-coordinate of point on sphere
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension zaxis(3), xaxis(3), yaxis(3)
c
      pi = 4.d0*datan(1.d0)
c
      call sph2cart (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
      call sph2cart (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1               xaxis(3))
      call cross (zaxis, xaxis, yaxis)
      xp = a*dcos(alpha)
      dxp = -a*dsin(alpha)
      yp = b*dsin(alpha)
      dyp = b*dcos(alpha)
      zp = dsqrt(1.d0 - xp**2 - yp**2)
      dzp = (-xp*dxp-yp*dyp)/zp
      dx = dxp*xaxis(1) + dyp*yaxis(1) + dzp*zaxis(1)
      dy = dxp*xaxis(2) + dyp*yaxis(2) + dzp*zaxis(2)
      dz = dxp*xaxis(3) + dyp*yaxis(3) + dzp*zaxis(3)
c
      return
      end ! dr_func   
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine d2r_func (alpha, a, b, th, phi, d2x, d2y, d2z)
c
c  *** DESCRIPTION :
c
c   Compute second derivative of ellipses on the unit sphere
c
c   *** INPUT PARAMETERS :
c
c   alpha     = argument used to parameterize the boundaries
c   a         = parameter for ellipses on sphere
c   b         = parameter for ellipses on sphere
c   th        = parameter for ellipses on sphere
c   phi       = parameter for ellipses on sphere
c
c   *** OUTPUT PARAMETERS :
c
c   d2x       = second derivative of x-coordinate of point on sphere
c   d2y       = second derivative of y-coordinate of point on sphere
c   d2z       = second derivative of z-coordinate of point on sphere
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension zaxis(3), xaxis(3), yaxis(3)
c
      pi = 4.d0*datan(1.d0)
c
      call sph2cart (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
      call sph2cart (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1               xaxis(3))
      call cross (zaxis, xaxis, yaxis)
      xp = a*dcos(alpha)
      dxp = -a*dsin(alpha)
      d2xp = -a*dcos(alpha)
      yp = b*dsin(alpha)
      dyp = b*dcos(alpha)
      d2yp = -b*dsin(alpha)
      zp = dsqrt(1.d0 - xp**2 - yp**2)
      dzp = (-xp*dxp-yp*dyp)/zp
      d2zp = (-dxp**2 - xp*d2xp - dyp**2 - yp*d2yp - dzp**2)/zp
      d2x = d2xp*xaxis(1) + d2yp*yaxis(1) + d2zp*zaxis(1)
      d2y = d2xp*xaxis(2) + d2yp*yaxis(2) + d2zp*zaxis(2)
      d2z = d2xp*xaxis(3) + d2yp*yaxis(3) + d2zp*zaxis(3)
c
      return
      end ! d2r_func   
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1                     dx, dy, dz, d2x, d2y, d2z, xn, yn, zn, dsda, 
     2                     arcLength, diag)
c
c  *** DESCRIPTION :
c
c   Compute points, derivatives, curvature, etc, and the diagonal term
c   of the second-kind integral equation for generic ellipses defined
c   on the boundary of the sphere
c
c   *** INPUT PARAMETERS :
c
c   k          = number of connected componenets in geometry
c   nd         = number of points per connected component
c   nbk        = total number of points k*nd
c   ak         = parameter for ellipses on sphere
c   bk         = parameter for ellipses on sphere
c   th_k       = parameter for ellipses on sphere
c   phi_k      = parameter for ellipses on sphere
c
c   *** OUTPUT PARAMETERS :
c
c   xs         = x-coordinate of the ellipses
c   ys         = y-coordinate of the ellipses
c   zs         = z-coordinate of the ellipses
c   dx         = first derivative of x-coordinate of the ellipses
c   dy         = first derivative of y-coordinate of the ellipses
c   dz         = first derivative of z-coordinate of the ellipses
c   d2x        = second derivative of x-coordinate of the ellipses
c   d2y        = second derivative of y-coordinate of the ellipses
c   d2z        = second derivative of z-coordinate of the ellipses
c   xn         = x-coordinate of the vector normal to the curve, but
c                tangent to the sphere
c   yn         = y-coordinate of the vector normal to the curve, but
c                tangent to the sphere
c   zn         = z-coordinate of the vector normal to the curve, but
c                tangent to the sphere
c   dsda       = arclength term
c   arcLength  = ds term of each component curve
c   diag       = diagonal term which depends on the curvature of the
c                ellipse
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension ak(k), bk(k), th_k(k), phi_k(k), r(3), t(3), pn(3),
     1          vn(3), diag(nbk), d2x(nbk), d2y(nbk), d2z(nbk)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          xn(nbk), yn(nbk), zn(nbk), dsda(nbk), arcLength(k)
c
      pi = 4.d0*datan(1.d0)
c
      dalph = 2.d0*pi/nd
      istart = 0
      do kbod = 1, k
        do i = 1, nd
          alpha = dalph*(i-1.d0)
c Construct points on the sphere
          call r_func (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                 phi_k(kbod), xs(istart+i), ys(istart+i), 
     2                 zs(istart+i))
c Construct derivatives of points on the sphere
          call dr_func (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                  phi_k(kbod), dx(istart+i), dy(istart+i), 
     2                  dz(istart+i))
c
c Construct normal to surface of sphere
          r(1) = xs(istart+i)
          r(2) = ys(istart+i)
          r(3) = zs(istart+i)
c
c Construct tangent to curve lying in plane of sphere
          t(1) = dx(istart+i)
          t(2) = dy(istart+i)
          t(3) = dz(istart+i)
          dsda(istart+i) = dsqrt((t(1))**2 + (t(2))**2 + (t(3))**2)
          do j = 1, 3
            t(j) = t(j)/dsda(istart+i)
          end do
c
c Construct normal to curve lying in plane of sphere
          call cross (t,r,vn)
          xn(istart+i) = vn(1)
          yn(istart+i) = vn(2)
          zn(istart+i) = vn(3)
c
c Construct the diagonal entry
          call d2r_func (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                    phi_k(kbod), pn(1), pn(2), pn(3))
          d2x(istart+i) = pn(1)
          d2y(istart+i) = pn(2)
          d2z(istart+i) = pn(3)
          pn(1) = pn(1)/(dsda(istart+i))**2
          pn(2) = pn(2)/(dsda(istart+i))**2
          pn(3) = pn(3)/(dsda(istart+i))**2
          call cross(pn,r,vn)
          call dot (t,vn,t_dot_vn)
          diag(istart+i) = -t_dot_vn*dsda(istart+i)/(4.d0*pi)
        end do
        istart = istart + nd
      end do

      istart = 0
      do kbod = 1, k
        arcLength(kbod) = 0.d0
        do i = 1, nd
          arcLength(kbod) = arcLength(kbod) + dsda(istart + i)*dalph
        enddo
        arcLength(kbod) = arcLength(kbod)/dble(nd)
        istart = istart + nd
      enddo
c find arclength of each componenet curve for doing near-singular
c integration


      is = 1
c make a Matlab file that can plot the geometry
      open (unit = 42, file = 'geo_3d.m')
      do kbod = 1, k
        call RS_3D_PLOT (xs(is),ys(is),zs(is), nd, 1, 42)
        is = is + nd
      end do
      close (42)
c

      return
      end ! MAKE_GEO 
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine findZone(x, y, z, ak, bk, th_k, phi_k, 
     1      nd, arcLength, izone)
c
c  *** DESCRIPTION :
c
c   Classify a point as either being far, intermediate, or near to
c   the ellipse given by the parameters ak, bk, th_k, phi_k
c
c   *** INPUT PARAMETERS :
c
c   x         = x-coordinate of point on sphere
c   y         = y-coordinate of point on sphere
c   z         = z-coordinate of point on sphere
c   ak        = parameter for ellipses on sphere
c   bk        = parameter for ellipses on sphere
c   th_k      = parameter for ellipses on sphere
c   phi_k     = parameter for ellipses on sphere
c   nd        = number of points on curve
c   arcLength = arclength term of curve
c
c   *** OUTPUT PARAMETERS :
c
c   izone     = flag which is 2 (far zone), 1 (intermediate zone),
c               or 0 (near zone)
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      twopi = 8.d0*datan(1.d0)
      dalph = twopi/dble(nd)

      dist2 = 1.d10

      do i = 1,nd
        alpha = dble(i-1)*dalph
        call r_func(alpha, ak, bk, th_k, phi_k, xs, ys, zs)
        dist2 = min(dist2,(xs-x)**2.d0 + (ys-y)**2.d0 + (zs-z)**2.d0)
      enddo

      if (dist2 > arcLength) then
c       in far zone
        izone = 2
      elseif (sqrt(dist2) > arcLength) then
c       in intermediate zone
        izone = 1
      else
c       in near zone
        izone = 0
      endif

      return
      end ! findZone
c      
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine IN_OR_OUT (x, y, z, eps, ak, bk, th_k, phi_k, itest)
c
c  *** DESCRIPTION :
c
c   Determine if a point is inside or outside of an ellipse.  Future
c   code should use the trick of computing the double-layer potential
c   of Laplace with a constant density function which is 1 on the
c   inside and 0 on the outside of the ellipses
c
c   *** INPUT PARAMETERS :
c
c   x       = x-coordinate of point on sphere
c   y       = y-coordinate of point on sphere
c   z       = z-coordinate of point on sphere
c   eps     = threshold used to account for numerical error
c   ak      = parameter for ellipses on sphere
c   bk      = parameter for ellipses on sphere
c   th_k    = parameter for ellipses on sphere
c   phi_k   = parameter for ellipses on sphere
c
c   *** OUTPUT PARAMETERS :
c
c   itest   = flag that determines if a point is inside or outside
c             a specific ellipse on the sphere
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension p(3), x_ax(3), y_ax(3), z_ax(3)
c
      pi = 4.d0*datan(1.d0)
c
      x2 = dsin(phi_k)*dcos(th_k)*x + 
     1      dsin(phi_k)*dsin(th_k)*y - dcos(phi_k)*z
      y2 = -dsin(th_k) * x + dcos(th_k) * y
      z2 = cos(phi_k)*cos(th_k)*x + 
     1      dcos(phi_k)*dsin(th_k)*y + dsin(phi_k)*z
      rad = dsqrt((x2/ak)**2.d0 + (y2/bk)**2.d0)

      if (rad<=1.0d0 .and. (z2 > 0.d0)) then 
c point is inside (or almost) the curve (outside the geometry)
        itest = 1
      elseif (rad<=1.0d0 + eps .and. (z2 > 0.d0)) then 
        itest = -1
c point is outside the curve (inside the geometry), but near
      else
        itest = 0
c point is outside the curve (inside the geometry) and well-separated
      end if


      return
      end ! IN_OR_OUT     
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SURFACE_GRID (k, nd, nbk, freq, nth, nphi, ak, bk, 
     1                   cx, cy, cz, th_k, phi_k, th_gr, phi_gr, 
     2                   arcLength,
     3                   nvort, x1Vort, x2Vort, x3Vort, vortK, 
     4                   x_gr, y_gr, z_gr, igrid, near_gr, uExact_gr)
c
c  *** DESCRIPTION :
c
c   Build a meshgrid structure in spherical coordinates and determine
c   which points are inside and which are outside of the geometry.
c   Write all of this information to dat files
c
c   *** INPUT PARAMETERS :
c
c   k         = number of connected componenets in geometry
c   nd        = number of points per connected component
c   nbk       = total number of points k*nd
c   freq      = constant in the pde (\Delta - freq^2)u = 0
c   nth       = number of points in theta direction for meshgrid
c   nphi      = number of points in phi direction for meshgrid
c   ak        = parameter for ellipses on sphere
c   bk        = parameter for ellipses on sphere
c   cx        = x-coordinate of center ellipses on sphere
c   cy        = y-coordinate of center ellipses on sphere
c   cz        = z-coordinate of center ellipses on sphere
c   th_k      = parameter for ellipses on sphere
c   phi_k     = parameter for ellipses on sphere
c   arcLength = arclength term of each component curve
c   nvort     = number of vorticies to form exact solution
c   x1Vort    = x-coordinate of the vorticies
c   x2Vort    = y-coordinate of the vorticies
c   x3Vort    = z-coordinate of the vorticies
c
c   *** OUTPUT PARAMETERS :
c
c   th_gr     = grid of azimuthal coordinates
c   phi_gr    = grid of z-uthal coordinates
c   x_gr      = x-coordinate of grid on sphere
c   y_gr      = y-coordinate of grid on sphere
c   z_gr      = z-coordinate of grid on sphere
c   igrid     = flag for if a point is inside or outside of the
c               geometry of interest
c   near_gr   = flag for each hole that determines if near-singular
c               integration is required or not
c   uExact_gr = exact solution formed using the fundamental solution
c               centered outside of the boundary components
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension ak(k), bk(k), cx(k), cy(k), cz(k), th_k(k), phi_k(k)
      dimension arcLength(k)
      dimension x1Vort(nvort), x2Vort(nvort), x3Vort(nvort)
      dimension vortK(nvort)
      dimension igrid(nth,nphi), near_gr(nth,nphi,k)
      dimension th_gr(nth,nphi), phi_gr(nth,nphi),
     1          x_gr(nth,nphi), y_gr(nth,nphi), z_gr(nth,nphi),
     2          uExact_gr(nth,nphi)
c
      pi = 4.d0*datan(1.d0)
c


      dth = 2.d0*pi/nth
      dphi = pi/(nphi-1)
      do i = 1, nth
        theta = (i-1)*dth
        do j = 1, nphi
          igrid(i,j) = 1
          phi = dble(j-1)*dphi - 0.5d0*pi
          th_gr(i,j) = theta
          phi_gr(i,j) = phi
c Cartesian coordinate of meshgrid on the sphere
          call sph2cart (th_gr(i,j), phi_gr(i,j), 1.d0, 
     1                   x_gr(i,j), y_gr(i,j), z_gr(i,j))
          in_out = 0
          do kbod = 1, k
c Check if it is inside any of the component curves
            eps = 2.d1 * 2.d0*pi*max(ak(kbod),bk(kbod))/nd
            eps = 1.0d10
c Hack for now to make sure we don't miss any intermediate points
            call IN_OR_OUT (x_gr(i,j), y_gr(i,j), z_gr(i,j), eps,
     1                       ak(kbod), bk(kbod), th_k(kbod), 
     2                       phi_k(kbod), itest)
c Check if point is possibly in the near or intermediate region
c Also, rule out the points outside the geometry here
            in_out = in_out + itest
            if (itest .eq. -1) then
c             point is potentially in near zone
              call findZone(x_gr(i,j), y_gr(i,j), z_gr(i,j), 
     1          ak(kbod),bk(kbod), th_k(kbod), phi_k(kbod), 
     2          nd, arcLength(kbod),near_gr(i,j,kbod))
            elseif (itest .eq. 0) then
c             point is definitely in far zone
              near_gr(i,j,kbod) = 2
            else
c             point is definitely outside geometry
              near_gr(i,j,kbod) = -1
            endif

            if (near_gr(i,j,kbod) .eq. -1) then
              igrid(i,j) = 0
            endif

          end do
        end do
      end do 

c Form the exact solution at all points that are inside the
c geometry
      do i = 1,nphi
        do j = 1, nphi
          uExact_gr(i,j) = 0.d0    
          if (igrid(i,j) .eq. 1) then 
            do ivort = 1, nvort
              dist2 = (x_gr(i,j) - x1Vort(ivort))**2.d0 + 
     1                (y_gr(i,j) - x2Vort(ivort))**2.d0 + 
     2                (z_gr(i,j) - x3Vort(ivort))**2.d0 
              uExact_gr(i,j) = uExact_gr(i,j) + 
     1            vortK(ivort)*yukawaSLP(freq,dist2)
            enddo
          endif
        enddo
      enddo


c     Set exact solution to be zero at the near points
c     Accurate evaulation at these points is not yet complete
      do i = 1,nphi
        do j = 1,nphi
          do kbod = 1,k
            if (near_gr(i,j,kbod) .eq. 0) then
              uExact_gr(i,j) = 1.d-5
            endif
          enddo
        enddo
      enddo



c Write the location of points in the geomtery and 
c meshgrid to dat files
      open (unit = 32, file = 'xgrid.dat')
      open (unit = 33, file = 'ygrid.dat')
      open (unit = 34, file = 'zgrid.dat')
      call dump (nth, nphi, x_gr, igrid, 1, 32)
      call dump (nth, nphi, y_gr, igrid, 1, 33)
      call dump (nth, nphi, z_gr, igrid, 1, 34)
      close (32)
      close (33)
      close (34)
c
      return ! SURFACE_GRID
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine getRHS (k, nd, nbk, nvort, freq, xs, ys, zs, 
     1    x1Vort, x2Vort, x3Vort, vortK, rhs)
c
c  *** DESCRIPTION :
c
c   Get the right-hand side which is the Dirichlet boundary condition
c
c   *** INPUT PARAMETERS :
c
c   k         = number of connected componenets in geometry
c   nd        = number of points per connected component
c   nbk       = total number of points k*nd
c   nvort     = number of vorticies to form exact solution
c   freq      = constant in the pde (\Delta - freq^2)u = 0
c   xs        = x-coordinate of the ellipses
c   ys        = y-coordinate of the ellipses
c   zs        = z-coordinate of the ellipses
c   x1Vort    = x-coordinate of the vorticies
c   x2Vort    = y-coordinate of the vorticies
c   x3Vort    = z-coordinate of the vorticies
c   vortK     = strength of the vorticies
c
c   *** OUTPUT PARAMETERS :
c
c   rhs       = Dirichlet boundary condition
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension xs(k*nd), ys(k*nd), zs(k*nd)
      dimension x1Vort(nvort), x2Vort(nvort), x3Vort(nvort)
      dimension vortK(nvort)
      dimension rhs(nbk)
c
      istart = 0
      do kbod = 1, k
        do j = 1, nd
          if (nvort .eq. 0) then
            if (kbod .le. 6) then
              rhs(istart+j) = 1.d0
            elseif (kbod .le. 12) then
              rhs(istart+j) = 1.d0
            elseif (kbod .le. 24) then
              rhs(istart+j) = 1.d0
            elseif (kbod .le. 30) then
              rhs(istart+j) = 1.d0
            else
              rhs(istart+j) = 1.d0
            endif
          else
            rhs(istart+j) = 0.d0
          endif
        end do
        istart = istart+nd
      end do
c     if there are no vorticies, then take a simple boundary
c     condition.  Otherwise, initialize the boundary condition
c     to zero

      istart = 0
      do kbod = 1,k
        do j = 1,nd
          do ivort = 1,nvort
            dist2 = (xs(istart+j) - x1Vort(ivort))**2.d0 + 
     1              (ys(istart+j) - x2Vort(ivort))**2.d0 + 
     2              (zs(istart+j) - x3Vort(ivort))**2.d0 
            rhs(istart+j) = rhs(istart+j) + vortK(ivort)*
     1            yukawaSLP(freq,dist2)
          enddo
        enddo
        istart = istart + nd
      enddo
c     Sum up the fundamental solution centered at various points
c     outside of the geometry.


c     Dump Dirichlet boundary condition to dat file
      open (unit = 24, file = 'rhs.dat')
      do i = 1, nbk
        write(24,'(e20.13,$)')(rhs(i))
        write (24,'(a)')  ''
      end do
      close (24) 
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine solveBIE(nd, k, nbk, rhs, density, rwork, 
     1                  lrwork, iwork, liwork, dsda, maxl)
c
c  *** DESCRIPTION :
c
c   Use GMRES to solve the second-kind integral equation
c
c   *** INPUT PARAMETERS :
c
c   nd        = number of points per connected component
c   k         = number of connected componenets in geometry
c   nbk       = total number of points k*nd
c   rhs       = Dirichlet boundary condition
c   rwork     = real workspace for gmres
c   lrwork    = length of real workspace for gmres
c   iwork     = integer workspace for gmres
c   liwork    = length of integer workspace for gmres
c   dsda      = arclength term
c   maxl      = maximum number of GMRES iterations before restart
c
c   *** OUTPUT PARAMETERS :
c
c   density   = density function
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      external matvecYukawa, msolve
c
c  System
      dimension rhs(nbk), density(nbk), dsda(nbk)
c
c  DGMRES work arrays
      dimension rwork(lrwork),iwork(liwork)
c
c  Timings
      real*4 timep(2), etime
c
c     parameters for DGMRES
      itol = 0
      tol = 1.0d-11
      isym = 0
      iwork(1) = maxl
      do i=2,liwork
        iwork(i) = 0
      enddo
c
c  Preconditioner flag
      iwork(4) = 0
c
c  Restart flag
      iwork(5) = -1
c
c  provide initial guess density
      do i=1,nbk
        density(i) = rhs(i)
      enddo
c
      t0 = etime(timep)
c  Solve linear system with GMRES
      call DGMRES (nbk, rhs, density, nelt, ia, ja, a, isym,
     1            matvecYukawa, msolve, itol, tol, itmax, iter, err,  
     1            ierr, 6, sb, sx, rwork, lrwork, iwork, 
     1            liwork, rw, iw)
      call PRIN2 (' after yukawa solve, err = *', err, 1)
      call PRINF ('  # GMRES ITERATIONS = *',iter,1)
      if (ierr.gt.2) then
        call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
        call PRINF ('  iwork = *',iwork,10)
        stop
      elseif (ierr.ge.0) then
        t1 = etime(timep)
        tsec = t1 - t0
        call PRIN2 (' time in solve = *', tsec, 1)
      end if

c Dump it out
      open (unit = 24, file = 'solution.dat')
      do i = 1, nbk
        write(24,'(e20.13,$)')(density(i))
        write (24,'(a)')  ''
      end do
      close (24) 
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine matvecYukawa (n, xx, yy, nelt, ia, ja, a, isym)
c
c  *** DESCRIPTION :
c
c  Required by DGMRES with this precise calling sequence.
c  We ignore most of the parameters except n, xx and yy
c  and call the fast multipole routine FASMVP to do the actual
c  work.
c
c   *** INPUT PARAMETERS :
c
c   n        = total number of points in linear system
c   xx       = input xx in yy = A*xx
c   nelt     = number of elements in matrix A (this doesn't matter to
c              us)
c   ia       = data structure for A (this doesn't matter to us)
c   ja       = data structure for A (this doesn't matter to us)
c   a        = matrix A in yy = A*xx.  In general, method is matrix
c              free, so this is never constructed
c   isym     = flag indicating if the problem is symmetric
c
c   *** OUTPUT PARAMETERS :
c
c   yy       = output yy in yy = A*xx
c
c***********************************************************************
c
      implicit real *8 (a-h,o-z)
c
      dimension xx(n), yy(n)
      parameter (kmax = 50, npmax = 512, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
c
c common blocks
      common /sys_size/ k, nd, nbk
      common /sphere_int/ xs, ys, zs, xn, yn, zn, dsda, diag, freq
c
c Geometry of holes
      dimension xs(nmax), ys(nmax), zs(nmax)
      dimension xn(nmax), yn(nmax), zn(nmax)
      dimension dsda(nmax), diag(nmax) 
c
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/nd

c     flag for deciding if we do Alpert corrections or not
      ialpert = 0
      call AlpertQuadrature(k,nd,xs,ys,zs,xn,yn,zn,
     1        dsda,freq,xx,ialpert,yy)

c     loop over target points
      do i = 1, nbk
c       loop over source points
c       START OF TRAPEZOID RULE
        do j = 1, nbk
          if (i .ne. j) then
c  distance squared between source and target
            dist2 = (xs(i) - xs(j))**2.d0 + 
     1              (ys(i) - ys(j))**2.d0 + 
     2              (zs(i) - zs(j))**2.d0 
c  dot product of difference vector with normal
            rdotn = (xs(i) - xs(j))*xn(j) + 
     1              (ys(i) - ys(j))*yn(j) + 
     2              (zs(i) - zs(j))*zn(j)

c  argument required by hypergeometric representation
            yy(i) = yy(i) + yukawaDLP(freq,dist2,rdotn)*
     1          xx(j)*dsda(j)*dalph
c           update solution
          endif
        end do !j
c       END OF TRAPEZOID RULE
        
c       jump term
        yy(i) = yy(i) + 5.d-1*xx(i)
        if (ialpert .eq. 0) then
c         diagonal term if not using Alpert's quadrature
          yy(i) = yy(i) + diag(i)*dalph*xx(i)
        endif
      end do !i
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      real *8 function yukawaDLP(freq,dist2,rdotn)
c
c  *** DESCRIPTION :
c
c   Compute the term in the integrand of the double-layer potential
c   This does not include the arclength or density function terms
c
c   *** INPUT PARAMETERS 
c
c   freq       = constant in the pde (\Delta - freq^2)u = 0
c   dist2      = distance squared between the source and target
c   rdotn      = dot product of the difference between the source
c                and target with the normal
c
c   *** OUTPUT PARAMETERS :
c
c   yukawaDLP  = the double-layer potential kernel
c
c***********************************************************************
c
      use someconstants
      use hyp_2f1_module
c
      implicit real*8 (a-h,o-z)
c
      complex *16 nu,a,b,c,hyp_2f1,cdlp,eye
c
      eye = (0.d0,1.d0)

c     parameter for using hypergeometric functions
      tau = dsqrt(4*freq**2.d0 - 1.d0)/2.d0
c     parameter if using associated Legendre P functions
      nu = -0.5d0 + eye*tau
c     parameters if using hypergeometric function 2F1
      a = nu + 1.d0
      b = -1.d0*nu;
      c = (1.d0,0.d0)
c     scaling that multiplies the fundamental solution so
c     that assymptotically it looks like 1/(2*Pi)*log||x-x_{0}||
      fundConst = 0.25d0/dcosh(pi*tau)

c     argument required by hypergeometric representation
      z = 1.d0 - 2.5d-1*dist2
      cdlp = (-1.d0 + 5.d-1*dist2) * hyp_2f1(a,b,c,dcmplx(z))
      cdlp = cdlp - hyp_2f1(a+1,b-1,c,dcmplx(z))
      cdlp = a * cdlp/((-1.d0 + 5.d-1*dist2)**2.d0 - 1.d0)
      cdlp = cdlp * rdotn
c     cdlp is the double-layer potential.  If it has an imaginary
c     component, there is a bug

c     multiply by the constant so that the jump is 1
c     across the boundary
      yukawaDLP = fundConst*dreal(cdlp)


      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      real *8 function yukawaSLP(freq,dist2)
c
c  *** DESCRIPTION :
c
c   Compute the term in the integrand of the single-layer potential
c   This does not include the arclength or density function terms
c
c   *** INPUT PARAMETERS 
c
c   freq       = constant in the pde (\Delta - freq^2)u = 0
c   dist2      = distance squared between the source and target
c   rdotn      = dot product of the difference between the source
c                and target with the normal
c
c   *** OUTPUT PARAMETERS :
c
c   yukawaSLP  = the single-layer potential kernel
c
c***********************************************************************
c
      use someconstants
      use hyp_2f1_module
c
      implicit real*8 (a-h,o-z)

      complex *16 nu,a,b,c,hyp_2f1,cdlp,eye
c
      eye = (0.d0,1.d0)

c     parameter for using hypergeometric functions
      tau = dsqrt(4*freq**2.d0 - 1.d0)/2.d0
c     parameter if using associated Legendre P functions
      nu = -0.5d0 + eye*tau
c     parameters if using hypergeometric function 2F1
      a = nu + 1.d0
      b = -1.d0*nu;
      c = (1.d0,0.d0)
c     scaling that multiplies the fundamental solution so
c     that assymptotically it looks like 1/(2*Pi)*log||x-x_{0}||
      fundConst = 0.25d0/dcosh(pi*tau)

c     argument required by hypergeometric representation
      z = 1.d0 - 2.5d-1*dist2
      yukawaSLP = fundConst*dreal(hyp_2f1(a,b,c,dcmplx(z)))


      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine AlpertQuadrature(nhole,nd,xs,ys,zs,xn,yn,zn,
     1        dsda,freq,xx,ialpert,yy)
c
c  *** DESCRIPTION :
c
c   Do modification to trapeziod rule that corresponds to Alpert's
c   quadrature rule.  If ialpert=0, this routine simple sets yy to
c   zero and then exits
c
c   *** INPUT PARAMETERS :
c
c   nhole   = number of connected componenets in geometry
c   nd      = number of points per connected component
c   xs      = x-coordinate of the ellipses
c   ys      = y-coordinate of the ellipses
c   zs      = z-coordinate of the ellipses
c   xn      = x-coordinate of the vector normal to the curve, but
c             tangent to the sphere
c   yn      = y-coordinate of the vector normal to the curve, but
c             tangent to the sphere
c   zn      = z-coordinate of the vector normal to the curve, but
c             tangent to the sphere
c   dsda    = arclength term
c   freq    = constant in the pde (\Delta - freq^2)u = 0
c   xx      = density function
c   ialpert = flag to determine if Alpert's rules are used or not
c
c   *** OUTPUT PARAMETERS :
c
c   yy      = layer potential due to the Alpert quadrature rules
c             This includes subtracting the terms around the singularity
c             which are added in by doing the usual trapezoid rule
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension xs(nhole*nd), ys(nhole*nd), zs(nhole*nd)
      dimension xn(nhole*nd), yn(nhole*nd), zn(nhole*nd)
      dimension dsda(nhole*nd)
      dimension xx(nhole*nd),yy(nhole*nd)

c Alpert quadrature rule nodes and weights
      dimension u(15),v(15),vv(30)

      complex *16, allocatable :: zfn(:),zInterp(:,:)
c     temporary space to store the functions that need
c     Fourier interpolation
      complex *16, allocatable :: wsave(:)
c     workspace for FFT and IFFT
      real *8, allocatable :: xsShift(:,:)
      real *8, allocatable :: ysShift(:,:)
      real *8, allocatable :: zsShift(:,:)
      real *8, allocatable :: xnShift(:,:)
      real *8, allocatable :: ynShift(:,:)
      real *8, allocatable :: znShift(:,:)
      real *8, allocatable :: dsdaShift(:,:)
      real *8, allocatable :: xxShift(:,:)
c     shifted versions of the different periodic variables
c
      complex *16 eye
c
      eye = (0.d0,1.d0)
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/dble(nd)

      do i = 1,nhole*nd
        yy(i) = 0.d0
      enddo
c     Initialize to zero

      if (ialpert .eq. 0) return
c     If we are not using Alpert, this the subroutine ends here

      call quad2(v,u,numquad,nbuffer,norder,6)
c     Get quadrature nodes, weights, and region around singularity
c     that is excluded

      allocate(zfn(nd),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(zInterp(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(wsave(4*nd+15),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(xsShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(ysShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(zsShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(xnShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(ynShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(znShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dsdaShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(xxShift(nd,2*numquad),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
c     Allocate memory to store all required functions at
c     the quadrature points using Fourier interpolation

      do i = 1,numquad
        vv(i) = v(i)
        vv(i+numquad) = -v(i)
      enddo
c     want to compute the forward and backward shift all
c     at once

      call DCFFTI(nd,wsave)
c     initialize work array wsave for doing FFT and IFFT

      istart = 0
      do ihole = 1,nhole
c       START OF FORMING PERIODIC FUNCTIONS ON THE SHIFTED GRID
c       FOR USING ALPERT QUADRATURE RULES
c       shift the first two components of the normal vector
        do i = 1,nd
          zfn(i) = xn(i + istart) + eye*yn(i + istart)
        enddo
        call fourierInterp(nd,zfn,2*numquad,vv,wsave,zInterp)
        do k = 1,2*numquad
          do i = 1,nd
            xnShift(i,k) = dreal(zInterp(i,k))
            ynShift(i,k) = dimag(zInterp(i,k))
          enddo
         enddo

c       shift the first two components of the position vector
        do i = 1,nd
          zfn(i) = xs(i + istart) + eye*ys(i + istart)
        enddo
        call fourierInterp(nd,zfn,2*numquad,vv,wsave,zInterp)
        do k = 1,2*numquad
          do i = 1,nd
            xsShift(i,k) = dreal(zInterp(i,k))
            ysShift(i,k) = dimag(zInterp(i,k))
          enddo
        enddo

c       shift the third component of the position and 
c       normal vector
        do i = 1,nd
          zfn(i) = zs(i + istart) + eye*zn(i + istart)
        enddo
        call fourierInterp(nd,zfn,2*numquad,vv,wsave,zInterp)
        do k = 1,2*numquad
          do i = 1,nd
            zsShift(i,k) = dreal(zInterp(i,k))
            znShift(i,k) = dimag(zInterp(i,k))
          enddo
        enddo

c       shift the arclength component and the density function
        do i = 1,nd
          zfn(i) = dsda(i + istart) + eye*xx(i + istart)
        enddo
        call fourierInterp(nd,zfn,2*numquad,vv,wsave,zInterp)
        do k = 1,2*numquad
          do i = 1,nd
            dsdaShift(i,k) = dreal(zInterp(i,k))
            xxShift(i,k) = dimag(zInterp(i,k))
          enddo
        enddo
c       END OF FORMING PERIODIC FUNCTIONS ON THE SHIFTED GRID
c       FOR USING ALPERT QUADRATURE RULES

c       START OF FIRST AND THIRD TERM IN ALPERT'S CORRECTION
        do i = 1,nd
          do k = 1,numquad
c  distance squared between source and Alpert quadrature point
c  forward from the target point
            dist2F = (xs(i+istart) - xsShift(i,k))**2.d0 + 
     1               (ys(i+istart) - ysShift(i,k))**2.d0 + 
     2               (zs(i+istart) - zsShift(i,k))**2.d0
c  distance squared between source and Alpert quadrature point
c  backward from the target point
            dist2B = (xs(i+istart) - xsShift(i,k+numquad))**2.d0 + 
     1               (ys(i+istart) - ysShift(i,k+numquad))**2.d0 + 
     2               (zs(i+istart) - zsShift(i,k+numquad))**2.d0
c  dot product of difference vector with normal at Alpert quadrature
c  point forward from the target point
            rdotnF = (xs(i+istart) - xsShift(i,k))*xnShift(i,k) + 
     1               (ys(i+istart) - ysShift(i,k))*ynShift(i,k) + 
     2               (zs(i+istart) - zsShift(i,k))*znShift(i,k)
c  dot product of difference vector with normal at Alpert quadrature
c  point backward from the target point
            rdotnB = (xs(i+istart) - xsShift(i,k+numquad))*
     1               xnShift(i,k+numquad) + 
     2               (ys(i+istart) - ysShift(i,k+numquad))*
     3               ynShift(i,k+numquad) + 
     4               (zs(i+istart) - zsShift(i,k+numquad))*
     5               znShift(i,k+numquad)

c  argument required by hypergeometric representation
            yy(i+istart) = yy(i+istart) + 
     1          u(k)*yukawaDLP(freq,dist2F,rdotnF)*
     2          xxShift(i,k)*dsdaShift(i,k)*dalph
c           first term update
            yy(i+istart) = yy(i+istart) + 
     1          u(k)*yukawaDLP(freq,dist2B,rdotnB)*
     2          xxShift(i,k+numquad)*dsdaShift(i,k+numquad)*dalph
c           third term update

          enddo
        enddo
c       END OF FIRST AND THIRD TERM IN ALPERT'S CORRECTION

c       START OF SECOND TERM IN ALPERT'S CORRECTION
        do i = 1,nd
          midSum = nd - 2*nbuffer + 1
          do k = midSum,nd-1
            index = mod(i+nbuffer+k,nd)
            if (index .eq. 0) index = nd
            if (index .ne. i) then
c  distance squared between source and target point
              dist2 = (xs(i+istart) - xs(index+istart))**2.d0 + 
     1                (ys(i+istart) - ys(index+istart))**2.d0 + 
     2                (zs(i+istart) - zs(index+istart))**2.d0
c  dot product of difference vector with normal vector 
              rdotn = (xs(i+istart) - xs(index+istart))*
     1                xn(index+istart) + 
     2                (ys(i+istart) - ys(index+istart))*
     3                yn(index+istart) + 
     4                (zs(i+istart) - zs(index+istart))*
     5                zn(index+istart)

              yy(i+istart) = yy(i+istart) - yukawaDLP(freq,dist2,rdotn)*
     1              xx(index+istart)*dsda(index+istart)*dalph 
            endif
          enddo
        enddo
c       END OF SECOND TERM IN ALPERT'S CORRECTION



        istart = istart + nd
      enddo !ihole

      deallocate(zfn)
      deallocate(zInterp)
      deallocate(wsave)
      deallocate(xsShift)
      deallocate(ysShift)
      deallocate(zsShift)
      deallocate(xnShift)
      deallocate(ynShift)
      deallocate(znShift)
      deallocate(dsdaShift)
      deallocate(xxShift)

      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine fourierUpsample(nIn,wsaveIn,zIn,nOut,wsaveOut,zOut)
c
c  *** DESCRIPTION :
c
c  Upsample a complex valued periodic function by padding the spectrum
c  with zeros
c
c   *** INPUT PARAMETERS :
c
c   nIn       = number of input points
c   wsaveIn   = precomputed work array for doing fft at input grid
c   zIn       = periodic function that needs to be interpolated
c   nOut      = number of output points
c   wsaveOut  = precomputed work array for doing ifft at output grid
c
c   *** OUTPUT PARAMETERS :
c
c   zOut      = output interpolated periodic function
c
c***********************************************************************
c
      implicit real *8 (a-h,o-z)
      complex *16 zIn(nIn)
      complex *16 wsaveIn(*)
      complex *16 wsaveOut(*)

      complex *16 zOut(nOut)

      do i = 1,nOut
        zOut(i) = (0.d0,0.d0)
      enddo

      call DCFFTF(nIn,zIn,wsaveIn)

      do i = 1,nIn/2
        zOut(i) = zIn(i)/dble(nIn)
        zOut(nOut - i + 1) = zIn(nIn - i + 1)/dble(nIn)
      enddo

      call DCFFTB(nOut,zOut,wsaveOut)

      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine fourierInterp(nd,zfn,nshifts,shift,wsave,zInterp)
c
c  *** DESCRIPTION :
c
c   Translate a complex valued function zfn by amounts in v so that
c   Alpert's quadrature can be used
c
c   *** INPUT PARAMETERS :
c
c   nd        = number of points per connected component
c   zfn       = periodic function that needs to be shifted
c   nshifts   = number of periodic shifts to be done
c   shift     = amounts periodic functions are to be shifted
c   wsave     = precomputed work array for doing fft and ifft
c
c   *** OUTPUT PARAMETERS :
c
c   zInterp   = shifted periodic function 
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      complex *16 zfn(nd)
      real *8 shift(nshifts) 
      complex *16 wsave(*)

      complex *16 zInterp(nd,nshifts)

      complex *16, allocatable :: alpha(:)
      complex *16 eye

      allocate (alpha(nd),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif

      twopi = 8.d0*datan(1.0d0)
      dx = twopi/dble(nd)
      eye = (0.d0,1.0d0)

      call DCFFTF(nd,zfn,wsave)
c     Fourier series of the density function

      do k=1,nshifts
        do j=1,nd
          if (j .le. nd/2) then
            alpha(j) = zfn(j)*exp(dble(j-1)*eye*shift(k)*dx)
          else
            alpha(j) = zfn(j)*exp(dble(j-nd-1)*eye*shift(k)*dx)
          endif
        enddo
        call DCFFTB(nd,alpha,wsave)
        do j = 1,nd
          zInterp(j,k) = alpha(j)/dble(nd)
        enddo
c       Compute xx on the grid centered forwards by distances        
      enddo

      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SOL_GRID (freq, nd, nhole, nbk, nth, nphi, density, 
     1        xs, ys, zs, xn, yn, zn, dsda, arcLength, 
     2        x_gr, y_gr, z_gr, igrid, near_gr, u_gr)
c
c  *** DESCRIPTION :
c
c   Evaluate the double-layer potential representation at 
c   a collection of points that are inside the geometry of interest.
c   Have to extend this to the double-layer potential
c
c   *** INPUT PARAMETERS :
c
c   freq      = constant in the pde (\Delta - freq^2)u = 0
c   nd        = number of points per connected component
c   k         = number of connected componenets in geometry
c   nbk       = total number of points k*nd
c   nth       = number of points in the azimuthal direction
c   nphi      = number of points in the z-uthal direction
c   density   = density function defined on the boundary of the ellipses
c   xs        = x-coordinate of the ellipses
c   ys        = y-coordinate of the ellipses
c   zs        = z-coordinate of the ellipses
c   xn        = x-coordinate of the vector normal to the curve, but
c               tangent to the sphere
c   yn        = y-coordinate of the vector normal to the curve, but
c               tangent to the sphere
c   zn        = z-coordinate of the vector normal to the curve, but
c               tangent to the sphere
c   dsda      = arclength term
c   arcLength = infintesimal arclength term
c   x_gr      = x-coordinate of grid on sphere
c   y_gr      = y-coordinate of grid on sphere
c   z_gr      = z-coordinate of grid on sphere
c   igrid     = flag to determine if a point is inside or outside
c   near_gr   = flag to determine if point is near, intermediate or far
c
c   *** OUTPUT PARAMETERS :
c
c   u_gr      = solution of the PDE at points in the geometry
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension xs(nbk), ys(nbk), zs(nbk)
      dimension xn(nbk), yn(nbk), zn(nbk)
      dimension density(nbk), dsda(nbk)
      dimension arcLength(nhole)
      dimension x_gr(nth,nphi), y_gr(nth,nphi), z_gr(nth,nphi)
      dimension igrid(nth,nphi), near_gr(nth,nphi,nhole)
      dimension u_gr(nth,nphi)
      real*4 timep(2), etime
      complex *16 eye
c
      complex *16, allocatable :: wsaveIn(:),wsaveOut(:)
      complex *16, allocatable :: zIn(:),zOut(:)
      real *8, allocatable :: xsUp(:),ysUp(:),zsUp(:)
      real *8, allocatable :: xnUp(:),ynUp(:),znUp(:)
      real *8, allocatable :: densityUp(:), dsdaUp(:)
c
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/nd
      eye = (0.d0,1.0d0)

      base2 = dlog(dble(nd))/dlog(2.d0)
      nup = 2**ceiling(base2/2.0d0)
c     upscaling factor

      allocate(zIn(nd),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(zOut(nup*nd),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(wsaveIn(4*nd+15),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(wsaveOut(4*nup*nd+15),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(xsUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(ysUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(zsUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(xnUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(ynUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(znUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dsdaUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(densityUp(nup*nbk),stat=ierr)
      if (ierr .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif

      call DCFFTI(nd,wsaveIn)
      call DCFFTI(nd*nup,wsaveOut)

      nhole = nbk/nd
      istartIn = 0
      istartOut = 0
      do ihole = 1,nhole

        do k = 1,nd
          zIn(k) = xs(istartIn + k) + eye*xn(istartIn + k)
        enddo
        call fourierUpsample(nd,wsaveIn,zIn,nd*nup,wsaveOut,zOut)
        do k = 1,nup*nd
          xsUp(istartOut + k) = dreal(zOut(k))
          xnUp(istartOut + k) = dimag(zOut(k))
        enddo

        do k = 1,nd
          zIn(k) = ys(istartIn + k) + eye*yn(istartIn + k)
        enddo
        call fourierUpsample(nd,wsaveIn,zIn,nd*nup,wsaveOut,zOut)
        do k = 1,nup*nd
          ysUp(istartOut + k) = dreal(zOut(k))
          ynUp(istartOut + k) = dimag(zOut(k))
        enddo

        do k = 1,nd
          zIn(k) = zs(istartIn + k) + eye*zn(istartIn + k)
        enddo
        call fourierUpsample(nd,wsaveIn,zIn,nd*nup,wsaveOut,zOut)
        do k = 1,nup*nd
          zsUp(istartOut + k) = dreal(zOut(k))
          znUp(istartOut + k) = dimag(zOut(k))
        enddo

        do k = 1,nd
          zIn(k) = density(istartIn + k) + eye*dsda(istartIn + k)
        enddo
        call fourierUpsample(nd,wsaveIn,zIn,nd*nup,wsaveOut,zOut)
        do k = 1,nup*nd
          densityUp(istartOut + k) = dreal(zOut(k))
          dsdaUp(istartOut + k) = dimag(zOut(k))
        enddo

        istartIn = istartIn + nd
        istartOut = istartOut + nd*nup
      enddo
c     Upsample the boundary and density function for handling near
c     points

c     Compute the double-layer potential at the target
c     points using hypergeometric functions.
      do i = 1,nth
        do j = 1,nphi
          u_gr(i,j) = 0.d0
        enddo
      enddo
c     initialize solution to be 0


      tbeg = etime(timep)
      istart = 0
      istartUp = 0
c     loop over source points


c      do i = 1,nth
c        do j = 1,nphi
c          if (near_gr(i,j,1) .eq. 2 .and. near_gr(i,j,2) .eq. 2) then
c            u_gr(i,j) = 3.d0
c          elseif (near_gr(i,j,1) .eq. 1 .or. near_gr(i,j,2) .eq. 1) then
c            u_gr(i,j) = 2.d0
c          elseif (near_gr(i,j,1) .eq. 0 .or. near_gr(i,j,2) .eq. 0) then
c            u_gr(i,j) = 1.d0
c          else
c            u_gr(i,j) = 0.d0
c          endif
c        enddo
c      enddo

      do ihole = 1,nhole
        do i = 1,nth
          do j = 1,nphi
            if (near_gr(i,j,ihole) .eq. 2 .and. igrid(i,j) .ne. 0) then
              do k = 1,nd
                dist2 = (x_gr(i,j) - xs(k + istart))**2.d0 + 
     1                  (y_gr(i,j) - ys(k + istart))**2.d0 + 
     2                  (z_gr(i,j) - zs(k + istart))**2.d0 
                rdotn = (x_gr(i,j) - xs(k + istart))*xn(k + istart) + 
     1                  (y_gr(i,j) - ys(k + istart))*yn(k + istart) + 
     2                  (z_gr(i,j) - zs(k + istart))*zn(k + istart)

                u_gr(i,j) = u_gr(i,j) + 
     1            yukawaDLP(freq,dist2,rdotn)*
     2            density(k + istart)*dsda(k + istart)*dalph
              enddo !k
            elseif (near_gr(i,j,ihole) .eq. 1 .and. 
     1            igrid(i,j) .ne. 0) then
              do k = 1,nd*nup
                dist2 = (x_gr(i,j) - xsUp(k + istartUp))**2.d0 + 
     1                  (y_gr(i,j) - ysUp(k + istartUp))**2.d0 + 
     2                  (z_gr(i,j) - zsUp(k + istartUp))**2.d0 
                rdotn = (x_gr(i,j) - xsUp(k + istartUp))*
     1                      xnUp(k + istartUp) + 
     2                  (y_gr(i,j) - ysUp(k + istartUp))*
     3                      ynUp(k + istartUp) + 
     4                  (z_gr(i,j) - zsUp(k + istartUp))*
     5                      znUp(k + istartUp)

                u_gr(i,j) = u_gr(i,j) + 
     1            yukawaDLP(freq,dist2,rdotn)*
     2            densityUp(k+istartUp)*dsdaUp(k+istartUp)*dalph/nup
              enddo !k
            elseif (near_gr(i,j,ihole) .eq. 0 .and. 
     1          igrid(i,j) .ne. 0) then
              u_gr(i,j) = 1.d-5

              call nearInterp(freq,nd*nup,
     1          x_gr(i,j),y_gr(i,j),z_gr(i,j),
     2          xsUp(istartUp+1),ysUp(istartUp+1),zsUp(istartUp+1),
     3          xnUp(istartUp+1),ynUp(istartUp+1),znUp(istartUp+1),
     4          densityUp(istartUp+1),dsdaUp(istartUp+1),
     5          arcLength(ihole), valNear)

              u_gr(i,j) = u_gr(i,j) + valNear

            else
c              u_gr(i,j) = 1.d0
              u_gr(i,j) = 0.d0
            endif
          enddo !j
        enddo !i
        istart = istart + nd
        istartUp = istartUp + nd*nup
      enddo !ihole
      tend = etime(timep)


      do i = 1,nth
        do j = 1,nphi
          do ihole = 1,nhole
            if (near_gr(i,j,ihole) .eq. 0) then
              u_gr(i,j) = 1.d-5
            endif
          enddo
        enddo
      enddo


      call PRIN2 (' TIME FOR GRID = *',tend-tbeg,1)
c
c  write solution on meshgrid to a dat file
      open (unit=43, file = 'ugrid.dat')
      call dump (nth, nphi, u_gr, igrid, 1, 43)
      close (43)

      return
      end ! SOL_GRID
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine nearInterp(freq,nd,x_gr,y_gr,z_gr,
     1      xs,ys,zs,xn,yn,zn,density,dsda,arcLength,valNear)
c
c  *** DESCRIPTION :
c
c    Do a one-dimensional Lagrange interplation to evaluate the
c    layer potential at (x_gr,y_gr,z_gr) which is in the near region
c    of (xs,ys,zs)
c
c   *** INPUT PARAMETERS :
c
c   freq     = constant in the pde (\Delta - freq^2)u = 0
c   nd       = total number of points in linear system
c   x_gr     = x-coordinate of the target point
c   y_gr     = y-coordinate of the target point
c   z_gr     = z-coordinate of the target point
c   xs       = x-coordinate of the source points
c   ys       = y-coordinate of the source points
c   zs       = z-coordinate of the source points
c   xn       = x-coordinate of the normal vector
c   yn       = y-coordinate of the normal vector
c   zn       = z-coordinate of the normal vector
c   density  = density function
c   dsda     = Jacobian
c
c   *** OUTPUT PARAMETERS :
c
c   valNear  = value of the double-layer potential due to
c              (xs,ys,zs)
c
c***********************************************************************
c
      implicit real *8 (a-h,o-z)
c
      dimension xs(nd),ys(nd),zs(nd)
      dimension xn(nd),yn(nd),zn(nd)
      dimension density(nd),dsda(nd)
c
      dimension xLag(10),yLag(10),zLag(10),valLag(10)

      twopi = 8.d0*datan(1.d0)
      dalph = twopi/dble(nd)

      distMin = 10.d0
      do k = 1,nd
        dist2 = (xs(k)-x_gr)**2.d0 + (ys(k)-y_gr)**2.d0 + 
     1      (zs(k) - z_gr)**2.d0
        if (dist2 .lt. distMin) then
          iCP = k
          distMin = dist2
        endif
      enddo
c     Find the location of the closest point

      dirx = x_gr - xs(iCP)
      diry = y_gr - ys(iCP)
      dirz = z_gr - ys(iCP)
      dirNorm = dsqrt(dirx**2.d0 + diry**2.d0 + dirz**2.d0)
      dirx = dirx/dirNorm*arcLength
      diry = diry/dirNorm*arcLength
      dirz = dirz/dirNorm*arcLength
c     Compute direction vector that joins the source and target point
c     Note that these points will not be on the sphere, but this should
c     not effect the accuracy

c     number of Lagrange interpolation points to use
      nLag = 6
      do k = 1,nLag
        xLag(k) = xs(iCP) + dble(k-1)*dirx
        yLag(k) = ys(iCP) + dble(k-1)*diry 
        zLag(k) = zs(iCP) + dble(k-1)*dirz 
        xNorm = dsqrt(xLag(k)**2.d0 + yLag(k)**2.d0 + 
     1        zLag(k)**2.d0)
        xLag(k) = xLag(k)/xNorm
        yLag(k) = yLag(k)/xNorm
        zLag(k) = zLag(k)/xNorm
c       Project Largrange point onto the sphere
c        print*,xLag(k),',',yLag(k),',',zLag(k)
      enddo

      valLag(1) = 0.d0
c     This will come from evaluting the layer potential later

      do n = 2,nLag
        valLag(n) = 0.d0
        do k = 1,nd
          dist2 = (xLag(n) - xs(k))**2.d0 + 
     1            (yLag(n) - ys(k))**2.d0 + 
     2            (zLag(n) - zs(k))**2.d0 
          rdotn = (xLag(n) - xs(k))*xn(k) + 
     1            (yLag(n) - ys(k))*yn(k) + 
     2            (zLag(n) - zs(k))*zn(k)

          valLag(n) = valLag(n) + 
     1        yukawaDLP(freq,dist2,rdotn)*
     2        density(k)*dsda(k)*dalph
        enddo !k
c        print*,xLag(n),yLag(n),zLag(n),valLag(n)

      enddo





      valNear = 0.d0


      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c
c  *** DESCRIPTION :
c
c   Another routine required by DGMRES. It allows the use of a
c   preconditioner.
c
c   *** INPUT PARAMETERS :
c
c   n        = total number of points in linear system
c   z        = right-hand side z of A*r = z
c   nelt     = number of elements in matrix A (this doesn't matter to
c              us)
c   ia       = data structure for A (this doesn't matter to us)
c   ja       = data structure for A (this doesn't matter to us)
c   a        = matrix A in yy = A*xx.  In general, method is matrix
c              free, so this is never constructed
c   isym     = flag indicating if the problem is symmetric
c   rwork    = real workspace array
c   iwork    = integer workspace array
c
c   *** OUTPUT PARAMETERS :
c
c   r        = solution r of A*r = z
c
c***********************************************************************
c
      implicit real *8 (a-h,o-z)
c
      dimension r(n), z(n)
c
      do i = 1,n
        r(i) = z(i)
      enddo
c     no preconditioning for now
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine dump (nx,ny,ugrid,igrid,ireal,ifile)
c
c  *** DESCRIPTION :
c
c   Write variables to dat files
c
c   *** INPUT PARAMETERS :
c
c   nx       = number of points in the first coordinate
c   ny       = number of points in the second coordinate
c   ugrid    = data to be written
c   igrid    = location of points that are in domain of interest
c   ireal    = flag to indicate of data is real or not
c   ifile    = file ID where information is dumped
c
c   *** OUTPUT PARAMETERS :
c
c   NONE
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension ugrid(nx,ny), igrid(nx,ny)
c
      do i = 1,nx
        do j = 1, ny
          if (ireal.eq.1) then 
            write(ifile,'(e20.13,$)')(ugrid(I,J))
            write (ifile,'(a)')  ''
          else
            write(ifile,'(i4,$)') (igrid(i,j))
            write (ifile,'(a)')  ''
          end if
        end do
      end do
c
c periodic
      do j = 1, ny
        if (ireal.eq.1) then 
          write(ifile,'(e20.13,$)')(ugrid(1,J))
          write (ifile,'(a)')  ''
        else
          write(ifile,'(i4,$)') (igrid(1,j))
          write (ifile,'(a)')  ''
        end if
      end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine CHECK_ERROR (nhole, nth, nphi, near_gr, u, uExact)
c
c  *** DESCRIPTION :
c
c   Compute the maximum absolute and relative errors of the
c   solution of the PDE.  Exact solution is found by using
c   the boundary conditions that come from taking the fundamental
c   solution centered inside of each hole.
c
c   *** INPUT PARAMETERS :
c
c   k       = number of connected componenets in geometry
c   nth     = number of points in the azimuthal direction
c   nphi    = number of points in the z-uthal direction
c   near_gr = want to check error in different regions independently
c   u       = numerical solution 
c   uExact  = exact solution found using vorticies centered
c             outside of the geometry
c
c   *** OUTPUT PARAMETERS :
c
c   NONE
c
c***********************************************************************
c
      implicit real *8 (a-h,o-z)
c
      dimension igrid(nth,nphi), near_gr(nth,nphi,nhole)
      dimension u(nth,nphi), uExact(nth,nphi)
c
      uExactFarMax = 0.d0
      uExactIntermediateMax = 0.d0
      uExactNearMax = 0.d0
      absErrorFar = 0.d0
      absErrorIntermediate = 0.d0
      absErrorNear = 0.d0
      do i = 1,nth
        do j = 1,nphi
          ifar = 1
          iint = 0
          inear = 0
c         flags to determine if a point is far, intermediate, or near         
          do k = 1,nhole
            if (near_gr(i,j,k) .eq. 1) then
              iint = 1
              ifar = 0
            elseif (near_gr(i,j,k) .eq. 0) then
              inear = 1
              ifar = 0
            endif
          enddo

          if (ifar .eq. 1) then
            uExactFarMax = max(uExactFarMax, dabs(uExact(i,j)))
            absErrorFar = max(absErrorFar,
     1          dabs(uExact(i,j) - u(i,j)))
          elseif (iint .eq. 1) then
            uExactIntermediateMax = max(uExactIntermediateMax, 
     1              dabs(uExact(i,j)))
            absErrorIntermediate = max(absErrorIntermediate,
     1          dabs(uExact(i,j) - u(i,j)))
          elseif (inear .eq. 1) then
            uExactNearMax = max(uExactNearMax, dabs(uExact(i,j)))
            absErrorNear = max(absErrorNear,
     1          dabs(uExact(i,j) - u(i,j)))
          endif

        enddo
      enddo

      relErrorFar = absErrorFar/uExactFarMax
      relErrorIntermediate = absErrorIntermediate/uExactIntermediateMax
      relErrorNear = absErrorNear/uExactNearMax


      if (uExactFarMax .gt. 0.d0) then
c        call PRIN2(' Absolute Error (Far) = *',absErrorFar, 1)
        call PRIN2(' Relative Error (Far)= *', relErrorFar, 1)
      else
        print*,'Far region is empty'
      endif
      if (uExactIntermediateMax .gt. 0.d0) then
c        call PRIN2(' Absolute Error (Intermediate) = *',
c     1      absErrorIntermediate, 1)
        call PRIN2(' Relative Error (Intermediate)= *',
     1      relErrorIntermediate, 1)
      else
        print*,'Intermediate region is empty'
      endif
c      if (uExactNearMax .gt. 0.d0) then
c        call PRIN2(' Absolute Error (Near) = *',absErrorNear, 1)
c        call PRIN2(' Relative Error (Near)= *', relErrorNear, 1)
c      else
c        print*,'Near region is empty'
c      endif

      return ! CHECK_ERROR
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine dump_movie (nx, ny, xgrid, ygrid, zgrid, 
     1          ugrid, uExact, ifile)
c
c  *** DESCRIPTION :
c
c   Write a matlab file that plots the sphere with the colormap
c   correpsonding to the solution of the PDE
c
c   *** INPUT PARAMETERS :
c
c   nx       = number of points in the first coordinate
c   ny       = number of points in the second coordinate
c   xgrid    = x-coordinate of grid on sphere
c   ygrid    = y-coordinate of grid on sphere
c   zgrid    = z-coordinate of grid on sphere
c   ugrid    = solution which is a function defined on the sphere
c   uExact   = exact solution which is a function defined on the sphere
c   ifile    = file ID where information is dumped
c
c   *** OUTPUT PARAMETERS :
c
c   NONE
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension xgrid(nx,ny), ygrid(nx,ny), zgrid(nx,ny)
      dimension ugrid(nx,ny), uExact(nx,ny)
c
      write (ifile,*) 'nx = ',nx,';'
      write (ifile,*) 'ny = ',ny,';'
      write (ifile,*) 'xgrid = zeros(nx+1,ny);'
      write (ifile,*) 'ygrid = zeros(nx+1,ny);'
      write (ifile,*) 'zgrid = zeros(nx+1,ny);'
      write (ifile,*) 'ugrid = zeros(nx+1,ny);'
      write (ifile,*) 'uExact = zeros(nx+1,ny);'
c  periodically write the x-coordinate of the meshgrid
      do i = 1,nx
        do j = 1,ny
          write(ifile,*) 'xgrid(',i,',',j,') = ',xgrid(i,j),';'
        enddo
      enddo
      do j = 1,ny
        write(ifile,*) 'xgrid(',nx+1,',',j,') = ',xgrid(1,j),';'
      enddo
c  periodically write the y-coordinate of the meshgrid
      do i = 1,nx
        do j = 1,ny
          write(ifile,*) 'ygrid(',i,',',j,') = ',ygrid(i,j),';'
        enddo
      enddo
      do j = 1,ny
        write(ifile,*) 'ygrid(',nx+1,',',j,') = ',ygrid(1,j),';'
      enddo
c  periodically write the z-coordinate of the meshgrid
      do i = 1,nx
        do j = 1,ny
          write(ifile,*) 'zgrid(',i,',',j,') = ',zgrid(i,j),';'
        enddo
      enddo
      do j = 1,ny
        write(ifile,*) 'zgrid(',nx+1,',',j,') = ',zgrid(1,j),';'
      enddo
c  periodically write the solution on the meshgrid
      do i = 1,nx
        do j = 1,ny
          write(ifile,*) 'ugrid(',i,',',j,') = ',ugrid(i,j),';'
        enddo
      enddo
      do j = 1,ny
        write(ifile,*) 'ugrid(',nx+1,',',j,') = ',ugrid(1,j),';'
      enddo
c  periodically write the exact solution on the meshgrid
      do i = 1,nx
        do j = 1,ny
          write(ifile,*) 'uExact(',i,',',j,') = ',uExact(i,j),';'
        enddo
      enddo
      do j = 1,ny
        write(ifile,*) 'uExact(',nx+1,',',j,') = ',uExact(1,j),';'
      enddo

c  write Matlab code that makes pretty surface plot
      write (ifile,*) 'figure(1); clf;'
      write (ifile,*) '   surf(xgrid,ygrid,zgrid,ugrid)'
      write (ifile,*) '   view([64,-4])'
      write (ifile,*) '   shading interp'
      write (ifile,*) '   colormap(jet)'
      write (ifile,*) '   lighting phong'
      write (ifile,*) '   material dull'
      write (ifile,*) '   camlight(''headlight'')'
      write (ifile,*) '   colorbar'
      write (ifile,*) '   hold on'
      write (ifile,*) '   axis([-1 1 -1 1 -1 1])'
      write (ifile,*) '   axis equal'
      write (ifile,*) '   axis off'
c
      return
      end

