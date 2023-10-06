*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine spheremean ( n, w, x, y, z, ax, ay, az)
      implicit none
      integer n
      double precision w(n), x(n), y(n), z(n), ax, ay, az
      intrinsic ABS, ATAN2, SQRT, COS, SIN, MAX
      EXTERNAL rwarn

*---  Compute the spherical weighted average of a set of input vectors.
*---  All inputs are assumed to be in a single hemisphere. This is a
*---  Fortran translation of the C++ version coded by Samuel R. Buss
*---  and available from the websites:
*---     https://mathweb.ucsd.edu/~sbuss/ResearchWeb/spheremean/
*---     https://web.archive.org/web/20230201151538/https://mathweb.ucsd.edu/~sbuss/ResearchWeb/spheremean/
*---  The code implements Algorithm A1 described in the following
*---  article by Samuel R. Buss and Jay Fillmore:
*---     "Spherical Averages and Applications to Spherical Splines and Interpolation." 
*---     ACM Transactions on Graphics 20 (2001) 95-126

      integer i, it, maxit
      double precision s, err, tol, anorm, theta, costheta, sintheta, 
     .                 lx, ly, lz, vx, vy, vz, dx, dy, dz,
     .                 axold, ayold, azold

      ax = 0.0d0
      ay = 0.0d0
      az = 0.0d0
      tol = 1e-8
      
      it = 0
      maxit = 10000

*---  Initial estimate for mean vector
      do i=1,n
         ax = ax + w(i)*x(i)
         ay = ay + w(i)*y(i)
         az = az + w(i)*z(i)
      enddo

      anorm = SQRT(ax*ax + ay*ay + az*az)

      if (anorm.eq.0.0d0) goto 200

      ax = ax / anorm
      ay = ay / anorm
      az = az / anorm

*---  Main loop to improve estimate
 100  dx = 0.0d0
      dy = 0.0d0
      dz = 0.0d0
      it = it + 1
      do i=1,n
*---  Compute the tangent vector from current mean vector to
*---  each input vector. Its length is equal to the spherical
*---  distance between the mean and the input
         costheta = ax*x(i) + ay*y(i) + az*z(i)
         vx = x(i) - costheta*ax
         vy = y(i) - costheta*ay
         vz = z(i) - costheta*az
         sintheta = SQRT(vx*vx + vy*vy + vz*vz)
         if (sintheta.eq.0.0d0) then
            lx = 0.0d0
            ly = 0.0d0
            lz = 0.0d0
         else
            theta = ATAN2(sintheta, costheta)
            s = theta / sintheta
            lx = s*vx
            ly = s*vy
            lz = s*vz
         endif
*---  Update the rotation vector
         dx = dx + w(i)*lx
         dy = dy + w(i)*ly
         dz = dz + w(i)*lz
      enddo
      
*---  Store the current mean vector estimate
      axold = ax
      ayold = ay
      azold = az
      
*---  Rotate current mean vector in the direction of (dx,dy,dz)
*---  to get new mean vector. The length of (dx,dy,dz) is the
*---  rotation angle. Note that the mean vector must be a unit
*---  vector and (dx,dy,dz) must be perpendicular to it.
      theta = SQRT(dx*dx + dy*dy + dz*dz)
      if (theta.gt.0) then
         dx = dx / theta
         dy = dy / theta
         dz = dz / theta
         costheta = COS(theta)
         sintheta = SIN(theta)
         ax = ax*costheta + dx*sintheta
         ay = ay*costheta + dy*sintheta
         az = az*costheta + dz*sintheta
      endif

*---  Break out of loop if rotation was small enough
      err = MAX(ABS(ax-axold), ABS(ay-ayold), ABS(az-azold))
      
      if (it.gt.maxit) then
         call rwarn("maximum iterations exceeded")
         goto 200
      endif

      if (err.gt.tol) goto 100

 200  continue
      END
