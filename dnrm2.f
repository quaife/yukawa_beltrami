      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
c     .. Scalar Arguments ..
      INTEGER INCX,N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION X(*)
c     ..
c
c  Purpose
c  =======
c
c  DNRM2 returns the euclidean norm of a vector via the function
c  name, so that
c
c     DNRM2 := sqrt( x'*x )
c
c  Further Details
c  ===============

c  -- This version written on 25-October-1982.
c     Modified on 14-October-1993 to inline the call to DLASSQ.
c     Sven Hammarling, Nag Ltd.
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
c     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
         NORM = ZERO
      ELSE IF (N.EQ.1) THEN
         NORM = ABS(X(1))
      ELSE
         SCALE = ZERO
         SSQ = ONE
c     The following loop is equivalent to this call to the LAPACK
c     auxiliary routine:
c     CALL DLASSQ( N, X, INCX, SCALE, SSQ )
c     
         DO 10 IX = 1,1 + (N-1)*INCX,INCX
            IF (X(IX).NE.ZERO) THEN
               ABSXI = ABS(X(IX))
               IF (SCALE.LT.ABSXI) THEN
                  SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                  SCALE = ABSXI
               ELSE
                  SSQ = SSQ + (ABSXI/SCALE)**2
               END IF
            END IF
 10      CONTINUE
         NORM = SCALE*SQRT(SSQ)
      END IF
c     
      DNRM2 = NORM
      RETURN
c
c     End of DNRM2.
c
      END
