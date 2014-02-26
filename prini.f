      SUBROUTINE PRINI(IP1,IQ1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /WRITING/ IP, IQ
      IP=IP1
      IQ=IQ1
      RETURN
      END 

C
C
C
      SUBROUTINE PRIN2(MES,A2,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 MES(*)
      DIMENSION A2(*)
      COMMON /WRITING/ IP, IQ
      CALL MESSPR(MES,IP,IQ)
      IF(IP.NE.0 .AND. N.NE.0) THEN
        WRITE(IP,1400)(A2(J),J=1,N)
        WRITE (IP,*)
      ENDIF
      IF(IQ.NE.0 .AND. N.NE.0) THEN
        WRITE(IQ,1400)(A2(J),J=1,N)
        WRITE (IQ,*)
      ENDIF
 1400 FORMAT(6(2X,E12.5))
      RETURN
      END
C
C
C
      SUBROUTINE PRIND2(MES,A2,N)
C     Same as PRIN2 except that format 
C     gives double precision
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 MES(*)
      DIMENSION A2(*)
      COMMON /WRITING/ IP, IQ
      CALL MESSPR(MES,IP,IQ)
      IF(IP.NE.0 .AND. N.NE.0) THEN
        WRITE(IP,1450)(A2(J),J=1,N)
        WRITE (IP,*)
      ENDIF
      IF(IQ.NE.0 .AND. N.NE.0) THEN
        WRITE(IQ,1450)(A2(J),J=1,N)
        WRITE (IQ,*)
      ENDIF
 1450 FORMAT(6(2X,E23.16))
      RETURN
      END
C
C
C
C
      SUBROUTINE PRINF(MES,IA,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 MES(*)
      DIMENSION IA(*)
      COMMON /WRITING/ IP, IQ
      CALL MESSPR(MES,IP,IQ)
      IF(IP.NE.0 .AND. N.NE.0) THEN
        WRITE(IP,1600)(IA(J),J=1,N)
      ENDIF
      IF(IQ.NE.0 .AND. N.NE.0) THEN
        WRITE(IQ,1600)(IA(J),J=1,N)
      ENDIF


 1600 FORMAT(10(1X,I7))
      RETURN
      END
C
C
c
c
      SUBROUTINE MESSPR(MES,IP,IQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER *1 MES(*),AST
      DATA AST/'*'/
C
C     DETERMINE THE LENGTH OF THE MESSAGE
C
      I=0
      DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
	I1=I
 1400 CONTINUE
 1600 CONTINUE
      IF ( (I1.NE.0) .AND. (IP.NE.0) ) THEN
        WRITE(IP,1800) (MES(I),I=1,I1)
      ENDIF
      IF ( (I1.NE.0) .AND. (IQ.NE.0) ) THEN
        WRITE(IQ,1800) (MES(I),I=1,I1)
      ENDIF
 1800 FORMAT(1X,80A1)
      RETURN
      END
c
