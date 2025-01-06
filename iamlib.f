C
C   Herbert M. Pickett, 24 Feb 1989
C   Adapted from subset of BLAS routines in LINPAK 
C
      FUNCTION IDAMAX(N,SX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C
      REAL*8  SX(0:*),SMAX
      INTEGER IDAMAX,N,INCX,I,M
      INTEGER*4 IX
C     EMA SX
C
      IDAMAX = 0
      M=N-1
      IF(M.LE.0) THEN
         IF(M.LT.0) IDAMAX=-1
      ELSE IF(INCX.NE.1) THEN
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
         SMAX = ABS(SX(0))
         IX = INCX
         DO 10 I = 1,M
             IF(ABS(SX(IX)).GT.SMAX) THEN
               IDAMAX = I
               SMAX = ABS(SX(IX))
             ENDIF   
             IX = IX + INCX
   10        CONTINUE
      ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
C
         SMAX = ABS(SX(0))
         DO 30 I = 1,M
            IF(ABS(SX(I)).GT.SMAX) THEN
              IDAMAX = I
              SMAX = ABS(SX(I))
            ENDIF  
   30       CONTINUE
      ENDIF
      IDAMAX=IDAMAX+1
      RETURN
      END
      FUNCTION DASUM(N,SX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C
      REAL*8  DASUM,SX(0:*)
      INTEGER I,INCX,N,M
      INTEGER*4 IX
C     EMA SX
C
      DASUM = 0
      M=N-1
      IF(M.GE.0) THEN
        IF(INCX.NE.1) THEN
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
          IX=0
          DO 10 I = 1,N
            DASUM = DASUM + ABS(SX(IX))
            IX = IX + INCX
   10     CONTINUE
        ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C       
          DO 30 I = 0,M
            DASUM = DASUM + ABS(SX(I))
   30     CONTINUE
        ENDIF
      ENDIF  
      RETURN
      END
      SUBROUTINE DAXPY(N,SA,SX,INCX,SY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C
      REAL*8  SX(0:*),SY(0:*),SA
      INTEGER I,INCX,INCY,N,M
      INTEGER*4 IX,IY
C     EMA SX,SY
C
      M=N-1
      IF(M.LT.0)RETURN
      IF(SA.EQ.0) RETURN
      IF(INCX.NE.1)GO TO 20
      IF(INCY.NE.1)GO TO 20
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
      DO 30 I = 0,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      RETURN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
   20 IX = 0
      IF(INCX.LT.0) IX = -M*INCX
      IY = 0
      IF(INCY.LT.0) IY = -M*INCY
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
      SUBROUTINE  DCOPY(N,SX,INCX,SY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C
      REAL*8  SX(0:*),SY(0:*),SA
      INTEGER I,INCX,INCY,N,M
      INTEGER*4 IX,IY
C     EMA SX,SY
C
      M=N-1
      IF(M.LT.0)RETURN
      IF(INCX.EQ.0) THEN
        SA=SX(0)
        IF(INCY.NE.1) THEN
C
C        CODE FOR FILL AND INCREMENT NOT EQUAL TO 1
C
          IY = 0
          IF(INCY.LT.0) IY = -M*INCY
          DO 10 I = 1,N
            SY(IY) = SA
            IY = IY + INCY
   10       CONTINUE
        ELSE
C
C        CODE FOR FILL AND INCREMENT EQUAL TO 1
C
C
          DO 20 I = 0,M
            SY(I) = SA
   20       CONTINUE
        ENDIF
      ELSE IF(INCX.NE.1.OR.INCY.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
        IX = 0
        IF(INCX.LT.0)IX = -M*INCX
        IY = 0
        IF(INCY.LT.0)IY = -M*INCY
        DO 30 I = 1,N
          SY(IY) = SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   30     CONTINUE
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
        DO 40 I = 0,M
          SY(I) = SX(I)
   40     CONTINUE
      ENDIF
      RETURN
      END
      FUNCTION DDOT(N,SX,INCX,SY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C
      REAL*8  DDOT,SX(0:*),SY(0:*)
      INTEGER I,INCX,INCY,N,M
      INTEGER*4 IX,IY
C     EMA SX,SY
C
      DDOT = 0
      M=N-1
      IF(M.LT.0)RETURN
      IF(INCX.NE.1)GO TO 20
      IF(INCY.NE.1)GO TO 20
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
      DO 30 I = 0,M
        DDOT = DDOT + SX(I)*SY(I)
   30 CONTINUE
      RETURN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
   20 IX = 0
      IF(INCX.LT.0)IX = -M*INCX
      IY = 0
      IF(INCY.LT.0)IY = -M*INCY
      DO 10 I = 1,N
        DDOT = DDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
      SUBROUTINE  DSCAL(N,SA,SX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C
      REAL*8  SA,SX(0:*)
      INTEGER I,INCX,N,M
      INTEGER*4 IX
C     EMA SX
C     
      M=N-1
      IF(M.LT.0)  RETURN
      IF(SA.EQ.1) RETURN
      IF(INCX.EQ.1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX=0
      DO 10 I = 1,N
        SX(IX) = SA*SX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 0,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE  DSWAP (N,SX,INCX,SY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C
      REAL*8  SX(0:*),SY(0:*),STEMP
      INTEGER I,INCX,INCY,N,M
      INTEGER*4 IX,IY
C     EMA SX,SY
C
      M=N-1
      IF(M.LT.0) RETURN
      IF(INCX.NE.1) GO TO 20
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
      DO 30 I = 0,M
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
   30 CONTINUE
      RETURN
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
   20 IX = 0
      IF(INCX.LT.0)IX = -M*INCX
      IY = 0
      IF(INCY.LT.0)IY = -M*INCY
      DO 10 I = 1,N
        STEMP = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
      SUBROUTINE  DROT (N,SX,INCX,SY,INCY,C,S)
C
C     APPLIES A PLANE ROTATION.
C
      REAL*8  SX(0:*),SY(0:*),C,S,STEMP
      INTEGER INCX,INCY,N,I,M
      INTEGER*4 IX,IY
C     EMA SX,SY
C
      M=N-1
      IF(M.LT.0)RETURN
      IF(INCX.NE.1) GO TO 20
      IF(INCY.NE.1) GO TO 20
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
      DO 30 I = 0,M
        STEMP = SX(I)*C + SY(I)*S
        SY(I) = SY(I)*C - SX(I)*S
        SX(I) = STEMP
   30 CONTINUE
C
      RETURN
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
   20 IX = 0
      IF(INCX.LT.0)IX = -M*INCX
      IY = 0
      IF(INCY.LT.0)IY = -M*INCY
      DO 10 I = 1,N
        STEMP  = SX(IX)*C + SY(IY)*S
        SY(IY) = SY(IY)*C - SX(IX)*S
        SX(IX) = STEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END

      SUBROUTINE EIGSRT(D,V,N,NP)
C     sort eigenvectors in ascending order
      implicit none
      INTEGER N,NP,I,J,K
      REAL*8 D(NP),V(NP,NP),P
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END


      SUBROUTINE HDIAG(NM,NX,Z,D,E,IERR)
C
C
      INTEGER   NM,NX,IERR
      REAL*8    Z(NM,*),D(*),E(*)
C      EMA Z,D,E
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     MODIFIED TO SPEED UP FOR SPARCE MATRICIES
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     THE SUBROUTINE THEN CALLS TRIAG WHICH DIAGONALIZES
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        NX IS THE ORDER OF THE MATRIX,
C
C        Z CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON CALL TO TRIAG -
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY,
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION.
C
C     ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ARBITRARY ORDER
C
C        Z CONTAINS THE EIGENVECTORS
C
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
      INTEGER I0,I1
      PARAMETER (I0=0,I1=1)
      REAL*8 MACHEP,F,G,H,T,BB,DDOT,ONE,ZERO
      PARAMETER (ONE=1.,ZERO=0.)
      INTEGER I,J,K,L,IZ,LZ,NT1,NT2,NDM,N
      MACHEP=1.E-30
      N=NX
      NDM=NM
C     ********** FOR I=N STEP -1 UNTIL 3 DO -- **********
      DO 300 I = N,3,-1
         L = I - 1
         J = L - 1
         CALL DCOPY(L,Z(I,1),NDM,D,I1)
         F=D(L)
         IZ=0
         H=0
C        GET SUM OF SQUARES OF ELEMENTS AND FIND FIRST NON-ZERO ELEMENT
         DO 220 K=J,1,-1
             G=D(K)*D(K)
             IF(G.GT.MACHEP) THEN
                 IZ=K
                 H=H+G
             ELSE
                 D(K)=0
             ENDIF
  220        CONTINUE
         E(I)=F
         IF(IZ.NE.0) THEN
             T=F*F
             IF(H.LT.MACHEP*T) IZ=0
         ENDIF
         IF(IZ.EQ.0) THEN
             D(L)=0
         ELSE
             G= SIGN( SQRT( H +T ), F)
             E(I)= -G
             F= F + G
             D(L)=F
             BB=SIGN( ONE, F ) / SQRT( F*G )
             LZ=I-IZ
             CALL DSCAL(LZ,BB,D(IZ),I1)
             CALL DCOPY(L,D,I1,Z(1,I),I1)
C
             DO 230 J= 1, IZ
               E(J)= DDOT(LZ,Z(IZ,J),I1,D(IZ),I1)
  230          CONTINUE                           !Herbers2024, Warning removal...
             F=-E(IZ)*D(IZ)
             DO 240 J = IZ+1, L
                NT1=J-IZ
                NT2=I-J
                T= DDOT(NT1,Z(J,IZ),NDM,D(IZ),I1)
     +            +DDOT(NT2,Z(J, J), I1,D( J),I1)
                E(J)=T
                F=F-T*D(J)
  240           CONTINUE                          !Herbers2024, Warning removal
C
C     ********** FORM REDUCED A **********
             DO 260 J = IZ, L
                H = -D(J)
                G =  F * H - E(J)
                K =  I-J
                IF(ABS(H).GT.MACHEP) CALL DAXPY(K,H,E(J),I1,Z(J,J),I1)
                CALL DAXPY(K,G,D(J),I1,Z(J,J),I1)
  260           CONTINUE
             DO 265 J = 1,IZ-1
                G = -E(J)
                CALL DAXPY(LZ,G,D(IZ),I1,Z(IZ,J),I1)
  265           CONTINUE
          ENDIF
  300     CONTINUE
C
      D(1)=Z(1,1)
      Z(1,1)=ONE
      IF(N.LE.1) THEN
         IERR=0
         RETURN
      ENDIF
      D(N)=0
      E(2)=Z(2,1)
C     ********** ACCUMULATION OF TRANSFORMATION MATRICES **********
      DO 500 L = 2, N
         I = L + 1
         J = L - 1
         E(J)=E(L)
         F=D(L)
         D(L)=Z(L,L)
         IF(F.GT.MACHEP) THEN
             DO 360 K=1,J
                G= -DDOT(J,Z(1,I),I1,Z(1,K),I1)
                IF(ABS(G).GT.MACHEP) CALL DAXPY(J,G,Z(1,I),I1
     +                                             ,Z(1,K),I1)
                Z(L,K)=   F * G
                Z(K,L)= - F * Z(K,I)
  360           CONTINUE
             Z(L,L) = ONE -F*F
         ELSE !By now rank-1 is not equal to scalar...2023
             Z(L,L) = ZERO                       !Sven2023
             CALL DCOPY(J,Z(L,L),I0,Z(1,L),I1)   !Sven2023
             CALL DCOPY(J,Z(L,L),I0,Z(L,1),NDM)  !Sven2023
             Z(L,L) = ONE                        !Sven2023
         ENDIF
  500    CONTINUE
C
      CALL TRIAG(NM,N,Z,D,E,IERR)
      RETURN
      END


      SUBROUTINE TRIAG(NM,NX,Z,D,E,IERR)
C
      INTEGER   NM,NX,IERR
      REAL*8    Z(0:*),D(*),E(*)
C      EMA Z,D,E
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE imtql2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        NX IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ARBITRARY ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      INTEGER I1
      PARAMETER (I1=1)
      REAL*8 MACHEP,F,B,P,G,R,C,S
      INTEGER*4 INDX1,INDX2
      INTEGER N,NN,ITR,M,I,NDM,L,MM
      MACHEP=1.E-15
      N=NX
C
      IERR = 0
      IF (N .LE. 1) RETURN
      NN=N-1
      NDM=NM
      E(N) = 0
      F = 0
      B = 0
C
      DO 240 L = 1, NN
        ITR = 31
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
        DO WHILE(ITR.NE.0)
         M=N
         DO 110 I = L, NN
            IF (ABS(E(I)).LE.MACHEP*(ABS(D(I))+ABS(D(I+1)))) THEN
               M=I
               GO TO 120
            ENDIF
  110       CONTINUE
C
  120    P=D(L)
         IF (M .EQ. L) GO TO 240
         ITR = ITR - 1
         IF (ITR .EQ. 0) THEN
C        ** SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER 30 ITERATIONS **
            IERR = L
            RETURN
         ENDIF
C     ********** FORM SHIFT **********
         G = E(L)
         G = ( D(L+1)-P )/( G + G )
         R = SQRT(1+G*G)
         G = D(M) - P +  E(L) / ( G +SIGN(R,G) )
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C = 1
         S = C
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         MM=M-1
         INDX2=NDM*MM
         DO 200 I = MM,L,-1
            B=E(I)*C
            F=E(I)*S
            IF(ABS(F).GE.ABS(G)) THEN
              C=G/F
              R=SQRT(1+C*C)
              E(I+1)=F*R
              S=1/R
              C=C*S
            ELSE
              S=F/G
              R=SQRT(1+S*S)
              E(I+1)=G*R
              C=1/R
              S=S*C
            ENDIF
            F=C*D(I)-S*B
            G=C*B-S*P
            R=D(I)+P
            P=C*F-S*G
            G=S*F+C*G
            D(I+1)=R-P
C     ********** FORM VECTOR **********
            INDX1=INDX2
            INDX2=INDX2-NDM
            CALL DROT(N,Z(INDX1),I1,Z(INDX2),I1,C,S)
C
  200    CONTINUE
C
         D(L) = P
         E(L) = G
         E(M) = 0
        ENDDO
  240 CONTINUE
      RETURN
      END
      SUBROUTINE HEIGSRT(D,VR,VI,N,NP)
C     sort hermitian eigenvectors in ascending order 
      implicit real*8 (a-h,o-z)
      DIMENSION D(NP),VR(NP,NP),VI(NP,NP)
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=VR(J,I)
            VR(J,I)=VR(J,K)
            VR(J,K)=P
            P=VI(J,I)
            VI(J,I)=VI(J,K)
            VI(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END
      SUBROUTINE HTRIB3(NM,N,A,TAU,M,ZR,ZI)
C
      INTEGER I,J,K,L,M,N,NM
      DOUBLE PRECISION A(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      DOUBLE PRECISION H,S,SI
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRID3.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION BY  HTRID3.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
      DO 60 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50    CONTINUE
   60 CONTINUE             !Herbers2024 warning removal.
C
      IF (N .EQ. 1) GO TO 200
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = A(I,I)
         IF (H .EQ. 0.0D0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0D0
            SI = 0.0D0
C
            DO 110 K = 1, L
               S = S + A(I,K) * ZR(K,J) - A(K,I) * ZI(K,J)
               SI = SI + A(I,K) * ZI(K,J) + A(K,I) * ZR(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * A(I,K) - SI * A(K,I)
               ZI(K,J) = ZI(K,J) - SI * A(I,K) + S * A(K,I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
C
C
C
      SUBROUTINE HTRID3(NM,N,A,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JM1,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N),TAU(2,N)
      DOUBLE PRECISION F,G,H,FI,GI,HH,SI,SCALE,PYTHAG
C
C     DIAGONALISIERUNG MIT EISPACK- ROUTINEN
C
C        IERR=0
C        CALL HTRID3 (NARZ,N,ZM1,D,E,E2,TAU)
C        CALL TQL2 (NARZ,N,D,E,ZM2,IERR)
C        IF (IERR.NE.0) THEN
C         WRITE (0,513) IERR
C  513    FORMAT ('FEHLER BEI EIGENWERT:',I5)
C         STOP
C        ENDIF
C        CALL HTRIB3 (NARZ,N,ZM1,TAU,NARZ,ZM2,ZM3)
C
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX, STORED AS
C     A SINGLE SQUARE ARRAY, TO A REAL SYMMETRIC TRIDIAGONAL MATRIX
C     USING UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE COMPLEX HERMITIAN INPUT
C          MATRIX.  THE REAL PARTS OF THE MATRIX ELEMENTS ARE STORED
C          IN THE FULL LOWER TRIANGLE OF A, AND THE IMAGINARY PARTS
C          ARE STORED IN THE TRANSPOSED POSITIONS OF THE STRICT UPPER
C          TRIANGLE OF A.  NO STORAGE IS REQUIRED FOR THE ZERO
C          IMAGINARY PARTS OF THE DIAGONAL ELEMENTS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
         SCALE = SCALE + DABS(A(I,K)) + DABS(A(K,I))
  120    CONTINUE      !Herbers 2024 warning removal
C
         IF (SCALE .NE. 0.0D0) GO TO 140
         TAU(1,L) = 1.0D0
         TAU(2,L) = 0.0D0
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            A(K,I) = A(K,I) / SCALE
            H = H + A(I,K) * A(I,K) + A(K,I) * A(K,I)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = DSQRT(H)
         E(I) = SCALE * G
         F = PYTHAG(A(I,L),A(L,I))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0D0) GO TO 160
         TAU(1,L) = (A(L,I) * TAU(2,I) - A(I,L) * TAU(1,I)) / F
         SI = (A(I,L) * TAU(2,I) + A(L,I) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0D0 + G / F
         A(I,L) = G * A(I,L)
         A(L,I) = G * A(L,I)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         A(I,L) = G
  170    F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
            GI = 0.0D0
            IF (J .EQ. 1) GO TO 190
            JM1 = J - 1
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, JM1
               G = G + A(J,K) * A(I,K) + A(K,J) * A(K,I)
               GI = GI - A(J,K) * A(K,I) + A(K,J) * A(I,K)
  180       CONTINUE
C
  190       G = G + A(J,J) * A(I,J)
            GI = GI - A(J,J) * A(J,I)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * A(I,K) - A(J,K) * A(K,I)
               GI = GI - A(K,J) * A(K,I) - A(J,K) * A(I,K)
  200       CONTINUE
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * A(I,J) - TAU(2,J) * A(J,I)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -A(J,I)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
            A(J,J) = A(J,J) - 2.0D0 * (F * G + FI * GI)
            IF (J .EQ. 1) GO TO 260
            JM1 = J - 1
C
            DO 250 K = 1, JM1
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
     X                         + FI * TAU(2,K) + GI * A(K,I)
               A(K,J) = A(K,J) - F * TAU(2,K) - G * A(K,I)
     X                         - FI * E(K) - GI * A(I,K)
  250       CONTINUE
C
  260    CONTINUE
C
  270    DO 280 K = 1, L
            A(I,K) = SCALE * A(I,K)
            A(K,I) = SCALE * A(K,I)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    D(I) = A(I,I)
         A(I,I) = SCALE * DSQRT(H)
  300 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
      E(I-1) = E(I)
  100 CONTINUE       !Herbers 2024 warning removal
C
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
         D(I) = D(I) - H
  140    CONTINUE
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END

