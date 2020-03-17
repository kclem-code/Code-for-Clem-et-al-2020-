      PROGRAM slope
c     This program calculates the slope (trend in time) of any input 1-D 
c     timeseries.  Returns the slope b and the two-tailed probability that the 
c     slope is zero.

      parameter (iu1 = 10, io1 = 20, io2 = 21)
      integer ihead (8),dof,nx,ny
      character * 256 file1,file2,file3
      real tval,cov,var,sb,var1,cnt
      double precision, dimension(:,:), allocatable :: sumx,sumy,sumx2
      double precision, dimension(:,:), allocatable :: sumy2,sumxy 
      real, dimension(:,:), allocatable :: indat,b,sig,nt

      namelist /param/file1, file2, file3

      data file2 /' '/, file1 /' '/, file3 /' '/

      read (5, param, end = 100)
 100  continue

      open (iu1, file = file1, form = 'unformatted',
     &     access = 'sequential')
      open (io1, file = file2, form = 'unformatted',
     &     access = 'sequential')
      open (io2, file = file3, form = 'unformatted',
     &     access = 'sequential')


      write (*, '(a10, a80)') '  >       ', file1
      write (*, '(a10, a80)') '  >       ', file2
      write (*, '(a10, a80)') '  >       ', file3

      read(iu1) ihead

      nx=ihead(5)
      ny=ihead(6)
      rewind(iu1)

      allocate(indat(nx,ny),b(nx,ny),sig(nx,ny))
      allocate(sumx(nx,ny),sumx2(nx,ny),sumy(nx,ny),sumy2(nx,ny))
      allocate(sumxy(nx,ny),nt(nx,ny))


    
      do i=1,nx
      do j=1,ny
      sumx(i,j)=0.
      sumx2(i,j)=0.
      sumy(i,j)=0.
      sumy2(i,j)=0.
      sumxy(i,j)=0.
      nt(i,j)=0.
      enddo
      enddo

 300  read (iu1,end = 910) ihead


      read (iu1) indat
      
      do i=1,nx
       do j=1,ny 
       if(indat(i,j).gt.-8e33)then
       nt(i,j)=nt(i,j)+ 1.0
       sumx2(i,j)=sumx2(i,j)+nt(i,j)*nt(i,j)
       sumx(i,j)=sumx(i,j)+nt(i,j)
       sumy(i,j)=sumy(i,j)+indat(i,j)
       sumy2(i,j)=sumy2(i,j)+indat(i,j)*indat(i,j)
       sumxy(i,j)=sumxy(i,j)+nt(i,j)*indat(i,j)      
       endif
       enddo
      enddo

      goto 300

 910  continue



      do i=1,nx
      do j=1,ny
      cnt=nt(i,j)
      if(cnt.gt.2)then
      dof=cnt-2
      cov=cnt*sumxy(i,j)-sumx(i,j)*sumy(i,j)
      var=cnt*sumx2(i,j)-sumx(i,j)*sumx(i,j)
      var1=sumx2(i,j)-1.0/cnt*sumx(i,j)*sumx(i,j)

      b(i,j)=cov/var
      sb=sqrt((1/(cnt*(cnt-2))*(cnt*sumy2(i,j)-sumy(i,j)*sumy(i,j)-
     & (cov*cov/var)))/var1)
      tval=b(i,j)/sb
      if(tval.ge.0)then
      sig(i,j)=1-PROBT1(tval,tval,dof,2)
      else
      sig(i,j)=PROBT1(tval,tval,dof,2)-1
      endif       
      else
      b(i,j)=-9e33
      sig(i,j)=-9e33
      endif

c      print*,'test...dof,nt,b,sb,tval,sig..... '
c     & ,dof,nt,b(i,j),sb,tval,sig(i,j)

      enddo
      enddo

      write (io1) ihead
      write (io1) b
      write (io2) ihead
      write (io2) sig

 


 905  stop 'regular end'
 
      END


c**************************************************************************
C
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C
C     FUNCTION PROBND1.
C
C     DATE:    14-NOV-1989.
C     AUTHOR:  M. Abramowitz and I.A. Stegun, 1984. In: "Pocketbook of
C              Mathematical Functions", pp 407-408 (26.2.1 and 26.2.17).
C              Programmed by B.D. Santer, MPI, Hamburg.
C
C     Calculates the area under a normal curve (mean=0.0, variance=1.0) to
C     the right of X. The accuracy is better than 7.5*10.**-8.
C
      FUNCTION PROBND1(X)
      IMPLICIT REAL(A-H,O-Z)
C
C     ** INPUT **
      REAL X            ! Input Z-score
C
C     ** OUTPUT **
      REAL PROBND1      ! Area under normal curve to right of X
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C
      B1 =  0.319381530
      B2 = -0.356563782
      B3 =  1.781477937
      B4 = -1.821255978
      B5 =  1.330274429
      P  =  0.2316419
      PI =  4.0*ATAN(1.0)
      T  =  1.0/(1.0+P*X)
C
      TERM1   = B1*T+B2*(T**2)+B3*(T**3)+B4*(T**4)+B5*(T**5)
      Z       = (1.0/SQRT(2.0*PI))*EXP(-X*X/2.0)
      PROBND1 = 0.0
      IF(X.GT.7.)GO TO 1
      PROBND1 = Z*TERM1
    1 CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C
C     FUNCTION PROBT1.
C
C     DATE:    14-NOV-1989.
C     AUTHOR:  M. Abramowitz and I.A. Stegun, 1964. In: "Handbook of
C              Mathematical Functions". Programmed by T.M.L. Wigley, CRU,
C              Norwich.
C
C     The output is either the one- or two-tailed test area, i.e., the area
C     under a Students curve (with N degrees of freedom) to the right of
C     ABS(X) (one-tailed test), or twice this area (two-tailed test).
C
      FUNCTION PROBT1(Y,X,N,ID)
      IMPLICIT REAL(A-H,O-Z)
C
C     ** INPUT **
      REAL Y              ! Input calculated t-value (unchanged on output)
      REAL X              ! Absolute value of Y
      INTEGER N           ! Degrees of freedom
      INTEGER ID          ! Identifier for one- or two-tailed test
C
C     ** OUTPUT **
      REAL PROBT1         ! Significance level (p-value) for t-value
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C
      X  = ABS(Y)
      DN = FLOAT(N)
      U  = X*(1.0-0.25/DN)/SQRT(1.0+X*X*0.5/DN)
C
C     Calculate probability for one- or two-tailed test.
C
      IF(ID.EQ.2) PROBT1 = 2.0*PROBND1(U)
      IF(ID.EQ.1) PROBT1 = PROBND1(U)
      RETURN
      END
C
c*****************************************************************************


