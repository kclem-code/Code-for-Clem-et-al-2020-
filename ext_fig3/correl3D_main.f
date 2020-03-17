      PROGRAM slope
c     This program calculates the correlation of an input 1-D 
c     timeseries (text file, just data) with a 2D spatial field.
c     Returns the correlation coefficient r and the 1 minus the two-tailed 
c     probability that the correlation is zero (based on a t-test)
c     Modified by R. L. Fogt of Ohio Univ on 2/20/13

      integer in1,in2,out1,out2
      parameter (in1 = 10, in2 = 20, out1 = 21, out2 = 22)
      integer ihead (8),dof,nx,ny,levfirst,nz
      character * 256 file1,file2,file3,file4
      real tval,val,cnt,cntlev
      double precision, dimension(:,:,:), allocatable :: sumx,sumy,sumx2
      double precision, dimension(:,:,:), allocatable :: sumy2,sumxy 
      double precision var,var1,cov
      real, dimension(:,:), allocatable :: indat
      real, dimension(:,:,:), allocatable :: nt
      real, dimension(:,:), allocatable :: r,sig
      integer, dimension(:), allocatable :: nlev

      namelist /param/file1, file2, file3,file4

      data file2 /' '/, file1 /' '/, file3 /' '/,file4 /'   '/

      read (5, param, end = 100)
 100  continue

      open (in1, file = file1, form = 'formatted')   !input timeseries
      open (in2, file = file2, form = 'unformatted',
     &     access = 'sequential') !input spatial field
      open (out1, file = file3, form = 'unformatted',
     &     access = 'sequential') !output correlation coefficient
      open (out2, file = file4, form = 'unformatted',
     &     access = 'sequential') !output significance (p-values)


      write (*, '(a10, a80)') '  >       ', file1
      write (*, '(a10, a80)') '  >       ', file2
      write (*, '(a10, a80)') '  >       ', file3
      write (*, '(a10, a80)') '  >       ', file4

      read(in2) ihead

      nx=ihead(5)
      ny=ihead(6)
 
      levfirst=ihead(2)
      
      cntlev=1  

      allocate(indat(nx,ny))

      read(in2) indat

 325  read(in2) ihead

      if(ihead(2).eq.levfirst)then

      goto 326
      else
      read(in2) indat
      cntlev=cntlev+1
      goto 325
      endif

 326  rewind(in2)

      nz = cntlev 
      allocate(nlev(nz))
      allocate(r(nx,ny),sig(nx,ny))
      allocate(sumx(nx,ny,nz),sumx2(nx,ny,nz),sumy(nx,ny,nz))  
      allocate(sumy2(nx,ny,nz),sumxy(nx,ny,nz),nt(nx,ny,nz))

      do i=1,nz
      read(in2) ihead
      nlev(i)=ihead(2)
      read(in2) indat
      print*,i,nlev(i)
      enddo


      rewind(in2)
    
      do i=1,nx
      do j=1,ny
      do k=1,nz
      sumx(i,j,k)=0.
      sumx2(i,j,k)=0.
      sumy(i,j,k)=0.
      sumy2(i,j,k)=0.
      sumxy(i,j,k)=0.
      nt(i,j,k)=0.
      enddo
      enddo
      enddo

      cntlev=1
 300  read (in2,end = 910) ihead

      
      read (in2) indat
      if(cntlev.eq.1)then
      read (in1,*) val
      endif 
      k=cntlev
c      print*,'check',val,indat(1,1)

      do i=1,nx
       do j=1,ny 
       if(indat(i,j).gt.-8e33)then
       nt(i,j,k)=nt(i,j,k)+ 1.0
       sumx2(i,j,k)=sumx2(i,j,k)+val*val
       sumx(i,j,k)=sumx(i,j,k)+val
       sumy(i,j,k)=sumy(i,j,k)+indat(i,j)
       sumy2(i,j,k)=sumy2(i,j,k)+indat(i,j)*indat(i,j)
       sumxy(i,j,k)=sumxy(i,j,k)+val*indat(i,j)      
       endif
       enddo
      enddo

      cntlev=cntlev+1
      if(cntlev.gt.nz)then
      cntlev=1
      endif

      goto 300

 910  continue

c      print*,'check 2',sumx(1,1),sumy(1,1),sumx2(1,1),sumy2(1,1)
c     + ,sumxy(1,1)

      do k=1,nz
      do i=1,nx
      do j=1,ny
      cnt=nt(i,j,k)
      if(cnt.gt.2)then
      dof=cnt-2
c      print*,i,j,sumx(i,j),sumy(i,j),sumx2(i,j),sumy2(i,j),sumxy(i,j)
      cov=sumxy(i,j,k)-1/cnt*sumx(i,j,k)*sumy(i,j,k)
      var=sqrt(sumx2(i,j,k)-1/cnt*sumx(i,j,k)*sumx(i,j,k))
      var1=sqrt(sumy2(i,j,k)-1/cnt*sumy(i,j,k)*sumy(i,j,k))

      r(i,j)=cov/var/var1
      tval=r(i,j)*sqrt(cnt-2)/sqrt(1-r(i,j)*r(i,j))
      if(tval.ge.0)then
      sig(i,j)=1-PROBT1(tval,tval,dof,2)
      else
      sig(i,j)=PROBT1(tval,tval,dof,2)-1
      endif
      else
      r(i,j)=-9e33
      sig(i,j)=-9e33
      endif

c      print*,'test...dof,r,tval,sig..... '
c     & ,dof,r(i,j),tval,sig(i,j)

      enddo
      enddo
      ihead(2)=nlev(k)

      write (out1) ihead
      write (out1) r
      write (out2) ihead
      write (out2) sig
      enddo
 


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


