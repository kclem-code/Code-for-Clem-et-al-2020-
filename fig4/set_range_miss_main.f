      program srv2grads
c
c Program uses statistical significance field to set grid cells 
c in the data field below a threshold (p value) to missing (-9e33).
c
c 
c
      parameter (iu1 = 10,iu2=12, io1=14, io2=16)
      parameter (iu3 = 18, iu4=20)
      integer ihead (8),idat
      character*256 file1, file2, file3,file4,file5,file6
      character cform * 16
      real fi1, undef,thold,miss
      real, dimension (:,:), allocatable :: in1,in2,out1,out2
      real, dimension (:,:), allocatable :: sig1,sig2
      parameter (undef = 9e12)

      namelist /param/file1,file2,file3,file4,file5,file6,thold
      data nlon/-99/, nlat/-99/, file1 /'  '/, file2 /'  '/,
     +file3/'   '/,file4/'    '/,thold/-99./

      read (5, param, end = 101)
 101  continue
      if (file1 .eq. ' ') stop 'No input file specified'
c      lfile1 = index (file1, ' ') - 1
c      if (lfile1 .eq. -1) lfile1 = mlfile
c      if (file2 .eq. ' ')
c     &     file2 = file1 (1:min (lfile1, mlfile - 4)) // '.dat'
c      lfile2 = index (file2, ' ') - 1
c      if (lfile2 .eq. -1) lfile2 = mlfile
      if(thold.eq.-99) stop 'No filter threshold specified'

c      write (*, *)
c      write (cform, '(a13, i2, a1)') '(x, a10, x, a', lfile1, ')'
c      write (*, cform) '  >       ', file1 (1:lfile1)
c      write (cform, '(a13, i2, a1)') '(x, a10, x, a', lfile2, ')'
c      write (*, cform) '  >       ', file2 (1:lfile2)
c      write (*, *)
c      write (*, '(x, a11)') 'Parameters:'
c      write (*, *) 'nlon      ', nlon
c      write (*, *) 'nlat      ', nlat
c      write (*,*) 'hallo'
      open (iu1, file = file1, form = 'unformatted',
     &     access = 'sequential', status='old')

      open (iu2, file = file2, form = 'unformatted',
     &     access = 'sequential',status='old')

      open (iu3, file = file3, form = 'unformatted',
     &     access = 'sequential', status='old')

      open (iu4, file = file4, form = 'unformatted',
     &     access = 'sequential',status='old')

      open (io1, file = file5, form = 'unformatted',
     &     access = 'sequential')

      open (io2, file = file6, form = 'unformatted',
     &     access = 'sequential')

      print*,'file 1...',file1
      print*,'file 2...',file2
      print*,'file 3...',file3
      print*,'file 4...',file4
      print*,'file 5...',file5
      print*,'file 6...',file6


      miss=-9e33
      
      read (iu1) ihead
      nlon=ihead (5)
      nlat=ihead (6)
      rewind (iu1)
      allocate (in1 (nlon,  nlat),in2(nlon,nlat),out1(nlon,nlat))
      allocate(out2(nlon,nlat),sig1(nlon,nlat),sig2(nlon,nlat))     

 100  read (iu1,end=901) ihead
      read (iu2) ihead

      read (iu1) in1

      read (iu2) in2

      read (iu3) ihead
      read (iu4) ihead
      read (iu3) sig1
      read (iu4) sig2

      do i=1,nlon
      do j=1,nlat
      if(abs(sig1(i,j)).ge.thold.or.abs(sig2(i,j)).ge.thold)then
      out1(i,j)=in1(i,j)
      out2(i,j)=in2(i,j)
      else
      out1(i,j)=miss
      out2(i,j)=miss
      endif
      enddo
      enddo

      write (io1) ihead
      write (io2) ihead
      write (io1) out1
      write (io2) out2

      go to 100
 900  stop 'error field is emty'
 901  print *, 'regular end'
      
      end




