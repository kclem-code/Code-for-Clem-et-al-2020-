      program get_ratio
c
c Program calculates the ratio of two files, only if the values
c are at a specified threshold of a third file.  
c Returns the ratio (first divided by second)
c 
c
      parameter (iu1 = 10,iu2=12,iu3=14,io1=16)
      integer ihead (8),idat
      character*256 file1, file2, file3,file4
      character cform * 16
      real fi1, undef,thold
      real, dimension (:,:), allocatable :: indat1,indat2
      real, dimension (:,:), allocatable :: val,ratio
      parameter (undef = 9e12)

      namelist /param/file1,file2,file3,file4,thold
      data nlon/-99/, nlat/-99/, file1 /'  '/, file2 /'  '/,
     + file3 /' '/,file4 /' '/,thold/-999/

      read (5, param, end = 101)
 101  continue
      if (file1 .eq. ' ') stop 'No input file1 specified'
      if (file2 .eq. ' ') stop 'No input file2 specified'
      if (file3 .eq. ' ') stop 'No input file3 specified'
      if (file4 .eq. ' ') stop 'No output file specified'
      if(thold.eq.-999) stop 'no threshold value specified'
c      lfile1 = index (file1, ' ') - 1
c      if (lfile1 .eq. -1) lfile1 = mlfile
c      if (file2 .eq. ' ')
c     &     file2 = file1 (1:min (lfile1, mlfile - 4)) // '.dat'
c      lfile2 = index (file2, ' ') - 1
c      if (lfile2 .eq. -1) lfile2 = mlfile

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
     &     access = 'sequential', status='old')
      open (iu3, file = file3, form = 'unformatted',
     &     access = 'sequential', status='old')


      open (io1, file = file4, form = 'unformatted',
     &     access = 'sequential')


      
      read (iu1) ihead
      nlon=ihead (5)
      nlat=ihead (6)
      rewind (iu1)
      allocate (indat1(nlon,nlat),indat2(nlon,nlat))
      allocate (ratio(nlon,nlat),val(nlon,nlat))

 100  read (iu1,end=901) ihead
      read (iu1) indat1

      read(iu2) ihead
      read(iu2) indat2
      read(iu3) ihead
      read(iu3) val

      do i=1,nlon
      do j=1,nlat
      if(abs(val(i,j)).ge.thold)then
      ratio(i,j)=indat1(i,j)/indat2(i,j)
      else
      ratio(i,j)=-9e33
      endif
      enddo
      enddo
      
      write (io1) ihead
      write (io1) ratio
      go to 100
 900  stop 'error field is emty'
 901  print *, 'regular end'
      
      end




