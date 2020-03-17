      program srv2grads2
c
c Program removes headers from service files to make them 
c readable by grads in binary format
c
c 
c
      parameter (iu1 = 10,io2=12)
      integer ihead (8),ilev
      character*256 file1, file2, file3
      character cform * 16
      real fi1, undef
      real, dimension (:,:), allocatable :: xs,ys
      parameter (undef = 9e12)

      namelist /param/file1,file2,file3
      data nlon/-99/, nlat/-99/, file1 /'  '/, file2 /'  '/

      read (5, param, end = 101)
 101  continue
      if (file1 .eq. ' ') stop 'No input file specified'
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

      open (io2, file = file2, form = 'unformatted',
     &     access = 'sequential')



      
      read (iu1) ihead
      print*,ihead,'testing headers'
      nlon=ihead (5)
      nlat=ihead (6)
      rewind (iu1)
      allocate (xs (nlon,  nlat))

 100  read (iu1,end=901) ihead
      read (iu1)((xs (jl, jk), jl = 1, nlon), 
     &     jk = 1,  nlat)
      ilev=ihead(2)

      write (*,*) 'removing headers for level',ilev,'...',nlon,nlat
      write (io2) xs
      go to 100
 900  stop 'error field is emty'
 901  print *, 'regular end'
      
      end




