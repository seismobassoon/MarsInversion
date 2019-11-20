subroutine invbyCG(nd,ata,atd,eps,x)
  ! modified from invbyCG Jun. 2009 Nobuaki Fuji
  !       
  !                       Nov. 2019 Nobuaki Fuji

  implicit none
  integer :: nd, ii
  real(kind(0d0)) :: ata(1:nd,1:nd)
  real(kind(0d0)) :: x(1:nd),r(1:nd),w(1:nd),z(1:nd),x0(1:nd),atd(1:nd)
  real(kind(0d0)) :: eps
  real(8) :: a, b, residual,initres,threshold

  x0 = 0
    
  r = atd - matmul(ata,x0)
  w = -r
  z = matmul(ata,w)
  a = dot_product(r,w) / dot_product(w,z)
  x = x0 +a*w
  b = 0

  initres=dot_product(atd,atd)
  threshold=eps*initrest
  
  do ii=1,nd
     r = r - a*z
     residual=dot_product(r,r)

     if(residual.le.threshold) exit
 
     b = dot_product(r,z)/dot_product (w,z)
     w = -r + b*w
     z = matmul(ata,w)
     
     !for pAAp calculation
     paap = dot_product(w,z)
     !end for pAAp calculation
     
     a = dot_product(r,w)/dot_product(w,z)
     x = x+a*w
     
  enddo

  

end subroutine invbyCG









subroutine pinput

  use parameters
  implicit none

  character(200) :: argv
  integer :: argc
  character(200) :: tmpfile,metafile
  character(200) :: obsfile
  character(200) :: dummy
  !integer, external :: getpid

  integer :: iloop,it

  call getarg(1,argv)
  metafile=argv
  call getarg(2,argv)
  workingDir=argv

  !write(tmpfile,"(Z5.5)") getpid()
  !tmpfile='tmpfileMarsInversion'//tmpfile
  tmpfile='tmpfileMarsInversion'
  open(unit=5,file=metafile,status='unknown')
  open(unit=1,file=tmpfile,status='unknown')
100 continue
  read(5,110) dummy
110 format(a200)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(5)

  open(unit=1,file=tmpfile,status='unknown')
  calculMode=0
  read(1,110) dummy
  if(dummy(1:4).eq.'test') calculMode=1
  read(1,*) dt
  read(1,*) tlen
  np=int(tlen/dt)
  read(1,*) npButterworth
  read(1,*) fmin
  read(1,*) fmax
  read(1,*) ntwin
  allocate(twin(1:4,1:ntwin))
  allocate(itwin(1:4,1:ntwin))
  do iloop=1,ntwin
     read(1,*) twin(1,iloop),twin(2,iloop),twin(3,iloop),twin(4,iloop)
  enddo
  itwin=int(twin/dt)
  read(1,*) tlenData
  npData = int(tlenData/dt)

  allocate(obsRaw(1:npData,1:3))
  allocate(obsFilt(1:npData,1:3))
  


  do iloop=1,3
     read(1,110) obsfile
     open(unit=10,file=obsfile,status='unknown')
     do it=1,npData
        read(10,*) obsRaw(it,iloop)
     enddo
     close(10)
  enddo
        

  read(1,*) nConfiguration
  allocate(filenames(nConfiguration))
  do iloop=1,nConfiguration
     read(1,110) filenames(iloop)
  enddo



  
  close(1)

end subroutine pinput


subroutine bwfilt(x,y,dt,n,irek,norder,f1,f2)
  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points
  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering
  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder<0: no starplots of transfer function and impulse response
  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter
  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter
  implicit none
  real(kind(0d0)), dimension(1) :: x,y
  real(kind(0d0)), dimension(10) :: a,b1,b2
  real(kind(0d0)) :: dt,f1,f2
  integer :: iunit,npoles,norder,irek,n,lx
  
  iunit = 3
  if (norder.ne.0) then
    npoles = iabs(norder)
    !determination of filter coefficients
    call bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    if (norder.ge.0) then
      !plot of transfer function and impuulse response
      lx = 100
      !filtering
    endif
  endif
  if (n.ne.0) then
    call rekurs(x,y,n,a,b1,b2,npoles,irek)
  endif

  return

end subroutine bwfilt

subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag.ne.0: forward and backward filtering
  implicit none
  real(kind(0d0)), dimension(10) :: z,z1,z2,a,b1,b2
  real(kind(0d0)) :: x1,x2
  integer :: ndat,npoles,iflag,n,i
  real(kind(0d0)), dimension(ndat) :: x,y
  
  !forward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = 1,ndat
    z(1) = a(1)*(x(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = x(n)
    do i = 1,npoles
      z2(i) = z1(i)
      z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo
  if (iflag.eq.0) then
    return
  endif
  !backward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = ndat,1,-1
    z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = y(n)
    do i = 1,npoles
       z2(i) = z1(i)
       z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo

  return

end subroutine rekurs

subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter
  implicit none
  real(kind(0d0)), dimension(10) :: a,b1,b2
  complex(kind(0d0)), dimension(20) :: s
  complex(kind(0d0)) :: t1,t2,p
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: f1,f2,dt,d2,w0,w1,w2,ssum,sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles
  
  if (npoles.gt.10) stop 'npoles greater than 10: STOP'
  d2 = 2.d0/dt
  w1 = d2*tan(2.d0*pi*f1/d2)
  w2 = d2*tan(2.d0*pi*f2/d2)
  w0 = 0.5*(w2-w1)
  i = 1
  npol2 = npoles/2+1
  do n = 1,npoles
    p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
    t1 = p*cmplx(w0,0.d0)
    t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
    s(i) = t1+t2
    s(i+1) = t1-t2
    i = i+2
  enddo 
  do n = 1,npoles
    ssum = 2*real(s(n))
    sprod = dble(s(n)*conjg(s(n)))
    fact1 = d2*d2-d2*ssum+sprod
    fact2 = 2.d0*(sprod-d2*d2)
    fact3 = d2*d2+d2*ssum+sprod
    a(n) = 2.d0*d2*w0/fact1
    b1(n) = fact2/fact1
    b2(n) = fact3/fact1
  enddo

  return

end subroutine bpcoeff
