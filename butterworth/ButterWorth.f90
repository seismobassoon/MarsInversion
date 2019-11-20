program mtInversion

!!!!!!
!
!   MT inversion for InSight data 
!
!                     Nov. 2019 Nobuaki Fuji @ IPGP
!  
!!!!!!


  use parameters
  implicit none

  integer :: mtcomp,jmtcomp
  integer ::icomp,iWindow,it
  integer :: iMovingWindow,iConfiguration
  real(kind(0d0)), allocatable :: taper(:)
  real(kind(0d0)), allocatable :: tmparray(:,:,:)
  real(kind(0d0)), allocatable :: GreenArray(:,:,:)
  real(kind(0d0)), allocatable :: obsArray(:,:)
  real(kind(0d0)), allocatable :: filtbefore(:),filtafter(:)
  real(kind(0d0)) :: xfwin
  real(kind(0d0)) :: ata(1:nmt,1:nmt),atd(1:nmt)
  real(kind(0d0)) :: mtInverted(1:nmt,1:npData-np+1,1:nConfiguration)
  real(kind(0d0)) :: misfitTaper(1:nmt,1:npData-np+1,1:nConfiguration)
  real(kind(0d0)) :: misfitRaw(1:nmt,1:npData-np+1,1:nConfiguration)
  character(200) :: synfile
  real(kind(0d0)) :: dummyFloat
  

110 format(a200)
  ! making taper function

  call pinput

  !print *, np,ntwin
  !print *, itwin(1:4,1)

  allocate(taper(1:np))

  taper=0.d0
  do iWindow=1,ntwin
     do it=itwin(1,iWindow),itwin(4,iWindow)
        if((it.gt.itwin(1,iWindow)).and.(it.lt.itwin(2,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwin(1,iWindow))/dble(itwin(2,iWindow)-itwin(1,iWindow)))
           taper(it)=xfwin*xfwin
        elseif((it.ge.itwin(2,iWindow)).and.(it.le.itwin(3,iWindow))) then
           taper(it)=1.d0
        elseif((it.gt.itwin(3,iWindow)).and.(it.lt.itwin(4,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwin(4,iWindow))/dble(itwin(3,iWindow)-itwin(4,iWindow)))
           taper(it)=xfwin*xfwin
        endif
     enddo
  enddo
   
  ! this is only for the use of confirmation
  !do it=1,np
  !   write(13,*) taper(it)
  !enddo
  !stop
           



  
  ! we apply the butterworth filter to the Green's functions


  allocate(tmparray(1:np,1:3,1:nmt))
  allocate(GreenArray(1:np,1:3,1:nmt))
  allocate(obsArray(1:np,1:3))
  allocate(filtbefore(1:np),filtafter(1:np))



  ! Grande boucle pour chaque configuration

  do iConfiguration=1,nConfiguration


     ! Lecture des fonctions de Greens a partir des fichiers dans le filenames(iConfiguration)

     open(unit=1,file=filenames(iConfiguration),status='unknown')
     do mtcomp=1,nmt
        read(1,110) synfile
        open(unit=10,file=synfile,status='unknown')
        do it=1,np
           read(10,*) dummyFloat,tmparray(it,1,mtcomp),tmparray(it,2,mtcomp),tmparray(it,3,mtcomp)
        enddo
        close(10)
     enddo
     close(1)
     


     ! Here we first filter Green's function as a whole and taper them
     
     do mtcomp=1,nmt
        do icomp=1,3
           filtbefore(1:np)=tmparray(1:np,icomp,mtcomp)
           call bwfilt(filtbefore,filtafter,dt,np,1,npButterworth,fmin,fmax)
           tmparray(1:np,icomp,mtcomp)=filtafter(1:np)
           GreenArray(1:np,icomp,mtcomp)=filtafter(1:np)*taper(1:np)
        enddo
     enddo
     

     open(11,file="synZ"//iConfiguration//".txt",status='unknown')
     open(12,file="synN"//iConfiguration//".txt",status='unknown')
     open(13,file="synE"//iConfiguration//".txt",status='unknown')
     do it=1,npData
        write(11,*) dble(it)*dt,obsRaw(it,1),obsFilt(it,1)
        write(12,*) dble(it)*dt,obsRaw(it,2),obsFilt(it,2)
        write(13,*) dble(it)*dt,obsRaw(it,3),obsFilt(it,3)
     enddo
     close(11)
     close(12)
     close(13)



     
     ! Construct AtA
     
     
     ata=0.d0
     do mtcomp=1,nmt
        do jmtcomp=1,mtcomp
           ata(mtcomp,jmtcomp)=sum(GreenArray(1:np,1:icomp,mtcomp)*GreenArray(1:np,1:icomp,jmtcomp))
        enddo
     enddo
     
     ! AtA is symmetric
     
     do mtcomp=1,nmt
        do jmtcomp=mtcomp,nmt
           ata(mtcomp,jmtcomp)=ata(jmtcomp,mtcomp)
        enddo
     enddo
     
     
     
     
     do mtcomp=1,nmt
        do jmtcomp=1,nmt
           print *, mtcomp,jmtcomp,ata(mtcomp,jmtcomp)
        enddo
     enddo
     
     
!!!!!!! NOW WE START MT INVERSION FOR EACH TIME iMovingWindow
     
     
     
     
     ! For the observed data we first filter the whole signal of a lenght of npData 
     
     do icomp=1,3
        call bwfilt(obsRaw(1:npData,icomp),obsFilt(1:npData,icomp),dt,npData,1,npButterworth,fmin,fmax)
     enddo
     
     ! Raw data and filtered data are written as fort.11-13 for references
     

     open(11,file="obsZtreated.txt",status='unknown')
     open(12,file="obsNtreated.txt",status='unknown')
     open(13,file="obsEtreated.txt",status='unknown')
     do it=1,npData
        write(11,*) dble(it)*dt,obsRaw(it,1),obsFilt(it,1)
        write(12,*) dble(it)*dt,obsRaw(it,2),obsFilt(it,2)
        write(13,*) dble(it)*dt,obsRaw(it,3),obsFilt(it,3)
     enddo
     close(11)
     close(12)
     close(13)

     
     
     ! Here we taper the observed function for each moving window of np
     
     do iMovingWindow=1,(npData-np+1)
        do icomp=1,3
           obsArray(1:np,icomp)=obsFilt(iMovingWindow:iMovingWindow+np-1,icomp)*taper(1:np)
        enddo
        atd=0.d0
        ! Atd construction
        
        do 
           



        
     enddo
     
     
    
     
  enddo
     
 
  
     
  

  
  
  
end program mtInversion
