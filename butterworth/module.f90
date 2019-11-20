module parameters 
  implicit none 
  
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: npButterworth
  real(kind(0d0)) :: fmin,fmax
  real(kind(0d0)) :: eps ! tolerance during MT inversion 
  integer, parameter :: nmt = 6
  integer :: np ! integer length of synthetics (to be considered)
  integer :: npData ! integer lenght of observed 
  integer :: ntStep
  real(kind(0d0)) :: dt, tlen, tlenData 
  real(kind(0d0)), allocatable :: twin(:,:)
  integer, allocatable :: itwin(:,:)
  integer :: ntwin
  character(200), allocatable :: filenames(:)
  real(kind(0d0)), allocatable :: obsRaw(:,:), obsFilt(:,:)
  integer :: calculMode ! 0=normal; 1=filter and stop 
  character(200) :: workingDir
  character(200) :: resultDir 
  character(200) :: inversionName
  integer :: nConfiguration ! The number of configurations (source location in 3D space, model)


end module parameters
