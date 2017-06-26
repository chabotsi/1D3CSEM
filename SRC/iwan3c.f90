subroutine Iwan3c(Tdomain,Enonli,it,TimeP,Evisco)

use sdomain
use nonlinear 
use stimeparam
use viscopara

implicit none
   
type(nonli), intent(INOUT),  target      ::  Enonli
type(domain),intent(INOUT),  target      ::  Tdomain
integer,     intent (IN)                 ::  it
type (time), intent (IN),    target      ::  TimeP
type(visco) , intent (INOUT), target     ::  Evisco


integer                    :: n 
integer                    :: i 
integer                    :: nummat 
real                       :: Emax 
real                       :: Kmax 
real                       :: Gmax 
real                       :: R     (Enonli%Nspr)
real                       :: CNinv (Enonli%Nspr-1)
real                       :: dsigma(6)
real                       :: deps  (6)
real                       :: Ni


integer, pointer                :: ngllx 
integer, pointer                :: nel 

real, dimension(:), pointer     :: Gm0
real, dimension(:), pointer     :: dplastic
real, dimension(:), pointer     :: S1 
real, dimension(:), pointer     :: F2

real, dimension(:,:), pointer   :: Stress
real, dimension(:,:), pointer   :: Sa1 



nel     =>    Tdomain%nelem

do n = 0,nel-1

  ngllx     =>  Tdomain%specel(n)%ngllx 
  nummat    =   Tdomain%specel(n)%Numdomain-1 
  Emax      =   Tdomain%specel(n)%Emodulus
  Gmax      =   Tdomain%specel(n)%Gmodulus
  Kmax      =   Tdomain%specel(n)%Kmodulus

  Stress    =>  Tdomain%specel(n)%Stress

  Gm0       =>  Enonli%specel(n)%Gm0
  Ni        =   Tdomain%specel(n)%Ni


  ! gref given model
  if (.NOT. Tdomain%Epressmod) then
    Gmax = Tdomain%specel(n)%Gmodulus
    Emax = Tdomain%specel(n)%Emodulus
    Kmax = Tdomain%specel(n)%Kmodulus

    ! Viscosity in
    if (Tdomain%rheovisco) then
      Gmax  = Evisco%specel(n)%visM(0)
      Emax  = 2.0* Gmax* (1.0+ Ni)
      Kmax  = Emax/ (3.0*(1.0- 2.0* Ni)) 
    endif  
  endif

  do i = 0, ngllx-1  

    dplastic => Tdomain%specel(n)%dplastic(i,:) 
    S1       => Enonli%specel(n)%hexaS1(i,:)
    Sa1      => Enonli%specel(n)%plastSa1(i,:,:)
    F2       => Enonli%specel(n)%plastF2(i,:)


    deps     = Tdomain%specel(n)%deps(i,:)


    if (.not. Tdomain%Epressmod) then
      R     = Enonli%R(:,nummat+1)
      CNinv = Enonli%CNinv(:,nummat+1) 
    else
      R     = Enonli%specel(n)%RSAT(:,i)
      CNinv = Enonli%specel(n)%CNinvSAT(:,i)

      Gmax  = Enonli%specel(n)%Gmod2(i)
      Emax  = 2.0* Gm0(i)* (1.0+ Ni)
      Kmax  = Emax/ (3.0*(1.0- 2.0* Ni)) 
    endif

    if (Tdomain%Epressmod .AND. .not. Tdomain%Sub_domain(nummat)%Eeffective) then
      Gmax  = Gm0(i)
      Emax  = 2.0* Gmax* (1.0+ Ni)
      Kmax  = Emax/ (3.0*(1.0- 2.0* Ni)) 
    endif


    call neoiwan(it,Enonli%Nspr,Emax,Gmax,Kmax,deps,dsigma,Sa1,F2,&
                    S1,R,CNinv,dplastic, Ni,Gm0(i),Enonli%specel(n)%aktifsur(i))


    ! Updating global variables
    if (.NOT. Tdomain%specel(n)%PML)    &
    Stress(i,:) = Stress(i,:)+ dsigma 


!     if (n==2 .AND. i==0)&
!     write(99,*)  Gmax, Emax* Ni/(1.0+ Ni)/(1- 2.0*Ni), Kmax


  enddo
enddo

end subroutine Iwan3c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine neoiwan(it,Nspr,Emax,Gmax,Kmax,deps,dSigma,Sa1,F2,S1,R,CNinv,&
                      dplastic,Ni,Ginitial,aktifsur)

implicit none


integer,  intent(IN)       ::  Nspr 
integer,  intent(IN)       ::  it 
real,     intent(IN)       ::  Emax 
real,     intent(IN)       ::  Gmax
real,     intent(IN)       ::  Kmax
real,     intent(INOUT)    ::  Sa1(Nspr,6)
real,     intent(INOUT)    ::  F2    (Nspr)
real,     intent(INOUT)    ::  R     (Nspr)
real,     intent(INOUT)    ::  CNinv (Nspr-1)
real,     intent(INOUT)    ::  dSigma(6)
real,     intent(INOUT)    ::  S1    (6)
real,     intent(IN)       ::  deps  (6)
real,     intent(INOUT)    ::  dplastic(6)
real,     intent(IN)       ::  Ni
real,     intent(IN)       ::  Ginitial
integer,  intent(INOUT)    ::  aktifsur


real :: depsm, dsigm
real :: S2(6), dS(6), F1(Nspr)

real :: dF(Nspr)
! real :: dF

real      ::  ss  (6)
real      ::  Ln(Nspr)
real      ::  de  (6)
real      ::  Kbm (6,6)
real      ::  Sm  (6,6)
real      ::  Af  (6,6)
real      ::  m1  (6,6)
real      ::  m2  (6,6)
real      ::  m3  (6,6)
real*8    ::  Esd (6,6), Ed(6,6)
integer   ::  j 
integer   ::  m 
integer   ::  k 
integer   ::  errorflag
real      ::  lambda
real      ::  depsvol
INTEGER   ::  aktif, start



depsm   = (deps(1) + deps(2)+deps(6))/ 3.0    
dsigm   = Kmax* depsm* 3.0


lambda  = Emax* Ni/(1.0+ Ni)/(1- 2.0*Ni)
depsvol = deps(1)+ deps(2)+ deps(6)


! Emax and Kmax are used only for principal directions 
! and held independent of saturation

de      = 0.0
if (it == 1) then

  call IWAN_elastic(Ginitial,Gmax,lambda,deps,dsigma) 

  dS    = dsigma
  dS(1) = dsigma(1)- dsigm
  dS(2) = dsigma(2)- dsigm
  dS(6) = dsigma(6)- dsigm   

  S1 = S1+ dS 

  Sa1 = 0.0
  dplastic  = 0.0
  aktifsur  = 0

  return
endif


de(1) = deps(1)- depsm
de(2) = deps(2)- depsm
de(3) = deps(3)/ 2.0
de(4) = deps(4)/ 2.0
de(5) = deps(5)/ 2.0
de(6) = deps(6)- depsm


if (aktifsur == 0) then
  call IWAN_surface(S1,Sa1(1,:),F2(1))
  call IWAN_dsurface(S1,Sa1(1,:),de,dF(1))
endif



aktif = 0
if (aktifsur .GT. 0) then
do j = 1,aktifsur 

  do m = 1,6
    Sa1(j,m) = S1(m)- R(j)/ sqrt(F2(j))*(S1(m)-Sa1(j,m))
  enddo


  call IWAN_surface(S1,Sa1(j,:),F2(j))
  call IWAN_dsurface(S1,Sa1(j,:),de,dF(j))


  if ( (dF(j) .GE. 0.0)  .AND. (F2(j) .GE. R(j)**2) ) &
  aktif = aktif+ 1
enddo
endif


if ( F2(1) .lt. R(1)**2  .AND. dF(1) .ge. 0.0)   then

  call IWAN_elastic(Ginitial,Gmax,lambda,deps,dsigma) 

  dS    = dsigma
  dS(1) = dsigma(1)- dsigm
  dS(2) = dsigma(2)- dsigm
  dS(6) = dsigma(6)- dsigm   

  S1 = S1+ dS 
  dplastic = 0.0

  return    
endif



! Ed computation
Ed    = 0.0
start = 1

Ed(1,1) = 0.5/Gmax
Ed(2,2) = 0.5/Gmax
Ed(3,3) = 0.5/Gmax
Ed(4,4) = 0.5/Gmax
Ed(5,5) = 0.5/Gmax
Ed(6,6) = 0.5/Gmax

if (aktif .GT. 0) &
call Ematris(start,Nspr,Ed,CNinv,S1,Sa1,F2,aktif)


do j = aktif+1, Nspr-1

  call IWAN_surface(S1,Sa1(j,:),F2(j))
  call IWAN_dsurface(S1,Sa1(j,:),de,dF(j))


  if ( (dF(j) .GE. 0.0)  .AND. (F2(j) .GE. R(j)**2) ) then 
    aktif = aktif+ 1
    start = aktif

    call Ematris(start,Nspr,Ed,CNinv,S1,Sa1,F2,aktif)
  else 
    EXIT
  endif
enddo


Esd   = 0.0
call FINDInv(Ed, Esd, 6, errorflag)
dS    = MATMUL(sngl(Esd),de)
S1 = S1+ dS 

aktifsur = max(aktif,1)



! Total Stresses
dsigma(1) = dS(1)+ dsigm
dsigma(2) = dS(2)+ dsigm
dsigma(3) = dS(3)
dsigma(4) = dS(4)
dsigma(5) = dS(5)
dsigma(6) = dS(6)+ dsigm


! Plastic strain increment array (gamma for shear not epsilon)
dplastic    = 0.0
dplastic(4) = deps(4)- dsigma(4)/ Gmax
dplastic(5) = deps(5)- dsigma(5)/ Gmax
dplastic(6) = deps(6)- dsigma(6)/ (lambda+ 2.0*Ginitial)


end subroutine neoiwan
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Ematris(start,Nspr,Ed,CNinv,S1,Sa1,F1,aktif)

integer, intent(IN)     :: start
integer, intent(IN)     :: Nspr
real*8,  INTENT(INOUT)  :: Ed(6,6)
real,    INTENT(IN)     :: CNinv(Nspr-1)
real,    INTENT(IN)     :: S1(6)
real,    INTENT(IN)     :: Sa1(Nspr,6)
real,    INTENT(IN)     :: F1 (Nspr)
integer, intent(IN)     :: aktif


integer :: j,m,k
real    :: ss(6)


ss(1) = 1.0
ss(2) = 1.0
ss(3) = 2.0
ss(4) = 2.0
ss(5) = 2.0
ss(6) = 1.0

j = start
do while (j .lt. aktif+1)
  do m = 1,6
    do k = 1,6    
      Ed(m,k) = Ed(m,k)+ CNinv(j)* ss(k)* (S1(m)-Sa1(j,m))&
                  *(S1(k)-Sa1(j,k))/ (2.0* F1(j))
    enddo
  enddo
  j = j+1
enddo

end subroutine Ematris
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Iwan surface(s) computation 

subroutine IWAN_surface(S,Sa,F)  
  
  real, intent(in) :: S(6)
  real, intent(in) :: Sa(6)
  real, intent(inout) :: F
 
  F = 0.0
  F = 0.5* ((S(1)-Sa(1))**2+ (S(2)-Sa(2))**2+ &
                  2.0* (S(3)-Sa(3))**2 + &
                  2.0* (S(4)-Sa(4))**2 + &
                  2.0* (S(5)-Sa(5))**2 &
                  + (S(6)-Sa(6))**2)

end subroutine IWAN_surface  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Iwan surface(s) movement computation 

subroutine IWAN_dsurface(S,Sa,de,dF) 
  
  real, intent(in) :: S(6)
  real, intent(in) :: Sa(6)
  real, intent(in) :: de(6)
  real, intent(inout) :: dF
 
  dF = 0.0
  dF = 0.5* ((S(1)-Sa(1))*de(1)+ (S(2)-Sa(2))*de(2)+ &
                  2d0* (S(3)-Sa(3))*de(3) + &
                  2d0* (S(4)-Sa(4))*de(4) + &
                  2d0* (S(5)-Sa(5))*de(5) &
                  + (S(6)-Sa(6))*de(6))

end subroutine IWAN_dsurface  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

subroutine IWAN_elastic (Ginitial,Gmax,lambda,deps,dsigma)  

  real, intent(in) :: Ginitial, Gmax, lambda
  real, intent(in) :: deps(6)
  real, intent(out):: dsigma(6)

  real :: depsvol


  depsvol = deps(1)+ deps(2)+ deps(6)

  dsigma(1) = 2.0*Ginitial*deps(1)+ lambda*depsvol
  dsigma(2) = 2.0*Ginitial*deps(2)+ lambda*depsvol
  dsigma(3) = Gmax * deps(3)
  dsigma(4) = Gmax * deps(4)
  dsigma(5) = Gmax * deps(5)
  dsigma(6) = 2.0*Ginitial*deps(6)+ lambda*depsvol


end subroutine IWAN_elastic




!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !Subroutine to find the inverse of a square matrix
        !Author : Louisda16th a.k.a Ashwith J. Rego
        !Reference : Algorithm has been well explained in:
        !http://math.uww.edu/~mcfarlat/inverse.htm
        !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

      SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
      IMPLICIT NONE
      !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k, l
        REAL*8 :: m
        REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix

!Augment input matrix with an identity matrix
        DO i = 1, n
      DO j = 1, 2*n
         IF (j <= n ) THEN
            augmatrix(i,j) = matrix(i,j)
         ELSE IF ((i+n) == j) THEN
            augmatrix(i,j) = 1.0D0
         Else
            augmatrix(i,j) = 0.0D0
         ENDIF
      END DO
        END DO

!Reduce augmented matrix to upper traingular form
        DO k =1, n-1
      IF (augmatrix(k,k) == 0.0D0) THEN
         FLAG = .FALSE.
         DO i = k+1, n
            IF (augmatrix(i,k) /= 0.0D0) THEN
               DO j = 1,2*n
                  augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
               END DO
               FLAG = .TRUE.
               EXIT
            ENDIF
            IF (FLAG .EQV. .FALSE.) THEN
               PRINT*, "Matrix is non - invertible"
               inverse = 0.0D0
               errorflag = -1
               return
            ENDIF
         END DO
      ENDIF
      DO j = k+1, n
         m = augmatrix(j,k)/augmatrix(k,k)
         DO i = k, 2*n
            augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
         END DO
      END DO
        END DO

!Test for invertibility
        DO i = 1, n
      IF (augmatrix(i,i) == 0.0D0) THEN
         PRINT*, "Matrix is non - invertible"
         inverse = 0.0D0
         errorflag = -1
         return
      ENDIF
        END DO

!Make diagonal elements as 1
        DO i = 1 , n
      m = augmatrix(i,i)
                DO j = i , (2 * n)
         augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO

!Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
         m = augmatrix(i,k+1)
         DO j = k, (2*n)
            augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
         END DO
      END DO
        END DO

!store answer
        DO i =1, n
      DO j = 1, n
         inverse(i,j) = augmatrix(i,j+n)
      END DO
        END DO
        errorflag = 0
        END SUBROUTINE FINDinv
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       