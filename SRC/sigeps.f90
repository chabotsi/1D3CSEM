subroutine sigeps(Tdomain,TimeP,Source,it,Evisco,Ewave,posi,Enonli)

use sdomain 
use stimeparam
use ssource
use viscopara
use inciwave
use nonlinear


implicit none

type (domain),  intent (INOUT), target    ::  Tdomain
type (time),    intent (IN),    target    ::  TimeP
type (sour),    intent (IN)               ::  Source
integer, intent (IN)                      ::  it
type(visco), intent(INOUT),     target    ::  Evisco
type(inci), intent(INOUT),      target    ::  Ewave
integer, intent (IN)                      ::  posi
type(nonli),    intent(INOUT),  target    ::  Enonli


integer                   ::  n
integer                   ::  i_domain
integer                   ::  j 
integer                   ::  xx
integer                   ::  tfactor
integer                   ::  k 
integer                   ::  nummat
integer                   ::  interpos

real                      ::  Dt 
real                      ::  coef1 
real                      ::  coef2
real                      ::  interpo
real                      ::  yenihiz
real                      ::  veloc1
real                      ::  veloc2

integer, pointer                  ::  nel 
integer, pointer                  ::  npt 
integer, pointer                  ::  ngllx
integer, pointer                  ::  ndom 

real, dimension(:), allocatable   ::  ni 
real, dimension(:), allocatable   ::  lambda
real, dimension(:), allocatable   ::  mu 

real, dimension(:), allocatable   ::  s1
real, dimension(:), allocatable   ::  Vxloc 
real, dimension(:), allocatable   ::  elVxloc 
real, dimension(:), allocatable   ::  els1
real, dimension(:), allocatable   ::  els2
real, dimension(:), allocatable   ::  s2 
real, dimension(:), allocatable   ::  zipto
real, dimension(:,:), allocatable ::  HT 

real, dimension(:), pointer       ::  DumpSx3
real, dimension(:), pointer       ::  DumpSx3YZ
real, dimension(:), pointer       ::  DumpSx3ZZ
real, dimension(:), pointer       ::  DumpSx3VOL


real, dimension(:,:), pointer     ::  zipt
real, dimension(:,:), pointer     ::  ziptYZ
real, dimension(:,:), pointer     ::  ziptZZ
real, dimension(:,:), pointer     ::  ziptVOL


real, dimension(:,:), pointer     ::  LocStress 
real, dimension(:,:), pointer     ::  Veloc 
real, dimension(:,:), pointer     ::  Accel
real, dimension(:,:), pointer     ::  gamma1
real, dimension(:,:), pointer     ::  gamma2



nel   =>   Tdomain%nelem

! STRAIN and/or STRESS CALCULATION

do n = 0,nel-1  
    
  ngllx       =>  Tdomain%specel(n)%ngllx
  LocStress   =>  Tdomain%specel(n)%Stress 
  Veloc       =>  Tdomain%specel(n)%Veloc
  Accel       =>  Tdomain%specel(n)%Accel 
  gamma1      =>  Tdomain%specel(n)%gamma1
  gamma2      =>  Tdomain%specel(n)%gamma2

  Dt      =   Tdomain%specel(n)%dt
  coef1   =   TimeP%beta/TimeP%gamm
  coef2   =   1.0/TimeP%gamm

  ! Updating previous strain
  gamma1  =   gamma2			
  gamma2  =   0.0

  
  do xx = 1,3

    allocate(HT(0:ngllx-1,0:ngllx-1))
    allocate(s1(0:ngllx-1))
    allocate(s2(0:ngllx-1)) 
    allocate(Vxloc(0:ngllx-1))
    allocate(elVxloc(0:ngllx-1))
    allocate(els1(0:ngllx-1))
    allocate(els2(0:ngllx-1))  

    allocate(ni(0:ngllx-1))
    allocate(lambda(0:ngllx-1))
    allocate(mu(0:ngllx-1))

    ni = Tdomain%specel(n)%Ni 
    mu = Tdomain%specel(n)%Mu 
    lambda = -2.0* Tdomain%specel(n)%Mu+ Tdomain%specel(n)%lambda2mu


    s1      = 0.0
    s2      = 0.0 
    HT      = 0.0 
    Vxloc   = 0.0 
    elVxloc = 0.0
    els1    = 0.0 
    els2    = 0.0
 
    i_domain  =   Tdomain%specel(n)%NumDomain- 1
    HT        =   Tdomain%Sub_domain(i_domain)%hTprime               
    Vxloc     =   Veloc(:,xx)+ Dt*(0.5- coef1)*(1.0- coef2)* Accel(:,xx)



! ABSORBING LAYER (PML) USED WITH INCIDENT WAVE (of 1 element of PML)
! input velocity multiplied by 2 (shared by upside and downside)
  if (Tdomain%specel(nel-1)%PML) then
  Ewave%sposi = 1
! if ( Tdomain%specel(n)%PML ) then
      if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
      ! Input signal implementation & interpolation                     
      elVxloc    =   0.0
      interpos   =   MOD((it- 1),Ewave%tfactor)
      interpo    =   float(interpos)/ float(Ewave%tfactor) 

      if (n == nel- Ewave%sposi-1  .AND. TimeP%rtime .le. Ewave%timeup)   then     

        veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
        veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
        yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
        elVxloc(ngllx-1)       =   yenihiz* 2.0    

        els1             =   MATMUL(hT,elVxloc)
        els2             =   Tdomain%specel(n)%InvGrad* els1
      else
        els2 = 0.0
      endif 
    else if (Ewave%tfactor == 1)    then
      elVxloc  =  0.0
      if ( n == nel- Ewave%sposi-1  .AND. TimeP%rtime .le. Ewave%timeup)   then

        yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
        elVxloc(ngllx-1) =   yenihiz* 2.0

        els1             =   MATMUL(hT,elVxloc)
        els2             =   Tdomain%specel(n)%InvGrad* els1
      else 
        els2 = 0.0
      endif 
    else if (.not. Ewave%imposource)    then
      els2 = 0.0
    endif  
endif





! RIGID BOUNDARY 
! (NORMALLY INCIDENT WAVE IS DEFINED ON THE LAST POINT, SO NO MULTIPLICATION BY 2)
! For some other tests, it is useful to define it on different elements
! thus, Ewave%sposi is still kept here !

if ( .NOT. Tdomain%specel(nel-1)%PML  .AND. .not. Tdomain%borehole) then
    Ewave%sposi = 0
    if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
      ! Input signal implementation & interpolation                     
      elVxloc    =   0.0
      interpos   =   MOD((it- 1),Ewave%tfactor)
      interpo    =   float(interpos)/ float(Ewave%tfactor) 

     if (n == nel- 1- Ewave%sposi  .AND. TimeP%rtime .le. Ewave%timeup)   then        

        veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
        veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
        yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
        elVxloc(ngllx-1) =   yenihiz


        els1             =   MATMUL(hT,elVxloc)
        els2             =   Tdomain%specel(n)%InvGrad* els1
      else
        els2 = 0.0
      endif 
    else if (Ewave%tfactor == 1)    then
      elVxloc  =  0.0
      if ( n == nel- 1- Ewave%sposi  .AND. TimeP%rtime .le. Ewave%timeup)   then     

        yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
        elVxloc(ngllx-1) =   yenihiz

        els1             =   MATMUL(hT,elVxloc)
        els2             =   Tdomain%specel(n)%InvGrad* els1
      else 
        els2 = 0.0
      endif 
    else if (.not. Ewave%imposource)    then
      els2 = 0.0
    endif  
endif






if (Tdomain%borehole) then
! Borehole uses always the last point of the model (bottom)
! Ewave%sposi is neglected
    yenihiz = 0.0
    if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
      ! Input signal implementation & interpolation                     
      interpos   =   MOD((it- 1),Ewave%tfactor)
      interpo    =   float(interpos)/ float(Ewave%tfactor) 
      if (n == nel- 1  .AND. TimeP%rtime .le. Ewave%timeup)   then        
        veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
        veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
        yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
      endif
    else if (Ewave%imposource .AND. Ewave%tfactor == 1)    then
      if ( n == nel- 1  .AND. TimeP%rtime .le. Ewave%timeup)   then     
        yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
      endif 
    endif 
    Tdomain%Vborehole(xx) = yenihiz 
    els2 = 0.0
endif




    ! Internal strain rate 
    s1 = MATMUL(hT,Vxloc)
    s2 = Tdomain%specel(n)%InvGrad* s1

    ! Strain rate contribution by input signal
    s2 = s2 + els2


    ! PML CONDITION
    if ( Tdomain%specel(n)%PML ) then
      if (xx .NE. 3) then
          LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3)+ &           
                                Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%Mu
      else
          LocStress(:,xx+3) = Tdomain%specel(n)%DumpSx1E* LocStress(:,xx+3)+ &           
                       Tdomain%specel(n)%DumpSx2E* Dt* s2* Tdomain%specel(n)%lambda2mu
      endif
    endif


    ! HERE IT BEGINS 
    if ( .not. Tdomain%specel(n)%PML ) then
      if (.not. Tdomain%rheovisco  .AND.  .not. Tdomain%rheononli) then

        ! ELASTICITY
        gamma2(:,xx+ 3) = s2

        if (xx .NE. 3) then
          LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3)+ &           
                                Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%Mu
        else
          LocStress(:,1) = Tdomain%specel(n)%DumpSx1* LocStress(:,1)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* lambda

          LocStress(:,2) = Tdomain%specel(n)%DumpSx1* LocStress(:,2)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* lambda

          LocStress(:,6) = Tdomain%specel(n)%DumpSx1* LocStress(:,6)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%lambda2mu
        endif

      else if ( Tdomain%rheovisco ) then
        ! Visco pointers & arrays
        DumpSx3   => Evisco%specel(n)%DumpSx3
        DumpSx3YZ => Evisco%specel(n)%DumpSx3YZ
        DumpSx3ZZ => Evisco%specel(n)%DumpSx3ZZ
        DumpSx3VOL=> Evisco%specel(n)%DumpSx3VOL

        zipt    => Evisco%specel(n)%viszipt
        ziptYZ  => Evisco%specel(n)%visziptYZ
        ziptZZ  => Evisco%specel(n)%visziptZZ
        ziptVOL => Evisco%specel(n)%visziptVOL

   	    allocate(zipto(0:ngllx-1)) 
        zipto   = 0.0  

        if (xx == 1) then
          zipto   = 0.0  
          ! Viscoelastic strain rate XZ
          do k = 0,ngllx-1
            do j = 1,8
              zipt(k,j) = zipt(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%visw(j)&
                            *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ zipt(k,j)
            enddo
            DumpSx3(k)  =  (s2(k)- zipto(k))	
          enddo  
        else if (xx == 2) then
          zipto   = 0.0   
          ! Viscoelastic strain rate YZ
          do k = 0,ngllx-1
            do j = 1,8
              ziptYZ(k,j) = ziptYZ(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%visw(j)&
                              *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ ziptYZ(k,j)
            enddo
            DumpSx3YZ(k)  =  (s2(k)- zipto(k))  
          enddo 
        else
          zipto   = 0.0       
          ! Viscoelastic strain rate ZZ
          do k = 0,ngllx-1
            do j = 1,8
              ziptZZ(k,j) = ziptZZ(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%visw(j)&
                              *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ ziptZZ(k,j)
            enddo
            DumpSx3ZZ(k)  =  (s2(k)- zipto(k)) 
          enddo   

          zipto   = 0.0       
          ! Viscoelastic strain rate VOLUMIC
          do k = 0,ngllx-1
            do j = 1,8
              ziptVOL(k,j) = ziptVOL(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%viswp(j)&
                              *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ ziptVOL(k,j)
            enddo
            DumpSx3VOL(k)  =  (s2(k)- zipto(k))  
          enddo   
        endif
              

        if (.not. Tdomain%rheononli ) then
          ! VISCOELASTICITY
          gamma2(:,xx+ 3)    = s2		                              

          if (xx == 1) then
            LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* Dt* Evisco%specel(n)%visM* DumpSx3 

          else if (xx == 2) then
            LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* Dt* Evisco%specel(n)%visM* DumpSx3YZ                     
          else

            ! These sigma_xx and sigma_yy formulas are my approximation given the fact that
            ! there is no formula for lambda and
            ! there are two different corrected strain rate matrices for zz.
            ! although they are never employed in any routine so never influence outputs
            ! sigma_xx = (lambda+2*mu)*eps_zz(P)-(2*mu)*eps_zz(S)


            LocStress(:,1)     = Tdomain%specel(n)%DumpSx1* LocStress(:,1) + &
                                    Tdomain%specel(n)%DumpSx2* Dt &
                                    *(Evisco%specel(n)%visMp* DumpSx3VOL- 2.0*Evisco%specel(n)%visM*DumpSx3ZZ)

            LocStress(:,2)     = Tdomain%specel(n)%DumpSx1* LocStress(:,2) + &
                                    Tdomain%specel(n)%DumpSx2* Dt &
                                    *(Evisco%specel(n)%visMp* DumpSx3VOL- 2.0*Evisco%specel(n)%visM*DumpSx3ZZ)

            LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* Dt* Evisco%specel(n)%visMp* DumpSx3VOL
          endif



        endif
                  
        if  (Tdomain%rheononli)   then
          ! ELASTOVISCOPLASTICITY
          if (xx == 1)&
          gamma2(:,xx+ 3) = DumpSx3 

          if (xx == 2)&
          gamma2(:,xx+ 3) = DumpSx3YZ 

          if (xx == 3)&
          gamma2(:,xx+ 3) = DumpSx3VOL 
        endif
        deallocate(zipto)                     
      else if (.not. Tdomain%rheovisco .AND. Tdomain%rheononli)   then
        ! ELASTOPLASTICITY
        gamma2(:,xx+ 3)    = s2
      endif 
    endif 
    
    deallocate(s1,s2,hT,Vxloc,elVxloc,els1,els2)            
    deallocate(ni,lambda,mu)

    ! Actual strain
    gamma2(:,xx+ 3) = gamma1(:,xx+ 3) + Dt* gamma2(:,xx+ 3) 
  enddo

  ! Strain increment
  Tdomain%specel(n)%deps = gamma2- gamma1


enddo 

end subroutine sigeps