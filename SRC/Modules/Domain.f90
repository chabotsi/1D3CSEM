module sdomain
! This module contains the derived type: computation domain, 
! with its elements, and element properties.

use selement
!use sfault

use subdomains

!use snode

!type :: DTN_1
!real :: Length, Fmax
!end type

!type :: TimeMultires
!integer :: node
!integer :: Fninf, Fnsup, Cn, iFi, iC
!end type

type :: PML_Conditions
logical :: Filtering, Left
integer :: n_pow
real :: A_pow,freq
end type

type :: node
integer :: num
real :: coord
integer, dimension(:), pointer :: elemconnect
 character (len =2) :: type  
end type

type :: domain

logical :: source

logical :: borehole
real, dimension (:), pointer :: Vborehole

logical :: rheovisco, rheononli , eliflogic
logical :: Epressmod

integer :: nelem, nnode ,npoints, nsubdomain
integer :: n_pml
integer :: nbdirichlet  , Ewlevel        
integer, dimension (:), pointer :: typemat,  nbelem

real :: dtmaximum      
real, dimension (:,:), pointer :: paramdiric
real, dimension (:,:), pointer :: dxmaxima

 character (len =1), dimension (:),allocatable :: dirichlet
 character (len =20) :: Title_simulation

type(element), dimension(:), pointer :: specel
type(node), dimension(:), pointer :: node
type (Subdomain), dimension (:), pointer :: Sub_Domain
type (PML_Conditions), dimension (:), pointer :: PMLCondition

end type

end module sdomain
