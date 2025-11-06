!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    GLOBAL VARIABLES (SYSTEM DIMENSIONS ETC.)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Global
!
  implicit none
  double precision, parameter                  :: pi = 3.141592653589793d0
!
  double precision                             :: ds  !1/average chain length 
  integer, parameter                           :: monomer_types=2

  double precision, dimension(0:monomer_types) :: monomer_volume ,aver_vol_fra
  double precision, dimension(monomer_types)   :: aver_frac

  integer                                      :: gridx,gridy, gridz, threads
  integer                                      :: ng !Number of chains
  double precision                             :: sizex, sizey, sizez, volume, dvol ,sig1 ,sig2 , init_rr
  double precision                             :: phibar ,tt , lll ,average_chain_length
  integer                                      :: every_out , seed , poly_stru ,blur ,polymer_ensemble
  integer, parameter                           :: mixing_dim = 10 ! To be defined in Module Global

  double precision                             :: kappa               ! Helfand compressibility parameter
  double precision                             :: zeta

  double precision                             :: chi(0:monomer_types,0:monomer_types) ! interaction parameters
  double precision                             :: chi_eff(monomer_types,monomer_types) ! interaction parameters
  double precision                             :: average_volume_fraction        ! average polymer volume fraction
  logical                                      :: compiled_with_openmp = .false.
  integer                                      ::  umbrella
  double precision                             ::  um_rr,um_phi

  integer                                      :: constraint
  double precision                             :: con_rr, asp_ratio , con_phi,con_shift,kapcon,con_const,con_ene
  integer, dimension(:,:,:), allocatable       :: con_mask
  double precision, dimension(:,:,:), allocatable    :: con_add

  double precision, dimension(:,:,:), allocatable  :: rad

end module Global
