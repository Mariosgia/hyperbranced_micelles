!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!        ROUTINES FOR FLORY HUGGINS TYPE INTERACTIONS IN SOLUTION
!             IN A MULTICOMPONENT SYSTEM
!             WITH INCOMPRESSIBLE IMPLICIT SOLVENT 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Interactions_FloryHuggins_Solution
!
  use Global
!
  implicit none
!
  private
!
  integer, parameter     :: components = monomer_types     ! solvent is component 0
!
!
  double precision, save :: small_number = 1.0d-3
!
  public :: Initialize_Interactions
  public :: Calculate_Mean_Fields
  public :: Interaction_Energy
  public :: Check_Volume_Fraction
!  public :: Init_tt

contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE INTERACTIONS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Interactions
!
    character*30 :: file_name
    integer      :: io_error
    integer      :: component, component2
!
    file_name='input_interactions'
    open(unit=20,file=file_name,status='old',action='read',iostat=io_error) 
!
    chi = 0.0d0
    if(io_error == 0) then
      read(20,*) kappa
      read(20,*)
      do component = 0, components
        do component2 = component+1, components
          read(20,*) chi(component,component2)
          chi(component2,component) = chi(component,component2)
        end do
      end do
      read(20,*)
      do component = 0, components
        read(20,*) monomer_volume(component)
      end do
    else
      print*, 'Error', io_error, ' while trying to open', file_name 
    end if
!
    do component = 1, components
      do component2 = 1,components
        chi_eff(component2, component) = chi(component2, component) &
           - chi(0,component) - chi(0,component2)
      end do
    end do
    zeta = 1.0d0/monomer_volume(0)
!
!  print*, " 1. Zeta , 2. ChiAS 3.ChiAB 4.ChiBS 5. ChiAAeff 6. ChiABeff 7. ChiBBeff"
!  print*,  zeta, chi(1,0) ,chi(1,2) , chi(2,0) , chi_eff(1,1) , chi_eff(1,2) , chi_eff(2,2) 


  end subroutine Initialize_Interactions
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE INTERACTION ENERGY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  double precision function Interaction_Energy(density)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
!

    double precision, dimension(:,:,:,:), allocatable   :: volume_fraction
    double precision, dimension(:,:,:), allocatable   :: solvent_volume_fraction
    double precision, dimension(:,:,:), allocatable   :: local_interaction_energy

    double precision :: field_threshold, energy_threshold
    integer          :: component, component2,allocstatus
!

    allocate(volume_fraction(gridx,gridy,gridz,components),stat=allocstatus)
    allocate(solvent_volume_fraction(gridx,gridy,gridz),stat=allocstatus)
    allocate(local_interaction_energy(gridx,gridy,gridz),stat=allocstatus)
    
    field_threshold  = - log(small_number)
    energy_threshold = (small_number*log(small_number)- small_number)
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do

    solvent_volume_fraction = 1.0d0 - sum(volume_fraction,dim=4)
    where (solvent_volume_fraction .ge. small_number)
      local_interaction_energy =  &
        ( solvent_volume_fraction*log(solvent_volume_fraction)- solvent_volume_fraction ) 
    else where
      local_interaction_energy = energy_threshold  &
          - field_threshold * (solvent_volume_fraction - small_number) &
          + 0.5d0*kappa*(solvent_volume_fraction - small_number)**2
    endwhere
    local_interaction_energy = local_interaction_energy * zeta

    do component = 1,components
      local_interaction_energy = local_interaction_energy &
           + chi(0,component) * volume_fraction(:,:,:,component)
      do component2 = 1, components
        local_interaction_energy = local_interaction_energy  &
           + 0.5d0*chi_eff(component,component2) &
              * volume_fraction(:,:,:,component)*volume_fraction(:,:,:,component2)
      end do
    end do

    interaction_energy = dvol*sum(local_interaction_energy)
 
    deallocate(volume_fraction,solvent_volume_fraction,local_interaction_energy)
    
  end function Interaction_Energy
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE MEAN FIELDS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Mean_Fields(density,field)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,gridz,components), intent(out):: field
!
    double precision, dimension(:,:,:,:), allocatable   :: volume_fraction
    double precision, dimension(:,:,:), allocatable   :: solvent_volume_fraction
    double precision, dimension(:,:,:), allocatable   :: field_0

    double precision :: field_threshold
    integer          :: component, component2,allocstatus


    allocate(volume_fraction(gridx,gridy,gridz,components),stat=allocstatus)
    allocate(solvent_volume_fraction(gridx,gridy,gridz),stat=allocstatus)
    allocate(field_0(gridx,gridy,gridz),stat=allocstatus)


    field_threshold = - log(small_number)
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do
!
    solvent_volume_fraction = 1.0d0 - sum(volume_fraction,dim=4)
    where (solvent_volume_fraction .ge. small_number)
      field_0 = - log(solvent_volume_fraction)
    else where
      field_0 = field_threshold - kappa * (solvent_volume_fraction-small_number)
    endwhere
    field_0 = field_0 * zeta
!
    do component = 1,components
      field(:,:,:,component) = field_0 + chi(0,component)
      do component2 = 1, components
        field(:,:,:,component) = field(:,:,:,component) &
          + chi_eff(component,component2)*volume_fraction(:,:,:,component2)
      end do
      field(:,:,:,component) = field(:,:,:,component) &
           * monomer_volume(component)
    end do

    
    if ( (con_rr .gt. 1e-4) .and. ( (constraint .eq. 1) .or. (constraint .eq. 3) ) ) then

        con_const=sum(con_mask*volume_fraction(:,:,:,1))*dvol
        con_const=kapcon*(con_const-con_phi)

!        print* , "Const field term is " , con_const
        where (con_mask .eq. 1 )
            field(:,:,:,1)=field(:,:,:,1) + con_const 
        end where

   end if



    if ( (con_rr .gt. 1e-4) .and. (constraint .eq. 2) ) then

        field(:,:,:,1)=field(:,:,:,1) + con_add 

   end if
    
    deallocate(volume_fraction,solvent_volume_fraction,field_0)
 
  end subroutine Calculate_Mean_Fields



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CHECK THAT SUM OF AVERAGE VOLUME FRACTIONS GIVES ONE.
!		(FOR USE IN THE CANONICAL ENSEMBLE)
!           THIN FILMS: COUNT ONLY ACCESSIBLE GRID POINTS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Check_Volume_Fraction(density)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density

!
    double precision, dimension(:,:,:,:), allocatable   :: volume_fraction


    integer           :: component,allocstatus
    
    allocate(volume_fraction(gridx,gridy,gridz,components),stat=allocstatus)

    aver_vol_fra(0)=1.0d0
    average_volume_fraction=0.0d0
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
      aver_vol_fra(component)=sum( volume_fraction(:,:,:,component) ) / dble(gridx*gridy*gridz)
      aver_vol_fra(0)=aver_vol_fra(0)-aver_vol_fra(component)
      average_volume_fraction=average_volume_fraction+aver_vol_fra(component)
    end do
!
    
!
!    write(30,*) 'Average total volume fraction is', average_volume_fraction
!    print*, 'Average total volume fraction is', average_volume_fraction

!    if ( abs(average_volume_fraction - 1.0d0) .gt. 1.0d-3) then

!      write(30,*) 'Caution: Differs noticeably from 1'
!      print*, 'Caution: Differs noticeably from 1'
!    end if
!
   deallocate(volume_fraction)

  end subroutine Check_Volume_Fraction



!  subroutine Init_tt()


!  tt=monomer_volume(1)*(-zeta*DLOG(1.0d0-phibar)+&
!                        phibar*(chi_eff(1,1)*frac(1)*frac(1)+2.0d0*chi_eff(1,2)*frac(1)*frac(2)+chi_eff(2,2)*frac(2)*frac(2))+&
!                        chi(1,0)*frac(1)+chi(2,0)*frac(2))


!  end subroutine Init_tt

end module Interactions_FloryHuggins_Solution
