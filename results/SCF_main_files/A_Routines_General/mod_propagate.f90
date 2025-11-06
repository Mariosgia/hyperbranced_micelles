!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           GENERAL PROPAGATOR SUBROUTINES
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Propagate
!
  use Global
!
  implicit none
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE STEP: CALCULATE EXP(-field*ds*0.5) 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Expfield(field,expfield)
!
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out)   :: expfield

!
      expfield = exp(-0.5d0*ds*field)
!
    return

  end subroutine Calculate_Expfield
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           PROPAGATE prop BY ONE STEP ds
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Propagate_Step(prop,expfield,th_num)
!
    use Fourier_fftw3            
!
    implicit none
!
    double precision, dimension(gridx,gridy,gridz), intent(in)    :: expfield
    integer, intent(in)                                           :: th_num
    double precision, dimension(gridx,gridy,gridz), intent(inout) :: prop
    double complex, dimension(:,:,:), allocatable           :: propfourier
!    double complex, dimension(gridx/2+1,gridy,gridz)   :: propfourier
    integer                         :: x,y,z ,iAllocStatus
    

     allocate(propfourier(gridx/2+1,gridy,gridz),stat=iAllocStatus) 

    prop = expfield * prop  

    call Real2Fourier(prop,propfourier,th_num)

    propfourier = explaplace * propfourier

    call Fourier2Real(prop,propfourier,th_num)

    prop = expfield * prop  

    deallocate(propfourier)
   return

  end subroutine Propagate_Step

end module Propagate
