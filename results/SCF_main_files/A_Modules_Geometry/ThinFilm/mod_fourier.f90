!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS FOR FOURIER TRANSFORMS - THIN FILM GEOMETRY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Fourier_fftw3
!
  use global
  implicit none

  private

  include "fftw3.f"

  integer*8           ::  sinplan,sinbackwardplan,ftplan,ftbackwardplan
  ! for fortran complex values are supported by fftw, !!!
  ! fftw's row major is warped to column major in Fortran!!!

  double precision, dimension(gridx, gridy) :: ftin 
  double complex, dimension(gridx/2+1, gridy) :: ftout
  double precision, dimension(gridz) :: sinin
  double precision, dimension(gridz) :: sinout
  double precision, dimension(gridx,gridy,gridz) :: out
  double precision, dimension(gridx, gridy) :: ftbackwardout 
  double complex, dimension(gridx/2+1, gridy) :: ftbackwardin

  double complex, save, dimension(gridx/2+1,gridy, gridz) :: laplace
  double complex, save, dimension(gridx/2+1,gridy, gridz) :: explaplace

  public :: Initialize_Fourier
  public :: End_Fourier
  public :: Fourier2Real
  public :: Real2Fourier
  public :: Calculate_Volume_Parameters
  public :: Calculate_Laplace

  public :: laplace
  public :: explaplace

contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    INITIALIZATION ROUTINE
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Fourier

    call dfftw_plan_dft_r2c_2d(ftplan,gridx,gridy, ftin, ftout,FFTW_PATIENT)
    call dfftw_plan_dft_c2r_2d(ftbackwardplan,gridx,gridy, ftbackwardin, ftbackwardout,FFTW_PATIENT)
    call dfftw_plan_many_r2r(sinplan,1,gridz,1,sinin,gridz,1,1,sinout,gridz,1,1, FFTW_RODFT00,FFTW_PATIENT)

  end subroutine Initialize_Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    CLOSING ROUTINE
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine End_Fourier

    call dfftw_destroy_plan(ftplan)
    call dfftw_destroy_plan(ftbackwardplan)
    call dfftw_destroy_plan(sinplan)

  end subroutine End_Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    FOURIER TRANSFORM (BACKWARD)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine fourier2real(prop,propfourier)

    integer                         :: x,y,z
    double precision, dimension(gridx,gridy,gridz), intent(out) :: prop
    double complex, dimension((gridx/2+1),gridy,gridz), intent(in) :: propfourier
    do z=1,gridz
      ftbackwardin=propfourier(:,:,z)
       call dfftw_execute(ftbackwardplan)
       out(:,:,z)=ftbackwardout
    end do
    do y=1,gridy
       do x=1,gridx
          sinin=out(x,y,:)
          call dfftw_execute(sinplan)              
          prop(x,y,:)=sinout
       end do
    end do
    return
  end subroutine fourier2real
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    FOURIER TRANSFORM (FORWARD)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine real2fourier(prop,propfourier)

    integer                         :: x,y,z
    double precision, dimension(gridx,gridy,gridz), intent(in) :: prop
    double complex, dimension((gridx/2+1),gridy,gridz), intent(inout) :: propfourier
  
    do y=1,gridy
       do x=1,gridx
          sinin=prop(x,y,:)
          call dfftw_execute(sinplan)             
          out(x,y,:)=sinout
       end do
    end do
    do z=1,gridz
       ftin=out(:,:,z)
       call dfftw_execute(ftplan)
       propfourier(:,:,z)=ftout/(2*(gridz+1)*gridx*gridy)
    end do
    return
  end subroutine real2fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    ROUTINE TO TEST FOURIER TRANSFORMS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine fourier_test
!
    use mtmod

    integer             :: x,y,z
    double precision, dimension(gridx,gridy,gridz)   :: testfield
    double precision, dimension(gridx,gridy,gridz)   :: testfieldnew
    double complex, dimension(gridx/2+1,gridy,gridz) :: testfourier

    do z=1,gridz
       do y=1,gridy
          do x=1,gridx
             testfield(x,y,z)=1.0d0*(1.0d0-2.0d0*grnd())
          end do
       end do
    end do
    do z=1,gridz
       do y=1,gridy
          write(10,'(8g12.4)') (testfield(x,y,z),x=1,8)
       end do
    end do
    call real2fourier(testfield,testfourier)
    call fourier2real(testfieldnew,testfourier)
    do z=1,gridz
       do y=1,gridy
	  write(12,'(8g12.4)') (testfieldnew(x,y,z),x=1,8)
       end do
    end do

   return

  end subroutine fourier_test
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE: CALCULATE VOLUME PARAMETERS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Volume_Parameters

    volume = sizex*sizey*sizez
    dvol   = volume/dble(gridx*gridy*(gridz+1)) 
    ! z=0 gridpoint is not stored, but still counts !

   return

  end subroutine Calculate_Volume_Parameters
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE: CALCULATE EXP(LAPLACE*DT) IN FOURIER SPACE 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Laplace

    integer :: x,y,z
    double precision :: qx,qy,qz
!
     do z=1,gridz !remember indexshift in fortran
       do y=1,gridy !remember indexshift in fortran
          do x=1,(gridx/2+1) !remember indexshift in fortran
             qx = 2.0d0*pi * dble(x-1) / sizex
             qy = 2.0d0*pi * dble(y-1) / sizey
             qz = pi * dble(z) / sizez
             !sintrafo has a shift in it (zeroth-component=0 isnt stored) 
             if (qy .gt. (gridy/2+1)) then
               qy = 2.0d0*pi * (dble(gridy)-dble(y-1)) / sizey
             end if
             laplace(x,y,z) = dcmplx(-(qx*qx+qy*qy+qz*qz), 0.0d0)
             explaplace(x,y,z) = dcmplx(exp(-ds*(qx*qx+qy*qy+qz*qz)), 0.0d0)
          end do
       end do
     end do

    return

  end subroutine Calculate_Laplace
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Fourier_fftw3
