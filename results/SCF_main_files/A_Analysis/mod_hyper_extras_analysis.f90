
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           EXTRA ROUTINES FOR HYPER COPOLYMERS 
!
!           XX.11.2023
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module Molecule_Hyper_Copolymer_Extras
   
   use Global
   use Molecule_Hyper_Copolymer_Analysis
   use Fourier_fftw3

   implicit none


   public :: Initialize_Density
   public :: Convolution_via_fft
   public :: Read_Topology_Data
   public :: Read_Weights 

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           FUNCTION FOR KERNEL USE IN THE GENERATION OF INIT DENSITIES
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!


  function kernel(dis,sig) result(otpt)
    
    double precision, intent (in)   :: dis , sig
    double precision                :: otpt 

    otpt=exp(-(dis**2)/(2*sig**2))

  end function



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINE FOR READING IN TOPOLOGY DATA
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



subroutine Read_Topology_Data(hyper_chains)

    type(hyper_copolymer), dimension(:), allocatable,  intent(out) :: hyper_chains 

    character(100)     :: filename
    character(100)     :: line 
    integer            :: i, j, k ,gen
    integer            :: num_columns=9
    integer            :: molecule_count ,num_segments ,pos, line_num, total_length 
    integer, dimension(9)  :: temp
    integer, dimension(:), allocatable  :: temp_ids  !Holds temporary info


    filename="./hyperbranched_data/topology_data.dat"
   
    ! Open the file
    open(unit=10, file=filename, status='old',action='read')

    molecule_count = 0
    read(10,'(A)',iostat=i) line !First and second lines are just comments
    read(10,'(A)',iostat=i) line 
    read(10,'(A)',iostat=i) line !Second line contains number of different kinds of molecules
    pos = INDEX(line, ',')
    read(line(pos+1:), *) ng 

!    print*, "Number of chains" , ng

    allocate(hyper_chains(ng))   
 
    do
      ! Read a line from the file
      read(10,'(A)',iostat=i) line
      if(i /= 0) exit 



      if (index(line, 'Mol-') == 1) then
        molecule_count = molecule_count + 1

        pos = INDEX(line, ',')
        
        read(line(pos+1:), *) num_segments 
    
        hyper_chains(molecule_count)%num_segments=num_segments




        allocate(hyper_chains(molecule_count)%molecule_data(num_segments, num_columns))
        line_num=0
      else
        line_num=line_num+1
        read(line, *) temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),temp(7),temp(8),temp(9)
        do k=1,num_columns
          hyper_chains(molecule_count)%molecule_data(line_num,k)=temp(k)
        end do
       end if    
    end do

    close(10)


    !Extra information about chain
    do i=1,ng
        hyper_chains(i)%gmax=MAXVAL(hyper_chains(i)%molecule_data(:,6))
        hyper_chains(i)%gmin=MINVAL(hyper_chains(i)%molecule_data(:,6))
        hyper_chains(i)%total_length=sum(hyper_chains(i)%molecule_data(:,2))  
        temp_ids=pack(hyper_chains(i)%molecule_data(:,1),hyper_chains(i)%molecule_data(:,6) .eq. 0)
        j=temp_ids(1)
        hyper_chains(i)%frac(1)=hyper_chains(i)%molecule_data(j,2)*1.0d0/hyper_chains(i)%total_length
        hyper_chains(i)%frac(2)=1.0d0-hyper_chains(i)%frac(1)
        !print*, hyper_chains(i)%frac(1)
    end do
    
        
    
    return


end subroutine Read_Topology_Data

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINE FOR READING IN ENSEMBLE DATA
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine Read_Weights(hyper_chains)

    type(hyper_copolymer), dimension(ng),  intent(inout) :: hyper_chains 

    character(100)     :: filename
    integer            :: io_error ,i


    filename="./hyperbranched_data/weights.dat"

    open(unit=10,file=filename,status='old',action='read',iostat=io_error)

      if(io_error ==0) then

        read(10,*)
        do i=1,ng 
            read(10,*) hyper_chains(i)%weight
        end do

      end if 
    close(unit=10)

end subroutine Read_Weights
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CONVOLUTION VIA FFT FOR BLURRING VOLUME FRACTIONS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 subroutine Convolution_via_fft(vol_frac,kk)
   

    double precision, dimension(gridx,gridy,gridz), intent(in)    :: kk
    double precision, dimension(gridx,gridy,gridz), intent(inout) :: vol_frac 
    double complex, dimension(gridx/2+1,gridy,gridz)   :: vol_frac_fourier, kk_fourier 

    call Real2Fourier(vol_frac,vol_frac_fourier,0) 
    call Real2Fourier(kk,kk_fourier,0) 

    vol_frac_fourier = vol_frac_fourier * kk_fourier  

    call Fourier2Real(vol_frac,vol_frac_fourier,0)
    
    vol_frac=vol_frac*volume

    return

 end subroutine Convolution_via_fft



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALISATION SCHEMES FOR A VARIETY OF GEOMETRIES
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Initialize_Density(geo,dx,dy,dz,init_stat,density)  
   
    
    double precision, dimension(gridx,gridy,gridz,monomer_types) , intent(out):: density
    integer                 , intent(out):: init_stat
    integer,  intent(in)              :: geo
    double precision,  intent(in)     :: dx,dy,dz  
    double precision , dimension(1:monomer_types) :: sigs
    double precision, dimension(gridx,gridy,gridz)              :: temp
    integer              :: x,y,z,mio_error ,i ,p ,q, r, cutoff , n 
    double precision     :: s ,cnx, cny, cnz , px, py, pz ,temp1,temp2, temp3 ,temp4
    double precision     :: pxy, pxz, x1,x2,x3,x4,y1,y2,z1,z2, alp, bet, theta ,d ,sig ,norm , incr 
    double precision     :: phia,phib, phiaver , a ,b ,c 
    double precision, dimension(gridx,gridy,gridz)    :: kk

    double precision, dimension(0:100,monomer_types)   ::  phirs
    double precision, dimension(0:100,monomer_types)   ::  rrs 
    integer         , dimension(monomer_types)         ::  nns
    

    call srand(seed)
    init_stat=0


    if ( polymer_ensemble .eq. 1) then
      phia=phibar*aver_frac(1)
      phib=phibar*aver_frac(2)
    else if (polymer_ensemble .eq. 2) then
      phia=0.5d0*aver_frac(1)
      phib=0.5d0*aver_frac(2)
    end if 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!CALCULATING INITIAL CONFIGURATIONS !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    cnx=sizex/2.0d0
    cny=sizey/2.0d0
    cnz=sizez/2.0d0
    
    incr=1.0d0

    pxy=sizex/sizey
    pxz=sizex/sizez

    
    
    if ( (geo .eq. 1) .or. (geo .eq. 4) ) then !ELLIPSOID
!        x1=(3.0d0/(4.0d0*pi)*pxy*pxz*copolymer%frac(1)*monomer_volume(1)*n)**(1.0d0/3.0d0)
!        x2=(3.0d0/(4.0d0*pi)*pxy*pxz*copolymer%frac(2)*monomer_volume(2)*n+x1**3)**(1.0d0/3.0d0)

        x1=init_rr*(0.5)**(1.0d0/3.0d0)
        x2=init_rr
        y1=x1/pxy
        z1=x1/pxz
        y2=x2/pxy
        z2=x2/pxz
        

        x1=x1*incr
        y1=y1*incr        
        z1=z1*incr

        x2=x2*incr
        y2=y2*incr        
        z2=z2*incr



        sig=(x2**2+y2**2+z2**2)**(0.5d0) 

        print*, "rad 1 : " , x1 , "rad 2 : " ,x2
        if (x2>=cnx) then
            print* , "Box is too small" 
            init_stat=1
            return
        end if
    end if


!    if (geo .eq. 2) then !SPHERICAL VESICLE
!        alp=2.0d0
!        bet=2.0d0

!        x1=(3.0d0/(4.0d0*pi*((alp*bet)**3-1.0d0))*copolymer%frac(2)*monomer_volume(2)*n)**(1.0d0/3.0d0)
!        x2=alp*x1 !Assumption
!        x3=bet*x2 !Assumption
!        x4=(3.0d0/(4.0d0*pi)*copolymer%frac(1)*monomer_volume(1)*n +1.0d0+(alp*bet)**3-alp**3)**(1.0d0/3.0d0)
!       
!        x1=x1*incr
!        x2=x2*incr
!        x3=x3*incr
!        x4=x4*incr


!        sig=x1*0.5d0

!        print*, "rad 0 : " , x1 , "rad 1 : " ,x2, "rad 2 : " , x3 , "rad 3 : " ,x4
!        if (x4>=cnx) then
!            print* , "Box is too small" 
!            init_stat=1
!            return
!        end if
!    end if

!    if (geo .eq. 3) then !JANUS PARTICLE 50-50
!        alp=2.0d0
!        bet=2.0d0
!!        theta=5*pi/4.0d0 !Angle of normal vector of the plane cutting the particle in the x-y plane
!        theta=3*pi/4.0d0
!        d=-(Cos(theta)*cnx+Sin(theta)*cny)

!        x1=(3.0d0/(4.0d0*pi)*monomer_volume(1)*n)**(1.0d0/3.0d0)
!        x1=x1*1.2d0

!        sig=x1*0.4d0


!        print*, "rad 0 : " , x1 
!        if (x4>=cnx) then
!            print* , "Box is too small" 
!            init_stat=1
!            return
!        end if
!    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!ASSIGNING 1S AND 0S TO EACH MONOMER TYPE!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
    density(:,:,:,:)=0.0d0

    do z = 1, gridz
     do y = 1, gridy
         do x = 1 , gridx

            if (geo .eq. -1 ) then
                c=rand()
                
                density(x,y,z,1)=c
               
                density(x,y,z,2)=1-c
                

            end if

            if (geo.eq.0) then !LAMELLAR
                
                px=(x*1.0d0/dble(gridx))*sizex
                temp1=abs((px-cnx)/init_rr)

                if ( temp1 <=1.0d0) then
                    density(x,y,z,1)=1.0d0
                end if

                if (( temp1 >1.0d0) .and. ( temp1 <=1.5d0)) then
                    density(x,y,z,2)=1.0d0
                end if
  
            else if ((geo .gt. 0) ) then 
                px=(x*1.0d0/dble(gridx))*sizex
                py=(y*1.0d0/dble(gridy))*sizey
                pz=(z*1.0d0/dble(gridz))*sizez

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!ELLIPSOID!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (geo .eq. 1) then 
                    temp1=( ((px-cnx)/x1)**2+((py-cny)/y1)**2+((pz-cnz)/z1)**2 )
                    temp2=( ((px-cnx)/x2)**2+((py-cny)/y2)**2+((pz-cnz)/z2)**2 )

                    if ( temp1 <=1.0d0) then

                        density(x,y,z,1)=1.0d0
                    end if
                    
                    if ( (temp2 <=1.0d0) .and. (temp1>1.0d0)) then
                        density(x,y,z,2)=1.0d0
                    end if

                    !if (temp2>1.0) then
                      !  density(x,y,z,1)=phia
                     !   density(x,y,z,2)=phib
                    !end if

    
    
                end if 


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!VESICLE!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (geo .eq. 2) then 
                    
                    temp1=(( (px-cnx)**2+((py-cny))**2+((pz-cnz))**2 )/x1 )**2 
                    temp2=(( (px-cnx)**2+((py-cny))**2+((pz-cnz))**2 )/x2 )**2  
                    temp3=(( (px-cnx)**2+((py-cny))**2+((pz-cnz))**2 )/x3 )**2 
                    temp4=(( (px-cnx)**2+((py-cny))**2+((pz-cnz))**2 )/x4 )**2
                    
     
                    if ( (temp2 <=1.0d0) .and. (temp1>=1.0d0)) then
                        density(x,y,z,1)=1.0d0
                    end if
    
                    if ( (temp3 <=1.0d0) .and. (temp2>1.0d0)) then
                        density(x,y,z,2)=1.0d0
                    end if
           
                    if ( (temp4 <=1.0d0) .and. (temp3>1.0d0)) then
                        density(x,y,z,1)=1.0d0
                    end if
                

               end if 

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!JANUS PARTICLE!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (geo .eq. 3) then 
                
                temp2=Cos(theta)*px+Sin(theta)*py+d
                if (temp1<=1.0d0) then
                     
                    
                    if ( temp2>=0.0d0) then
!                        density(x,y,z,1)=1.0d0*Exp(-temp1/(2*sig**2))
                         density(x,y,z,1)=1.0d0
                    else
!                        density(x,y,z,2)=1.0d0*Exp(-temp1/(2*sig**2))0
                         density(x,y,z,2)=1.0d0
                    end if

                end if

               end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!CYLINDER!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (geo .eq. 4) then 
                    temp1=( ((px-cnx)/x1)**2+((py-cny)/y1)**2 )
                    temp2=( ((px-cnx)/x2)**2+((py-cny)/y2)**2 )

                    if ( temp1 <=1.0d0) then

                        density(x,y,z,1)=1.0d0
                    end if
                    
                    if ( (temp2 <=1.0d0) ) then
                        density(x,y,z,2)=0.5d0
                    end if

                end if 



            end if 
           end do
         end do
      end do
    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!NOW BLURRING CONFIGURATION WITH KERNEL!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   


    if ((geo .ge. 0) .and. (blur .eq. 1)) then 
        sigs(1)=sig1*sig
        sigs(2)=sig2*sig

        do i=1,monomer_types
          cutoff=nint(3.0d0*sigs(i)/(dx**2+dy**2+dz**2)**(0.5d0)) !cutoff at 3 sigma

          if ( (cutoff .gt. gridx/2) .or. ((cutoff .gt. gridy/2) .or. (cutoff .gt. gridz/2))) then
            print*, "WARNING CUTOFF USED FOR INIT CONFIGURATION IS TOO BIG USE SMALLER SIG1 AND/0R SIG2"
          end if


          print*, "Cutoff was set to : " , cutoff , " and sigma is " , sig

          if (cutoff .le. 1) then 
            print*, " Cutoff radius is too low so we will skip blurring"
            exit
          end if

          kk=0
          do p = -cutoff, cutoff
              do q = -cutoff, cutoff
                  do r = -cutoff, cutoff
                      
                      temp1=(p**2+q**2+r**2)**(0.5d0)
                      kk(p,q,r) = kernel(temp1,sigs(i))
                  end do
              end do
          end do

          norm = sum(kk)*dvol
          kk = kk / norm

          call Convolution_via_fft(density(:,:,:,i),kk(:,:,:))
        end do

        
      
    end if
    



    !Dividing phi by monomer volume to get density

    density(:,:,:,1)=density(:,:,:,1) / monomer_volume(1) 
    density(:,:,:,2)=density(:,:,:,2) / monomer_volume(2) 



   close (unit=112)


  print*, "DENSITY INITIALIZATION DONE"
  end subroutine Initialize_Density

end module Molecule_Hyper_Copolymer_Extras

