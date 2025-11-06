!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       ANALYSIS FILE THAT WORKS THE SAME WAY AS MAIN ALTHOUGH JUST USES THE LAST CONFIGURATION
!           	     TO CALCULATE CERTAIN QUANTITIES THAT ARE NOT EXPORTED IN THE MAIN
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_hyper
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

  use Global
  use Molecule_Hyper_Copolymer_Analysis
  use Molecule_Hyper_Copolymer_Extras
  use Interactions_FloryHuggins_Solution
  use Iterate
  USE OMP_LIB
!
  implicit none


  double precision    ::  number_of_chains
  double precision    ::  energy , field_energy, chains_energy
  integer             ::  itermax_scf    ! maximum number of iterations
  integer             ::  iterations     ! actual number of iterations
  double precision    ::  accuracy_scf   ! desired accuracy in iterations
  double precision    ::  accuracy       ! actual accuracy
  double precision    ::  lambda         ! mixing parameters
  integer             ::  mixing_type    ! type of mixing
  integer             ::  stat           ! status of converged solution
  integer             ::  geo            ! Init geometry
  real                :: start, fin
!      -------------  1: simple mixing, 2: lambda mixing, 3: anderson mixing ----------

  type(hyper_copolymer), dimension(:), allocatable       :: hyper_chains !Creates a number of chain objects

  double precision, dimension(:,:,:,:), allocatable ::density,newdensity,dfield, newfield ,field
  double precision, dimension(:,:,:,:), allocatable ::saddledensity,saddlefield,temp_dens
  integer, dimension(:,:,:,:), allocatable          :: umbrella_mask
  double precision, dimension(:,:),allocatable :: var,dvar
  double precision, dimension(:,:,:), allocatable :: temp_arr

  logical          :: qexist
  integer          :: qread,every_show
  integer          :: allocStatus,init_status
  character*24     :: file_name ,x1
  character*80     :: output_str0 ,output_str1,output_str2
  character(len=8) :: fmtout ! format descriptor
  integer          :: polymer_label, monomer_type_A, monomer_type_B ,grad_data
  double precision :: dx, dy, dz, energy_homogeneous  , cx, cy, cz

  integer          :: io_error, x,y,z,iter,iter_status, i ,u,j,k,id ,xc,yc,zc,con_temp,gen,added_index,arg_count,arg_val

  double precision ::  total_rho , dummy 

  integer, dimension(:), allocatable :: temp_ids,temp_max
  character(len=100) :: arg_char

  double precision ::cnx,cny,cnz,dis,temp_val,radx,rady,radz , phiori ,max_point, mean_dis,norm



arg_val=0
arg_count= command_argument_count()
if (arg_count .eq. 1) then
    call get_command_argument(1, arg_char)  ! Retrieve command-line argument
    READ (arg_char, *) arg_val
end if

print*, " Argumen value is " , arg_val

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READING PARAMETERS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

  file_name='input_general'
  call read_input_parameters(file_name)

!!!!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!    ALLOCATING ARRAYS

!!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   allocate(density(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(newdensity(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(field(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(newfield(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(dfield(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(saddledensity(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(saddlefield(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(var(gridx*gridy*gridz*monomer_types,0:mixing_dim),stat=allocStatus)
   allocate(dvar(gridx*gridy*gridz*monomer_types,0:mixing_dim),stat=allocStatus)
   allocate(umbrella_mask(gridx,gridy,gridz,monomer_types),stat=allocStatus)
   allocate(temp_arr(gridx,gridy,gridz),stat=allocStatus)
   allocate(con_mask(gridx,gridy,gridz),stat=allocStatus)
   allocate(rad(gridx,gridy,gridz),stat=allocStatus) 
   allocate(temp_max(gridx),stat=allocStatus) 
   allocate(temp_dens(gridx,gridy,gridz,monomer_types),stat=allocStatus)


!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!	     INITIALIZATIONS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
! ------ GENERAL INITIALIZATIONS 
!


  dx = sizex/dble(gridx)
  dy = sizey/dble(gridy)
  dz = sizez/dble(gridz)


  call Read_Topology_Data(hyper_chains)
  
  call Read_Weights(hyper_chains)  

  call Initialize_Interactions

  average_chain_length=sum(hyper_chains(:)%weight*hyper_chains(:)%total_length)
  ds=1.0d0/(average_chain_length)
  print*, "ds" , ds
    

  aver_frac(1)=0.0d0
  

  aver_frac(1)=sum(hyper_chains(:)%weight*hyper_chains(:)%total_length*hyper_chains(:)%frac(1))/average_chain_length
  aver_frac(2)=1.0d0-aver_frac(1) 

  print* , aver_frac(1) ,aver_frac(2)


  lll=phibar/(monomer_volume(1))

  call Initialize



 !$ compiled_with_openmp = .true.

if (compiled_with_openmp) then
  write(*,*) 'OpenMP used...'
    
else
  write(*,*) 'Openmp not used...'
end if

  if (compiled_with_openmp) then
    if ( (mod(ng,threads)) .ne. 0) then
        print*, "ng has to be a multiple of threads for faster implemantation"
    end if
  end if
    


!!!!! ------ INITIALIZE POLYMERS

    polymer_label = 1
    monomer_type_A = 1
    monomer_type_B = 2

    do i=1,ng
        call Initialize_Hyper(hyper_chains(i), i,&
                     (/monomer_type_A, monomer_type_B/),&
                     polymer_ensemble)
    end do



! ------ INITIALIZE INTERACTIONS
!


  con_temp=constraint
  !Creating array with radius valuse from the center
  do z = 1, gridz
    do y = 1, gridy
      do x = 1 , gridx 
         rad(x,y,z)=( ((x-gridx/2)*dx)**2 +((y-gridy/2)*dy)**2+((z-gridz/2)*dz)**2  )**0.5d0
      end do
    end do
  end do



  con_temp=constraint

  if (constraint .ne. 0) then
  !Creating array with radius valuse from the center
    do z = 1, gridz
      do y = 1, gridy
        do x = 1 , gridx 
           rad(x,y,z)=( ((x-gridx/2)*dx)**2 +((y-gridy/2)*dy)**2+((z-gridz/2)*dz)**2  )**0.5d0
        end do
      end do
    end do

    !Need to read con densities for constraint potential to be the same as in the converged stata
    open (unit=33, file='con_densities.in', status='old', action='read', iostat=io_error) 
        if(io_error /=0) then                   
          write(*,*) 'Error', io_error , ' while trying to open con_densities.in'
          stop
        end if
        read (33,*) density
        print*, "Read con densities"
    close (unit=33)



  !Creating array with mask where constraint potential will apply
    if ( (constraint .eq. 1) .or. (constraint .eq. 3)) then
      if (qread .ne. 1) then
          print*, "Consrtaint 1 uses the initial volume of the core and that relates to the initial configuration"
          stop 
      end if

      if (constraint.eq. 3) then 
        temp_val=COUNT(density(:,:,:,1) >= 0.5d0)*dvol !volume of micelle core
        temp_val=(temp_val*3.0/4.0/pi)**(1.0/3.0)
        print*, " radius micelle : " , temp_val
        rady=temp_val/asp_ratio**(1.0/3.0)
        radz=rady
        radx=asp_ratio*rady
      else
        radx=asp_ratio*con_rr
        rady=con_rr
        radz=con_rr
      end if

      con_mask=0
      do z = 1, gridz
        do y = 1, gridy
          do x = 1 , gridx 
             temp_val=( ((x-gridx/2)*dx)**2/radx**2 +((y-gridy/2)*dy)**2/rady**2+((z-gridz/2)*dz)**2/radz**2 )**0.5d0

             if (temp_val .le. 1) then
                con_mask(x,y,z)=1
             end if

          end do
        end do
      end do
      
      print*, "Volume ellipsoid : " , sum(con_mask)*dvol   
      print*, "Volume core : " , COUNT(density(:,:,:,1) >= 0.5d0)*dvol
   
      if (constraint.eq. 3) then 
        phiori=SUM(density(:,:,:,1), MASK = density(:,:,:,1) >= 0.5d0)*monomer_volume(1)*dvol
      else
        phiori=sum(con_mask*density(:,:,:,1))*monomer_volume(1)*dvol
      end if

      con_phi=phiori+con_shift

      if ( abs(con_rr) .le. 1e-4) kapcon=0.0

      print*
      print*, "Phi average of solvophobic part in con_rr and shifted phi and kapcon : " , phiori , con_phi , kapcon
      print*, "Radx Rady Radz : " , radx , rady , radz
      print*


   end if 



    if (constraint .eq. 2) then
    con_add(:,:,:)=0.0d0
    do x=1,gridx

        if ( ( abs(((x-gridx/2)*dx)) .le. con_rr ) ) then
            con_add(x,gridy/2,gridz/2)=-kapcon 
        end if       
    end do

   end if 

  end if

!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	   RUNNING QREAD CASES
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
if (qread.eq.1) then
!
    open (unit=33, file='cnf_densities.in', status='old', action='read', iostat=io_error)
    if(io_error /=0) then                   
      write(*,*) 'Error', io_error , ' while trying to open cnf_densities.in'
    end if
    read (33,*) density
    print*, "Read densities"
    close (unit=33)

end if



if (qread.ne.0) then
 call Calculate_Mean_Fields(density,field)
 call Check_Volume_Fraction(density(:,:,:,:))
print*, " Initial Vol. fraction. average. : " , average_volume_fraction

end if




constraint=0 !Need to calculate mus without the constraint potential
call Calculate_Homogeneous_Energy()
constraint=con_temp

   


    

call Calculate_Densities(field,newdensity,hyper_chains,1)

call Calculate_Mean_Fields(newdensity,newfield)


dfield = newfield - field
call  Check_Accuracy(density,newdensity,accuracy)

call Check_Volume_Fraction(newdensity(:,:,:,:))

call Calculate_Energy(field,density,hyper_chains)                                        ! FOR TESTING



!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	                ANALYSIS FILES
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!




saddledensity=density
saddlefield=field

call Calculate_Centre(saddledensity)



call system("mkdir terminal_distributions")
call system("mkdir terminal_gen_index_info")

do i=1,ng
  

  write (x1,*) i
  output_str1='./terminal_gen_index_info/mol_id_'//TRIM(ADJUSTL(x1))//'.dat' 
  output_str0='./terminal_distributions/mol_id_'//TRIM(ADJUSTL(x1))//'.dat'

  open(unit=222,file=output_str1, status='replace',action='write',iostat=io_error)
    write (222, *) "# 1.Target node 2.Generation 3.Max point 4.Mean distance"


  temp_ids=pack(hyper_chains(i)%molecule_data(:,1),hyper_chains(i)%molecule_data(:,3) .eq. 0)   !Choosing terminal segments
  temp_arr=0


  do k=1,size(temp_ids) !going over each segment
    id=temp_ids(k) !segment id of polymer i
    added_index=hyper_chains(i)%molecule_data(id,5)
    gen=hyper_chains(i)%molecule_data(id,6)

    if ( added_index .ge. 92) then !Monomers with indices below 92 are part of the linear and core parts which we want to ignore

      temp_max=MAXLOC(hyper_chains(i)%node_distribution(xc:,yc,zc,id))
      max_point=((xc+temp_max(1))*dx)-cx  !x-coord of maximum point of node distribution
      
      call Integrate(hyper_chains(i)%node_distribution(xc:,yc,zc,id),rad(xc:,yc,zc),dx,gridx-xc,norm)

      call Integrate(( rad(xc:,yc,zc)*hyper_chains(i)%node_distribution(xc:,yc,zc,id) ),rad(xc:,yc,zc),dx,gridx-xc,mean_dis)

      mean_dis=mean_dis/norm

      
!      print*, "MAX LOC" , temp_max(1), size(temp_max),max_point
      write (222, *) added_index, gen , max_point, mean_dis 
    end if

    temp_arr=temp_arr+hyper_chains(i)%node_distribution(:,:,:,id)*exp(hyper_chains(i)%mu) 
  end do 

  
  open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
        write (111, *) size(temp_ids) !Number of terminal segments
        do x = 1  , gridx
            write (111, *) cx-(x*dx), temp_arr(x,yc,zc)    
        end do
  close(unit=111)

  close(unit=222)


end do


call system("mkdir new_weigths")
output_str0='./new_weigths/weights.dat'
output_str1='./rhos.dat'
output_str2='./micelle_rhos.dat'
open(unit=333,file=output_str2, status='replace',action='write',iostat=io_error)
open(unit=222,file=output_str1, status='replace',action='write',iostat=io_error)
 open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
    write (333, *) "#Micelle Rhos"
    write (222, *) "#Rhos"
    write (111, *) "#Weights n_i/n (ratio of number of chains) of each polymer found in the topology file in the same order"
    temp_val=sum(hyper_chains(:)%rho /hyper_chains(:)%total_length )
    do i = 1  , ng
        temp_dens=(exp(hyper_chains(i)%mu))*(hyper_chains(i)%Qsum*hyper_chains(i)%distribution)   
        write (333,'(F18.16)') sum(temp_dens(:,:,:,1),mask=saddledensity(:,:,:,1)>=0.05)*dvol/volume
        write (222,'(F18.16)') hyper_chains(i)%rho
        write (111,'(F18.16)') (hyper_chains(i)%rho /hyper_chains(i)%total_length)/temp_val
    end do
 close(unit=111)
 close(unit=222)
 close(unit=333)

!! Assumed to have a canonical ensemble
if (polymer_ensemble .eq. 1) then


  if (arg_val .eq. 1) then
    call system("mkdir dens_xy")
    do i=1,ng
      write (x1,*) i
      output_str0='./dens_xy/xy_dens_iter_mol_id_'//TRIM(ADJUSTL(x1))//'.dat'
        open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
          do y=1,gridy
            do x=1,gridx
            write(111, *) x*dx, '   ' , y*dy, '  ', &
                (volume*hyper_chains(i)%weighted_rho*sum(hyper_chains(i)%distribution(x,y,:,j))*dz,j=1,monomer_types)
            end do
          write(111, *) 
        end do 
      close(unit=111)
    end do
  end if



  if (arg_val .eq. 2) then
    ! Assuming centre of micelle is in the centre of the box
    xc=gridx/2
    yc=gridy/2
    zc=gridz/2

    call system("mkdir volume_fractions_x_and_y")
    !! #1. X distance from centre 2. Y distance from centre 3. Solvophobic vol frac in x dir 4. Solvophilic vol frac in x dir 
     !! 5. Solvophobic vol frac in y dir 6. Solvophilic vol frac in y dir 
    do i=1,ng
      write (x1,*) i
      output_str0='./volume_fractions_x_and_y/mol_id_'//TRIM(ADJUSTL(x1))//'.dat'
        open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
          do j=0,gridx/2
            write(111, *) (j)*dx, '   ' , (j)*dy, '  ', &
             (monomer_volume(k)*volume*hyper_chains(i)%weighted_rho*hyper_chains(i)%distribution(xc+j,yc,zc,k),k=1,monomer_types),&
             (monomer_volume(k)*volume*hyper_chains(i)%weighted_rho*hyper_chains(i)%distribution(xc,yc+j,zc,k),k=1,monomer_types)
          end do
        close(unit=111)
    end do
  end if


end if

      

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   FINISH (CLOSE FILES AND DEALLOCATE ETC.)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  close (unit=30)

    
  deallocate(density,newdensity,dfield,saddledensity,saddlefield,temp_dens)
  deallocate(var,dvar,temp_arr,con_mask,temp_max)



  call finish
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   END OF MAIN PROGRAM. NOW INNER SUBROUTINES
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO READ INPUT PARAMETERS  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine read_input_parameters(file_name)
    
    character*24, intent(in) :: file_name
    integer                  :: io_error

    open(unit=20,file=file_name,status='old',action='read',iostat=io_error)

      if(io_error ==0) then

        read(20,*) sizex
        read(20,*) sizey
        read(20,*) sizez
        read(20,*) gridx
        read(20,*) gridy
        read(20,*) gridz
        read(20,*) threads
        read(20,*)
        read(20,*) phibar
        read(20,*)
        read(20,*) itermax_scf
        read(20,*) every_out
        read(20,*) every_show
        read(20,*) accuracy_scf
        read(20,*) mixing_type
        read(20,*) lambda
        read(20,*) 
        read(20,*) qread  ! qread =  1: read from cnf_densities.in
                          ! qread =  0: read from cnf_fields.in
                          ! qread = -1: setup from program
        read(20,*) 
        read(20,*) geo
        read(20,*) blur
        read(20,*) sig1
        read(20,*) sig2
        read(20,*) init_rr
        read(20,*) seed

        read(20,*) 
        read(20,*) polymer_ensemble
        read(20,*) 
        read(20,*) umbrella
        read(20,*) um_rr
        read(20,*) um_phi
        read(20,*) 
        read(20,*) constraint
        read(20,*) kapcon
        read(20,*) con_rr
        read(20,*) con_shift
        read(20,*) asp_ratio

      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if

    close(unit=20)
    
  
    return

  end subroutine read_input_parameters
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO CALCULATE DENSITIES  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Densities(field,density,hyper_chains,node_calc)
!
   use Propagate

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
    integer                                                     , intent(in)  :: node_calc !1 calculates node distribution, produces large files
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density
    type(hyper_copolymer), dimension(ng), intent(inout)       :: hyper_chains !Creates a number of chain objects


    double precision, dimension(:,:,:,:), allocatable :: expfield

    integer :: monomer_type,i,allocStatus ,th_num

    allocate(expfield(gridx,gridy,gridz,monomer_types),stat=allocStatus)


    call Calculate_Expfield(field,expfield)

    th_num=0

    !$OMP PARALLEL PRIVATE(th_num)  SHARED(expfield,hyper_chains,node_calc)
    !$OMP DO 
    DO i=1,ng
        !$ th_num=omp_get_thread_num()

        call Propagate_Hyper(hyper_chains(i),expfield,th_num,node_calc)

    END DO
    !$OMP END DO 
    !$OMP END PARALLEL


    call Calculate_Macro(hyper_chains,density)

    
    if ( (umbrella .eq. 1) .or. (umbrella .eq. 2) .or. (umbrella .eq. 3).or. (umbrella .eq. 4)) then
        
      where (umbrella_mask(:,:,:,1) .eq. 1) 
        where (density(:,:,:,1) .lt. um_phi/monomer_volume(1)) 
            density(:,:,:,1)=um_phi/monomer_volume(1)
        end where
      end where
        
    end if

    deallocate(expfield)

    return

  end subroutine Calculate_Densities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO CALCULATE ENERGY  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          SUBROUTINE TO CALCULATE ENERGY  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Energy(field,density,hyper_chains)

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: density
    type(hyper_copolymer), dimension(ng), intent(in)                         :: hyper_chains !Creates a number of chain objects




    field_energy = dvol*sum(density*field)
    chains_energy =0.0d0
    total_rho=0.0d0
    do i=1,ng
        chains_energy =chains_energy+ hyper_chains(i)%energy 
        total_rho=total_rho+hyper_chains(i)%rho 
    end do

    energy = chains_energy - field_energy + interaction_energy(density)          


   
!     print*,NEW_LINE('A'), " E " , energy , " IE ", interaction_energy(density),& 
!                             " FE ", field_energy,  "CE ",  chains_energy              

    return

  end subroutine Calculate_Energy

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO CALCULATE HOMOGENEOUS ENERGY  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Calculate_Homogeneous_Energy()
!
    double precision, dimension(gridx,gridy,gridz,monomer_types)        :: field_homo, density_homogeneous,new_density
    double precision, dimension(ng)                                     :: inte ,qsum_homo
    double precision                                                    :: polymer_energy_homogeneous, field_energy_homogeneous
    
    do i=1,monomer_types
        density_homogeneous(:,:,:,i)=phibar*aver_frac(i)
    end do



    call Calculate_Mean_Fields(density_homogeneous,field_homo)

    call Calculate_Densities(field_homo,new_density,hyper_chains,0) !newdensity not needed

    
    do i=1,ng
        hyper_chains(i)%Qsum_homo=hyper_chains(i)%Qsum
        hyper_chains(i)%mu= DLOG(phibar*hyper_chains(i)%weight*volume&
                                                    /(hyper_chains(i)%Qsum_homo)) !For canonical this is updated
!          print*, "Chemical potential of chain", i ,hyper_chains(i)%mu
    end do

    if (polymer_ensemble .eq. 2) then
        call Calculate_Macro(hyper_chains,density_homogeneous) !To update the energies according to mu
    end if

    Call Calculate_Energy(field_homo,density_homogeneous,hyper_chains)
    
    energy_homogeneous=energy

    return

  end subroutine Calculate_Homogeneous_Energy

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO CALCULATE CENTRE OF MASS OF POLYMER AND THE DENSITY ACROSS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Centre(density)

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: density
    double precision, dimension(3)                                           :: cc
    double precision                   :: mass ,xhi,mi,theta,xhi_bar,mi_bar,theta_bar, temp, temp1
    integer                            :: x,y,z ,i

    mass=sum(density)*dvol

    !cx calculation
    xhi_bar=0
    mi_bar=0
    theta_bar=0
    do x=1,gridx
        theta=x/gridx*2*pi
        xhi=cos(theta)
        mi=sin(theta)

        temp=sum(density(x,:,:,:))
        xhi_bar=xhi_bar+xhi*temp*dy*dz
        mi_bar=mi_bar+mi*temp*dy*dz
    end do
    xhi_bar=xhi_bar/mass
    mi_bar=mi_bar/mass
    theta_bar=DATAN(mi_bar/xhi_bar)+pi
    cx=sizex*theta_bar/(2*pi)


    !cy calculation
    xhi_bar=0
    mi_bar=0
    theta_bar=0
    do y=1,gridy
        theta=y/gridy*2*pi
        xhi=cos(theta)
        mi=sin(theta)

        temp=sum(density(:,y,:,:))
        xhi_bar=xhi_bar+xhi*temp*dx*dz
        mi_bar=mi_bar+mi*temp*dx*dz
    end do

    xhi_bar=xhi_bar/mass
    mi_bar=mi_bar/mass
    theta_bar=DATAN(mi_bar/xhi_bar)+pi
    cy=sizey*theta_bar/(2*pi)


    !cz calculation
    xhi_bar=0
    mi_bar=0
    theta_bar=0
    do z=1,gridz
        theta=z/gridz*2*pi
        xhi=cos(theta)
        mi=sin(theta)

        temp=sum(density(:,:,z,:))
        xhi_bar=xhi_bar+xhi*temp*dy*dx
        mi_bar=mi_bar+mi*temp*dy*dx
    end do

    xhi_bar=xhi_bar/mass
    mi_bar=mi_bar/mass

    theta_bar=DATAN(mi_bar/xhi_bar)+pi

    cz=sizez*theta_bar/(2*pi)


    !!Finding the centre indices
    temp1=10**8

    do z=1,gridz
        do y=1,gridy
            do x=1,gridx    
                
                temp=(cx-(x*dx))**2+(cy-(y*dy))**2+(cz-(z*dz))**2             
                
                if (temp .lt. temp1) then
                    xc=x
                    yc=y
                    zc=z
                    temp1=temp
                end if
                    
                
            end do
        end do    
    end do
    
    print*, "centre x,y,z", cx, cy, cz ,xc,yc,zc


    !!!Assuming radial symmetry 
!    open(unit=111,file='data/x_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

!        do x = 1  , gridx

!            write (111, *) cx-(x*dx), density(x,yc,zc,1)*monomer_volume(1)  &
!                                             , density(x,yc,zc,2)*monomer_volume(2)      
!        end do
!    close(unit=111)

!    open(unit=111,file='data/y_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

!        do y = 1  , gridy

!            write (111, *) cy-(y*dy),  density(xc,y,zc,1)*monomer_volume(1)  &
!                                             , density(xc,y,zc,2)*monomer_volume(2)     
!        end do
!    close(unit=111)


!    open(unit=111,file='data/z_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

!        do z = 1  , gridz

!            write (111, *) cz-(z*dz), density(xc,yc,z,1)*monomer_volume(1)  &
!                                             , density(xc,yc,z,2)*monomer_volume(2)      
!        end do
!    close(unit=111)


    return

  end subroutine Calculate_Centre
!


!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO MIX FIELDS IN SCF ITERATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mix_Fields(mixing_type,field,dfield,lambda,iter_status,xgrid,ygrid,zgrid)


     integer, intent(in)          :: mixing_type
     integer, intent(in)          :: xgrid,ygrid,zgrid
     double precision, intent(in) :: lambda
     integer, intent(inout)       :: iter_status
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: dfield
     integer :: x,y,z,n,monomer_type,allocStatus

     integer :: nvar 
    
     nvar=xgrid*ygrid*zgrid*monomer_types


      n=0
      do monomer_type = 1,monomer_types 
        do z = 1,zgrid
         do y = 1,ygrid
          do x = 1,xgrid
            n=n+1
             var(n,0)  = field(x,y,z,monomer_type)
             dvar(n,0) = dfield(x,y,z,monomer_type)
          end do
         end do
        end do
      end do

      mixing: select case (mixing_type)
           case (1)
              call mixing_simple(var,dvar,nvar,lambda,iter_status)
           case (2)
              call mixing_lambda(var,dvar,nvar,lambda,iter_status)
           case (3)
              call mixing_anderson(var,dvar,nvar,lambda,iter_status)
           case default
              print*, 'Invalid mixing option in find_saddle_point!'
              stop
      end select mixing

      n=0
      do monomer_type = 1,monomer_types 
       do z = 1,zgrid
         do y = 1,ygrid
          do x = 1,xgrid
            n=n+1
             field(x,y,z,monomer_type) = var(n,0)  
             dfield(x,y,z,monomer_type) = dvar(n,0) 
          end do
         end do
        end do
      end do

      iter_status = iter_status + 1


    return

    
  end subroutine Mix_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE ACHIEVED ACCURACY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Accuracy(density,newdensity,accuracy)

    integer :: x,y,z
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) ::    density
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: newdensity
    double precision, intent(out)       :: accuracy

    accuracy=sum((density-newdensity)*(density-newdensity))*dvol

    return
  end subroutine check_accuracy

 subroutine Integrate(arr,rr,dr,gridrr,outvalue)

    integer, intent(in)                                  :: gridrr
    double precision, dimension(0:gridrr), intent(in)    :: arr
    double precision, dimension(0:gridrr), intent(in)    :: rr
    double precision, intent(in)                         :: dr
    double precision, intent(out)                        :: outvalue
    integer                                              :: r
    
    
    
    outvalue=0  

  !Simpsons rule
  do r=1,gridrr-1
      if (mod(r,2) .eq. 0 ) then
        outvalue=outvalue+2.0d0*(rr(r))**2*arr(r)
      else
        outvalue=outvalue+4.0d0*(rr(r))**2*arr(r)
      end if
  end do
  
  outvalue=outvalue+(gridrr*dr)**2*arr(gridrr)   
  outvalue=outvalue/3.0d0

  outvalue=4.0*pi*outvalue*dr


  end subroutine Integrate

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!


!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_hyper
