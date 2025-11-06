!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM 
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
  use Molecule_Hyper_Copolymer
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
  double precision, dimension(:,:,:,:), allocatable ::saddledensity,saddlefield
  integer, dimension(:,:,:,:), allocatable          :: umbrella_mask
  double precision, dimension(:,:),allocatable :: var,dvar


  logical          :: qexist
  integer          :: qread,every_show
  integer          :: allocStatus,init_status
  character*24     :: file_name ,x1
  character*80     :: output_str0 ,output_str1
  character(len=8) :: fmtout ! format descriptor
  integer          :: polymer_label, monomer_type_A, monomer_type_B ,grad_data
  double precision :: dx, dy, dz, energy_homogeneous

  integer          :: io_error, x,y,z,iter,iter_status, i ,u,con_temp

  double precision ::  total_rho , dummy ,temp_val,radx,rady,radz



  double precision ::cnx,cny,cnz,dis ,  phiori

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
   allocate(con_mask(gridx,gridy,gridz),stat=allocStatus) 
   allocate(con_add(gridx,gridy,gridz),stat=allocStatus) 
   allocate(rad(gridx,gridy,gridz),stat=allocStatus) 

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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Calculating ds for all chains !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  average_chain_length=sum(hyper_chains(:)%weight*hyper_chains(:)%total_length)
  ds=1.0d0/(average_chain_length)
  print*, "ds" , ds
    

  aver_frac(1)=0.0d0
  aver_frac(1)=sum(hyper_chains(:)%weight*hyper_chains(:)%total_length*hyper_chains(:)%frac(1))/average_chain_length
  aver_frac(2)=1.0d0-aver_frac(1) 
  print* , aver_frac(1) ,aver_frac(2)
  lll=phibar/(monomer_volume(1))

  call Initialize




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! MAKING SURE OPENMP IS USED !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$ compiled_with_openmp = .true.
  if (compiled_with_openmp) then
  write(*,*) 'OpenMP used...'
    
  else
  write(*,*) 'Openmp not used...'
  end if

  if (compiled_with_openmp) then
    if ( (mod(ng,threads)) .ne. 0) then
        print*
        print*, "It is better if the number of molecules is a multiple of the number of threads for faster implementation "
        print*
    end if
    if ( threads .gt. ng) then
        print*
        print*, "Warning using more threads than there are polymers"
        print*
        stop
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


    !We propagate once for the homogeneous case to find Qsum_homo which is used for both canonical and grand







! ------ INITIALIZE INTERACTIONS
!


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

  else if (qread.eq.0) then
!
    open (unit=34, file='cnf_fields.in', status='old', action='read', iostat=io_error)
    if(io_error /=0) then                   
      write(*,*) 'Error', io_error , ' while trying to open cnf_fields.in'
    end if
    read (34,*) field
    print*, "IO_error" ,io_error
    close (unit=34)
!
  else if (qread.eq.-1) then
    density = 0.0d0

    call Initialize_Density(geo,dx,dy,dz,init_status,density)
    
    
    if (init_status .eq. 1 ) then
        print*, "Failed to initialize correctly"

        open(unit=63,file='data/details.dat', status='replace',action='write',iostat=io_error)
        write(63,"(a,F10.6)") '# 1. Status 2. Accuracy 3. Energy 4. Concentration  5. Aver. vol. frac.'
        write(63,*)   "3"  , "1000000"  , "1000000"  , "1000000"  , "1000000"
        close(unit=63) 
 
        stop
    end if

    

    

   else
     print*, 'Invalid value of qread'
     stop
!
  end if

  if (qread.ne.0) then
     call Calculate_Mean_Fields(density,field)
     
      
    call Check_Volume_Fraction(density(:,:,:,:))
    print*, " Initial Vol. fraction. average. : " , average_volume_fraction
 
  end if
!
call system("mkdir data")
open(unit=111,file='data/init_volume_fractions.dat', status='replace',action='write',iostat=io_error)
do z = 1, gridz
 do y = 1, gridy
   do x = 1 , gridx 
     write (111, *) x*dx, ' ', y*dy, ' ', z*dz,' ', density(x,y,z,1)  &
                                                    , ' ', density(x,y,z,2)
   end do
 end do
end do

 close (unit=111)
 

call system("mkdir OUTPUT_DENS")


  con_temp=constraint

!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	   INIT CONTRAINTS WHICH USE INITIAL DENSITY TO SET PARAMETERS
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (constraint .ne. 0) then
  !Creating array with radius valuse from the center
    do z = 1, gridz
      do y = 1, gridy
        do x = 1 , gridx 
           rad(x,y,z)=( ((x-gridx/2)*dx)**2 +((y-gridy/2)*dy)**2+((z-gridz/2)*dz)**2  )**0.5d0
        end do
      end do
    end do

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
!!	   CALCULATING THE CHEMICAL POTENTIALS (IN ENSEMBLE 2) AND HOMOGENEOUS ENERGY
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
constraint=0 !Need to calculate mus without the constraint potential
call Calculate_Homogeneous_Energy()
constraint=con_temp







!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	   IF UMBRELLA IS SET TO !=0 THEN THE VOLUME FRACTION IN CERTAIN AREAS ARE SET TO A MINIMUM
!!     OR ABOVE SO THAT CONVERGENCE IS ASSISTED. ONLY USED TO OBTAIN AN INITIAL STATE.
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (umbrella .ge. 1) then
    
    umbrella_mask=0
    cnx=sizex/2.0
    cny=sizey/2.0
    cnz=sizez/2.0

     do z=1,gridz
       do y=1,gridy
         do x=1,gridx
            dis=(cnx-x*dx)**2+(cny-y*dy)**2+(cnz-z*dz)**2

            if (umbrella .eq. 3) dis=(cnx-x*dx)**2+(cny-y*dy)**2
            if (umbrella .eq. 4) dis=(cnx-x*dx)**2
            dis=dis**0.5d0
            if (( dis .le. um_rr) .and. ( umbrella .eq. 1)) then !Micelle core

             umbrella_mask(x,y,z,1)=1
!             print*, "x y z ", x*dx,y*dy,z*dz 
    
            end if

            if (( dis .gt. um_rr) .and. ( dis .lt. 1.5*um_rr) .and. ( umbrella .eq. 2)) then !Vesicle

             umbrella_mask(x,y,z,1)=1
    
            end if

            if ( ( dis .lt. um_rr).and.( umbrella .eq. 3)) then !Cylindrical micelle to be used in 2d

             umbrella_mask(x,y,z,1)=1
    
            end if

            if ( ( dis .lt. um_rr).and.( umbrella .eq. 4)) then !Lammela micelle to be used in 1d

             umbrella_mask(x,y,z,1)=1
    
            end if

         end do
       end do
    end do

    print*
    print*, "Total points used for umbrella mask : " ,sum(umbrella_mask)
    print*
end if

!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	   RUNNING ITERATION LOOP
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!


  iter_status = 1

  do iter = 1,itermax_scf
    

    call Calculate_Densities(field,newdensity,hyper_chains,0)

    call Calculate_Mean_Fields(newdensity,newfield)


    dfield = newfield - field
    call  Check_Accuracy(density,newdensity,accuracy)
!
    call Check_Volume_Fraction(newdensity(:,:,:,:))

    call Calculate_Energy(field,density,hyper_chains)                                        ! FOR TESTING

!    call cpu_time(fin)

    if (((iter/every_show)*every_show) .eq. iter) then                         ! FOR TESTING
     WRITE(*, '(A, I7, A, F17.7,A, F17.7,A, F17.15,A, F12.9,A, F12.9,A, F12.9,A, F16.8)') &
             " Iter ", iter, "  E  ", energy, " HE ", energy_homogeneous ," Accu. ", accuracy,"  rho  ", total_rho ,  &
             "  V1  ",aver_vol_fra(1),  "  V2  " ,aver_vol_fra(2),&
             "  Aver. Chem  ", sum(hyper_chains(:)%mu*hyper_chains(:)%weight)      

    end if                                                      ! FOR TESTING                                         ! FOR TESTING

    

    !OUTPUT FILES FOR TESTING
    if (((iter/every_out)*every_out).eq. iter) then                         ! FOR TESTING

      write (x1,*) iter

        if ( ((iter/(5000))*(5000) ).eq. iter) then
    
            output_str0='./OUTPUT_DENS/densities_iter_'//TRIM(ADJUSTL(x1))//'.dat'
            open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
                 write(111, *) density
            close(unit=111)
           
        end if

        output_str0='./OUTPUT_DENS/xy_dens_iter_'//TRIM(ADJUSTL(x1))//'.dat'
        open(unit=111,file=output_str0, status='replace',action='write',iostat=io_error)
        do y=1,gridy
          do x=1,gridx
          write(111, *) x*dx, '   ' , y*dy, '  ', (sum(density(x,y,:,i))*dz,i=1,monomer_types)
          end do
        write(111, *) 
       end do 
       close(unit=111)



    end if   
    
    if ( .not.(accuracy>accuracy_scf) .and. (iter .ge. 50) ) then     !double negation to deal with NaNs
      exit
    end if

   
    if ( (total_rho .ge. 0.3) .and. (iter .ge. 100)) then !COULD BE REMOVED BUT HELPED CUT TIME FOR UNECCESSARY STATES
        print*, " rho is bigger than 0.3 "
        stat=0
        stop
    end if





    Call Mix_Fields(mixing_type,field,dfield,lambda,iter_status,gridx,gridy,gridz)
    density=newdensity
!
  end do


!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!
!!	                OUTPUT FILES
!!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!


  iterations = iter

  density = newdensity
  call Calculate_Energy(field,density,hyper_chains)

  saddledensity = newdensity
  saddlefield   = field

  if (.not.(accuracy < accuracy_scf)) then
    stat=1 !Did not converge
  else
    stat=0 ! Converged
    
  end if
  
  call Check_Volume_Fraction(saddledensity(:,:,:,:)) 


 if (polymer_ensemble .eq. 1) then
        dummy=phibar*volume
  else if ( polymer_ensemble .ge. 2) then
        dummy=volume
  end if 


  if (constraint .eq. 0) then 
    con_ene=0 
  else if ((constraint .eq. 1) .or. (constraint .eq. 3)) then
      
      temp_val=sum(con_mask*saddledensity(:,:,:,1))*dvol*monomer_volume(1)
      con_ene=0.5d0*kapcon*(temp_val-con_phi)**2

  end if

  open(unit=63,file='data/details.dat', status='replace',action='write',iostat=io_error)

    write(63,"(a,F10.6)") '# 1. Status 2. Accuracy 3. Energy 4. Int_ene 5. Polymers ene 6.Field ene '// &
                            '7. Energy/volume 8. Concentration  9. Sizex 10. Volume 11.Ensemble' // &
                    '12.Homo_Energy 13.V1 14.V2 15.NE 16.Aver. mu 17. Phibar 18. Edif 19. Con_ene'
    write(63,*)   stat  , accuracy  , energy  , interaction_energy(saddledensity), chains_energy, &
                    field_energy, energy/volume, total_rho  , sizex, volume, polymer_ensemble, energy_homogeneous, &
                    aver_vol_fra(1), aver_vol_fra(2), energy/(dummy) , sum(hyper_chains(:)%mu*hyper_chains(:)%weight)&
                    , phibar, (energy-energy_homogeneous)/dummy , con_ene
  close(unit=63)    



  if (.not.(accuracy < accuracy_scf)) then  ! DOUBLE NEGATION DEALS WITH NaNs
    print*, 'Caution: SCF saddle point not found within ',iterations,' iteration steps!'
    print*, 'Accuracy: ', accuracy
    write(30,*) '   Saddle point not found after', iterations, 'iterations'
    write(30,*) '   Accuracy: ', accuracy, ' Energy: ', energy, &
                    ' C: ', total_rho
                     

    goto 100
  end if

  write(30,*) '   Saddle point found after', iterations, 'iterations'
  write(30,*) '   Accuracy: ', accuracy, ' Energy: ', energy, &
           ' C: ', total_rho
  write(30,*) 
!

  print*, '   Saddle point found after', iterations, 'iterations'
  print*, '   Accuracy: ', accuracy, ' Energy: ', energy, &
           ' C: ', total_rho
 

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF SADDLE POINT RESULTS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
100 continue
!
      temp_val=sum(hyper_chains(:)%rho /hyper_chains(:)%total_length )
      open(unit=63,file='data/chain_details.dat', status='replace',action='write',iostat=io_error)
          do i=1,ng
          write(63, '(I7,F25.16,F25.16,F25.16)') i, hyper_chains(i)%mu , hyper_chains(i)%frac(1) ,& 
                          (hyper_chains(i)%rho /hyper_chains(i)%total_length)/temp_val
          end do
      close(unit=63)

!      call Calculate_Densities(saddlefield,newdensity,hyper_chains,1)

!      call system("mkdir OUTPUT_NODE_DISTRIBUTIONS")

!      do i=1,ng
!          write (x1,*) i
!          output_str0='./OUTPUT_NODE_DISTRIBUTIONS/node_distributions_mol_id_'//TRIM(ADJUSTL(x1))//'.out'

!            open(newunit=u,file=output_str0,access='stream',status='REPLACE',form='unformatted')       

!            WRITE(u) gridx, gridy, gridz, hyper_chains(i)%num_segments+1, hyper_chains(i)%node_distribution


!            CLOSE(u)

!      end do







!    write (x1,*) 1
!    output_str0='./OUTPUT_NODE_DISTRIBUTIONS/node_distributions_mol_id_'//TRIM(ADJUSTL(x1))//'.out'
!    open(newunit=u,file=output_str0,access='stream',status='REPLACE',form='unformatted')
!    write(u) hyper_chains(1)%node_distribution
!    close(u)





      call Calculate_Centre(saddledensity)

      open(unit=63,file='data/volume_fractions.dat', status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
             do x=1,gridx
              write(63, *) x*dx, '   ', y*dy, '   ' , z*dz , '   ' , saddledensity(x,y,z,1)*monomer_volume(1), &
                                            '   '  , saddledensity(x,y,z,2)*monomer_volume(2)
             end do
           end do
        end do
      close(unit=63)

    open(unit=111,file='data/final_volume_fractions_phirs.dat', status='replace',action='write',iostat=io_error)
        do x = gridx/2  , gridx

        write (111, *) abs(sizex/2.0d0-(x/dble(gridx)*sizex)), saddledensity(x,gridy/2,gridz/2,1)*monomer_volume(1)  &
                                             , saddledensity(x,gridy/2,gridz/2,2)*monomer_volume(2)      
        end do
!

      open(unit=33, file='data/cnf_densities.out', status = 'replace', action='write', iostat=io_error)
        write (33,*) saddledensity
      close (unit=33)
     ! open(unit=34, file='cnf_fields.out', status = 'replace', action='write', iostat=io_error)
      !  write (34,*) saddlefield
     ! close (unit=34)
!   
     open(unit=63,file='data/dens_xy.dat', status='replace',action='write',iostat=io_error)
    do y=1,gridy
       do x=1,gridx
          write(63, *) x*dx, '   ' , y*dy, '  ', (sum(saddledensity(x,y,:,i))*dz,i=1,monomer_types)
       end do
       write(63, *) 
    end do
    close(unit=63)
    !
    open(unit=63,file='data/dens_xz.dat', status='replace',action='write',iostat=io_error)
    do z=1,gridz
       do x=1,gridx
          write(63, *) x*dx, '   ' , z*dz, (sum(saddledensity(x,:,z,i))*dy,i=1,monomer_types)
       end do
       write(63, *) 
    end do
    close(unit=63)
    !
    open(unit=63,file='data/dens_yz.dat', status='replace',action='write',iostat=io_error)
    do z=1,gridz
       do y=1,gridy
          write(63, *) y*dy, '   ' , z*dz, (sum(saddledensity(:,y,z,i))*dx,i=1,monomer_types)
       end do
       write(63, *) 
    end do
    close(unit=63)

      

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   FINISH (CLOSE FILES AND DEALLOCATE ETC.)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  close (unit=30)

    
  deallocate(density,newdensity,dfield,saddledensity,saddlefield)
  deallocate(var,dvar,rad,con_mask,con_add)



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
!

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
         ! print*, "Chemical potential of chain", i ,hyper_chains(i)%mu
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
    double precision                   :: mass ,xhi,mi,theta,xhi_bar,mi_bar,theta_bar, temp, temp1, cx,cy,cz
    integer                            :: x,y,z ,i,xc,yc,zc

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
    if (gridz .eq. 1) then
      cz=1
      zc=1
    end if
    if (gridy .eq. 1) then
      cy=1
      yc=1
    end if

!    cx=sizex/2.0
    print*, "centre x,y,z", cx, cy, cz ,xc,yc,zc
    !!!Assuming radial symmetry 
    open(unit=111,file='data/x_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

        do x = 1  , gridx

            write (111, *) cx-(x*dx), density(x,yc,zc,1)*monomer_volume(1)  &
                                             , density(x,yc,zc,2)*monomer_volume(2)      
        end do
    close(unit=111)

    if (gridy .gt. 1) then
      open(unit=111,file='data/y_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

          do y = 1  , gridy

              write (111, *) cy-(y*dy),  density(xc,y,zc,1)*monomer_volume(1)  &
                                               , density(xc,y,zc,2)*monomer_volume(2)     
          end do
      close(unit=111)
    end if
    if (gridz .gt. 1) then
      open(unit=111,file='data/z_centered_radial_volume_fractions.dat', status='replace',action='write',iostat=io_error)

          do z = 1  , gridz

              write (111, *) cz-(z*dz), density(xc,yc,z,1)*monomer_volume(1)  &
                                               , density(xc,yc,z,2)*monomer_volume(2)      
          end do
      close(unit=111)
    end if

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
