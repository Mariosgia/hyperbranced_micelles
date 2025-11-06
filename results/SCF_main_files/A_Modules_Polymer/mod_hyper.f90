!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR HYPERBRANCHED COPOLYMERS 
!
!           XX.11.2023
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Every polymer molecule type is given a "Fortran type" with its 
!   connectivity being recorded in molecule_data and its prelevance 
!   being recorded in weights. 
!   
!   In "Propagate_Hyper" each block of the polymer is then propagated first  
!   starting with the backward propagator "qqdag" and then forward propagator "qq"
!   The density, and other quanitites are then calculated in "Calculate_Macro"
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module Molecule_Hyper_Copolymer
!
  use Global
  use Propagate
  USE OMP_LIB
!
  implicit none
!


  private
  

  type hyper_copolymer
    integer           :: num_segments
    integer           :: total_length
    double precision  :: weight !Proportion of chains in ensemble n_i/n_tot.
    integer           :: gmax,gmin
    integer           :: label
    integer           :: ensemble         ! canonical:1, grand canonical:2
    double precision  :: energy           ! free energy in SCF field
    double precision  :: rho
    double precision  :: mu
    double precision  :: weighted_rho !rho_i*v^*^
    double precision  :: frac(monomer_types)
    double precision                                                :: Qsum ,Qsum_homo
    integer, dimension(:,:), allocatable                            :: molecule_data !Allocated in mod_hyper_extras.f90
    double precision, dimension(:,:,:,:), allocatable               :: distribution
    REAL*4          , dimension(:,:,:,:), allocatable               :: node_distribution !
    integer         , dimension(monomer_types)                      :: monomer_type

  end type hyper_copolymer


  type propagator

    double precision, dimension(:,:,:,:), allocatable     :: qq
    double precision, dimension(:,:,:,:), allocatable     :: qqdag

  end type propagator
!
  public :: hyper_copolymer
  public :: Initialize_Hyper
  public :: Propagate_Hyper
  public :: Calculate_Macro
  public :: Test_Propagate


contains

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES COPOLYMER
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!


  subroutine Initialize_Hyper(copolymer,label,monomer_type,ensemble)
!
   
   integer, intent(in)                                   :: label, ensemble
   integer, dimension(monomer_types), intent(in)         :: monomer_type

   
!
   type(hyper_copolymer), intent(inout):: copolymer
!
   character*15     :: string_ensemble
   integer          :: iAllocStatus , i ,j , n ,io_error
   double precision :: x

   allocate(copolymer%distribution(gridx,gridy,gridz,monomer_types),stat=iAllocStatus)
!
   if (ensemble.eq.1) then
      string_ensemble = 'canonical'
   else if (ensemble.eq.2) then
      string_ensemble = 'grand canonical'
   else
      print*, 'Invalid Ensemble'
   end if
!
   copolymer%label            = label
   copolymer%monomer_type     = monomer_type
   copolymer%ensemble         = ensemble

   copolymer%weighted_rho     = lll*copolymer%weight !Only for canonical
   

!   print*, "copolymer ensemble is: ", copolymer%ensemble


  end subroutine Initialize_Hyper

 subroutine Test_Propagate()

    implicit none

    double precision, dimension(gridx,gridy,gridz,monomer_types)             :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types)             :: expfield
    double precision, dimension(:,:,:,:), allocatable   :: qq, qqdag
    double precision, dimension(gridx,gridy,gridz)      :: prop
    integer :: iAllocStatus , i, io_error

    integer :: s,s0,mono,block,monomer_type
!
   allocate(qq(gridx,gridy,gridz,0:100),stat=iAllocStatus)
   allocate(qqdag(gridx,gridy,gridz,0:100),stat=iAllocStatus)

    ds=0.01d0


    do i=1,gridx
        field(i,:,:,:)=sin(i*1.0d0)
    end do

    call Calculate_Expfield(field,expfield)

    print*, "Field A", expfield(32,32,32,1) 
    print*, "Field B", expfield(32,32,32,2) 
    print*, sizex,sizey,sizez

    s0 = 0
    prop=1.0d0
    qq(:,:,:,s0)=prop  
   
    
    do i=0,1
      if (i.eq. 0) then
        monomer_type = 1
      else
        monomer_type = 2
      end if
      s0 = s0 + 50*i
      do s = s0+1, s0+50
         call Propagate_Step(prop,expfield(:,:,:,monomer_type),0)
         qq(:,:,:,s)=prop
            print*, "forward prop" , qq(32,32,32,s), "type", monomer_type , "s" , s
!         print*, "SUm prop" , sum( qq(:,:,:,s))*dvol, "block", block
      end do
    end do
!

    s0 = 0
    prop=1.0d0
    qqdag(:,:,:,s0)=prop

    do i=0,1
      if (i.eq. 0) then
        monomer_type = 2
      else
        monomer_type = 1
      end if
      s0 = s0 + 50*i
      do s = s0+1, s0+50
         call Propagate_Step(prop,expfield(:,:,:,monomer_type),0)
         qqdag(:,:,:,s)=prop
            print*, "backward prop" , qqdag(32,32,32,s), "type", monomer_type,"s" , s
!         print*, "SUm prop" , sum( qq(:,:,:,s))*dvol, "block", block
      end do
    end do

    deallocate(qq,qqdag)

  end subroutine Test_Propagate

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE PROPAGATORS FOR EACH HYPER CHAIN 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Propagate_Hyper(copolymer,expfield,th_num,node_calc)
!
    implicit none

    type(hyper_copolymer), intent(inout)  :: copolymer

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: expfield
    integer , intent(in)                                                     :: th_num
    integer , intent(in)                                                     :: node_calc   


    double precision, dimension(:,:,:), allocatable       :: prop 
    integer, dimension(:), allocatable                    :: sel_ids ,tar_ids, src_ids, temp_ids, cp_ids  


    type(propagator), dimension(:), allocatable       :: propagators 

    integer :: iAllocStatus ,io_error

    integer :: s,s0,mono,monomer_type,typee,sd , i, j, k
    integer :: id,segment_length,gen , species,src,tar,tp_id,tp_id0,tp_len,similar_id,cp_id

     !Creating set of propagators for each linear segment in the copolymer
     allocate(propagators(copolymer%num_segments),stat=iAllocStatus) 
     allocate(prop(gridx,gridy,gridz),stat=iAllocStatus)

     !Creating set of qq and qqdag for each linear segment in the copolymer

     do i=1,copolymer%num_segments
        id=copolymer%molecule_data(i,1)
        segment_length=copolymer%molecule_data(id,2)
        allocate(propagators(id)%qq(gridx,gridy,gridz,0:segment_length),stat=iAllocStatus)
        allocate(propagators(id)%qqdag(gridx,gridy,gridz,0:segment_length),stat=iAllocStatus)
     end do

     copolymer%distribution=0.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  QQDAG CALCULATION  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !First we find all the terminal groups and propagate qqdag for segments that have similar id=0 and 
    sel_ids=pack(copolymer%molecule_data(:,1),&
                (copolymer%molecule_data(:,3) .eq. 0) .and. (copolymer%molecule_data(:,8) .eq. 0))                                                                              

    

    do i=1,size(sel_ids)

        id=sel_ids(i) 
        mono=copolymer%molecule_data(id,2) !Segment length

        prop=1.0d0
        propagators(id)%qqdag(:,:,:,0)=prop !Initial condition for terminal groups

        species=copolymer%molecule_data(id,7)

        do sd=1,mono !sd here is actually N-s     
           call Propagate_Step(prop,expfield(:,:,:,species),th_num)
           propagators(id)%qqdag(:,:,:,sd)=prop
        end do
     

    end do

    !Copying qqdags of the rest of the terminal segments 
    cp_ids=pack(copolymer%molecule_data(:,1),&
            (copolymer%molecule_data(:,3) .eq. 0) .and. (copolymer%molecule_data(:,8) .ne. 0)&
             .and. (copolymer%molecule_data(:,9) .ne. 2))
    do j=1,size(cp_ids)
        cp_id=cp_ids(j) 
        similar_id=copolymer%molecule_data(cp_id,8)
        mono=copolymer%molecule_data(cp_id,2) !Segment length
        do sd=0,mono !sd here is actually N-s     
           propagators(cp_id)%qqdag(:,:,:,sd)=propagators(similar_id)%qqdag(:,:,:,sd)
        end do
    end do
    



    !Second we propagate the internal segments from the largest generation up to and including g=gmin doing 
    gen=copolymer%gmax-1
    do while (gen .ge. copolymer%gmin) 

        !Internal segments or stem ids of generation "gen" 
        src_ids=pack(copolymer%molecule_data(:,1),&
        ((copolymer%molecule_data(:,3) .ne. 0) .and. (copolymer%molecule_data(:,6) .eq. gen )&
          .and. (copolymer%molecule_data(:,8) .eq. 0) )) 

        
        do i=1,size(src_ids)

            id=src_ids(i) 
            mono=copolymer%molecule_data(id,2) !Segment length
            species=copolymer%molecule_data(id,7)

            !Setting up initial condition         
            tar=copolymer%molecule_data(id,5) !target node
            !temp_ids need to have the same number as the target node since we are calculating 
            temp_ids=pack(copolymer%molecule_data(:,1),copolymer%molecule_data(:,4) .eq. tar)
            prop=1.0d0

            
            do j=1,size(temp_ids)
                   
                tp_id=temp_ids(j)
                tp_len=copolymer%molecule_data(tp_id,2)
                prop=prop*propagators(tp_id)%qqdag(:,:,:,tp_len)
            end do
            
            propagators(id)%qqdag(:,:,:,0)=prop !Initial condition 


            do sd=1,mono !sd here is actually N-s           
                call Propagate_Step(prop,expfield(:,:,:,species),th_num)
                propagators(id)%qqdag(:,:,:,sd)=prop
            end do

        end do


        !Copying qqdags for the rest of the generation 
        cp_ids=pack(copolymer%molecule_data(:,1),&
        ((copolymer%molecule_data(:,3) .ne. 0) .and. (copolymer%molecule_data(:,6) .eq. gen )&
          .and. (copolymer%molecule_data(:,8) .ne. 0) .and. ((copolymer%molecule_data(:,9) .ne. 2)) )) 
        do j=1,size(cp_ids)
            cp_id=cp_ids(j) 
            similar_id=copolymer%molecule_data(cp_id,8)
            mono=copolymer%molecule_data(cp_id,2) !Segment length
            do sd=0,mono !sd here is actually N-s     
               propagators(cp_id)%qqdag(:,:,:,sd)=propagators(similar_id)%qqdag(:,:,:,sd)
            end do
        end do

        gen=gen-1
    end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   QQ CALCULATION  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !First we propagate qq for g=gmin case
    gen=copolymer%gmin
    temp_ids=pack(copolymer%molecule_data(:,1),copolymer%molecule_data(:,6) .eq. gen)
    id=temp_ids(1)
    species=copolymer%molecule_data(id,7)
    mono=copolymer%molecule_data(id,2) !Segment length
    prop=1.0d0
    propagators(id)%qq(:,:,:,0)=prop !Initial condition for g=gmin

    do s=1,mono      
       call Propagate_Step(prop,expfield(:,:,:,species),th_num)
       propagators(id)%qq(:,:,:,s)=prop
    end do


    !Second we propagate the remaining segments from gen=gmin+1 up to and including g=gmax 
    do while (gen .lt. copolymer%gmax) 

        !Internal segments or stem ids of generation "gen" 
        tar_ids=pack(copolymer%molecule_data(:,1),&
                 ((copolymer%molecule_data(:,6) .eq. gen ) .and. (copolymer%molecule_data(:,3) .ne. 0 )))

        !Propagating "gen+1" segments
        do i=1,size(tar_ids)

            id=tar_ids(i)
            tar=copolymer%molecule_data(id,5) !target node
!            src_ids=pack(copolymer%molecule_data(:,1),&
!                    (copolymer%molecule_data(:,4) .eq. tar)) !Ids of segments connected to target of previous segment
            src_ids=pack(copolymer%molecule_data(:,1),&
                    ((copolymer%molecule_data(:,4) .eq. tar) .and. ((copolymer%molecule_data(:,8) .eq. 0 )&
                    .or. (copolymer%molecule_data(:,9) .ne. 2 )) )) !Ids of segments connected to target of previous segment

        
            do j=1,size(src_ids)

                tp_id0=src_ids(j) 
                mono=copolymer%molecule_data(tp_id0,2) !Segment length
                species=copolymer%molecule_data(tp_id0,7)

                !Setting up initial condition         
                !Finding out which other segments have source nodes the same as tp_id0 but not tp_id0  
                temp_ids=pack(copolymer%molecule_data(:,1),&
                    ((copolymer%molecule_data(:,4) .eq. tar) .and. (copolymer%molecule_data(:,1) .ne. tp_id0)))

                prop=propagators(id)%qq(:,:,:,copolymer%molecule_data(id,2)) !Include incoming segment from lower gen


                do k=1,size(temp_ids)
                       
                    tp_id=temp_ids(k)
                    tp_len=copolymer%molecule_data(tp_id,2)
                    prop=prop*propagators(tp_id)%qqdag(:,:,:,tp_len)

                end do
               

                propagators(tp_id0)%qq(:,:,:,0)=prop !Initial condition 

                do s=1,mono !sd here is actually N-s           
                    call Propagate_Step(prop,expfield(:,:,:,species),th_num)
                    propagators(tp_id0)%qq(:,:,:,s)=prop

                end do

            end do
        end do

        !Copying the rest of the segments qq
        cp_ids=pack(copolymer%molecule_data(:,1),&
                    ((copolymer%molecule_data(:,6) .eq. gen+1 ) .and. (copolymer%molecule_data(:,8) .ne. 0 )&
                    .and. (copolymer%molecule_data(:,9) .ne. 3 ) )) 
        do j=1,size(cp_ids)
            cp_id=cp_ids(j) 
            similar_id=copolymer%molecule_data(cp_id,8)
            mono=copolymer%molecule_data(cp_id,2) !Segment length
            do sd=0,mono !sd here is actually N-s     
               propagators(cp_id)%qq(:,:,:,sd)=propagators(similar_id)%qq(:,:,:,sd)
            end do
        end do




        gen=gen+1
    end do  !End of while loop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   CALCULATION OF QSUM  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gen=copolymer%gmin
temp_ids=pack(copolymer%molecule_data(:,1),copolymer%molecule_data(:,6) .eq. gen)
id=temp_ids(1)
copolymer%Qsum=dvol * sum(propagators(id)%qqdag(:,:,:,copolymer%molecule_data(id,2)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   CALCULATION OF DISTRIBUTIONS    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,copolymer%num_segments
    
    id=i
    mono=copolymer%molecule_data(id,2)
    src=copolymer%molecule_data(id,4)
    tar=copolymer%molecule_data(id,5)
    gen=copolymer%molecule_data(id,6)
    

    species=copolymer%molecule_data(id,7)
    
    !Either this works 
!     copolymer%distribution(:,:,:,species) =copolymer%distribution(:,:,:,species) &
!       +0.5d0*propagators(id)%qq(:,:,:,0)*propagators(id)%qqdag(:,:,:,mono)
!    do s = 1, mono-1
!         copolymer%distribution(:,:,:,species) =copolymer%distribution(:,:,:,species) &
!       +propagators(id)%qq(:,:,:,s)*propagators(id)%qqdag(:,:,:,mono-s)
!    end do

!      copolymer%distribution(:,:,:,species) =copolymer%distribution(:,:,:,species) &
!       +0.5d0*propagators(id)%qq(:,:,:,mono)*propagators(id)%qqdag(:,:,:,0)

    !Or this works 
    do s = 0, mono-1
         copolymer%distribution(:,:,:,species) =copolymer%distribution(:,:,:,species) &
       +propagators(id)%qq(:,:,:,s)*propagators(id)%qqdag(:,:,:,mono-s)
    end do
end do



copolymer%distribution = copolymer%distribution*ds   
copolymer%distribution= copolymer%distribution/copolymer%Qsum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   CALCULATION OF NODE DISTRIBUTION  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ( node_calc .eq. 1) then
    call Calculate_Node(copolymer,propagators)
end if

IF (ALLOCATED(tar_ids)) THEN
    DEALLOCATE(tar_ids)
END IF

IF (ALLOCATED(sel_ids)) THEN
    DEALLOCATE(sel_ids)
END IF

IF (ALLOCATED(src_ids)) THEN
    DEALLOCATE(src_ids)
END IF

deallocate(  temp_ids, prop ,propagators,cp_ids)


return

  end subroutine Propagate_Hyper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  subroutine Calculate_Macro(hyper_chains,density)
    
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
    type(hyper_copolymer), dimension(ng),  intent(inout) :: hyper_chains 
    integer                                              :: i,j

    density=0.0d0
    if (polymer_ensemble .eq. 1) then
        do i=1,ng
!            print*, "PS sum", sum(hyper_chains(i)%distribution)*dvol
            density=density+volume*hyper_chains(i)%weighted_rho*hyper_chains(i)%distribution
        end do
        
        do i=1,ng
              hyper_chains(i)%rho = sum(hyper_chains(i)%weighted_rho*hyper_chains(i)%distribution)*dvol  
              hyper_chains(i)%energy = -volume*phibar*hyper_chains(i)%weight*&
                                        (DLOG( hyper_chains(i)%Qsum / (volume*phibar*hyper_chains(i)%weight))+1.0d0)
              hyper_chains(i)%mu = DLOG(hyper_chains(i)%rho*volume/hyper_chains(i)%Qsum)
        end do
        

    else if (polymer_ensemble .eq. 2) then 
        
        do i=1,ng
            density=density&
                    +(exp(hyper_chains(i)%mu))&
                    *(hyper_chains(i)%Qsum*hyper_chains(i)%distribution)   
        end do
        
        do j=1,monomer_types
        density(:,:,:,j)=density(:,:,:,j)/monomer_volume(j)
        end do

        do i=1,ng
                
              hyper_chains(i)%rho = (phibar*hyper_chains(i)%weight)*(hyper_chains(i)%Qsum/hyper_chains(i)%Qsum_homo)&
                                        *sum(hyper_chains(i)%distribution )*dvol/monomer_volume(1) !Because we assumed same monomer volume for both species

              hyper_chains(i)%energy = -hyper_chains(i)%Qsum*(exp(hyper_chains(i)%mu)) !Because we assumed same monomer volume for both species
                                       
        end do

    end if
    
    
  end subroutine Calculate_Macro



  subroutine Calculate_Node(copolymer,propagators)

    type(hyper_copolymer), intent(inout)                        :: copolymer
    type(propagator), dimension(:), allocatable,  intent(in)     :: propagators 
    
    REAL*4  :: norm
    integer :: iAllocStatus, i , segment_length, gen, id
    integer, dimension(:), allocatable                    ::  temp_ids  !Holds temporary linear segment ids

    allocate(copolymer%node_distribution(gridx,gridy,gridz,0:copolymer%num_segments),stat=iAllocStatus)


    
    gen=copolymer%gmin
    temp_ids=pack(copolymer%molecule_data(:,1),copolymer%molecule_data(:,6) .eq. gen)
    id=temp_ids(1)
    segment_length=copolymer%molecule_data(id,2)
    copolymer%node_distribution(:,:,:,0)=propagators(id)%qq(:,:,:,0)*propagators(id)%qqdag(:,:,:,segment_length) !Setting the first element as the start of gen=gmin segment
    norm=sum(copolymer%node_distribution(:,:,:,0))*dvol
    copolymer%node_distribution(:,:,:,0)=copolymer%node_distribution(:,:,:,0)/norm

    do i=1,copolymer%num_segments
        segment_length=copolymer%molecule_data(i,2)
        copolymer%node_distribution(:,:,:,i)=propagators(i)%qq(:,:,:,segment_length)*propagators(i)%qqdag(:,:,:,0)
        norm=sum(copolymer%node_distribution(:,:,:,i))*dvol
        copolymer%node_distribution(:,:,:,i)=copolymer%node_distribution(:,:,:,i)/norm
    end do 


  end subroutine Calculate_Node



end module Molecule_Hyper_Copolymer
