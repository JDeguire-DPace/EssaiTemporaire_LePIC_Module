subroutine read_reactions(sig,sig_Er,sig_list,sig_Eex,ncol_mx,sig_type, &
     npt_mx,rname,col_info,scol_rank,scol_info,sigv_mx,ni0,ntype,n_neu,mpi_rank)
!     ===================================================================
!     VERSION:         0.3
!     LAST MOD:      Nov/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      I am not particularly amateur of goto's,
!                    (they can be quite confusing) but for this 
!                    subroutine I make an exception ...
!
!     NOTES:         (1) Structure of col_info(reaction ind, info)
!                     info= 1 -> n1= # of reactants
!                     info= 2 -> n2= # of byproducts
!                     info= 3 to 2+n1 -> reactant ptypes
!                     info= 2+n1 to 2+n1+n2 -> byproduct ptypes  
!                     info= 2+n1+n2+1 -> reaction type
!                                        1= COLLISION
!                                        2= IONIZATION
!                                        3= EXCITATION
!                                        4= CHARGEEXCHANGE
!                                        5= DISSOCIATION
!                     info= 2+n1+n2+2 -> # of byproduct electrons 
!                     p_ncol is the number of reactions per ptype
!
!                    (2) Structure of sig_list(ptype, icol)
!                     Reactions are listed per incident particles (not 
!                     targets)
!                     icol= cross section number, max is p_ncol(ptype)
!                     sig_list(ptype, icol)= index of reaction in input
!                     file
!
!                    (3) Structure of scol_info(reaction ind, info)
!                     info=1 -> incoming particle type 
!                     info=2 -> reflected particle type (off wall)
!                     info=3 -> RN
!                     info=4 -> RE
!                     p_nscol is the number of wall reactions per ptype
!
!                    (4) Structure of scol_rank(ptype, icol)
!                     Reactions are listed per particle types incident
!                     on wall
!                     icol= number of wall reaction, max is p_nscol(ptype)
!                     scol_rank(ptype, icol)= rank of reaction in list
!                     order is from smallest to largest.
!     -------------------------------------------------------------------
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: mpi_rank,ipt,i,j,jmax,npt_mx,icol,ncol_mx,sig_npt(ncol_mx), &
       sig_type(ncol_mx),sig_list(npart,ncol_mx), &
       ptype,ntype,flag_pa,flag_re,i_re,r_sy,l_sy, &
       col_info(ncol_mx,10),n_re,n_by,ind_label,iscol, &
       scol_rank(npart,ncol_mx),indx(ncol_mx),n_neu
  real(kind=8):: sig_Er(npt_mx),sig(npt_mx,ncol_mx),sig_Eex(ncol_mx,2), &
       scol_info(ncol_mx,4),sort_arr(ncol_mx),sig_scale(2),ni0(npart), &
       sigv_mx(npart,ncol_mx),Ps
  real(kind=8), allocatable:: sig_tmp(:,:,:)
  character:: name*4,rname*20,rlabel*50,tname(2)*6,rtype*20
  parameter ( jmax=10000 )

  allocate ( sig_tmp(npt_mx,ncol_mx,2) )

  !
  ! Open external file
  !
  open(10,file='input_dir/'//rname)
  open(11,file='DATA/particles.out')

  !
  ! Initialize variables & arrays
  !
  do icol=1,ncol_mx
     sig_npt(icol)= 0
     sig_type(icol)= 0
     do ipt=1,npt_mx
        sig_Er(ipt)= 0.d0
        sig(ipt,icol)= 0.d0
        do j=1,2 
           sig_tmp(ipt,icol,j)= 0.d0
        enddo
     enddo
  enddo

  icol=0
  flag_pa=0
  flag_re=0
  col_info=0
  p_ncol=0
  iscol=0
  p_nscol=0
  sort_arr=0.d0
  indx=0
  tag_neg=0
  tag_neu=0
  tag_beam=0
  n_neu=0

  ! Set electron parameters
  pname(1)='[e]'
  charge(1)= -qe
  mass(1)= 9.10938188d-31
  ni0(1)=1.d0
  ntype=1

  !
  ! Scan through file
  !
  do j=1,jmax

19   read(10,*,END=101,ERR=100) name

     !
     !  Find section type
     !

     ! Ions
     if( name.eq.'IONS' .or. name.eq.'Ions' .or. name.eq.'ions' ) then
        flag_pa=1  
        ! Write info on file
        if(mpi_rank.eq.0) then 
           write(11,*) ' '
           write(11,*) 'Ions:'
        endif
        goto 20
     endif

     ! Neutrals
     if( name.eq.'NEUT' .or. name.eq.'Neut' .or. name.eq.'neut' ) then
        flag_pa=1  
        ! Write info on file
        if(mpi_rank.eq.0) then 
           write(11,*) ' '
           write(11,*) 'Neutrals:'
        endif
        goto 20
     endif

     ! Reactions
     if( name.eq.'REAC' .or. name.eq.'Reac' .or. name.eq.'reac' ) then
        if (flag_pa.eq.0) goto 101 ! Error
        ! Add another reaction
        icol= icol + 1
        ! Warning
        if( mpi_rank.eq.0 .and. icol.gt.ncol_mx ) then
           print*, 'icol > ncol_mx, please correct.'
           print*, 'Abort simulation ...'           
           call stop_calculation
        endif
        ! Store collision type
        sig_type(icol)=1 
        goto 21
     endif

     ! Surface reactions
     if( name.eq.'SURF' .or. name.eq.'Surf' .or. name.eq.'surf' ) then
        if (flag_pa.eq.0) goto 101 ! Error
        goto 22
     endif

     goto 19 

     !
     ! Section dedicated to the definition of ions
     ! 

     ! Allows space for some text
20   do while(name.ne.'----')
        read(10,*,END=101,ERR=100) name
     enddo

200  read(10,FMT='(A6)',END=101,ERR=100,ADVANCE='NO') tname(1)

     if(tname(1)(1:2).eq.'--') then 
        ! Warning
        if( mpi_rank.eq.0 .and. ntype.gt.npart ) then
           print*, 'Warning ntype>npart, please correct ...'
           call stop_calculation
        endif
        goto 99 ! jump to next section 
     else
        ntype= ntype + 1 
        pname(ntype)=tname(1)
        if(mpi_rank.eq.0) write(11,*) pname(ntype) ! write particle name on file
        if(pname(ntype).eq.'[H-]') tag_neg= ntype
        if(pname(ntype).eq.'[H]') tag_neu= ntype
        if(pname(ntype).eq.'[eb]') tag_beam= ntype
        read(10,*,END=101,ERR=100) mass(ntype),charge(ntype), &
             Ti(ntype),ni0(ntype)
        if( mpi_rank.eq.0 .and. ni0(ntype).gt.1.d0 ) then
           print*, 'Warning (n/n0)>1 for ', pname(ntype)
           print*, 'Please correct ...'
           call stop_calculation
        endif
        charge(ntype)= charge(ntype)*qe
        mass(ntype)= mass(ntype)*amu
        if(charge(ntype).eq.0.d0) then 
           n_neu= n_neu + 1
           ni0(ntype)= ni0(ntype)*ngas
           ! Neutrals have negative masses
           mass(ntype)= -ABS(mass(ntype))
        endif
        goto 200 ! check for another particle
     endif

     !
     ! Section dedicated to reactions
     !     
21   read(10,*,END=101,ERR=100) rlabel ! Read reaction label

     ! Write info on file
     if(flag_re.eq.0) then
        if(mpi_rank.eq.0) then
           write(11,*) ' '
           write(11,*) 'Reactions:'
        endif
        flag_re=1
     endif

     ! Extract reaction characteristics
     ind_label=2
     l_sy=0
     r_sy=0
     do i_re=1,50 ! Scan reaction label

        ! Extract particles name "[.]"
        if( rlabel(i_re:i_re).eq.'[' ) l_sy=i_re
        if( rlabel(i_re:i_re).eq.']' ) r_sy=i_re
        if( rlabel(i_re:i_re).eq.'>' ) then 
           n_re=ind_label-2
           ! Save # of reactants
           col_info(icol,ind_nre)=n_re
        endif

        ! Find ptype
        if(r_sy.gt.l_sy) then 
           do ptype=1,(ntype+1) ! Loop over particles
              if( mpi_rank.eq.0 .and. ptype.eq.(ntype+1) ) then
                 write(*,500) icol 
500              format(' Warning: unknown particle in reaction #',1x,i2)
                 print*, 'Please correct ...'
                 call stop_calculation
              endif
              if( rlabel(l_sy:r_sy).eq.pname(ptype) ) exit
           enddo
           ! Save reactant/byproduct ptype
           ind_label= ind_label + 1 
           col_info(icol,ind_label)=ptype
           l_sy=0
           r_sy=0
        endif

     enddo
     
     ! Save # of byproducts
     n_by= ind_label - n_re - 2
     col_info(icol,ind_nby)=n_by


     ! Read information about energy exchange between species
     sig_scale=0.d0
     read(10,*,END=101,ERR=100) sig_Eex(icol,ind_Eth), &
          sig_Eex(icol,ind_dE)
     read(10,*,END=101,ERR=100) sig_scale(1),sig_scale(2)
     read(10,*,END=101,ERR=100) rtype
        
     sig_Eex(icol,ind_Eth)= sig_Eex(icol,ind_Eth)*sig_scale(1)

     if( rtype.eq.'COLLISION' .or. rtype.eq.'collision' ) &
          col_info(icol,ind_nby+1+n_re+n_by)=1
     if( rtype.eq.'IONIZATION' .or. rtype.eq.'ionization' ) &
          col_info(icol,ind_nby+1+n_re+n_by)=2
     if( rtype.eq.'EXCITATION' .or. rtype.eq.'excitation' ) &
          col_info(icol,ind_nby+1+n_re+n_by)=3
     if( rtype.eq.'CHARGEEXCHANGE' .or. rtype.eq.'chargeexchange' ) &
          col_info(icol,ind_nby+1+n_re+n_by)=4
     if( rtype.eq.'DISSOCIATION' .or. rtype.eq.'dissociation' ) &
          col_info(icol,ind_nby+1+n_re+n_by)=5
        
     ! Warning
     if( mpi_rank.eq.0 .and. col_info(icol,ind_nby+1+n_re+n_by).eq.0 ) then
        print*, 'Warning: unknown reaction type in reaction #',icol
        print*, 'Please correct ...'
        call stop_calculation
     endif

     ! Store number of byproduct electrons
     do i_re=ind_nby+1+n_re,ind_nby+n_re+n_by
        if(col_info(icol,i_re).eq.1) then
           col_info(icol,ind_nby+1+n_re+n_by+1)= &
                col_info(icol,ind_nby+1+n_re+n_by+1) + 1
        endif
     enddo
        
     ! Allows space for some text
     do while(name.ne.'----')
        read(10,*,END=101,ERR=100) name
     enddo

     if(mpi_rank.eq.0) write(11,'(i3,1x,a30,1x,f8.2)') icol,rlabel,sig_Eex(icol,ind_Eth) ! write reaction on file 

     !
     ! Read data points
     !
     do ipt=1,npt_mx


        ! Warning
        if( mpi_rank.eq.0 .and. ipt.eq.npt_mx ) then
           print*, 'Maximum number of data points in cross section reached'
           print*, 'Abort simulation ...'
           call stop_calculation
        endif
        
        ! Check for end of array string        
        read(10,FMT='(A)',ADVANCE='NO') name
        if(name.eq.'----') then
           read(10,FMT='(A)') name ! reread once more 
           exit ! exit loop
        else
           backspace(10)
        endif

        ! Add another data point
        sig_npt(icol)= sig_npt(icol) + 1

        read(10,*) sig_tmp(ipt,icol,1), sig_tmp(ipt,icol,2)
        sig_tmp(ipt,icol,1)= sig_tmp(ipt,icol,1)*sig_scale(1)
        sig_tmp(ipt,icol,2)= sig_tmp(ipt,icol,2)*sig_scale(2)

     enddo

     goto 99 ! Jump to next section 

     !
     ! Section dedicated to surface reactions
     !     
22   do while(name.ne.'----') ! Allows space for some text
        read(10,*,END=101,ERR=100) name
     enddo

     ! Write info on file
     if(mpi_rank.eq.0) then
        write(11,*) ' '
        write(11,*) 'Surface reactions:'
     endif

     ! Loop over surface reactions
220  iscol= iscol + 1 
     read(10,FMT='(A6)',END=101,ERR=100,ADVANCE='NO') tname(1)

     if(tname(1)(1:2).eq.'--') then 

        ! Sort cross-sections
        nscol= iscol - 1
        call indexx(nscol,sort_arr,indx)

        ! Save surface reaction index by particle type and rank
        do iscol=1,nscol
           ptype= scol_info(indx(iscol),1)
           p_nscol(ptype)= p_nscol(ptype) + 1
           scol_rank(ptype,p_nscol(ptype))= indx(iscol)
        enddo
        
        ! Jump to next section 
        goto 99 

     else

        ! read surface reaction characteristics
        read(10,*,END=101,ERR=100) tname(2),scol_info(iscol,3), & ! RN
             scol_info(iscol,4) ! RE
     
        ! write surface reaction on file
        if(mpi_rank.eq.0) write(11,*) tname(1),tname(2) 

        ! Find particle labels
        do i=1,2 ! loop over incident & byproduct particles
           do ptype=1,(ntype+1) ! loop over particles
              if( mpi_rank.eq.0 .and. ptype.eq.(ntype+1) ) then
                 write(*,501) iscol 
501              format(' Warning: unknown particle found in surface reaction #',1x,i2)
                 print*, 'Please correct ...'
                 call stop_calculation
              endif
              if( tname(i).eq.pname(ptype) ) exit
           enddo
           ! Save particle type (incident=1 and byproduct=2)
           scol_info(iscol,i)= ptype 
        enddo
       
        ! Temporary store RN
        sort_arr(iscol)= scol_info(iscol,3)

        goto 220 ! check for another surface reaction

     endif

99   continue
  enddo

  if(mpi_rank.eq.0) then
     print*, '# of collisions > maximum allowed, please correct'
     print*, 'Abort calculation ...'
     call stop_calculation
  endif
  
100 continue
  !
  ! Error
  !
  if(mpi_rank.eq.0) then
     print*, 'An error occured while reading cross section #',icol
     print*, 'Abort calculation ...'
     call stop_calculation
  endif

101 continue
  !
  ! Warnings & messages
  !
  if( mpi_rank.eq.0 .and. flag_pa.eq.0 ) then
     print*, 'Warning: "PARTICLE" section must be placed first, please correct ...'
     call stop_calculation
  endif

  if(mpi_rank.eq.0) print*, 'Gas chemistry read correctly (end of file reached)'
  ncol=icol

  if(ngas.eq.0.d0) ncol=0
  if(ncol.eq.0) then
     if(mpi_rank.eq.0) print*, 'Collision-less mode is set'
  else
     if(mpi_rank.eq.0) print*, 'Total number of reactions:',ncol
  endif


  !
  ! Calculate null collision frequency & rank of collisions (sorting)
  !
  call ordering(sig,sig_Er,sig_tmp,sig_npt,sig_list,col_info, &
       sig_Eex,sigv_mx,ncol_mx,npt_mx,ntype,mpi_rank)


  !
  ! Deallocate sig_tmp array
  !
  deallocate ( sig_tmp )

  !
  ! Check if wall collision probabilities does not exceed one
  !
  do ptype=1,ntype
     Ps= 0.d0
     do icol=1,p_nscol(ptype)
        Ps= Ps + scol_info(scol_rank(ptype,icol),3)
     enddo
     ! Warning
     if( mpi_rank.eq.0 .and. Ps.gt.1.d0 ) then
        print*, 'Total wall collision probability greater than 1 for ',&
             pname(ptype)
        print*, 'Please correct ...'
        call stop_calculation
     endif
  enddo

  !
  ! Close input file
  !
  close(10)
  close(11)

  return
end subroutine read_reactions


subroutine ordering(sig,sig_Er,sig_tmp,sig_npt,sig_list,col_info, &
     sig_Eex,sigv_mx,ncol_mx,npt_mx,ntype,mpi_rank)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Nov/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:        /
!
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: ipt,ipt_new,npt,npt_mx,icol,ncol_mx,sig_npt(ncol_mx), &
       sig_list(npart,ncol_mx),ntype,ptype,length_sort,mpi_rank
  integer:: col_info(ncol_mx,10),rtype1,rtype2,ind_col(ntype)
  real(kind=8):: sig_tmp(npt_mx,ncol_mx,2),sig_Er(npt_mx), &
       sig(npt_mx,ncol_mx),sigv_mx(npart,ncol_mx),mu,vr,&
       sig_Eex(ncol_mx,2)
  real(kind=8):: Ek_L,Ek_R,sig_L,sig_R,sig_p,Ekp,sort_Ek_sav,Eth
  real(kind=8), allocatable :: sort_Ek(:),sort_Ek_tmp(:)
  integer, allocatable :: indx(:)
  character:: pnum*1

  !
  ! Initialize variables & arrays
  !
  sigv_mx=0.d0
  sig_list=0
  rtype1=0
  rtype2=0


  !
  ! Concatenate and sort kinetic energies in cross-sections 
  !
  length_sort= SUM(sig_npt(1:ncol))
  allocate ( sort_Ek(length_sort), &
       sort_Ek_tmp(length_sort), &
       indx(length_sort) )

  ! Concatenate
  ipt_new=0
  do icol=1,ncol
     do ipt=1,sig_npt(icol)
        ipt_new= ipt_new + 1 
        sort_Ek(ipt_new)= sig_tmp(ipt,icol,1)
     enddo
  enddo

  ! Sort
  call indexx(length_sort,sort_Ek,indx)           

  ipt_new=0
  sort_Ek_sav= 1.d10
  do ipt=1,length_sort
     if( sort_Ek(indx(ipt)).ne.sort_Ek_sav ) then
        ipt_new= ipt_new + 1 
        sort_Ek_tmp(ipt_new)= sort_Ek(indx(ipt))
        sort_Ek_sav= sort_Ek_tmp(ipt_new)
     endif
  enddo

  ! Update number of points in cross-section
  sig_npt_mx= ipt_new

  ! Warning (only to reassure anxious people ...)
  if( mpi_rank.eq.0 .and. sig_npt_mx.ge.npt_mx ) then
     print*, 'sig_npt_mx > npt_mx, please correct'
     print*, 'sig_npt_mx=',sig_npt_mx
     print*, 'Abort calculation ...'
     call stop_calculation
  endif

  ! Array storing kinetic energies
  sig_Er(1:sig_npt_mx)= sort_Ek_tmp(1:sig_npt_mx)

  deallocate( sort_Ek, sort_Ek_tmp )

  !
  ! Extrapolate other cross sections to finest energy grid  
  !
  do icol=1,ncol

     ! Get array length
     npt=sig_npt(icol)
        
     ! Threshold energy
     Eth=sig_Eex(icol,ind_Eth)

     do ipt=1,sig_npt_mx

        ! Get stored kinetic energy
        Ekp= sig_Er(ipt)

        ! Scan energy from cross-section arrays
        do ipt_new=1,npt

           Ek_L= sig_tmp(ipt_new,icol,1)
           Ek_R= sig_tmp(ipt_new+1,icol,1)
           sig_L= sig_tmp(ipt_new,icol,2)
           sig_R= sig_tmp(ipt_new+1,icol,2)           

           if(Ekp.lt.Ek_L) then 
              sig_p= sig_L
              ! Exit loop
              exit
           endif

           if( Ekp.ge.Ek_L .and. Ekp.le.Ek_R ) then 
              ! Perform linear interpolation
              sig_p = sig_L + ( Ekp - Ek_L )* & ! Calculate sigma
                   ( sig_R - sig_L )/( Ek_R - Ek_L )
              ! Exit loop
              exit
           endif

        enddo

        ! Warning 
        if( sig_L.eq.0.d0 .and. sig_R.gt.0.d0) then 
           if(Eth.gt.Ek_R) then
              print'(" Eth greater than threshold energy in reaction #",1x,i3)',icol
              print*, 'Please correct...'
              call stop_calculation
           endif
        endif

        ! Save interpolated cross-section
        sig(ipt,icol)= sig_p

     enddo
  enddo
 
  !
  ! Calculate (sig*v)_max & rank of cross-sections
  !
  ind_col=0
  do icol=1,ncol

     ! List cross sections per incident particles
     rtype1= col_info(icol,ind_nby+1)
     ind_col(rtype1)= ind_col(rtype1) + 1 
     sig_list(rtype1,ind_col(rtype1))= icol

     do ipt=1,sig_npt_mx

        ! Reduced mass
        rtype1= col_info(icol,ind_nby+1)
        rtype2= col_info(icol,ind_nby+2)
        mu=ABS(mass(rtype1))*ABS(mass(rtype2))/ &
             (ABS(mass(rtype1))+ABS(mass(rtype2)))

        ! Calculate relative velocity (energy in eV)
        vr= dsqrt(2.d0*sig_Er(ipt)*qe/mu)
  
        ! Calculate ( sig*v_r )_max
        if(sig_Er(ipt).ge.1.d0) & ! Cut off at 1. eV
             sigv_mx(rtype1,icol)= MAX( sigv_mx(rtype1,icol), vr*sig(ipt,icol) )

     enddo

  enddo

  ! Save number of reactions per particle specie
  p_ncol(1:ntype)= ind_col(1:ntype)


  ! Save cross sections per common reactant
  if(mpi_rank.eq.0) then
     do ptype=1,ntype
        write (pnum,'(i1)'),ptype
        open(10,file='DATA/reactions.'//pnum)
        do icol=1,p_ncol(ptype)
           do ipt=1,sig_npt_mx
              write(10,'(1x,f12.4,1x,es12.4)') sig_Er(ipt), &
                   sig(ipt,sig_list(ptype,icol))  
           enddo
           write(10,*) ' '
        enddo
        close(10)
     enddo
  endif

  return
end subroutine ordering