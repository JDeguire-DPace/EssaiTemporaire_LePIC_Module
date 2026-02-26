!
! Counters and seed number
!
  integer:: npart,np_cell,flag_pdf, &
       flag_restart,n_cell,nbak,n_rg, &
       flag_pbc,cnt_seg,tag_neg,num_grd,&
       tag_neu,flag_nmn,flag_heat,flag_inj,&
       opt_inj,flag_B_pos,omp_rank_max,mpi_rank_max,&
       flag_die,flag_thr,ix_pl,iz_pl,ix_thr,flag_convP,&
       ig_die(5),flag_pbcz,tag_beam,every,flag_bak,&
       flag_circxh,flag_ahp,flag_gridB,flag_RFpot
  parameter (npart=10)

!
! Collisions
!
  integer:: ncol,p_ncol(npart),sig_npt_mx,nscol,p_nscol(npart)
  real(kind=8):: np_mx(npart),nu_max(npart),nu_uplim(npart)

!
! Particle arrays & flags
!
  integer:: np_loss,P_w,ind_nre,ind_nby,ind_Eth, &
            ind_dE,cnt_neg
  parameter ( np_loss=1, P_w=2, ind_nre=1, &
       ind_nby=2, ind_Eth=1, ind_dE=2 )
  real(kind=8):: charge(npart),mass(npart),vt0(npart),Ti(npart)

!
! Magnetic field
!
  integer:: nB
  real(kind=8):: B_scale(5),B0(5),dL(5),x0(5),y0(5),z0(5)
  character:: B_file(5)*1,B_name(5)*20,B_info(5)*2

!
! Particle and simulation info
!
  integer:: ixl_pow,ixr_pow,ns_heat,ns_coll,nm_rg,n_holes,flag_grd,ind_g,&
            ptype_pdf,ns_flx,flag_B,nhy,nhz,ns_inj,ixg,ix_PE,igrid_sec,dir_sec,&
            n_cath,plt_src
  parameter ( nm_rg=5 )
  integer:: ixl_rg(nm_rg),ixr_rg(nm_rg),iyl_rg,iyr_rg,izl_rg,izr_rg,cnt_plt
  real(kind=8):: Pabs,n0,ngas,Nm(npart),dt,xmax,ymax,zmax,x_load, &
            nu_h,Lgy,Lgz,Sg,Lhy,Lhz,xg1,Bmax,nudt,R_ahp
  real(kind=8):: k_eps0,jne,RNeta,THm,x_thr,y_thr,tseq,tseq_init,&
  tseq_final,gam_sec,zg_sec(2)
  character:: pname(npart)*6

!
! Common blocks
!
  common /part_info_int/ np_cell,ncol,sig_npt_mx, &
       flag_pdf,n_cell,nbak,flag_restart,ixl_pow,ixr_pow, &
       flag_B,ns_heat,ns_coll,p_ncol,nscol,p_nscol,&
       n_rg,flag_pbc,n_holes,flag_grd,ind_g,nB,&
       cnt_seg,tag_neg,cnt_neg,tag_neu,flag_nmn, &
       num_grd,ptype_pdf,ns_flx,nhy,nhz,ixl_rg,ns_inj,&
       ixr_rg,iyl_rg,iyr_rg,izl_rg,izr_rg,flag_heat,&
       flag_inj,opt_inj,flag_B_pos,ixg,ix_PE,omp_rank_max,mpi_rank_max,&
       flag_die,flag_thr,ix_pl,iz_pl,ix_thr,flag_convP,ig_die,&
       flag_pbcz,tag_beam,every,cnt_plt,flag_bak,igrid_sec,dir_sec,&
       n_cath,flag_circxh,flag_ahp,flag_gridB,plt_src,flag_RFpot
  common /part_info_dp_1/ charge,mass,Ti,vt0,k_eps0,B_scale, &
          B0,dL,x0,y0,z0,jne,RNeta,THm,x_thr,y_thr,tseq,tseq_init,&
       tseq_final,gam_sec,zg_sec
  common /part_info_dp_2/ Pabs,n0,ngas,Nm,dt,xmax,ymax,zmax,x_load, &
       nu_h,Lgy,Lgz,Sg,Lhy,Lhz,xg1,Bmax,nudt,R_ahp
  common /part_info_dp_3/ np_mx,nu_max,nu_uplim
  common /part_info_ch/ pname,B_file,B_name,B_info
