# Fortran routine headers

- **subroutine** — `subroutine read_Bfield_map(Bi,n,name,namlen,B_dir,scaling,mpi_rank)`  
  *File:* `Bfield.f90`  —  *Line:* 1

- **subroutine** — `subroutine gaussian_Bfield(Bi,n,h,B0,x0,Lx,B_dir)`  
  *File:* `Bfield.f90`  —  *Line:* 44

- **subroutine** — `subroutine PG_current(Bi,n,h,B0,B_dir)`  
  *File:* `Bfield.f90`  —  *Line:* 70

- **subroutine** — `subroutine EE_magnets(Bi,n,h,B0,x0,z0,d)`  
  *File:* `Bfield.f90`  —  *Line:* 116

- **subroutine** — `subroutine Bfield_Hall_thruster(Bi,n,h,B0,x0,a)`  
  *File:* `Bfield.f90`  —  *Line:* 148

- **function** — `function Bx_magnet(x,y,a)`  
  *File:* `Bfield.f90`  —  *Line:* 201

- **function** — `function By_magnet(x,y,a)`  
  *File:* `Bfield.f90`  —  *Line:* 217

- **subroutine** — `subroutine find_B_dir(B_info,B_dir)`  
  *File:* `Bfield.f90`  —  *Line:* 233

- **subroutine** — `subroutine calc_rho(n,np,rhs_dom,bcnd,ntype,nproc,nproc_mpi,mpi_rank)`  
  *File:* `calc_rho.f90`  —  *Line:* 1

- **subroutine** — `subroutine dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)`  
  *File:* `calc_rho.f90`  —  *Line:* 114

- **subroutine** — `subroutine collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex, ncol_mx,npt_mx,cnt_col,P_loss,sour_xy,sour_xz,Plist,pl_max,sigv_mx, col_info,np_red,flag_dead,flag_cex,nproc,np_tot,iseed,nproc_mpi,mpi_rank, ss2D_xy,ss2D_xz,ss2D_yz,flag_diag)`  
  *File:* `collisions.f90`  —  *Line:* 1

- **subroutine** — `subroutine collision_OMP(vxp,n,h,ntype,nmax,sig,sig_Er, sig_list,sig_Eex,ncol_mx,npt_mx,cnt_col,P_loss,sour_xy,sour_xz, col_info,np_red,flag_dead,nproc,np_tot,iseed,Plist,pl_max,np_add,Nc, ptype,sav_ttype,iproc,np_tot_rg,nu_max_OMP,err_coll,flag_cex, ss2D_xy,ss2D_xz,ss2D_yz,flag_diag)`  
  *File:* `collisions.f90`  —  *Line:* 202

- **subroutine** — `SUBROUTINE scatter(vx1,vy1,vz1,vx,vy,vz,costheta,phi)`  
  *File:* `collisions.f90`  —  *Line:* 1010

- **subroutine** — `subroutine calc_Efield(n,h,phi,Ei,bcnd)`  
  *File:* `Efield.f90`  —  *Line:* 1

- **subroutine** — `subroutine generate_boundary(u,n,h,bcnd,V,ngrid,dtype,xl_pow,xr_pow,mpi_rank)`  
  *File:* `generate_boundary.f90`  —  *Line:* 1

- **subroutine** — `SUBROUTINE indexx(n,arr,indx)`  
  *File:* `indexx.f90`  —  *Line:* 1

- **subroutine** — `subroutine load_part(n,h,bcnd,np,vxp,ntype,nmax,kq, ni0,np_tot,nproc,iseed,sum_dEk,Nh,mpi_rank,nproc_mpi)`  
  *File:* `load_part.f90`  —  *Line:* 1

- **subroutine** — `subroutine load_part_OMP(n,h,bcnd,np,vxp,ntype,nmax,kq, ni0,np_tot,nproc,iseed,sum_dEk,Nh,nproc_mpi,iproc)`  
  *File:* `load_part.f90`  —  *Line:* 60

- **subroutine** — `subroutine load_gauss(vx,vy,vt,rnd)`  
  *File:* `load_part.f90`  —  *Line:* 235

- **subroutine** — `subroutine introduction`  
  *File:* `main.f90`  —  *Line:* 1542

- **subroutine** — `subroutine mg(u,b,bcnd,h,res,n,omega,ng,eps,k,ktot,rank,nproc, e2,e4,e8,e16,e32,e64,e128,e256,e512,e1024,e2048, r2,r4,r8,r16,r32,r64,r128,r256,r512,r1024,r2048, bcnd2,bcnd4,bcnd8,bcnd16,bcnd32,bcnd64,bcnd128, bcnd256,bcnd512,bcnd1024,bcnd2048)`  
  *File:* `mg.f90`  —  *Line:* 1

- **subroutine** — `subroutine restriction(u,b,h,bcnd,bcnd2h,r2h,e2h,n,n3,rank,nproc)`  
  *File:* `mg.f90`  —  *Line:* 248

- **subroutine** — `subroutine getres(uorg,rhsorg,h,ir,jr,kr,n,r,flag)`  
  *File:* `mg.f90`  —  *Line:* 808

- **subroutine** — `subroutine prolongation(u,e2h,bcnd,n,n3,nproc)`  
  *File:* `mg.f90`  —  *Line:* 849

- **subroutine** — `subroutine part_mover(n,h,Ei,Bi,p_mac,P_loss,vxp, bcnd,nmax,ntype,ngrid,flag_dead,nproc,np_tot, iproc,ptype,flag_cex,cnt_cex,cnt_dead,beam_div,sum_q_xz, sum_q_yz,dtype,iseed,n_B,h_B)`  
  *File:* `part_expmover.f90`  —  *Line:* 1

- **subroutine** — `subroutine eheating(h,vxp,nmax,ntype,nproc,iseed,np_tot,vt,iproc,P_loss)`  
  *File:* `part_expmover.f90`  —  *Line:* 471

- **subroutine** — `subroutine charge_deposition(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk, Nh,iproc,ptype)`  
  *File:* `part_expmover.f90`  —  *Line:* 557

- **subroutine** — `subroutine part_flux_injection(istep,vxp,n,h,ntype,nmax,bcnd,sour_xy, sour_xz,sour_fx_yz,nproc,iseed,flag_cex,np_tot,N_flx,nproc_mpi, mpi_rank,nsav,P_loss)`  
  *File:* `part_flux_injection.f90`  —  *Line:* 1

- **subroutine** — `subroutine load_flux_OMP(vxp,n,h,ntype,nmax,bcnd,sour_xy,sour_xz,sour_fx_yz, nproc,iseed,flag_cex,np_tot,Nh,ptype,iproc,P_loss)`  
  *File:* `part_flux_injection.f90`  —  *Line:* 105

- **subroutine** — `subroutine part_injection(n,h,bcnd,vxp,sour_xy,sour_xz,sour_yz,ntype,nmax, I_inj,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,np_tot,nproc,iseed, nproc_mpi,N_inj,iproc,P_loss,phi,ni0)`  
  *File:* `part_injection.f90`  —  *Line:* 1

- **subroutine** — `subroutine shifted_maxwellian_flux(v,vb,vt,fmax,iseed)`  
  *File:* `part_injection.f90`  —  *Line:* 267

- **subroutine** — `subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)`  
  *File:* `pdesolver.f90`  —  *Line:* 1

- **function** — `function ran2(irand)`  
  *File:* `ran2.f90`  —  *Line:* 1

- **subroutine** — `subroutine read_input(n,tmax,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow, nsav,eps,omega,kt,rname,ngrid,ng,xl_rg,xr_rg,yl_rg,yr_rg,zl_rg,zr_rg, I_inj,Ca,mpi_rank,flag_avg3D,np_dup,n_B,phi0_RF,f0_RF, phi1_RF,f1_RF)`  
  *File:* `read_input.f90`  —  *Line:* 1

- **subroutine** — `subroutine read_reactions(sig,sig_Er,sig_list,sig_Eex,ncol_mx,sig_type, npt_mx,rname,col_info,scol_rank,scol_info,sigv_mx,ni0,ntype,n_neu,mpi_rank)`  
  *File:* `read_reactions.f90`  —  *Line:* 1

- **subroutine** — `subroutine ordering(sig,sig_Er,sig_tmp,sig_npt,sig_list,col_info, sig_Eex,sigv_mx,ncol_mx,npt_mx,ntype,mpi_rank)`  
  *File:* `read_reactions.f90`  —  *Line:* 477

- **subroutine** — `subroutine restart(n,h,np,vxp,nmax,ntype,kq,time,nproc,iseed, np_tot,sum_dEk,Nh,flag_cex,mpi_rank,nproc_mpi,np_dup)`  
  *File:* `restart.f90`  —  *Line:* 1

- **subroutine** — `subroutine restart_OMP(n,h,np,vxp,nmax,ptype,ntype,kq, iproc,nproc,np_tot,sum_dEk,Nh)`  
  *File:* `restart.f90`  —  *Line:* 195

- **subroutine** — `subroutine sor_rb(u,b,h,bcnd,res,n,n3,omega,eps,ksor,kmg,dig,rank,nproc)`  
  *File:* `sors.f90`  —  *Line:* 1

- **subroutine** — `subroutine jacobi(u,b,h,bcnd,res,n,n3,eps,ksor,kmg,dig,rank,nproc)`  
  *File:* `sors.f90`  —  *Line:* 538

- **subroutine** — `subroutine part_sorting_OMP(vxp,sorting,nmax,ntype, n,h,Plist,pl_max,nproc,np_tot,iproc,ptype)`  
  *File:* `sorting.f90`  —  *Line:* 1

- **subroutine** — `subroutine part_sorting(vxp,sorting,nmax,ntype,n,h,Plist_thd,pl_max, nproc,np_tot,ptype)`  
  *File:* `sorting.f90`  —  *Line:* 97

- **subroutine** — `subroutine part_sorting_dom(istep,vxp,nmax,ntype,n,h,pl_max,nproc,np_tot)`  
  *File:* `unused.f90`  —  *Line:* 1

- **subroutine** — `subroutine part_moments(n,h,vxp,nmax,ntype,kq,nproc,np_tot, iproc,ptype,p_mts_xy,p_mts_xz,p_mts_yz,n_mts)`  
  *File:* `utils.f90`  —  *Line:* 1

- **subroutine** — `subroutine calc_avg(n,h,np,p_mts_xy,p_mts_xz,p_mts_yz,data_pavg_xy, data_pavg_xz,data_pavg_yz,sour_xy,sour_xz,sour_yz,ntype,n_mts, ss2D_xy,ss2D_xz,ss2D_yz,nproc,nproc_mpi)`  
  *File:* `utils.f90`  —  *Line:* 226

- **function** — `FUNCTION MSTIMER()`  
  *File:* `utils.f90`  —  *Line:* 503

- **subroutine** — `subroutine stop_calculation`  
  *File:* `utils.f90`  —  *Line:* 529

- **subroutine** — `subroutine write_data(it,time,n,h,p_mac,P_loss,phi_avg_xy, phi_avg_xz,phi_avg_yz,data_pavg_xy,data_pavg_xz,data_pavg_yz, ntype,ngrid,cnt_col,ncol_mx,sig_list,nproc,cnt_avg,mpi_rank, nproc_mpi,Vgrd,avg3D,flag_avg3D,I_inj,nEf)`  
  *File:* `write_data.f90`  —  *Line:* 1

