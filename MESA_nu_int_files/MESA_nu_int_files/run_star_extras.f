 ! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
!======Changed======
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         S% other_net_get => neu_other_net_get
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      
      subroutine neu_other_net_get(  &
            id, k, net_handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, theta_e_for_graboske_et_al,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)
         use net_lib, only: net_get
         use net_def, only: Net_Info, Net_General_Info
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell
         integer, intent(in) :: net_handle
         logical, intent(in) :: just_dxdt
         type (Net_Info), pointer:: n
         type (Net_General_Info), pointer:: g
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened ! if true. use given rate_screened
         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after subtract reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos)       
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)            
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
         integer, intent(in) :: screening_mode
         real(dp), intent(in)  :: theta_e_for_graboske_et_al
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         integer, intent(out) :: ierr ! 0 means okay

         type (star_info),pointer :: s
         real(dp) :: neu_reac_info(4)
         real(dp) :: rad, age
         logical :: neu_reac

         neu_reac = .true.

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! OOPS
            write(*,*) 'failed in call star_ptr in neu_other_net_get'
            return
         end if
         rad = s% r(k) * 1d-6
         age = s% star_age * secyer

         s% report_ierr = .true.
         s% stop_for_NaNs = .true.
         call net_neu_reacy(rad, age, g, n, neu_reac_info, ierr)

         call net_get( &
            net_handle, neu_reac_info, just_dxdt, neu_reac, n, num_isos, &
            num_reactions, x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, &
            reuse_rate_screened, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, &
            d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
            screening_mode, theta_e_for_graboske_et_al, &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)

      contains

      subroutine net_neu_reacy(radius, time, g, n, neu_reacy_info, ierr)
         real(dp), intent(in) :: radius, time
         real(dp), intent(out) :: neu_reacy_info(4)
         integer, intent(out) :: ierr
         type (Net_Info), pointer :: n
         type (Net_General_Info), pointer :: g

         real(dp) :: nu_L, nu_bar_L, nu_E, nu_bar_E, &
         mass_diff, rate_p_to_n, rate_n_to_p
         real(dp) :: nu_L_base, nu_bar_L_base, nu_L_0, nu_bar_L_0, &
         nu_L_timescale, nu_bar_L_timescale
         logical :: dbg

         ierr = 0
         dbg = .false.
         mass_diff = 1.293 ! mass difference between prot and neut in MeV

! THESE PARAMETERS ARE BASED ON THE PROTONEUTRON STAR OF THE CCSN MODEL.
! TO GET ACCURATE RESULTS THEY WOULD NEED TO DEPEND ON THE NS EOS.
! CURRENTLY THEY ARE JUST PLACEHOLDERS TO SHOW THE CODE WORKS.

         nu_L_base = s% x_ctrl(1)
         nu_bar_L_base = s% x_ctrl(2)
         nu_L_0 = s% x_ctrl(3)
         nu_bar_L_0 = s% x_ctrl(4)
         nu_L_timescale = s% x_ctrl(5)
         nu_bar_L_timescale = s% x_ctrl(6)

         ! (anti-)neutrino luminosity in 10^51 ergs/sec
         nu_L = nu_L_base+nu_L_0*EXP(-time/nu_L_timescale)
         nu_bar_L = nu_bar_L_base+nu_bar_L_0*EXP(-time/nu_bar_L_timescale)

         ! average (anti-)neutrino energy in MeV
         nu_E = s% x_ctrl(7)
         nu_bar_E = s% x_ctrl(8)

! THIS STUFF CALCULATES RATES OF CHARGED CURRENT NU REACTIONS

         ! The following formulas are Eqs. 65a,b from Qian & Woosley (1996)
         rate_p_to_n = 4.83*nu_bar_L*(nu_bar_E-2*mass_diff+1.2*(mass_diff* &
         mass_diff/nu_bar_E))*(1/(radius*radius))
         rate_n_to_p = 4.83*nu_L*(nu_E+2*mass_diff+1.2*(mass_diff* &
         mass_diff/nu_E))*(1/(radius*radius))

         if (dbg) then
            write(*,*) 'rate_p_to_n is:', rate_p_to_n, 'at radius,k:', radius, k
            write(*,*) 'rate_n_to_p is:', rate_n_to_p, 'at radius,k:', radius, k
            write(*,*) 'nu_L is:', nu_L, 'at time:', time
            write(*,*) 'nu_bar_l is:', nu_bar_L, 'at time:', time
         end if

         ! HERE WE PLACE THE REACTION RATES AND Q'S IN AN ARRAY WHICH WE PASS
         ! DOWN TO THE NEU_REAC_PLACE ROUTINE.
         neu_reacy_info(1) = rate_p_to_n
         neu_reacy_info(2) = nu_E+mass_diff
         neu_reacy_info(3) = rate_n_to_p
         neu_reacy_info(4) = nu_bar_E-mass_diff

      end subroutine net_neu_reacy

      end subroutine neu_other_net_get

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols=0
      end subroutine how_many_extra_history_header_items
      
      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra history header item
      !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items
      
      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info

!==================
      
      end module run_star_extras
