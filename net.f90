! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module net

      use star_private_def
      use const_def

      implicit none

      private
      public :: set_net, do_net, do1_net, do_micro_change_net, &
         get_screening_mode, default_set_which_rates, default_set_rate_factors, &
         default_set_op_mono_factors


      contains


      subroutine do_net(s, nzlo, nzhi, reuse_given_rates, ierr)
         use net_lib, only: net_work_size
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: reuse_given_rates
         integer, intent(out) :: ierr

         logical, parameter :: use_omp = .true.
         integer :: k, op_err, net_lwork, j, jj, cnt, kmax
         real(dp) :: abs_e, abs_e_dm, abs_e_limit, max_abs_e_dm, dm_limit, e_limit
         integer, pointer :: ks(:)
         logical, parameter :: only_dlnT = .false.
         logical :: okay

         include 'formats'

         ierr = 0

         if (s% eps_nuc_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) then
            do k = nzlo, nzhi
               s% eps_nuc(k) = 0d0
               s% d_epsnuc_dlnd(k) = 0d0
               s% d_epsnuc_dlnT(k) = 0d0
               s% d_epsnuc_dx(:,k) = 0d0
               s% eps_nuc_categories(:,k) = 0d0
               s% dxdt_nuc(:,k) =  0d0
               s% dxdt_dRho(:,k) =  0d0
               s% dxdt_dT(:,k) =  0d0
               s% d_dxdt_dx(:,:,k) =  0d0
               s% eps_nuc_neu_total(k) = 0d0
            end do
            return
         end if

         net_lwork = net_work_size(s% net_handle, ierr)

         if (nzlo == nzhi) then
            call do1_net( &
               s, nzlo, s% species, reuse_given_rates, &
               s% num_reactions, net_lwork, ierr)
            return
         end if

         if (use_omp) then
            okay = .true.
!$OMP PARALLEL DO PRIVATE(k,op_err)
            do k = nzlo, nzhi
               if (.not. okay) cycle
               op_err = 0
               call do1_net( &
                  s, k, s% species, reuse_given_rates, &
                  s% num_reactions, net_lwork, op_err)
               if (op_err /= 0) okay = .false.
            end do
!$OMP END PARALLEL DO
            if (.not. okay) ierr = -1
         else
            do k = nzlo, nzhi
               call do1_net( &
                  s, k, s% species, reuse_given_rates, s% num_reactions, net_lwork, ierr)
               if (ierr /= 0) exit
            end do
         end if

      end subroutine do_net


      subroutine do1_net( &
            s, k, species, reuse_given_rates, num_reactions, net_lwork, ierr)
         use rates_def, only: std_reaction_Qs, std_reaction_neuQs, i_rate
         use net_def, only: Net_Info
         use net_lib, only: net_get
         use chem_def, only: chem_isos, category_name, i_ni56_co56, i_co56_fe56
         use eos_def, only : i_eta
         use utils_lib,only: realloc_double, realloc_double3
         type (star_info), pointer :: s
         integer, intent(in) :: k, species, num_reactions, net_lwork
         logical, intent(in) :: reuse_given_rates
         integer, intent(out) :: ierr

         integer :: i, j, kk, screening_mode
         real(dp) :: log10_rho, log10_T, T, alfa, beta, &
            d_eps_nuc_dRho, d_eps_nuc_dT, cat_factor, tau_gamma
         real(dp), target :: net_work_ary(net_lwork)
         real(dp), pointer :: net_work(:)
         type (Net_Info), target :: net_info_target
         type (Net_Info), pointer :: netinfo

         character (len=100) :: message
!=======Changed======
         real(dp) :: neu_reac_info(4)
         integer :: w
!====================
         real(dp), pointer :: reaction_neuQs(:)
         integer :: sz
         real(dp) :: eps_nuc_factor
         logical :: clipped_T

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         net_work => net_work_ary
         netinfo => net_info_target

         s% eps_nuc(k) = 0d0
         s% d_epsnuc_dlnd(k) = 0d0
         s% d_epsnuc_dlnT(k) = 0d0
         s% d_epsnuc_dx(:,k) = 0d0
         s% eps_nuc_categories(:,k) = 0d0
         s% dxdt_nuc(:,k) =  0d0
         s% dxdt_dRho(:,k) =  0d0
         s% dxdt_dT(:,k) =  0d0
         s% d_dxdt_dx(:,:,k) =  0d0
         s% eps_nuc_neu_total(k) = 0d0

         if ((s% eps_nuc_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) .or. &
              s% abar(k) > s% max_abar_for_burning) then
            return
         end if

         log10_rho = s% lnd(k)/ln10
         log10_T = s% lnT(k)/ln10
         T = s% T(k)
         clipped_T = (s% max_logT_for_net > 0 .and. log10_T > s% max_logT_for_net)
         if (clipped_T) then
            T = s% max_logT_for_net
            log10_T = log10_cr(T)
         end if

         screening_mode = get_screening_mode(s,ierr)
         if (ierr /= 0) then
            write(*,*) 'unknown string for screening_mode: ' // trim(s% screening_mode)
            stop 'do1_net'
            return
         end if

         if (s% reaction_neuQs_factor /= 1d0) then
            sz = size(std_reaction_neuQs,dim=1)
            allocate(reaction_neuQs(sz))
            do j=1,sz
               reaction_neuQs(j) = std_reaction_neuQs(j)*s% reaction_neuQs_factor
            end do
         else
            reaction_neuQs => std_reaction_neuQs
         end if
         
         if (s% use_other_net_get) then
            call s% other_net_get( &
               s% id, k, &
               s% net_handle, .false., netinfo, species, num_reactions, s% xa(1:species,k), &
               T, log10_T, s% rho(k), log10_Rho, &
               s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
               s% eta(k), s% d_eos_dlnd(i_eta,k), s% d_eos_dlnT(i_eta,k), &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, reaction_neuQs, reuse_given_rates, .false., &
               s% eps_nuc(k), d_eps_nuc_dRho, d_eps_nuc_dT, s% d_epsnuc_dx(:,k), &
               s% dxdt_nuc(:,k), s% dxdt_dRho(:,k), s% dxdt_dT(:,k), s% d_dxdt_dx(:,:,k), &
               screening_mode, s% theta_e(k), &
               s% eps_nuc_categories(:,k), &
               s% eps_nuc_neu_total(k), net_lwork, net_work, ierr)
         else
!=======changed below, added zero array and second false as dummy args====
            do w=1, 4
               neu_reac_info(w) = 0d0
            end do
            call net_get( &
               s% net_handle, neu_reac_info, .false., .false., netinfo, &
               species, num_reactions, s% xa(1:species,k), &
               T, log10_T, s% rho(k), log10_Rho, &
               s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
               s% eta(k), s% d_eos_dlnd(i_eta,k), s% d_eos_dlnT(i_eta,k), &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, reaction_neuQs, reuse_given_rates, .false., &
               s% eps_nuc(k), d_eps_nuc_dRho, d_eps_nuc_dT, s% d_epsnuc_dx(:,k), &
               s% dxdt_nuc(:,k), s% dxdt_dRho(:,k), s% dxdt_dT(:,k), s% d_dxdt_dx(:,:,k), &
               screening_mode, s% theta_e(k), &
               s% eps_nuc_categories(:,k), &
               s% eps_nuc_neu_total(k), net_lwork, net_work, ierr)
         end if
     
         if (clipped_T) then
            d_eps_nuc_dT = 0
            s% dxdt_dT(1:species,k) = 0
         end if

         if (s% nonlocal_NiCo_kap_gamma > 0d0 .and. &
               .not. s% nonlocal_NiCo_decay_heat) then
            tau_gamma = 0
            do kk = 1, k
               tau_gamma = tau_gamma + s% dm(kk)/(4*pi*s% rmid(kk)*s% rmid(kk))
            end do
            tau_gamma = tau_gamma*s% nonlocal_NiCo_kap_gamma
            s% eps_nuc(k) = s% eps_nuc(k)*(1d0 - exp_cr(-tau_gamma))
         end if
         
         if (abs(s% eps_nuc(k)) > s% max_abs_eps_nuc) then
            s% eps_nuc(k) = sign(s% max_abs_eps_nuc, s% eps_nuc(k))
            d_eps_nuc_dRho = 0d0
            d_eps_nuc_dT = 0d0
            s% d_epsnuc_dx(:,k) = 0d0
         end if

         if (s% reaction_neuQs_factor /= 1d0) deallocate(reaction_neuQs)

         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*)
               write(*,*) 'do1_net: net_get failure for cell ', k
               !return
               call show_stuff(s,k,net_lwork,net_work)
            end if
            if (is_bad_num(s% eps_nuc(k))) then
               if (s% stop_for_NaNs) then
                  if (is_nan(s% eps_nuc(k))) then
                     write(*,2) 'eps_nuc', k, s% eps_nuc(k)
                     stop 'do1_net'
                  end if
               end if
            end if
            return
         end if

         if (is_bad_num(s% eps_nuc(k))) then
            if (s% stop_for_NaNs) then
               if (is_nan(s% eps_nuc(k))) then
                  write(*,2) 'eps_nuc', k, s% eps_nuc(k)
                  stop 'do1_net'
               end if
            end if
            ierr = -1
            return
         end if

         s% d_epsnuc_dlnd(k) = d_eps_nuc_dRho*s% rho(k)
         s% d_epsnuc_dlnT(k) = d_eps_nuc_dT*s% T(k)

         eps_nuc_factor = s% eps_nuc_factor
         if (eps_nuc_factor /= 1d0) then
            s% eps_nuc(k) = s% eps_nuc(k)*eps_nuc_factor
            s% d_epsnuc_dlnd(k) = s% d_epsnuc_dlnd(k)*eps_nuc_factor
            s% d_epsnuc_dlnT(k) = s% d_epsnuc_dlnT(k)*eps_nuc_factor
            s% d_epsnuc_dx(:,k) = s% d_epsnuc_dx(:,k)*eps_nuc_factor
            s% eps_nuc_categories(:,k) = s% eps_nuc_categories(:,k)*eps_nuc_factor
         end if

         if (s% dxdt_nuc_factor /= 1d0) then
            s% dxdt_nuc(:,k) = s% dxdt_nuc(:,k)*s% dxdt_nuc_factor
            s% dxdt_dRho(:,k) = s% dxdt_dRho(:,k)*s% dxdt_nuc_factor
            s% dxdt_dT(:,k) = s% dxdt_dT(:,k)*s% dxdt_nuc_factor
            s% d_dxdt_dx(:,:,k) = s% d_dxdt_dx(:,:,k)*s% dxdt_nuc_factor
         end if

         if (is_bad_num(s% eps_nuc(k))) then
            write(*,*) 'k', k
            write(*,1) 's% eps_nuc(k)', s% eps_nuc(k)
            ierr = -1
            call show_stuff(s,k,net_lwork,net_work)
            write(*,*) '(is_bad_num(s% eps_nuc(k)))'
            write(*,*) 'failed in do1_net'
            if (s% stop_for_NaNs) then
               if (is_nan(s% eps_nuc(k))) then
                  stop 'do1_net'
               end if
            end if
            return
         end if
         
         if (k == -1) then
            write(*,*)
            do j=1,species
               write(*,3) 'd_epsnuc_dx ' // trim(chem_isos% name(s% chem_id(j))), &
                  s% newton_iter, k, s% d_epsnuc_dx(j, k)
            end do
            call show_stuff(s,k,net_lwork,net_work)
         end if

         if (s% model_number == -1) then
            write(*,5) 'eps_nuc', k, s% newton_iter, s% model_number, s% newton_adjust_iter, &
                        s% eps_nuc(k)
         end if

         !if (k == 3123 .and. s% model_number == 2001 .and. s% newton_iter == 1) then
         if (.false.) then
            write(*,*)
            call show_stuff(s,k,net_lwork,net_work)
            write(*,1) 's% eps_nuc(k)', s% eps_nuc(k)
            write(*,1) 's% d_epsnuc_dlnd(k)', s% d_epsnuc_dlnd(k)
            write(*,1) 's% d_epsnuc_dlnT(k)', s% d_epsnuc_dlnT(k)
            write(*,*)
            write(*,*) 'do1_net'
            stop
            !ierr = -1
         end if

         if (.false.) call show_stuff(s,k,net_lwork,net_work)

      end subroutine do1_net



      subroutine show_stuff(s,k,lwork,work)
         use chem_def
         use rates_def
         use net_lib, only: get_reaction_id_table_ptr, get_net_rate_ptrs
         use num_lib, only: qsort
         type (star_info), pointer :: s
         integer, intent(in) :: k, lwork
         real(dp), pointer :: work(:)

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: i, j, ierr, species, num_reactions
         real(dp) :: log10_Rho, log10_T
         real(dp), pointer :: v(:)
         integer, pointer :: index(:)
         real(dp), pointer, dimension(:) :: &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho

         include 'formats'

         logical, parameter :: do_sort = .true.

         ierr = 0
         species = s% species
         num_reactions = s% num_reactions
         log10_T = s% lnT(k)/ln10
         log10_Rho = s% lnd(k)/ln10

         call get_net_rate_ptrs(s% net_handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_rate_ptrs'
            call mesa_error(__FILE__,__LINE__)
         end if

         call get_reaction_id_table_ptr(s% net_handle, reaction_id, ierr)
         if (ierr /= 0) return

         do i=1,num_reactions
            if (s% rate_factors(i) /= 1d0) then
               write(*,2) 'rate factor ' // trim(reaction_Name(reaction_id(i))), &
                  s% rate_factors(i)
            end if
         end do

         write(*,2) 'k', k
         write(*,*)
         write(*,*) 'net_name ', trim(s% net_name)
         write(*,*) 'species', species
         i = max(species, num_reactions)
         allocate(v(i), index(i))
         write(*,*)
         if (.true.) then
            write(*, *)
            if (do_sort) then
               do j=1,num_reactions
                  v(j) = abs(rate_raw(j))
               end do
               call qsort(index, num_reactions, v)
            else
               do j=1,num_reactions
                  index(j) = j
               end do
            end if

            write(*,*) 'reaction rate_raw'
            do i=1,num_reactions
               j = index(num_reactions+1-i)
               write(*,2) trim(reaction_Name(reaction_id(j))), k, rate_raw(j)
            end do
         end if

         if (.false.) then
            write(*,*)
            write(*,*) 'screened rates'
            do j=1,num_reactions
               write(*,3) 'screened rate ' // trim(reaction_Name(reaction_id(j))), &
                  j, k, rate_screened(j)
            end do
         end if

         if (.true.) then
            write(*,*)
            do j=1,species
               write(*,2) 'dxdt ' // trim(chem_isos% name(s% chem_id(j))), k, s% dxdt_nuc(j, k)
            end do
         end if
         write(*,*)

         if (.false.) then
            write(*,*)
            do j=1,species
               write(*,2) 'dt*dxdt ' // trim(chem_isos% name(s% chem_id(j))), k, &
                  s% dt * s% dxdt_nuc(j, k)
            end do
         end if

         if (.true.) then
            if (do_sort) then
               do j=1,species
                  v(j) = s% xa(j,k)
               end do
               call qsort(index, species, v)
            else
               do j=1,num_reactions
                  index(j) = j
               end do
            end if
            write(*,*)
            do i=1,species
               j = index(species+1-i)
               if (.true. .or. s% xa(j,k) > 1d-9) &
                  write(*,1) 'xin(net_iso(i' // &
                     trim(chem_isos% name(s% chem_id(j))) // '))= ', s% xa(j,k)
            end do
         end if

         write(*,*)
         write(*,1) 'logT =', log10_T
         write(*,1) 'logRho =', log10_Rho
         write(*,1) 'eta =', s% eta(k)
         write(*,*)
         write(*,1) 'T =', s% T(k)
         write(*,1) 'rho =', s% rho(k)
         write(*,1) 'abar =', s% abar(k)
         write(*,1) 'zbar =', s% zbar(k)
         write(*,1) 'z2bar =', s% z2bar(k)
         write(*,1) 'ye =', s% ye(k)
         write(*,*) 'screening_mode = ' // trim(s% screening_mode)
         write(*,1) 'theta_e =', s% theta_e(k)

         return
         !stop 'do1_net'

         if (.false.) then
            write(*,*)
            do j=1,num_categories
               write(*,2) trim(category_name(j)), k, s% eps_nuc_categories(j, k)
            end do
         end if
         if (.true.) then
            write(*, *)
            write(*,*) 'raw rates'
            do j=1,num_reactions
               write(*,2) 'raw rate ' // trim(reaction_Name(reaction_id(j))), &
                  k, rate_raw(j)
            end do
         end if
         if (.false.) then
            write(*, *)
            write(*,*) 'raw rates dlnT'
            do j=1,num_reactions
               write(*,2) 'raw rate dlnT ' // trim(reaction_Name(reaction_id(j))), &
                  k, rate_raw_dT(j)*s% T(k)
            end do
            write(*, *)
         end if



         !return


         if (.false.) then
            write(*, *)
            write(*,*) 'screened rates dlnT'
            do j=1,num_reactions
               write(*,2) 'screened rate dlnT ' // trim(reaction_Name(reaction_id(j))), &
                  k, rate_screened_dT(j)*s% T(k)
            end do
         end if
         if (.false.) then
            write(*, *)
            write(*,*) 'screened rates dlnRho'
            do j=1,num_reactions
               write(*,2) 'screened rate dlnRho ' // trim(reaction_Name(reaction_id(j))), &
                  k, rate_screened_dRho(j)*s% rho(k)
            end do
         end if
         if (.false.) then
            write(*,*)
            do j=1,species
               write(*,2) 'dxdt_dlnRho ' // &
                  trim(chem_isos% name(s% chem_id(j))), k, s% dxdt_dRho(j, k)*s% Rho(k)
            end do
         end if
         if (.false.) then
            write(*,*)
            do j=1,species
               write(*,2) 'dxdt_dlnT ' // &
                  trim(chem_isos% name(s% chem_id(j))), k, s% dxdt_dT(j, k)*s% T(k)
            end do
         end if
         write(*,*) 'X'
         write(*,*)
         write(*,2) 'sum(s% xa(1:species,k))', k, sum(s% xa(1:species,k))
         write(*,2) '1 - sum(s% xa(1:species,k))', k, 1 - sum(s% xa(1:species,k))
         !do j=1,species
         !   write(*,1) trim(chem_isos% name(s% chem_id(j))), s% xa(j,k)
         !end do
         write(*,*)
         write(*,2) 'nnuc = ', species
         write(*,*)
         do j=1,species
            write(*,'(a)') '      j' // trim(chem_isos% name(s% chem_id(j))) // ' = ' // &
               'get_nuclide_index_in_set("' // trim(chem_isos% name(s% chem_id(j))) // '", set)'
         end do
         write(*,*)
         do j=1,species
            write(*,'(a)',advance='no') 'j' // trim(chem_isos% name(s% chem_id(j))) // ', '
         end do
         write(*,*)
         write(*,*)
         do j=1,species
            write(*,1) 'x(j' // trim(chem_isos% name(s% chem_id(j))) // ')= ', s% xa(j,k)
         end do
         write(*,*)
         write(*,*)
         do j=1,species
            write(*,1) 'xin(net_iso(i' // trim(chem_isos% name(s% chem_id(j))) // '))= ', s% xa(j,k)
         end do
         write(*,*)
         write(*,*)
         write(*,1) 'T =', s% T(k)
         write(*,1) 'logT =', log10_T
         write(*,1) 'rho =', s% rho(k)
         write(*,1) 'logRho =', log10_Rho
         write(*,1) 'abar =', s% abar(k)
         write(*,1) 'zbar =', s% zbar(k)
         write(*,1) 'z2bar =', s% z2bar(k)
         write(*,1) 'ye =', s% ye(k)
         write(*,1) 'eta =', s% eta(k)
         write(*,*) 'screening_mode = ' // trim(s% screening_mode)
         write(*,1) 'theta_e =', s% theta_e(k)
         write(*,*)
      end subroutine show_stuff


      subroutine do1_zero_net_vars(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% eps_nuc(k) = 0
         s% eps_nuc_categories(:,k) = 0
         s% d_epsnuc_dlnd(k) = 0
         s% d_epsnuc_dlnT(k) = 0
         s% d_epsnuc_dx(:,k) = 0
         s% eps_nuc_neu_total(k) = 0
         s% dxdt_nuc(:,k) = 0
         s% dxdt_dRho(:,k) = 0
         s% dxdt_dT(:,k) = 0
         s% d_dxdt_dx(:,:,k) = 0
      end subroutine do1_zero_net_vars


      integer function get_screening_mode(s,ierr)
         use rates_lib, only: screening_option
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (s% screening_mode_value >= 0) then
            get_screening_mode = s% screening_mode_value
            return
         end if
         get_screening_mode = screening_option(s% screening_mode, ierr)
         if (ierr /= 0) return
         s% screening_mode_value = get_screening_mode
      end function get_screening_mode


      subroutine do_micro_change_net(s, new_net_name, ierr)
         use net_def
         type (star_info), pointer :: s
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr
         ierr = 0
         s% net_name = new_net_name
         call set_net(s, new_net_name, ierr)
      end subroutine do_micro_change_net


      subroutine set_net(s, new_net_name, ierr)
         use net_lib
         use utils_lib, only: realloc_double
         use alloc, only: update_nvar_allocs, set_chem_names
         use chem_def, only: ih1, ihe4
         use rates_def
         type (star_info), pointer :: s
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr

         integer :: i, ir
         integer :: old_num_reactions, old_nvar_chem, old_species
         integer, parameter :: num_lowT_rates = 10
         integer, pointer :: net_reaction_ptr(:)

         include 'formats'

         old_num_reactions = s% num_reactions

         if (s% net_handle /= 0) call free_net_handle(s% net_handle)

         s% net_handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in alloc_net_handle'
            return
         end if
         call net_ptr(s% net_handle, s% net_rq, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_ptr'
            return
         end if

         call net_tables(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_tables'
            return
         end if

         old_species = s% species
         s% species = net_num_isos(s% net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_num_isos'
            return
         end if

         s% num_reactions = net_num_reactions(s% net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_num_reactions'
            return
         end if

         old_nvar_chem = s% nvar_chem
         s% nvar_chem = s% species
         call update_nvar_allocs(s, s% nvar_hydro, old_nvar_chem, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in update_nvar_allocs'
            return
         end if

         call get_chem_id_table_ptr(s% net_handle, s% chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in get_chem_id_table_ptr'
            return
         end if

         call get_net_iso_table_ptr(s% net_handle, s% net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in get_net_iso_table_ptr'
            return
         end if
         
         if (s% net_iso(ih1) == 0 .or. s% net_iso(ihe4) == 0) then
            write(*,*) 'mesa/star requires both h1 and he4 in net isotopes'
            write(*,*) 'but they are not included in ' // trim(new_net_name)
            ierr = -1
            return
         end if

         if (associated(s% xa_removed)) deallocate(s% xa_removed)
         allocate(s% xa_removed(s% species))

         call set_chem_names(s)

         call s% set_rate_factors(s% id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in s% set_rate_factors'
            return
         end if

         if (associated(s% op_mono_factors)) deallocate(s% op_mono_factors)
         allocate(s% op_mono_factors(s% species))

         call s% set_op_mono_factors(s% id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in s% set_op_mono_factors'
            return
         end if

      end subroutine set_net


      subroutine net_tables(s, ierr)
         use net_lib ! setup net
         use rates_lib
         use rates_def, only: rates_reaction_id_max
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0

         call net_start_def(s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_start_def'
            return
         end if

         if (len_trim(s% net_name) == 0) then
            write(*,*) 'missing net_name -- please set it and try again'
            ierr = -1
            return
         end if

         call read_net_file(s% net_name, s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in read_net_file ' // trim(s% net_name)
            return
         end if

         call net_finish_def(s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_finish_def'
            return
         end if

         if (associated(s% rate_factors)) deallocate(s% rate_factors)
         allocate(s% rate_factors(rates_reaction_id_max))

         call s% set_rate_factors(s% id, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in set_rate_factors'
            return
         end if

         if (associated(s% which_rates)) deallocate(s% which_rates)
         allocate(s% which_rates(rates_reaction_id_max))

         call s% set_which_rates(s% id, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in set_which_rates'
            return
         end if

         call net_set_which_rates(s% net_handle, s% which_rates, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_set_which_rates'
            return
         end if

         call net_set_logTcut(s% net_handle, s% net_logTcut_lo, s% net_logTcut_lim, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_set_logTcut'
            return
         end if

         call net_set_fe56ec_fake_factor( &
            s% net_handle, s% fe56ec_fake_factor, s% min_T_for_fe56ec_fake_factor, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_set_fe56ec_fake_factor'
            return
         end if

         call net_setup_tables( &
            s% net_handle, rates_cache_suffix_for_star, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'failed in net_setup_tables'
            return
         end if

      end subroutine net_tables


      subroutine default_set_which_rates(id, ierr)
         use rates_def, only: rates_NACRE_if_available
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% which_rates(:) = rates_NACRE_if_available
      end subroutine default_set_which_rates


      subroutine default_set_rate_factors(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% rate_factors(:) = 1
      end subroutine default_set_rate_factors


      subroutine default_set_op_mono_factors(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% op_mono_factors(:) = 1
      end subroutine default_set_op_mono_factors


      end module net

