! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module net_eval
      
      use const_def
      use crlibm_lib, only: log10_cr, exp10_cr, exp_cr
      use chem_def
      use chem_lib, only: get_mass_excess
      use net_def, only: Net_General_Info, Net_Info
      
      implicit none
      
      
      contains
      
!========Changed below, added radius and neu_reac arg================
      subroutine eval_net( & 
            neu_reac_info, n, g, rates_only, just_dxdt, neu_reac, &
            num_isos, num_reactions, num_wk_reactions, &
            x, atemp, alogtemp, arho, alogrho, &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, theta_e_for_graboske_et_al, &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, ierr)
         use net_initialize, only: &
            set_rate_ptrs, setup_net_info, set_ptrs_for_approx21
         use net_approx21, only: approx21_nrat
         use net_approx21_plus_co56, only: approx21_plus_co56_nrat
         use net_approx21_plus_fe53_fe55_co56, only: &
            approx21_plus_fe53_fe55_co56_nrat
         use net_screen
         use net_derivs
!========Changed================
         use net_neu_reac
!===============================

         type (Net_Info), pointer:: n
         type (Net_General_Info), pointer :: g
         logical, intent(in) :: rates_only, just_dxdt
!========Changed================
         logical, intent(in) :: neu_reac
!         real(dp), intent(in) :: radius, age
         real(dp), intent(in) :: neu_reac_info(4)
!===============================
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions, num_wk_reactions
         real(dp), intent(in)  :: x(:)
         real(dp), intent(in)  :: atemp, alogtemp
         real(dp), intent(in)  :: arho, alogrho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! eta and derivatives
         real(dp), intent(in) :: rate_factors(:)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened
         real(dp), intent(out) :: eps_nuc ! ergs/gram/second from burning 
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) 
         real(dp), intent(inout) :: dxdt(:)
         real(dp), intent(inout) :: d_dxdt_dRho(:)
         real(dp), intent(inout) :: d_dxdt_dT(:)
         real(dp), intent(inout) :: d_dxdt_dx(:,:)
         real(dp), intent(inout) :: eps_nuc_categories(:)
         real(dp), intent(out) :: eps_neu_total
         integer, intent(in) :: screening_mode
         real(dp), intent(in)  :: theta_e_for_graboske_et_al
         integer, intent(in) :: lwork
         real(dp), pointer :: work(:) ! (lwork)
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs ! ignore if null
         logical, pointer :: from_weaklib(:) ! ignore if null
         logical, intent(in) :: symbolic
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         real(dp), target, dimension(num_rvs, num_isos) :: screen_h1, screen_he4
         integer, parameter :: max_z_for_cache = 14
         real(dp), target :: graboske_cache(3, max_z_for_cache, max_z_for_cache)
         real(qp), target :: dydt_a(num_rvs*num_isos)
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(dp), target :: mion_a(num_isos)
         real(dp), pointer :: mion(:)
         real(dp) :: enuc, temp, logtemp, T9, rho, logrho, total, prev, curr, prev_T
         real(dp) :: btemp, bden, eps_total, Ys, sum_dxdt, compare, Z_plus_N
         real(qp) :: eps_nuc_MeV(num_rvs)
         integer :: ci, i, j, ir, weak_id, h1, iwork, approx21_num_rates
         integer, pointer :: chem_id(:)
         integer :: time0, time1
         logical :: doing_timing
         
         ! for approx21
         real(dp), pointer :: dfdy(:,:)
         real(dp), dimension(:), pointer :: &
            dratdumdy1, dratdumdy2, d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho
         real(dp) :: &
            deps_total_dRho, deps_total_dT, &
            deps_neu_dT, deps_neu_dRho, fII
               
         real(dp) :: mev2gr
         
         logical, parameter :: dbg = .false.
         !logical, parameter :: dbg = .true.
         
         include 'formats.dek'

         if (dbg) write(*,*) 'enter eval_net'
         
         doing_timing = g% doing_timing
         if (doing_timing) then
            call system_clock(time0)
            g% doing_timing = .false.
         end if

         ierr = 0
         
         dydt(1:num_rvs,1:num_isos) => dydt_a(1:num_rvs*num_isos)
         chem_id => g% chem_id

         eps_nuc = 0

         temp = atemp; logtemp = alogtemp; rho = arho; logrho = alogrho
         call get_T_rho_args(temp, logtemp, rho, logrho, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_T_rho_args'
            return
         end if
         
         if (logtemp < rattab_tlo) then ! clip to table so can eval beta decays
            logtemp = rattab_tlo
            temp = exp10_cr(logtemp)
         end if
         T9 = temp*1d-9
         
         n% screen_h1 => screen_h1
         n% screen_he4 => screen_he4
         if (.not. reuse_rate_screened) then 
            screen_h1(:,:) = 0
            screen_he4(:,:) = 0
         end if
         
         n% graboske_cache => graboske_cache
         n% reaction_Qs => reaction_Qs
         n% reaction_neuQs => reaction_neuQs
         n% eps_neu_total = 0
         n% weak_rate_factor = weak_rate_factor
         
         if (dbg) write(*,*) 'call set_rate_ptrs'
         call set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            iwork, ierr) ! iwork is number of entries in work used for rates
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in set_ptrs_in_work'
            return
         end if

         if (dbg) write(*,*) 'call setup_net_info'
         call setup_net_info( &
            g, n, eps_nuc_categories,  &
            screening_mode, theta_e_for_graboske_et_al, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            reuse_rate_raw, reuse_rate_screened, &
            iwork, ierr) ! iwork updated for amount now used in work
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in setup_net_info'
            return
         end if

         if (g% doing_approx21) then
            if (g% add_fe53_fe55_co56_to_approx21) then
               approx21_num_rates = approx21_plus_fe53_fe55_co56_nrat
            else if (g% add_co56_to_approx21) then
               approx21_num_rates = approx21_plus_co56_nrat
            else
               approx21_num_rates = approx21_nrat
            end if
         else
            approx21_num_rates = -1
         end if
         
         if (g% doing_approx21) then
            call set_ptrs_for_approx21( &
               g% add_fe53_fe55_co56_to_approx21, g% add_co56_to_approx21, &
               iwork, work, dfdy, dratdumdy1, dratdumdy2, &
               d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho)
            mion => mion_a
            mev2gr = 1d6*ev2erg/(clight*clight)
            do i=1,num_isos
                mion(i) = get_mass_excess(chem_isos,g% chem_id(i))*mev2gr
            end do
         end if
         
         if (.not. g% net_has_been_defined) then
            ierr = -1
            if (dbg) write(*,*) 'failed (.not. g% net_has_been_defined)'
            return
         end if
         
         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            time0 = time1
         end if

         if (dbg) write(*,*) 'call set_molar_abundances'
         call set_molar_abundances(g, num_isos, x, n% y, dbg, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in set_molar_abundances'
            return
         end if
         
         if (num_wk_reactions > 0 .and. (.not. reuse_rate_screened)) then
            if ((.not. reuse_rate_raw) .or. &
                 T9 > T9_weaklib_full_off_hi_Z .or. &
                 T9 > T9_weaklib_full_off) then
               ! at high T, weak rates depend on Ye so need to redo
               ! even if okay to reuse raw rates for strong reactions.
               if (dbg) write(*,*) 'call get_weaklib_rates'
               call get_weaklib_rates(ierr)
               if (ierr /= 0) then
                  if (dbg) write(*,*) 'failed in get_weaklib_rates'
                  return
               end if
            end if
         end if
         
         if (associated(actual_Qs) .and. associated(actual_neuQs)) then
            do i = 1, g% num_reactions
               ir = g% reaction_id(i)
               from_weaklib(i) = .false.
               actual_Qs(i) = n% reaction_Qs(ir)
               actual_neuQs(i) = n% reaction_neuQs(ir)
               weak_id = g% weak_reaction_index(i)
               if (weak_id > 0) then
                  if (g% weaklib_ids(weak_id) > 0) then
                     from_weaklib(i) = .true.
                     actual_Qs(i) = n% Q(weak_id)
                     actual_neuQs(i) = n% Qneu(weak_id)
                  end if
               end if
            end do
         end if

         n% d_eps_nuc_dy(:) = 0

         ! limit range of temperatures and densities
         btemp = min(rattab_temp_hi, max(temp, rattab_temp_lo))
         bden = min(1d11, max(rho, 1d-10))
         
         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            time0 = time1
         end if

         if (dbg) write(*,*) 'call get_rates_with_screening'
         call get_rates_with_screening(ierr)
         if (dbg) write(*,*) 'done get_rates_with_screening'
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_rates_with_screening'
            return
         end if
         
!==========Changed==========         
         if (neu_reac) then
            call net_neu_place(neu_reac_info, g, n, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in net_neu_place in net_eval'
               return
            end if
         end if
!==========================

         if (rates_only) return

         d_eps_nuc_dT = 0
         d_eps_nuc_dRho = 0
         d_eps_nuc_dx(:) = 0
         dxdt(:) = 0
         d_dxdt_dRho(:) = 0
         d_dxdt_dT(:) = 0
         if (.not. just_dxdt) d_dxdt_dx(:,:) = 0
         eps_nuc_categories(:) = 0
         eps_neu_total = 0
         
         if (g% doing_approx21) then
            if (g% add_fe53_fe55_co56_to_approx21) then
               call eval_net_approx21_plus_fe53_fe55_co56(ierr)
            else if (g% add_co56_to_approx21) then
               call eval_net_approx21_plus_co56(ierr)
            else
               call eval_net_approx21(ierr)
            end if
            if (ierr /= 0) return                
            return            
         end if         
     
         if (dbg) write(*,*) 'call get_derivs'
         call get_derivs(  &
             n, dydt, eps_nuc_MeV, eta, ye, &
             logtemp, btemp, bden, abar, zbar,  &
             reuse_rate_screened, num_reactions, rate_factors, &
             symbolic, just_dxdt, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_derivs'
            return
         end if
                  
         if (symbolic) then
            do j=1, num_isos
               do i=1, num_isos
                  d_dxdt_dx(i,j) = n% d_dydt_dy(i,j)
               end do
            end do
            return
         end if
         
         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_derivs = g% clock_net_derivs + (time1 - time0)
            time0 = time1
         end if

         ! convert the eps_nuc_categories
         do i=1,num_categories
            n% eps_nuc_categories(i) = Qconv*n% eps_nuc_categories(i)
         end do

         ! store the results
         do i=1,num_isos
            ci = chem_id(i)
            dxdt(i) = chem_isos% Z_plus_N(ci)*dydt(i_rate, i)
         end do
         
         if (.not. just_dxdt) call store_partials
   
         eps_nuc = eps_nuc_MeV(i_rate)*Qconv 
         d_eps_nuc_dT = eps_nuc_MeV(i_rate_dT)*Qconv 
         d_eps_nuc_dRho = eps_nuc_MeV(i_rate_dRho)*Qconv 
         
         eps_neu_total = n% eps_neu_total*Qconv

         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            g% doing_timing = .true.
         end if
         
         contains

         subroutine eval_net_approx21(ierr)
            use net_approx21, only: &
               approx21_special_reactions, approx21_dydt, approx21_d_epsneu_dy, &
               approx21_eval_PPII_fraction, approx21_eps_info, &
               approx21_dfdy, approx21_dfdT_dfdRho
            integer, intent(out) :: ierr
#include "net_eval_approx21_proc.inc"
         end subroutine eval_net_approx21

         subroutine eval_net_approx21_plus_co56(ierr)
            use net_approx21_plus_co56, only: &
               approx21_special_reactions, approx21_dydt, approx21_d_epsneu_dy, &
               approx21_eval_PPII_fraction, approx21_eps_info, &
               approx21_dfdy, approx21_dfdT_dfdRho
            integer, intent(out) :: ierr
#include "net_eval_approx21_proc.inc"
         end subroutine eval_net_approx21_plus_co56

         subroutine eval_net_approx21_plus_fe53_fe55_co56(ierr)
            use net_approx21_plus_fe53_fe55_co56, only: &
               approx21_special_reactions, approx21_dydt, approx21_d_epsneu_dy, &
               approx21_eval_PPII_fraction, approx21_eps_info, &
               approx21_dfdy, approx21_dfdT_dfdRho
            integer, intent(out) :: ierr
#include "net_eval_approx21_proc.inc"
         end subroutine eval_net_approx21_plus_fe53_fe55_co56

         subroutine get_approx21_eps_info( &
               dydt1, rate_screened, do_eps_nuc_categories, eps_total, eps_neu_total, ierr)
            use net_approx21, only: approx21_eps_info
            real(dp), intent(in), dimension(:) :: dydt1, rate_screened
            logical, intent(in) :: do_eps_nuc_categories
            real(dp), intent(out) :: eps_total, eps_neu_total
            integer, intent(out) :: ierr
            call approx21_eps_info( &
               n% y, mion, dydt1, rate_screened, fII, &               
               reaction_Qs(irpp_to_he3), reaction_neuQs(irpp_to_he3), & 
               reaction_Qs(ir_he3_he3_to_h1_h1_he4), &
               reaction_Qs(ir34_pp2), reaction_neuQs(ir34_pp2), & 
               reaction_Qs(ir34_pp3), reaction_neuQs(ir34_pp3), & 
               reaction_Qs(irc12_to_n14), reaction_neuQs(irc12_to_n14), & 
               reaction_Qs(irn14_to_c12), reaction_neuQs(irn14_to_c12), & 
               reaction_Qs(iro16_to_n14), reaction_neuQs(iro16_to_n14), & 
               
               reaction_Qs(irprot_to_neut), reaction_neuQs(irprot_to_neut), & 
               reaction_Qs(irneut_to_prot), reaction_neuQs(irneut_to_prot), & 
               reaction_Qs(irni56ec_to_co56), reaction_neuQs(irni56ec_to_co56), & 
               reaction_Qs(irco56ec_to_fe56), reaction_neuQs(irco56ec_to_fe56), & 
               
               reaction_Qs(irn14ag_lite), &
               reaction_Qs(ir_he4_he4_he4_to_c12), &
               reaction_Qs(ir1212), &
               reaction_Qs(ir1216_to_mg24), reaction_Qs(ir1216_to_si28), &
               reaction_Qs(ir1616a), reaction_Qs(ir1616g), &
               reaction_Qs(ir_ne20_ag_mg24), &
               reaction_Qs(ir_mg24_ag_si28), &
               reaction_Qs(ir_si28_ag_s32), &
               reaction_Qs(ir_s32_ag_ar36), &
               reaction_Qs(ir_ar36_ag_ca40), &
               reaction_Qs(ir_ca40_ag_ti44), &
               reaction_Qs(ir_ti44_ag_cr48), &
               reaction_Qs(ir_cr48_ag_fe52), &
               reaction_Qs(ir_fe52_ag_ni56), &       
               reaction_Qs(ir_fe52_ng_fe53), &       
               reaction_Qs(ir_fe53_ng_fe54), &       
               reaction_Qs(ir_fe54_ng_fe55), &       
               reaction_Qs(ir_fe55_ng_fe56), &                              
               reaction_Qs(irfe52neut_to_fe54), &               
               reaction_Qs(irfe52aprot_to_fe54), &               
               reaction_Qs(irfe54ng_to_fe56), &               
               reaction_Qs(irfe54aprot_to_fe56), &               
               reaction_Qs(irfe52aprot_to_ni56), &               
               reaction_Qs(irfe54prot_to_ni56), &
!=========Changed==============
               reaction_Qs(ir_prot_wk_neut), &
               reaction_Qs(ir_neut_wk_prot), &
!==============================
               eps_total, eps_neu_total, do_eps_nuc_categories, &
               n% eps_nuc_categories(ipp), &
               n% eps_nuc_categories(icno), &
               n% eps_nuc_categories(i_burn_n), &
               n% eps_nuc_categories(i3alf), &
               n% eps_nuc_categories(icc), &
               n% eps_nuc_categories(ico), &
               n% eps_nuc_categories(ioo), &
               n% eps_nuc_categories(i_burn_ne), &
               n% eps_nuc_categories(i_burn_mg), &
               n% eps_nuc_categories(i_burn_si), &
               n% eps_nuc_categories(i_burn_s), &
               n% eps_nuc_categories(i_burn_ar), &
               n% eps_nuc_categories(i_burn_ca), &
               n% eps_nuc_categories(i_burn_ti), &
               n% eps_nuc_categories(i_burn_cr), &
               n% eps_nuc_categories(i_burn_fe), &
               n% eps_nuc_categories(iphoto), &
               n% eps_nuc_categories(i_ni56_co56), &
               n% eps_nuc_categories(i_co56_fe56), &
               .false., ierr)
         end subroutine get_approx21_eps_info
            
         subroutine get_approx21_plus_co56_eps_info( &
               dydt1, rate_screened, do_eps_nuc_categories, eps_total, eps_neu_total, ierr)
            use net_approx21_plus_co56, only: approx21_eps_info
            real(dp), intent(in), dimension(:) :: dydt1, rate_screened
            logical, intent(in) :: do_eps_nuc_categories
            real(dp), intent(out) :: eps_total, eps_neu_total
            integer, intent(out) :: ierr
            
            call approx21_eps_info( &
               n% y, mion, dydt1, rate_screened, fII, &               
               reaction_Qs(irpp_to_he3), reaction_neuQs(irpp_to_he3), & 
               reaction_Qs(ir_he3_he3_to_h1_h1_he4), &
               reaction_Qs(ir34_pp2), reaction_neuQs(ir34_pp2), & 
               reaction_Qs(ir34_pp3), reaction_neuQs(ir34_pp3), & 
               reaction_Qs(irc12_to_n14), reaction_neuQs(irc12_to_n14), & 
               reaction_Qs(irn14_to_c12), reaction_neuQs(irn14_to_c12), & 
               reaction_Qs(iro16_to_n14), reaction_neuQs(iro16_to_n14), & 
               
               reaction_Qs(irprot_to_neut), reaction_neuQs(irprot_to_neut), & 
               reaction_Qs(irneut_to_prot), reaction_neuQs(irneut_to_prot), & 
               reaction_Qs(irni56ec_to_co56), reaction_neuQs(irni56ec_to_co56), & 
               reaction_Qs(irco56ec_to_fe56), reaction_neuQs(irco56ec_to_fe56), & 
               
               reaction_Qs(irn14ag_lite), &
               reaction_Qs(ir_he4_he4_he4_to_c12), &
               reaction_Qs(ir1212), &
               reaction_Qs(ir1216_to_mg24), reaction_Qs(ir1216_to_si28), &
               reaction_Qs(ir1616a), reaction_Qs(ir1616g), &
               reaction_Qs(ir_ne20_ag_mg24), &
               reaction_Qs(ir_mg24_ag_si28), &
               reaction_Qs(ir_si28_ag_s32), &
               reaction_Qs(ir_s32_ag_ar36), &
               reaction_Qs(ir_ar36_ag_ca40), &
               reaction_Qs(ir_ca40_ag_ti44), &
               reaction_Qs(ir_ti44_ag_cr48), &
               reaction_Qs(ir_cr48_ag_fe52), &
               reaction_Qs(ir_fe52_ag_ni56), &               
               reaction_Qs(ir_fe52_ng_fe53), &       
               reaction_Qs(ir_fe53_ng_fe54), &       
               reaction_Qs(ir_fe54_ng_fe55), &       
               reaction_Qs(ir_fe55_ng_fe56), &                              
               reaction_Qs(irfe52neut_to_fe54), &               
               reaction_Qs(irfe52aprot_to_fe54), &               
               reaction_Qs(irfe54ng_to_fe56), &               
               reaction_Qs(irfe54aprot_to_fe56), &               
               reaction_Qs(irfe52aprot_to_ni56), &               
               reaction_Qs(irfe54prot_to_ni56), &
!======Changed=====
               reaction_Qs(ir_prot_wk_neut), &
               reaction_Qs(ir_neut_wk_prot), &
!==================
               eps_total, eps_neu_total, do_eps_nuc_categories, &
               n% eps_nuc_categories(ipp), &
               n% eps_nuc_categories(icno), &
               n% eps_nuc_categories(i_burn_n), &
               n% eps_nuc_categories(i3alf), &
               n% eps_nuc_categories(icc), &
               n% eps_nuc_categories(ico), &
               n% eps_nuc_categories(ioo), &
               n% eps_nuc_categories(i_burn_ne), &
               n% eps_nuc_categories(i_burn_mg), &
               n% eps_nuc_categories(i_burn_si), &
               n% eps_nuc_categories(i_burn_s), &
               n% eps_nuc_categories(i_burn_ar), &
               n% eps_nuc_categories(i_burn_ca), &
               n% eps_nuc_categories(i_burn_ti), &
               n% eps_nuc_categories(i_burn_cr), &
               n% eps_nuc_categories(i_burn_fe), &
               n% eps_nuc_categories(iphoto), &
               n% eps_nuc_categories(i_ni56_co56), &
               n% eps_nuc_categories(i_co56_fe56), &
               .false., ierr)
         end subroutine get_approx21_plus_co56_eps_info
            
         subroutine get_approx21_plus_fe53_fe55_co56_eps_info( &
               dydt1, rate_screened, do_eps_nuc_categories, eps_total, eps_neu_total, ierr)
            use net_approx21_plus_fe53_fe55_co56, only: approx21_eps_info
            real(dp), intent(in), dimension(:) :: dydt1, rate_screened
            logical, intent(in) :: do_eps_nuc_categories
            real(dp), intent(out) :: eps_total, eps_neu_total
            integer, intent(out) :: ierr
            call approx21_eps_info( &
               n% y, mion, dydt1, rate_screened, fII, &               
               reaction_Qs(irpp_to_he3), reaction_neuQs(irpp_to_he3), & 
               reaction_Qs(ir_he3_he3_to_h1_h1_he4), &
               reaction_Qs(ir34_pp2), reaction_neuQs(ir34_pp2), & 
               reaction_Qs(ir34_pp3), reaction_neuQs(ir34_pp3), & 
               reaction_Qs(irc12_to_n14), reaction_neuQs(irc12_to_n14), & 
               reaction_Qs(irn14_to_c12), reaction_neuQs(irn14_to_c12), & 
               reaction_Qs(iro16_to_n14), reaction_neuQs(iro16_to_n14), & 
               
               reaction_Qs(irprot_to_neut), reaction_neuQs(irprot_to_neut), & 
               reaction_Qs(irneut_to_prot), reaction_neuQs(irneut_to_prot), & 
               reaction_Qs(irni56ec_to_co56), reaction_neuQs(irni56ec_to_co56), & 
               reaction_Qs(irco56ec_to_fe56), reaction_neuQs(irco56ec_to_fe56), & 
               
               reaction_Qs(irn14ag_lite), &
               reaction_Qs(ir_he4_he4_he4_to_c12), &
               reaction_Qs(ir1212), &
               reaction_Qs(ir1216_to_mg24), reaction_Qs(ir1216_to_si28), &
               reaction_Qs(ir1616a), reaction_Qs(ir1616g), &
               reaction_Qs(ir_ne20_ag_mg24), &
               reaction_Qs(ir_mg24_ag_si28), &
               reaction_Qs(ir_si28_ag_s32), &
               reaction_Qs(ir_s32_ag_ar36), &
               reaction_Qs(ir_ar36_ag_ca40), &
               reaction_Qs(ir_ca40_ag_ti44), &
               reaction_Qs(ir_ti44_ag_cr48), &
               reaction_Qs(ir_cr48_ag_fe52), &
               reaction_Qs(ir_fe52_ag_ni56), &               
               reaction_Qs(ir_fe52_ng_fe53), &       
               reaction_Qs(ir_fe53_ng_fe54), &       
               reaction_Qs(ir_fe54_ng_fe55), &       
               reaction_Qs(ir_fe55_ng_fe56), &                              
               reaction_Qs(irfe52neut_to_fe54), &               
               reaction_Qs(irfe52aprot_to_fe54), &               
               reaction_Qs(irfe54ng_to_fe56), &               
               reaction_Qs(irfe54aprot_to_fe56), &               
               reaction_Qs(irfe52aprot_to_ni56), &               
               reaction_Qs(irfe54prot_to_ni56), &               
!======Changed=====
               reaction_Qs(ir_prot_wk_neut), &
               reaction_Qs(ir_neut_wk_prot), &
!==================
               eps_total, eps_neu_total, do_eps_nuc_categories, &
               n% eps_nuc_categories(ipp), &
               n% eps_nuc_categories(icno), &
               n% eps_nuc_categories(i_burn_n), &
               n% eps_nuc_categories(i3alf), &
               n% eps_nuc_categories(icc), &
               n% eps_nuc_categories(ico), &
               n% eps_nuc_categories(ioo), &
               n% eps_nuc_categories(i_burn_ne), &
               n% eps_nuc_categories(i_burn_mg), &
               n% eps_nuc_categories(i_burn_si), &
               n% eps_nuc_categories(i_burn_s), &
               n% eps_nuc_categories(i_burn_ar), &
               n% eps_nuc_categories(i_burn_ca), &
               n% eps_nuc_categories(i_burn_ti), &
               n% eps_nuc_categories(i_burn_cr), &
               n% eps_nuc_categories(i_burn_fe), &
               n% eps_nuc_categories(iphoto), &
               n% eps_nuc_categories(i_ni56_co56), &
               n% eps_nuc_categories(i_co56_fe56), &
               .false., ierr)
         end subroutine get_approx21_plus_fe53_fe55_co56_eps_info
         
         subroutine store_partials
            integer :: i, j
            do i=1,num_isos
               ci = chem_id(i)
               Z_plus_N = dble(chem_isos% Z_plus_N(ci))
               d_eps_nuc_dx(i) = Qconv*n% d_eps_nuc_dy(i)/Z_plus_N
               dxdt(i) = Z_plus_N*dydt(i_rate, i)
               d_dxdt_dRho(i) = Z_plus_N*dydt(i_rate_dRho, i)
               d_dxdt_dT(i) = Z_plus_N*dydt(i_rate_dT, i)
               do j=1, num_isos
                  d_dxdt_dx(i,j) = &
                     n% d_dydt_dy(i,j)*Z_plus_N/chem_isos% Z_plus_N(chem_id(j))
               end do
            end do
         end subroutine store_partials
         
         subroutine get_rates_with_screening(ierr)
            use rates_def, only: reaction_inputs
            use rates_lib, only: eval_using_rate_tables
            use net_approx21, only: approx21_nrat
            use net_approx21_plus_co56, only: approx21_plus_co56_nrat
            
            integer, intent(out) :: ierr
            
            integer :: i, num
            real(dp) :: f
            logical :: okay
            
            include 'formats.dek'

            do i=1,num_reactions
               if (g% reaction_id(i) <= 0) then
                  write(*,2) 'g% reaction_id(i)', i, g% reaction_id(i)
                  stop 'get_rates_with_screening'
               end if
            end do
            
            if (.not. reuse_rate_raw) then ! get the raw reaction rates
               if (dbg) write(*,*) 'call eval_using_rate_tables'
               call eval_using_rate_tables( &
                  num_reactions, g% reaction_id, g% rate_table, g% rattab_f1, nrattab,  &
                  ye, logtemp, btemp, bden, rate_factors, g% logttab, &
                  rate_raw, rate_raw_dT, rate_raw_dRho, ierr) 
               if (ierr /= 0) then
                  if (dbg) write(*,*) 'ierr from eval_using_rate_tables'
                  return
               end if
               
               if (doing_timing) then
                  call system_clock(time1)
                  g% clock_net_rate_tables = g% clock_net_rate_tables + (time1 - time0)
                  time0 = time1
               end if
            end if

            if (g% doing_approx21) then
               if (g% add_fe53_fe55_co56_to_approx21) then
                  call approx21_plus_fe53_fe55_co56_rates(ierr)
               else if (g% add_co56_to_approx21) then
                  call approx21_plus_co56_rates(ierr)
               else
                  call approx21_rates(ierr)
               end if
               if (ierr /= 0) return            
            end if
            
            if (.not. reuse_rate_screened) then ! get the screened reaction rates
               ! get the reaction rates including screening factors
               if (dbg) write(*,*) 'call screen_net with init=.false.'
               call screen_net( &
                  g, num_isos, n% y, btemp, bden, logtemp, logrho, .false.,  &
                  rate_raw, rate_raw_dT, rate_raw_dRho, &
                  rate_screened, rate_screened_dT, rate_screened_dRho, &
                  n% screening_mode, n% theta_e_for_graboske_et_al, n% graboske_cache, &
                  screen_h1, screen_he4, zbar, abar, z2bar, ye, ierr)
               if (dbg) write(*,*) 'done screen_net with init=.false.'
               if (ierr /= 0) return
               if (g% doing_approx21) then
                  if (g% add_fe53_fe55_co56_to_approx21) then
                     num = approx21_plus_fe53_fe55_co56_nrat
                  else if (g% add_co56_to_approx21) then
                     num = approx21_plus_co56_nrat
                  else
                     num = approx21_nrat
                  end if
                  do i=num_reactions+1,num
                     rate_screened(i) = rate_raw(i)
                     rate_screened_dT(i) = rate_raw_dT(i)
                     rate_screened_dRho(i) = rate_raw_dRho(i)
                  end do
                  do i=1,num
                     dratdumdy1(i) = 0d0
                     dratdumdy2(i) = 0d0
                  end do           
               end if
            end if
            
            if (doing_timing) then
               call system_clock(time1)
               g% clock_net_screen = g% clock_net_screen + (time1 - time0)
               time0 = time1
            end if
            
         end subroutine get_rates_with_screening 

         subroutine approx21_rates(ierr)
            use net_approx21, only: &
               approx21_pa_pg_fractions, approx21_weak_rates
            integer, intent(out) :: ierr
            ierr = 0
            call approx21_pa_pg_fractions( &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
            if (ierr /= 0) return            
            call approx21_weak_rates( &
               n% y, rate_raw, rate_raw_dT, rate_raw_dRho, &
               btemp, bden, ye, eta, zbar, &
               weak_rate_factor, reuse_rate_screened, ierr)
            if (ierr /= 0) return            
         end subroutine approx21_rates

         subroutine approx21_plus_co56_rates(ierr)
            use net_approx21_plus_co56, only: &
               approx21_pa_pg_fractions, approx21_weak_rates
            integer, intent(out) :: ierr
            ierr = 0
            call approx21_pa_pg_fractions( &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
            if (ierr /= 0) return            
            call approx21_weak_rates( &
               n% y, rate_raw, rate_raw_dT, rate_raw_dRho, &
               btemp, bden, ye, eta, zbar, &
               weak_rate_factor, reuse_rate_screened, ierr)
            if (ierr /= 0) return            
         end subroutine approx21_plus_co56_rates

         subroutine approx21_plus_fe53_fe55_co56_rates(ierr)
            use net_approx21_plus_fe53_fe55_co56, only: &
               approx21_pa_pg_fractions, approx21_weak_rates
            integer, intent(out) :: ierr
            ierr = 0
            call approx21_pa_pg_fractions( &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
            if (ierr /= 0) return            
            call approx21_weak_rates( &
               n% y, rate_raw, rate_raw_dT, rate_raw_dRho, &
               btemp, bden, ye, eta, zbar, &
               weak_rate_factor, reuse_rate_screened, ierr)
            if (ierr /= 0) return            
         end subroutine approx21_plus_fe53_fe55_co56_rates

         subroutine get_weaklib_rates(ierr)
            use rates_def, only : Coulomb_Info
            use rates_lib, only: eval_weak_reaction_info, coulomb_set_context
            use net_def, only: other_kind

            type (Coulomb_Info), target :: cc_info
            type (Coulomb_Info), pointer :: cc
            real(dp) :: iso_z(num_isos)

            integer, intent(out) :: ierr
            integer :: i, j, id, ir
            include 'formats.dek'

            ! before getting the weaklib rates, the Coulomb_Info
            ! structure must be populated.  the ecapture routines need
            ! to know some local quantities (functions of the density,
            ! temperature, and composition), to calculate Coulomb
            ! corrections to the rates

            ierr = 0
            cc => cc_info

            do i=1,num_isos
               iso_z(i) = chem_isos% Z(n% g% chem_id(i))
            end do

            call coulomb_set_context( &
               cc, temp, rho, logtemp, logrho, zbar, abar, z2bar,  &
               n% theta_e_for_graboske_et_al, num_isos, n% y, iso_z)
            
            call eval_weak_reaction_info( &
               num_wk_reactions, &
               g% weaklib_ids(1:num_wk_reactions), &
               g% reaction_id_for_weak_reactions(1:num_wk_reactions), &
               cc, temp*1d-9, ye*rho, zbar, eta, d_eta_dlnT, d_eta_dlnRho, &
               n% lambda, n% dlambda_dlnT, n% dlambda_dlnRho, &
               n% Q, n% dQ_dlnT, n% dQ_dlnRho, &
               n% Qneu, n% dQneu_dlnT, n% dQneu_dlnRho, &
               ierr)
            if (weak_rate_factor < 1d0) then
               do i=1,num_wk_reactions
                  n% lambda(i) = weak_rate_factor*n% lambda(i)
                  n% dlambda_dlnT(i) = weak_rate_factor*n% dlambda_dlnT(i)
                  n% dlambda_dlnRho(i) = weak_rate_factor*n% dlambda_dlnRho(i)
               end do
            end if          
            if (doing_timing) then
               call system_clock(time1)
               g% clock_net_weak_rates = g% clock_net_weak_rates + (time1 - time0)
               time0 = time1
            end if
         end subroutine get_weaklib_rates
      
      end subroutine eval_net
         
         
      subroutine get_T_limit_factor( &
            temp, lnT, T_lo, T_hi, lnT_lo, lnT_hi, &
            min_ln_factor, min_factor, reuse_rate_screened, &
            factor, d_factor_dT)
         real(dp), intent(in) ::  &
            temp, lnT, T_lo, T_hi, lnT_lo, lnT_hi, &
            min_ln_factor, min_factor
         logical, intent(in) :: reuse_rate_screened
         real(dp), intent(out) :: &
            factor, d_factor_dT
         real(dp) :: ln_factor, d_ln_factor_dlnT
         factor = 1d0
         d_factor_dT = 0d0
         if (reuse_rate_screened) return
         if (temp <= T_lo) return
         if (temp >= T_hi) then
            factor = min_factor
            return
         end if
         ln_factor = min_ln_factor*(lnT - lnT_lo)/(lnT_hi - lnT_lo)
         d_ln_factor_dlnT = min_ln_factor/(lnT_hi - lnT_lo)
         factor = exp_cr(ln_factor)
         d_factor_dT = d_ln_factor_dlnT*factor/temp
      end subroutine get_T_limit_factor

         
      subroutine set_molar_abundances(g, num_isos, x, y, dbg, ierr)
         type (Net_General_Info), pointer :: g
         integer, intent(in) :: num_isos
         real(dp), intent(in) :: x(:)
         real(dp), intent(inout) :: y(:)
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         
         real(dp) :: sum
         integer :: i, ci
         character (len=256) :: message
         include 'formats.dek'
         sum = 0
         do i = 1, g% num_isos
            sum = sum + x(i)
            ci = g% chem_id(i)
            if (ci <= 0) then
               write(*,*) 'problem in set_molar_abundances'
               write(*,*) 'i', i
               write(*,*) 'g% num_isos', g% num_isos
               write(*,*) 'g% chem_id(i)', g% chem_id(i)
               stop 'set_molar_abundances' 
            end if
            y(i) = min(1d0, max(x(i), 0d0)) / chem_isos% Z_plus_N(ci)
         enddo
         
         
         
         return      ! let it go even with bad xsum.
         
         
         
   
         if (abs(sum - 1d0) > 1d-2) then
            ierr = -1
            if (dbg) then
               do i = 1, g% num_isos
                  ci = g% chem_id(i)
                  write(*,2) chem_isos% name(ci), i, x(i)
               end do
               write(*,1) 'abs(sum - 1d0)', abs(sum - 1d0)
            end if
            return
         end if
      
      end subroutine set_molar_abundances

      
      subroutine get_T_rho_args(temp, logtemp, rho, logrho, info)
         real(dp), intent(inout) :: temp, logtemp ! log10 of temp
         real(dp), intent(inout) :: rho, logrho ! log10 of rho
         integer, intent(out) :: info
         info = 0
         if (temp == arg_not_provided .and. logtemp == arg_not_provided) then
            info = -2
            return
         end if
         if (logtemp == arg_not_provided) logtemp = log10_cr(temp)
         if (temp == arg_not_provided) temp = exp10_cr(logtemp)
         if (temp <= 0) then
            info = -1
            return
         end if
         if (rho == arg_not_provided .and. logrho == arg_not_provided) then
            info = -3
            return
         end if
         if (logrho == arg_not_provided) logrho = log10_cr(rho)
         if (rho == arg_not_provided) rho = exp10_cr(logrho)
         if (rho <= 0) then
            info = -1
            return
         end if
      end subroutine get_T_rho_args
      
      
      subroutine do_clean_up_fractions(nzlo, nzhi, species, nz, xa, max_sum_abs, xsum_tol, ierr)
         integer, intent(in) :: nzlo, nzhi, species, nz
         real(dp), intent(inout) :: xa(:,:) ! (species, nz)
         real(dp), intent(in) :: max_sum_abs, xsum_tol
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
         if (nzlo == nzhi) then
            call do_clean1(species, xa(1: species, nzlo), nzlo, max_sum_abs, xsum_tol, ierr)
            return
         end if         
!x$OMP  PARALLEL DO PRIVATE(k, op_err)
         do k = nzlo, nzhi
            op_err = 0
            call do_clean1(species, xa(1: species, k), k, max_sum_abs, xsum_tol, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!x$OMP  END PARALLEL DO
      end subroutine do_clean_up_fractions
      

      subroutine do_clean1(species, xa, k, max_sum_abs, xsum_tol, ierr)
         use utils_lib
         integer, intent(in) :: species, k
         real(dp), intent(inout) :: xa(:) ! (species)
         real(dp), intent(in) :: max_sum_abs, xsum_tol
         integer, intent(out) :: ierr
         integer :: j
         real(dp) :: xsum
         real(dp), parameter :: tiny_x = 1d-99
         character (len=256) :: message
         if (max_sum_abs > 1) then ! check for crazy values
            xsum = sum(abs(xa(1: species)))
            if (is_bad(xsum) .or. xsum > max_sum_abs) then
               ierr = -1
               return
            end if
         end if
         ierr = 0
         do j = 1, species
            if (xa(j) < tiny_x) xa(j) = tiny_x
            if (xa(j) > 1) xa(j) = 1
         end do         
         xsum = sum(xa(1: species))         
         if (abs(xsum-1) > xsum_tol) then
            ierr = -1
            return
         end if
         xa(1: species) = xa(1: species)/xsum
      end subroutine do_clean1
      

      end module net_eval

