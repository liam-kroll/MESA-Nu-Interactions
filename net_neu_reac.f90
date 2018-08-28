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
!===========Changed==============      
      module net_neu_reac

      use rates_def
      use net_def, only: Net_Info, Net_General_Info

      contains

      subroutine net_neu_place(neu_reac_info, g, n, ierr)
         real(dp), intent(in) :: neu_reac_info(4)
         integer, intent(out) :: ierr
         type (Net_Info), pointer :: n
         type (Net_General_Info), pointer :: g

         integer :: w, reac_in

         ierr = 0

!PLACEMENT OF CHARGED CURRENT NU REACTION RATES INTO NET INFO 

         do w = 1, g% num_reactions
            reac_in = g% reaction_id(w)
            if (reac_in == ir_prot_wk_neut) then
               n% rate_screened(w:w) = neu_reac_info(1)
               n% rate_screened_dT(w:w) = 0.0d0
               n% rate_screened_dRho(w:w) = 0.0d0
               n% reaction_Qs(reac_in) = neu_reac_info(2)
               n% reaction_neuQs(reac_in) = 0.0d0
            end if
            if (reac_in == ir_neut_wk_prot) then
               n% rate_screened(w:w) = neu_reac_info(3)
               n% rate_screened_dT(w:w) = 0.0d0
               n% rate_screened_dRho(w:w) = 0.0d0
               n% reaction_Qs(reac_in) = neu_reac_info(4)
               n% reaction_neuQs(reac_in) = 0.0d0
            end if
         end do
!END OF PLACEMENT OF CHARGED CURRENT NU REACTION RATES INTO NET INFO 

      end subroutine net_neu_place

      end module net_neu_reac
!========================================
