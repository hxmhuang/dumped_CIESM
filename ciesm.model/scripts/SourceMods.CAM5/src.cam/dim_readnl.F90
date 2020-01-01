      module dim_readnl

      use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl 
      use abortutils,       only: endrun

       implicit none
       private
       save

       public :: dim_intr_readnl
       REAL(r8),public, DIMENSION(14) :: nbandsw
       REAL(r8),public, DIMENSION(16) :: nbandlw

       CONTAINS

      subroutine dim_intr_readnl
      IMPLICIT REAL (A-H,O-Z),INTEGER (I-N)

        DO N=1,14
           nbandsw(N) = N
        ENDDO

        DO N=1,16
           nbandlw(N) = N
        ENDDO

      end subroutine dim_intr_readnl
      end module dim_readnl    
