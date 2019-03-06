module constants

implicit none

complex*16       , public, parameter   :: cmplx_i = (0.0,1.0), cmplx_0 = (0.0,0.0)
real*8          , public, parameter    :: twopi = 2*acos(-1.0), kboltz = 8.6173303505050d-5
integer          , public, parameter   :: nspin = 2
character(len=50), public, parameter   :: hamname     = 'hf_ham.in' ,   &
                                          socname     = 'hf_soc.in' ,   &
                                          xcname      = 'hf_xc.in',     &
                                          coulombname = 'hf_coulomb.in',&
                                          elecname    = 'hf_mag.in',    &
                                          kptname     = 'hf_kpt.in',    &
                                          lattname    = 'hf_struct.in', &
                                          pathname    = 'hf_path.in'   
                                              

end module constants
