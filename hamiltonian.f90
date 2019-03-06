module hamilton

  ! Case for spinorb = .true.
  complex*16, allocatable, save, public :: hf_potsoc(:,:,:,:,:), ham_fullsoc(:,:,:)
  ! Case for spinorb = .false.
  complex*16, allocatable, save, public :: hf_pot(:,:,:,:), ham_full(:,:,:,:)

  public :: init_hf_pot
  public :: set_hf_pot
  public :: set_full_ham

  contains

  subroutine init_hf_pot()

    use param,          only : u_mat, norb, nsite, ntot, nkpt, spinorb
    use dens,           only : init_dens, init_denssoc
    use constants,      only : nspin, cmplx_0

    integer    :: ispin, ispin1, isite, iorb, iorb1, iorb2, iorb3

    if (.not.spinorb) then

       allocate(hf_pot(norb,norb,nsite,nspin))
       allocate(ham_full(norb*nsite,norb*nsite,nspin,nkpt))

       hf_pot = cmplx_0
       do isite = 1, nsite
          do iorb = 1, norb      
             do iorb1 = 1, norb
                do iorb2 = 1, norb
                   do iorb3 = 1, norb
                      hf_pot(iorb,iorb1,isite,1) = hf_pot(iorb,iorb1,isite,1) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                     u_mat(iorb,iorb3,iorb2,iorb1,isite))*init_dens(iorb2,iorb3,isite,1) + &
                                                     u_mat(iorb,iorb1,iorb2,iorb3,isite)*init_dens(iorb2,iorb3,isite,2) 

                      hf_pot(iorb,iorb1,isite,2) = hf_pot(iorb,iorb1,isite,2) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                     u_mat(iorb,iorb3,iorb2,iorb1,isite))*init_dens(iorb2,iorb3,isite,2) + &
                                                     u_mat(iorb,iorb1,iorb2,iorb3,isite)*init_dens(iorb2,iorb3,isite,1)
                   end do
                end do
             end do
          end do
       end do

    else if (spinorb) then

       allocate(hf_potsoc(norb,norb,nsite,nspin,nspin))
       allocate(ham_fullsoc(ntot,ntot,nkpt))

       hf_potsoc = cmplx_0
       do isite = 1, nsite
          do iorb = 1, norb
             do iorb1 = 1, norb
                do iorb2 = 1, norb
                   do iorb3 = 1, norb
                      hf_potsoc(iorb,iorb1,isite,1,1) = hf_potsoc(iorb,iorb1,isite,1,1) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                        u_mat(iorb,iorb3,iorb2,iorb1,isite))*init_denssoc(iorb2,iorb3,isite,1,1) + &
                                                        u_mat(iorb,iorb1,iorb2,iorb3,isite)*init_denssoc(iorb2,iorb3,isite,2,2)

                      hf_potsoc(iorb,iorb1,isite,2,2) = hf_potsoc(iorb,iorb1,isite,2,2) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                        u_mat(iorb,iorb3,iorb2,iorb1,isite))*init_denssoc(iorb2,iorb3,isite,2,2) + &
                                                        u_mat(iorb,iorb1,iorb2,iorb3,isite)*init_denssoc(iorb2,iorb3,isite,1,1)

                      hf_potsoc(iorb,iorb1,isite,1,2) = hf_potsoc(iorb,iorb1,isite,1,2) - (u_mat(iorb,iorb3,iorb2,iorb1,isite)*    &
                                                        init_denssoc(iorb2,iorb3,isite,2,1))

                      hf_potsoc(iorb,iorb1,isite,2,1) = hf_potsoc(iorb,iorb1,isite,2,1) - (u_mat(iorb,iorb3,iorb2,iorb1,isite)*    &
                                                        init_denssoc(iorb2,iorb3,isite,1,2))
                   end do
                end do
             end do
          end do
       end do

    end if

  end subroutine init_hf_pot

  subroutine set_hf_pot(rho)

    use param,          only : u_mat, norb, nsite 
    use constants,      only : nspin

    complex*16, intent(in)  :: rho(norb,norb,nsite,nspin)
    integer                 :: ispin, ispin1, isite, iorb, iorb1, iorb2, iorb3

    hf_pot = cmplx_0
    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             do iorb2 = 1, norb
                do iorb3 = 1, norb
                   hf_pot(iorb,iorb1,isite,1) = hf_pot(iorb,iorb1,isite,1) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                u_mat(iorb,iorb3,iorb2,iorb1,isite))*rho(iorb2,iorb3,isite,1) +       &
                                                u_mat(iorb,iorb1,iorb2,iorb3,isite)*rho(iorb2,iorb3,isite,2)
 
                   hf_pot(iorb,iorb1,isite,2) = hf_pot(iorb,iorb1,isite,2) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                u_mat(iorb,iorb3,iorb2,iorb1,isite))*rho(iorb2,iorb3,isite,2) +       &
                                                u_mat(iorb,iorb1,iorb2,iorb3,isite)*rho(iorb2,iorb3,isite,1)
                end do
             end do
          end do
       end do
    end do
  
  end subroutine set_hf_pot

  subroutine set_hf_pot_soc(rho)

    use param,          only : u_mat, norb, nsite
    use constants,      only : nspin

    complex*16, intent(in)  :: rho(norb,norb,nsite,nspin,nspin)
    integer                 :: ispin, ispin1, isite, iorb, iorb1, iorb2, iorb3

    hf_potsoc = cmplx_0
    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             do iorb2 = 1, norb
                do iorb3 = 1, norb
                   hf_potsoc(iorb,iorb1,isite,1,1) = hf_potsoc(iorb,iorb1,isite,1,1) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                     u_mat(iorb,iorb3,iorb2,iorb1,isite))*rho(iorb2,iorb3,isite,1,1) +       &
                                                     u_mat(iorb,iorb1,iorb2,iorb3,isite)*rho(iorb2,iorb3,isite,2,2)

                   hf_potsoc(iorb,iorb1,isite,2,2) = hf_potsoc(iorb,iorb1,isite,2,2) + (u_mat(iorb,iorb1,iorb2,iorb3,isite)   - &
                                                     u_mat(iorb,iorb3,iorb2,iorb1,isite))*rho(iorb2,iorb3,isite,2,2) +       &
                                                     u_mat(iorb,iorb1,iorb2,iorb3,isite)*rho(iorb2,iorb3,isite,1,1)

                   hf_potsoc(iorb,iorb1,isite,1,2) = hf_potsoc(iorb,iorb1,isite,1,2) - (u_mat(iorb,iorb3,iorb2,iorb1,isite)*    &
                                                     rho(iorb2,iorb3,isite,2,1))

                   hf_potsoc(iorb,iorb1,isite,2,1) = hf_potsoc(iorb,iorb1,isite,2,1) - (u_mat(iorb,iorb3,iorb2,iorb1,isite)*    &
                                                     rho(iorb2,iorb3,isite,1,2))
                end do
             end do
          end do
       end do
    end do

  end subroutine set_hf_pot_soc


  subroutine set_full_ham()

    use param,          only : ntot, nkpt, norb, nsite, ham_k, ham_ksoc, spinorb
    use constants,      only : nspin, cmplx_0

    integer    :: ikpt, ispin, ispin1, isite, isite1

    if (.not.spinorb) then

       !Sum kinetic part to the hamiltonian
       ham_full = cmplx_0
       do ikpt = 1, nkpt
          do isite = 1, nsite
             do isite1 = 1, nsite
                ham_full((isite-1)*norb+1:isite*norb,(isite1-1)*norb+1:isite1*norb,1,ikpt) = ham_k(:,:,isite,isite1,ikpt,1)
                ham_full((isite-1)*norb+1:isite*norb,(isite1-1)*norb+1:isite1*norb,2,ikpt) = ham_k(:,:,isite,isite1,ikpt,2)
             end do        
          end do
       end do   

       !Sum hf_potential
       do ikpt = 1, nkpt
          do isite = 1, nsite
             ham_full((isite-1)*norb+1:isite*norb,(isite-1)*norb+1:isite*norb,1,ikpt) =   &
             ham_full((isite-1)*norb+1:isite*norb,(isite-1)*norb+1:isite*norb,1,ikpt) + hf_pot(:,:,isite,1)
             ham_full((isite-1)*norb+1:isite*norb,(isite-1)*norb+1:isite*norb,2,ikpt) =   &
             ham_full((isite-1)*norb+1:isite*norb,(isite-1)*norb+1:isite*norb,2,ikpt) + hf_pot(:,:,isite,2)
          end do
       end do 

    else if (spinorb) then

       !Sum kinetic part to the hamiltonian
       ham_fullsoc = cmplx_0
       do ikpt = 1, nkpt
          do ispin = 1, nspin
             do ispin1 = 1, nspin
                do isite = 1, nsite
                   do isite1 = 1, nsite
                      ham_fullsoc((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                                  (isite1-1)*norb+1+(ispin1-1)*norb*nsite:isite1*norb+(ispin1-1)*norb*nsite,ikpt) = &
                                  ham_ksoc(:,:,isite,isite1,ikpt,ispin,ispin1)
                   end do
                end do
             end do
          end do
       end do

       !Sum hf_potential
       do ikpt = 1, nkpt
          do ispin = 1, nspin
             do ispin1 = 1, nspin
                do isite = 1, nsite
                   ham_fullsoc((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                               (isite-1)*norb+1+(ispin1-1)*norb*nsite:isite*norb+(ispin1-1)*norb*nsite,ikpt) =   &
                   ham_fullsoc((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                               (isite-1)*norb+1+(ispin1-1)*norb*nsite:isite*norb+(ispin1-1)*norb*nsite,ikpt) +   &
                   hf_potsoc(:,:,isite,ispin,ispin1)
                end do
             end do
          end do
       end do

    end if
         
  end subroutine set_full_ham

end module hamilton

