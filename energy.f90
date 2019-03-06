module energies

  real*8, save, public     :: etot

  contains

  subroutine evaluate_energy(eigenval,pot,density)

    use param,         only : ntot, norb, nsite, nkpt
    use constants,     only : nspin
    use dens,          only : efermi

    integer                :: itot, ikpt, ispin, ispin1, iorb, iorb1, isite
    real*8, intent(in)     :: eigenval(norb*nsite,nspin,nkpt)
    complex*16, intent(in) :: pot(norb,norb,nsite,nspin), density(norb,norb,nsite,nspin)

    etot = 0.0
    do itot = 1, norb*nsite
       do ispin = 1, nspin
          do ikpt = 1, nkpt
             if (eigenval(itot,ispin,ikpt) > efermi) then
                cycle
             end if

             etot = etot + eigenval(itot,ispin,ikpt)/real(nkpt)
          end do
       end do
    end do

    do ispin = 1, nspin
       do isite = 1, nsite
          do iorb = 1, norb
             do iorb1 = 1, norb
                etot = etot - 0.5*pot(iorb,iorb1,isite,ispin)*density(iorb,iorb1,isite,ispin)
             end do
          end do
       end do
    end do

  end subroutine evaluate_energy

  subroutine evaluate_energy_soc(eigenval,pot,density)

    use param,         only : ntot, norb, nsite, nkpt
    use constants,     only : nspin
    use dens,          only : efermi

    integer                :: itot, ikpt, ispin, ispin1, iorb, iorb1, isite
    real*8, intent(in)     :: eigenval(ntot,nkpt)
    complex*16, intent(in) :: pot(norb,norb,nsite,nspin,nspin), density(norb,norb,nsite,nspin,nspin)
    
    etot = 0.0
    do itot = 1, ntot
       do ikpt = 1, nkpt
          if (eigenval(itot,ikpt) > efermi) then
             cycle
          end if
          
          etot = etot + eigenval(itot,ikpt)/real(nkpt)

       end do
    end do

    do ispin = 1, nspin
       do ispin1 = 1, nspin
          do isite = 1, nsite
             do iorb = 1, norb
                do iorb1 = 1, norb
                   etot = etot - 0.5*pot(iorb,iorb1,isite,ispin,ispin1)*density(iorb,iorb1,isite,ispin,ispin1)
                end do
             end do
          end do
       end do
    end do

  end subroutine evaluate_energy_soc

end module energies
