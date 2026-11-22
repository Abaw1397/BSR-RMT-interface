!======================================================================
!     UTILITY       BSR_SPLINEWAVES
!
!               C O P Y R I G H T -- 2026
!
!     Written by:   Aaron T. Bondy and Juan C. del Valle
!                   email: bondy11u@gmail.com, jcdvaller@gmail.com 
!
!     Original version: Juan C. del Valle, 2024.
!
!======================================================================
!   based on rsol.nnn files, it creates the Splinewaves file
!   This routine works in both LS and jK coupling
!======================================================================
!
!   INPUT FILES:
!
!      rsol.nnn      - BSR R-matrix eigenvectors in B-spline representation.
!                      Length of each vector is nhm = kch*ns + kcp.
!
!   OUTPUT FILES:
!
!      Splinewaves   - unformatted channel-only file in the RMT style.
!                      Each state contains kch*(ns-1) entries:
!                      in each channel, the coefficients 2:ns are kept.
!
!      bsr_Splinewaves.log - contains dimension information and separate
!                      norm diagnostics for channels and perturbers.
!
!      bsr_Splinewaves-states.log - contains diagnostics for each eigenvector,
!                      including energy and dominant perturber coefficient.
!                     
!   NOTES:
!
!      From BSR, the last kcp entries of each rsol eigenvector are perturber coefficients.
!      The log files in this routine use them, but they are not written to Splinewaves.
!      Splinewaves stores only the channel part of each eigenvector.
!      The perturbers affect energies, dipole matrix elements, etc. that are
!         passed to RMT.
!
!   USAGE:
!
!      bsr_Splinewaves nlsp=Number_of_Partial_Waves
!======================================================================
      Program bsr_Splinewaves

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AF, AF_all = 'Splinewaves'
      Integer, allocatable :: nchn(:), nstk(:), kcps(:), nhms(:), nsvals(:)
      Real(8), allocatable :: eval(:), vfull(:), vred(:)
      Real(8), allocatable :: pmin(:), pmax(:), psum(:)
      Real(8), allocatable :: fmin(:), fmax(:), fsum(:)
      Real(8), allocatable :: cmin(:), cmax(:), csum(:)
      Integer :: idiag, ip, ipmax
      Real(8) :: pert_coeff

      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'Convert BSR rsol.nnn files to RMT Splinewaves'
        write(*,'(/a)') 'Call as: bsr_Splinewaves nlsp=..'
        write(*,'(/a)')
        Stop
      end if

      Call Read_iarg('nlsp',nlsp)
      if(nlsp.le.0) Stop 'bsr_Splinewaves: nlsp <= 0'

      iout = 2
      open(iout,file=AF_all,form='UNFORMATTED')

      idiag = 4
      open(idiag,file='bsr_Splinewaves-states.log',status='REPLACE')
      write(idiag,'(a)') 'bsr_Splinewaves state-by-state diagnostics'
      write(idiag,'(a)') &
       'Norms and fractions below are not probabilities.'// &
       ' They are obtained from squared eigenvector coefficients without the overlap matrix.'
      write(idiag,'(a)') 'pert_coeff_fraction can be used to identify states with'// &
       ' perturber coefficients that are large compared to the channel coefficients.'
      write(idiag,'(a)') &
       ' ilsp   state     kcp               energy_au'// &
       '      channel_coeff_norm         pert_coeff_norm'// &
       '      pert_coeff_fraction   dominant_pert_index'// &
       '      dominant_pert_coeff'

      Allocate(nchn(nlsp), nstk(nlsp), kcps(nlsp), nhms(nlsp), nsvals(nlsp))
      Allocate(pmin(nlsp), pmax(nlsp), psum(nlsp))
      Allocate(fmin(nlsp), fmax(nlsp), fsum(nlsp))
      Allocate(cmin(nlsp), cmax(nlsp), csum(nlsp))

      pmin = 0.d0; pmax = 0.d0; psum = 0.d0
      fmin = 0.d0; fmax = 0.d0; fsum = 0.d0
      cmin = 0.d0; cmax = 0.d0; csum = 0.d0

      Do ilsp = 1,nlsp

       write(AF,'(a,i3.3)') 'rsol.',ilsp
       Call Check_file(AF)
       inp=1; open(inp,file=AF,form='UNFORMATTED')

       read(inp) nhm,khm,kch,kcp,ns
       Allocate(eval(khm))
       read(inp) eval

       if(ns.le.1) Stop 'bsr_Splinewaves: ns <= 1'
       if(ilsp.gt.1.and.ns.ne.nsvals(1)) &
        Stop 'bsr_Splinewaves: inconsistent ns across rsol files'
       if(nhm.ne.kch*ns+kcp) &
        Stop 'bsr_Splinewaves: inconsistent rsol dimensions'

       Allocate(vfull(nhm))
       Allocate(vred(kch*(ns-1)))

       if(kcp.gt.0) then
        pmin(ilsp) = huge(0.d0)
        fmin(ilsp) = huge(0.d0)
       end if
       cmin(ilsp) = huge(0.d0)

       Do j = 1,khm
        read(inp) vfull

        io = 0
!       There are ns coefficients in each channel in rsol and ns-1 in Splinewaves.
        Do ich = 1,kch
         i1 = (ich-1)*ns + 2
         i2 = ich*ns
         vred(io+1:io+ns-1) = vfull(i1:i2)
         io = io + ns - 1
        End do

        chan_norm = sum(vfull(1:kch*ns)**2)
        pert_norm = 0.d0
        if(kcp.gt.0) pert_norm = sum(vfull(kch*ns+1:nhm)**2)
        denom = chan_norm + pert_norm
        frac = 0.d0
        if(denom.gt.0.d0) frac = pert_norm / denom

        ipmax = 0
        pert_coeff = 0.d0
        Do ip = 1,kcp
         if(abs(vfull(kch*ns+ip)).gt.abs(pert_coeff)) then
          ipmax = ip
          pert_coeff = vfull(kch*ns+ip)
         end if
        End do

        write(idiag,'(3i8,4es24.14,i8,es24.14)') &
         ilsp,j,kcp,eval(j),chan_norm,pert_norm,frac,ipmax,pert_coeff

        cmin(ilsp) = min(cmin(ilsp),chan_norm)
        cmax(ilsp) = max(cmax(ilsp),chan_norm)
        csum(ilsp) = csum(ilsp) + chan_norm
        if(kcp.gt.0) then
         pmin(ilsp) = min(pmin(ilsp),pert_norm)
         pmax(ilsp) = max(pmax(ilsp),pert_norm)
         psum(ilsp) = psum(ilsp) + pert_norm
         fmin(ilsp) = min(fmin(ilsp),frac)
         fmax(ilsp) = max(fmax(ilsp),frac)
         fsum(ilsp) = fsum(ilsp) + frac
        end if

        write(iout) (vred(i),i=1,kch*(ns-1))
       End do

       close(inp)
       Deallocate(eval)
       Deallocate(vfull)
       Deallocate(vred)

       nchn(ilsp) = kch
       nstk(ilsp) = khm
       kcps(ilsp) = kcp
       nhms(ilsp) = nhm
       nsvals(ilsp) = ns

       if(kcp.eq.0) then
        pmin(ilsp) = 0.d0
        fmin(ilsp) = 0.d0
       end if

      End do

      close(iout)
      close(idiag)

      ipri = 3; open(ipri,file='bsr_Splinewaves.log')

      write(ipri,'(a)') 'bsr_Splinewaves'
      write(ipri,'(a)') 'RMT receives only channel information: kch*(ns-1) entries per state/eigenvector.'
      write(ipri,'(a)') 'Full rsol vectors need to be read; but the perturber tail (the last kcp entries) is omitted deliberately.'
      write(ipri,'(a)')

      write(ipri,'(a)') 'nlsp (inast), number of partial waves:'
      write(ipri,'(i5)') nlsp

      write(ipri,'(a)') 'raw rsol ns:'
      write(ipri,'(i5)') nsvals(1)

      write(ipri,'(a)') 'RMT Splinewaves contains ndim entries per channel:'
      write(ipri,'(i5)') nsvals(1)-1

      nchmx = maxval(nchn)
      write(ipri,'(a)') 'nchmx, maximum number of channels:'
      write(ipri,'(i5)') nchmx

      nstmx = maxval(nstk)
      write(ipri,'(a)') 'nstmx, maximum number of RM solutions:'
      write(ipri,'(i5)') nstmx

      write(ipri,'(a)') 'nchn(1:nlsp) - number of channels for each partial wave:'
      write(ipri,'(10i10)') nchn

      write(ipri,'(a)') 'nstk(1:nlsp) - number of RM solutions for each partial wave:'
      write(ipri,'(10i10)') nstk

      write(ipri,'(a)') 'nhm(1:nlsp) - BSR rsol vector length:'
      write(ipri,'(10i10)') nhms

      write(ipri,'(a)') 'kcp(1:nlsp) - number of perturber coefficients in rsol:'
      write(ipri,'(10i10)') kcps

      write(ipri,'(/a)') 'per-partial-wave norm diagnostics:'
      write(ipri,'(a)') &
       ' ilsp   kch   kcp   khm      min_chan      max_chan      avg_chan'// &
       '      min_pert      max_pert      avg_pert      min_pfrac'// &
       '     max_pfrac     avg_pfrac'
      Do ilsp = 1,nlsp
       avg_chan = csum(ilsp) / max(1,nstk(ilsp))
       avg_pert = psum(ilsp) / max(1,nstk(ilsp))
       avg_frac = fsum(ilsp) / max(1,nstk(ilsp))
       write(ipri,'(4i6,9es14.5)') ilsp,nchn(ilsp),kcps(ilsp),nstk(ilsp), &
        cmin(ilsp),cmax(ilsp),avg_chan, &
        pmin(ilsp),pmax(ilsp),avg_pert, &
        fmin(ilsp),fmax(ilsp),avg_frac
      End do

      close(ipri)

      End Program bsr_Splinewaves
