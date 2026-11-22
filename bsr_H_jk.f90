!======================================================================
!     UTILITY       BSR_H_JK
!
!               C O P Y R I G H T -- 2026
!
!     Written by:   Aaron T. Bondy and Juan C. del Valle
!                   email: bondy11u@gmail.com, jcdvaller@gmail.com
!
!     Original LS version: Juan C. del Valle, 2024.
!
!======================================================================
!   This routine merges the h.nnn files generated in BSR into the H file
!   used by RMT for jK coupling, in the RMI format.
!======================================================================
!
!   INPUT FILES:
!
!      Each h.nnn file contains Hamiltonian data for one or more symmetry
!      blocks. With the setting folder=1, this program reads 001/h.001,
!      002/h.002, etc.
!
!   OUTPUT FILES:
!
!      H is a sequential unformatted file that has five header records,
!      after which there are six records for each symmetry block.
!
!      bsr_H_jk_v3_candidate_02.log
!      The log file lists the channel labels and pole energy order in
!      each of these blocks.
!
!   NOTES:
!
!      Within each channel, the orbital angular momentum of the outer
!      electron is stored in L2P and the intermediate angular momentum K
!      is stored as 2K+1 in jkch. Here, we write L2P to the H file.
!      L2P and jkch are documented in the log.
!
!   USAGE:
!
!      bsr_H_jk ih1=1 ih2=2
!======================================================================

      Implicit none

      Integer :: ih1=1,ih2=999,folder=0,reverse_poles=-1,rmt_nfdm=18
      Integer :: in=1,iout=2,ipri=3,ih,last_file,nfiles
      Integer :: nelc,nz,lrang2_current,lamax,nast
      Integer :: nelc_ref,nz_ref,lamax_ref,nast_ref,lrang2
      Integer :: lrgl,nspn,npty,nchan,mnp2,more0,more,nstat
      Integer :: i,j,k
      Real(8) :: ra,bsto,ra_ref,bsto_ref,rmt_deltar=0.d0,max_rmt_deltar
      Real(8), allocatable :: enat(:),enat_ref(:),coeff(:,:)
      Real(8), allocatable :: value(:),wmat(:,:),cf(:,:,:)
      Integer, allocatable :: lat(:),lat_ref(:),isat(:),isat_ref(:)
      Integer, allocatable :: ipat(:),ipat_ref(:),nconat(:)
      Integer, allocatable :: l2p(:),jkch(:)
      Character(80) :: af,arg
      Logical :: exists,first_file,ascending,descending,reverse_this

      Call get_command_argument(1,arg)
      if(trim(arg).eq.'?') then
       write(*,'(/a)') 'bsr_H_jk merges the BSR h.nnn files into one H file for RMT in jK coupling.'
       write(*,'(/a)') 'Example: bsr_H_jk ih1=1 ih2=2'
       write(*,'(/a)') 'ih1 and ih2 are the first and last file numbers to check. By default they are 1 and 999.'
       write(*,'(a)')  'Use folder=1 if files are stored as 001/h.001, 002/h.002, etc.'
       write(*,'(/a)') 'The default setting reverse_poles=-1 inspects the pole energies within each symmetry block'
       write(*,'(a)')  'and writes them in increasing order. Setting reverse_poles=0 keeps the input order,'
       write(*,'(a)')  'while reverse_poles=1 reverses it. Each pole energy stays paired with its surface'
       write(*,'(a)')  'amplitudes, as required by RMT.'
       write(*,'(/a)') 'Optional: provide the intended RMT outer region grid rmt_deltar to ensure'
       write(*,'(a)')  'rmt_nfdm × rmt_deltar < R-matrix radius. The default is rmt_nfdm=18 unless'
       write(*,'(a)')  'a different value is supplied.'
       Stop
      end if

      Call Read_iarg('ih1',ih1)
      Call Read_iarg('ih2',ih2)
      Call Read_iarg('folder',folder)
      Call Read_iarg('reverse_poles',reverse_poles)
      Call Read_rarg('rmt_deltar',rmt_deltar)
      Call Read_iarg('rmt_nfdm',rmt_nfdm)
      if(ih1.le.0) Stop 'bsr_H: ih1 <= 0'
      if(ih2.lt.ih1) Stop 'bsr_H: ih2 < ih1'
      if(reverse_poles.lt.-1 .or. reverse_poles.gt.1) &
       Stop 'bsr_H_jk: reverse_poles must be -1, 0 or 1'
      if(rmt_nfdm.le.0) Stop 'bsr_H_jk: rmt_nfdm must be positive'
      if(rmt_deltar.lt.0.d0) Stop 'bsr_H_jk: rmt_deltar must not be negative'

      Open(ipri,file='bsr_H_jk_v3_candidate_02.log',status='replace')
      write(ipri,'(a)') 'bsr_H_jk'
      write(ipri,*)

!     On the first pass, the last h.nnn file is located and the largest
!     absolute value of LRANG2 across the input files is identified.
!     This value is needed in the common H header.

      lrang2=0; last_file=0; nfiles=0
      Do ih=ih1,ih2
       write(af,'(a,i3.3)') 'h.',ih
       if(folder.ne.0) write(af,'(i3.3,a,i3.3)') ih,'/h.',ih
       inquire(file=trim(af),exist=exists)
       if(.not.exists) Cycle

       Open(in,file=trim(af),form='unformatted',status='old')
       read(in) nelc,nz,lrang2_current,lamax,nast,ra,bsto
       Close(in)
       lrang2=max(lrang2,abs(lrang2_current))
       last_file=ih
       nfiles=nfiles+1
      End do
      if(nfiles.eq.0) Stop 'bsr_H: no input h.nnn files found'

      write(ipri,'(a,i0)') 'number of h files: ',nfiles
      write(ipri,'(a,i0)') 'LRANG2 written to H: ',lrang2
      write(ipri,'(a,i0)') 'reverse_poles setting: ',reverse_poles
      write(ipri,*)

      Open(iout,file='H',form='unformatted',status='replace')
      first_file=.true.

!     The second pass writes the common header only once, followed by
!     the symmetry blocks from every h.nnn file in sequence.

      Do ih=ih1,last_file
       write(af,'(a,i3.3)') 'h.',ih
       if(folder.ne.0) write(af,'(i3.3,a,i3.3)') ih,'/h.',ih
       inquire(file=trim(af),exist=exists)
       if(.not.exists) Cycle

       Open(in,file=trim(af),form='unformatted',status='old')
       read(in) nelc,nz,lrang2_current,lamax,nast,ra,bsto

       Allocate(enat(nast),lat(nast),isat(nast),ipat(nast))
       read(in) enat
       read(in) lat
!      Within BSR, ISAT (0 for the jK targets here) and the target parity
!      IPAT are stored in the same record. Both are read to ensure
!      consistent target data in the provided h.nnn files. Only ISAT is
!      written to this record in H.
       read(in) isat,ipat
       Allocate(coeff(3,abs(lrang2_current)))
       read(in) coeff

       if(first_file) then
        nelc_ref=nelc; nz_ref=nz; lamax_ref=lamax; nast_ref=nast
        ra_ref=ra; bsto_ref=bsto
        Allocate(enat_ref(nast),lat_ref(nast),isat_ref(nast),ipat_ref(nast))
        enat_ref=enat; lat_ref=lat; isat_ref=isat; ipat_ref=ipat

!       The innermost point of the RMT finite-difference stencil is at
!       r = RA − nfdm × DeltaR, where RA is the R-matrix radius.
!       r > 0 is required.
        max_rmt_deltar=ra_ref/real(rmt_nfdm,8)
        write(*,'(a,f12.6)') 'bsr_H_jk: R-matrix radius = ',ra_ref
        write(*,'(a,i0,a,f12.6)') 'bsr_H_jk: for rmt_nfdm=',rmt_nfdm, &
             ', DeltaR must be smaller than ',max_rmt_deltar
        write(ipri,'(a,a)') 'header from ',trim(af)
        write(ipri,'(a,5(1x,i0),2(1x,f12.6))') 'global:', &
             nelc_ref,nz_ref,lrang2,lamax_ref,nast_ref,ra_ref,bsto_ref
        write(ipri,'(a,1000(1x,es15.7))') 'target energies:', &
             (enat_ref(i),i=1,nast_ref)
        write(ipri,'(a,1000(1x,i0))') 'target angular labels:', &
             (lat_ref(i),i=1,nast_ref)
        write(ipri,'(a,1000(1x,i0))') 'target spin labels:', &
             (isat_ref(i),i=1,nast_ref)
        write(ipri,'(a,1000(1x,i0))') 'target parity labels read:', &
             (ipat_ref(i),i=1,nast_ref)
        write(ipri,'(a,f12.6)') 'R-matrix radius:',ra_ref
        write(ipri,'(a,i0,a,f12.6)') 'RMT grid guard for nfdm=',rmt_nfdm, &
             ', DeltaR < ',max_rmt_deltar
        write(ipri,*)
        if(rmt_deltar.gt.0.d0 .and. &
           real(rmt_nfdm,8)*rmt_deltar.ge.ra_ref) then
         write(*,'(a,f12.6)') 'requested RMT DeltaR = ',rmt_deltar
         Stop 'bsr_H_jk: RMT inner stencil lies outside the B-spline box'
        end if

!       Write the five header records read by RMT's parse_H_file.
        write(iout) nelc_ref,nz_ref,lrang2,lamax_ref,nast_ref,ra_ref,bsto_ref
        write(iout) enat_ref
        write(iout) lat_ref
        write(iout) isat_ref
!       The last of the five header records contains three entries
!       representing the Buttle correction for each orbital angular
!       momentum l. These entries are all set to 0 since BSR does not use
!       Buttle corrections here, though RMT needs to read these entries.
        write(iout) ((0.d0,k=1,3),j=1,lrang2)
        first_file=.false.
       else
        if(nelc.ne.nelc_ref) Stop 'bsr_H: inconsistent NELC'
        if(nz.ne.nz_ref) Stop 'bsr_H: inconsistent NZ'
        if(nast.ne.nast_ref) Stop 'bsr_H: inconsistent NAST'
        if(lamax.ne.lamax_ref) Stop 'bsr_H: inconsistent LAMAX'
        if(abs(ra-ra_ref).gt.100*epsilon(ra)*max(1.d0,abs(ra_ref))) &
         Stop 'bsr_H: inconsistent RA'
        if(abs(bsto-bsto_ref).gt.100*epsilon(bsto)*max(1.d0,abs(bsto_ref))) &
         Stop 'bsr_H: inconsistent BSTO'
        if(any(abs(enat-enat_ref).gt.100*epsilon(enat)* &
               max(1.d0,abs(enat_ref)))) &
         Stop 'bsr_H: inconsistent target energies'
        if(any(lat.ne.lat_ref)) Stop 'bsr_H: inconsistent target angular momenta'
        if(any(isat.ne.isat_ref)) Stop 'bsr_H: inconsistent target spins'
        if(any(ipat.ne.ipat_ref)) Stop 'bsr_H: inconsistent target parities'
       end if

       Deallocate(enat,lat,isat,ipat,coeff)

!      An h.nnn file can contain more than one symmetry block.

  100  Continue
       read(in) lrgl,nspn,npty,nchan,mnp2,more0
       nstat=abs(mnp2)
       if(nchan.le.0 .or. nstat.le.0) Stop 'bsr_H: invalid block dimensions'

       Allocate(nconat(nast),l2p(nchan),jkch(nchan))
       Allocate(cf(nchan,nchan,lamax),value(nstat),wmat(nchan,nstat))
       read(in) nconat
       read(in) l2p,jkch
       read(in) cf
       read(in) value
       read(in) wmat

       if(sum(nconat).ne.nchan) &
        Stop 'bsr_H_jk: SUM(NCONAT) does not equal NCHAN'
       if(any(l2p.lt.0) .or. any(l2p.ge.lrang2)) &
        Stop 'bsr_H_jk: channel l is outside the global LRANG2 range'

!      The pole energies are considered separately in each symmetry block.
!      By default, if they are in decreasing order, we reverse the order.
!      If the energies are not in increasing or decreasing order,
!      the program will stop.
       ascending=.true.; descending=.true.
       Do k=2,nstat
        if(value(k).lt.value(k-1)) ascending=.false.
        if(value(k).gt.value(k-1)) descending=.false.
       End do
       if(.not.ascending .and. .not.descending) &
        Stop 'bsr_H_jk: eigenvalues are not monotonic'
       if(reverse_poles.eq.-1) then
        reverse_this=descending .and. .not.ascending
       else
        reverse_this=reverse_poles.ne.0
       end if

!      MORE is set to 1 if another block, in the current h.nnn file or a
!      subsequent file, remains to be read. It is set to 0 once the
!      final block is reached.
       if(more0.gt.0 .or. ih.lt.last_file) then
        more=1
       else
        more=0
       end if

       write(ipri,'(a,a)') 'file ',trim(af)
       write(ipri,'(a,6(1x,i0))') 'block written:', &
            lrgl,nspn,npty,nchan,nstat,more
       write(ipri,'(a,l1)') '  reverse_this: ',reverse_this
       write(ipri,'(a,2(1x,es15.7))') '  eigen range input:', &
            value(1),value(nstat)
       write(ipri,'(a,1000(1x,i0))') '  NCONAT:', &
            (nconat(i),i=1,nast)
       write(ipri,'(a,1000(1x,i0))') '  L2P   :', &
            (l2p(i),i=1,nchan)
       write(ipri,'(a,1000(1x,i0))') '  jkch  :', &
            (jkch(i),i=1,nchan)
       write(ipri,*)

!      Write the six records for this symmetry block.
       write(iout) lrgl,nspn,npty,nchan,nstat,more
       write(iout) nconat
       write(iout) l2p
       write(iout) cf

!      Each column of WMAT contains the surface amplitudes for one pole.
!      If the pole energies are reversed, so are the columns of WMAT in
!      the same order.
       if(reverse_this) then
        write(iout) (value(nstat-k+1),k=1,nstat)
        write(iout) ((wmat(i,nstat-k+1),i=1,nchan),k=1,nstat)
       else
        write(iout) value
        write(iout) wmat
       end if

       Deallocate(nconat,l2p,jkch,cf,value,wmat)
       if(more0.gt.0) go to 100
       Close(in)
      End do

      Close(iout)
      write(*,'(a,i0,a)') 'bsr_H: merged ',nfiles,' files into H'
      write(ipri,'(a,i0,a)') 'bsr_H_jk_v3: merged ',nfiles,' files into H'
      Close(ipri)

      End
