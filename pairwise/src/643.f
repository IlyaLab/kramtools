C      ALGORITHM 643, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 19, NO. 4, DECEMBER, 1993, PP. 484-488.
c-----------------------------------------------------------------------
c  Name:       FEXACT
c
c  Purpose:    Computes Fisher's exact test probabilities and a hybrid
c              approximation to Fisher exact test probabilities for a
c              contingency table using the network algorithm.
c
c  Usage:      CALL FEXACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
c                          EMIN, PRT, PRE)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     NCOL   - The number of columns in the table.  (Input)
c     TABLE  - NROW by NCOL matrix containing the contingency table.
c              (Input)
c     LDTABL - Leading dimension of TABLE exactly as specified in the
c              dimension statement in the calling program.  (Input)
c     EXPECT - Expected value used in the hybrid algorithm for
c              deciding when to use asymptotic theory probabilities.
c              (Input)
c              If EXPECT .LE. 0.0 then asymptotic theory probabilities
c              are not used and Fisher exact test probabilities are
c              computed.  Otherwise, if PERCNT or more of the cells in
c              the remaining table have estimated expected values of
c              EXPECT or more, with no remaining cell having expected
c              value less than EMIN, then asymptotic chi-squared
c              probabilities are used.  See the algorithm section of the
c              manual document for details.  Use EXPECT = 5.0 to obtain
c              the 'Cochran' condition.
c     PERCNT - Percentage of remaining cells that must have estimated
c              expected  values greater than EXPECT before asymptotic
c              probabilities can be used.  (Input)
c              See argument EXPECT for details.  Use PERCNT = 80.0 to
c              obtain the 'Cochran' condition.
c     EMIN   - Minimum cell estimated expected value allowed for
c              asymptotic chi-squared probabilities to be used.  (Input)
c              See argument EXPECT for details.  Use EMIN = 1.0 to
c              obtain the 'Cochran' condition.
c     PRT    - Probability of the observed table for fixed marginal
c              totals.  (Output)
c     PRE    - Table p-value.  (Output)
c              PRE is the probability of a more extreme table, where
c              'extreme' is in a probabilistic sense.
c              If EXPECT .LT. 0 then the Fisher exact probability
c              is returned.  Otherwise, an approximation to the
c              Fisher exact probability is computed based upon
c              asymptotic chi-squared probabilities for ``large''
c              table expected values.  The user defines ``large''
c              through the arguments EXPECT, PERCNT, and EMIN.
c
c  Remarks:
c  1. For many problems one megabyte or more of workspace can be 
c     required.  If the environment supports it, the user should begin 
c     by increasing the workspace used to 200,000 units. 
c
c  2. In FEXACT, LDSTP = 30*LDKEY.  The proportion of table space used 
c     by STP may be changed by changing the line MULT = 30 below to 
c     another value.
c
c  3. FEXACT may be converted to single precision by setting IREAL = 3,
c     and converting all DOUBLE PRECISION specifications (except the 
c     specifications for RWRK, IWRK, and DWRK) to REAL.  This will 
c     require changing the names and specifications of the intrinsic
c     functions ALOG, AMAX1, AMIN1, EXP, and REAL.  In addition, the
c     machine specific constants will need to be changed, and the name
c     DWRK will need to be changed to RWRK in the call to F2XACT.
c
c  4. Machine specific constants are specified and documented in F2XACT.
c     A missing value code is specified in both FEXACT and F2XACT.
c
c  5. Although not a restriction, is is not generally practical to call
c     this routine with large tables which are not sparse and in
c     which the 'hybrid' algorithm has little effect.  For example,
c     although it is feasible to compute exact probabilities for the
c     table
c            1 8 5 4 4 2 2
c            5 3 3 4 3 1 0
c           10 1 4 0 0 0 0,
c     computing exact probabilities for a similar table which has been
c     enlarged by the addition of an extra row (or column) may not be
c     feasible.
c-----------------------------------------------------------------------
      subroutine fexact (nrow, ncol, table, ldtabl, expect, percnt,
     &                   emin, prt, pre,
     &                   iwkmax, rwrk, dwrk, iwrk )
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, ldtabl, iwkmax
      double precision expect, percnt, emin, prt, pre, table(ldtabl,*)
c***********************************************************************
c                                  Workspace is now provided by caller.
c***********************************************************************
      real       rwrk(iwkmax)
      double precision dwrk(iwkmax/2)
      integer    iwrk(iwkmax/2)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, i1, i10, i2, i3, i3a, i3b, i3c, i4, i5, i6, i7,
     &           i8, i9, i9a, iiwk, ireal, irwk, iwkpt,
     &           j, k, kk, ldkey, ldstp, mult, nco, nro,
     &           ntot, numb
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  max0
      integer    max0
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f2xact
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   iwork
      integer    iwork
c***********************************************************************
c                                  To increase the length of the table
c                                  of paste path lengths relative to the
c                                  length of the hash table,  increase
c                                  MULT 
c***********************************************************************
      mult   = 30
c***********************************************************************
c                                  Set IREAL = 4 for DOUBLE PRECISION
c                                  Set IREAL = 3 for SINGLE PRECISION
c***********************************************************************
      ireal  = 4
c***********************************************************************
c                                  AMISS is a missing value indicator
c                                  which is returned when the
c                                  probability is not defined.
c***********************************************************************
      amiss = -12345.0d0
c
      iwkpt  = 1
c
      if (nrow .gt. ldtabl) then
         call prterr (1, 'NROW must be less than or equal to '//
     &               'LDTABL.')
      end if
      ntot = 0
      do 20  i=1, nrow
         do 10  j=1, ncol
            if (table(i,j) .lt. 0) then
               call prterr (2, 'All elements of TABLE must '//
     &                     'be positive.')
            end if
            ntot = ntot + table(i,j)
   10    continue
   20 continue
      if (ntot .eq. 0) then
         call prterr (3, 'All elements of TABLE are zero.  '//
     &               'PRT and PRE are set to missing values '//
     &               '(NaN, not a number).')
         prt = amiss
         pre = amiss
         go to 9000
      end if
c
      nco = max0(nrow,ncol)
      nro = nrow + ncol - nco
      k   = nrow + ncol + 1
      kk  = k*max0(nrow,ncol)
c
      i1   = iwork(iwkmax,iwkpt,ntot+1,ireal)
      i2   = iwork(iwkmax,iwkpt,nco,2)
      i3   = iwork(iwkmax,iwkpt,nco,2)
      i3a  = iwork(iwkmax,iwkpt,nco,2)
      i3b  = iwork(iwkmax,iwkpt,nro,2)
      i3c  = iwork(iwkmax,iwkpt,nro,2)
      iiwk = iwork(iwkmax,iwkpt,max0(5*k+2*kk,800+7*max0(nrow,ncol)),2)
      irwk = iwork(iwkmax,iwkpt,max0(400+max0(nrow,ncol)+1,k),ireal)
c                                  Double precision
      if (ireal .eq. 4) then
         numb  = 18 + 10*mult
         ldkey = (iwkmax-iwkpt+1)/numb
      else
c                                  Real workspace
         numb  = 12 + 8*mult
         ldkey = (iwkmax-iwkpt+1)/numb
      end if
c
      ldstp = mult*ldkey
      i4    = iwork(iwkmax,iwkpt,2*ldkey,2)
      i5    = iwork(iwkmax,iwkpt,2*ldkey,2)
      i6    = iwork(iwkmax,iwkpt,2*ldstp,ireal)
      i7    = iwork(iwkmax,iwkpt,6*ldstp,2)
      i8    = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i9    = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i9a   = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i10   = iwork(iwkmax,iwkpt,2*ldkey,2)
c***********************************************************************
c                                  To convert to double precision,
c                                  change RWRK to WWRK in the next CALL
c***********************************************************************
c
      call f2xact (nrow, ncol, table, ldtabl, expect, percnt, emin,
     &             prt, pre, dwrk(i1), iwrk(i2), iwrk(i3), iwrk(i3a),
     &             iwrk(i3b), iwrk(i3c), iwrk(i4), ldkey, iwrk(i5),
     &             dwrk(i6), ldstp, iwrk(i7), dwrk(i8), dwrk(i9),
     &             dwrk(i9a), iwrk(i10), iwrk(iiwk), dwrk(irwk))
c
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F2XACT
c
c  Purpose:    Computes Fisher's exact test for a contingency table,
c              routine with workspace variables specified.
c
c  Usage:      CALL F2XACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
c                          EMIN, PRT, PRE, FACT, ICO, IRO, KYY, IDIF,
c                          IRN, KEY, LDKEY, IPOIN, STP, LDSTP, IFRQ,
c                          DLP, DSP, TM, KEY2, IWK, RWK)
c-----------------------------------------------------------------------
      subroutine f2xact (nrow, ncol, table, ldtabl, expect, percnt,
     &                   emin, prt, pre, fact, ico, iro, kyy, idif,
     &                   irn, key, ldkey, ipoin, stp, ldstp, ifrq,
     &                   dlp, dsp, tm, key2, iwk, rwk)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, ldtabl, ldkey, ldstp, ico(*), iro(*),
     &           kyy(*), idif(*), irn(*), key(*), ipoin(*), ifrq(*),
     &           key2(*), iwk(*)
      double precision expect, percnt, emin, prt, pre, table(ldtabl,*),
     &           fact(0:*), stp(*), dlp(*), dsp(*), tm(*), rwk(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, i31, i310, i311, i32, i33, i34, i35, i36, i37,
     &           i38, i39, i41, i42, i43, i44, i45, i46, i47, i48,
     &           iflag, ifreq, ii, ikkey, ikstp, ikstp2, ipn, ipo,
     &           itmp, itop, itp, j, jkey, jstp, jstp2, jstp3, jstp4,
     &           k, k1, kb, kd, kmax, ks, kval, last, n, ncell, nco,
     &           nrb, nro, nro2, ntot, ifault, imax
      double precision dd, ddf, df, drn, dro, dspt, emn, obs, obs2, 
     &           obs3, pastp, pv, tmp, tol
      logical    chisq, ipsh
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog, dmax1, dmin1, dexp, max0, min0, mod, nint, dble
      integer    max0, min0, mod, nint
      double precision dlog, dmax1, dmin1, dexp, dble
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f3xact, f4xact, f5xact, f6xact, f7xact, isort
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   f9xact, gammds
      double precision f9xact, gammds
c***********************************************************************
c                                  IMAX is the largest representable
c                                  integer on the machine
c***********************************************************************
      data imax/2147483647/
c***********************************************************************
c                                  AMISS is a missing value indicator
c                                  which is returned when the
c                                  probability is not defined.
c***********************************************************************
      data amiss/-12345.0d0/
c***********************************************************************
c                                  TOL is chosen as the square root of
c                                  the smallest relative spacing
c***********************************************************************
      data tol/3.45254d-07/
c***********************************************************************
c                                  EMX is a large positive value used 
c                                  in comparing expected values
c***********************************************************************
      data emx/1.0d30/
c                                  Initialize KEY array
      do 10  i=1, 2*ldkey
         key(i)  = -9999
         key2(i) = -9999
   10 continue
c                                  Initialize parameters
      pre  = 0.0
      itop = 0
      if (expect .gt. 0.0d0) then
         emn = emin
      else
         emn = emx
      end if
c                                  Initialize pointers for workspace
      k = max0(nrow,ncol)
c                                  f3xact
      i31  = 1
      i32  = i31 + k
      i33  = i32 + k
      i34  = i33 + k
      i35  = i34 + k
      i36  = i35 + k
      i37  = i36 + k
      i38  = i37 + k
      i39  = i38 + 400
      i310 = 1
      i311 = 401
c                                  f4xact
      k   = nrow + ncol + 1
      i41 = 1
      i42 = i41 + k
      i43 = i42 + k
      i44 = i43 + k
      i45 = i44 + k
      i46 = i45 + k
      i47 = i46 + k*max0(nrow,ncol)
      i48 = 1
c                                  Check table dimensions
      if (nrow .gt. ldtabl) then
         call prterr (1, 'NROW must be less than or equal to '//
     &               'LDTABL.')
      end if
      if (ncol .le. 1) then
         call prterr (4, 'NCOL must be greater than 1.0.')
      end if
c                                  Compute row marginals and total
      ntot = 0
      do 30  i=1, nrow
         iro(i) = 0
         do 20  j=1, ncol
            if (table(i,j) .lt. -0.0001d0) then
               call prterr (2, 'All elements of TABLE must be '//
     &                     'positive.')
            end if
            iro(i) = iro(i) + nint(table(i,j))
            ntot   = ntot + nint(table(i,j))
   20    continue
   30 continue
c
      if (ntot .eq. 0) then
         call prterr (3, 'All elements of TABLE are zero.  '//
     &               'PRT and PRE are set to missing values '//
     &               '(NaN, not a number).')
         prt = amiss
         pre = amiss
         go to 9000
      end if
c                                  Column marginals
      do 50  i=1, ncol
         ico(i) = 0
         do 40  j=1, nrow
            ico(i) = ico(i) + nint(table(j,i))
   40    continue
   50 continue
c                                  sort
      call isort (nrow, iro)
      call isort (ncol, ico)
c                                  Determine row and column marginals
c
      if (nrow .gt. ncol) then
         nro = ncol
         nco = nrow
c                                  Interchange row and column marginals
         do 60  i=1, nrow
            itmp = iro(i)
            if (i .le. ncol) iro(i) = ico(i)
            ico(i) = itmp
   60    continue
      else
         nro = nrow
         nco = ncol
      end if
c
c                                  Get multiplers for stack
      kyy(1) = 1
      do 70  i=2, nro
c                                  Hash table multipliers
         if (iro(i-1)+1 .le. imax/kyy(i-1)) then
            kyy(i) = kyy(i-1)*(iro(i-1)+1)
            j      = j/kyy(i-1)
         else
            call prterr (5, 'The hash table key cannot be computed'//
     &                  ' because the largest key is larger than the'//
     &                  ' largest representable integer.  The '//
     &                  'algorithm cannot proceed.')
         end if
   70 continue
c                                  Maximum product
      if (iro(nro-1)+1 .le. imax/kyy(nro-1)) then
         kmax = (iro(nro)+1)*kyy(nro-1)
      else
         call prterr (5, 'The hash table key cannot be computed'//
     &               ' because the largest key is larger than the'//
     &               ' largest representable integer.  The '//
     &               'algorithm cannot proceed.')
         go to 9000
      end if
c                                  Compute log factorials
      fact(0) = 0.0d0
      fact(1) = 0.0d0
      fact(2) = dlog(2.0d0)
      do 80  i=3, ntot, 2
         fact(i) = fact(i-1) + dlog(dble(i))
         j       = i + 1
         if (j .le. ntot) fact(j) = fact(i) + fact(2) + fact(j/2) -
     &       fact(j/2-1)
   80 continue
c                                  Compute observed path length: OBS
      obs  = tol
      ntot = 0
      do 100  j=1, nco
         dd = 0.0
         do 90  i=1, nro
            if (nrow .le. ncol) then
               dd   = dd + fact(nint(table(i,j)))
               ntot = ntot + nint(table(i,j))
            else
               dd   = dd + fact(nint(table(j,i)))
               ntot = ntot + nint(table(j,i))
            end if
   90    continue
         obs = obs + fact(ico(j)) - dd
  100 continue
c                                  Denominator of observed table: DRO
      dro = f9xact(nro,ntot,iro,fact)
      prt = dexp(obs-dro)
c                                  Initialize pointers
      k        = nco
      last     = ldkey + 1
      jkey     = ldkey + 1
      jstp     = ldstp + 1
      jstp2    = 3*ldstp + 1
      jstp3    = 4*ldstp + 1
      jstp4    = 5*ldstp + 1
      ikkey    = 0
      ikstp    = 0
      ikstp2   = 2*ldstp
      ipo      = 1
      ipoin(1) = 1
      stp(1)   = 0.0
      ifrq(1)  = 1
      ifrq(ikstp2+1) = -1
c
  110 kb = nco - k + 1
      ks   = 0
      n    = ico(kb)
      kd   = nro + 1
      kmax = nro
c                                  IDIF is the difference in going to th
c                                  daughter
      do 120  i=1, nro
         idif(i) = 0
  120 continue
c                                  Generate the first daughter
  130 kd = kd - 1
      ntot     = min0(n,iro(kd))
      idif(kd) = ntot
      if (idif(kmax) .eq. 0) kmax = kmax - 1
      n = n - ntot
      if (n.gt.0 .and. kd.ne.1) go to 130
      if (n .ne. 0) go to 310
c
      k1   = k - 1
      n    = ico(kb)
      ntot = 0
      do 140  i=kb + 1, nco
         ntot = ntot + ico(i)
  140 continue
c                                  Arc to daughter length=ICO(KB)
  150 do 160  i=1, nro
         irn(i) = iro(i) - idif(i)
  160 continue
c                                  Sort irn
      if (k1 .gt. 1) then
         if (nro .eq. 2) then
            if (irn(1) .gt. irn(2)) then
               ii     = irn(1)
               irn(1) = irn(2)
               irn(2) = ii
            end if
         else if (nro .eq. 3) then
            ii = irn(1)
            if (ii .gt. irn(3)) then
               if (ii .gt. irn(2)) then
                  if (irn(2) .gt. irn(3)) then
                     irn(1) = irn(3)
                     irn(3) = ii
                  else
                     irn(1) = irn(2)
                     irn(2) = irn(3)
                     irn(3) = ii
                  end if
               else
                  irn(1) = irn(3)
                  irn(3) = irn(2)
                  irn(2) = ii
               end if
            else if (ii .gt. irn(2)) then
               irn(1) = irn(2)
               irn(2) = ii
            else if (irn(2) .gt. irn(3)) then
               ii     = irn(2)
               irn(2) = irn(3)
               irn(3) = ii
            end if
         else
            do 180  j=2, nro
               i  = j - 1
               ii = irn(j)
  170          if (ii .lt. irn(i)) then
                  irn(i+1) = irn(i)
                  i        = i - 1
                  if (i .gt. 0) go to 170
               end if
               irn(i+1) = ii
  180       continue
         end if
c                                  Adjust start for zero
         do 190  i=1, nro
            if (irn(i) .ne. 0) go to 200
  190    continue
  200    nrb = i
         nro2 = nro - i + 1
      else
         nrb  = 1
         nro2 = nro
      end if
c                                  Some table values
      ddf = f9xact(nro,n,idif,fact)
      drn = f9xact(nro2,ntot,irn(nrb),fact) - dro + ddf
c                                  Get hash value
      if (k1 .gt. 1) then
         kval = irn(1) + irn(2)*kyy(2)
         do 210  i=3, nro
            kval = kval + irn(i)*kyy(i)
  210    continue
c                                  Get hash table entry
         i = mod(kval,2*ldkey) + 1
c                                  Search for unused location
         do 220  itp=i, 2*ldkey
            ii = key2(itp)
            if (ii .eq. kval) then
               go to 240
            else if (ii .lt. 0) then
               key2(itp) = kval
               dlp(itp)  = 1.0d0
               dsp(itp)  = 1.0d0
               go to 240
            end if
  220    continue
c
         do 230  itp=1, i - 1
            ii = key2(itp)
            if (ii .eq. kval) then
               go to 240
            else if (ii .lt. 0) then
               key2(itp) = kval
               dlp(itp)  = 1.0
               go to 240
            end if
  230    continue
c
         call prterr (6, 'LDKEY is too small.  It is not possible to '//
     &               'give thevalue of LDKEY required, but you could '//
     &               'try doubling LDKEY (and possibly LDSTP).')
      end if
c
  240 ipsh = .true.
c                                  Recover pastp
      ipn   = ipoin(ipo+ikkey)
      pastp = stp(ipn+ikstp)
      ifreq = ifrq(ipn+ikstp)
c                                  Compute shortest and longest path
      if (k1 .gt. 1) then
         obs2 = obs - fact(ico(kb+1)) - fact(ico(kb+2)) - ddf
         do 250  i=3, k1
            obs2 = obs2 - fact(ico(kb+i))
  250    continue
c
         if (dlp(itp) .gt. 0.0d0) then
            dspt = obs - obs2 - ddf
c                                  Compute longest path
            dlp(itp) = 0.0d0
            call f3xact (nro2, irn(nrb), k1, ico(kb+1), dlp(itp),
     &                   ntot, fact, iwk(i31), iwk(i32), iwk(i33),
     &                   iwk(i34), iwk(i35), iwk(i36), iwk(i37),
     &                   iwk(i38), iwk(i39), rwk(i310), rwk(i311), tol)
            dlp(itp) = dmin1(0.0d0,dlp(itp))
c                                  Compute shortest path
            dsp(itp) = dspt
            call f4xact (nro2, irn(nrb), k1, ico(kb+1), dsp(itp),
     &                   fact, iwk(i47), iwk(i41), iwk(i42), iwk(i43),
     &                   iwk(i44), iwk(i45), iwk(i46), rwk(i48), tol)
            dsp(itp) = dmin1(0.0d0,dsp(itp)-dspt)
c                                  Use chi-squared approximation?
            if (dble(irn(nrb)*ico(kb+1))/dble(ntot) .gt. emn) then
               ncell = 0.0
               do 270  i=1, nro2
                  do 260  j=1, k1
                     if (irn(nrb+i-1)*ico(kb+j) .ge. ntot*expect) then
                        ncell = ncell + 1
                     end if
  260             continue
  270          continue
               if (ncell*100 .ge. k1*nro2*percnt) then
                  tmp = 0.0
                  do 280  i=1, nro2
                     tmp = tmp + fact(irn(nrb+i-1)) -
     &                     fact(irn(nrb+i-1)-1)
  280             continue
                  tmp = tmp*(k1-1)
                  do 290  j=1, k1
                     tmp = tmp + (nro2-1)*(fact(ico(kb+j))-fact(ico(kb+
     &                     j)-1))
  290             continue
                  df      = (nro2-1)*(k1-1)
                  tmp     = tmp + df*1.83787706640934548356065947281d0
                  tmp     = tmp - (nro2*k1-1)*(fact(ntot)-fact(ntot-1))
                  tm(itp) = -2.0d0*(obs-dro) - tmp
               else
c                                  tm(itp) set to a flag value
                  tm(itp) = -9876.0d0
               end if
            else
               tm(itp) = -9876.0d0
            end if
         end if
         obs3 = obs2 - dlp(itp)
         obs2 = obs2 - dsp(itp)
         if (tm(itp) .eq. -9876.0d0) then
            chisq = .false.
         else
            chisq = .true.
            tmp   = tm(itp)
         end if
      else
         obs2 = obs - drn - dro
         obs3 = obs2
      end if
c                                  Process node with new PASTP
  300 if (pastp .le. obs3) then
c                                  Update pre
         pre = pre + dble(ifreq)*dexp(pastp+drn)
c
      else if (pastp .lt. obs2) then
         if (chisq) then
            df  = (nro2-1)*(k1-1)
            pv  = 1.0 - gammds(dmax1(0.0d0,tmp+2.0d0*(pastp+drn))/
     &            2.0d0,df/2.0d0,ifault)
            pre = pre + dble(ifreq)*dexp(pastp+drn)*pv
         else
c                                  Put daughter on queue
            call f5xact (pastp+ddf, tol, kval, key(jkey), ldkey,
     &                   ipoin(jkey), stp(jstp), ldstp, ifrq(jstp),
     &                   ifrq(jstp2), ifrq(jstp3), ifrq(jstp4), ifreq,
     &                   itop, ipsh)
            ipsh = .false.
         end if
      end if
c                                  Get next PASTP on chain
      ipn = ifrq(ipn+ikstp2)
      if (ipn .gt. 0) then
         pastp = stp(ipn+ikstp)
         ifreq = ifrq(ipn+ikstp)
         go to 300
      end if
c                                  Generate a new daughter node
      call f7xact (kmax, iro, idif, kd, ks, iflag)
      if (iflag .ne. 1) go to 150
c                                  Go get a new mother from stage K
  310 iflag = 1
      call f6xact (nro, iro, iflag, kyy, key(ikkey+1), ldkey, last,
     &             ipo)
c                                  Update pointers
      if (iflag .eq. 3) then
         k      = k - 1
         itop   = 0
         ikkey  = jkey - 1
         ikstp  = jstp - 1
         ikstp2 = jstp2 - 1
         jkey   = ldkey - jkey + 2
         jstp   = ldstp - jstp + 2
         jstp2  = 2*ldstp + jstp
         do 320  i=1, 2*ldkey
            key2(i) = -9999
  320    continue
         if (k .ge. 2) go to 310
      else
         go to 110
      end if
c
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F3XACT
c
c  Purpose:    Computes the shortest path length for a given table.
c
c  Usage:      CALL F3XACT (NROW, IROW, NCOL, ICOL, DLP, MM, FACT, ICO,
c                          IRO, IT, LB, NR, NT, NU, ITC, IST, STV, ALEN,
c                          TOL)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     IROW   - Vector of length NROW containing the row sums for the
c              table.  (Input)
c     NCOL   - The number of columns in the table.  (Input)
c     ICOL   - Vector of length K containing the column sums for the
c              table.  (Input)
c     DLP    - The longest path for the table.  (Output)
c     MM     - The total count in the table.  (Output)
c     FACT   - Vector containing the logarithms of factorials.  (Input)
c     ICO    - Work vector of length MAX(NROW,NCOL).
c     IRO    - Work vector of length MAX(NROW,NCOL).
c     IT     - Work vector of length MAX(NROW,NCOL).
c     LB     - Work vector of length MAX(NROW,NCOL).
c     NR     - Work vector of length MAX(NROW,NCOL).
c     NT     - Work vector of length MAX(NROW,NCOL).
c     NU     - Work vector of length MAX(NROW,NCOL).
c     ITC    - Work vector of length 400.
c     IST    - Work vector of length 400.
c     STV    - Work vector of length 400.
c     ALEN   - Work vector of length MAX(NROW,NCOL).
c     TOL    - Tolerance.  (Input)
c-----------------------------------------------------------------------
      subroutine f3xact (nrow, irow, ncol, icol, dlp, mm, fact, ico,
     &                   iro, it, lb, nr, nt, nu, itc, ist, stv, alen,
     &                   tol)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, mm, irow(*), icol(*), ico(*), iro(*),
     &           it(*), lb(*), nr(*), nt(*), nu(*), itc(*), ist(*)
      double precision dlp, tol, fact(0:*), stv(*), alen(0:*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ic1, ic2, ii, ipn, irl, itp, k, key, ks, kyy, lev,
     &           n11, n12, nc1, nc1s, nco, nct, nn, nn1, nr1, nro, nrt
      double precision v, val, vmn
      logical    xmin
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    ldst, nitc, nst
      save       ldst, nitc, nst
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  dmin1, int, mod, dble
      integer    int, mod
      double precision dmin1, dble
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f10act, isort
c
      data ldst/200/, nst/0/, nitc/0/
c
      do 10  i=0, ncol
         alen(i) = 0.0
   10 continue
      do 20  i=1, 400
         ist(i) = -1
   20 continue
c                                  nrow is 1
      if (nrow .le. 1) then
         if (nrow .gt. 0) then
            dlp = dlp - fact(icol(1))
            do 30  i=2, ncol
               dlp = dlp - fact(icol(i))
   30       continue
         end if
         go to 9000
      end if
c                                  ncol is 1
      if (ncol .le. 1) then
         if (ncol .gt. 0) then
            dlp = dlp - fact(irow(1)) - fact(irow(2))
            do 40  i=3, nrow
               dlp = dlp - fact(irow(i))
   40       continue
         end if
         go to 9000
      end if
c                                  2 by 2 table
      if (nrow*ncol .eq. 4) then
         n11 = (irow(1)+1)*(icol(1)+1)/(mm+2)
         n12 = irow(1) - n11
         dlp = dlp - fact(n11) - fact(n12) - fact(icol(1)-n11) -
     &         fact(icol(2)-n12)
         go to 9000
      end if
c                                  Test for optimal table
      val  = 0.0
      xmin = .false.
      if (irow(nrow) .le. irow(1)+ncol) then
         call f10act (nrow, irow, ncol, icol, val, xmin, fact, lb, nu,
     &                nr)
      end if
      if (.not.xmin) then
         if (icol(ncol) .le. icol(1)+nrow) then
            call f10act (ncol, icol, nrow, irow, val, xmin, fact, lb,
     &                   nu, nr)
         end if
      end if
c
      if (xmin) then
         dlp = dlp - val
         go to 9000
      end if
c                                  Setup for dynamic programming
      nn = mm
c                                  Minimize ncol
      if (nrow .ge. ncol) then
         nro = nrow
         nco = ncol
c
         do 50  i=1, nrow
            iro(i) = irow(i)
   50    continue
c
         ico(1) = icol(1)
         nt(1)  = nn - ico(1)
         do 60  i=2, ncol
            ico(i) = icol(i)
            nt(i)  = nt(i-1) - ico(i)
   60    continue
      else
         nro = ncol
         nco = nrow
c
         ico(1) = irow(1)
         nt(1)  = nn - ico(1)
         do 70  i=2, nrow
            ico(i) = irow(i)
            nt(i)  = nt(i-1) - ico(i)
   70    continue
c
         do 80  i=1, ncol
            iro(i) = icol(i)
   80    continue
      end if
c                                  Initialize pointers
      vmn  = 1.0d10
      nc1s = nco - 1
      irl  = 1
      ks   = 0
      k    = ldst
      kyy  = ico(nco) + 1
      go to 100
c                                  Test for optimality
   90 xmin = .false.
      if (iro(nro) .le. iro(irl)+nco) then
         call f10act (nro, iro(irl), nco, ico, val, xmin, fact, lb,
     &                nu, nr)
      end if
      if (.not.xmin) then
         if (ico(nco) .le. ico(1)+nro) then
            call f10act (nco, ico, nro, iro(irl), val, xmin, fact, lb,
     &                   nu, nr)
         end if
      end if
c
      if (xmin) then
         if (val .lt. vmn) vmn = val
         go to 200
      end if
c                                  Setup to generate new node
  100 lev = 1
      nr1   = nro - 1
      nrt   = iro(irl)
      nct   = ico(1)
      lb(1) = int(dble((nrt+1)*(nct+1))/dble(nn+nr1*nc1s+1)-tol) - 1
      nu(1) = int(dble((nrt+nc1s)*(nct+nr1))/dble(nn+nr1+nc1s)) -
     &        lb(1) + 1
      nr(1) = nrt - lb(1)
c                                  Generate a node
  110 nu(lev) = nu(lev) - 1
      if (nu(lev) .eq. 0) then
         if (lev .eq. 1) go to 200
         lev = lev - 1
         go to 110
      end if
      lb(lev) = lb(lev) + 1
      nr(lev) = nr(lev) - 1
  120 alen(lev) = alen(lev-1) + fact(lb(lev))
      if (lev .lt. nc1s) then
         nn1     = nt(lev)
         nrt     = nr(lev)
         lev     = lev + 1
         nc1     = nco - lev
         nct     = ico(lev)
         lb(lev) = dble((nrt+1)*(nct+1))/dble(nn1+nr1*nc1+1) - tol
         nu(lev) = dble((nrt+nc1)*(nct+nr1))/dble(nn1+nr1+nc1) -
     &             lb(lev) + 1
         nr(lev) = nrt - lb(lev)
         go to 120
      end if
      alen(nco) = alen(lev) + fact(nr(lev))
      lb(nco)   = nr(lev)
c
      v = val + alen(nco)
      if (nro .eq. 2) then
c                                  Only 1 row left
         v = v + fact(ico(1)-lb(1)) + fact(ico(2)-lb(2))
         do 130  i=3, nco
            v = v + fact(ico(i)-lb(i))
  130    continue
         if (v .lt. vmn) vmn = v
      else if (nro.eq.3 .and. nco.eq.2) then
c                                  3 rows and 2 columns
         nn1 = nn - iro(irl) + 2
         ic1 = ico(1) - lb(1)
         ic2 = ico(2) - lb(2)
         n11 = (iro(irl+1)+1)*(ic1+1)/nn1
         n12 = iro(irl+1) - n11
         v   = v + fact(n11) + fact(n12) + fact(ic1-n11) +
     &         fact(ic2-n12)
         if (v .lt. vmn) vmn = v
      else
c                                  Column marginals are new node
         do 140  i=1, nco
            it(i) = ico(i) - lb(i)
  140    continue
c                                  Sort column marginals
         if (nco .eq. 2) then
            if (it(1) .gt. it(2)) then
               ii    = it(1)
               it(1) = it(2)
               it(2) = ii
            end if
         else if (nco .eq. 3) then
            ii = it(1)
            if (ii .gt. it(3)) then
               if (ii .gt. it(2)) then
                  if (it(2) .gt. it(3)) then
                     it(1) = it(3)
                     it(3) = ii
                  else
                     it(1) = it(2)
                     it(2) = it(3)
                     it(3) = ii
                  end if
               else
                  it(1) = it(3)
                  it(3) = it(2)
                  it(2) = ii
               end if
            else if (ii .gt. it(2)) then
               it(1) = it(2)
               it(2) = ii
            else if (it(2) .gt. it(3)) then
               ii    = it(2)
               it(2) = it(3)
               it(3) = ii
            end if
         else
            call isort (nco, it)
         end if
c                                  Compute hash value
         key = it(1)*kyy + it(2)
         do 150  i=3, nco
            key = it(i) + key*kyy
  150    continue
c                                  Table index
         ipn = mod(key,ldst) + 1
c                                  Find empty position
         ii = ks + ipn
         do 160  itp=ipn, ldst
            if (ist(ii) .lt. 0) then
               go to 180
            else if (ist(ii) .eq. key) then
               go to 190
            end if
            ii = ii + 1
  160    continue
c
         ii = ks + 1
         do 170  itp=1, ipn - 1
            if (ist(ii) .lt. 0) then
               go to 180
            else if (ist(ii) .eq. key) then
               go to 190
            end if
            ii = ii + 1
  170    continue
c
         call prterr (30, 'Stack length exceeded in f3xact.'//
     &               '  This problem should not occur.')
c                                  Push onto stack
  180    ist(ii) = key
         stv(ii) = v
         nst     = nst + 1
         ii      = nst + ks
         itc(ii) = itp
         go to 110
c                                  Marginals already on stack
  190    stv(ii) = dmin1(v,stv(ii))
      end if
      go to 110
c                                  Pop item from stack
  200 if (nitc .gt. 0) then
c                                  Stack index
         itp      = itc(nitc+k) + k
         nitc     = nitc - 1
         val      = stv(itp)
         key      = ist(itp)
         ist(itp) = -1
c                                  Compute marginals
         do 210  i=nco, 2, -1
            ico(i) = mod(key,kyy)
            key    = key/kyy
  210    continue
         ico(1) = key
c                                  Set up nt array
         nt(1) = nn - ico(1)
         do 220  i=2, nco
            nt(i) = nt(i-1) - ico(i)
  220    continue
         go to 90
c
      else if (nro.gt.2 .and. nst.gt.0) then
c                                  Go to next level
         nitc = nst
         nst  = 0
         k    = ks
         ks   = ldst - ks
         nn   = nn - iro(irl)
         irl  = irl + 1
         nro  = nro - 1
         go to 200
      end if
c
      dlp = dlp - vmn
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F4XACT
c
c  Purpose:    Computes the longest path length for a given table.
c
c  Usage:      CALL F4XACT (NROW, IROW, NCOL, ICOL, DSP, FACT, ICSTK,
c                          NCSTK, LSTK, MSTK, NSTK, NRSTK, IRSTK, YSTK,
c                          TOL)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     IROW   - Vector of length NROW containing the row sums for the
c              table.  (Input)
c     NCOL   - The number of columns in the table.  (Input)
c     ICOL   - Vector of length K containing the column sums for the
c              table.  (Input)
c     DSP    - The shortest path for the table.  (Output)
c     FACT   - Vector containing the logarithms of factorials.  (Input)
c     ICSTK  - NCOL by NROW+NCOL+1 work array.
c     NCSTK  - Work vector of length NROW+NCOL+1.
c     LSTK   - Work vector of length NROW+NCOL+1.
c     MSTK   - Work vector of length NROW+NCOL+1.
c     NSTK   - Work vector of length NROW+NCOL+1.
c     NRSTK  - Work vector of length NROW+NCOL+1.
c     IRSTK  - NROW by MAX(NROW,NCOL) work array.
c     YSTK   - Work vector of length NROW+NCOL+1.
c     TOL    - Tolerance.  (Input)
c-----------------------------------------------------------------------
      subroutine f4xact (nrow, irow, ncol, icol, dsp, fact, icstk,
     &                   ncstk, lstk, mstk, nstk, nrstk, irstk, ystk,
     &                   tol)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, irow(*), icol(*), icstk(ncol,*),
     &           ncstk(*), lstk(*), mstk(*), nstk(*), nrstk(*),
     &           irstk(nrow,*)
      double precision dsp, tol, fact(0:*), ystk(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ic1, ict, ir1, irt, istk, j, k, l, m, mn, n, nco,
     &           nro
      double precision amx, y
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   f11act, f8xact
c                                  Take care of the easy cases firstkt
      if (nrow .eq. 1) then
         do 10  i=1, ncol
            dsp = dsp - fact(icol(i))
   10    continue
         go to 9000
      end if
c
      if (ncol .eq. 1) then
         do 20  i=1, nrow
            dsp = dsp - fact(irow(i))
   20    continue
         go to 9000
      end if
c
      if (nrow*ncol .eq. 4) then
         if (irow(2) .le. icol(2)) then
            dsp = dsp - fact(irow(2)) - fact(icol(1)) -
     &            fact(icol(2)-irow(2))
         else
            dsp = dsp - fact(icol(2)) - fact(irow(1)) -
     &            fact(irow(2)-icol(2))
         end if
         go to 9000
      end if
c                                  initialization before loop
      do 30  i=1, nrow
         irstk(i,1) = irow(nrow-i+1)
   30 continue
c
      do 40  j=1, ncol
         icstk(j,1) = icol(ncol-j+1)
   40 continue
c
      nro      = nrow
      nco      = ncol
      nrstk(1) = nro
      ncstk(1) = nco
      ystk(1)  = 0.0
      y        = 0.0
      istk     = 1
      l        = 1
      amx      = 0.0
c
   50 ir1 = irstk(1,istk)
      ic1 = icstk(1,istk)
      if (ir1 .gt. ic1) then
         if (nro .ge. nco) then
            m = nco - 1
            n = 2
         else
            m = nro
            n = 1
         end if
      else if (ir1 .lt. ic1) then
         if (nro .le. nco) then
            m = nro - 1
            n = 1
         else
            m = nco
            n = 2
         end if
      else
         if (nro .le. nco) then
            m = nro - 1
            n = 1
         else
            m = nco - 1
            n = 2
         end if
      end if
c
   60 if (n .eq. 1) then
         i = l
         j = 1
      else
         i = 1
         j = l
      end if
c
      irt = irstk(i,istk)
      ict = icstk(j,istk)
      mn  = irt
      if (mn .gt. ict) mn = ict
      y = y + fact(mn)
      if (irt .eq. ict) then
         nro = nro - 1
         nco = nco - 1
         call f11act (irstk(1,istk), i, nro, irstk(1,istk+1))
         call f11act (icstk(1,istk), j, nco, icstk(1,istk+1))
      else if (irt .gt. ict) then
         nco = nco - 1
         call f11act (icstk(1,istk), j, nco, icstk(1,istk+1))
         call f8xact (irstk(1,istk), irt-ict, i, nro, irstk(1,istk+1))
      else
         nro = nro - 1
         call f11act (irstk(1,istk), i, nro, irstk(1,istk+1))
         call f8xact (icstk(1,istk), ict-irt, j, nco, icstk(1,istk+1))
      end if
c
      if (nro .eq. 1) then
         do 70  k=1, nco
            y = y + fact(icstk(k,istk+1))
   70    continue
         go to 90
      end if
c
      if (nco .eq. 1) then
         do 80  k=1, nro
            y = y + fact(irstk(k,istk+1))
   80    continue
         go to 90
      end if
c
      lstk(istk)  = l
      mstk(istk)  = m
      nstk(istk)  = n
      istk        = istk + 1
      nrstk(istk) = nro
      ncstk(istk) = nco
      ystk(istk)  = y
      l           = 1
      go to 50
c
   90 if (y .gt. amx) then
         amx = y
         if (dsp-amx .le. tol) then
            dsp = 0.0
            go to 9000
         end if
      end if
c
  100 istk = istk - 1
      if (istk .eq. 0) then
         dsp = dsp - amx
         if (dsp-amx .le. tol) dsp = 0.0
         go to 9000
      end if
      l = lstk(istk) + 1
c
  110 if (l .gt. mstk(istk)) go to 100
      n   = nstk(istk)
      nro = nrstk(istk)
      nco = ncstk(istk)
      y   = ystk(istk)
      if (n .eq. 1) then
         if (irstk(l,istk) .lt. irstk(l-1,istk)) go to 60
      else if (n .eq. 2) then
         if (icstk(l,istk) .lt. icstk(l-1,istk)) go to 60
      end if
c
      l = l + 1
      go to 110
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F5XACT
c
c  Purpose:    Put node on stack in network algorithm.
c
c  Usage:      CALL F5XACT (PASTP, TOL, KVAL, KEY, LDKEY, IPOIN, STP,
c                          LDSTP, IFRQ, NPOIN, NR, NL, IFREQ, ITOP,
c                          IPSH)
c
c  Arguments:
c     PASTP  - The past path length.  (Input)
c     TOL    - Tolerance for equivalence of past path lengths.  (Input)
c     KVAL   - Key value.  (Input)
c     KEY    - Vector of length LDKEY containing the key values.
c              (Input/output)
c     LDKEY  - Length of vector KEY.  (Input)
c     IPOIN  - Vector of length LDKEY pointing to the linked list
c              of past path lengths.  (Input/output)
c     STP    - Vector of length LSDTP containing the linked lists
c              of past path lengths.  (Input/output)
c     LDSTP  - Length of vector STP.  (Input)
c     IFRQ   - Vector of length LDSTP containing the past path
c              frequencies.  (Input/output)
c     NPOIN  - Vector of length LDSTP containing the pointers to
c              the next past path length.  (Input/output)
c     NR     - Vector of length LDSTP containing the right object
c              pointers in the tree of past path lengths.
c              (Input/output)
c     NL     - Vector of length LDSTP containing the left object
c              pointers in the tree of past path lengths.
c              (Input/output)
c     IFREQ  - Frequency of the current path length.  (Input)
c     ITOP   - Pointer to the top of STP.  (Input)
c     IPSH   - Option parameter.  (Input)
c              If IPSH is true, the past path length is found in the
c              table KEY.  Otherwise the location of the past path
c              length is assumed known and to have been found in
c              a previous call.
c-----------------------------------------------------------------------
      subroutine f5xact (pastp, tol, kval, key, ldkey, ipoin, stp,
     &                   ldstp, ifrq, npoin, nr, nl, ifreq, itop, ipsh)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    kval, ldkey, ldstp, ifreq, itop, key(*), ipoin(*),
     &           ifrq(*), npoin(*), nr(*), nl(*)
      double precision pastp, tol, stp(*)
      logical    ipsh
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    ipn, ird, itmp
      double precision test1, test2
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      integer    itp
      save       itp
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  mod
      integer    mod
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr
c
      if (ipsh) then
c                                  Convert KVAL to integer in range
c                                  1, ..., LDKEY.
         ird = mod(kval,ldkey) + 1
c                                  Search for an unused location
         do 10  itp=ird, ldkey
            if (key(itp) .eq. kval) go to 40
            if (key(itp) .lt. 0) go to 30
   10    continue
c
         do 20  itp=1, ird - 1
            if (key(itp) .eq. kval) go to 40
            if (key(itp) .lt. 0) go to 30
   20    continue
c                                  Return if KEY array is full
         call prterr(6, 'LDKEY is too small for this problem.  It is '//
     &               'not possible to estimate the value of LDKEY '//
     &               'required, but twice the current value may be '//
     &               'sufficient.')
c                                  Update KEY
   30    key(itp) = kval
         itop       = itop + 1
         ipoin(itp) = itop
c                                  Return if STP array full
         if (itop .gt. ldstp) then
            call prterr(7, 'LDSTP is too small for this problem.  It '//
     &                  'is not possible to estimate the value of '//
     &                  'LDSTP required, but twice the current value '//
     &                  'may be sufficient.')
         end if
c                                  Update STP, etc.
         npoin(itop) = -1
         nr(itop)    = -1
         nl(itop)    = -1
         stp(itop)   = pastp
         ifrq(itop)  = ifreq
         go to 9000
      end if
c                                  Find location, if any, of pastp
   40 ipn = ipoin(itp)
      test1 = pastp - tol
      test2 = pastp + tol
c
   50 if (stp(ipn) .lt. test1) then
         ipn = nl(ipn)
         if (ipn .gt. 0) go to 50
      else if (stp(ipn) .gt. test2) then
         ipn = nr(ipn)
         if (ipn .gt. 0) go to 50
      else
         ifrq(ipn) = ifrq(ipn) + ifreq
         go to 9000
      end if
c                                  Return if STP array full
      itop = itop + 1
      if (itop .gt. ldstp) then
         call prterr(7, 'LDSTP is too small for this problem.  It is '//
     &               'not possible to estimate the value of LDSTP '//
     &               'rerquired, but twice the current value may be '//
     &               'sufficient.')
         go to 9000
      end if
c                                  Find location to add value
      ipn  = ipoin(itp)
      itmp = ipn
   60 if (stp(ipn) .lt. test1) then
         itmp = ipn
         ipn  = nl(ipn)
         if (ipn .gt. 0) then
            go to 60
         else
            nl(itmp) = itop
         end if
      else if (stp(ipn) .gt. test2) then
         itmp = ipn
         ipn  = nr(ipn)
         if (ipn .gt. 0) then
            go to 60
         else
            nr(itmp) = itop
         end if
      end if
c                                  Update STP, etc.
      npoin(itop) = npoin(itmp)
      npoin(itmp) = itop
      stp(itop)   = pastp
      ifrq(itop)  = ifreq
      nl(itop)    = -1
      nr(itop)    = -1
c
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F6XACT
c
c  Purpose:    Pop a node off the stack.
c
c  Usage:      CALL F6XACT (NROW, IROW, IFLAG, KYY, KEY, LDKEY, LAST,
c                          IPN)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     IROW   - Vector of length nrow containing the row sums on output.
c              (Output)
c     IFLAG  - Set to 3 if there are no additional nodes to process.
c              (Output)
c     KYY    - Constant mutlipliers used in forming the hash table key.
c              (Input)
c     KEY    - Vector of length LDKEY containing the hash table keys.
c              (Input/output)
c     LDKEY  - Length of vector KEY.  (Input)
c     LAST   - Index of the last key popped off the stack.
c              (Input/output)
c     IPN    - Pointer to the linked list of past path lengths.
c              (Output)
c-----------------------------------------------------------------------
      subroutine f6xact (nrow, irow, iflag, kyy, key, ldkey, last, ipn)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, iflag, ldkey, last, ipn, irow(*), kyy(*), key(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    j, kval
c                                  SPECIFICATIONS FOR SAVE VARIABLES
c
   10 last = last + 1
      if (last .le. ldkey) then
         if (key(last) .lt. 0) go to 10
c                                  Get KVAL from the stack
         kval      = key(last)
         key(last) = -9999
         do 20  j=nrow, 2, -1
            irow(j) = kval/kyy(j)
            kval    = kval - irow(j)*kyy(j)
   20    continue
         irow(1) = kval
         ipn     = last
      else
         last  = 0
         iflag = 3
      end if
      return
      end
c-----------------------------------------------------------------------
c  Name:       F7XACT
c
c  Purpose:    Generate the new nodes for given marinal totals.
c
c  Usage:      CALL F7XACT (NROW, IMAX, IDIF, K, KS, IFLAG)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     IMAX   - The row marginal totals.  (Input)
c     IDIF   - The column counts for the new column.  (Input/output)
c     K      - Indicator for the row to decrement.  (Input/output)
c     KS     - Indicator for the row to increment.  (Input/output)
c     IFLAG  - Status indicator.  (Output)
c              If IFLAG is zero, a new table was generated.  For
c              IFLAG = 1, no additional tables could be generated.
c-----------------------------------------------------------------------
      subroutine f7xact (nrow, imax, idif, k, ks, iflag)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, k, ks, iflag, imax(*), idif(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, k1, m, mm
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  min0
      integer    min0
c
      iflag = 0
c                                  Find node which can be
c                                  incremented, ks
      if (ks .eq. 0) then
   10    ks = ks + 1
         if (idif(ks) .eq. imax(ks)) go to 10
      end if
c                                 Find node to decrement (>ks)
   20 if (idif(k).gt.0 .and. k.gt.ks) then
         idif(k) = idif(k) - 1
   30    k = k - 1
         if (imax(k) .eq. 0) go to 30
         m = k
c                                 Find node to increment (>=ks)
   40    if (idif(m) .ge. imax(m)) then
            m = m - 1
            go to 40
         end if
         idif(m) = idif(m) + 1
c                                 Change ks
         if (m .eq. ks) then
            if (idif(m) .eq. imax(m)) ks = k
         end if
      else
c                                 Check for finish
   50    do 60  k1=k + 1, nrow
            if (idif(k1) .gt. 0) go to 70
   60    continue
         iflag = 1
         go to 9000
c                                 Reallocate counts
   70    mm = 1
         do 80  i=1, k
            mm      = mm + idif(i)
            idif(i) = 0
   80    continue
         k = k1
   90    k = k - 1
         m       = min0(mm,imax(k))
         idif(k) = m
         mm      = mm - m
         if (mm.gt.0 .and. k.ne.1) go to 90
c                                 Check that all counts
c                                 reallocated
         if (mm .gt. 0) then
            if (k1 .ne. nrow) then
               k = k1
               go to 50
            end if
            iflag = 1
            go to 9000
         end if
c                                 Get ks
         idif(k1) = idif(k1) - 1
         ks       = 0
  100    ks = ks + 1
         if (ks .gt. k) go to 9000
         if (idif(ks) .ge. imax(ks)) go to 100
      end if
c
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       F8XACT
c
c  Purpose:    Routine for reducing a vector when there is a zero
c              element.
c
c  Usage:      CALL F8XACT (IROW, IS, I1, IZERO, NEW)
c
c  Arguments:
c     IROW   - Vector containing the row counts.  (Input)
c     IS     - Indicator.  (Input)
c     I1     - Indicator.  (Input)
c     IZERO  - Position of the zero.  (Input)
c     NEW    - Vector of new row counts.  (Output)
c-----------------------------------------------------------------------
      subroutine f8xact (irow, is, i1, izero, new)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    is, i1, izero, irow(*), new(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i
c
      do 10  i=1, i1 - 1
         new(i) = irow(i)
   10 continue
c
      do 20  i=i1, izero - 1
         if (is .ge. irow(i+1)) go to 30
         new(i) = irow(i+1)
   20 continue
c
      i = izero
   30 new(i) = is
   40 i = i + 1
      if (i .gt. izero) return
      new(i) = irow(i)
      go to 40
      end
c-----------------------------------------------------------------------
c  Name:       F9XACT
c
c  Purpose:    Computes the log of a multinomial coefficient.
c
c  Usage:      F9XACT(N, MM, IR, FACT)
c
c  Arguments:
c     N      - Length of IR.  (Input)
c     MM     - Number for factorial in numerator.  (Input)
c     IR     - Vector of length N containing the numebers for the
c              denominator of the factorial.  (Input)
c     FACT   - Table of log factorials.  (Input)
c     F9XACT  - The log of the multinomal coefficient.  (Output)
c-----------------------------------------------------------------------
      double precision function f9xact (n, mm, ir, fact)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    n, mm, ir(*)
      double precision fact(0:*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    k
c
      f9xact = fact(mm)
      do 10  k=1, n
         f9xact = f9xact - fact(ir(k))
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c  Name:       F10ACT
c
c  Purpose:    Computes the shortest path length for special tables.
c
c  Usage:      CALL F10ACT (NROW, IROW, NCOL, ICOL, VAL, XMIN, FACT, ND,
c                          NE, M)
c
c  Arguments:
c     NROW   - The number of rows in the table.  (Input)
c     IROW   - Vector of length NROW containing the row totals.  (Input)
c     NCOL   - The number of columns in the table.  (Input)
c     ICO    - Vector of length NCOL containing the column totals.
c              (Input)
c     VAL    - The shortest path.  (Output)
c     XMIN   - Set to true if shortest path obtained.  (Output)
c     FACT   - Vector containing the logarithms of factorials.
c              (Input)
c     ND     - Workspace vector of length NROW.
c     NE     - Workspace vector of length NCOL.
c     M      - Workspace vector of length NCOL.
c
c  Chapter:    STAT/LIBRARY Categorical and Discrete Data Analysis
c-----------------------------------------------------------------------
      subroutine f10act (nrow, irow, ncol, icol, val, xmin, fact, nd,
     &                   ne, m)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, irow(*), icol(*), nd(*), ne(*), m(*)
      double precision val, fact(0:*)
      logical    xmin
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, is, ix, nrw1
c
      do 10  i=1, nrow - 1
         nd(i) = 0
   10 continue
c
      is    = icol(1)/nrow
      ne(1) = is
      ix    = icol(1) - nrow*is
      m(1)  = ix
      if (ix .ne. 0) nd(ix) = nd(ix) + 1
c
      do 20  i=2, ncol
         ix    = icol(i)/nrow
         ne(i) = ix
         is    = is + ix
         ix    = icol(i) - nrow*ix
         m(i)  = ix
         if (ix .ne. 0) nd(ix) = nd(ix) + 1
   20 continue
c
      do 30  i=nrow - 2, 1, -1
         nd(i) = nd(i) + nd(i+1)
   30 continue
c
      ix   = 0
      nrw1 = nrow + 1
      do 40  i=nrow, 2, -1
         ix = ix + is + nd(nrw1-i) - irow(i)
         if (ix .lt. 0) return
   40 continue
c
      do 50  i=1, ncol
         ix  = ne(i)
         is  = m(i)
         val = val + is*fact(ix+1) + (nrow-is)*fact(ix)
   50 continue
      xmin = .true.
c
      return
      end
c-----------------------------------------------------------------------
c  Name:       F11ACT
c
c  Purpose:    Routine for revising row totals.
c
c  Usage:      CALL F11ACT (IROW, I1, I2, NEW)
c
c  Arguments:
c     IROW   - Vector containing the row totals.  (Input)
c     I1     - Indicator.  (Input)
c     I2     - Indicator.  (Input)
c     NEW    - Vector containing the row totals.  (Input)
c-----------------------------------------------------------------------
      subroutine f11act (irow, i1, i2, new)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    i1, i2, irow(*), new(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i
c
      do 10  i=1, i1 - 1
         new(i) = irow(i)
   10 continue
c
      do 20  i=i1, i2
         new(i) = irow(i+1)
   20 continue
c
      return
      end
c-----------------------------------------------------------------------
c  Name:       ERPRT
c
c  Purpose:    Print an error message and stop.
c
c  Usage:      CALL ERPRT (ICODE, MES)
c
c  Arguments:
c     ICODE  - Integer code for the error message.  (Input)
c     MES    - Character string containing the error message.  (Input)
c-----------------------------------------------------------------------
      subroutine prterr (icode, mes)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    icode
      character  mes*(*)
c
      write (*,*) 'FEXACT ERROR: ', icode, ' ', mes
      stop
      end
c-----------------------------------------------------------------------
c  Name:       IWORK
c
c  Purpose:    Routine for allocating workspace.
c
c  Usage:      IWORK (IWKMAX, IWKPT, NUMBER, ITYPE)
c
c  Arguments:
c     IWKMAX - Maximum length of workspace.  (Input)
c     IWKPT  - Amount of workspace currently allocated.  (Input/output)
c     NUMBER - Number of elements of workspace desired.  (Input)
c     ITYPE  - Worspace type.  (Input)
c              ITYPE  TYPE
c                2    Integer
c                3    Real
c                4    Double Precision
c     IWORK  - Index in RWRK, DWRK, or IWRK of the beginning of the
c              first element in the workspace array.  (Output)
c-----------------------------------------------------------------------
      integer function iwork (iwkmax, iwkpt, number, itype)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    iwkmax, iwkpt, number, itype
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  mod
      integer    mod
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr
c
      iwork = iwkpt
      if (itype.eq.2 .or. itype.eq.3) then
         iwkpt = iwkpt + number
      else
         if (mod(iwork,2) .ne. 0) iwork = iwork + 1
         iwkpt = iwkpt + 2*number
         iwork = iwork/2
      end if
      if (iwkpt .gt. iwkmax+1) then
         call prterr (40, 'Out of workspace.')
      end if
      return
      end
c-----------------------------------------------------------------------
c  Name:       ISORT
c
c  Purpose:    Shell sort for an integer vector.
c
c  Usage:      CALL ISORT (N, IX)
c
c  Arguments:
c     N      - Lenth of vector IX.  (Input)
c     IX     - Vector to be sorted.  (Input/output)
c-----------------------------------------------------------------------
      subroutine isort (n, ix)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    n, ix(*)
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ikey, il(10), it, iu(10), j, kl, ku, m
c                                  SPECIFICATIONS FOR SUBROUTINES
      external   prterr
c                                  Sort IX
      m = 1
      i = 1
      j = n
   10 if (i .ge. j) go to 40
      kl   = i
      ku   = j
      ikey = i
      j    = j + 1
c                                  Find element in first half
   20 i = i + 1
      if (i .lt. j) then
         if (ix(ikey) .gt. ix(i)) go to 20
      end if
c                                  Find element in second half
   30 j = j - 1
      if (ix(j) .gt. ix(ikey)) go to 30
c                                  Interchange
      if (i .lt. j) then
         it    = ix(i)
         ix(i) = ix(j)
         ix(j) = it
         go to 20
      end if
      it       = ix(ikey)
      ix(ikey) = ix(j)
      ix(j)    = it
c                                  Save upper and lower subscripts of
c                                  the array yet to be sorted
      if (m .lt. 11) then
         if (j-kl .lt. ku-j) then
            il(m) = j + 1
            iu(m) = ku
            i     = kl
            j     = j - 1
         else
            il(m) = kl
            iu(m) = j - 1
            i     = j + 1
            j     = ku
         end if
         m = m + 1
         go to 10
      else
         call prterr (20, 'This should never occur.')
      end if
c                                  Use another segment
   40 m = m - 1
      if (m .eq. 0) go to 9000
      i = il(m)
      j = iu(m)
      go to 10
c
 9000 return
      end
c-----------------------------------------------------------------------
c  Name:       GAMMDS
c
c  Purpose:    Cumulative distribution for the gamma distribution.
c
c  Usage:      PGAMMA (Q, ALPHA,IFAULT)
c
c  Arguments:
c     Q      - Value at which the distribution is desired.  (Input)
c     ALPHA  - Parameter in the gamma distribution.  (Input)
c     IFAULT - Error indicator.  (Output)
c               IFAULT  DEFINITION
c                 0     No error
c                 1     An argument is misspecified.
c                 2     A numerical error has occurred.
c     PGAMMA - The cdf for the gamma distribution with parameter alpha
c              evaluated at Q.  (Output)
c-----------------------------------------------------------------------
c
c       Algorithm AS 147 APPL. Statist. (1980) VOL. 29, P. 113
c
c       Computes the incomplete gamma integral for positive
c       parameters Y, P using and infinite series.
c
      double precision function gammds (y, p, ifault)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    ifault
      double precision y, p
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer    ifail
      double precision a, c, f
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      double precision e, one, zero
      save       e, one, zero
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog, dexp
      double precision dlog, dexp
c                                  SPECIFICATIONS FOR FUNCTIONS
      external   alogam
      double precision alogam
      double precision zexp, zlog
c
      data e, zero, one/1.0d-6, 0.0d0, 1.0d0/
c
      zexp(a) = dexp(a)
      zlog(a) = dlog(a)
c
c       Checks for the admissibility of arguments and value of F
c
      ifault = 1
      gammds = zero
      if (y.le.zero .or. p.le.zero) return
      ifault = 2
c
c       ALOGAM is natural log of gamma function
c       no need to test ifail as an error is impossible
c
      f = zexp(p*zlog(y)-alogam(p+one,ifail)-y)
      if (f .eq. zero) return
      ifault = 0
c
c       Series begins
c
      c      = one
      gammds = one
      a      = p
   10 a = a + one
      c      = c*y/a
      gammds = gammds + c
      if (c/gammds .gt. e) go to 10
      gammds = gammds*f
      return
      end
c-----------------------------------------------------------------------
c  Name:       ALOGAM
c
c  Purpose:    Value of the log-gamma function.
c
c  Usage:      ALOGAM (X, IFAULT)
c
c  Arguments:
c     X      - Value at which the log-gamma function is to be evaluated.
c              (Input)
c     IFAULT  - Error indicator.  (Output)
c               IFAULT  DEFINITION
c                 0     No error
c                 1     X .LT. 0
c     ALGAMA - The value of the log-gamma function at XX.  (Output)
c-----------------------------------------------------------------------
c
c        Algorithm ACM 291, Comm. ACM. (1966) Vol. 9, P. 684
c
c        Evaluates natural logarithm of gamma(x)
c        for X greater than zero.
c
      double precision function alogam (x, ifault)
c                                  SPECIFICATIONS FOR ARGUMENTS
      integer    ifault
      double precision x
c                                  SPECIFICATIONS FOR LOCAL VARIABLES
      double precision f, y, z
c                                  SPECIFICATIONS FOR SAVE VARIABLES
      double precision a1, a2, a3, a4, a5, half, one, seven, zero
      save       a1, a2, a3, a4, a5, half, one, seven, zero
c                                  SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog
      double precision dlog
      double precision zlog
c
c        The following constants are dlog(2PI)/2,
c        half, zero, one, seven
c
      data a1, a2, a3, a4, a5/0.918938533204673d0, 0.000595238095238d0,
     &     0.000793650793651d0, 0.002777777777778d0, 
     &     0.083333333333333d0/
      data half, zero, one, seven/0.5d0, 0.0d0, 1.0d0, 7.0d0/
c
      zlog(f) = dlog(f)
c
      alogam = zero
      ifault = 1
      if (x .lt. zero) return
      ifault = 0
      y      = x
      f      = zero
      if (y .ge. seven) go to 30
      f = y
   10 y = y + one
      if (y .ge. seven) go to 20
      f = f*y
      go to 10
   20 f = -zlog(f)
   30 z = one/(y*y)
      alogam = f + (y-half)*zlog(y) - y + a1 + (((-a2*z+a3)*z-a4)*z+a5)
     &         /y
      return
      end

