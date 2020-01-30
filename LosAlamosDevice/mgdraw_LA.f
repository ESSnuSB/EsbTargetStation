*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2013      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created on   01 March 1990   by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*     Last change   12-Nov-13      by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(QUEMGD)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
*
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
c$$$      WRITE (IODRAW) NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
c$$$     &               SNGL (WTRACK)
c$$$      WRITE (IODRAW) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
c$$$     &                 SNGL (ZTRACK (I)), I = 0, NTRACK ),
c$$$     &               ( SNGL (DTRACK (I)), I = 1, MTRACK ),
c$$$     &                 SNGL (CTRACK)
*  +-------------------------------------------------------------------*
*  |  Quenching is activated
      IF ( LQEMGD ) THEN
         IF ( MTRACK .GT. 0 ) THEN
            RULLL  = ZERZER
            CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
c$$$            WRITE (IODRAW) ( ( SNGL (DTQUEN (I,JBK)), I = 1, MTRACK ),
c$$$     &                         JBK = 1, NQEMGD )
         END IF
      END IF
*  |  End of quenching
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
*
c$$$ 100     FORMAT(I9,3X,I3,3X,F10.4,3X,F10.4,3X,F10.4,3X,F10.4,3X,F10.4,
c$$$     &        3X,F10.4,3X,F10.4)
c$$$*---------------------------------------
c$$$*    Detection plane 1, right behind LA:      
c$$$*---------------------------------------      
c$$$      IF(MREG.EQ.4.and.NEWREG.EQ.5) THEN
c$$$         OPEN (UNIT=31, file= 'DetectionPlane.txt')
c$$$         WRITE(31,100) NCASE, JTRACK, SNGL(XSCO),SNGL(YSCO),SNGL(ZSCO),
c$$$     &           SNGL(CXTRCK),SNGL(CYTRCK),SNGL(CZTRCK),SNGL(PTRACK)
c$$$       END  IF
*---------------------------------------
*    Detection plane 2, end of tunnel:
*---------------------------------------      
      IF(MREG.EQ.12.and.NEWREG.EQ.13) THEN
*-----------------------------------------------
*    NUE
*-----------------------------------------------         
         IF(JTRACK.EQ.5) THEN
            OPEN (UNIT=51, file= 'LA_end_tunnel_nue.txt')
            WRITE(51,*) NCASE,SNGL(XSCO),SNGL(YSCO),SNGL(CXTRCK),
     &           SNGL(CYTRCK),SNGL(CZTRCK),SNGL(PTRACK),SNGL(WTRACK)
*-----------------------------------------------
*    NUEBAR
*-----------------------------------------------
         ELSEIF(JTRACK.EQ.6) THEN    
            OPEN (UNIT=52, file= 'LA_end_tunnel_nuebar.txt')
            WRITE(52,*) NCASE,SNGL(XSCO),SNGL(YSCO),SNGL(CXTRCK),
     &           SNGL(CYTRCK),SNGL(CZTRCK),SNGL(PTRACK),SNGL(WTRACK)
*-----------------------------------------------
*    NUMU
*-----------------------------------------------
         ELSEIF(JTRACK.EQ.27) THEN    
            OPEN (UNIT=53, file= 'LA_end_tunnel_numu.txt')
            WRITE(53,*) NCASE,SNGL(XSCO),SNGL(YSCO),SNGL(CXTRCK),
     &           SNGL(CYTRCK),SNGL(CZTRCK),SNGL(PTRACK),SNGL(WTRACK)
*-----------------------------------------------
*    NUMUBAR
*-----------------------------------------------
         ELSEIF(JTRACK.EQ.28) THEN    
            OPEN (UNIT=54, file= 'LA_end_tunnel_numubar.txt')
            WRITE(54,*) NCASE,SNGL(XSCO),SNGL(YSCO),SNGL(CXTRCK),
     &           SNGL(CYTRCK),SNGL(CZTRCK),SNGL(PTRACK),SNGL(WTRACK)
         END IF
       END  IF
*      
*
      RETURN
*
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
*
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
c$$$      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
c$$$      WRITE (IODRAW)  SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
*  +-------------------------------------------------------------------*
*  |  Quenching is activated : calculate quenching factor
*  |  and store quenched energy in DTQUEN(1, jbk)
      IF ( LQEMGD ) THEN
         RULLL = RULL
         CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
c$$$         WRITE (IODRAW) ( SNGL (DTQUEN(1, JBK)), JBK = 1, NQEMGD )
      END IF
*  |  end quenching
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
c$$$      WRITE (IODRAW) -NCASE, NPFLKA, NSTMAX, SNGL (TKESUM),
c$$$     &                SNGL (WEIPRI)
*======================================================================*
*======================================================================*
*======================================================================*
c$$$      OPEN (UNIT=70, file='/home/majol640/ESSnuSB/collector/Los_Alamos/s
c$$$     &imulations/mina/test_20190121/Output/SourceSODRAW.txt')
c$$$      WRITE(70,*)  SNGL(XFLK(1)), SNGL(YFLK(1)), SNGL(ZFLK(1)),
c$$$     &                 SNGL(TXFLK(1)),SNGL(TYFLK(1)),SNGL(TZFLK(1)),
c$$$     &                 SNGL(PMOFLK(1)), SNGL(TKEFLK(1))
*======================================================================*
*======================================================================*
*======================================================================*
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope: it works only for 1 source particle on
*  |  the stack for the time being
      IF ( ILOFLK (NPFLKA) .GE. 100000 .AND. LRADDC (NPFLKA) ) THEN
         IARES  = MOD ( ILOFLK (NPFLKA), 100000  )  / 100
         IZRES  = MOD ( ILOFLK (NPFLKA), 10000000 ) / 100000
         IISRES = ILOFLK (NPFLKA) / 10000000
         IONID  = ILOFLK (NPFLKA)
c$$$         WRITE (IODRAW) ( IONID,SNGL(-TKEFLK(I)),
c$$$     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c$$$     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c$$$     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c$$$     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
      ELSE IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
         IONID = ILOFLK (NPFLKA)
         CALL DCDION ( IONID )
c$$$         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-IONID)),
c$$$     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c$$$     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c$$$     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c$$$     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: ???
      ELSE IF ( ILOFLK (NPFLKA) .LT. -6 ) THEN
c$$$         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-ILOFLK(NPFLKA))),
c$$$     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c$$$     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c$$$     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c$$$     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
c$$$         WRITE (IODRAW) ( ILOFLK(I), SNGL (TKEFLK(I)+AM(ILOFLK(I))),
c$$$     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
c$$$     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
c$$$     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
c$$$     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*             110: decay products                                      *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*             237: mu pair     production secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=npflka          *
*                                                                      *
*======================================================================*
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
* No output by default:
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

