*$ CREATE MAGFLD.FOR
*COPY MAGFLD
*
*===magfld=============================================================*
*
      SUBROUTINE MAGFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1988-2010      by Alberto Fasso` & Alfredo Ferrari *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     Created  in     1988         by    Alberto Fasso`                *
*                                                                      *
*                                                                      *
*     Last change on 06-Nov-10     by    Alfredo Ferrari               *
*                                                                      *
*     Input variables:                                                 *
*            x,y,z = current position                                  *
*            nreg  = current region                                    *
*     Output variables:                                                *
*            btx,bty,btz = cosines of the magn. field vector           *
*            B = magnetic field intensity (Tesla)                      *
*            idisc = set to 1 if the particle has to be discarded      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CMEMFL)'
      INCLUDE '(CSMCRY)'
*
*
*----------------------------------------------------------------------*
*     Parameters Section *
*      
        INTEGER i, j, NR                
        PARAMETER (NR = 851)
        DOUBLE PRECISION BR(NR),RB(NR),PI,R,dr,theta,thetar,phi
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* LFIRST Section  *
* Reads the magnetic field map and saves values in respective vectors *
* RB(i) - Discrete position vector in radial direction *
* BR(i) - Magnetic field magnitude at radial position *
*
        LOGICAL LFIRST
        SAVE LFIRST
        DATA LFIRST / .TRUE. / 


        IF (LFIRST) THEN
*
*----------------------------------------------------------------------*

        OPEN(UNIT = 40, FILE = '/data/WORK/ESSnuSB/collector/fluka_run/2
     &0190320/input/LA_Bfield/RDATA.txt')
        OPEN(UNIT = 41, FILE = '/data/WORK/ESSnuSB/collector/fluka_run/2
     &0190320/input/LA_Bfield/BRDATA.txt')
*
        DO i = 1, NR
        READ(40, *) RB(i)     
        READ(41, *) BR(i)     
        END DO
        CLOSE(40 41)
        
        LFIRST = .FALSE.
        ENDIF
*----------------------------------------------------------------------*
*
        IDISC = 0
*----------------------------------------------------------------------*
        pi=4.0D0*DATAN(1.0D0)
*     
        B=0.0D0
        BTX=0.0D0
        BTY=0.0D0
        BTZ=1.0D0
*        
*     Finds closest representation of simulated particle
c$$$* in the position vector RB(i) and thereafter returns the index of the 
c$$$* found position
*
        dr = RB(NR)/(NR-1)
        R = SQRT(X**2 + Y**2)
        theta = DATAN2(Y,X)*180/PI
        thetar= DATAN2(Y,X)
        j=NINT(R/dr)+1
*
        IF(J.GE.NR) THEN
           j=NR
        ELSE IF(J.LE.0) THEN
           j=1
           R=0.0D0
        ENDIF
*
*     
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
***   As of 2019-03-08:
***   Region in the core: VoCore, region No. 2
***   Region in between Vanes: VoLA, region No. 3
***   Region between Vanes and outer shell: VoCirc, region No. 4
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
       
*     Find in which magnetic segment the particle is 
        IF(J.EQ.1) THEN
           THETA=0.0D0
           B=0.0D0
        ELSE
           B=Br(j)
           theta=theta
        ENDIF
        IF(theta.LT.0) THEN
           theta=theta+360
        ENDIF
*
        phi=0.0D0
*     If in between Vanes (region 3: VoLA),
*        choose the direction cosines of that segment
        IF(NREG.EQ.3) THEN
           phi=9
            IF(theta.GE.10.and.theta.LE.35) THEN
              phi=pi/8
           ELSE IF(theta.GE.55.and.theta.LE.80) THEN
              phi=3*pi/8
           ELSE IF(theta.GE.100.and.theta.LE.125) THEN
              phi=5*pi/8
           ELSE IF(theta.GE.145.and.theta.LE.170) THEN
              phi=7*pi/8
           ELSE IF(theta.GE.190.and.theta.LE.215) THEN
              phi=9*pi/8
           ELSE IF(theta.GE.235.and.theta.LE.260) THEN
              phi=11*pi/8
           ELSE IF(theta.GE.280.and.theta.LE.305) THEN
              phi=13*pi/8
           ELSE IF(theta.GE.325.and.theta.LE.350) THEN
              phi=15*pi/8
           ELSE
              B=0.0D0
           ENDIF
*     
        ELSE IF(NREG.EQ.4) THEN
           phi=thetar
        ELSE IF(NREG.EQ.2) THEN
           phi=thetar
        ELSE
           phi=thetar
           B=0.0D0
        ENDIF
        BTX=-sin(phi)
        BTY= cos(phi)
        BTZ=0.0D0
*        BTZ=sqrt(ONEONE-BTX**2-BTY**2)
**** This seems to crash:
*        BTZ=sqrt(1-BTX**2-BTY**2)  
**** Does it become imaginary?!?
*
c$$$        OPEN (UNIT=80, file='/home/majol640/ESSnuSB/collector/Los_Alamos
c$$$     &/simulations/mina/test_20190117/Output/test.txt')
c$$$        WRITE(80,*) theta,phi,R,j,NREG,BTX,BTY,BTZ,B
*     
*----------------------------------------------------------------------*
c$$$*  +-------------------------------------------------------------------*
c$$$*  |  Earth geomagnetic field:
c$$$      IF ( LGMFLD ) THEN
c$$$         CALL GEOFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )
c$$$         RETURN
c$$$      END IF
c$$$*  |
c$$$*  +-------------------------------------------------------------------*
*
*
*
*
      RETURN
*=== End of subroutine Magfld =========================================*
      END

