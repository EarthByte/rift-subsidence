#!/usr/bin/env python
# coding: utf-8
#
#**********************************************************************
#
#             CALCULATING THEORETICAL SUBSIDENCE CURVES
#
#                     NICKY WHITE OCTOBER 1986
#                                 AUGUST 1989
#
#                     Paul Bellingham Jan 1996
#
#***********************************************************************
#
#   Modified 5/4/90 to have realistic later stretching events
#   at large BETAs.
#
#   PTB 22/1/96 Modified so that you can have a variable crustal
#   thickness and then program calculates the lithospheric
#   thickness isostatically balanced against a MOR. This is now
#   finished.
#  
#   RDM 06/2004 removed all plotting functions and compile for Linux/OSX
#   with g77.  Also fixed various mistakes in how input is read in.
#
#   02/2020 michael.chin@sydney.edu.au rewrote in Python.
#
#   Output:  4 column file with Time (my) after initial rifting,
#   tectonic water-loaded basement subsidence (km),
#   heat flow (mW/m2) and strain rate (per billion year)
#
#   Syn-rift subsidence is calculated based on given time period of
#   rifting and constant beta value.  For each selected beta value
#   a 4-column set of output data will be written to the file
#   Each selected beta value will result.
#

#You need python3 to run this program!

#This file is translated from gforisos.f, which is the reason of fortran style programming.

import math

#
# read in a floating point number from stdin
#
def get_float_from_stdin(msg):
    while(1):
        tmp = input(msg)
        #print(tmp)
        if(len(tmp)==0):
            return None
        
        try:
            tmp=float(tmp)
        except:
            print('Please type in a floating point number.')
            continue
        return tmp

#
# read in an integer number from stdin
#
def get_int_from_stdin(msg,values=[]):
    while(1):
        tmp = input(msg)
        #print(tmp)
        if(len(tmp)==0):
            return None
        
        try:
            tmp=int(tmp)
        except:
            print('Please type in an interger')
            continue
        
        if values:
            if tmp in values:
                return tmp
            else:
                continue
        else:
            return tmp

#
# read in "yes or no" response from user
#
def get_yes_or_no_from_stdin(msg, no_default=False):
    while(1):
        tmp = input(msg)
        #print(tmp)
        if(len(tmp)==0 and not no_default):
            return None
        
        if tmp.lower() in ['y','yes']: return 'y'
        elif tmp.lower() in ['n','no']: return 'n'
        else: 
            print('Reply y/yes or n/no')
            continue
        
        
#PARAMETER
NDIM1=4
NDIM2=10000
RAD=0.0174533
PI=3.1415927
RHOM=3.35
RHOC=2.78
RHOW=1.03
TASTH=1333.0
DIFFUS=0.804E-6
TEXPAN=3.28E-5
TAU=62.8
TAUS=43.28
IHFLW=1
CONDTY=38.544
TCO=7.E+3
DW=2.5E+3

#INITVR
BETMIN=1.5
BETMAX=2.0
BETINC=0.5
BETA2=1.1
BETA3=1.1
ILOAD=1
ICOMP=0
TBEG1=160.
TEND1=100.
TBEG2=60.
TEND2=50.
TBEG3=20.
TEND3=10.
TSTOP=0.
ITYPE1=1
ITYPE2=1
ITYPE3=1
IRIFT=1
RHOS=2.6
PHI0=0.6
CONST=2.0
CYCLES = 1.5
FRAC = 0.8
TL = 125000
TC = 30000

# Nondimensionalised parametres
TIMSC = DTBEG1 = DTEND1 = DTBEG2 = DTEND2 = DTBEG3 = DTEND3 = DTSTOP = DTST = DTAUS = 0.0

# constants for subsidence calculations
RHOA = CONBET = CONTRM = CONB1 = CONB2 = GP1 = GP2 = GP3 = 0.0

#ITYPE = 4 parameters ??
AMPL = WAVEN = 1.

SINIT = TRINIT = 0.

K=1 # fixme: this is a counter. should not be here as global. 

DATA=[ [ None for y in range( NDIM2 ) ] 
         for x in range( NDIM1 ) ] 
TEM=[ [ None for y in range( NDIM2 ) ] 
         for x in range( NDIM1 ) ] 

PLABEL = 'my_label'

def calculate_TL():
    #
    #
    #   This bit takes input of a continental crustal
    #   thickness and gives a lithospheric thickness that
    #   is balanced isostatically against a MOR
    #
    #   This assumes that the system starts at sea level.
    #

    #   Theory is in Bown 1993 thesis p.125 and is checked
    #
    global TL
    AINT = ((RHOM*TEXPAN*TASTH)/2)
    BINT = ((RHOW*DW + RHOC*TCO) - (DW+TCO)*RHOM*(1-TEXPAN*TASTH) - TC*(RHOC-RHOM))
    CINT = (((TC**2)*TASTH*TEXPAN)*(RHOM-RHOC))/2
    #
    #  SOLVING THE QUADRATIC
    #
    TL = (-BINT + (math.sqrt((BINT**2)-(4.*AINT*CINT))))/(2.*AINT)
    if TL <= 20.E+3 or TL>=200.E+3:
        TL = (-BINT - (math.sqrt((BINT**2)-(4.*AINT*CINT))))/(2.*AINT)

    TL = -TL
    if TL<=TC: 
        print('ERROR: The TL is not supposed to be less than TC! Check the equation and numbers.')

    if TL >= TC:
        TLM = TL/1000
    print(f' Lithospheric thickness is {TLM:.6f} kilometres')

#    
# Collect parameters from stdin (keyboard, screen)    
#    
def collect_parameters_from_stdin():
    print('                ********************************************')
    print('                *                                          *')
    print('                *    Calculating subsidence curves as a    *')
    print('                *     function of time and stretching      *')
    print('                *                                          *')
    print('                ********************************************')

    print('Input parameters required to run SUBSID.')
    print('This program is about to collect parameters on screen.')
    print('Alternatively, you can save the parameters in a json file and')
    print('run the program with commmand "python3 RiftSubsidence.py parameters.json".')
    print('An example parameters.json has been prepared for your convenience.')
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    print('Stretching periods:')
    print(' ~~~~~~~~~~~~~~~~~~')
    
    global BETMIN,BETMAX,BETINC,BETA2,BETA3,ILOAD,ICOMP,TBEG1,TEND1,TBEG2,TEND2,TBEG3,TEND3,TSTOP,ITYPE1
    global ITYPE2,ITYPE3,IRIFT,RHOS,PHI0,CONST,CYCLES,FRAC,TL,TC

    global PLABEL 
    tmp = input(f'Label ({PLABEL}): ')
    if tmp: PLABEL = tmp
    
    tmp = get_int_from_stdin(f' 1, 2, or 3 periods of rifting ({IRIFT}): ', [1,2,3])
    if tmp: IRIFT = tmp

    tmp = get_float_from_stdin(f'Initial Cont. crustal thickness (m)({TC:.0f}): ')
    if tmp: TC = tmp

    tmp = get_float_from_stdin(f' start of first rifting period ({TBEG1:.1f}): ')
    if tmp: TBEG1 = tmp

    tmp = get_float_from_stdin(f' end of first rifting period ({TEND1:.1f}): ')
    if tmp: TEND1 = tmp

    tmp = get_float_from_stdin(f' minimum value of BETA ({BETMIN:.1f}): ')
    if tmp: BETMIN = tmp

    tmp = get_float_from_stdin(f' maximum value of BETA ({BETMAX:.1f}): ')
    if tmp: BETMAX = tmp

    tmp = get_float_from_stdin(f' increment in BETA ({BETINC:.1f}): ')
    if tmp: BETINC = tmp

    tmp = get_int_from_stdin(f' strain rate for first rift period ({ITYPE1}): ')
    if tmp: ITYPE1 = tmp

    if ITYPE1 == 4:
        tmp = get_float_from_stdin(f'  CYCLES ({CYCLES:.1f}): ')
        if tmp: CYCLES = tmp
        tmp = get_float_from_stdin(f'  FRAC ({FRAC:.2f}): ')
        if tmp: FRAC = tmp

    if IRIFT >= 2:
        tmp = get_float_from_stdin(f' start of second rifting period ({TBEG2:.1f}): ')
        if tmp: TBEG2 = tmp

        tmp = get_float_from_stdin(f' end of second rifting period ({TEND2:.1f}): ')
        if tmp: TEND2 = tmp

        tmp = get_float_from_stdin(f' BETA for second rifting phase ({BETA2:.1f}): ')
        if tmp: BETA2 = tmp

        #tmp = get_int_from_stdin(f' strain rate for second rift period ({ITYPE2}): ')
        #if tmp: ITYPE2 = tmp

    if IRIFT >= 3:
        tmp = get_float_from_stdin(f' start of third rifting period ({TBEG3:.1f}): ')
        if tmp: TBEG3 = tmp

        tmp = get_float_from_stdin(f' end of third rifting period ({TEND3:.1f}): ')
        if tmp: TEND3 = tmp

        tmp = get_float_from_stdin(f' BETA for third rifting phase ({BETA3:.1f}): ')
        if tmp: BETA3 = tmp

        #tmp = get_int_from_stdin(f' strain rate for third rift period ({ITYPE3}): ')
        #if tmp: ITYPE3 = tmp

    tmp = get_float_from_stdin(f'  final time ({TSTOP:.1f}): ')
    if tmp: TSTOP = tmp

    tmp = get_int_from_stdin(f' loading (0=air, 1=water, 2=sediments) ({ILOAD}): ', [0,1,2])
    if tmp:  ILOAD = tmp

    if ILOAD == 2:
        xxx='Y' if ICOMP else 'N'
        tmp = get_yes_or_no_from_stdin(f' Is compaction required? ({xxx})')
        if tmp:
            ICOMP = 0 if tmp=='n' else 1

        if not ICOMP:
            tmp = get_float_from_stdin(f' Sediment density ({RHOS:.1f}): ')
            if tmp: RHOS = tmp
        else:
            tmp = get_float_from_stdin(f' Solid sediment density ({RHOS:.1f}): ')
            if tmp: RHOS = tmp

            while(1):
                tmp = get_float_from_stdin(f' Initial porosity ({PHI0*100}%): ')
                if tmp: PHI0 = tmp/100.
                if PHI0 < 0 or PHI0 > 0.9:
                    print(' ** Reply in the range 0 to 90')
                    continue
                else:
                    break

            tmp = get_float_from_stdin(f' Depth constant, LAMBDA ({CONST}): ')
            if tmp: CONST = tmp

    print('Chosen parameters:')
    print(' ~~~~~~~~~~~~~~~~~')
    print(f'number of rift periods = {IRIFT},')
    print(f'min. BETA = {BETMIN}, max. BETA = {BETMAX}, Increment in BETA = {BETINC}')
    print(f'lithos. thick. = {TL:.0f} m,  cont. crust = {TC} m')
    print(f'start of first rift = {TBEG1}, end of first rift = {TEND1}')

    #print('Strain rate = {ITYPE1}')
    if IRIFT >=2:
        print('BETA2 = {BETA2}')
        print('start second rift = {TBEG2}, end second rift = {TEND2}')
        #print('Strain rate = {ITYPE2}')

    if IRIFT >=3:
        print('BETA3 = {BETA3}')
        print('start third rift = {TBEG3}, end third rift = {TEND3}')
        #print('Strain rate = {ITYPE2}')

    print(f'Final time = {TSTOP}')
    print(f'Type of loading = {ILOAD}')
    if ILOAD == 2 and ICOMP == 0:
        print(f'Sediment density = {RHOS}')
    elif ILOAD == 2 and ICOMP == 1:
        print(f'Solid sediment density = {RHOS}')
        print(f'Initial porosity = {PHI0}, Depth constant = {CONST}')

    
#
#  Make all parameters dimensionless.
#
def nondimensionalise_parameters():
    global TIMSC, DTBEG1, DTEND1, DTBEG2, DTEND2, DTBEG3, DTEND3, DTSTOP, DTST, DTAUS
          
    TIMSC = (TL**2)/(DIFFUS*3.15E+13)
    if IRIFT >= 1:
        DTBEG1 = (-TBEG1+TBEG1)/TIMSC
        DTEND1 = (-TEND1+TBEG1)/TIMSC
    
    if IRIFT >= 2:    
        DTBEG2 = (-TBEG2+TBEG1)/TIMSC
        DTEND2 = (-TEND2+TBEG1)/TIMSC
        
    if IRIFT >= 3:
        DTBEG3 = (-TBEG3+TBEG1)/TIMSC
        DTEND3 = (-TEND3+TBEG1)/TIMSC
        
    DTSTOP = (-TSTOP+TBEG1)/TIMSC
    DTST = 1./TIMSC
    DTAUS = TAUS/TIMSC

    
#    
# Set up constants for subsidence calculations  
#
def setup_constants():
    global RHOA, CONBET, CONTRM, CONB1, CONB2
    RHOA = RHOM*(1.0 - (TEXPAN*TASTH))
    if ICOMP == 0:
        if ILOAD == 2:
            CONBET = ((RHOM-RHOC)*TC)/(RHOA-RHOS)
            CONTRM = (TEXPAN*RHOM*TL*TASTH)/(RHOA-RHOS)
        elif ILOAD == 1:
            CONBET = ((RHOM-RHOC)*TC)/(RHOA-RHOW)
            CONTRM = (TEXPAN*RHOM*TL*TASTH)/(RHOA-RHOW)
        elif ILOAD == 0:
            CONBET = ((RHOM-RHOC)*TC)/(RHOA)
            CONTRM = (TEXPAN*RHOM*TL*TASTH)/(RHOA)
       
    else:
        CONBET = ((RHOM-RHOC)*TC)/(RHOA)
        CONTRM = (TEXPAN*RHOM*TL*TASTH)/(RHOA)
     
    CONB1 = (RHOA-RHOW)/(RHOA-RHOS)
    CONB2 = (RHOS-RHOW)/(RHOA-RHOS)

    #
    #   Set GP, dless magnitude of the
    #   vertical velocity gradient across
    #   a depth TL for a second and third
    #   rift period, as appropriate.
    #
    global GP2, GP3 
    if IRIFT >= 2:
        if ITYPE2 == 1:
            GP2 = math.log(BETA2)/(DTEND2-DTBEG2)
        elif ITYPE2 == 2:
            GP2 = math.log(BETA2)/(1.-math.exp(-(DTEND2-DTBEG2)/DTAUS))
        elif ITYPE2 == 3:
            GP2 = math.log(BETA2)/(math.exp((DTEND2-DTBEG2)/DTAUS)-1.)
       
    elif IRIFT >= 3:
        if ITYPE3 == 1:
            GP3 = math.log(BETA3)/(DTEND3-DTBEG3)
        elif ITYPE3 == 2:
            GP3 = math.log(BETA3)/(1.-math.exp(-(DTEND3-DTBEG3)/DTAUS))
        elif ITYPE3 == 3:
            GP3 = math.log(BETA3)/(math.exp((DTEND3-DTBEG3)/DTAUS)-1.)
 

#
#   Set dless constants for sinusoidally varying strain rate
#
def set_dless_constants(BETA1):
    print(f'\nBeta1 = {BETA1} .... ')    
    #
    #   Set GP1, dless magnitude of the
    #   vertical velocity gradient across
    #   a depth TL for main rift period.
    #
    global GP1, AMPL, WAVEN
    
    if ITYPE1 == 1:
        GP1 = math.log(BETA1)/DTEND1
    elif ITYPE1 == 2:
        GP1 = math.log(BETA1)/(1.-math.exp(-DTEND1/DTAUS))
    elif ITYPE1 == 3:
        GP1 = math.log(BETA1)/(math.exp(DTEND1/DTAUS)-1.)
    elif ITYPE1 == 4:
        GP1 = math.log(BETA1)/DTEND1
        CG = GP1
        AMPL = FRAC*math.log(BETA1)/DTEND1
        
        while(1):
            WAVEL = (DTEND1)/(CYCLES + (GP1/(PI*AMPL)))
            WAVEN = (2.*PI)/WAVEL
            ALOB = (2.*PI*DTEND1)/WAVEL
            BLOB = (AMPL*WAVEL)/(2.*PI*DTEND1)
            CLOB = WAVEL/(2.*PI*AMPL*DTEND1)
            GP2 = CG  + BLOB*(math.cos(ALOB-((2.*GP1)/AMPL)) - 1.) + ((GP1**2)*CLOB)
            ERROR = math.fabs(GP2-GP1)
            print(f' GP1={GP1}, GP2={GP2} ')
            GP1 = GP2
            if ERROR > 1.E-6:
                continue
            else:
                break
                
        print(f' ALOB = {ALOB}, WAVEL = {WAVEL}, AMPL = {AMPL}')
        print(f' DTEND1 = {DTEND1}, TIMSC = {TIMSC}')
        A = AMPL/WAVEN
        B = WAVEN*DTEND1
        C = (2.*GP1)/AMPL
        D = (GP1**2)/(2.*WAVEN*AMPL)
        BTEST = math.exp((GP1*DTEND1) + (A*(1. - math.cos(B-C))) + D)
        print(f' GP1 = {GP1}, BETA = {BTEST} ')

#
#
#
def TIMSTR(G,DTEND,DTSTOP,DTAUS,CONBET,CONTRM,
        DTBEGX,DTBEGY,DTST,ITYPE,IHFLW,KRR,IRIFT,
        BETA1,CONB1,CONB2):
#
#   calls SUBCAL
#
#   solves the dimensionless equation governing the one
#   dimensional advection and diffusion of heat in a
#   plane layer using a finite difference scheme.
#   If KRR=0, solves for one rifting episode
#   (IRIFT=1) orthe first rifting episode(IRIFT=2).
#   If KRR=1, solves for the second rifting episode.
#
#
# DTST IS THE DIMENSIONLESS TIME INTERVAL GOVERNING THE SAVE FREQUENCY
# FOR TIMSAV AND HFLSAV.
# DTEND DLESS TIME AT END OF STRETCHING
# DTSTOP DLESS TIME AT END OF CALCULATION
#

    global K
    
    print(f'\nRift Period {KRR+1} ....')
    print('Calculating temperature structure....\n')
    
    JSUBS = 0
    DZ=0.05
    JN=20
    JNP1=JN+1
    
    STOT  = 0. #???

    # CONDUCTIVE STABILITY LIMIT ON DT
    DTC=DZ**2/2.

    # ADVECTIVE STABILITY LIMIT ON DT
    # Make sure that biggest G is used (i.e. one near beginning).
    GMAX = AMPL + G
    DTCFL = DZ/GMAX
    DT = min(DTC,DTCFL)
    DTI = min(DT,DTST)

    # INITIALISE PARAMETERS AND TEMPERATURE TO CONDUCTIVE PROFILE
    #TEM(I,1) AND TEM(I,JNP1) ARE NOT CHANGED BY THE CALCULATION
    if KRR == 0: T = 0.
    if KRR > 0:  T = DTBEGX
    DT = DTI

    #DT = DT #?
    JGI=1
    JGF=2
    JFLAG=0
    IADV=1
    if KRR == 0:
        DATA[0][K-1] = T
        DATA[1][K-1] = 0.0
    
    if IHFLW == 1:
        if KRR == 0:
            DATA[2][K-1] = 1.
            if ITYPE == 1:
                DATA[3][K-1] = G
            elif ITYPE == 2:
                DATA[3][K-1] = G*math.exp(-T/DTAUS)/DTAUS
            elif ITYPE == 3:
                DATA[3][K-1] = G*math.exp(T/DTAUS)/DTAUS
            elif ITYPE == 4:
                #DATA[3][K-1] = AMPL*math.sin(WAVEN*T) + G 
                if T < G/(AMPL*WAVEN):
                    DATA[3][K-1] =  T*AMPL*WAVEN
                elif T >= G/(AMPL*WAVEN) and T < DTEND-(G/(AMPL*WAVEN)):
                    DATA[3][K-1] = AMPL*math.sin(WAVEN*(T-(G/(WAVEN*AMPL)))) + G
                elif T >= DTEND-(G/(AMPL*WAVEN)):
                    DATA[3][K-1] = (DTEND-T)*WAVEN*AMPL
            
    if KRR == 0:
        for J in range(1,JNP1+1):
            TEM[JGI-1][J-1]=1.-DZ*(J-1)
            TEM[JGF-1][J-1]=1.-DZ*(J-1)

#
# ADVANCE TIME BY DT AND CALCULATE TEMPERATURE
#
    while(1):
        for J in range(2, JN+1):
            DELSQ = (TEM[JGI-1][J]-2.*TEM[JGI-1][J-1]+TEM[JGI-1][J-2])/DZ**2
            ADV=0.
            if IADV == 1 and ITYPE == 1:
                ADV=(G*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)
            elif IADV == 1 and ITYPE == 2:
                if KRR == 0:
                    ADV=(((G*math.exp(-T/DTAUS))/DTAUS)*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)
                else:
                    ADV=(((G*math.exp(-(T-DTBEGX)/DTAUS))/DTAUS)*(1.-(J-1)*DZ)*(TEM[JGI][J]-TEM[JGI][J-2]))/(2.*DZ)
            elif IADV == 1 and ITYPE == 3:
                if KRR == 0:
                    ADV=(((G*math.exp(T/DTAUS))/DTAUS)*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)
                else:
                    ADV=(((G*math.exp((T-DTBEGX)/DTAUS))/DTAUS)*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)
            elif IADV ==1 and ITYPE == 4:
                if KRR == 0:
                    if T <= G/(AMPL*WAVEN):
                        GVECT =  T*AMPL*WAVEN
                    elif T > G/(AMPL*WAVEN) and T <= DTEND-(G/(AMPL*WAVEN)):
                        GVECT = AMPL*math.sin(WAVEN*(T-(G/(WAVEN*AMPL)))) + G
                    elif T > DTEND-(G/(AMPL*WAVEN)):
                        GVECT = (DTEND-T)*WAVEN*AMPL
                    ADV=((GVECT)*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)
                else:
                    if T <= G/(AMPL*WAVEN):
                        GVECT =  (T-DTBEGX)*AMPL*WAVEN
                    elif T > G/(AMPL*WAVEN) and T <= DTEND-(G/(AMPL*WAVEN)):
                        GVECT = AMPL*math.sin(WAVEN*(T-DTBEGX-(G/(WAVEN*AMPL)))) + G
                    elif T > DTEND-(G/(AMPL*WAVEN)):
                        GVECT = (DTEND-(T-DTBEGX))*WAVEN*AMPL
                    ADV=((GVECT)*(1.-(J-1)*DZ)*(TEM[JGI-1][J]-TEM[JGI-1][J-2]))/(2.*DZ)

            TEM[JGF-1][J-1]=TEM[JGI-1][J-1]+DT*(DELSQ-ADV)

        T=T+DT

    #
    # RESTORE VALUE OF DT
    #
        DT=DTI
        if JGI == 1:
            JGI=2
            JGF=1
        elif JGI == 2:
            JGI=1
            JGF=2

    #
    #   Use dless temp to calculate subsidence.
    #
        DTHALT = DTSTOP
        if KRR == 0 and IRIFT == 2: DTHALT = DTBEGX
        if KRR == 0 and IRIFT == 3: DTHALT = DTBEGX
        if KRR == 1 and IRIFT == 3: DTHALT = DTBEGY
        if(T >= DATA[0][K-1]+DTST) or (T == DTHALT) or (T == DTEND):
            K=K+1
            DATA[0][K-1] = T
            SUBCAL(T,IADV,CONBET,CONTRM,CONB1,CONB2,JNP1,DTEND,G,K,
                JGI,DZ,DTAUS,ITYPE,KRR,IRIFT,DTBEGX,BETA1,JSUBS)

            #
            #   Calculate heat flux as appropriate..
            #
            if IHFLW == 1:
                DATA[2][K-1] = TEM[JGI-1][JN-1]/DZ
            #
            #   Set strain rate
            #
                if IADV == 0:
                    DATA[3][K-1] = 0.
                elif IADV == 1:
                    if ITYPE == 1:
                        DATA[3][K-1] = G
                    elif ITYPE == 2:
                        DATA[3][K-1] = G*math.exp(-T/DTAUS)/DTAUS
                    elif ITYPE == 3:
                        DATA[3][K-1] = G*math.exp(T/DTAUS)/DTAUS
                    elif ITYPE == 4:
                        #DATA(4,K) = AMPL*SIN(WAVEN*T) + G
                        if T <= G/(AMPL*WAVEN):
                            DATA[3][K-1] =  T*AMPL*WAVEN
                        elif T > G/(AMPL*WAVEN) and T <= DTEND-(G/(AMPL*WAVEN)):
                            DATA[3][K-1] = AMPL*math.sin(WAVEN*(T-(G/(WAVEN*AMPL)))) + G
                        elif T > DTEND-(G/(AMPL*WAVEN)):
                            DATA[3][K-1] = (DTEND-T)*WAVEN*AMPL

        #
        # WILL NEXT STEP GO BEYOND DTE WITH IADV=1? IF NOT CONTINUE
        # IF SO REDUCE DT
        #       
        if (T+DT) <= DTEND:
            continue
        elif (T+DT) > DTEND and JFLAG == 0:
            JFLAG=1
            DT = DTEND-T
            continue

        #
        # NEXT STEP PURELY CONDUCTIVE
        #
        elif (T+DT) < DTHALT and JFLAG == 1:
            IADV=0
            continue
        #
        # NEXT STEP LAST STEP IF DT.NE.0.
        #
        elif (T+DT) >= DTHALT and JFLAG != 2:
            JFLAG=2
            DT=DTHALT-T
            if DT != 0: continue
        elif JFLAG == 2:
            #
            #   Remember  total subsidence  if 2 or 3 stage rifting.
            #
            if IRIFT > 1: STOT = DATA[1][K-1]
        break #break the while(1) loop

#
#
#
def SUBCAL(T,IADV,CONBET,CONTRM,CONB1,CONB2,JNP1,DTEND,G,K,
        JGI,DZ,DTAUS,ITYPE,KRR,IRIFT,DTBEGX,BETA1,JSUBS):
#
#   Calls no subroutines
#
#   Calculates the subsidence at time T by
#   integrating dimensionless temperature
#   over dless depth using the trapezoidal rule.
#   During the stretching interval (IADV=1)
#   a crustal thinning factor must be accounted for.
#   After this period (IADV=0) then subsidence
#   is solely due to thermal contraction.
#   If KRR=1 (second stage rift), allow
#   for fact of previously thinned crust.
#

    global TRINIT, SINIT
    R = 0.
    P = 0.
    if IADV == 1:
        #
        #   If KRR=0 then
        #   Calculate the stretching factor after a time T, BETAT
        #   If KRR=1,2 then
        #   calculate the stretching factor after a time (T-DTBEGX)
        #
        if KRR == 0: WWW = T
        if KRR > 0: WWW = T-DTBEGX
        if ITYPE == 1: 
            BETAT = math.exp(G*WWW)
        elif ITYPE == 2: 
            BETAT = math.exp(G*(1.-math.exp(-WWW/DTAUS)))
        elif ITYPE == 3:
            BETAT = math.exp(G*(math.exp(WWW/DTAUS) - 1.))
        elif ITYPE == 4:
        #BETAT = EXP((G*WWW) - (AMPL/WAVEN)*(COS(WAVEN*WWW)-1.))
            if WWW <= G/(AMPL*WAVEN):
                BETAT =  math.exp(0.5*WAVEN*AMPL*(WWW**2))
            elif WWW > G/(AMPL*WAVEN) and WWW <= (DTEND-(G/(AMPL*WAVEN))):
                if JSUBS == 0:
                    P = G/(AMPL*WAVEN)
                    JSUBS = 1

                Q = WWW - P
                BETAT = math.exp((0.5*WAVEN*AMPL*(P**2)) + (G*Q) - ((AMPL/WAVEN)*(math.cos(WAVEN*(Q)) - 1.)))
            elif WWW > (DTEND-(G/(AMPL*WAVEN))):
                if JSUBS == 1:
                    R = DTEND-(G/(AMPL*WAVEN))
                    JSUBS = 2

                S = WWW - R
                BETAT = math.exp((0.5*WAVEN*AMPL*(P**2)) + (G*(R-P)) - ((AMPL/WAVEN)*(math.cos(WAVEN*(R-P)) - 1.)) + 
			(WAVEN*AMPL*(((G*S)/(WAVEN*AMPL)) - (0.5*(S**2)))))

        TRAP = 0.0
        for I in range(2,JNP1+1):
            S = 0.5*DZ
            TRAP = S*(TEM[JGI-1][I-2] + TEM[JGI-1][I-1]) + TRAP

        if T == DTEND: TRINIT = TRAP
        if KRR == 0:
            SUBS = CONBET*(1.-(1./BETAT)) + CONTRM*(0.5 - TRAP)
            #print(f'CONBET = {round(CONBET,4):<12} \t CONTRM = {round(CONTRM,4):<12}')
            print(f'T = {round(T,4):<12} \t Q = {round(0.5-TRAP,4):<12} \t BETA = {round(BETAT,4):<12}')
        else:
            SUBS = CONBET*(2.-(1./BETAT)-(1./BETA1)) + CONTRM*(0.5 - TRAP)
            
        DATA[1][K-1] = SUBS*1.E-3
        if T == DTEND: SINIT = DATA[1][K-1]
    else: #IADV != 1
        TRAP = 0.0
        for I in range(2, JNP1+1):
            S = 0.5*DZ
            TRAP = S*(TEM[JGI-1][I-2] + TEM[JGI-1][I-1]) + TRAP

        SUBS = CONTRM*(TRINIT - TRAP)
        DATA[1][K-1] = SUBS*1.E-3 + SINIT
   
#
#
#
def COMPAC(NPTS):
#
#   Loads curves with compacted sediment
#   Calls no subroutines
#

    F = [None]*30
    ALOB = RHOA/(RHOA-RHOS)
    DLOB = RHOA - RHOS
    ELOB = PHI0*CONST*(RHOS-RHOW)
    FLOB = PHI0*(RHOS-RHOW)
    
    for I in range(1, NPTS+1):
        if DATA[0][I-1] is None: continue
        
        J=1
        F[0] = DATA[1][I-1]*ALOB
        #
        #   Newton-Raphson iteration
        #
        while(1):
            X1 = F[J-1]*DLOB - DATA[1][I-1]*RHOA + ELOB*(1.-math.exp(-F[J-1]/CONST))
            X2 = DLOB + FLOB*math.exp(-F[J-1]/CONST)
            F[J] = F[J-1] - (X1/X2)

            ERROR = math.fabs(F[J]-F[J-1])
            #
            #   Accurate to 1 cm
            #
            if ERROR > 1.E-5:
                J=J+1
                if J>= 30:
                    break
                else:
                    continue
            break
        DATA[1][I-1] = F[J]


#
#----------------------------------------------------------------------
#      
def WRITE(NPTS):
    filename = PLABEL.replace(" ", "_")+'.dat'
    with open(filename, 'w+') as file: 
        #file.write(PLABEL+'\n')
        #file.write(str(NPTS)+'\n')
        for I in range(1, NPTS+1):
            if DATA[0][I-1] is None: DATA[0][I-1] = -999999.
            if DATA[1][I-1] is None: DATA[1][I-1] = -999999.
            if DATA[2][I-1] is None: DATA[2][I-1] = -999999.
            if DATA[3][I-1] is None: DATA[3][I-1] = -999999.
            print(f'{round(DATA[0][I-1],4):<12}\t {round(DATA[1][I-1],4):<12}\t'+ 
                  f'{round(DATA[2][I-1],4):<12}\t {round(DATA[3][I-1],4):<12}')
            file.write(f'{round(DATA[0][I-1],4):<12}\t\t'+
                       f'{round(DATA[1][I-1],4):<12}\t\t'+
                       f'{round(DATA[2][I-1],4):<12}\t\t'+
                       f'{round(DATA[3][I-1],4):<12}\n')
    print(f'The data has been saved in {filename}.')
    print(f'Data length: {NPTS}')

    
def read_cfg_file(cfg_file):
    import json
    
    global BETMIN,BETMAX,BETINC,BETA2,BETA3,ILOAD,ICOMP,TBEG1,TEND1,TBEG2,TEND2,TBEG3,TEND3,TSTOP,ITYPE1
    global ITYPE2,ITYPE3,IRIFT,RHOS,PHI0,CONST,CYCLES,FRAC,TL,TC,PLABEL 
    
    with open(cfg_file, 'r') as f:
        cfg = json.load(f)
        PLABEL = cfg["PLABEL"]
        IRIFT = cfg["IRIFT"] 
        TC = cfg["TC"]
        TBEG1 = cfg["TBEG1"]
        TEND1 = cfg["TEND1"]
        BETMIN = cfg["BETMIN"]
        BETMAX = cfg["BETMAX"]
        BETINC = cfg["BETINC"]
        ITYPE1 = cfg["ITYPE1"]
        CYCLES = cfg["CYCLES"]
        FRAC = cfg["FRAC"]
        TBEG2 = cfg["TBEG2"]
        TEND2 = cfg["TEND2"]
        BETA2 = cfg["BETA2"]
        ITYPE2 = cfg["ITYPE2"]
        TBEG3 = cfg["TBEG3"]
        TEND3 = cfg["TEND3"]
        BETA3 = cfg["BETA3"]
        ITYPE3 = cfg["ITYPE3"]
        TSTOP = cfg["TSTOP"]
        ILOAD = cfg["ILOAD"]
        ICOMP = cfg["ICOMP"]
        RHOS = cfg["RHOS"]
        PHI0 = cfg["PHI0"]
        CONST = cfg["CONST"]

#save the parameters for reference
def write_out_parameters():
    import time
    timestr = time.strftime("%Y%m%d_%H%M%S")
    with open(PLABEL.replace(" ", "_")+'_parameters_'+timestr+'.json', 'w+') as file: 
        file.write("{\n")
        file.write(f"\t\"PLABEL\" : \"{PLABEL}\",\n")
        file.write(f"\t\"IRIFT\" : {IRIFT},\n")
        file.write(f"\t\"TC\" : {TC},\n")
        file.write(f"\t\"TL\" : {TL},\n")
        file.write(f"\t\"TBEG1\" : {TBEG1},\n")
        file.write(f"\t\"TEND1\" : {TEND1},\n")
        file.write(f"\t\"BETMIN\" : {BETMIN},\n")
        file.write(f"\t\"BETMAX\" : {BETMAX},\n")
        file.write(f"\t\"BETINC\" : {BETINC},\n")
        file.write(f"\t\"ITYPE1\" : {ITYPE1},\n")
        file.write(f"\t\"CYCLES\" : {CYCLES},\n")
        file.write(f"\t\"FRAC\" : {FRAC},\n")
        file.write(f"\t\"TBEG2\" : {TBEG2},\n")
        file.write(f"\t\"TEND2\" : {TEND2},\n")
        file.write(f"\t\"BETA2\" : {BETA2},\n")
        file.write(f"\t\"ITYPE2\" : {ITYPE2},\n")
        file.write(f"\t\"TBEG3\" : {TBEG3},\n")
        file.write(f"\t\"TEND3\" : {TEND3},\n")
        file.write(f"\t\"BETA3\" : {BETA3},\n")
        file.write(f"\t\"ITYPE3\" : {ITYPE3},\n")
        file.write(f"\t\"TSTOP\" : {TSTOP},\n")
        file.write(f"\t\"ILOAD\" : {ILOAD},\n")
        file.write(f"\t\"ICOMP\" : {ICOMP},\n")
        file.write(f"\t\"RHOS\" : {RHOS},\n")
        file.write(f"\t\"PHI0\" : {PHI0},\n")
        file.write(f"\t\"CONST\" : {CONST}\n")
        file.write("}\n")

SAVE_PARAMS=False
        
#
# Main function
#
def main(cfg_file=None):
    global SAVE_PARAMS
    if not cfg_file:
        #collect parameters from standard input(keyboard)
        while(1):
            collect_parameters_from_stdin()
            tmp = get_yes_or_no_from_stdin(f' Alter parameters? ', True)
            if tmp and tmp=='n':
                break
            else:
                continue
        SAVE_PARAMS=True
    else:
        read_cfg_file(cfg_file)
    
    run()#run the model
   
def check_cfg():
    if IRIFT >= 1:
        assert TBEG1 > TEND1
        
    if IRIFT >= 2:   
        assert TBEG2 > TEND2
        
    if IRIFT >= 3:
        assert TBEG3 > TEND3
    
def run():
    global K, AMPL, WAVEN, SINIT, TRINIT, DATA, TEM
    K = 1
    AMPL = 1 
    WAVEN = 1.
    SINIT = 0.
    TRINIT = 0.
    
    check_cfg()
    
    DATA=[ [ None for y in range( NDIM2 ) ] 
             for x in range( NDIM1 ) ] 
    TEM=[ [ None for y in range( NDIM2 ) ] 
             for x in range( NDIM1 ) ] 
    
    calculate_TL()
    
    nondimensionalise_parameters()
    #print(DTAUS, DTSTOP, DTST, TIMSC)
    
    setup_constants()
    #print(RHOA, CONBET, CONTRM, CONB1, CONB2, GP2, GP3)

    BETA1 = BETMIN
    
    while(1):
        if IRIFT == 3:
            BETA12 = BETA1*BETA2
            
        set_dless_constants(BETA1)

        KRR = 0
        TIMSTR(GP1,DTEND1,DTSTOP,DTAUS,CONBET,CONTRM,
            DTBEG2,DTBEG3,DTST,ITYPE1,IHFLW,KRR,IRIFT,
            BETA1,CONB1,CONB2) #fixme: this function call needs improvement

        if IRIFT >= 2:
            KRR = 1
            TIMSTR(GP2,DTEND2,DTSTOP,DTAUS,CONBET,CONTRM,
                DTBEG2,DTBEG3,DTST,ITYPE2,IHFLW,KRR,IRIFT,
                BETA1,CONB1,CONB2)

        if IRIFT >= 3:
            KRR = 2
            TIMSTR(GP3,DTEND3,DTSTOP,DTAUS,CONBET,CONTRM,
                DTBEG3,DTBEG3,DTST,ITYPE3,IHFLW,KRR,IRIFT,
                BETA12,CONB1,CONB2) 
            
        BETA1 = BETA1 + BETINC
        if BETA1 <= BETMAX:
            K = K+1
            DATA[0][K-1] = None #???
            DATA[1][K-1] = None #???
            K = K + 1
            continue
      
        break
        
    # Dimensionalise time and Heat Flux (as appropriate).
    for I in range(1, K+1):
        if DATA[0][I-1] is not None:
            DATA[0][I-1] = -DATA[0][I-1]*TIMSC  + TBEG1
            if IHFLW == 1:
                DATA[2][I-1] = DATA[2][I-1]*CONDTY
                DATA[3][I-1] = DATA[3][I-1]/TIMSC

    #Compaction if required.
    if ICOMP == 1:
        COMPAC(K)
 
    #save the data
    WRITE(K)
    
    #save the parameters for reference
    if SAVE_PARAMS:
        write_out_parameters()
    
    #print(DATA[1])
    print('done!')
  
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()








