C***********************************************************
C      RANDOM NUMBER GENERATOR: HP11C RANDOM ROUTINES      *
C      COUPLED WITH ANOTHER RANDOM ROUTINE, ONE CALCUALTES *
C      THE RANDOM SEED FOR THE OTHER                       *
C***********************************************************
        SUBROUTINE RANDM(NSEED,RAN)
           implicit double precision(a-h,o-z)
           DOUBLE PRECISION X,RA,RAN
           INTEGER A,P,NSEED,B15,B16,XHI,XALO,LEFTLO,FHI,K,C,B,N
999        CONTINUE
             RA = 997.0*.5284163
             FRAN = RA - AINT(RA)
             NSEED = INT(NSEED*FRAN)
        DATA A/16807/,B15/32768/,B16/65536/,P/2147483647/
             XHI = NSEED/B16
             XALO = (NSEED-XHI*B16)*A
             LEFTLO = XALO/B16
             FHI = XHI*A + LEFTLO
             K = FHI/B15
       NSEED = (((XALO - LEFTLO*B16) - P) + (FHI - K*B15)*B16) + K
          RA = FLOAT(NSEED)*4.656612875E-10
          ran = abs(ra)
C*************** END RANDOM NUMBER SUBROUTINE ********************
       return
       end
