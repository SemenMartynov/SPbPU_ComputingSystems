/* Livermore Loops coded in C        La//test File Modification  20 Oct 92,
 *  by Tim Peters, Kendall Square Res. Corp. tim@ksr.com, ksr!tim@uunet.uu.net
 *     SUBROUTINE KERNEL( TK)  replaces the Fortran routine in LFK //Test program.
 ************************************************************************
 *                                                                      *
 *            KERNEL     executes 24 samples of "C" computation         *
 *                                                                      *
 *                TK(1) - total cpu time to execute only the 24 kernels.*
 *                TK(2) - total Flops executed by the 24 Kernels        *
 *                                                                      *
 ************************************************************************
 *                                                                      *
 *     L. L. N. L.   " C "   K E R N E L S:   M F L O P S               *
 *                                                                      *
 *     These kernels measure   " C "   numerical computation            *
 *     rates for  a  spectrum  of  cpu-limited computational            *
 *     structures or benchmarks.   Mathematical  through-put            *
 *     is measured  in  units  of millions of floating-point            *
 *     operations executed per second, called Megaflops/sec.            *
 *                                                                      *
 *     Fonzi's Law: There is not now and there never will be a language *
 *                  in which it is the least bit difficult to write     *
 *                  bad programs.                                       *
 *                                                    F.H.MCMAHON  1972 *
 ************************************************************************
 *Originally from  Greg Astfalk, AT&T, P.O.Box 900, Princeton, NJ. 08540*
 *               by way of Frank McMahon (LLNL).                        *
 *                                                                      *
 *                               REFERENCE                              *
 *                                                                      *
 *              F.H.McMahon,   The Livermore Fortran Kernels:           *
 *              A Computer //Test Of The Numerical Performance Range,     *
 *              Lawrence Livermore National Laboratory,                 *
 *              Livermore, California, UCRL-53745, December 1986.       *
 *                                                                      *
 *       from:  National Technical Information Service                  *
 *              U.S. Department of Commerce                             *
 *              5285 Port Royal Road                                    *
 *              Springfield, VA.  22161                                 *
 *                                                                      *
 *    Changes made to correct many array subscripting problems,         *
 *      make more readable (added #define's), include the original      *
 *      FORTRAN versions of the runs as comments, and make more         *
 *      portable by Kelly O'Hair (LLNL) and Chuck Rasbold (LLNL).       *
 *                                                                      *
 ************************************************************************
 */

#include <stdio.h>
#include <math.h>

/* Define the structs or COMMON BLOCKS */
struct {
    long Mk;
    long Ik;
    long Im;
    long Ml;
    long Il;
    long Mruns;
    long Nruns;
    long Jr;
    long Npfs[47][3][8];
} ALPHA ;
#define mk     ALPHA.Mk
#define ik     ALPHA.Ik
#define im     ALPHA.Im
#define ml     ALPHA.Ml
#define il     ALPHA.Il
#define mruns  ALPHA.Mruns;
#define nruns  ALPHA.Nruns;
#define jr     ALPHA.Jr
#define npfs   ALPHA.Npfs

struct {
    double Tic;
    double Times[47][3][8];
    double See[3][8][3][5];
    double Terrs[47][3][8];
    double Csums[47][3][8];
    double Fopn[47][3][8];
    double Dos[47][3][8];
} BETA ;
#define tic     BETA.Tic
#define times   BETA.Times
#define see     BETA.See
#define terrs   BETA.Terrs
#define csums   BETA.Csums
#define fopn    BETA.Fopn
#define dos     BETA.Dos

struct {
    long Ion;
    long J5;
    long K2;
    long K3;
    long MULTI;
    long Laps;
    long Loop;
    long M;
    long Kr;
    long It;
    long N13h;
    long Ibuf;
    long Npass;
    long Nfail;
    long N;
    long N1;
    long N2;
    long N13;
    long N213;
    long N813;
    long N14;
    long N16;
    long N416;
    long N21;
    long Nt1;
    long Nt2;
} SPACES ;
#define  ion    SPACES.Ion
#define  j5     SPACES.J5
#define  k2     SPACES.K2
#define  k3     SPACES.K3
#define  multi  SPACES.MULTI
#define  laps   SPACES.Laps
#define  loop   SPACES.Loop
#define  m      SPACES.M
#define  kr     SPACES.Kr
#define  it     SPACES.It
#define  n13h   SPACES.N13h
#define  ibuf   SPACES.Ibuf
#define  npass  SPACES.Npass
#define  nfail  SPACES.Nfail
#define  n      SPACES.N
#define  n1     SPACES.N1
#define  n2     SPACES.N2
#define  n13    SPACES.N13
#define  n213   SPACES.N213
#define  n813   SPACES.N813
#define  n14    SPACES.N14
#define  n16    SPACES.N16
#define  n416   SPACES.N416
#define  n21    SPACES.N21
#define  nt1    SPACES.Nt1
#define  nt2    SPACES.Nt2

struct {
    double A11;
    double A12;
    double A13;
    double A21;
    double A22;
    double A23;
    double A31;
    double A32;
    double A33;
    double Ar;
    double Br;
    double C0;
    double Cr;
    double Di;
    double Dk;
    double Dm22;
    double Dm23;
    double Dm24;
    double Dm25;
    double Dm26;
    double Dm27;
    double Dm28;
    double Dn;
    double E3;
    double E6;
    double Expmax;
    double Flx;
    double Q;
    double Qa;
    double R;
    double Ri;
    double S;
    double Scale;
    double Sig;
    double Stb5;
    double T;
    double Xnc;
    double Xnei;
    double Xnm;
} SPACER ;
#define  a11     SPACER.A11
#define  a12     SPACER.A12
#define  a13     SPACER.A13
#define  a21     SPACER.A21
#define  a22     SPACER.A22
#define  a23     SPACER.A23
#define  a31     SPACER.A31
#define  a32     SPACER.A32
#define  a33     SPACER.A33
#define  ar      SPACER.Ar
#define  br      SPACER.Br
#define  c0      SPACER.C0
#define  cr      SPACER.Cr
#define  di      SPACER.Di
#define  dk      SPACER.Dk
#define  dm22    SPACER.Dm22
#define  dm23    SPACER.Dm23
#define  dm24    SPACER.Dm24
#define  dm25    SPACER.Dm25
#define  dm26    SPACER.Dm26
#define  dm27    SPACER.Dm27
#define  dm28    SPACER.Dm28
#define  dn      SPACER.Dn
#define  e3      SPACER.E3
#define  e6      SPACER.E6
#define  expmax  SPACER.Expmax
#define  flx     SPACER.Flx
#define  q       SPACER.Q
#define  qa      SPACER.Qa
#define  r       SPACER.R
#define  ri      SPACER.Ri
#define  s       SPACER.S
#define  scale   SPACER.Scale
#define  sig     SPACER.Sig
#define  stb5    SPACER.Stb5
#define  t       SPACER.T
#define  xnc     SPACER.Xnc
#define  xnei    SPACER.Xnei
#define  xnm     SPACER.Xnm

struct {
    double Time[47];
    double Csum[47];
    double Ww[47];
    double Wt[47];
    double Ticks;
    double Fr[9];
    double Terr1[47];
    double Sumw[7];
    double Start;
    double Skale[47];
    double Bias[47];
    double Ws[95];
    double Total[47];
    double Flopn[47];
    long Iq[7];
    long Npf;
    long Npfs1[47];
} SPACE0 ;
#define  time    SPACE0.Time
#define  csum    SPACE0.Csum
#define  ww      SPACE0.Ww
#define  wt      SPACE0.Wt
#define  ticks   SPACE0.Ticks
#define  fr      SPACE0.Fr
#define  terr1   SPACE0.Terr1
#define  sumw    SPACE0.Sumw
#define  start   SPACE0.Start
#define  skale   SPACE0.Skale
#define  bias    SPACE0.Bias
#define  ws      SPACE0.Ws
#define  total   SPACE0.Total
#define  flopn   SPACE0.Flopn
#define  iq      SPACE0.Iq
#define  npf     SPACE0.Npf
#define  npfs1   SPACE0.Npfs1

struct {
    double Wtp[3];
    long Mult[3];
    long Ispan[3][47];
    long Ipass[3][47];
} SPACEI ;
#define wtp    SPACEI.Wtp
#define mult   SPACEI.Mult
#define ispan  SPACEI.Ispan
#define ipass  SPACEI.Ipass

struct {
    long E[96];
    long F[96];
    long Ix[1001];
    long Ir[1001];
    long Zone[300];
} ISPACE ;
#define e    ISPACE.E
#define f    ISPACE.F
#define ix   ISPACE.Ix
#define ir   ISPACE.Ir
#define zone ISPACE.Zone

struct {
    double U[1001];
    double V[1001];
    double W[1001];
    double X[1001];
    double Y[1001];
    double Z[1001];
    double G[1001];
    double Du1[101];
    double Du2[101];
    double Du3[101];
    double Grd[1001];
    double Dex[1001];
    double Xi[1001];
    double Ex[1001];
    double Ex1[1001];
    double Dex1[1001];
    double Vx[1001];
    double Xx[1001];
    double Rx[1001];
    double Rh[2048];
    double Vsp[101];
    double Vstp[101];
    double Vxne[101];
    double Vxnd[101];
    double Ve3[101];
    double Vlr[101];
    double Vlin[101];
    double B5[101];
    double Plan[300];
    double D[300];
    double Sa[101];
    double Sb[101];
} SPACE1 ;
#define  u    SPACE1.U
#define  v    SPACE1.V
#define  w    SPACE1.W
#define  x    SPACE1.X
#define  y    SPACE1.Y
#define  z    SPACE1.Z
#define  g    SPACE1.G
#define  du1  SPACE1.Du1
#define  du2  SPACE1.Du2
#define  du3  SPACE1.Du3
#define  grd  SPACE1.Grd
#define  dex  SPACE1.Dex
#define  xi   SPACE1.Xi
#define  ex   SPACE1.Ex
#define  ex1  SPACE1.Ex1
#define  dex1 SPACE1.Dex1
#define  vx   SPACE1.Vx
#define  xx   SPACE1.Xx
#define  rx   SPACE1.Rx
#define  rh   SPACE1.Rh
#define  vsp  SPACE1.Vsp
#define  vstp SPACE1.Vstp
#define  vxne SPACE1.Vxne
#define  vxnd SPACE1.Vxnd
#define  ve3  SPACE1.Ve3
#define  vlr  SPACE1.Vlr
#define  vlin SPACE1.Vlin
#define  b5   SPACE1.B5
#define  plan SPACE1.Plan
#define  d    SPACE1.D
#define  sa   SPACE1.Sa
#define  sb   SPACE1.Sb

struct {
    double P[512][4];
    double Px[101][25];
    double Cx[101][25];
    double Vy[25][101];
    double Vh[7][101];
    double Vf[7][101];
    double Vg[7][101];
    double Vs[7][101];
    double Za[7][101];
    double Zp[7][101];
    double Zq[7][101];
    double Zr[7][101];
    double Zm[7][101];
    double Zb[7][101];
    double Zu[7][101];
    double Zv[7][101];
    double Zz[7][101];
    double B[64][64];
    double C[64][64];
    double H[64][64];
    double U1[2][101][5];
    double U2[2][101][5];
    double U3[2][101][5];
} SPACE2 ;
#define  p       SPACE2.P
#define  px      SPACE2.Px
#define  cx      SPACE2.Cx
#define  vy      SPACE2.Vy
#define  vh      SPACE2.Vh
#define  vf      SPACE2.Vf
#define  vg      SPACE2.Vg
#define  vs      SPACE2.Vs
#define  za      SPACE2.Za
#define  zp      SPACE2.Zp
#define  zq      SPACE2.Zq
#define  zr      SPACE2.Zr
#define  zm      SPACE2.Zm
#define  zb      SPACE2.Zb
#define  zu      SPACE2.Zu
#define  zv      SPACE2.Zv
#define  zz      SPACE2.Zz
#define  b       SPACE2.B
#define  c       SPACE2.C
#define  h       SPACE2.H
#define  u1      SPACE2.U1
#define  u2      SPACE2.U2
#define  u3      SPACE2.U3

/* KERNEL routine */

int main(int argc, char** argv)
{

#pragma nodyneqv
#pragma o=i

    long argument , k , l , ipnt , ipntp , i;
    long lw , j , nl1 , nl2 , kx , ky , ip , kn;
    long i1 , j1 , i2 , j2 , nz , ink , jn , kb5i;
    long ii , lb , j4 , ng;
    double tmp , temp, sum, som;
    char name[8];

    /*
     *******************************************************************
     *   Kernel 1 -- hydro fragment
     *******************************************************************
     *       DO 1 L = 1,Loop
     *       DO 1 k = 1,n
     *  1       X(k)= Q + Y(k)*(R*ZX(k+10) + T*ZX(k+11))
     */

    for ( l=1 ; l<=loop ; l++ ) {
        for ( k=0 ; k<n ; k++ ) {
            x[k] = q + y[k]*( r*z[k+10] + t*z[k+11] );
        }
    }
    argument = 1;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 2 -- ICCG excerpt (Incomplete Cholesky Conjugate Gradient)
     *******************************************************************
     *    DO 200  L= 1,Loop
     *        II= n
     *     IPNTP= 0
     *222   IPNT= IPNTP
     *     IPNTP= IPNTP+II
     *        II= II/2
     *         i= IPNTP
     CDIR$ IVDEP
     *    DO 2 k= IPNT+2,IPNTP,2
     *         i= i+1
     *  2   X(i)= X(k) - V(k)*X(k-1) - V(k+1)*X(k+1)
     *        IF( II.GT.1) GO TO 222
     *200 CONTINUE
     */

    for ( l=1 ; l<=loop ; l++ ) {
        ii = n;
        ipntp = 0;
        do {
            ipnt = ipntp;
            ipntp += ii;
            ii /= 2;
            i = ipntp - 1;
#pragma nohazard
            for ( k=ipnt+1 ; k<ipntp ; k=k+2 ) {
                i++;
                x[i] = x[k] - v[k  ]*x[k-1] - v[k+1]*x[k+1];
            }
        } while ( ii>0 );
    }
    argument = 2;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 3 -- inner product
     *******************************************************************
     *    DO 3 L= 1,Loop
     *         Q= 0.0
     *    DO 3 k= 1,n
     *  3      Q= Q + Z(k)*X(k)
     */

    for ( l=1 ; l<=loop ; l++ ) {
        q = 0.0;
        for ( k=0 ; k<n ; k++ ) {
            q += z[k]*x[k];
        }
    }
    argument = 3;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 4 -- banded linear equations
     *******************************************************************
     *            m= (1001-7)/2
     *    DO 444  L= 1,Loop
     *    DO 444  k= 7,1001,m
     *           lw= k-6
     *         temp= X(k-1)
     CDIR$ IVDEP
     *    DO   4  j= 5,n,5
     *       temp  = temp   - XZ(lw)*Y(j)
     *  4        lw= lw+1
     *       X(k-1)= Y(5)*temp
     *444 CONTINUE
     */

    m = ( 1001-7 )/2;
    for ( l=1 ; l<=loop ; l++ ) {
        for ( k=6 ; k<1001 ; k=k+m ) {
            lw = k - 6;
            temp = x[k-1];
#pragma nohazard
            for ( j=4 ; j<n ; j=j+5 ) {
                temp -= x[lw]*y[j];
                lw++;
            }
            x[k-1] = y[4]*temp;
        }
    }
    argument = 4;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 5 -- tri-diagonal elimination, below diagonal
     *******************************************************************
     *    DO 5 L = 1,Loop
     *    DO 5 i = 2,n
     *  5    X(i)= Z(i)*(Y(i) - X(i-1))
     */

    for ( l=1 ; l<=loop ; l++ ) {
        for ( i=1 ; i<n ; i++ ) {
            x[i] = z[i]*( y[i] - x[i-1] );
        }
    }
    argument = 5;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 6 -- general linear recurrence equations
     *******************************************************************
     *    DO  6  L= 1,Loop
     *    DO  6  i= 2,n
     *    DO  6  k= 1,i-1
     *        W(i)= W(i)  + B(i,k) * W(i-k)
     *  6 CONTINUE
     */

    for ( l=1 ; l<=loop ; l++ ) {
        for ( i=1 ; i<n ; i++ ) {
            for ( k=0 ; k<i ; k++ ) {
                w[i] += b[k][i] * w[(i-k)-1];
            }
        }
    }
    argument = 6;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 7 -- equation of state fragment
     *******************************************************************
     *    DO 7 L= 1,Loop
     *    DO 7 k= 1,n
     *      X(k)=     U(k  ) + R*( Z(k  ) + R*Y(k  )) +
     *   .        T*( U(k+3) + R*( U(k+2) + R*U(k+1)) +
     *   .        T*( U(k+6) + R*( U(k+5) + R*U(k+4))))
     *  7 CONTINUE
     */

    for ( l=1 ; l<=loop ; l++ ) {
#pragma nohazard
        for ( k=0 ; k<n ; k++ ) {
            x[k] = u[k] + r*( z[k] + r*y[k] ) +
                   t*( u[k+3] + r*( u[k+2] + r*u[k+1] ) +
                      t*( u[k+6] + r*( u[k+5] + r*u[k+4] ) ) );
        }
    }
    argument = 7;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 8 -- ADI integration
     *******************************************************************
     *    DO  8      L = 1,Loop
     *             nl1 = 1
     *             nl2 = 2
     *    DO  8     kx = 2,3
     CDIR$ IVDEP
     *    DO  8     ky = 2,n
     *          DU1(ky)=U1(kx,ky+1,nl1)  -  U1(kx,ky-1,nl1)
     *          DU2(ky)=U2(kx,ky+1,nl1)  -  U2(kx,ky-1,nl1)
     *          DU3(ky)=U3(kx,ky+1,nl1)  -  U3(kx,ky-1,nl1)
     *    U1(kx,ky,nl2)=U1(kx,ky,nl1) +A11*DU1(ky) +A12*DU2(ky) +A13*DU3(ky)
     *   .       + SIG*(U1(kx+1,ky,nl1) -2.*U1(kx,ky,nl1) +U1(kx-1,ky,nl1))
     *    U2(kx,ky,nl2)=U2(kx,ky,nl1) +A21*DU1(ky) +A22*DU2(ky) +A23*DU3(ky)
     *   .       + SIG*(U2(kx+1,ky,nl1) -2.*U2(kx,ky,nl1) +U2(kx-1,ky,nl1))
     *    U3(kx,ky,nl2)=U3(kx,ky,nl1) +A31*DU1(ky) +A32*DU2(ky) +A33*DU3(ky)
     *   .       + SIG*(U3(kx+1,ky,nl1) -2.*U3(kx,ky,nl1) +U3(kx-1,ky,nl1))
     *  8 CONTINUE
     */

    for ( l=1 ; l<=loop ; l++ ) {
        nl1 = 0;
        nl2 = 1;
        for ( kx=1 ; kx<3 ; kx++ ){
#pragma nohazard
            for ( ky=1 ; ky<n ; ky++ ) {
               du1[ky] = u1[nl1][ky+1][kx] - u1[nl1][ky-1][kx];
               du2[ky] = u2[nl1][ky+1][kx] - u2[nl1][ky-1][kx];
               du3[ky] = u3[nl1][ky+1][kx] - u3[nl1][ky-1][kx];
               u1[nl2][ky][kx]=
                  u1[nl1][ky][kx]+a11*du1[ky]+a12*du2[ky]+a13*du3[ky] + sig*
                     (u1[nl1][ky][kx+1]-2.0*u1[nl1][ky][kx]+u1[nl1][ky][kx-1]);
               u2[nl2][ky][kx]=
                  u2[nl1][ky][kx]+a21*du1[ky]+a22*du2[ky]+a23*du3[ky] + sig*
                     (u2[nl1][ky][kx+1]-2.0*u2[nl1][ky][kx]+u2[nl1][ky][kx-1]);
               u3[nl2][ky][kx]=
                  u3[nl1][ky][kx]+a31*du1[ky]+a32*du2[ky]+a33*du3[ky] + sig*
                     (u3[nl1][ky][kx+1]-2.0*u3[nl1][ky][kx]+u3[nl1][ky][kx-1]);
            }
        }
    }
    argument = 8;
    //TEST( &argument );

    /*
     *******************************************************************
     *   Kernel 9 -- integrate predictors
     *******************************************************************
     *    DO 9  L = 1,Loop
     *    DO 9  i = 1,n
     *    PX( 1,i)= DM28*PX(13,i) + DM27*PX(12,i) + DM26*PX(11,i) +
     *   .          DM25*PX(10,i) + DM24*PX( 9,i) + DM23*PX( 8,i) +
     *   .          DM22*PX( 7,i) +  C0*(PX( 5,i) +      PX( 6,i))+ PX( 3,i)
     *  9 CONTINUE
     */

    for ( l=1 ; l<=loop ; l++ ) {
        for ( i=0 ; i<n ; i++ ) {
            px[i][0] = dm28*px[i][12] + dm27*px[i][11] + dm26*px[i][10] +
                       dm25*px[i][ 9] + dm24*px[i][ 8] + dm23*px[i][ 7] +
                       dm22*px[i][ 6] + c0*( px[i][ 4] + px[i][ 5]) + px[i][ 2];
        }
    }
    argument = 9;
    //TEST( &argument );

    return;
}
