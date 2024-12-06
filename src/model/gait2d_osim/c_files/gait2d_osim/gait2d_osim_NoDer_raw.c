/* The name of this program is gait2d_osim/gait2d_osim_NoDer_raw.c */
/* Created by Autolev 4.1 on Thu Nov 28 15:08:29 2024 */

#include <ctype.h> 
#include <math.h>  
#include <stdarg.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#define  _NAN      9.99999999999999E+305

void     eqns1     (void);
void     output    (FILE *Fptr[] );
void     readf     (FILE *Fp, double *next, ...);
void     pgets     (FILE *Fp, double *x);
void     writef    (FILE *Fp, char format[], ...);

double   par__airdrag,par__back_x,par__back_y,par__bodyweight,par__gravity_x,par__gravity_y,
  par__Lankle_y,par__Lfemur_CMx,par__Lfemur_CMy,par__Lfemur_Izz,par__Lfemur_M,
  par__Lfoot_CMx,par__Lfoot_CMy,par__Lfoot_Izz,par__Lfoot_M,par__Lhip_x,par__Lhip_y,
  par__Lknee_x1,par__Lknee_x2,par__Lknee_x3,par__Lknee_x4,par__Lknee_x5,par__Lknee_y1,
  par__Lknee_y2,par__Lknee_y3,par__Lknee_y4,par__Lknee_y5,par__Ltibia_CMx,par__Ltibia_CMy,
  par__Ltibia_Izz,par__Ltibia_M,par__pelvis_CMx,par__pelvis_CMy,par__pelvis_Izz,
  par__pelvis_M,par__Rankle_y,par__Rfemur_CMx,par__Rfemur_CMy,par__Rfemur_Izz,
  par__Rfemur_M,par__Rfoot_CMx,par__Rfoot_CMy,par__Rfoot_Izz,par__Rfoot_M,par__Rhip_x,
  par__Rhip_y,par__Rknee_x1,par__Rknee_x2,par__Rknee_x3,par__Rknee_x4,par__Rknee_x5,
  par__Rknee_y1,par__Rknee_y2,par__Rknee_y3,par__Rknee_y4,par__Rknee_y5,par__Rtibia_CMx,
  par__Rtibia_CMy,par__Rtibia_Izz,par__Rtibia_M,par__torso_CMx,par__torso_CMy,
  par__torso_Izz,par__torso_M,par__wind,G1,G2,G3,G4,G5,G6,q1,q2,q3,q4,q5,q6,q7,
  q8,q9,q1p,q2p,q3p,q4p,q5p,q6p,q7p,q8p,q9p,q1pp,q2pp,q3pp,q4pp,q5pp,q6pp,q7pp,
  q8pp,q9pp;
double   s,sx,sy;
double   Pi,DEGtoRAD,RADtoDEG,z[545],fk[48],dfk_dq[48][9],fkdot[48],dfkdot_dq[
  48][9],f[9],Encode[5];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait2d_osim/gait2d_osim_NoDer_raw.in");
  else sprintf(fileName, "gait2d_osim/gait2d_osim_NoDer_raw.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&par__airdrag,&par__back_x,&par__back_y,&par__bodyweight,&par__gravity_x,&
  par__gravity_y,&par__Lankle_y,&par__Lfemur_CMx,&par__Lfemur_CMy,&par__Lfemur_Izz,&
  par__Lfemur_M,&par__Lfoot_CMx,&par__Lfoot_CMy,&par__Lfoot_Izz,&par__Lfoot_M,&
  par__Lhip_x,&par__Lhip_y,&par__Lknee_x1,&par__Lknee_x2,&par__Lknee_x3,&par__Lknee_x4,&
  par__Lknee_x5,&par__Lknee_y1,&par__Lknee_y2,&par__Lknee_y3,&par__Lknee_y4,&par__Lknee_y5,&
  par__Ltibia_CMx,&par__Ltibia_CMy,&par__Ltibia_Izz,&par__Ltibia_M,&par__pelvis_CMx,&
  par__pelvis_CMy,&par__pelvis_Izz,&par__pelvis_M,&par__Rankle_y,&par__Rfemur_CMx,&
  par__Rfemur_CMy,&par__Rfemur_Izz,&par__Rfemur_M,&par__Rfoot_CMx,&par__Rfoot_CMy,&
  par__Rfoot_Izz,&par__Rfoot_M,&par__Rhip_x,&par__Rhip_y,&par__Rknee_x1,&par__Rknee_x2,&
  par__Rknee_x3,&par__Rknee_x4,&par__Rknee_x5,&par__Rknee_y1,&par__Rknee_y2,&par__Rknee_y3,&
  par__Rknee_y4,&par__Rknee_y5,&par__Rtibia_CMx,&par__Rtibia_CMy,&par__Rtibia_Izz,&
  par__Rtibia_M,&par__torso_CMx,&par__torso_CMy,&par__torso_Izz,&par__torso_M,&
  par__wind,&G1,&G2,&G3,&G4,&G5,&G6,&q1,&q2,&q3,&q4,&q5,&q6,&q7,&q8,&q9,&q1p,&
  q2p,&q3p,&q4p,&q5p,&q6p,&q7p,&q8p,&q9p,&q1pp,&q2pp,&q3pp,&q4pp,&q5pp,&q6pp,&
  q7pp,&q8pp,&q9pp,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file gait2d_osim/gait2d_osim_NoDer_raw.in" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (void)
{

/* Evaluate constants */
  z[1] = cos(q1);
  z[2] = sin(q1);
  z[3] = cos(q4);
  z[4] = sin(q4);
  z[5] = z[1]*z[3] - z[2]*z[4];
  z[6] = -z[1]*z[4] - z[2]*z[3];
  z[7] = z[1]*z[4] + z[2]*z[3];
  z[8] = par__Rhip_x*z[4] - par__Rhip_y*z[3];
  z[9] = par__Rhip_x*z[3] + par__Rhip_y*z[4];
  z[10] = z[8] - par__Rfemur_CMy;
  z[11] = par__Rfemur_CMx + z[9];
  z[12] = par__Rknee_x5 + par__Rknee_x4*q5 + par__Rknee_x1*pow(q5,4) + par__Rknee_x2*
  pow(q5,3) + par__Rknee_x3*pow(q5,2);
  z[13] = par__Rknee_y5 + par__Rknee_y4*q5 + par__Rknee_y1*pow(q5,4) + par__Rknee_y2*
  pow(q5,3) + par__Rknee_y3*pow(q5,2);
  z[14] = z[8] - z[13];
  z[15] = z[9] + z[12];
  z[16] = cos(q5);
  z[17] = sin(q5);
  z[18] = z[5]*z[16] + z[6]*z[17];
  z[19] = z[6]*z[16] - z[5]*z[17];
  z[20] = z[5]*z[17] + z[7]*z[16];
  z[21] = z[5]*z[16] - z[7]*z[17];
  z[22] = z[12]*z[17] - z[13]*z[16];
  z[23] = z[14]*z[16] + z[15]*z[17];
  z[24] = z[12]*z[16] + z[13]*z[17];
  z[25] = z[15]*z[16] - z[14]*z[17];
  z[26] = z[22] - par__Rtibia_CMy;
  z[27] = z[23] - par__Rtibia_CMy;
  z[28] = par__Rtibia_CMx + z[24];
  z[29] = par__Rtibia_CMx + z[25];
  z[30] = z[22] - par__Rankle_y;
  z[31] = z[23] - par__Rankle_y;
  z[32] = cos(q6);
  z[33] = sin(q6);
  z[34] = z[18]*z[32] + z[19]*z[33];
  z[35] = z[19]*z[32] - z[18]*z[33];
  z[36] = z[20]*z[32] + z[21]*z[33];
  z[37] = z[21]*z[32] - z[20]*z[33];
  z[38] = z[24]*z[33] + z[30]*z[32];
  z[39] = z[25]*z[33] + z[31]*z[32];
  z[40] = par__Rankle_y*z[32];
  z[41] = z[24]*z[32] - z[30]*z[33];
  z[42] = z[25]*z[32] - z[31]*z[33];
  z[43] = par__Rankle_y*z[33];
  z[44] = z[38] - par__Rfoot_CMy;
  z[45] = z[39] - par__Rfoot_CMy;
  z[46] = -par__Rfoot_CMy - z[40];
  z[47] = par__Rfoot_CMx + z[41];
  z[48] = par__Rfoot_CMx + z[42];
  z[49] = par__Rfoot_CMx + z[43];
  z[50] = cos(q7);
  z[51] = sin(q7);
  z[52] = z[1]*z[50] - z[2]*z[51];
  z[53] = -z[1]*z[51] - z[2]*z[50];
  z[54] = z[1]*z[51] + z[2]*z[50];
  z[55] = par__Lhip_x*z[51] - par__Lhip_y*z[50];
  z[56] = par__Lhip_x*z[50] + par__Lhip_y*z[51];
  z[57] = z[55] - par__Lfemur_CMy;
  z[58] = par__Lfemur_CMx + z[56];
  z[59] = par__Lknee_x5 + par__Lknee_x4*q8 + par__Lknee_x1*pow(q8,4) + par__Lknee_x2*
  pow(q8,3) + par__Lknee_x3*pow(q8,2);
  z[60] = par__Lknee_y5 + par__Lknee_y4*q8 + par__Lknee_y1*pow(q8,4) + par__Lknee_y2*
  pow(q8,3) + par__Lknee_y3*pow(q8,2);
  z[61] = z[55] - z[60];
  z[62] = z[56] + z[59];
  z[63] = cos(q8);
  z[64] = sin(q8);
  z[65] = z[52]*z[63] + z[53]*z[64];
  z[66] = z[53]*z[63] - z[52]*z[64];
  z[67] = z[52]*z[64] + z[54]*z[63];
  z[68] = z[52]*z[63] - z[54]*z[64];
  z[69] = z[59]*z[64] - z[60]*z[63];
  z[70] = z[61]*z[63] + z[62]*z[64];
  z[71] = z[59]*z[63] + z[60]*z[64];
  z[72] = z[62]*z[63] - z[61]*z[64];
  z[73] = z[69] - par__Ltibia_CMy;
  z[74] = z[70] - par__Ltibia_CMy;
  z[75] = par__Ltibia_CMx + z[71];
  z[76] = par__Ltibia_CMx + z[72];
  z[77] = z[69] - par__Lankle_y;
  z[78] = z[70] - par__Lankle_y;
  z[79] = cos(q9);
  z[80] = sin(q9);
  z[81] = z[65]*z[79] + z[66]*z[80];
  z[82] = z[66]*z[79] - z[65]*z[80];
  z[83] = z[67]*z[79] + z[68]*z[80];
  z[84] = z[68]*z[79] - z[67]*z[80];
  z[85] = z[71]*z[80] + z[77]*z[79];
  z[86] = z[72]*z[80] + z[78]*z[79];
  z[87] = par__Lankle_y*z[79];
  z[88] = z[71]*z[79] - z[77]*z[80];
  z[89] = z[72]*z[79] - z[78]*z[80];
  z[90] = par__Lankle_y*z[80];
  z[91] = z[85] - par__Lfoot_CMy;
  z[92] = z[86] - par__Lfoot_CMy;
  z[93] = -par__Lfoot_CMy - z[87];
  z[94] = par__Lfoot_CMx + z[88];
  z[95] = par__Lfoot_CMx + z[89];
  z[96] = par__Lfoot_CMx + z[90];
  z[97] = -par__back_y - par__torso_CMy;
  z[98] = par__back_x + par__torso_CMx;
  z[99] = z[2]*z[4] - z[1]*z[3];
  z[100] = par__Rknee_x4 + 2*par__Rknee_x3*q5 + 3*par__Rknee_x2*pow(q5,2) + 4*
  par__Rknee_x1*pow(q5,3);
  z[101] = par__Rknee_y4 + 2*par__Rknee_y3*q5 + 3*par__Rknee_y2*pow(q5,2) + 4*
  par__Rknee_y1*pow(q5,3);
  z[102] = z[6]*z[16] + z[17]*z[99];
  z[103] = z[16]*z[99] - z[6]*z[17];
  z[104] = -z[5]*z[16] - z[6]*z[17];
  z[105] = -z[5]*z[17] - z[7]*z[16];
  z[106] = z[32]*z[102] + z[33]*z[103];
  z[107] = z[19]*z[32] + z[33]*z[104];
  z[108] = z[32]*z[103] - z[33]*z[102];
  z[109] = z[32]*z[104] - z[19]*z[33];
  z[110] = -z[18]*z[32] - z[19]*z[33];
  z[111] = z[21]*z[32] + z[33]*z[105];
  z[112] = z[32]*z[105] - z[21]*z[33];
  z[113] = -z[20]*z[32] - z[21]*z[33];
  z[114] = z[2]*z[51] - z[1]*z[50];
  z[115] = par__Lknee_x4 + 2*par__Lknee_x3*q8 + 3*par__Lknee_x2*pow(q8,2) + 4*
  par__Lknee_x1*pow(q8,3);
  z[116] = par__Lknee_y4 + 2*par__Lknee_y3*q8 + 3*par__Lknee_y2*pow(q8,2) + 4*
  par__Lknee_y1*pow(q8,3);
  z[117] = z[53]*z[63] + z[64]*z[114];
  z[118] = z[63]*z[114] - z[53]*z[64];
  z[119] = -z[52]*z[63] - z[53]*z[64];
  z[120] = -z[52]*z[64] - z[54]*z[63];
  z[121] = z[79]*z[117] + z[80]*z[118];
  z[122] = z[66]*z[79] + z[80]*z[119];
  z[123] = z[79]*z[118] - z[80]*z[117];
  z[124] = z[79]*z[119] - z[66]*z[80];
  z[125] = -z[65]*z[79] - z[66]*z[80];
  z[126] = z[68]*z[79] + z[80]*z[120];
  z[127] = z[79]*z[120] - z[68]*z[80];
  z[128] = -z[67]*z[79] - z[68]*z[80];
  z[129] = -par__Rhip_x*z[2] - par__Rhip_y*z[1];
  z[130] = par__Rhip_x*z[1] - par__Rhip_y*z[2];
  z[131] = z[6]*z[12] + z[13]*z[99] - par__Rhip_x*z[2] - par__Rhip_y*z[1];
  z[132] = z[6]*z[12] + z[13]*z[99];
  z[133] = z[5]*z[100] + z[6]*z[101];
  z[134] = par__Rhip_x*z[1] + z[5]*z[12] + z[6]*z[13] - par__Rhip_y*z[2];
  z[135] = z[5]*z[12] + z[6]*z[13];
  z[136] = z[5]*z[101] + z[7]*z[100];
  z[137] = par__Rankle_y*z[103] + z[6]*z[12] + z[13]*z[99] - par__Rhip_x*z[2] - 
  par__Rhip_y*z[1];
  z[138] = par__Rankle_y*z[103] + z[6]*z[12] + z[13]*z[99];
  z[139] = par__Rankle_y*z[104] + z[5]*z[100] + z[6]*z[101];
  z[140] = par__Rankle_y*z[19] + par__Rhip_x*z[1] + z[5]*z[12] + z[6]*z[13] - 
  par__Rhip_y*z[2];
  z[141] = par__Rankle_y*z[19] + z[5]*z[12] + z[6]*z[13];
  z[142] = par__Rankle_y*z[105] + z[5]*z[101] + z[7]*z[100];
  z[143] = -par__Lhip_x*z[2] - par__Lhip_y*z[1];
  z[144] = par__Lhip_x*z[1] - par__Lhip_y*z[2];
  z[145] = z[53]*z[59] + z[60]*z[114] - par__Lhip_x*z[2] - par__Lhip_y*z[1];
  z[146] = z[53]*z[59] + z[60]*z[114];
  z[147] = z[52]*z[115] + z[53]*z[116];
  z[148] = par__Lhip_x*z[1] + z[52]*z[59] + z[53]*z[60] - par__Lhip_y*z[2];
  z[149] = z[52]*z[59] + z[53]*z[60];
  z[150] = z[52]*z[116] + z[54]*z[115];
  z[151] = par__Lankle_y*z[118] + z[53]*z[59] + z[60]*z[114] - par__Lhip_x*
  z[2] - par__Lhip_y*z[1];
  z[152] = par__Lankle_y*z[118] + z[53]*z[59] + z[60]*z[114];
  z[153] = par__Lankle_y*z[119] + z[52]*z[115] + z[53]*z[116];
  z[154] = par__Lankle_y*z[66] + par__Lhip_x*z[1] + z[52]*z[59] + z[53]*z[60] - 
  par__Lhip_y*z[2];
  z[155] = par__Lankle_y*z[66] + z[52]*z[59] + z[53]*z[60];
  z[156] = par__Lankle_y*z[120] + z[52]*z[116] + z[54]*z[115];
  z[157] = -par__back_x*z[2] - par__back_y*z[1];
  z[158] = par__back_x*z[1] - par__back_y*z[2];
  z[159] = z[2]*q1p;
  z[160] = z[1]*q1p;
  z[161] = q2p - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[162] = q3p + par__Rhip_x*z[1]*q1p - par__Rhip_y*z[2]*q1p;
  z[163] = z[6]*(q1p+q4p);
  z[164] = z[99]*(q1p+q4p);
  z[165] = z[5]*(q1p+q4p);
  z[166] = q2p + z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + 
  z[13]*z[99]*(q1p+q4p) - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[167] = q3p + par__Rhip_x*z[1]*q1p + z[5]*z[101]*q5p + z[7]*z[100]*q5p + 
  z[5]*z[12]*(q1p+q4p) + z[6]*z[13]*(q1p+q4p) - par__Rhip_y*z[2]*q1p;
  z[168] = z[19]*q5p + z[102]*q1p + z[102]*q4p;
  z[169] = z[103]*q1p + z[103]*q4p + z[104]*q5p;
  z[170] = z[18]*q1p + z[18]*q4p + z[21]*q5p;
  z[171] = z[19]*q1p + z[19]*q4p + z[105]*q5p;
  z[172] = q2p + z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + 
  z[13]*z[99]*(q1p+q4p) + par__Rankle_y*(z[103]*q1p+z[103]*q4p+z[104]*q5p) - 
  par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[173] = q3p + par__Rhip_x*z[1]*q1p + z[5]*z[101]*q5p + z[7]*z[100]*q5p + 
  z[5]*z[12]*(q1p+q4p) + z[6]*z[13]*(q1p+q4p) + par__Rankle_y*(z[19]*q1p+
  z[19]*q4p+z[105]*q5p) - par__Rhip_y*z[2]*q1p;
  z[174] = z[35]*q6p + z[106]*q1p + z[106]*q4p + z[107]*q5p;
  z[175] = z[108]*q1p + z[108]*q4p + z[109]*q5p + z[110]*q6p;
  z[176] = z[34]*q1p + z[34]*q4p + z[37]*q6p + z[111]*q5p;
  z[177] = z[35]*q1p + z[35]*q4p + z[112]*q5p + z[113]*q6p;
  z[178] = q2p - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[179] = q3p + par__Lhip_x*z[1]*q1p - par__Lhip_y*z[2]*q1p;
  z[180] = z[53]*(q1p+q7p);
  z[181] = z[114]*(q1p+q7p);
  z[182] = z[52]*(q1p+q7p);
  z[183] = q2p + z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p) - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[184] = q3p + par__Lhip_x*z[1]*q1p + z[52]*z[116]*q8p + z[54]*z[115]*q8p + 
  z[52]*z[59]*(q1p+q7p) + z[53]*z[60]*(q1p+q7p) - par__Lhip_y*z[2]*q1p;
  z[185] = z[66]*q8p + z[117]*q1p + z[117]*q7p;
  z[186] = z[118]*q1p + z[118]*q7p + z[119]*q8p;
  z[187] = z[65]*q1p + z[65]*q7p + z[68]*q8p;
  z[188] = z[66]*q1p + z[66]*q7p + z[120]*q8p;
  z[189] = q2p + z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p) + par__Lankle_y*(z[118]*q1p+z[118]*q7p+z[119]*q8p) - 
  par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[190] = q3p + par__Lhip_x*z[1]*q1p + z[52]*z[116]*q8p + z[54]*z[115]*q8p + 
  z[52]*z[59]*(q1p+q7p) + z[53]*z[60]*(q1p+q7p) + par__Lankle_y*(z[66]*q1p+
  z[66]*q7p+z[120]*q8p) - par__Lhip_y*z[2]*q1p;
  z[191] = z[82]*q9p + z[121]*q1p + z[121]*q7p + z[122]*q8p;
  z[192] = z[123]*q1p + z[123]*q7p + z[124]*q8p + z[125]*q9p;
  z[193] = z[81]*q1p + z[81]*q7p + z[84]*q9p + z[126]*q8p;
  z[194] = z[82]*q1p + z[82]*q7p + z[127]*q8p + z[128]*q9p;
  z[195] = q2p - par__back_x*z[2]*q1p - par__back_y*z[1]*q1p;
  z[196] = q3p + par__back_x*z[1]*q1p - par__back_y*z[2]*q1p;
  z[197] = (par__Rhip_x*z[1]-par__Rhip_y*z[2])*q1p;
  z[198] = (par__Rhip_x*z[2]+par__Rhip_y*z[1])*q1p;
  z[199] = par__Rhip_y*z[2]*q1p + z[6]*z[100]*q5p + z[99]*z[101]*q5p + z[7]*
  z[13]*(q1p+q4p) + z[12]*z[99]*(q1p+q4p) - par__Rhip_x*z[1]*q1p;
  z[200] = z[6]*z[100]*q5p + z[99]*z[101]*q5p + z[7]*z[13]*(q1p+q4p) + z[12]*
  z[99]*(q1p+q4p);
  z[201] = 2*par__Rknee_x3 + 6*par__Rknee_x2*q5 + 12*par__Rknee_x1*pow(q5,2);
  z[202] = 2*par__Rknee_y3 + 6*par__Rknee_y2*q5 + 12*par__Rknee_y1*pow(q5,2);
  z[203] = z[5]*z[201]*q5p + z[6]*z[202]*q5p + z[6]*z[100]*(q1p+q4p) + z[99]*
  z[101]*(q1p+q4p);
  z[204] = z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[99]*(q1p+q4p) - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[205] = z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[99]*(q1p+q4p);
  z[206] = z[5]*z[202]*q5p + z[7]*z[201]*q5p + z[5]*z[100]*(q1p+q4p) + z[6]*
  z[101]*(q1p+q4p);
  z[207] = z[7]*z[17] + z[16]*z[99];
  z[208] = z[103]*q5p + z[207]*q1p + z[207]*q4p;
  z[209] = z[7]*z[16] - z[17]*z[99];
  z[210] = -z[6]*z[16] - z[17]*z[99];
  z[211] = z[209]*q1p + z[209]*q4p + z[210]*q5p;
  z[212] = z[5]*z[17] - z[6]*z[16];
  z[213] = z[210]*q1p + z[210]*q4p + z[212]*q5p;
  z[214] = z[7]*z[17] - z[5]*z[16];
  z[215] = z[104]*q1p + z[104]*q4p + z[214]*q5p;
  z[216] = par__Rhip_y*z[2]*q1p + z[6]*z[100]*q5p + z[99]*z[101]*q5p + z[7]*
  z[13]*(q1p+q4p) + z[12]*z[99]*(q1p+q4p) + par__Rankle_y*(z[209]*q1p+z[209]*
  q4p+z[210]*q5p) - par__Rhip_x*z[1]*q1p;
  z[217] = z[6]*z[100]*q5p + z[99]*z[101]*q5p + z[7]*z[13]*(q1p+q4p) + z[12]*
  z[99]*(q1p+q4p) + par__Rankle_y*(z[209]*q1p+z[209]*q4p+z[210]*q5p);
  z[218] = z[5]*z[201]*q5p + z[6]*z[202]*q5p + z[6]*z[100]*(q1p+q4p) + z[99]*
  z[101]*(q1p+q4p) + par__Rankle_y*(z[210]*q1p+z[210]*q4p+z[212]*q5p);
  z[219] = z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[99]*(q1p+q4p) + par__Rankle_y*(z[103]*q1p+z[103]*q4p+z[104]*q5p) - par__Rhip_x*
  z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[220] = z[5]*z[100]*q5p + z[6]*z[101]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[99]*(q1p+q4p) + par__Rankle_y*(z[103]*q1p+z[103]*q4p+z[104]*q5p);
  z[221] = z[5]*z[202]*q5p + z[7]*z[201]*q5p + z[5]*z[100]*(q1p+q4p) + z[6]*
  z[101]*(q1p+q4p) + par__Rankle_y*(z[104]*q1p+z[104]*q4p+z[214]*q5p);
  z[222] = z[32]*z[207] + z[33]*z[209];
  z[223] = z[32]*z[103] + z[33]*z[210];
  z[224] = z[108]*q6p + z[222]*q1p + z[222]*q4p + z[223]*q5p;
  z[225] = z[32]*z[104] + z[33]*z[212];
  z[226] = z[109]*q6p + z[223]*q1p + z[223]*q4p + z[225]*q5p;
  z[227] = z[32]*z[209] - z[33]*z[207];
  z[228] = z[32]*z[210] - z[33]*z[103];
  z[229] = -z[32]*z[102] - z[33]*z[103];
  z[230] = z[227]*q1p + z[227]*q4p + z[228]*q5p + z[229]*q6p;
  z[231] = z[32]*z[212] - z[33]*z[104];
  z[232] = -z[19]*z[32] - z[33]*z[104];
  z[233] = z[228]*q1p + z[228]*q4p + z[231]*q5p + z[232]*q6p;
  z[234] = z[18]*z[33] - z[19]*z[32];
  z[235] = z[229]*q1p + z[229]*q4p + z[232]*q5p + z[234]*q6p;
  z[236] = z[32]*z[105] + z[33]*z[214];
  z[237] = z[107]*q1p + z[107]*q4p + z[112]*q6p + z[236]*q5p;
  z[238] = z[32]*z[214] - z[33]*z[105];
  z[239] = -z[21]*z[32] - z[33]*z[105];
  z[240] = z[109]*q1p + z[109]*q4p + z[238]*q5p + z[239]*q6p;
  z[241] = z[20]*z[33] - z[21]*z[32];
  z[242] = z[110]*q1p + z[110]*q4p + z[239]*q5p + z[241]*q6p;
  z[243] = (par__Lhip_x*z[1]-par__Lhip_y*z[2])*q1p;
  z[244] = (par__Lhip_x*z[2]+par__Lhip_y*z[1])*q1p;
  z[245] = par__Lhip_y*z[2]*q1p + z[53]*z[115]*q8p + z[114]*z[116]*q8p + 
  z[54]*z[60]*(q1p+q7p) + z[59]*z[114]*(q1p+q7p) - par__Lhip_x*z[1]*q1p;
  z[246] = z[53]*z[115]*q8p + z[114]*z[116]*q8p + z[54]*z[60]*(q1p+q7p) + 
  z[59]*z[114]*(q1p+q7p);
  z[247] = 2*par__Lknee_x3 + 6*par__Lknee_x2*q8 + 12*par__Lknee_x1*pow(q8,2);
  z[248] = 2*par__Lknee_y3 + 6*par__Lknee_y2*q8 + 12*par__Lknee_y1*pow(q8,2);
  z[249] = z[52]*z[247]*q8p + z[53]*z[248]*q8p + z[53]*z[115]*(q1p+q7p) + 
  z[114]*z[116]*(q1p+q7p);
  z[250] = z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p) - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[251] = z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p);
  z[252] = z[52]*z[248]*q8p + z[54]*z[247]*q8p + z[52]*z[115]*(q1p+q7p) + 
  z[53]*z[116]*(q1p+q7p);
  z[253] = z[54]*z[64] + z[63]*z[114];
  z[254] = z[118]*q8p + z[253]*q1p + z[253]*q7p;
  z[255] = z[54]*z[63] - z[64]*z[114];
  z[256] = -z[53]*z[63] - z[64]*z[114];
  z[257] = z[255]*q1p + z[255]*q7p + z[256]*q8p;
  z[258] = z[52]*z[64] - z[53]*z[63];
  z[259] = z[256]*q1p + z[256]*q7p + z[258]*q8p;
  z[260] = z[54]*z[64] - z[52]*z[63];
  z[261] = z[119]*q1p + z[119]*q7p + z[260]*q8p;
  z[262] = par__Lhip_y*z[2]*q1p + z[53]*z[115]*q8p + z[114]*z[116]*q8p + 
  z[54]*z[60]*(q1p+q7p) + z[59]*z[114]*(q1p+q7p) + par__Lankle_y*(z[255]*q1p+
  z[255]*q7p+z[256]*q8p) - par__Lhip_x*z[1]*q1p;
  z[263] = z[53]*z[115]*q8p + z[114]*z[116]*q8p + z[54]*z[60]*(q1p+q7p) + 
  z[59]*z[114]*(q1p+q7p) + par__Lankle_y*(z[255]*q1p+z[255]*q7p+z[256]*q8p);
  z[264] = z[52]*z[247]*q8p + z[53]*z[248]*q8p + z[53]*z[115]*(q1p+q7p) + 
  z[114]*z[116]*(q1p+q7p) + par__Lankle_y*(z[256]*q1p+z[256]*q7p+z[258]*q8p);
  z[265] = z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p) + par__Lankle_y*(z[118]*q1p+z[118]*q7p+z[119]*q8p) - 
  par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[266] = z[52]*z[115]*q8p + z[53]*z[116]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[114]*(q1p+q7p) + par__Lankle_y*(z[118]*q1p+z[118]*q7p+z[119]*q8p);
  z[267] = z[52]*z[248]*q8p + z[54]*z[247]*q8p + z[52]*z[115]*(q1p+q7p) + 
  z[53]*z[116]*(q1p+q7p) + par__Lankle_y*(z[119]*q1p+z[119]*q7p+z[260]*q8p);
  z[268] = z[79]*z[253] + z[80]*z[255];
  z[269] = z[79]*z[118] + z[80]*z[256];
  z[270] = z[123]*q9p + z[268]*q1p + z[268]*q7p + z[269]*q8p;
  z[271] = z[79]*z[119] + z[80]*z[258];
  z[272] = z[124]*q9p + z[269]*q1p + z[269]*q7p + z[271]*q8p;
  z[273] = z[79]*z[255] - z[80]*z[253];
  z[274] = z[79]*z[256] - z[80]*z[118];
  z[275] = -z[79]*z[117] - z[80]*z[118];
  z[276] = z[273]*q1p + z[273]*q7p + z[274]*q8p + z[275]*q9p;
  z[277] = z[79]*z[258] - z[80]*z[119];
  z[278] = -z[66]*z[79] - z[80]*z[119];
  z[279] = z[274]*q1p + z[274]*q7p + z[277]*q8p + z[278]*q9p;
  z[280] = z[65]*z[80] - z[66]*z[79];
  z[281] = z[275]*q1p + z[275]*q7p + z[278]*q8p + z[280]*q9p;
  z[282] = z[79]*z[120] + z[80]*z[260];
  z[283] = z[122]*q1p + z[122]*q7p + z[127]*q9p + z[282]*q8p;
  z[284] = z[79]*z[260] - z[80]*z[120];
  z[285] = -z[68]*z[79] - z[80]*z[120];
  z[286] = z[124]*q1p + z[124]*q7p + z[284]*q8p + z[285]*q9p;
  z[287] = z[67]*z[80] - z[68]*z[79];
  z[288] = z[125]*q1p + z[125]*q7p + z[285]*q8p + z[287]*q9p;
  z[289] = (par__back_x*z[1]-par__back_y*z[2])*q1p;
  z[290] = (par__back_x*z[2]+par__back_y*z[1])*q1p;
  z[291] = z[7]*(q1p+q4p);
  z[292] = z[54]*(q1p+q7p);
  z[293] = par__gravity_x*par__Lfemur_M;
  z[294] = par__gravity_y*par__Lfemur_M;
  z[295] = par__gravity_x*par__Lfoot_M;
  z[296] = par__gravity_y*par__Lfoot_M;
  z[297] = par__gravity_x*par__Ltibia_M;
  z[298] = par__gravity_y*par__Ltibia_M;
  z[299] = par__gravity_x*par__pelvis_M;
  z[300] = par__gravity_y*par__pelvis_M;
  z[301] = par__gravity_x*par__Rfemur_M;
  z[302] = par__gravity_y*par__Rfemur_M;
  z[303] = par__gravity_x*par__Rfoot_M;
  z[304] = par__gravity_y*par__Rfoot_M;
  z[305] = par__gravity_x*par__Rtibia_M;
  z[306] = par__gravity_y*par__Rtibia_M;
  z[307] = par__gravity_x*par__torso_M;
  z[308] = par__gravity_y*par__torso_M;
  sx = -par__wind - z[1]*(par__pelvis_CMy*q1p-z[1]*q2p) - z[2]*(par__pelvis_CMx*
  q1p-z[2]*q2p);
  sy = z[1]*(par__pelvis_CMx*q1p+z[1]*q3p) - z[2]*(par__pelvis_CMy*q1p-z[2]*
  q3p);
  s = pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[309] = par__airdrag*par__wind;
  z[310] = z[309]*s;
  z[311] = par__airdrag*s*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*q3p);
  z[312] = par__airdrag*s*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*q2p);
  z[313] = z[299] + z[310];
  z[321] = par__Rfoot_CMy + z[44];
  z[322] = par__Rfoot_CMy + z[45];
  z[323] = par__Rfoot_CMy + z[46];
  z[324] = z[47] - par__Rfoot_CMx;
  z[325] = z[48] - par__Rfoot_CMx;
  z[326] = z[49] - par__Rfoot_CMx;
  z[327] = par__bodyweight*G1;
  z[328] = par__bodyweight*G2;
  z[330] = par__bodyweight*G3;
  z[338] = par__Lfoot_CMy + z[91];
  z[339] = par__Lfoot_CMy + z[92];
  z[340] = par__Lfoot_CMy + z[93];
  z[341] = z[94] - par__Lfoot_CMx;
  z[342] = z[95] - par__Lfoot_CMx;
  z[343] = z[96] - par__Lfoot_CMx;
  z[344] = par__bodyweight*G4;
  z[345] = par__bodyweight*G5;
  z[348] = par__bodyweight*G6;
  z[350] = par__pelvis_CMx*z[300];
  z[351] = z[97]*z[307];
  z[352] = z[97]*z[308];
  z[353] = z[98]*z[308];
  z[354] = par__pelvis_CMy*z[300];
  z[355] = z[98]*z[307];
  z[356] = z[330] + z[348] + z[350]*z[1] + z[351]*z[1] + z[352]*z[2] + z[353]*
  z[1] + q3*z[327] + q3*z[344] + par__Lankle_y*z[65]*z[344] + par__Lankle_y*
  z[67]*z[345] + par__Lhip_x*z[2]*z[344] + par__Lhip_y*z[1]*z[344] + par__Lhip_y*
  z[2]*z[345] + par__Rankle_y*z[18]*z[327] + par__Rankle_y*z[20]*z[328] + 
  par__Rhip_x*z[2]*z[327] + par__Rhip_y*z[1]*z[327] + par__Rhip_y*z[2]*z[328] + 
  z[293]*z[52]*z[57] + z[293]*z[53]*z[58] + z[294]*z[52]*z[58] + z[294]*z[54]*
  z[57] + z[295]*z[81]*z[92] + z[295]*z[82]*z[95] + z[296]*z[83]*z[92] + 
  z[296]*z[84]*z[95] + z[297]*z[65]*z[74] + z[297]*z[66]*z[76] + z[298]*z[67]*
  z[74] + z[298]*z[68]*z[76] + z[301]*z[5]*z[10] + z[301]*z[6]*z[11] + z[302]*
  z[5]*z[11] + z[302]*z[7]*z[10] + z[303]*z[34]*z[45] + z[303]*z[35]*z[48] + 
  z[304]*z[36]*z[45] + z[304]*z[37]*z[48] + z[305]*z[18]*z[27] + z[305]*z[19]*
  z[29] + z[306]*z[20]*z[27] + z[306]*z[21]*z[29] + z[5]*z[13]*z[327] + z[7]*
  z[13]*z[328] + z[34]*z[322]*z[327] + z[35]*z[325]*z[327] + z[36]*z[322]*
  z[328] + z[37]*z[325]*z[328] + z[52]*z[60]*z[344] + z[54]*z[60]*z[345] + 
  z[81]*z[339]*z[344] + z[82]*z[342]*z[344] + z[83]*z[339]*z[345] + z[84]*
  z[342]*z[345] - z[354]*z[2] - z[355]*z[2] - q2*z[328] - q2*z[345] - par__Lhip_x*
  z[1]*z[345] - par__pelvis_CMx*z[2]*z[313] - par__pelvis_CMy*z[1]*z[313] - 
  par__Rhip_x*z[1]*z[328] - z[5]*z[12]*z[328] - z[6]*z[12]*z[327] - z[52]*
  z[59]*z[345] - z[53]*z[59]*z[344] - par__pelvis_CMx*z[312] - par__pelvis_CMy*
  z[311];
  z[357] = z[293]*pow(z[52],2) + z[293]*pow(z[53],2) + z[294]*z[52]*z[53] + 
  z[294]*z[52]*z[54] + z[295]*pow(z[81],2) + z[295]*pow(z[82],2) + z[296]*
  z[81]*z[83] + z[296]*z[82]*z[84] + z[297]*pow(z[65],2) + z[297]*pow(z[66],2) + 
  z[298]*z[65]*z[67] + z[298]*z[66]*z[68] + z[301]*pow(z[5],2) + z[301]*pow(
  z[6],2) + z[302]*z[5]*z[6] + z[302]*z[5]*z[7] + z[303]*pow(z[34],2) + 
  z[303]*pow(z[35],2) + z[304]*z[34]*z[36] + z[304]*z[35]*z[37] + z[305]*pow(
  z[18],2) + z[305]*pow(z[19],2) + z[306]*z[18]*z[20] + z[306]*z[19]*z[21] + 
  z[34]*z[36]*z[328] + z[35]*z[37]*z[328] + z[81]*z[83]*z[345] + z[82]*z[84]*
  z[345] + z[313]*pow(z[1],2) + z[313]*pow(z[2],2) + z[327]*pow(z[34],2) + 
  z[327]*pow(z[35],2) + z[344]*pow(z[81],2) + z[344]*pow(z[82],2) + z[307]*(
  pow(z[1],2)+pow(z[2],2)) + z[1]*z[311] + z[2]*z[312];
  z[358] = z[293]*z[52]*z[53] + z[293]*z[52]*z[54] + z[294]*pow(z[52],2) + 
  z[294]*pow(z[54],2) + z[295]*z[81]*z[83] + z[295]*z[82]*z[84] + z[296]*pow(
  z[83],2) + z[296]*pow(z[84],2) + z[297]*z[65]*z[67] + z[297]*z[66]*z[68] + 
  z[298]*pow(z[67],2) + z[298]*pow(z[68],2) + z[300]*pow(z[1],2) + z[300]*
  pow(z[2],2) + z[301]*z[5]*z[6] + z[301]*z[5]*z[7] + z[302]*pow(z[5],2) + 
  z[302]*pow(z[7],2) + z[303]*z[34]*z[36] + z[303]*z[35]*z[37] + z[304]*pow(
  z[36],2) + z[304]*pow(z[37],2) + z[305]*z[18]*z[20] + z[305]*z[19]*z[21] + 
  z[306]*pow(z[20],2) + z[306]*pow(z[21],2) + z[34]*z[36]*z[327] + z[35]*
  z[37]*z[327] + z[81]*z[83]*z[344] + z[82]*z[84]*z[344] + z[328]*pow(z[36],2) + 
  z[328]*pow(z[37],2) + z[345]*pow(z[83],2) + z[345]*pow(z[84],2) + z[308]*(
  pow(z[1],2)+pow(z[2],2)) + z[2]*z[311] - z[1]*z[312];
  z[359] = par__Rfemur_CMx*z[301];
  z[360] = par__Rfemur_CMx*z[302];
  z[361] = par__Rfemur_CMy*z[301];
  z[362] = par__Rfemur_CMy*z[302];
  z[363] = z[330] + z[359]*z[6] + z[360]*z[5] + q3*z[327] + par__Rankle_y*
  z[18]*z[327] + par__Rankle_y*z[20]*z[328] + par__Rhip_x*z[2]*z[327] + 
  par__Rhip_y*z[1]*z[327] + par__Rhip_y*z[2]*z[328] + z[303]*z[34]*z[44] + 
  z[303]*z[35]*z[47] + z[304]*z[36]*z[44] + z[304]*z[37]*z[47] + z[305]*z[18]*
  z[26] + z[305]*z[19]*z[28] + z[306]*z[20]*z[26] + z[306]*z[21]*z[28] + z[5]*
  z[13]*z[327] + z[7]*z[13]*z[328] + z[34]*z[321]*z[327] + z[35]*z[324]*
  z[327] + z[36]*z[321]*z[328] + z[37]*z[324]*z[328] - z[361]*z[5] - z[362]*
  z[7] - q2*z[328] - par__Rhip_x*z[1]*z[328] - z[5]*z[12]*z[328] - z[6]*z[12]*
  z[327];
  z[364] = par__Rtibia_CMx*z[305];
  z[365] = par__Rtibia_CMx*z[306];
  z[366] = par__Rtibia_CMy*z[305];
  z[367] = par__Rtibia_CMy*z[306];
  z[368] = z[330] + z[364]*z[19] + z[365]*z[21] + q3*z[327] + par__Rankle_y*
  z[18]*z[327] + par__Rankle_y*z[20]*z[328] + par__Rhip_x*z[2]*z[327] + 
  par__Rhip_y*z[1]*z[327] + par__Rhip_y*z[2]*z[328] + z[303]*z[34]*z[46] + 
  z[303]*z[35]*z[49] + z[304]*z[36]*z[46] + z[304]*z[37]*z[49] + z[5]*z[13]*
  z[327] + z[7]*z[13]*z[328] + z[34]*z[323]*z[327] + z[35]*z[326]*z[327] + 
  z[36]*z[323]*z[328] + z[37]*z[326]*z[328] - z[366]*z[18] - z[367]*z[20] - 
  q2*z[328] - par__Rhip_x*z[1]*z[328] - z[5]*z[12]*z[328] - z[6]*z[12]*z[327];
  z[369] = par__Rfoot_CMx*z[303];
  z[370] = par__Rfoot_CMx*z[304];
  z[371] = par__Rfoot_CMy*z[303];
  z[372] = par__Rfoot_CMy*z[304];
  z[373] = z[330] + z[369]*z[35] + z[370]*z[37] + q3*z[327] + par__Rankle_y*
  z[18]*z[327] + par__Rankle_y*z[20]*z[328] + par__Rhip_x*z[2]*z[327] + 
  par__Rhip_y*z[1]*z[327] + par__Rhip_y*z[2]*z[328] + z[5]*z[13]*z[327] + 
  z[7]*z[13]*z[328] - z[371]*z[34] - z[372]*z[36] - q2*z[328] - par__Rhip_x*
  z[1]*z[328] - z[5]*z[12]*z[328] - z[6]*z[12]*z[327];
  z[374] = par__Lfemur_CMx*z[293];
  z[375] = par__Lfemur_CMx*z[294];
  z[376] = par__Lfemur_CMy*z[293];
  z[377] = par__Lfemur_CMy*z[294];
  z[378] = z[348] + z[374]*z[53] + z[375]*z[52] + q3*z[344] + par__Lankle_y*
  z[65]*z[344] + par__Lankle_y*z[67]*z[345] + par__Lhip_x*z[2]*z[344] + 
  par__Lhip_y*z[1]*z[344] + par__Lhip_y*z[2]*z[345] + z[295]*z[81]*z[91] + 
  z[295]*z[82]*z[94] + z[296]*z[83]*z[91] + z[296]*z[84]*z[94] + z[297]*z[65]*
  z[73] + z[297]*z[66]*z[75] + z[298]*z[67]*z[73] + z[298]*z[68]*z[75] + 
  z[52]*z[60]*z[344] + z[54]*z[60]*z[345] + z[81]*z[338]*z[344] + z[82]*
  z[341]*z[344] + z[83]*z[338]*z[345] + z[84]*z[341]*z[345] - z[376]*z[52] - 
  z[377]*z[54] - q2*z[345] - par__Lhip_x*z[1]*z[345] - z[52]*z[59]*z[345] - 
  z[53]*z[59]*z[344];
  z[379] = par__Ltibia_CMx*z[297];
  z[380] = par__Ltibia_CMx*z[298];
  z[381] = par__Ltibia_CMy*z[297];
  z[382] = par__Ltibia_CMy*z[298];
  z[383] = z[348] + z[379]*z[66] + z[380]*z[68] + q3*z[344] + par__Lankle_y*
  z[65]*z[344] + par__Lankle_y*z[67]*z[345] + par__Lhip_x*z[2]*z[344] + 
  par__Lhip_y*z[1]*z[344] + par__Lhip_y*z[2]*z[345] + z[295]*z[81]*z[93] + 
  z[295]*z[82]*z[96] + z[296]*z[83]*z[93] + z[296]*z[84]*z[96] + z[52]*z[60]*
  z[344] + z[54]*z[60]*z[345] + z[81]*z[340]*z[344] + z[82]*z[343]*z[344] + 
  z[83]*z[340]*z[345] + z[84]*z[343]*z[345] - z[381]*z[65] - z[382]*z[67] - 
  q2*z[345] - par__Lhip_x*z[1]*z[345] - z[52]*z[59]*z[345] - z[53]*z[59]*
  z[344];
  z[384] = par__Lfoot_CMx*z[295];
  z[385] = par__Lfoot_CMx*z[296];
  z[386] = par__Lfoot_CMy*z[295];
  z[387] = par__Lfoot_CMy*z[296];
  z[388] = z[348] + z[384]*z[82] + z[385]*z[84] + q3*z[344] + par__Lankle_y*
  z[65]*z[344] + par__Lankle_y*z[67]*z[345] + par__Lhip_x*z[2]*z[344] + 
  par__Lhip_y*z[1]*z[344] + par__Lhip_y*z[2]*z[345] + z[52]*z[60]*z[344] + 
  z[54]*z[60]*z[345] - z[386]*z[81] - z[387]*z[83] - q2*z[345] - par__Lhip_x*
  z[1]*z[345] - z[52]*z[59]*z[345] - z[53]*z[59]*z[344];
  z[389] = q1p + q7p;
  z[390] = q1p + q7p + q8p + q9p;
  z[391] = q1p + q7p + q8p;
  z[392] = q1p + q4p;
  z[393] = q1p + q4p + q5p + q6p;
  z[394] = q1p + q4p + q5p;
  z[395] = z[52]*q2p + z[54]*q3p + z[57]*q1p - par__Lfemur_CMy*q7p;
  z[396] = par__Lfemur_CMx*q7p + z[52]*q3p + z[53]*q2p + z[58]*q1p;
  z[397] = -z[1]*z[51]*q1p - z[1]*z[51]*q7p - z[2]*z[50]*q1p - z[2]*z[50]*q7p;
  z[398] = z[1]*z[50]*q1p + z[1]*z[50]*q7p - z[2]*z[51]*q1p - z[2]*z[51]*q7p;
  z[399] = (par__Lhip_x*z[50]+par__Lhip_y*z[51])*q7p;
  z[400] = q1p*z[399] + q2p*z[397] + q3p*z[398];
  z[401] = z[2]*z[51]*q1p + z[2]*z[51]*q7p - z[1]*z[50]*q1p - z[1]*z[50]*q7p;
  z[402] = (par__Lhip_x*z[51]-par__Lhip_y*z[50])*q7p;
  z[403] = q2p*z[401] + q3p*z[397] - q1p*z[402];
  z[404] = z[400] - z[389]*z[396];
  z[405] = z[403] + z[389]*z[395];
  z[406] = z[81]*q2p + z[83]*q3p + z[91]*q7p + z[92]*q1p + z[93]*q8p - par__Lfoot_CMy*
  q9p;
  z[407] = par__Lfoot_CMx*q9p + z[82]*q2p + z[84]*q3p + z[94]*q7p + z[95]*q1p + 
  z[96]*q8p;
  z[408] = z[53]*z[63]*q8p + z[63]*z[397] + z[64]*z[401] - z[52]*z[64]*q8p;
  z[409] = z[63]*z[401] - z[52]*z[63]*q8p - z[53]*z[64]*q8p - z[64]*z[397];
  z[410] = z[66]*z[79]*q9p + z[79]*z[408] + z[80]*z[409] - z[65]*z[80]*q9p;
  z[411] = z[52]*z[63]*q8p + z[63]*z[398] + z[64]*z[397] - z[54]*z[64]*q8p;
  z[412] = z[63]*z[397] - z[52]*z[64]*q8p - z[54]*z[63]*q8p - z[64]*z[398];
  z[413] = z[68]*z[79]*q9p + z[79]*z[411] + z[80]*z[412] - z[67]*z[80]*q9p;
  z[414] = (z[59]*z[64]-z[60]*z[63]-z[63]*z[115]-z[64]*z[116])*q8p;
  z[415] = (z[63]*z[116]-z[59]*z[63]-z[60]*z[64]-z[64]*z[115])*q8p;
  z[416] = z[71]*z[79]*q9p - z[77]*z[80]*q9p - z[79]*z[415] - z[80]*z[414];
  z[417] = z[115]*q8p - z[402];
  z[418] = z[399] - z[116]*q8p;
  z[419] = z[63]*z[417] - z[61]*z[63]*q8p - z[62]*z[64]*q8p - z[64]*z[418];
  z[420] = z[62]*z[63]*q8p + z[63]*z[418] + z[64]*z[417] - z[61]*z[64]*q8p;
  z[421] = z[72]*z[79]*q9p + z[79]*z[420] + z[80]*z[419] - z[78]*z[80]*q9p;
  z[422] = par__Lankle_y*z[80]*q9p;
  z[423] = q1p*z[421] + q2p*z[410] + q3p*z[413] + q7p*z[416] + q8p*z[422];
  z[424] = z[79]*z[409] - z[65]*z[79]*q9p - z[66]*z[80]*q9p - z[80]*z[408];
  z[425] = z[79]*z[412] - z[67]*z[79]*q9p - z[68]*z[80]*q9p - z[80]*z[411];
  z[426] = z[80]*z[415] - z[71]*z[80]*q9p - z[77]*z[79]*q9p - z[79]*z[414];
  z[427] = z[79]*z[419] - z[72]*z[80]*q9p - z[78]*z[79]*q9p - z[80]*z[420];
  z[428] = par__Lankle_y*z[79]*q9p;
  z[429] = q1p*z[427] + q2p*z[424] + q3p*z[425] + q7p*z[426] + q8p*z[428];
  z[430] = z[423] - z[390]*z[407];
  z[431] = z[429] + z[390]*z[406];
  z[432] = z[65]*q2p + z[67]*q3p + z[73]*q7p + z[74]*q1p - par__Ltibia_CMy*
  q8p;
  z[433] = par__Ltibia_CMx*q8p + z[66]*q2p + z[68]*q3p + z[75]*q7p + z[76]*
  q1p;
  z[434] = q1p*z[420] + q2p*z[408] + q3p*z[411] - q7p*z[415];
  z[435] = q1p*z[419] + q2p*z[409] + q3p*z[412] - q7p*z[414];
  z[436] = z[434] - z[391]*z[433];
  z[437] = z[435] + z[391]*z[432];
  z[438] = z[1]*q2p + z[2]*q3p - par__pelvis_CMy*q1p;
  z[439] = par__pelvis_CMx*q1p + z[1]*q3p - z[2]*q2p;
  z[440] = q1p*(z[1]*q3p-z[2]*q2p);
  z[441] = q1p*(z[1]*q2p+z[2]*q3p);
  z[442] = z[440] - q1p*z[439];
  z[443] = q1p*z[438] - z[441];
  z[444] = z[5]*q2p + z[7]*q3p + z[10]*q1p - par__Rfemur_CMy*q4p;
  z[445] = par__Rfemur_CMx*q4p + z[5]*q3p + z[6]*q2p + z[11]*q1p;
  z[446] = -z[1]*z[4]*q1p - z[1]*z[4]*q4p - z[2]*z[3]*q1p - z[2]*z[3]*q4p;
  z[447] = z[1]*z[3]*q1p + z[1]*z[3]*q4p - z[2]*z[4]*q1p - z[2]*z[4]*q4p;
  z[448] = (par__Rhip_x*z[3]+par__Rhip_y*z[4])*q4p;
  z[449] = q1p*z[448] + q2p*z[446] + q3p*z[447];
  z[450] = z[2]*z[4]*q1p + z[2]*z[4]*q4p - z[1]*z[3]*q1p - z[1]*z[3]*q4p;
  z[451] = (par__Rhip_x*z[4]-par__Rhip_y*z[3])*q4p;
  z[452] = q2p*z[450] + q3p*z[446] - q1p*z[451];
  z[453] = z[449] - z[392]*z[445];
  z[454] = z[452] + z[392]*z[444];
  z[455] = z[34]*q2p + z[36]*q3p + z[44]*q4p + z[45]*q1p + z[46]*q5p - par__Rfoot_CMy*
  q6p;
  z[456] = par__Rfoot_CMx*q6p + z[35]*q2p + z[37]*q3p + z[47]*q4p + z[48]*q1p + 
  z[49]*q5p;
  z[457] = z[6]*z[16]*q5p + z[16]*z[446] + z[17]*z[450] - z[5]*z[17]*q5p;
  z[458] = z[16]*z[450] - z[5]*z[16]*q5p - z[6]*z[17]*q5p - z[17]*z[446];
  z[459] = z[19]*z[32]*q6p + z[32]*z[457] + z[33]*z[458] - z[18]*z[33]*q6p;
  z[460] = z[5]*z[16]*q5p + z[16]*z[447] + z[17]*z[446] - z[7]*z[17]*q5p;
  z[461] = z[16]*z[446] - z[5]*z[17]*q5p - z[7]*z[16]*q5p - z[17]*z[447];
  z[462] = z[21]*z[32]*q6p + z[32]*z[460] + z[33]*z[461] - z[20]*z[33]*q6p;
  z[463] = (z[12]*z[17]-z[13]*z[16]-z[16]*z[100]-z[17]*z[101])*q5p;
  z[464] = (z[16]*z[101]-z[12]*z[16]-z[13]*z[17]-z[17]*z[100])*q5p;
  z[465] = z[24]*z[32]*q6p - z[30]*z[33]*q6p - z[32]*z[464] - z[33]*z[463];
  z[466] = z[100]*q5p - z[451];
  z[467] = z[448] - z[101]*q5p;
  z[468] = z[16]*z[466] - z[14]*z[16]*q5p - z[15]*z[17]*q5p - z[17]*z[467];
  z[469] = z[15]*z[16]*q5p + z[16]*z[467] + z[17]*z[466] - z[14]*z[17]*q5p;
  z[470] = z[25]*z[32]*q6p + z[32]*z[469] + z[33]*z[468] - z[31]*z[33]*q6p;
  z[471] = par__Rankle_y*z[33]*q6p;
  z[472] = q1p*z[470] + q2p*z[459] + q3p*z[462] + q4p*z[465] + q5p*z[471];
  z[473] = z[32]*z[458] - z[18]*z[32]*q6p - z[19]*z[33]*q6p - z[33]*z[457];
  z[474] = z[32]*z[461] - z[20]*z[32]*q6p - z[21]*z[33]*q6p - z[33]*z[460];
  z[475] = z[33]*z[464] - z[24]*z[33]*q6p - z[30]*z[32]*q6p - z[32]*z[463];
  z[476] = z[32]*z[468] - z[25]*z[33]*q6p - z[31]*z[32]*q6p - z[33]*z[469];
  z[477] = par__Rankle_y*z[32]*q6p;
  z[478] = q1p*z[476] + q2p*z[473] + q3p*z[474] + q4p*z[475] + q5p*z[477];
  z[479] = z[472] - z[393]*z[456];
  z[480] = z[478] + z[393]*z[455];
  z[481] = z[18]*q2p + z[20]*q3p + z[26]*q4p + z[27]*q1p - par__Rtibia_CMy*
  q5p;
  z[482] = par__Rtibia_CMx*q5p + z[19]*q2p + z[21]*q3p + z[28]*q4p + z[29]*
  q1p;
  z[483] = q1p*z[469] + q2p*z[457] + q3p*z[460] - q4p*z[464];
  z[484] = q1p*z[468] + q2p*z[458] + q3p*z[461] - q4p*z[463];
  z[485] = z[483] - z[394]*z[482];
  z[486] = z[484] + z[394]*z[481];
  z[487] = z[97]*q1p + z[1]*q2p + z[2]*q3p;
  z[488] = z[98]*q1p + z[1]*q3p - z[2]*q2p;
  z[489] = z[440] - q1p*z[488];
  z[490] = q1p*z[487] - z[441];
  z[491] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz + par__pelvis_Izz + 
  par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz + par__torso_Izz + 
  par__pelvis_M*(pow(par__pelvis_CMx,2)+pow(par__pelvis_CMy,2)) + par__torso_M*(
  pow(z[97],2)+pow(z[98],2));
  z[492] = z[491] + par__Lfemur_M*(pow(z[57],2)+pow(z[58],2)) + par__Lfoot_M*(
  pow(z[92],2)+pow(z[95],2)) + par__Ltibia_M*(pow(z[74],2)+pow(z[76],2)) + 
  par__Rfemur_M*(pow(z[10],2)+pow(z[11],2)) + par__Rfoot_M*(pow(z[45],2)+pow(
  z[48],2)) + par__Rtibia_M*(pow(z[27],2)+pow(z[29],2));
  z[493] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz;
  z[494] = z[493] + par__Lfoot_M*(z[91]*z[92]+z[94]*z[95]) + par__Ltibia_M*(
  z[73]*z[74]+z[75]*z[76]) + par__Lfemur_M*(par__Lfemur_CMx*z[58]-par__Lfemur_CMy*
  z[57]);
  z[495] = par__Lfoot_Izz + par__Ltibia_Izz;
  z[496] = z[495] + par__Lfoot_M*(z[92]*z[93]+z[95]*z[96]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[76]-par__Ltibia_CMy*z[74]);
  z[497] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[95]-par__Lfoot_CMy*
  z[92]);
  z[498] = par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz;
  z[499] = z[498] + par__Rfoot_M*(z[44]*z[45]+z[47]*z[48]) + par__Rtibia_M*(
  z[26]*z[27]+z[28]*z[29]) + par__Rfemur_M*(par__Rfemur_CMx*z[11]-par__Rfemur_CMy*
  z[10]);
  z[500] = par__Rfoot_Izz + par__Rtibia_Izz;
  z[501] = z[500] + par__Rfoot_M*(z[45]*z[46]+z[48]*z[49]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[29]-par__Rtibia_CMy*z[27]);
  z[502] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[48]-par__Rfoot_CMy*
  z[45]);
  z[503] = par__Lfemur_M*(z[52]*z[58]+z[54]*z[57]) + par__Lfoot_M*(z[83]*
  z[92]+z[84]*z[95]) + par__Ltibia_M*(z[67]*z[74]+z[68]*z[76]) + par__Rfemur_M*(
  z[5]*z[11]+z[7]*z[10]) + par__Rfoot_M*(z[36]*z[45]+z[37]*z[48]) + par__Rtibia_M*(
  z[20]*z[27]+z[21]*z[29]) + par__torso_M*(z[97]*z[2]+z[98]*z[1]) + par__pelvis_M*(
  par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2]);
  z[504] = par__Lfemur_M*(z[52]*z[57]+z[53]*z[58]) + par__Lfoot_M*(z[81]*
  z[92]+z[82]*z[95]) + par__Ltibia_M*(z[65]*z[74]+z[66]*z[76]) + par__Rfemur_M*(
  z[5]*z[10]+z[6]*z[11]) + par__Rfoot_M*(z[34]*z[45]+z[35]*z[48]) + par__Rtibia_M*(
  z[18]*z[27]+z[19]*z[29]) + par__torso_M*(z[97]*z[1]-z[98]*z[2]) - par__pelvis_M*(
  par__pelvis_CMx*z[2]+par__pelvis_CMy*z[1]);
  z[505] = par__Lfemur_M*(z[57]*z[404]+z[58]*z[405]) + par__Lfoot_M*(z[92]*
  z[430]+z[95]*z[431]) + par__Ltibia_M*(z[74]*z[436]+z[76]*z[437]) + par__Rfemur_M*(
  z[10]*z[453]+z[11]*z[454]) + par__Rfoot_M*(z[45]*z[479]+z[48]*z[480]) + 
  par__Rtibia_M*(z[27]*z[485]+z[29]*z[486]) + par__torso_M*(z[97]*z[489]+
  z[98]*z[490]) + par__pelvis_M*(par__pelvis_CMx*z[443]-par__pelvis_CMy*
  z[442]);
  z[506] = par__Lfoot_M*(z[81]*z[91]+z[82]*z[94]) + par__Ltibia_M*(z[65]*
  z[73]+z[66]*z[75]) + par__Lfemur_M*(par__Lfemur_CMx*z[53]-par__Lfemur_CMy*
  z[52]);
  z[507] = par__Lfemur_M*z[52]*(z[53]+z[54]) + par__Rfemur_M*z[5]*(z[6]+z[7]) + 
  par__Lfoot_M*(z[81]*z[83]+z[82]*z[84]) + par__Ltibia_M*(z[65]*z[67]+z[66]*
  z[68]) + par__Rfoot_M*(z[34]*z[36]+z[35]*z[37]) + par__Rtibia_M*(z[18]*
  z[20]+z[19]*z[21]);
  z[508] = par__Lfemur_M*(pow(z[52],2)+pow(z[53],2)) + par__Lfoot_M*(pow(
  z[81],2)+pow(z[82],2)) + par__Ltibia_M*(pow(z[65],2)+pow(z[66],2)) + par__pelvis_M*(
  pow(z[1],2)+pow(z[2],2)) + par__Rfemur_M*(pow(z[5],2)+pow(z[6],2)) + par__Rfoot_M*(
  pow(z[34],2)+pow(z[35],2)) + par__Rtibia_M*(pow(z[18],2)+pow(z[19],2)) + 
  par__torso_M*(pow(z[1],2)+pow(z[2],2));
  z[509] = par__Lfoot_M*(par__Lfoot_CMx*z[82]-par__Lfoot_CMy*z[81]);
  z[510] = par__Lfoot_M*(z[81]*z[93]+z[82]*z[96]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[66]-par__Ltibia_CMy*z[65]);
  z[511] = par__Rfoot_M*(z[34]*z[44]+z[35]*z[47]) + par__Rtibia_M*(z[18]*
  z[26]+z[19]*z[28]) + par__Rfemur_M*(par__Rfemur_CMx*z[6]-par__Rfemur_CMy*
  z[5]);
  z[512] = par__Rfoot_M*(par__Rfoot_CMx*z[35]-par__Rfoot_CMy*z[34]);
  z[513] = par__Rfoot_M*(z[34]*z[46]+z[35]*z[49]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[19]-par__Rtibia_CMy*z[18]);
  z[514] = par__Lfemur_M*(z[52]*z[404]+z[53]*z[405]) + par__Lfoot_M*(z[81]*
  z[430]+z[82]*z[431]) + par__Ltibia_M*(z[65]*z[436]+z[66]*z[437]) + par__Rfemur_M*(
  z[5]*z[453]+z[6]*z[454]) + par__Rfoot_M*(z[34]*z[479]+z[35]*z[480]) + 
  par__Rtibia_M*(z[18]*z[485]+z[19]*z[486]) + par__pelvis_M*(z[1]*z[442]-z[2]*
  z[443]) + par__torso_M*(z[1]*z[489]-z[2]*z[490]);
  z[515] = par__Lfoot_M*(z[83]*z[91]+z[84]*z[94]) + par__Ltibia_M*(z[67]*
  z[73]+z[68]*z[75]) + par__Lfemur_M*(par__Lfemur_CMx*z[52]-par__Lfemur_CMy*
  z[54]);
  z[516] = par__Lfemur_M*(pow(z[52],2)+pow(z[54],2)) + par__Lfoot_M*(pow(
  z[83],2)+pow(z[84],2)) + par__Ltibia_M*(pow(z[67],2)+pow(z[68],2)) + par__pelvis_M*(
  pow(z[1],2)+pow(z[2],2)) + par__Rfemur_M*(pow(z[5],2)+pow(z[7],2)) + par__Rfoot_M*(
  pow(z[36],2)+pow(z[37],2)) + par__Rtibia_M*(pow(z[20],2)+pow(z[21],2)) + 
  par__torso_M*(pow(z[1],2)+pow(z[2],2));
  z[517] = par__Lfoot_M*(par__Lfoot_CMx*z[84]-par__Lfoot_CMy*z[83]);
  z[518] = par__Lfoot_M*(z[83]*z[93]+z[84]*z[96]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[68]-par__Ltibia_CMy*z[67]);
  z[519] = par__Rfoot_M*(z[36]*z[44]+z[37]*z[47]) + par__Rtibia_M*(z[20]*
  z[26]+z[21]*z[28]) + par__Rfemur_M*(par__Rfemur_CMx*z[5]-par__Rfemur_CMy*
  z[7]);
  z[520] = par__Rfoot_M*(par__Rfoot_CMx*z[37]-par__Rfoot_CMy*z[36]);
  z[521] = par__Rfoot_M*(z[36]*z[46]+z[37]*z[49]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[21]-par__Rtibia_CMy*z[20]);
  z[522] = par__Lfemur_M*(z[52]*z[405]+z[54]*z[404]) + par__Lfoot_M*(z[83]*
  z[430]+z[84]*z[431]) + par__Ltibia_M*(z[67]*z[436]+z[68]*z[437]) + par__pelvis_M*(
  z[1]*z[443]+z[2]*z[442]) + par__Rfemur_M*(z[5]*z[454]+z[7]*z[453]) + par__Rfoot_M*(
  z[36]*z[479]+z[37]*z[480]) + par__Rtibia_M*(z[20]*z[485]+z[21]*z[486]) + 
  par__torso_M*(z[1]*z[490]+z[2]*z[489]);
  z[523] = par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz + par__Rfemur_M*(
  pow(par__Rfemur_CMx,2)+pow(par__Rfemur_CMy,2));
  z[524] = z[523] + par__Rfoot_M*(pow(z[44],2)+pow(z[47],2)) + par__Rtibia_M*(
  pow(z[26],2)+pow(z[28],2));
  z[525] = z[500] + par__Rfoot_M*(z[44]*z[46]+z[47]*z[49]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[28]-par__Rtibia_CMy*z[26]);
  z[526] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[47]-par__Rfoot_CMy*
  z[44]);
  z[527] = par__Rfoot_M*(z[44]*z[479]+z[47]*z[480]) + par__Rtibia_M*(z[26]*
  z[485]+z[28]*z[486]) + par__Rfemur_M*(par__Rfemur_CMx*z[454]-par__Rfemur_CMy*
  z[453]);
  z[528] = par__Rfoot_Izz + par__Rtibia_Izz + par__Rtibia_M*(pow(par__Rtibia_CMx,
  2)+pow(par__Rtibia_CMy,2));
  z[529] = z[528] + par__Rfoot_M*(pow(z[46],2)+pow(z[49],2));
  z[530] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[49]-par__Rfoot_CMy*
  z[46]);
  z[531] = par__Rfoot_M*(z[46]*z[479]+z[49]*z[480]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[486]-par__Rtibia_CMy*z[485]);
  z[532] = par__Rfoot_Izz + par__Rfoot_M*(pow(par__Rfoot_CMx,2)+pow(par__Rfoot_CMy,
  2));
  z[533] = par__Rfoot_M*(par__Rfoot_CMx*z[480]-par__Rfoot_CMy*z[479]);
  z[534] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz + par__Lfemur_M*(
  pow(par__Lfemur_CMx,2)+pow(par__Lfemur_CMy,2));
  z[535] = z[534] + par__Lfoot_M*(pow(z[91],2)+pow(z[94],2)) + par__Ltibia_M*(
  pow(z[73],2)+pow(z[75],2));
  z[536] = z[495] + par__Lfoot_M*(z[91]*z[93]+z[94]*z[96]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[75]-par__Ltibia_CMy*z[73]);
  z[537] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[94]-par__Lfoot_CMy*
  z[91]);
  z[538] = par__Lfoot_M*(z[91]*z[430]+z[94]*z[431]) + par__Ltibia_M*(z[73]*
  z[436]+z[75]*z[437]) + par__Lfemur_M*(par__Lfemur_CMx*z[405]-par__Lfemur_CMy*
  z[404]);
  z[539] = par__Lfoot_Izz + par__Ltibia_Izz + par__Ltibia_M*(pow(par__Ltibia_CMx,
  2)+pow(par__Ltibia_CMy,2));
  z[540] = z[539] + par__Lfoot_M*(pow(z[93],2)+pow(z[96],2));
  z[541] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[96]-par__Lfoot_CMy*
  z[93]);
  z[542] = par__Lfoot_M*(z[93]*z[430]+z[96]*z[431]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[437]-par__Ltibia_CMy*z[436]);
  z[543] = par__Lfoot_Izz + par__Lfoot_M*(pow(par__Lfoot_CMx,2)+pow(par__Lfoot_CMy,
  2));
  z[544] = par__Lfoot_M*(par__Lfoot_CMx*z[431]-par__Lfoot_CMy*z[430]);


/* Evaluate output quantities */
}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[] )
{
int      i1;

/* Write output to screen and to output file(s) */
  fk[0] = q2;
  fk[1] = q3;
  fk[2] = z[1];
  fk[3] = -z[2];
  fk[4] = z[2];
  fk[5] = z[1];
  fk[6] = q2 + par__Rhip_x*z[1] - par__Rhip_y*z[2];
  fk[7] = q3 + par__Rhip_x*z[2] + par__Rhip_y*z[1];
  fk[8] = z[5];
  fk[9] = z[6];
  fk[10] = z[7];
  fk[11] = z[5];
  fk[12] = q2 + par__Rhip_x*z[1] + z[5]*z[12] + z[6]*z[13] - par__Rhip_y*z[2];
  fk[13] = q3 + par__Rhip_x*z[2] + par__Rhip_y*z[1] + z[5]*z[13] + z[7]*z[12];
  fk[14] = z[18];
  fk[15] = z[19];
  fk[16] = z[20];
  fk[17] = z[21];
  fk[18] = q2 + par__Rankle_y*z[19] + par__Rhip_x*z[1] + z[5]*z[12] + z[6]*
  z[13] - par__Rhip_y*z[2];
  fk[19] = q3 + par__Rankle_y*z[21] + par__Rhip_x*z[2] + par__Rhip_y*z[1] + 
  z[5]*z[13] + z[7]*z[12];
  fk[20] = z[34];
  fk[21] = z[35];
  fk[22] = z[36];
  fk[23] = z[37];
  fk[24] = q2 + par__Lhip_x*z[1] - par__Lhip_y*z[2];
  fk[25] = q3 + par__Lhip_x*z[2] + par__Lhip_y*z[1];
  fk[26] = z[52];
  fk[27] = z[53];
  fk[28] = z[54];
  fk[29] = z[52];
  fk[30] = q2 + par__Lhip_x*z[1] + z[52]*z[59] + z[53]*z[60] - par__Lhip_y*
  z[2];
  fk[31] = q3 + par__Lhip_x*z[2] + par__Lhip_y*z[1] + z[52]*z[60] + z[54]*
  z[59];
  fk[32] = z[65];
  fk[33] = z[66];
  fk[34] = z[67];
  fk[35] = z[68];
  fk[36] = q2 + par__Lankle_y*z[66] + par__Lhip_x*z[1] + z[52]*z[59] + z[53]*
  z[60] - par__Lhip_y*z[2];
  fk[37] = q3 + par__Lankle_y*z[68] + par__Lhip_x*z[2] + par__Lhip_y*z[1] + 
  z[52]*z[60] + z[54]*z[59];
  fk[38] = z[81];
  fk[39] = z[82];
  fk[40] = z[83];
  fk[41] = z[84];
  fk[42] = q2 + par__back_x*z[1] - par__back_y*z[2];
  fk[43] = q3 + par__back_x*z[2] + par__back_y*z[1];
  fk[44] = z[1];
  fk[45] = -z[2];
  fk[46] = z[2];
  fk[47] = z[1];
  dfk_dq[0][0] = 0;
  dfk_dq[0][1] = 1;
  dfk_dq[0][2] = 0;
  dfk_dq[0][3] = 0;
  dfk_dq[0][4] = 0;
  dfk_dq[0][5] = 0;
  dfk_dq[0][6] = 0;
  dfk_dq[0][7] = 0;
  dfk_dq[0][8] = 0;
  dfk_dq[1][0] = 0;
  dfk_dq[1][1] = 0;
  dfk_dq[1][2] = 1;
  dfk_dq[1][3] = 0;
  dfk_dq[1][4] = 0;
  dfk_dq[1][5] = 0;
  dfk_dq[1][6] = 0;
  dfk_dq[1][7] = 0;
  dfk_dq[1][8] = 0;
  dfk_dq[2][0] = -z[2];
  dfk_dq[2][1] = 0;
  dfk_dq[2][2] = 0;
  dfk_dq[2][3] = 0;
  dfk_dq[2][4] = 0;
  dfk_dq[2][5] = 0;
  dfk_dq[2][6] = 0;
  dfk_dq[2][7] = 0;
  dfk_dq[2][8] = 0;
  dfk_dq[3][0] = -z[1];
  dfk_dq[3][1] = 0;
  dfk_dq[3][2] = 0;
  dfk_dq[3][3] = 0;
  dfk_dq[3][4] = 0;
  dfk_dq[3][5] = 0;
  dfk_dq[3][6] = 0;
  dfk_dq[3][7] = 0;
  dfk_dq[3][8] = 0;
  dfk_dq[4][0] = z[1];
  dfk_dq[4][1] = 0;
  dfk_dq[4][2] = 0;
  dfk_dq[4][3] = 0;
  dfk_dq[4][4] = 0;
  dfk_dq[4][5] = 0;
  dfk_dq[4][6] = 0;
  dfk_dq[4][7] = 0;
  dfk_dq[4][8] = 0;
  dfk_dq[5][0] = -z[2];
  dfk_dq[5][1] = 0;
  dfk_dq[5][2] = 0;
  dfk_dq[5][3] = 0;
  dfk_dq[5][4] = 0;
  dfk_dq[5][5] = 0;
  dfk_dq[5][6] = 0;
  dfk_dq[5][7] = 0;
  dfk_dq[5][8] = 0;
  dfk_dq[6][0] = z[129];
  dfk_dq[6][1] = 1;
  dfk_dq[6][2] = 0;
  dfk_dq[6][3] = 0;
  dfk_dq[6][4] = 0;
  dfk_dq[6][5] = 0;
  dfk_dq[6][6] = 0;
  dfk_dq[6][7] = 0;
  dfk_dq[6][8] = 0;
  dfk_dq[7][0] = z[130];
  dfk_dq[7][1] = 0;
  dfk_dq[7][2] = 1;
  dfk_dq[7][3] = 0;
  dfk_dq[7][4] = 0;
  dfk_dq[7][5] = 0;
  dfk_dq[7][6] = 0;
  dfk_dq[7][7] = 0;
  dfk_dq[7][8] = 0;
  dfk_dq[8][0] = z[6];
  dfk_dq[8][1] = 0;
  dfk_dq[8][2] = 0;
  dfk_dq[8][3] = z[6];
  dfk_dq[8][4] = 0;
  dfk_dq[8][5] = 0;
  dfk_dq[8][6] = 0;
  dfk_dq[8][7] = 0;
  dfk_dq[8][8] = 0;
  dfk_dq[9][0] = z[99];
  dfk_dq[9][1] = 0;
  dfk_dq[9][2] = 0;
  dfk_dq[9][3] = z[99];
  dfk_dq[9][4] = 0;
  dfk_dq[9][5] = 0;
  dfk_dq[9][6] = 0;
  dfk_dq[9][7] = 0;
  dfk_dq[9][8] = 0;
  dfk_dq[10][0] = z[5];
  dfk_dq[10][1] = 0;
  dfk_dq[10][2] = 0;
  dfk_dq[10][3] = z[5];
  dfk_dq[10][4] = 0;
  dfk_dq[10][5] = 0;
  dfk_dq[10][6] = 0;
  dfk_dq[10][7] = 0;
  dfk_dq[10][8] = 0;
  dfk_dq[11][0] = z[6];
  dfk_dq[11][1] = 0;
  dfk_dq[11][2] = 0;
  dfk_dq[11][3] = z[6];
  dfk_dq[11][4] = 0;
  dfk_dq[11][5] = 0;
  dfk_dq[11][6] = 0;
  dfk_dq[11][7] = 0;
  dfk_dq[11][8] = 0;
  dfk_dq[12][0] = z[131];
  dfk_dq[12][1] = 1;
  dfk_dq[12][2] = 0;
  dfk_dq[12][3] = z[132];
  dfk_dq[12][4] = z[133];
  dfk_dq[12][5] = 0;
  dfk_dq[12][6] = 0;
  dfk_dq[12][7] = 0;
  dfk_dq[12][8] = 0;
  dfk_dq[13][0] = z[134];
  dfk_dq[13][1] = 0;
  dfk_dq[13][2] = 1;
  dfk_dq[13][3] = z[135];
  dfk_dq[13][4] = z[136];
  dfk_dq[13][5] = 0;
  dfk_dq[13][6] = 0;
  dfk_dq[13][7] = 0;
  dfk_dq[13][8] = 0;
  dfk_dq[14][0] = z[102];
  dfk_dq[14][1] = 0;
  dfk_dq[14][2] = 0;
  dfk_dq[14][3] = z[102];
  dfk_dq[14][4] = z[19];
  dfk_dq[14][5] = 0;
  dfk_dq[14][6] = 0;
  dfk_dq[14][7] = 0;
  dfk_dq[14][8] = 0;
  dfk_dq[15][0] = z[103];
  dfk_dq[15][1] = 0;
  dfk_dq[15][2] = 0;
  dfk_dq[15][3] = z[103];
  dfk_dq[15][4] = z[104];
  dfk_dq[15][5] = 0;
  dfk_dq[15][6] = 0;
  dfk_dq[15][7] = 0;
  dfk_dq[15][8] = 0;
  dfk_dq[16][0] = z[18];
  dfk_dq[16][1] = 0;
  dfk_dq[16][2] = 0;
  dfk_dq[16][3] = z[18];
  dfk_dq[16][4] = z[21];
  dfk_dq[16][5] = 0;
  dfk_dq[16][6] = 0;
  dfk_dq[16][7] = 0;
  dfk_dq[16][8] = 0;
  dfk_dq[17][0] = z[19];
  dfk_dq[17][1] = 0;
  dfk_dq[17][2] = 0;
  dfk_dq[17][3] = z[19];
  dfk_dq[17][4] = z[105];
  dfk_dq[17][5] = 0;
  dfk_dq[17][6] = 0;
  dfk_dq[17][7] = 0;
  dfk_dq[17][8] = 0;
  dfk_dq[18][0] = z[137];
  dfk_dq[18][1] = 1;
  dfk_dq[18][2] = 0;
  dfk_dq[18][3] = z[138];
  dfk_dq[18][4] = z[139];
  dfk_dq[18][5] = 0;
  dfk_dq[18][6] = 0;
  dfk_dq[18][7] = 0;
  dfk_dq[18][8] = 0;
  dfk_dq[19][0] = z[140];
  dfk_dq[19][1] = 0;
  dfk_dq[19][2] = 1;
  dfk_dq[19][3] = z[141];
  dfk_dq[19][4] = z[142];
  dfk_dq[19][5] = 0;
  dfk_dq[19][6] = 0;
  dfk_dq[19][7] = 0;
  dfk_dq[19][8] = 0;
  dfk_dq[20][0] = z[106];
  dfk_dq[20][1] = 0;
  dfk_dq[20][2] = 0;
  dfk_dq[20][3] = z[106];
  dfk_dq[20][4] = z[107];
  dfk_dq[20][5] = z[35];
  dfk_dq[20][6] = 0;
  dfk_dq[20][7] = 0;
  dfk_dq[20][8] = 0;
  dfk_dq[21][0] = z[108];
  dfk_dq[21][1] = 0;
  dfk_dq[21][2] = 0;
  dfk_dq[21][3] = z[108];
  dfk_dq[21][4] = z[109];
  dfk_dq[21][5] = z[110];
  dfk_dq[21][6] = 0;
  dfk_dq[21][7] = 0;
  dfk_dq[21][8] = 0;
  dfk_dq[22][0] = z[34];
  dfk_dq[22][1] = 0;
  dfk_dq[22][2] = 0;
  dfk_dq[22][3] = z[34];
  dfk_dq[22][4] = z[111];
  dfk_dq[22][5] = z[37];
  dfk_dq[22][6] = 0;
  dfk_dq[22][7] = 0;
  dfk_dq[22][8] = 0;
  dfk_dq[23][0] = z[35];
  dfk_dq[23][1] = 0;
  dfk_dq[23][2] = 0;
  dfk_dq[23][3] = z[35];
  dfk_dq[23][4] = z[112];
  dfk_dq[23][5] = z[113];
  dfk_dq[23][6] = 0;
  dfk_dq[23][7] = 0;
  dfk_dq[23][8] = 0;
  dfk_dq[24][0] = z[143];
  dfk_dq[24][1] = 1;
  dfk_dq[24][2] = 0;
  dfk_dq[24][3] = 0;
  dfk_dq[24][4] = 0;
  dfk_dq[24][5] = 0;
  dfk_dq[24][6] = 0;
  dfk_dq[24][7] = 0;
  dfk_dq[24][8] = 0;
  dfk_dq[25][0] = z[144];
  dfk_dq[25][1] = 0;
  dfk_dq[25][2] = 1;
  dfk_dq[25][3] = 0;
  dfk_dq[25][4] = 0;
  dfk_dq[25][5] = 0;
  dfk_dq[25][6] = 0;
  dfk_dq[25][7] = 0;
  dfk_dq[25][8] = 0;
  dfk_dq[26][0] = z[53];
  dfk_dq[26][1] = 0;
  dfk_dq[26][2] = 0;
  dfk_dq[26][3] = 0;
  dfk_dq[26][4] = 0;
  dfk_dq[26][5] = 0;
  dfk_dq[26][6] = z[53];
  dfk_dq[26][7] = 0;
  dfk_dq[26][8] = 0;
  dfk_dq[27][0] = z[114];
  dfk_dq[27][1] = 0;
  dfk_dq[27][2] = 0;
  dfk_dq[27][3] = 0;
  dfk_dq[27][4] = 0;
  dfk_dq[27][5] = 0;
  dfk_dq[27][6] = z[114];
  dfk_dq[27][7] = 0;
  dfk_dq[27][8] = 0;
  dfk_dq[28][0] = z[52];
  dfk_dq[28][1] = 0;
  dfk_dq[28][2] = 0;
  dfk_dq[28][3] = 0;
  dfk_dq[28][4] = 0;
  dfk_dq[28][5] = 0;
  dfk_dq[28][6] = z[52];
  dfk_dq[28][7] = 0;
  dfk_dq[28][8] = 0;
  dfk_dq[29][0] = z[53];
  dfk_dq[29][1] = 0;
  dfk_dq[29][2] = 0;
  dfk_dq[29][3] = 0;
  dfk_dq[29][4] = 0;
  dfk_dq[29][5] = 0;
  dfk_dq[29][6] = z[53];
  dfk_dq[29][7] = 0;
  dfk_dq[29][8] = 0;
  dfk_dq[30][0] = z[145];
  dfk_dq[30][1] = 1;
  dfk_dq[30][2] = 0;
  dfk_dq[30][3] = 0;
  dfk_dq[30][4] = 0;
  dfk_dq[30][5] = 0;
  dfk_dq[30][6] = z[146];
  dfk_dq[30][7] = z[147];
  dfk_dq[30][8] = 0;
  dfk_dq[31][0] = z[148];
  dfk_dq[31][1] = 0;
  dfk_dq[31][2] = 1;
  dfk_dq[31][3] = 0;
  dfk_dq[31][4] = 0;
  dfk_dq[31][5] = 0;
  dfk_dq[31][6] = z[149];
  dfk_dq[31][7] = z[150];
  dfk_dq[31][8] = 0;
  dfk_dq[32][0] = z[117];
  dfk_dq[32][1] = 0;
  dfk_dq[32][2] = 0;
  dfk_dq[32][3] = 0;
  dfk_dq[32][4] = 0;
  dfk_dq[32][5] = 0;
  dfk_dq[32][6] = z[117];
  dfk_dq[32][7] = z[66];
  dfk_dq[32][8] = 0;
  dfk_dq[33][0] = z[118];
  dfk_dq[33][1] = 0;
  dfk_dq[33][2] = 0;
  dfk_dq[33][3] = 0;
  dfk_dq[33][4] = 0;
  dfk_dq[33][5] = 0;
  dfk_dq[33][6] = z[118];
  dfk_dq[33][7] = z[119];
  dfk_dq[33][8] = 0;
  dfk_dq[34][0] = z[65];
  dfk_dq[34][1] = 0;
  dfk_dq[34][2] = 0;
  dfk_dq[34][3] = 0;
  dfk_dq[34][4] = 0;
  dfk_dq[34][5] = 0;
  dfk_dq[34][6] = z[65];
  dfk_dq[34][7] = z[68];
  dfk_dq[34][8] = 0;
  dfk_dq[35][0] = z[66];
  dfk_dq[35][1] = 0;
  dfk_dq[35][2] = 0;
  dfk_dq[35][3] = 0;
  dfk_dq[35][4] = 0;
  dfk_dq[35][5] = 0;
  dfk_dq[35][6] = z[66];
  dfk_dq[35][7] = z[120];
  dfk_dq[35][8] = 0;
  dfk_dq[36][0] = z[151];
  dfk_dq[36][1] = 1;
  dfk_dq[36][2] = 0;
  dfk_dq[36][3] = 0;
  dfk_dq[36][4] = 0;
  dfk_dq[36][5] = 0;
  dfk_dq[36][6] = z[152];
  dfk_dq[36][7] = z[153];
  dfk_dq[36][8] = 0;
  dfk_dq[37][0] = z[154];
  dfk_dq[37][1] = 0;
  dfk_dq[37][2] = 1;
  dfk_dq[37][3] = 0;
  dfk_dq[37][4] = 0;
  dfk_dq[37][5] = 0;
  dfk_dq[37][6] = z[155];
  dfk_dq[37][7] = z[156];
  dfk_dq[37][8] = 0;
  dfk_dq[38][0] = z[121];
  dfk_dq[38][1] = 0;
  dfk_dq[38][2] = 0;
  dfk_dq[38][3] = 0;
  dfk_dq[38][4] = 0;
  dfk_dq[38][5] = 0;
  dfk_dq[38][6] = z[121];
  dfk_dq[38][7] = z[122];
  dfk_dq[38][8] = z[82];
  dfk_dq[39][0] = z[123];
  dfk_dq[39][1] = 0;
  dfk_dq[39][2] = 0;
  dfk_dq[39][3] = 0;
  dfk_dq[39][4] = 0;
  dfk_dq[39][5] = 0;
  dfk_dq[39][6] = z[123];
  dfk_dq[39][7] = z[124];
  dfk_dq[39][8] = z[125];
  dfk_dq[40][0] = z[81];
  dfk_dq[40][1] = 0;
  dfk_dq[40][2] = 0;
  dfk_dq[40][3] = 0;
  dfk_dq[40][4] = 0;
  dfk_dq[40][5] = 0;
  dfk_dq[40][6] = z[81];
  dfk_dq[40][7] = z[126];
  dfk_dq[40][8] = z[84];
  dfk_dq[41][0] = z[82];
  dfk_dq[41][1] = 0;
  dfk_dq[41][2] = 0;
  dfk_dq[41][3] = 0;
  dfk_dq[41][4] = 0;
  dfk_dq[41][5] = 0;
  dfk_dq[41][6] = z[82];
  dfk_dq[41][7] = z[127];
  dfk_dq[41][8] = z[128];
  dfk_dq[42][0] = z[157];
  dfk_dq[42][1] = 1;
  dfk_dq[42][2] = 0;
  dfk_dq[42][3] = 0;
  dfk_dq[42][4] = 0;
  dfk_dq[42][5] = 0;
  dfk_dq[42][6] = 0;
  dfk_dq[42][7] = 0;
  dfk_dq[42][8] = 0;
  dfk_dq[43][0] = z[158];
  dfk_dq[43][1] = 0;
  dfk_dq[43][2] = 1;
  dfk_dq[43][3] = 0;
  dfk_dq[43][4] = 0;
  dfk_dq[43][5] = 0;
  dfk_dq[43][6] = 0;
  dfk_dq[43][7] = 0;
  dfk_dq[43][8] = 0;
  dfk_dq[44][0] = -z[2];
  dfk_dq[44][1] = 0;
  dfk_dq[44][2] = 0;
  dfk_dq[44][3] = 0;
  dfk_dq[44][4] = 0;
  dfk_dq[44][5] = 0;
  dfk_dq[44][6] = 0;
  dfk_dq[44][7] = 0;
  dfk_dq[44][8] = 0;
  dfk_dq[45][0] = -z[1];
  dfk_dq[45][1] = 0;
  dfk_dq[45][2] = 0;
  dfk_dq[45][3] = 0;
  dfk_dq[45][4] = 0;
  dfk_dq[45][5] = 0;
  dfk_dq[45][6] = 0;
  dfk_dq[45][7] = 0;
  dfk_dq[45][8] = 0;
  dfk_dq[46][0] = z[1];
  dfk_dq[46][1] = 0;
  dfk_dq[46][2] = 0;
  dfk_dq[46][3] = 0;
  dfk_dq[46][4] = 0;
  dfk_dq[46][5] = 0;
  dfk_dq[46][6] = 0;
  dfk_dq[46][7] = 0;
  dfk_dq[46][8] = 0;
  dfk_dq[47][0] = -z[2];
  dfk_dq[47][1] = 0;
  dfk_dq[47][2] = 0;
  dfk_dq[47][3] = 0;
  dfk_dq[47][4] = 0;
  dfk_dq[47][5] = 0;
  dfk_dq[47][6] = 0;
  dfk_dq[47][7] = 0;
  dfk_dq[47][8] = 0;
  fkdot[0] = q2p;
  fkdot[1] = q3p;
  fkdot[2] = -z[159];
  fkdot[3] = -z[160];
  fkdot[4] = z[160];
  fkdot[5] = -z[159];
  fkdot[6] = z[161];
  fkdot[7] = z[162];
  fkdot[8] = z[163];
  fkdot[9] = z[164];
  fkdot[10] = z[165];
  fkdot[11] = z[163];
  fkdot[12] = z[166];
  fkdot[13] = z[167];
  fkdot[14] = z[168];
  fkdot[15] = z[169];
  fkdot[16] = z[170];
  fkdot[17] = z[171];
  fkdot[18] = z[172];
  fkdot[19] = z[173];
  fkdot[20] = z[174];
  fkdot[21] = z[175];
  fkdot[22] = z[176];
  fkdot[23] = z[177];
  fkdot[24] = z[178];
  fkdot[25] = z[179];
  fkdot[26] = z[180];
  fkdot[27] = z[181];
  fkdot[28] = z[182];
  fkdot[29] = z[180];
  fkdot[30] = z[183];
  fkdot[31] = z[184];
  fkdot[32] = z[185];
  fkdot[33] = z[186];
  fkdot[34] = z[187];
  fkdot[35] = z[188];
  fkdot[36] = z[189];
  fkdot[37] = z[190];
  fkdot[38] = z[191];
  fkdot[39] = z[192];
  fkdot[40] = z[193];
  fkdot[41] = z[194];
  fkdot[42] = z[195];
  fkdot[43] = z[196];
  fkdot[44] = -z[159];
  fkdot[45] = -z[160];
  fkdot[46] = z[160];
  fkdot[47] = -z[159];
  dfkdot_dq[0][0] = 0;
  dfkdot_dq[0][1] = 0;
  dfkdot_dq[0][2] = 0;
  dfkdot_dq[0][3] = 0;
  dfkdot_dq[0][4] = 0;
  dfkdot_dq[0][5] = 0;
  dfkdot_dq[0][6] = 0;
  dfkdot_dq[0][7] = 0;
  dfkdot_dq[0][8] = 0;
  dfkdot_dq[1][0] = 0;
  dfkdot_dq[1][1] = 0;
  dfkdot_dq[1][2] = 0;
  dfkdot_dq[1][3] = 0;
  dfkdot_dq[1][4] = 0;
  dfkdot_dq[1][5] = 0;
  dfkdot_dq[1][6] = 0;
  dfkdot_dq[1][7] = 0;
  dfkdot_dq[1][8] = 0;
  dfkdot_dq[2][0] = -z[160];
  dfkdot_dq[2][1] = 0;
  dfkdot_dq[2][2] = 0;
  dfkdot_dq[2][3] = 0;
  dfkdot_dq[2][4] = 0;
  dfkdot_dq[2][5] = 0;
  dfkdot_dq[2][6] = 0;
  dfkdot_dq[2][7] = 0;
  dfkdot_dq[2][8] = 0;
  dfkdot_dq[3][0] = z[159];
  dfkdot_dq[3][1] = 0;
  dfkdot_dq[3][2] = 0;
  dfkdot_dq[3][3] = 0;
  dfkdot_dq[3][4] = 0;
  dfkdot_dq[3][5] = 0;
  dfkdot_dq[3][6] = 0;
  dfkdot_dq[3][7] = 0;
  dfkdot_dq[3][8] = 0;
  dfkdot_dq[4][0] = -z[159];
  dfkdot_dq[4][1] = 0;
  dfkdot_dq[4][2] = 0;
  dfkdot_dq[4][3] = 0;
  dfkdot_dq[4][4] = 0;
  dfkdot_dq[4][5] = 0;
  dfkdot_dq[4][6] = 0;
  dfkdot_dq[4][7] = 0;
  dfkdot_dq[4][8] = 0;
  dfkdot_dq[5][0] = -z[160];
  dfkdot_dq[5][1] = 0;
  dfkdot_dq[5][2] = 0;
  dfkdot_dq[5][3] = 0;
  dfkdot_dq[5][4] = 0;
  dfkdot_dq[5][5] = 0;
  dfkdot_dq[5][6] = 0;
  dfkdot_dq[5][7] = 0;
  dfkdot_dq[5][8] = 0;
  dfkdot_dq[6][0] = -z[197];
  dfkdot_dq[6][1] = 0;
  dfkdot_dq[6][2] = 0;
  dfkdot_dq[6][3] = 0;
  dfkdot_dq[6][4] = 0;
  dfkdot_dq[6][5] = 0;
  dfkdot_dq[6][6] = 0;
  dfkdot_dq[6][7] = 0;
  dfkdot_dq[6][8] = 0;
  dfkdot_dq[7][0] = -z[198];
  dfkdot_dq[7][1] = 0;
  dfkdot_dq[7][2] = 0;
  dfkdot_dq[7][3] = 0;
  dfkdot_dq[7][4] = 0;
  dfkdot_dq[7][5] = 0;
  dfkdot_dq[7][6] = 0;
  dfkdot_dq[7][7] = 0;
  dfkdot_dq[7][8] = 0;
  dfkdot_dq[8][0] = z[164];
  dfkdot_dq[8][1] = 0;
  dfkdot_dq[8][2] = 0;
  dfkdot_dq[8][3] = z[164];
  dfkdot_dq[8][4] = 0;
  dfkdot_dq[8][5] = 0;
  dfkdot_dq[8][6] = 0;
  dfkdot_dq[8][7] = 0;
  dfkdot_dq[8][8] = 0;
  dfkdot_dq[9][0] = z[291];
  dfkdot_dq[9][1] = 0;
  dfkdot_dq[9][2] = 0;
  dfkdot_dq[9][3] = z[291];
  dfkdot_dq[9][4] = 0;
  dfkdot_dq[9][5] = 0;
  dfkdot_dq[9][6] = 0;
  dfkdot_dq[9][7] = 0;
  dfkdot_dq[9][8] = 0;
  dfkdot_dq[10][0] = z[163];
  dfkdot_dq[10][1] = 0;
  dfkdot_dq[10][2] = 0;
  dfkdot_dq[10][3] = z[163];
  dfkdot_dq[10][4] = 0;
  dfkdot_dq[10][5] = 0;
  dfkdot_dq[10][6] = 0;
  dfkdot_dq[10][7] = 0;
  dfkdot_dq[10][8] = 0;
  dfkdot_dq[11][0] = z[164];
  dfkdot_dq[11][1] = 0;
  dfkdot_dq[11][2] = 0;
  dfkdot_dq[11][3] = z[164];
  dfkdot_dq[11][4] = 0;
  dfkdot_dq[11][5] = 0;
  dfkdot_dq[11][6] = 0;
  dfkdot_dq[11][7] = 0;
  dfkdot_dq[11][8] = 0;
  dfkdot_dq[12][0] = z[199];
  dfkdot_dq[12][1] = 0;
  dfkdot_dq[12][2] = 0;
  dfkdot_dq[12][3] = z[200];
  dfkdot_dq[12][4] = z[203];
  dfkdot_dq[12][5] = 0;
  dfkdot_dq[12][6] = 0;
  dfkdot_dq[12][7] = 0;
  dfkdot_dq[12][8] = 0;
  dfkdot_dq[13][0] = z[204];
  dfkdot_dq[13][1] = 0;
  dfkdot_dq[13][2] = 0;
  dfkdot_dq[13][3] = z[205];
  dfkdot_dq[13][4] = z[206];
  dfkdot_dq[13][5] = 0;
  dfkdot_dq[13][6] = 0;
  dfkdot_dq[13][7] = 0;
  dfkdot_dq[13][8] = 0;
  dfkdot_dq[14][0] = z[208];
  dfkdot_dq[14][1] = 0;
  dfkdot_dq[14][2] = 0;
  dfkdot_dq[14][3] = z[208];
  dfkdot_dq[14][4] = z[169];
  dfkdot_dq[14][5] = 0;
  dfkdot_dq[14][6] = 0;
  dfkdot_dq[14][7] = 0;
  dfkdot_dq[14][8] = 0;
  dfkdot_dq[15][0] = z[211];
  dfkdot_dq[15][1] = 0;
  dfkdot_dq[15][2] = 0;
  dfkdot_dq[15][3] = z[211];
  dfkdot_dq[15][4] = z[213];
  dfkdot_dq[15][5] = 0;
  dfkdot_dq[15][6] = 0;
  dfkdot_dq[15][7] = 0;
  dfkdot_dq[15][8] = 0;
  dfkdot_dq[16][0] = z[168];
  dfkdot_dq[16][1] = 0;
  dfkdot_dq[16][2] = 0;
  dfkdot_dq[16][3] = z[168];
  dfkdot_dq[16][4] = z[171];
  dfkdot_dq[16][5] = 0;
  dfkdot_dq[16][6] = 0;
  dfkdot_dq[16][7] = 0;
  dfkdot_dq[16][8] = 0;
  dfkdot_dq[17][0] = z[169];
  dfkdot_dq[17][1] = 0;
  dfkdot_dq[17][2] = 0;
  dfkdot_dq[17][3] = z[169];
  dfkdot_dq[17][4] = z[215];
  dfkdot_dq[17][5] = 0;
  dfkdot_dq[17][6] = 0;
  dfkdot_dq[17][7] = 0;
  dfkdot_dq[17][8] = 0;
  dfkdot_dq[18][0] = z[216];
  dfkdot_dq[18][1] = 0;
  dfkdot_dq[18][2] = 0;
  dfkdot_dq[18][3] = z[217];
  dfkdot_dq[18][4] = z[218];
  dfkdot_dq[18][5] = 0;
  dfkdot_dq[18][6] = 0;
  dfkdot_dq[18][7] = 0;
  dfkdot_dq[18][8] = 0;
  dfkdot_dq[19][0] = z[219];
  dfkdot_dq[19][1] = 0;
  dfkdot_dq[19][2] = 0;
  dfkdot_dq[19][3] = z[220];
  dfkdot_dq[19][4] = z[221];
  dfkdot_dq[19][5] = 0;
  dfkdot_dq[19][6] = 0;
  dfkdot_dq[19][7] = 0;
  dfkdot_dq[19][8] = 0;
  dfkdot_dq[20][0] = z[224];
  dfkdot_dq[20][1] = 0;
  dfkdot_dq[20][2] = 0;
  dfkdot_dq[20][3] = z[224];
  dfkdot_dq[20][4] = z[226];
  dfkdot_dq[20][5] = z[175];
  dfkdot_dq[20][6] = 0;
  dfkdot_dq[20][7] = 0;
  dfkdot_dq[20][8] = 0;
  dfkdot_dq[21][0] = z[230];
  dfkdot_dq[21][1] = 0;
  dfkdot_dq[21][2] = 0;
  dfkdot_dq[21][3] = z[230];
  dfkdot_dq[21][4] = z[233];
  dfkdot_dq[21][5] = z[235];
  dfkdot_dq[21][6] = 0;
  dfkdot_dq[21][7] = 0;
  dfkdot_dq[21][8] = 0;
  dfkdot_dq[22][0] = z[174];
  dfkdot_dq[22][1] = 0;
  dfkdot_dq[22][2] = 0;
  dfkdot_dq[22][3] = z[174];
  dfkdot_dq[22][4] = z[237];
  dfkdot_dq[22][5] = z[177];
  dfkdot_dq[22][6] = 0;
  dfkdot_dq[22][7] = 0;
  dfkdot_dq[22][8] = 0;
  dfkdot_dq[23][0] = z[175];
  dfkdot_dq[23][1] = 0;
  dfkdot_dq[23][2] = 0;
  dfkdot_dq[23][3] = z[175];
  dfkdot_dq[23][4] = z[240];
  dfkdot_dq[23][5] = z[242];
  dfkdot_dq[23][6] = 0;
  dfkdot_dq[23][7] = 0;
  dfkdot_dq[23][8] = 0;
  dfkdot_dq[24][0] = -z[243];
  dfkdot_dq[24][1] = 0;
  dfkdot_dq[24][2] = 0;
  dfkdot_dq[24][3] = 0;
  dfkdot_dq[24][4] = 0;
  dfkdot_dq[24][5] = 0;
  dfkdot_dq[24][6] = 0;
  dfkdot_dq[24][7] = 0;
  dfkdot_dq[24][8] = 0;
  dfkdot_dq[25][0] = -z[244];
  dfkdot_dq[25][1] = 0;
  dfkdot_dq[25][2] = 0;
  dfkdot_dq[25][3] = 0;
  dfkdot_dq[25][4] = 0;
  dfkdot_dq[25][5] = 0;
  dfkdot_dq[25][6] = 0;
  dfkdot_dq[25][7] = 0;
  dfkdot_dq[25][8] = 0;
  dfkdot_dq[26][0] = z[181];
  dfkdot_dq[26][1] = 0;
  dfkdot_dq[26][2] = 0;
  dfkdot_dq[26][3] = 0;
  dfkdot_dq[26][4] = 0;
  dfkdot_dq[26][5] = 0;
  dfkdot_dq[26][6] = z[181];
  dfkdot_dq[26][7] = 0;
  dfkdot_dq[26][8] = 0;
  dfkdot_dq[27][0] = z[292];
  dfkdot_dq[27][1] = 0;
  dfkdot_dq[27][2] = 0;
  dfkdot_dq[27][3] = 0;
  dfkdot_dq[27][4] = 0;
  dfkdot_dq[27][5] = 0;
  dfkdot_dq[27][6] = z[292];
  dfkdot_dq[27][7] = 0;
  dfkdot_dq[27][8] = 0;
  dfkdot_dq[28][0] = z[180];
  dfkdot_dq[28][1] = 0;
  dfkdot_dq[28][2] = 0;
  dfkdot_dq[28][3] = 0;
  dfkdot_dq[28][4] = 0;
  dfkdot_dq[28][5] = 0;
  dfkdot_dq[28][6] = z[180];
  dfkdot_dq[28][7] = 0;
  dfkdot_dq[28][8] = 0;
  dfkdot_dq[29][0] = z[181];
  dfkdot_dq[29][1] = 0;
  dfkdot_dq[29][2] = 0;
  dfkdot_dq[29][3] = 0;
  dfkdot_dq[29][4] = 0;
  dfkdot_dq[29][5] = 0;
  dfkdot_dq[29][6] = z[181];
  dfkdot_dq[29][7] = 0;
  dfkdot_dq[29][8] = 0;
  dfkdot_dq[30][0] = z[245];
  dfkdot_dq[30][1] = 0;
  dfkdot_dq[30][2] = 0;
  dfkdot_dq[30][3] = 0;
  dfkdot_dq[30][4] = 0;
  dfkdot_dq[30][5] = 0;
  dfkdot_dq[30][6] = z[246];
  dfkdot_dq[30][7] = z[249];
  dfkdot_dq[30][8] = 0;
  dfkdot_dq[31][0] = z[250];
  dfkdot_dq[31][1] = 0;
  dfkdot_dq[31][2] = 0;
  dfkdot_dq[31][3] = 0;
  dfkdot_dq[31][4] = 0;
  dfkdot_dq[31][5] = 0;
  dfkdot_dq[31][6] = z[251];
  dfkdot_dq[31][7] = z[252];
  dfkdot_dq[31][8] = 0;
  dfkdot_dq[32][0] = z[254];
  dfkdot_dq[32][1] = 0;
  dfkdot_dq[32][2] = 0;
  dfkdot_dq[32][3] = 0;
  dfkdot_dq[32][4] = 0;
  dfkdot_dq[32][5] = 0;
  dfkdot_dq[32][6] = z[254];
  dfkdot_dq[32][7] = z[186];
  dfkdot_dq[32][8] = 0;
  dfkdot_dq[33][0] = z[257];
  dfkdot_dq[33][1] = 0;
  dfkdot_dq[33][2] = 0;
  dfkdot_dq[33][3] = 0;
  dfkdot_dq[33][4] = 0;
  dfkdot_dq[33][5] = 0;
  dfkdot_dq[33][6] = z[257];
  dfkdot_dq[33][7] = z[259];
  dfkdot_dq[33][8] = 0;
  dfkdot_dq[34][0] = z[185];
  dfkdot_dq[34][1] = 0;
  dfkdot_dq[34][2] = 0;
  dfkdot_dq[34][3] = 0;
  dfkdot_dq[34][4] = 0;
  dfkdot_dq[34][5] = 0;
  dfkdot_dq[34][6] = z[185];
  dfkdot_dq[34][7] = z[188];
  dfkdot_dq[34][8] = 0;
  dfkdot_dq[35][0] = z[186];
  dfkdot_dq[35][1] = 0;
  dfkdot_dq[35][2] = 0;
  dfkdot_dq[35][3] = 0;
  dfkdot_dq[35][4] = 0;
  dfkdot_dq[35][5] = 0;
  dfkdot_dq[35][6] = z[186];
  dfkdot_dq[35][7] = z[261];
  dfkdot_dq[35][8] = 0;
  dfkdot_dq[36][0] = z[262];
  dfkdot_dq[36][1] = 0;
  dfkdot_dq[36][2] = 0;
  dfkdot_dq[36][3] = 0;
  dfkdot_dq[36][4] = 0;
  dfkdot_dq[36][5] = 0;
  dfkdot_dq[36][6] = z[263];
  dfkdot_dq[36][7] = z[264];
  dfkdot_dq[36][8] = 0;
  dfkdot_dq[37][0] = z[265];
  dfkdot_dq[37][1] = 0;
  dfkdot_dq[37][2] = 0;
  dfkdot_dq[37][3] = 0;
  dfkdot_dq[37][4] = 0;
  dfkdot_dq[37][5] = 0;
  dfkdot_dq[37][6] = z[266];
  dfkdot_dq[37][7] = z[267];
  dfkdot_dq[37][8] = 0;
  dfkdot_dq[38][0] = z[270];
  dfkdot_dq[38][1] = 0;
  dfkdot_dq[38][2] = 0;
  dfkdot_dq[38][3] = 0;
  dfkdot_dq[38][4] = 0;
  dfkdot_dq[38][5] = 0;
  dfkdot_dq[38][6] = z[270];
  dfkdot_dq[38][7] = z[272];
  dfkdot_dq[38][8] = z[192];
  dfkdot_dq[39][0] = z[276];
  dfkdot_dq[39][1] = 0;
  dfkdot_dq[39][2] = 0;
  dfkdot_dq[39][3] = 0;
  dfkdot_dq[39][4] = 0;
  dfkdot_dq[39][5] = 0;
  dfkdot_dq[39][6] = z[276];
  dfkdot_dq[39][7] = z[279];
  dfkdot_dq[39][8] = z[281];
  dfkdot_dq[40][0] = z[191];
  dfkdot_dq[40][1] = 0;
  dfkdot_dq[40][2] = 0;
  dfkdot_dq[40][3] = 0;
  dfkdot_dq[40][4] = 0;
  dfkdot_dq[40][5] = 0;
  dfkdot_dq[40][6] = z[191];
  dfkdot_dq[40][7] = z[283];
  dfkdot_dq[40][8] = z[194];
  dfkdot_dq[41][0] = z[192];
  dfkdot_dq[41][1] = 0;
  dfkdot_dq[41][2] = 0;
  dfkdot_dq[41][3] = 0;
  dfkdot_dq[41][4] = 0;
  dfkdot_dq[41][5] = 0;
  dfkdot_dq[41][6] = z[192];
  dfkdot_dq[41][7] = z[286];
  dfkdot_dq[41][8] = z[288];
  dfkdot_dq[42][0] = -z[289];
  dfkdot_dq[42][1] = 0;
  dfkdot_dq[42][2] = 0;
  dfkdot_dq[42][3] = 0;
  dfkdot_dq[42][4] = 0;
  dfkdot_dq[42][5] = 0;
  dfkdot_dq[42][6] = 0;
  dfkdot_dq[42][7] = 0;
  dfkdot_dq[42][8] = 0;
  dfkdot_dq[43][0] = -z[290];
  dfkdot_dq[43][1] = 0;
  dfkdot_dq[43][2] = 0;
  dfkdot_dq[43][3] = 0;
  dfkdot_dq[43][4] = 0;
  dfkdot_dq[43][5] = 0;
  dfkdot_dq[43][6] = 0;
  dfkdot_dq[43][7] = 0;
  dfkdot_dq[43][8] = 0;
  dfkdot_dq[44][0] = -z[160];
  dfkdot_dq[44][1] = 0;
  dfkdot_dq[44][2] = 0;
  dfkdot_dq[44][3] = 0;
  dfkdot_dq[44][4] = 0;
  dfkdot_dq[44][5] = 0;
  dfkdot_dq[44][6] = 0;
  dfkdot_dq[44][7] = 0;
  dfkdot_dq[44][8] = 0;
  dfkdot_dq[45][0] = z[159];
  dfkdot_dq[45][1] = 0;
  dfkdot_dq[45][2] = 0;
  dfkdot_dq[45][3] = 0;
  dfkdot_dq[45][4] = 0;
  dfkdot_dq[45][5] = 0;
  dfkdot_dq[45][6] = 0;
  dfkdot_dq[45][7] = 0;
  dfkdot_dq[45][8] = 0;
  dfkdot_dq[46][0] = -z[159];
  dfkdot_dq[46][1] = 0;
  dfkdot_dq[46][2] = 0;
  dfkdot_dq[46][3] = 0;
  dfkdot_dq[46][4] = 0;
  dfkdot_dq[46][5] = 0;
  dfkdot_dq[46][6] = 0;
  dfkdot_dq[46][7] = 0;
  dfkdot_dq[46][8] = 0;
  dfkdot_dq[47][0] = -z[160];
  dfkdot_dq[47][1] = 0;
  dfkdot_dq[47][2] = 0;
  dfkdot_dq[47][3] = 0;
  dfkdot_dq[47][4] = 0;
  dfkdot_dq[47][5] = 0;
  dfkdot_dq[47][6] = 0;
  dfkdot_dq[47][7] = 0;
  dfkdot_dq[47][8] = 0;
  f[0] = (z[356]-z[505]-z[492]*q1pp-z[494]*q7pp-z[496]*q8pp-z[497]*q9pp-
  z[499]*q4pp-z[501]*q5pp-z[502]*q6pp-z[503]*q3pp-z[504]*q2pp)/par__bodyweight;
  f[1] = (z[357]-z[514]-z[504]*q1pp-z[506]*q7pp-z[507]*q3pp-z[508]*q2pp-
  z[509]*q9pp-z[510]*q8pp-z[511]*q4pp-z[512]*q6pp-z[513]*q5pp)/par__bodyweight;
  f[2] = (z[358]-z[522]-z[503]*q1pp-z[507]*q2pp-z[515]*q7pp-z[516]*q3pp-
  z[517]*q9pp-z[518]*q8pp-z[519]*q4pp-z[520]*q6pp-z[521]*q5pp)/par__bodyweight;
  f[3] = (z[363]-z[527]-z[499]*q1pp-z[511]*q2pp-z[519]*q3pp-z[524]*q4pp-
  z[525]*q5pp-z[526]*q6pp)/par__bodyweight;
  f[4] = (z[368]-z[531]-z[501]*q1pp-z[513]*q2pp-z[521]*q3pp-z[525]*q4pp-
  z[529]*q5pp-z[530]*q6pp)/par__bodyweight;
  f[5] = (z[373]-z[533]-z[532]*q6pp-z[502]*q1pp-z[512]*q2pp-z[520]*q3pp-
  z[526]*q4pp-z[530]*q5pp)/par__bodyweight;
  f[6] = (z[378]-z[538]-z[494]*q1pp-z[506]*q2pp-z[515]*q3pp-z[535]*q7pp-
  z[536]*q8pp-z[537]*q9pp)/par__bodyweight;
  f[7] = (z[383]-z[542]-z[496]*q1pp-z[510]*q2pp-z[518]*q3pp-z[536]*q7pp-
  z[540]*q8pp-z[541]*q9pp)/par__bodyweight;
  f[8] = (z[388]-z[544]-z[543]*q9pp-z[497]*q1pp-z[509]*q2pp-z[517]*q3pp-
  z[537]*q7pp-z[541]*q8pp)/par__bodyweight;

  Encode[0] = 0.0;
  Encode[1] = 0.0;
  Encode[2] = 0.0;
  Encode[3] = 0.0;

  Encode[4] = 0.0;

}


/*................................... READF ................................*/
void  readf( FILE *Fp, double *next, ... )
{
va_list args;                                    /* Variable argument list  */
for( va_start(args,next);  next;  next=va_arg(args,double *) )
  pgets(Fp,next);
va_end(args);                          /* Help function make normal return  */
pgets(Fp,NULL);                        /* Always get a newline at the end   */
}


/*................................... PGETS ................................*/
void  pgets( FILE *Fp, double *x )
{
static long  lineNumber = 0;
char         line[256];

lineNumber++;
if( !fgets(line,256,Fp) )
{
   printf("\n Unable to read line %ld of input file."
          "\n Premature end of file found while reading input file\n", lineNumber);
   exit(0);
}
if( !x ) return;
if( strlen(line) >= 60 )
{
   char *endOfNumber;
   *x = strtod(line+59,&endOfNumber);
   while( isspace(*endOfNumber) )  endOfNumber++;
   if( !*endOfNumber ) return;
}
printf("\n An error occured while reading line %ld of the input file."
       "\n The program was expecting to find one double precision number"
       "\n beginning with the 60th character in the following line:\n\n%s\n",
           lineNumber, line );
exit(0);
}


/* .................................. WRITEF .............................. */
void  writef( FILE *Fp, char format[], ... )
{
va_list  args;                                   /* Variable argument list  */
double   next;                                   /* Current place in list   */

va_start(args,format);                           /* args start after format */
while( (next=va_arg(args,double)) != _NAN )  fprintf(Fp, format, next);
va_end(args);                                    /* End of variable list    */
fprintf(Fp, "\n");                               /* End with newline        */
}


