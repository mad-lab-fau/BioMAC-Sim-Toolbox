/* The name of this program is gait2d_osim/gait2d_osim_raw.c */
/* Created by Autolev 4.1 on Thu Nov 28 15:08:30 2024 */

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
double   Pi,DEGtoRAD,RADtoDEG,z[1268],fk[48],dfk_dq[48][9],fkdot[48],dfkdot_dq[
  48][9],f[9],df_dq[9][9],df_dqd[9][9],df_dqdd[9][9],df_dG[9][6],Encode[9];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait2d_osim/gait2d_osim_raw.in");
  else sprintf(fileName, "gait2d_osim/gait2d_osim_raw.%d", iloop);
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
puts( "\n Input is in the file gait2d_osim/gait2d_osim_raw.in" );
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
  z[545] = (par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2])*q1p;
  z[546] = (par__pelvis_CMx*z[2]+par__pelvis_CMy*z[1])*q1p;
  z[547] = (sx*z[545]+sy*z[546])/pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[548] = par__airdrag*(s*(z[1]*q2p+z[2]*q3p)+(par__pelvis_CMx*q1p+z[1]*q3p-
  z[2]*q2p)*z[547]);
  z[549] = par__airdrag*(s*(z[1]*q3p-z[2]*q2p)+(par__pelvis_CMy*q1p-z[1]*q2p-
  z[2]*q3p)*z[547]);
  z[550] = par__pelvis_CMx*z[309];
  z[551] = par__pelvis_CMy*z[309];
  z[552] = z[352]*z[1] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*z[117]*
  z[344] + par__Lhip_x*z[1]*z[344] + par__Lhip_x*z[2]*z[345] + par__Lhip_y*
  z[1]*z[345] + par__pelvis_CMy*z[2]*z[313] + par__Rankle_y*z[18]*z[328] + 
  par__Rankle_y*z[102]*z[327] + par__Rhip_x*z[1]*z[327] + par__Rhip_x*z[2]*
  z[328] + par__Rhip_y*z[1]*z[328] + z[293]*z[53]*z[57] + z[293]*z[58]*z[114] + 
  z[294]*z[52]*z[57] + z[294]*z[53]*z[58] + z[295]*z[92]*z[121] + z[295]*
  z[95]*z[123] + z[296]*z[81]*z[92] + z[296]*z[82]*z[95] + z[297]*z[74]*
  z[117] + z[297]*z[76]*z[118] + z[298]*z[65]*z[74] + z[298]*z[66]*z[76] + 
  z[301]*z[6]*z[10] + z[301]*z[11]*z[99] + z[302]*z[5]*z[10] + z[302]*z[6]*
  z[11] + z[303]*z[45]*z[106] + z[303]*z[48]*z[108] + z[304]*z[34]*z[45] + 
  z[304]*z[35]*z[48] + z[305]*z[27]*z[102] + z[305]*z[29]*z[103] + z[306]*
  z[18]*z[27] + z[306]*z[19]*z[29] + z[5]*z[13]*z[328] + z[6]*z[13]*z[327] + 
  z[34]*z[322]*z[328] + z[35]*z[325]*z[328] + z[52]*z[60]*z[345] + z[53]*
  z[60]*z[344] + z[81]*z[339]*z[345] + z[82]*z[342]*z[345] + z[106]*z[322]*
  z[327] + z[108]*z[325]*z[327] + z[121]*z[339]*z[344] + z[123]*z[342]*z[344] + 
  par__pelvis_CMx*z[548] + par__pelvis_CMy*z[549] + z[550]*z[2]*z[547] + 
  z[551]*z[1]*z[547] - z[350]*z[2] - z[351]*z[2] - z[353]*z[2] - z[354]*z[1] - 
  z[355]*z[1] - par__Lhip_y*z[2]*z[344] - par__pelvis_CMx*z[1]*z[313] - 
  par__Rhip_y*z[2]*z[327] - z[6]*z[12]*z[328] - z[12]*z[99]*z[327] - z[53]*
  z[59]*z[345] - z[59]*z[114]*z[344];
  z[553] = q2p*z[401] + q3p*z[397];
  z[554] = z[53]*q3p + z[114]*q2p;
  z[555] = z[1]*z[51]*q1p + z[1]*z[51]*q7p + z[2]*z[50]*q1p + z[2]*z[50]*q7p;
  z[556] = q2p*z[555] + q3p*z[401];
  z[557] = z[52]*q3p + z[53]*q2p;
  z[558] = z[63]*z[114]*q8p + z[63]*z[401] + z[64]*z[555] - z[53]*z[64]*q8p;
  z[559] = z[63]*z[555] - z[53]*z[63]*q8p - z[64]*z[114]*q8p - z[64]*z[401];
  z[560] = z[79]*z[118]*q9p + z[79]*z[558] + z[80]*z[559] - z[80]*z[117]*q9p;
  z[561] = q2p*z[560] + q3p*z[410];
  z[562] = z[82]*q3p + z[123]*q2p;
  z[563] = z[79]*z[559] - z[79]*z[117]*q9p - z[80]*z[118]*q9p - z[80]*z[558];
  z[564] = q2p*z[563] + q3p*z[424];
  z[565] = z[81]*q3p + z[121]*q2p;
  z[566] = q2p*z[558] + q3p*z[408];
  z[567] = z[66]*q3p + z[118]*q2p;
  z[568] = q2p*z[559] + q3p*z[409];
  z[569] = z[65]*q3p + z[117]*q2p;
  z[570] = q2p*z[450] + q3p*z[446];
  z[571] = z[6]*q3p + z[99]*q2p;
  z[572] = z[1]*z[4]*q1p + z[1]*z[4]*q4p + z[2]*z[3]*q1p + z[2]*z[3]*q4p;
  z[573] = q2p*z[572] + q3p*z[450];
  z[574] = z[5]*q3p + z[6]*q2p;
  z[575] = z[16]*z[99]*q5p + z[16]*z[450] + z[17]*z[572] - z[6]*z[17]*q5p;
  z[576] = z[16]*z[572] - z[6]*z[16]*q5p - z[17]*z[99]*q5p - z[17]*z[450];
  z[577] = z[32]*z[103]*q6p + z[32]*z[575] + z[33]*z[576] - z[33]*z[102]*q6p;
  z[578] = q2p*z[577] + q3p*z[459];
  z[579] = z[35]*q3p + z[108]*q2p;
  z[580] = z[32]*z[576] - z[32]*z[102]*q6p - z[33]*z[103]*q6p - z[33]*z[575];
  z[581] = q2p*z[580] + q3p*z[473];
  z[582] = z[34]*q3p + z[106]*q2p;
  z[583] = q2p*z[575] + q3p*z[457];
  z[584] = z[19]*q3p + z[103]*q2p;
  z[585] = q2p*z[576] + q3p*z[458];
  z[586] = z[18]*q3p + z[102]*q2p;
  z[587] = -z[1]*q2p - z[2]*q3p;
  z[588] = z[1]*q3p - z[2]*q2p;
  z[589] = par__Lfemur_M*(z[58]*(z[556]+z[389]*z[557])+z[57]*(z[553]-z[389]*
  z[554])) + par__Lfoot_M*(z[95]*(z[564]+z[390]*z[565])+z[92]*(z[561]-z[390]*
  z[562])) + par__Ltibia_M*(z[76]*(z[568]+z[391]*z[569])+z[74]*(z[566]-z[391]*
  z[567])) + par__Rfemur_M*(z[11]*(z[573]+z[392]*z[574])+z[10]*(z[570]-z[392]*
  z[571])) + par__Rfoot_M*(z[48]*(z[581]+z[393]*z[582])+z[45]*(z[578]-z[393]*
  z[579])) + par__Rtibia_M*(z[29]*(z[585]+z[394]*z[586])+z[27]*(z[583]-z[394]*
  z[584])) + par__pelvis_M*(par__pelvis_CMy*(z[441]+q1p*z[587])-par__pelvis_CMx*(
  z[440]-q1p*z[588])) - par__torso_M*(z[97]*(z[441]+q1p*z[587])+z[98]*(z[440]-
  q1p*z[588]));
  z[590] = par__Lfemur_M*(z[53]*z[57]+z[58]*z[114]) + par__Lfoot_M*(z[92]*
  z[121]+z[95]*z[123]) + par__Ltibia_M*(z[74]*z[117]+z[76]*z[118]) + par__Rfemur_M*(
  z[6]*z[10]+z[11]*z[99]) + par__Rfoot_M*(z[45]*z[106]+z[48]*z[108]) + par__Rtibia_M*(
  z[27]*z[102]+z[29]*z[103]) - par__torso_M*(z[97]*z[2]+z[98]*z[1]) - par__pelvis_M*(
  par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2]);
  z[591] = par__Rhip_y*z[3] - par__Rhip_x*z[4];
  z[592] = z[16]*z[591] - z[9]*z[17];
  z[593] = z[9]*z[16] + z[17]*z[591];
  z[594] = z[32]*z[593] + z[33]*z[592];
  z[595] = z[32]*z[592] - z[33]*z[593];
  z[596] = par__Rankle_y*z[18]*z[328] + par__Rankle_y*z[102]*z[327] + z[301]*
  z[5]*z[9] + z[301]*z[6]*z[10] + z[301]*z[6]*z[591] + z[301]*z[11]*z[99] + 
  z[302]*z[5]*z[10] + z[302]*z[5]*z[591] + z[302]*z[6]*z[11] + z[302]*z[7]*
  z[9] + z[303]*z[34]*z[594] + z[303]*z[35]*z[595] + z[303]*z[45]*z[106] + 
  z[303]*z[48]*z[108] + z[304]*z[34]*z[45] + z[304]*z[35]*z[48] + z[304]*
  z[36]*z[594] + z[304]*z[37]*z[595] + z[305]*z[18]*z[593] + z[305]*z[19]*
  z[592] + z[305]*z[27]*z[102] + z[305]*z[29]*z[103] + z[306]*z[18]*z[27] + 
  z[306]*z[19]*z[29] + z[306]*z[20]*z[593] + z[306]*z[21]*z[592] + z[5]*z[13]*
  z[328] + z[6]*z[13]*z[327] + z[34]*z[322]*z[328] + z[34]*z[327]*z[594] + 
  z[35]*z[325]*z[328] + z[35]*z[327]*z[595] + z[36]*z[328]*z[594] + z[37]*
  z[328]*z[595] + z[106]*z[322]*z[327] + z[108]*z[325]*z[327] - z[6]*z[12]*
  z[328] - z[12]*z[99]*z[327];
  z[597] = z[6]*q3p + z[99]*q2p + z[591]*q1p;
  z[598] = q2p*z[572] + q3p*z[450] - q1p*z[448];
  z[599] = z[5]*q3p + z[6]*q2p + z[9]*q1p;
  z[600] = z[16]*z[591]*q5p - z[9]*z[17]*q5p - z[16]*z[451] - z[17]*z[448];
  z[601] = z[17]*z[451] - z[9]*z[16]*q5p - z[17]*z[591]*q5p - z[16]*z[448];
  z[602] = z[32]*z[592]*q6p + z[32]*z[600] + z[33]*z[601] - z[33]*z[593]*q6p;
  z[603] = q1p*z[602] + q2p*z[577] + q3p*z[459];
  z[604] = z[35]*q3p + z[108]*q2p + z[595]*q1p;
  z[605] = z[32]*z[601] - z[32]*z[593]*q6p - z[33]*z[592]*q6p - z[33]*z[600];
  z[606] = q1p*z[605] + q2p*z[580] + q3p*z[473];
  z[607] = z[34]*q3p + z[106]*q2p + z[594]*q1p;
  z[608] = q1p*z[600] + q2p*z[575] + q3p*z[457];
  z[609] = z[19]*q3p + z[103]*q2p + z[592]*q1p;
  z[610] = q1p*z[601] + q2p*z[576] + q3p*z[458];
  z[611] = z[18]*q3p + z[102]*q2p + z[593]*q1p;
  z[612] = par__Rfemur_M*(z[9]*z[453]+z[591]*z[454]+z[11]*(z[598]+z[392]*
  z[599])+z[10]*(z[452]-z[392]*z[597])) + par__Rfoot_M*(z[594]*z[479]+z[595]*
  z[480]+z[48]*(z[606]+z[393]*z[607])+z[45]*(z[603]-z[393]*z[604])) + par__Rtibia_M*(
  z[592]*z[486]+z[593]*z[485]+z[29]*(z[610]+z[394]*z[611])+z[27]*(z[608]-
  z[394]*z[609]));
  z[613] = 2*par__Rfemur_M*(z[9]*z[10]+z[11]*z[591]) + 2*par__Rfoot_M*(z[45]*
  z[594]+z[48]*z[595]) + 2*par__Rtibia_M*(z[27]*z[593]+z[29]*z[592]);
  z[614] = par__Rfoot_M*(z[44]*z[594]+z[47]*z[595]) + par__Rtibia_M*(z[26]*
  z[593]+z[28]*z[592]) + par__Rfemur_M*(par__Rfemur_CMx*z[591]-par__Rfemur_CMy*
  z[9]);
  z[615] = par__Rfoot_M*(z[46]*z[594]+z[49]*z[595]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[592]-par__Rtibia_CMy*z[593]);
  z[616] = par__Rfoot_M*(par__Rfoot_CMx*z[595]-par__Rfoot_CMy*z[594]);
  z[617] = par__Rfemur_M*(z[5]*z[10]+z[5]*z[591]+z[6]*z[11]+z[7]*z[9]) + 
  par__Rfoot_M*(z[34]*z[45]+z[35]*z[48]+z[36]*z[594]+z[37]*z[595]) + par__Rtibia_M*(
  z[18]*z[27]+z[19]*z[29]+z[20]*z[593]+z[21]*z[592]);
  z[618] = par__Rfemur_M*(z[5]*z[9]+z[6]*z[10]+z[6]*z[591]+z[11]*z[99]) + 
  par__Rfoot_M*(z[34]*z[594]+z[35]*z[595]+z[45]*z[106]+z[48]*z[108]) + par__Rtibia_M*(
  z[18]*z[593]+z[19]*z[592]+z[27]*z[102]+z[29]*z[103]);
  z[619] = z[16]*z[100] + z[17]*z[101] - z[14]*z[16] - z[15]*z[17];
  z[620] = z[15]*z[16] + z[17]*z[100] - z[14]*z[17] - z[16]*z[101];
  z[621] = z[32]*z[620] + z[33]*z[619];
  z[622] = z[32]*z[619] - z[33]*z[620];
  z[623] = par__Rankle_y*z[19]*z[327] + par__Rankle_y*z[21]*z[328] + z[303]*
  z[34]*z[621] + z[303]*z[35]*z[622] + z[303]*z[45]*z[107] + z[303]*z[48]*
  z[109] + z[304]*z[36]*z[621] + z[304]*z[37]*z[622] + z[304]*z[45]*z[111] + 
  z[304]*z[48]*z[112] + z[305]*z[18]*z[620] + z[305]*z[19]*z[27] + z[305]*
  z[19]*z[619] + z[305]*z[29]*z[104] + z[306]*z[20]*z[620] + z[306]*z[21]*
  z[27] + z[306]*z[21]*z[619] + z[306]*z[29]*z[105] + z[5]*z[101]*z[327] + 
  z[7]*z[101]*z[328] + z[34]*z[327]*z[621] + z[35]*z[327]*z[622] + z[36]*
  z[328]*z[621] + z[37]*z[328]*z[622] + z[107]*z[322]*z[327] + z[109]*z[325]*
  z[327] + z[111]*z[322]*z[328] + z[112]*z[325]*z[328] - z[5]*z[100]*z[328] - 
  z[6]*z[100]*z[327];
  z[624] = z[16]*z[100]*q5p + z[17]*z[101]*q5p + z[17]*z[201]*q5p + z[16]*
  z[466] - z[14]*z[16]*q5p - z[15]*z[17]*q5p - z[16]*z[202]*q5p - z[17]*
  z[467];
  z[625] = z[14]*z[17]*q5p + z[16]*z[101]*q5p + z[16]*z[201]*q5p + z[17]*
  z[202]*q5p - z[15]*z[16]*q5p - z[17]*z[100]*q5p - z[16]*z[467] - z[17]*
  z[466];
  z[626] = z[32]*z[619]*q6p + z[32]*z[624] + z[33]*z[625] - z[33]*z[620]*q6p;
  z[627] = z[5]*z[17]*q5p - z[6]*z[16]*q5p - z[16]*z[446] - z[17]*z[450];
  z[628] = z[32]*z[104]*q6p + z[32]*z[458] + z[33]*z[627] - z[19]*z[33]*q6p;
  z[629] = z[7]*z[17]*q5p - z[5]*z[16]*q5p - z[16]*z[447] - z[17]*z[446];
  z[630] = z[32]*z[105]*q6p + z[32]*z[461] + z[33]*z[629] - z[21]*z[33]*q6p;
  z[631] = z[13]*z[16] + z[16]*z[100] + z[17]*z[101] - z[12]*z[17];
  z[632] = z[12]*z[16] + z[13]*z[17] + z[17]*z[100] - z[16]*z[101];
  z[633] = (z[12]*z[17]+z[16]*z[202]-2*z[16]*z[100]-2*z[17]*z[101]-z[13]*
  z[16]-z[17]*z[201])*q5p;
  z[634] = (z[12]*z[16]+z[13]*z[17]+2*z[17]*z[100]-2*z[16]*z[101]-z[16]*
  z[201]-z[17]*z[202])*q5p;
  z[635] = z[32]*z[631]*q6p - z[33]*z[632]*q6p - z[32]*z[633] - z[33]*z[634];
  z[636] = q1p*z[626] + q2p*z[628] + q3p*z[630] + q4p*z[635];
  z[637] = z[32]*z[631] - z[33]*z[632];
  z[638] = z[109]*q2p + z[112]*q3p + z[622]*q1p + z[637]*q4p;
  z[639] = z[32]*z[625] - z[32]*z[620]*q6p - z[33]*z[619]*q6p - z[33]*z[624];
  z[640] = z[32]*z[627] - z[19]*z[32]*q6p - z[33]*z[104]*q6p - z[33]*z[458];
  z[641] = z[32]*z[629] - z[21]*z[32]*q6p - z[33]*z[105]*q6p - z[33]*z[461];
  z[642] = z[33]*z[633] - z[32]*z[632]*q6p - z[33]*z[631]*q6p - z[32]*z[634];
  z[643] = q1p*z[639] + q2p*z[640] + q3p*z[641] + q4p*z[642];
  z[644] = z[32]*z[632] + z[33]*z[631];
  z[645] = z[107]*q2p + z[111]*q3p + z[621]*q1p + z[644]*q4p;
  z[646] = q1p*z[624] + q2p*z[458] + q3p*z[461] - q4p*z[633];
  z[647] = z[104]*q2p + z[105]*q3p + z[619]*q1p + z[631]*q4p;
  z[648] = q1p*z[625] + q2p*z[627] + q3p*z[629] - q4p*z[634];
  z[649] = z[19]*q2p + z[21]*q3p + z[620]*q1p + z[632]*q4p;
  z[650] = par__Rfoot_M*(z[621]*z[479]+z[622]*z[480]+z[48]*(z[643]+z[393]*
  z[645])+z[45]*(z[636]-z[393]*z[638])) + par__Rtibia_M*(z[619]*z[486]+z[620]*
  z[485]+z[29]*(z[648]+z[394]*z[649])+z[27]*(z[646]-z[394]*z[647]));
  z[651] = 2*par__Rfoot_M*(z[45]*z[621]+z[48]*z[622]) + 2*par__Rtibia_M*(
  z[27]*z[620]+z[29]*z[619]);
  z[652] = par__Rfoot_M*(z[44]*z[621]+z[45]*z[644]+z[47]*z[622]+z[48]*z[637]) + 
  par__Rtibia_M*(z[26]*z[620]+z[27]*z[632]+z[28]*z[619]+z[29]*z[631]);
  z[653] = par__Rfoot_M*(z[46]*z[621]+z[49]*z[622]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[619]-par__Rtibia_CMy*z[620]);
  z[654] = par__Rfoot_M*(par__Rfoot_CMx*z[622]-par__Rfoot_CMy*z[621]);
  z[655] = par__Rfoot_M*(z[36]*z[621]+z[37]*z[622]+z[45]*z[111]+z[48]*z[112]) + 
  par__Rtibia_M*(z[20]*z[620]+z[21]*z[27]+z[21]*z[619]+z[29]*z[105]);
  z[656] = par__Rfoot_M*(z[34]*z[621]+z[35]*z[622]+z[45]*z[107]+z[48]*z[109]) + 
  par__Rtibia_M*(z[18]*z[620]+z[19]*z[27]+z[19]*z[619]+z[29]*z[104]);
  z[657] = -z[25]*z[33] - z[31]*z[32];
  z[658] = z[303]*z[34]*z[42] + z[303]*z[35]*z[45] + z[303]*z[35]*z[657] + 
  z[303]*z[48]*z[110] + z[304]*z[36]*z[42] + z[304]*z[37]*z[45] + z[304]*
  z[37]*z[657] + z[304]*z[48]*z[113] + z[34]*z[42]*z[327] + z[35]*z[322]*
  z[327] + z[35]*z[327]*z[657] + z[36]*z[42]*z[328] + z[37]*z[322]*z[328] + 
  z[37]*z[328]*z[657] + z[110]*z[325]*z[327] + z[113]*z[325]*z[328];
  z[659] = par__Rankle_y*z[32]*q5p*q6p + q1p*z[476] + q2p*z[473] + q3p*z[474] + 
  q4p*z[475];
  z[660] = -z[24]*z[33] - z[30]*z[32];
  z[661] = z[110]*q2p + z[113]*q3p + z[657]*q1p + z[660]*q4p + par__Rankle_y*
  z[32]*q5p;
  z[662] = z[31]*z[33]*q6p - z[25]*z[32]*q6p - z[32]*z[469] - z[33]*z[468];
  z[663] = z[18]*z[33]*q6p - z[19]*z[32]*q6p - z[32]*z[457] - z[33]*z[458];
  z[664] = z[20]*z[33]*q6p - z[21]*z[32]*q6p - z[32]*z[460] - z[33]*z[461];
  z[665] = z[30]*z[33]*q6p + z[32]*z[464] + z[33]*z[463] - z[24]*z[32]*q6p;
  z[666] = q1p*z[662] + q2p*z[663] + q3p*z[664] + q4p*z[665] - par__Rankle_y*
  z[33]*q5p*q6p;
  z[667] = z[35]*q2p + z[37]*q3p + z[41]*q4p + z[42]*q1p + par__Rankle_y*
  z[33]*q5p;
  z[668] = par__Rfoot_M*(z[42]*z[479]+z[657]*z[480]+z[48]*(z[666]+z[393]*
  z[667])+z[45]*(z[659]-z[393]*z[661]));
  z[669] = par__Rfoot_M*(z[42]*z[45]+z[48]*z[657]);
  z[670] = par__Rfoot_M*(z[41]*z[45]+z[42]*z[44]+z[47]*z[657]+z[48]*z[660]);
  z[671] = par__Rfoot_M*(z[42]*z[46]+z[49]*z[657]+par__Rankle_y*z[32]*z[48]+
  par__Rankle_y*z[33]*z[45]);
  z[672] = par__Rfoot_M*(par__Rfoot_CMx*z[657]-par__Rfoot_CMy*z[42]);
  z[673] = par__Rfoot_M*(z[36]*z[42]+z[37]*z[45]+z[37]*z[657]+z[48]*z[113]);
  z[674] = par__Rfoot_M*(z[34]*z[42]+z[35]*z[45]+z[35]*z[657]+z[48]*z[110]);
  z[675] = par__Lhip_y*z[50] - par__Lhip_x*z[51];
  z[676] = z[63]*z[675] - z[56]*z[64];
  z[677] = z[56]*z[63] + z[64]*z[675];
  z[678] = z[79]*z[677] + z[80]*z[676];
  z[679] = z[79]*z[676] - z[80]*z[677];
  z[680] = par__Lankle_y*z[65]*z[345] + par__Lankle_y*z[117]*z[344] + z[293]*
  z[52]*z[56] + z[293]*z[53]*z[57] + z[293]*z[53]*z[675] + z[293]*z[58]*
  z[114] + z[294]*z[52]*z[57] + z[294]*z[52]*z[675] + z[294]*z[53]*z[58] + 
  z[294]*z[54]*z[56] + z[295]*z[81]*z[678] + z[295]*z[82]*z[679] + z[295]*
  z[92]*z[121] + z[295]*z[95]*z[123] + z[296]*z[81]*z[92] + z[296]*z[82]*
  z[95] + z[296]*z[83]*z[678] + z[296]*z[84]*z[679] + z[297]*z[65]*z[677] + 
  z[297]*z[66]*z[676] + z[297]*z[74]*z[117] + z[297]*z[76]*z[118] + z[298]*
  z[65]*z[74] + z[298]*z[66]*z[76] + z[298]*z[67]*z[677] + z[298]*z[68]*
  z[676] + z[52]*z[60]*z[345] + z[53]*z[60]*z[344] + z[81]*z[339]*z[345] + 
  z[81]*z[344]*z[678] + z[82]*z[342]*z[345] + z[82]*z[344]*z[679] + z[83]*
  z[345]*z[678] + z[84]*z[345]*z[679] + z[121]*z[339]*z[344] + z[123]*z[342]*
  z[344] - z[53]*z[59]*z[345] - z[59]*z[114]*z[344];
  z[681] = z[53]*q3p + z[114]*q2p + z[675]*q1p;
  z[682] = q2p*z[555] + q3p*z[401] - q1p*z[399];
  z[683] = z[52]*q3p + z[53]*q2p + z[56]*q1p;
  z[684] = z[63]*z[675]*q8p - z[56]*z[64]*q8p - z[63]*z[402] - z[64]*z[399];
  z[685] = z[64]*z[402] - z[56]*z[63]*q8p - z[64]*z[675]*q8p - z[63]*z[399];
  z[686] = z[79]*z[676]*q9p + z[79]*z[684] + z[80]*z[685] - z[80]*z[677]*q9p;
  z[687] = q1p*z[686] + q2p*z[560] + q3p*z[410];
  z[688] = z[82]*q3p + z[123]*q2p + z[679]*q1p;
  z[689] = z[79]*z[685] - z[79]*z[677]*q9p - z[80]*z[676]*q9p - z[80]*z[684];
  z[690] = q1p*z[689] + q2p*z[563] + q3p*z[424];
  z[691] = z[81]*q3p + z[121]*q2p + z[678]*q1p;
  z[692] = q1p*z[684] + q2p*z[558] + q3p*z[408];
  z[693] = z[66]*q3p + z[118]*q2p + z[676]*q1p;
  z[694] = q1p*z[685] + q2p*z[559] + q3p*z[409];
  z[695] = z[65]*q3p + z[117]*q2p + z[677]*q1p;
  z[696] = par__Lfemur_M*(z[56]*z[404]+z[675]*z[405]+z[58]*(z[682]+z[389]*
  z[683])+z[57]*(z[403]-z[389]*z[681])) + par__Lfoot_M*(z[678]*z[430]+z[679]*
  z[431]+z[95]*(z[690]+z[390]*z[691])+z[92]*(z[687]-z[390]*z[688])) + par__Ltibia_M*(
  z[676]*z[437]+z[677]*z[436]+z[76]*(z[694]+z[391]*z[695])+z[74]*(z[692]-
  z[391]*z[693]));
  z[697] = 2*par__Lfemur_M*(z[56]*z[57]+z[58]*z[675]) + 2*par__Lfoot_M*(z[92]*
  z[678]+z[95]*z[679]) + 2*par__Ltibia_M*(z[74]*z[677]+z[76]*z[676]);
  z[698] = par__Lfoot_M*(z[91]*z[678]+z[94]*z[679]) + par__Ltibia_M*(z[73]*
  z[677]+z[75]*z[676]) + par__Lfemur_M*(par__Lfemur_CMx*z[675]-par__Lfemur_CMy*
  z[56]);
  z[699] = par__Lfoot_M*(z[93]*z[678]+z[96]*z[679]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[676]-par__Ltibia_CMy*z[677]);
  z[700] = par__Lfoot_M*(par__Lfoot_CMx*z[679]-par__Lfoot_CMy*z[678]);
  z[701] = par__Lfemur_M*(z[52]*z[57]+z[52]*z[675]+z[53]*z[58]+z[54]*z[56]) + 
  par__Lfoot_M*(z[81]*z[92]+z[82]*z[95]+z[83]*z[678]+z[84]*z[679]) + par__Ltibia_M*(
  z[65]*z[74]+z[66]*z[76]+z[67]*z[677]+z[68]*z[676]);
  z[702] = par__Lfemur_M*(z[52]*z[56]+z[53]*z[57]+z[53]*z[675]+z[58]*z[114]) + 
  par__Lfoot_M*(z[81]*z[678]+z[82]*z[679]+z[92]*z[121]+z[95]*z[123]) + par__Ltibia_M*(
  z[65]*z[677]+z[66]*z[676]+z[74]*z[117]+z[76]*z[118]);
  z[703] = z[63]*z[115] + z[64]*z[116] - z[61]*z[63] - z[62]*z[64];
  z[704] = z[62]*z[63] + z[64]*z[115] - z[61]*z[64] - z[63]*z[116];
  z[705] = z[79]*z[704] + z[80]*z[703];
  z[706] = z[79]*z[703] - z[80]*z[704];
  z[707] = par__Lankle_y*z[66]*z[344] + par__Lankle_y*z[68]*z[345] + z[295]*
  z[81]*z[705] + z[295]*z[82]*z[706] + z[295]*z[92]*z[122] + z[295]*z[95]*
  z[124] + z[296]*z[83]*z[705] + z[296]*z[84]*z[706] + z[296]*z[92]*z[126] + 
  z[296]*z[95]*z[127] + z[297]*z[65]*z[704] + z[297]*z[66]*z[74] + z[297]*
  z[66]*z[703] + z[297]*z[76]*z[119] + z[298]*z[67]*z[704] + z[298]*z[68]*
  z[74] + z[298]*z[68]*z[703] + z[298]*z[76]*z[120] + z[52]*z[116]*z[344] + 
  z[54]*z[116]*z[345] + z[81]*z[344]*z[705] + z[82]*z[344]*z[706] + z[83]*
  z[345]*z[705] + z[84]*z[345]*z[706] + z[122]*z[339]*z[344] + z[124]*z[342]*
  z[344] + z[126]*z[339]*z[345] + z[127]*z[342]*z[345] - z[52]*z[115]*z[345] - 
  z[53]*z[115]*z[344];
  z[708] = z[63]*z[115]*q8p + z[64]*z[116]*q8p + z[64]*z[247]*q8p + z[63]*
  z[417] - z[61]*z[63]*q8p - z[62]*z[64]*q8p - z[63]*z[248]*q8p - z[64]*
  z[418];
  z[709] = z[61]*z[64]*q8p + z[63]*z[116]*q8p + z[63]*z[247]*q8p + z[64]*
  z[248]*q8p - z[62]*z[63]*q8p - z[64]*z[115]*q8p - z[63]*z[418] - z[64]*
  z[417];
  z[710] = z[79]*z[703]*q9p + z[79]*z[708] + z[80]*z[709] - z[80]*z[704]*q9p;
  z[711] = z[52]*z[64]*q8p - z[53]*z[63]*q8p - z[63]*z[397] - z[64]*z[401];
  z[712] = z[79]*z[119]*q9p + z[79]*z[409] + z[80]*z[711] - z[66]*z[80]*q9p;
  z[713] = z[54]*z[64]*q8p - z[52]*z[63]*q8p - z[63]*z[398] - z[64]*z[397];
  z[714] = z[79]*z[120]*q9p + z[79]*z[412] + z[80]*z[713] - z[68]*z[80]*q9p;
  z[715] = z[60]*z[63] + z[63]*z[115] + z[64]*z[116] - z[59]*z[64];
  z[716] = z[59]*z[63] + z[60]*z[64] + z[64]*z[115] - z[63]*z[116];
  z[717] = (z[59]*z[64]+z[63]*z[248]-2*z[63]*z[115]-2*z[64]*z[116]-z[60]*
  z[63]-z[64]*z[247])*q8p;
  z[718] = (z[59]*z[63]+z[60]*z[64]+2*z[64]*z[115]-2*z[63]*z[116]-z[63]*
  z[247]-z[64]*z[248])*q8p;
  z[719] = z[79]*z[715]*q9p - z[80]*z[716]*q9p - z[79]*z[717] - z[80]*z[718];
  z[720] = q1p*z[710] + q2p*z[712] + q3p*z[714] + q7p*z[719];
  z[721] = z[79]*z[715] - z[80]*z[716];
  z[722] = z[124]*q2p + z[127]*q3p + z[706]*q1p + z[721]*q7p;
  z[723] = z[79]*z[709] - z[79]*z[704]*q9p - z[80]*z[703]*q9p - z[80]*z[708];
  z[724] = z[79]*z[711] - z[66]*z[79]*q9p - z[80]*z[119]*q9p - z[80]*z[409];
  z[725] = z[79]*z[713] - z[68]*z[79]*q9p - z[80]*z[120]*q9p - z[80]*z[412];
  z[726] = z[80]*z[717] - z[79]*z[716]*q9p - z[80]*z[715]*q9p - z[79]*z[718];
  z[727] = q1p*z[723] + q2p*z[724] + q3p*z[725] + q7p*z[726];
  z[728] = z[79]*z[716] + z[80]*z[715];
  z[729] = z[122]*q2p + z[126]*q3p + z[705]*q1p + z[728]*q7p;
  z[730] = q1p*z[708] + q2p*z[409] + q3p*z[412] - q7p*z[717];
  z[731] = z[119]*q2p + z[120]*q3p + z[703]*q1p + z[715]*q7p;
  z[732] = q1p*z[709] + q2p*z[711] + q3p*z[713] - q7p*z[718];
  z[733] = z[66]*q2p + z[68]*q3p + z[704]*q1p + z[716]*q7p;
  z[734] = par__Lfoot_M*(z[705]*z[430]+z[706]*z[431]+z[95]*(z[727]+z[390]*
  z[729])+z[92]*(z[720]-z[390]*z[722])) + par__Ltibia_M*(z[703]*z[437]+z[704]*
  z[436]+z[76]*(z[732]+z[391]*z[733])+z[74]*(z[730]-z[391]*z[731]));
  z[735] = 2*par__Lfoot_M*(z[92]*z[705]+z[95]*z[706]) + 2*par__Ltibia_M*(
  z[74]*z[704]+z[76]*z[703]);
  z[736] = par__Lfoot_M*(z[91]*z[705]+z[92]*z[728]+z[94]*z[706]+z[95]*z[721]) + 
  par__Ltibia_M*(z[73]*z[704]+z[74]*z[716]+z[75]*z[703]+z[76]*z[715]);
  z[737] = par__Lfoot_M*(z[93]*z[705]+z[96]*z[706]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[703]-par__Ltibia_CMy*z[704]);
  z[738] = par__Lfoot_M*(par__Lfoot_CMx*z[706]-par__Lfoot_CMy*z[705]);
  z[739] = par__Lfoot_M*(z[83]*z[705]+z[84]*z[706]+z[92]*z[126]+z[95]*z[127]) + 
  par__Ltibia_M*(z[67]*z[704]+z[68]*z[74]+z[68]*z[703]+z[76]*z[120]);
  z[740] = par__Lfoot_M*(z[81]*z[705]+z[82]*z[706]+z[92]*z[122]+z[95]*z[124]) + 
  par__Ltibia_M*(z[65]*z[704]+z[66]*z[74]+z[66]*z[703]+z[76]*z[119]);
  z[741] = -z[72]*z[80] - z[78]*z[79];
  z[742] = z[295]*z[81]*z[89] + z[295]*z[82]*z[92] + z[295]*z[82]*z[741] + 
  z[295]*z[95]*z[125] + z[296]*z[83]*z[89] + z[296]*z[84]*z[92] + z[296]*
  z[84]*z[741] + z[296]*z[95]*z[128] + z[81]*z[89]*z[344] + z[82]*z[339]*
  z[344] + z[82]*z[344]*z[741] + z[83]*z[89]*z[345] + z[84]*z[339]*z[345] + 
  z[84]*z[345]*z[741] + z[125]*z[342]*z[344] + z[128]*z[342]*z[345];
  z[743] = par__Lankle_y*z[79]*q8p*q9p + q1p*z[427] + q2p*z[424] + q3p*z[425] + 
  q7p*z[426];
  z[744] = -z[71]*z[80] - z[77]*z[79];
  z[745] = z[125]*q2p + z[128]*q3p + z[741]*q1p + z[744]*q7p + par__Lankle_y*
  z[79]*q8p;
  z[746] = z[78]*z[80]*q9p - z[72]*z[79]*q9p - z[79]*z[420] - z[80]*z[419];
  z[747] = z[65]*z[80]*q9p - z[66]*z[79]*q9p - z[79]*z[408] - z[80]*z[409];
  z[748] = z[67]*z[80]*q9p - z[68]*z[79]*q9p - z[79]*z[411] - z[80]*z[412];
  z[749] = z[77]*z[80]*q9p + z[79]*z[415] + z[80]*z[414] - z[71]*z[79]*q9p;
  z[750] = q1p*z[746] + q2p*z[747] + q3p*z[748] + q7p*z[749] - par__Lankle_y*
  z[80]*q8p*q9p;
  z[751] = z[82]*q2p + z[84]*q3p + z[88]*q7p + z[89]*q1p + par__Lankle_y*
  z[80]*q8p;
  z[752] = par__Lfoot_M*(z[89]*z[430]+z[741]*z[431]+z[95]*(z[750]+z[390]*
  z[751])+z[92]*(z[743]-z[390]*z[745]));
  z[753] = par__Lfoot_M*(z[89]*z[92]+z[95]*z[741]);
  z[754] = par__Lfoot_M*(z[88]*z[92]+z[89]*z[91]+z[94]*z[741]+z[95]*z[744]);
  z[755] = par__Lfoot_M*(z[89]*z[93]+z[96]*z[741]+par__Lankle_y*z[79]*z[95]+
  par__Lankle_y*z[80]*z[92]);
  z[756] = par__Lfoot_M*(par__Lfoot_CMx*z[741]-par__Lfoot_CMy*z[89]);
  z[757] = par__Lfoot_M*(z[83]*z[89]+z[84]*z[92]+z[84]*z[741]+z[95]*z[128]);
  z[758] = par__Lfoot_M*(z[81]*z[89]+z[82]*z[92]+z[82]*z[741]+z[95]*z[125]);
  z[759] = z[294]*z[52]*z[114] + z[294]*pow(z[52],2) + z[294]*z[53]*z[54] + 
  z[294]*pow(z[53],2) + z[296]*pow(z[81],2) + z[296]*pow(z[82],2) + z[296]*
  z[83]*z[121] + z[296]*z[84]*z[123] + z[298]*pow(z[65],2) + z[298]*pow(z[66],
  2) + z[298]*z[67]*z[117] + z[298]*z[68]*z[118] + z[302]*z[5]*z[99] + z[302]*
  pow(z[5],2) + z[302]*z[6]*z[7] + z[302]*pow(z[6],2) + z[304]*pow(z[34],2) + 
  z[304]*pow(z[35],2) + z[304]*z[36]*z[106] + z[304]*z[37]*z[108] + z[306]*
  pow(z[18],2) + z[306]*pow(z[19],2) + z[306]*z[20]*z[102] + z[306]*z[21]*
  z[103] + z[36]*z[106]*z[328] + z[37]*z[108]*z[328] + z[83]*z[121]*z[345] + 
  z[84]*z[123]*z[345] + z[328]*pow(z[34],2) + z[328]*pow(z[35],2) + z[345]*
  pow(z[81],2) + z[345]*pow(z[82],2) + 2*z[293]*z[52]*z[53] + 2*z[293]*z[53]*
  z[114] + 2*z[295]*z[81]*z[121] + 2*z[295]*z[82]*z[123] + 2*z[297]*z[65]*
  z[117] + 2*z[297]*z[66]*z[118] + 2*z[301]*z[5]*z[6] + 2*z[301]*z[6]*z[99] + 
  2*z[303]*z[34]*z[106] + 2*z[303]*z[35]*z[108] + 2*z[305]*z[18]*z[102] + 2*
  z[305]*z[19]*z[103] + 2*z[34]*z[106]*z[327] + 2*z[35]*z[108]*z[327] + 2*
  z[81]*z[121]*z[344] + 2*z[82]*z[123]*z[344] + z[1]*z[312] - z[1]*z[549] - 
  z[2]*z[311] - z[2]*z[548] - z[309]*pow(z[1],2)*z[547] - z[309]*pow(z[2],2)*
  z[547];
  z[760] = par__Lfemur_M*(z[53]*z[404]+z[114]*z[405]+z[53]*(z[556]+z[389]*
  z[557])+z[52]*(z[553]-z[389]*z[554])) + par__Lfoot_M*(z[121]*z[430]+z[123]*
  z[431]+z[82]*(z[564]+z[390]*z[565])+z[81]*(z[561]-z[390]*z[562])) + par__Ltibia_M*(
  z[117]*z[436]+z[118]*z[437]+z[66]*(z[568]+z[391]*z[569])+z[65]*(z[566]-
  z[391]*z[567])) + par__Rfemur_M*(z[6]*z[453]+z[99]*z[454]+z[6]*(z[573]+
  z[392]*z[574])+z[5]*(z[570]-z[392]*z[571])) + par__Rfoot_M*(z[106]*z[479]+
  z[108]*z[480]+z[35]*(z[581]+z[393]*z[582])+z[34]*(z[578]-z[393]*z[579])) + 
  par__Rtibia_M*(z[102]*z[485]+z[103]*z[486]+z[19]*(z[585]+z[394]*z[586])+
  z[18]*(z[583]-z[394]*z[584])) - par__pelvis_M*(z[1]*z[443]+z[2]*z[442]+z[1]*(
  z[441]+q1p*z[587])-z[2]*(z[440]-q1p*z[588])) - par__torso_M*(z[1]*z[490]+
  z[2]*z[489]+z[1]*(z[441]+q1p*z[587])-z[2]*(z[440]-q1p*z[588]));
  z[761] = par__Lfoot_M*(z[91]*z[121]+z[94]*z[123]) + par__Ltibia_M*(z[73]*
  z[117]+z[75]*z[118]) + par__Lfemur_M*(par__Lfemur_CMx*z[114]-par__Lfemur_CMy*
  z[53]);
  z[762] = par__Lfemur_M*z[52]*(z[52]+z[114]) + par__Lfemur_M*z[53]*(z[53]+
  z[54]) + par__Rfemur_M*z[5]*(z[5]+z[99]) + par__Rfemur_M*z[6]*(z[6]+z[7]) + 
  par__Lfoot_M*(pow(z[81],2)+pow(z[82],2)+z[83]*z[121]+z[84]*z[123]) + par__Ltibia_M*(
  pow(z[65],2)+pow(z[66],2)+z[67]*z[117]+z[68]*z[118]) + par__Rfoot_M*(pow(
  z[34],2)+pow(z[35],2)+z[36]*z[106]+z[37]*z[108]) + par__Rtibia_M*(pow(z[18],
  2)+pow(z[19],2)+z[20]*z[102]+z[21]*z[103]);
  z[763] = 2*par__Lfemur_M*z[53]*(z[52]+z[114]) + 2*par__Rfemur_M*z[6]*(z[5]+
  z[99]) + 2*par__Lfoot_M*(z[81]*z[121]+z[82]*z[123]) + 2*par__Ltibia_M*(
  z[65]*z[117]+z[66]*z[118]) + 2*par__Rfoot_M*(z[34]*z[106]+z[35]*z[108]) + 2*
  par__Rtibia_M*(z[18]*z[102]+z[19]*z[103]);
  z[764] = par__Lfoot_M*(par__Lfoot_CMx*z[123]-par__Lfoot_CMy*z[121]);
  z[765] = par__Lfoot_M*(z[93]*z[121]+z[96]*z[123]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[118]-par__Ltibia_CMy*z[117]);
  z[766] = par__Rfoot_M*(z[44]*z[106]+z[47]*z[108]) + par__Rtibia_M*(z[26]*
  z[102]+z[28]*z[103]) + par__Rfemur_M*(par__Rfemur_CMx*z[99]-par__Rfemur_CMy*
  z[6]);
  z[767] = par__Rfoot_M*(par__Rfoot_CMx*z[108]-par__Rfoot_CMy*z[106]);
  z[768] = par__Rfoot_M*(z[46]*z[106]+z[49]*z[108]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[103]-par__Rtibia_CMy*z[102]);
  z[769] = z[302]*z[5]*z[99] + z[302]*pow(z[5],2) + z[302]*z[6]*z[7] + z[302]*
  pow(z[6],2) + z[304]*pow(z[34],2) + z[304]*pow(z[35],2) + z[304]*z[36]*
  z[106] + z[304]*z[37]*z[108] + z[306]*pow(z[18],2) + z[306]*pow(z[19],2) + 
  z[306]*z[20]*z[102] + z[306]*z[21]*z[103] + z[36]*z[106]*z[328] + z[37]*
  z[108]*z[328] + z[328]*pow(z[34],2) + z[328]*pow(z[35],2) + 2*z[301]*z[5]*
  z[6] + 2*z[301]*z[6]*z[99] + 2*z[303]*z[34]*z[106] + 2*z[303]*z[35]*z[108] + 
  2*z[305]*z[18]*z[102] + 2*z[305]*z[19]*z[103] + 2*z[34]*z[106]*z[327] + 2*
  z[35]*z[108]*z[327];
  z[770] = par__Rfemur_M*(z[6]*z[453]+z[99]*z[454]+z[6]*(z[598]+z[392]*z[599])+
  z[5]*(z[452]-z[392]*z[597])) + par__Rfoot_M*(z[106]*z[479]+z[108]*z[480]+
  z[35]*(z[606]+z[393]*z[607])+z[34]*(z[603]-z[393]*z[604])) + par__Rtibia_M*(
  z[102]*z[485]+z[103]*z[486]+z[19]*(z[610]+z[394]*z[611])+z[18]*(z[608]-
  z[394]*z[609]));
  z[771] = par__Rfemur_M*z[5]*(z[5]+z[99]) + par__Rfemur_M*z[6]*(z[6]+z[7]) + 
  par__Rfoot_M*(pow(z[34],2)+pow(z[35],2)+z[36]*z[106]+z[37]*z[108]) + par__Rtibia_M*(
  pow(z[18],2)+pow(z[19],2)+z[20]*z[102]+z[21]*z[103]);
  z[772] = 2*par__Rfemur_M*z[6]*(z[5]+z[99]) + 2*par__Rfoot_M*(z[34]*z[106]+
  z[35]*z[108]) + 2*par__Rtibia_M*(z[18]*z[102]+z[19]*z[103]);
  z[773] = z[304]*z[34]*z[111] + z[304]*z[35]*z[112] + z[304]*z[36]*z[107] + 
  z[304]*z[37]*z[109] + z[306]*z[18]*z[21] + z[306]*z[19]*z[20] + z[306]*
  z[19]*z[105] + z[306]*z[21]*z[104] + z[34]*z[111]*z[328] + z[35]*z[112]*
  z[328] + z[36]*z[107]*z[328] + z[37]*z[109]*z[328] + 2*z[303]*z[34]*z[107] + 
  2*z[303]*z[35]*z[109] + 2*z[305]*z[18]*z[19] + 2*z[305]*z[19]*z[104] + 2*
  z[34]*z[107]*z[327] + 2*z[35]*z[109]*z[327];
  z[774] = par__Rfoot_M*(z[107]*z[479]+z[109]*z[480]+z[35]*(z[643]+z[393]*
  z[645])+z[34]*(z[636]-z[393]*z[638])) + par__Rtibia_M*(z[19]*z[485]+z[104]*
  z[486]+z[19]*(z[648]+z[394]*z[649])+z[18]*(z[646]-z[394]*z[647]));
  z[775] = par__Rfoot_M*(z[34]*z[111]+z[35]*z[112]+z[36]*z[107]+z[37]*z[109]) + 
  par__Rtibia_M*(z[18]*z[21]+z[19]*z[20]+z[19]*z[105]+z[21]*z[104]);
  z[776] = 2*par__Rtibia_M*z[19]*(z[18]+z[104]) + 2*par__Rfoot_M*(z[34]*
  z[107]+z[35]*z[109]);
  z[777] = par__Rfoot_M*(z[34]*z[644]+z[35]*z[637]+z[44]*z[107]+z[47]*z[109]) + 
  par__Rtibia_M*(z[18]*z[632]+z[19]*z[26]+z[19]*z[631]+z[28]*z[104]);
  z[778] = par__Rfoot_M*(par__Rfoot_CMx*z[109]-par__Rfoot_CMy*z[107]);
  z[779] = par__Rfoot_M*(z[46]*z[107]+z[49]*z[109]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[104]-par__Rtibia_CMy*z[19]);
  z[780] = z[304]*z[34]*z[37] + z[304]*z[35]*z[36] + z[304]*z[35]*z[113] + 
  z[304]*z[37]*z[110] + z[34]*z[37]*z[328] + z[35]*z[36]*z[328] + z[35]*
  z[113]*z[328] + z[37]*z[110]*z[328] + 2*z[303]*z[34]*z[35] + 2*z[303]*z[35]*
  z[110] + 2*z[34]*z[35]*z[327] + 2*z[35]*z[110]*z[327];
  z[781] = par__Rfoot_M*(z[35]*z[479]+z[110]*z[480]+z[35]*(z[666]+z[393]*
  z[667])+z[34]*(z[659]-z[393]*z[661]));
  z[782] = par__Rfoot_M*(z[34]*z[37]+z[35]*z[36]+z[35]*z[113]+z[37]*z[110]);
  z[783] = par__Rfoot_M*z[35]*(z[34]+z[110]);
  z[784] = par__Rfoot_M*(z[34]*z[41]+z[35]*z[44]+z[35]*z[660]+z[47]*z[110]);
  z[785] = par__Rfoot_M*(par__Rfoot_CMx*z[110]-par__Rfoot_CMy*z[35]);
  z[786] = par__Rfoot_M*(z[35]*z[46]+z[49]*z[110]+par__Rankle_y*z[32]*z[35]+
  par__Rankle_y*z[33]*z[34]);
  z[787] = z[294]*z[52]*z[114] + z[294]*pow(z[52],2) + z[294]*z[53]*z[54] + 
  z[294]*pow(z[53],2) + z[296]*pow(z[81],2) + z[296]*pow(z[82],2) + z[296]*
  z[83]*z[121] + z[296]*z[84]*z[123] + z[298]*pow(z[65],2) + z[298]*pow(z[66],
  2) + z[298]*z[67]*z[117] + z[298]*z[68]*z[118] + z[83]*z[121]*z[345] + 
  z[84]*z[123]*z[345] + z[345]*pow(z[81],2) + z[345]*pow(z[82],2) + 2*z[293]*
  z[52]*z[53] + 2*z[293]*z[53]*z[114] + 2*z[295]*z[81]*z[121] + 2*z[295]*
  z[82]*z[123] + 2*z[297]*z[65]*z[117] + 2*z[297]*z[66]*z[118] + 2*z[81]*
  z[121]*z[344] + 2*z[82]*z[123]*z[344];
  z[788] = par__Lfemur_M*(z[53]*z[404]+z[114]*z[405]+z[53]*(z[682]+z[389]*
  z[683])+z[52]*(z[403]-z[389]*z[681])) + par__Lfoot_M*(z[121]*z[430]+z[123]*
  z[431]+z[82]*(z[690]+z[390]*z[691])+z[81]*(z[687]-z[390]*z[688])) + par__Ltibia_M*(
  z[117]*z[436]+z[118]*z[437]+z[66]*(z[694]+z[391]*z[695])+z[65]*(z[692]-
  z[391]*z[693]));
  z[789] = par__Lfemur_M*z[52]*(z[52]+z[114]) + par__Lfemur_M*z[53]*(z[53]+
  z[54]) + par__Lfoot_M*(pow(z[81],2)+pow(z[82],2)+z[83]*z[121]+z[84]*z[123]) + 
  par__Ltibia_M*(pow(z[65],2)+pow(z[66],2)+z[67]*z[117]+z[68]*z[118]);
  z[790] = 2*par__Lfemur_M*z[53]*(z[52]+z[114]) + 2*par__Lfoot_M*(z[81]*
  z[121]+z[82]*z[123]) + 2*par__Ltibia_M*(z[65]*z[117]+z[66]*z[118]);
  z[791] = z[296]*z[81]*z[126] + z[296]*z[82]*z[127] + z[296]*z[83]*z[122] + 
  z[296]*z[84]*z[124] + z[298]*z[65]*z[68] + z[298]*z[66]*z[67] + z[298]*
  z[66]*z[120] + z[298]*z[68]*z[119] + z[81]*z[126]*z[345] + z[82]*z[127]*
  z[345] + z[83]*z[122]*z[345] + z[84]*z[124]*z[345] + 2*z[295]*z[81]*z[122] + 
  2*z[295]*z[82]*z[124] + 2*z[297]*z[65]*z[66] + 2*z[297]*z[66]*z[119] + 2*
  z[81]*z[122]*z[344] + 2*z[82]*z[124]*z[344];
  z[792] = par__Lfoot_M*(z[122]*z[430]+z[124]*z[431]+z[82]*(z[727]+z[390]*
  z[729])+z[81]*(z[720]-z[390]*z[722])) + par__Ltibia_M*(z[66]*z[436]+z[119]*
  z[437]+z[66]*(z[732]+z[391]*z[733])+z[65]*(z[730]-z[391]*z[731]));
  z[793] = par__Lfoot_M*(z[81]*z[728]+z[82]*z[721]+z[91]*z[122]+z[94]*z[124]) + 
  par__Ltibia_M*(z[65]*z[716]+z[66]*z[73]+z[66]*z[715]+z[75]*z[119]);
  z[794] = par__Lfoot_M*(z[81]*z[126]+z[82]*z[127]+z[83]*z[122]+z[84]*z[124]) + 
  par__Ltibia_M*(z[65]*z[68]+z[66]*z[67]+z[66]*z[120]+z[68]*z[119]);
  z[795] = 2*par__Ltibia_M*z[66]*(z[65]+z[119]) + 2*par__Lfoot_M*(z[81]*
  z[122]+z[82]*z[124]);
  z[796] = par__Lfoot_M*(par__Lfoot_CMx*z[124]-par__Lfoot_CMy*z[122]);
  z[797] = par__Lfoot_M*(z[93]*z[122]+z[96]*z[124]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[119]-par__Ltibia_CMy*z[66]);
  z[798] = z[296]*z[81]*z[84] + z[296]*z[82]*z[83] + z[296]*z[82]*z[128] + 
  z[296]*z[84]*z[125] + z[81]*z[84]*z[345] + z[82]*z[83]*z[345] + z[82]*
  z[128]*z[345] + z[84]*z[125]*z[345] + 2*z[295]*z[81]*z[82] + 2*z[295]*z[82]*
  z[125] + 2*z[81]*z[82]*z[344] + 2*z[82]*z[125]*z[344];
  z[799] = par__Lfoot_M*(z[82]*z[430]+z[125]*z[431]+z[82]*(z[750]+z[390]*
  z[751])+z[81]*(z[743]-z[390]*z[745]));
  z[800] = par__Lfoot_M*(z[81]*z[88]+z[82]*z[91]+z[82]*z[744]+z[94]*z[125]);
  z[801] = par__Lfoot_M*(z[81]*z[84]+z[82]*z[83]+z[82]*z[128]+z[84]*z[125]);
  z[802] = par__Lfoot_M*z[82]*(z[81]+z[125]);
  z[803] = par__Lfoot_M*(par__Lfoot_CMx*z[125]-par__Lfoot_CMy*z[82]);
  z[804] = par__Lfoot_M*(z[82]*z[93]+z[96]*z[125]+par__Lankle_y*z[79]*z[82]+
  par__Lankle_y*z[80]*z[81]);
  z[805] = z[293]*z[52]*z[114] + z[293]*pow(z[52],2) + z[293]*z[53]*z[54] + 
  z[293]*pow(z[53],2) + z[295]*pow(z[81],2) + z[295]*pow(z[82],2) + z[295]*
  z[83]*z[121] + z[295]*z[84]*z[123] + z[297]*pow(z[65],2) + z[297]*pow(z[66],
  2) + z[297]*z[67]*z[117] + z[297]*z[68]*z[118] + z[301]*z[5]*z[99] + z[301]*
  pow(z[5],2) + z[301]*z[6]*z[7] + z[301]*pow(z[6],2) + z[303]*pow(z[34],2) + 
  z[303]*pow(z[35],2) + z[303]*z[36]*z[106] + z[303]*z[37]*z[108] + z[305]*
  pow(z[18],2) + z[305]*pow(z[19],2) + z[305]*z[20]*z[102] + z[305]*z[21]*
  z[103] + z[36]*z[106]*z[327] + z[37]*z[108]*z[327] + z[83]*z[121]*z[344] + 
  z[84]*z[123]*z[344] + z[327]*pow(z[34],2) + z[327]*pow(z[35],2) + z[344]*
  pow(z[81],2) + z[344]*pow(z[82],2) + 2*z[294]*z[52]*z[53] + 2*z[294]*z[52]*
  z[54] + 2*z[296]*z[81]*z[83] + 2*z[296]*z[82]*z[84] + 2*z[298]*z[65]*z[67] + 
  2*z[298]*z[66]*z[68] + 2*z[302]*z[5]*z[6] + 2*z[302]*z[5]*z[7] + 2*z[304]*
  z[34]*z[36] + 2*z[304]*z[35]*z[37] + 2*z[306]*z[18]*z[20] + 2*z[306]*z[19]*
  z[21] + 2*z[34]*z[36]*z[328] + 2*z[35]*z[37]*z[328] + 2*z[81]*z[83]*z[345] + 
  2*z[82]*z[84]*z[345] + z[1]*z[311] + z[1]*z[548] + z[2]*z[312] - z[2]*
  z[549];
  z[806] = par__Lfemur_M*(z[52]*z[404]+z[53]*z[405]+z[52]*(z[556]+z[389]*
  z[557])+z[54]*(z[553]-z[389]*z[554])) + par__Lfoot_M*(z[81]*z[430]+z[82]*
  z[431]+z[84]*(z[564]+z[390]*z[565])+z[83]*(z[561]-z[390]*z[562])) + par__Ltibia_M*(
  z[65]*z[436]+z[66]*z[437]+z[68]*(z[568]+z[391]*z[569])+z[67]*(z[566]-z[391]*
  z[567])) + par__Rfemur_M*(z[5]*z[453]+z[6]*z[454]+z[5]*(z[573]+z[392]*
  z[574])+z[7]*(z[570]-z[392]*z[571])) + par__Rfoot_M*(z[34]*z[479]+z[35]*
  z[480]+z[37]*(z[581]+z[393]*z[582])+z[36]*(z[578]-z[393]*z[579])) + par__Rtibia_M*(
  z[18]*z[485]+z[19]*z[486]+z[21]*(z[585]+z[394]*z[586])+z[20]*(z[583]-z[394]*
  z[584])) + par__pelvis_M*(z[1]*z[442]-z[2]*z[443]-z[2]*(z[441]+q1p*z[587])-
  z[1]*(z[440]-q1p*z[588])) + par__torso_M*(z[1]*z[489]-z[2]*z[490]-z[2]*(
  z[441]+q1p*z[587])-z[1]*(z[440]-q1p*z[588]));
  z[807] = 2*par__Lfemur_M*z[52]*(z[53]+z[54]) + 2*par__Rfemur_M*z[5]*(z[6]+
  z[7]) + 2*par__Lfoot_M*(z[81]*z[83]+z[82]*z[84]) + 2*par__Ltibia_M*(z[65]*
  z[67]+z[66]*z[68]) + 2*par__Rfoot_M*(z[34]*z[36]+z[35]*z[37]) + 2*par__Rtibia_M*(
  z[18]*z[20]+z[19]*z[21]);
  z[808] = z[301]*z[5]*z[99] + z[301]*pow(z[5],2) + z[301]*z[6]*z[7] + z[301]*
  pow(z[6],2) + z[303]*pow(z[34],2) + z[303]*pow(z[35],2) + z[303]*z[36]*
  z[106] + z[303]*z[37]*z[108] + z[305]*pow(z[18],2) + z[305]*pow(z[19],2) + 
  z[305]*z[20]*z[102] + z[305]*z[21]*z[103] + z[36]*z[106]*z[327] + z[37]*
  z[108]*z[327] + z[327]*pow(z[34],2) + z[327]*pow(z[35],2) + 2*z[302]*z[5]*
  z[6] + 2*z[302]*z[5]*z[7] + 2*z[304]*z[34]*z[36] + 2*z[304]*z[35]*z[37] + 2*
  z[306]*z[18]*z[20] + 2*z[306]*z[19]*z[21] + 2*z[34]*z[36]*z[328] + 2*z[35]*
  z[37]*z[328];
  z[809] = par__Rfemur_M*(z[5]*z[453]+z[6]*z[454]+z[5]*(z[598]+z[392]*z[599])+
  z[7]*(z[452]-z[392]*z[597])) + par__Rfoot_M*(z[34]*z[479]+z[35]*z[480]+
  z[37]*(z[606]+z[393]*z[607])+z[36]*(z[603]-z[393]*z[604])) + par__Rtibia_M*(
  z[18]*z[485]+z[19]*z[486]+z[21]*(z[610]+z[394]*z[611])+z[20]*(z[608]-z[394]*
  z[609]));
  z[810] = 2*par__Rfemur_M*z[5]*(z[6]+z[7]) + 2*par__Rfoot_M*(z[34]*z[36]+
  z[35]*z[37]) + 2*par__Rtibia_M*(z[18]*z[20]+z[19]*z[21]);
  z[811] = z[303]*z[34]*z[111] + z[303]*z[35]*z[112] + z[303]*z[36]*z[107] + 
  z[303]*z[37]*z[109] + z[305]*z[18]*z[21] + z[305]*z[19]*z[20] + z[305]*
  z[19]*z[105] + z[305]*z[21]*z[104] + z[34]*z[111]*z[327] + z[35]*z[112]*
  z[327] + z[36]*z[107]*z[327] + z[37]*z[109]*z[327] + 2*z[304]*z[36]*z[111] + 
  2*z[304]*z[37]*z[112] + 2*z[306]*z[20]*z[21] + 2*z[306]*z[21]*z[105] + 2*
  z[36]*z[111]*z[328] + 2*z[37]*z[112]*z[328];
  z[812] = par__Rfoot_M*(z[111]*z[479]+z[112]*z[480]+z[37]*(z[643]+z[393]*
  z[645])+z[36]*(z[636]-z[393]*z[638])) + par__Rtibia_M*(z[21]*z[485]+z[105]*
  z[486]+z[21]*(z[648]+z[394]*z[649])+z[20]*(z[646]-z[394]*z[647]));
  z[813] = 2*par__Rtibia_M*z[21]*(z[20]+z[105]) + 2*par__Rfoot_M*(z[36]*
  z[111]+z[37]*z[112]);
  z[814] = par__Rfoot_M*(z[36]*z[644]+z[37]*z[637]+z[44]*z[111]+z[47]*z[112]) + 
  par__Rtibia_M*(z[20]*z[632]+z[21]*z[26]+z[21]*z[631]+z[28]*z[105]);
  z[815] = par__Rfoot_M*(par__Rfoot_CMx*z[112]-par__Rfoot_CMy*z[111]);
  z[816] = par__Rfoot_M*(z[46]*z[111]+z[49]*z[112]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[105]-par__Rtibia_CMy*z[21]);
  z[817] = z[303]*z[34]*z[37] + z[303]*z[35]*z[36] + z[303]*z[35]*z[113] + 
  z[303]*z[37]*z[110] + z[34]*z[37]*z[327] + z[35]*z[36]*z[327] + z[35]*
  z[113]*z[327] + z[37]*z[110]*z[327] + 2*z[304]*z[36]*z[37] + 2*z[304]*z[37]*
  z[113] + 2*z[36]*z[37]*z[328] + 2*z[37]*z[113]*z[328];
  z[818] = par__Rfoot_M*(z[37]*z[479]+z[113]*z[480]+z[37]*(z[666]+z[393]*
  z[667])+z[36]*(z[659]-z[393]*z[661]));
  z[819] = par__Rfoot_M*z[37]*(z[36]+z[113]);
  z[820] = par__Rfoot_M*(z[36]*z[41]+z[37]*z[44]+z[37]*z[660]+z[47]*z[113]);
  z[821] = par__Rfoot_M*(par__Rfoot_CMx*z[113]-par__Rfoot_CMy*z[37]);
  z[822] = par__Rfoot_M*(z[37]*z[46]+z[49]*z[113]+par__Rankle_y*z[32]*z[37]+
  par__Rankle_y*z[33]*z[36]);
  z[823] = z[293]*z[52]*z[114] + z[293]*pow(z[52],2) + z[293]*z[53]*z[54] + 
  z[293]*pow(z[53],2) + z[295]*pow(z[81],2) + z[295]*pow(z[82],2) + z[295]*
  z[83]*z[121] + z[295]*z[84]*z[123] + z[297]*pow(z[65],2) + z[297]*pow(z[66],
  2) + z[297]*z[67]*z[117] + z[297]*z[68]*z[118] + z[83]*z[121]*z[344] + 
  z[84]*z[123]*z[344] + z[344]*pow(z[81],2) + z[344]*pow(z[82],2) + 2*z[294]*
  z[52]*z[53] + 2*z[294]*z[52]*z[54] + 2*z[296]*z[81]*z[83] + 2*z[296]*z[82]*
  z[84] + 2*z[298]*z[65]*z[67] + 2*z[298]*z[66]*z[68] + 2*z[81]*z[83]*z[345] + 
  2*z[82]*z[84]*z[345];
  z[824] = par__Lfemur_M*(z[52]*z[404]+z[53]*z[405]+z[52]*(z[682]+z[389]*
  z[683])+z[54]*(z[403]-z[389]*z[681])) + par__Lfoot_M*(z[81]*z[430]+z[82]*
  z[431]+z[84]*(z[690]+z[390]*z[691])+z[83]*(z[687]-z[390]*z[688])) + par__Ltibia_M*(
  z[65]*z[436]+z[66]*z[437]+z[68]*(z[694]+z[391]*z[695])+z[67]*(z[692]-z[391]*
  z[693]));
  z[825] = 2*par__Lfemur_M*z[52]*(z[53]+z[54]) + 2*par__Lfoot_M*(z[81]*z[83]+
  z[82]*z[84]) + 2*par__Ltibia_M*(z[65]*z[67]+z[66]*z[68]);
  z[826] = z[295]*z[81]*z[126] + z[295]*z[82]*z[127] + z[295]*z[83]*z[122] + 
  z[295]*z[84]*z[124] + z[297]*z[65]*z[68] + z[297]*z[66]*z[67] + z[297]*
  z[66]*z[120] + z[297]*z[68]*z[119] + z[81]*z[126]*z[344] + z[82]*z[127]*
  z[344] + z[83]*z[122]*z[344] + z[84]*z[124]*z[344] + 2*z[296]*z[83]*z[126] + 
  2*z[296]*z[84]*z[127] + 2*z[298]*z[67]*z[68] + 2*z[298]*z[68]*z[120] + 2*
  z[83]*z[126]*z[345] + 2*z[84]*z[127]*z[345];
  z[827] = par__Lfoot_M*(z[126]*z[430]+z[127]*z[431]+z[84]*(z[727]+z[390]*
  z[729])+z[83]*(z[720]-z[390]*z[722])) + par__Ltibia_M*(z[68]*z[436]+z[120]*
  z[437]+z[68]*(z[732]+z[391]*z[733])+z[67]*(z[730]-z[391]*z[731]));
  z[828] = par__Lfoot_M*(z[83]*z[728]+z[84]*z[721]+z[91]*z[126]+z[94]*z[127]) + 
  par__Ltibia_M*(z[67]*z[716]+z[68]*z[73]+z[68]*z[715]+z[75]*z[120]);
  z[829] = 2*par__Ltibia_M*z[68]*(z[67]+z[120]) + 2*par__Lfoot_M*(z[83]*
  z[126]+z[84]*z[127]);
  z[830] = par__Lfoot_M*(par__Lfoot_CMx*z[127]-par__Lfoot_CMy*z[126]);
  z[831] = par__Lfoot_M*(z[93]*z[126]+z[96]*z[127]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[120]-par__Ltibia_CMy*z[68]);
  z[832] = z[295]*z[81]*z[84] + z[295]*z[82]*z[83] + z[295]*z[82]*z[128] + 
  z[295]*z[84]*z[125] + z[81]*z[84]*z[344] + z[82]*z[83]*z[344] + z[82]*
  z[128]*z[344] + z[84]*z[125]*z[344] + 2*z[296]*z[83]*z[84] + 2*z[296]*z[84]*
  z[128] + 2*z[83]*z[84]*z[345] + 2*z[84]*z[128]*z[345];
  z[833] = par__Lfoot_M*(z[84]*z[430]+z[128]*z[431]+z[84]*(z[750]+z[390]*
  z[751])+z[83]*(z[743]-z[390]*z[745]));
  z[834] = par__Lfoot_M*(z[83]*z[88]+z[84]*z[91]+z[84]*z[744]+z[94]*z[128]);
  z[835] = par__Lfoot_M*z[84]*(z[83]+z[128]);
  z[836] = par__Lfoot_M*(par__Lfoot_CMx*z[128]-par__Lfoot_CMy*z[84]);
  z[837] = par__Lfoot_M*(z[84]*z[93]+z[96]*z[128]+par__Lankle_y*z[79]*z[84]+
  par__Lankle_y*z[80]*z[83]);
  z[838] = z[359]*z[99] + z[360]*z[6] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + par__Rhip_x*z[1]*z[327] + par__Rhip_x*z[2]*z[328] + par__Rhip_y*
  z[1]*z[328] + z[303]*z[44]*z[106] + z[303]*z[47]*z[108] + z[304]*z[34]*
  z[44] + z[304]*z[35]*z[47] + z[305]*z[26]*z[102] + z[305]*z[28]*z[103] + 
  z[306]*z[18]*z[26] + z[306]*z[19]*z[28] + z[5]*z[13]*z[328] + z[6]*z[13]*
  z[327] + z[34]*z[321]*z[328] + z[35]*z[324]*z[328] + z[106]*z[321]*z[327] + 
  z[108]*z[324]*z[327] - z[361]*z[6] - z[362]*z[5] - par__Rhip_y*z[2]*z[327] - 
  z[6]*z[12]*z[328] - z[12]*z[99]*z[327];
  z[839] = par__Rfoot_M*(z[47]*(z[581]+z[393]*z[582])+z[44]*(z[578]-z[393]*
  z[579])) + par__Rtibia_M*(z[28]*(z[585]+z[394]*z[586])+z[26]*(z[583]-z[394]*
  z[584])) + par__Rfemur_M*(par__Rfemur_CMx*(z[573]+z[392]*z[574])-par__Rfemur_CMy*(
  z[570]-z[392]*z[571]));
  z[840] = z[359]*z[99] + z[360]*z[6] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + z[303]*z[44]*z[106] + z[303]*z[47]*z[108] + z[304]*z[34]*
  z[44] + z[304]*z[35]*z[47] + z[305]*z[26]*z[102] + z[305]*z[28]*z[103] + 
  z[306]*z[18]*z[26] + z[306]*z[19]*z[28] + z[5]*z[13]*z[328] + z[6]*z[13]*
  z[327] + z[34]*z[321]*z[328] + z[35]*z[324]*z[328] + z[106]*z[321]*z[327] + 
  z[108]*z[324]*z[327] - z[361]*z[6] - z[362]*z[5] - z[6]*z[12]*z[328] - 
  z[12]*z[99]*z[327];
  z[841] = par__Rfoot_M*(z[47]*(z[606]+z[393]*z[607])+z[44]*(z[603]-z[393]*
  z[604])) + par__Rtibia_M*(z[28]*(z[610]+z[394]*z[611])+z[26]*(z[608]-z[394]*
  z[609])) + par__Rfemur_M*(par__Rfemur_CMx*(z[598]+z[392]*z[599])-par__Rfemur_CMy*(
  z[452]-z[392]*z[597]));
  z[842] = par__Rankle_y*z[19]*z[327] + par__Rankle_y*z[21]*z[328] + z[303]*
  z[34]*z[644] + z[303]*z[35]*z[637] + z[303]*z[44]*z[107] + z[303]*z[47]*
  z[109] + z[304]*z[36]*z[644] + z[304]*z[37]*z[637] + z[304]*z[44]*z[111] + 
  z[304]*z[47]*z[112] + z[305]*z[18]*z[632] + z[305]*z[19]*z[26] + z[305]*
  z[19]*z[631] + z[305]*z[28]*z[104] + z[306]*z[20]*z[632] + z[306]*z[21]*
  z[26] + z[306]*z[21]*z[631] + z[306]*z[28]*z[105] + z[5]*z[101]*z[327] + 
  z[7]*z[101]*z[328] + z[34]*z[327]*z[644] + z[35]*z[327]*z[637] + z[36]*
  z[328]*z[644] + z[37]*z[328]*z[637] + z[107]*z[321]*z[327] + z[109]*z[324]*
  z[327] + z[111]*z[321]*z[328] + z[112]*z[324]*z[328] - z[5]*z[100]*z[328] - 
  z[6]*z[100]*z[327];
  z[843] = par__Rfoot_M*(z[637]*z[480]+z[644]*z[479]+z[47]*(z[643]+z[393]*
  z[645])+z[44]*(z[636]-z[393]*z[638])) + par__Rtibia_M*(z[631]*z[486]+z[632]*
  z[485]+z[28]*(z[648]+z[394]*z[649])+z[26]*(z[646]-z[394]*z[647]));
  z[844] = 2*par__Rfoot_M*(z[44]*z[644]+z[47]*z[637]) + 2*par__Rtibia_M*(
  z[26]*z[632]+z[28]*z[631]);
  z[845] = par__Rfoot_M*(z[46]*z[644]+z[49]*z[637]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[631]-par__Rtibia_CMy*z[632]);
  z[846] = par__Rfoot_M*(par__Rfoot_CMx*z[637]-par__Rfoot_CMy*z[644]);
  z[847] = z[303]*z[34]*z[41] + z[303]*z[35]*z[44] + z[303]*z[35]*z[660] + 
  z[303]*z[47]*z[110] + z[304]*z[36]*z[41] + z[304]*z[37]*z[44] + z[304]*
  z[37]*z[660] + z[304]*z[47]*z[113] + z[34]*z[41]*z[327] + z[35]*z[321]*
  z[327] + z[35]*z[327]*z[660] + z[36]*z[41]*z[328] + z[37]*z[321]*z[328] + 
  z[37]*z[328]*z[660] + z[110]*z[324]*z[327] + z[113]*z[324]*z[328];
  z[848] = par__Rfoot_M*(z[41]*z[479]+z[660]*z[480]+z[47]*(z[666]+z[393]*
  z[667])+z[44]*(z[659]-z[393]*z[661]));
  z[849] = par__Rfoot_M*(z[41]*z[44]+z[47]*z[660]);
  z[850] = par__Rfoot_M*(z[41]*z[46]+z[49]*z[660]+par__Rankle_y*z[32]*z[47]+
  par__Rankle_y*z[33]*z[44]);
  z[851] = par__Rfoot_M*(par__Rfoot_CMx*z[660]-par__Rfoot_CMy*z[41]);
  z[852] = z[364]*z[103] + z[365]*z[19] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + par__Rhip_x*z[1]*z[327] + par__Rhip_x*z[2]*z[328] + par__Rhip_y*
  z[1]*z[328] + z[303]*z[46]*z[106] + z[303]*z[49]*z[108] + z[304]*z[34]*
  z[46] + z[304]*z[35]*z[49] + z[5]*z[13]*z[328] + z[6]*z[13]*z[327] + z[34]*
  z[323]*z[328] + z[35]*z[326]*z[328] + z[106]*z[323]*z[327] + z[108]*z[326]*
  z[327] - z[366]*z[102] - z[367]*z[18] - par__Rhip_y*z[2]*z[327] - z[6]*
  z[12]*z[328] - z[12]*z[99]*z[327];
  z[853] = par__Rfoot_M*(z[49]*(z[581]+z[393]*z[582])+z[46]*(z[578]-z[393]*
  z[579])) + par__Rtibia_M*(par__Rtibia_CMx*(z[585]+z[394]*z[586])-par__Rtibia_CMy*(
  z[583]-z[394]*z[584]));
  z[854] = z[364]*z[103] + z[365]*z[19] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + z[303]*z[46]*z[106] + z[303]*z[49]*z[108] + z[304]*z[34]*
  z[46] + z[304]*z[35]*z[49] + z[5]*z[13]*z[328] + z[6]*z[13]*z[327] + z[34]*
  z[323]*z[328] + z[35]*z[326]*z[328] + z[106]*z[323]*z[327] + z[108]*z[326]*
  z[327] - z[366]*z[102] - z[367]*z[18] - z[6]*z[12]*z[328] - z[12]*z[99]*
  z[327];
  z[855] = par__Rfoot_M*(z[49]*(z[606]+z[393]*z[607])+z[46]*(z[603]-z[393]*
  z[604])) + par__Rtibia_M*(par__Rtibia_CMx*(z[610]+z[394]*z[611])-par__Rtibia_CMy*(
  z[608]-z[394]*z[609]));
  z[856] = z[364]*z[104] + z[365]*z[105] + par__Rankle_y*z[19]*z[327] + 
  par__Rankle_y*z[21]*z[328] + z[303]*z[46]*z[107] + z[303]*z[49]*z[109] + 
  z[304]*z[46]*z[111] + z[304]*z[49]*z[112] + z[5]*z[101]*z[327] + z[7]*
  z[101]*z[328] + z[107]*z[323]*z[327] + z[109]*z[326]*z[327] + z[111]*z[323]*
  z[328] + z[112]*z[326]*z[328] - z[366]*z[19] - z[367]*z[21] - z[5]*z[100]*
  z[328] - z[6]*z[100]*z[327];
  z[857] = par__Rfoot_M*(z[49]*(z[643]+z[393]*z[645])+z[46]*(z[636]-z[393]*
  z[638])) + par__Rtibia_M*(par__Rtibia_CMx*(z[648]+z[394]*z[649])-par__Rtibia_CMy*(
  z[646]-z[394]*z[647]));
  z[858] = par__Rankle_y*z[303];
  z[859] = par__Rankle_y*z[304];
  z[860] = z[303]*z[35]*z[46] + z[303]*z[49]*z[110] + z[304]*z[37]*z[46] + 
  z[304]*z[49]*z[113] + z[858]*z[32]*z[35] + z[858]*z[33]*z[34] + z[859]*
  z[32]*z[37] + z[859]*z[33]*z[36] + z[35]*z[323]*z[327] + z[37]*z[323]*
  z[328] + z[110]*z[326]*z[327] + z[113]*z[326]*z[328] + par__Rankle_y*z[32]*
  z[35]*z[327] + par__Rankle_y*z[32]*z[37]*z[328] + par__Rankle_y*z[33]*z[34]*
  z[327] + par__Rankle_y*z[33]*z[36]*z[328];
  z[861] = par__Rfoot_M*(par__Rankle_y*z[32]*z[480]+par__Rankle_y*z[33]*
  z[479]+z[49]*(z[666]+z[393]*z[667])+z[46]*(z[659]-z[393]*z[661]));
  z[862] = par__Rankle_y*par__Rfoot_M;
  z[863] = z[862]*(z[32]*z[49]+z[33]*z[46]);
  z[864] = z[862]*(par__Rfoot_CMx*z[32]-par__Rfoot_CMy*z[33]);
  z[865] = z[369]*z[108] + z[370]*z[35] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + par__Rhip_x*z[1]*z[327] + par__Rhip_x*z[2]*z[328] + par__Rhip_y*
  z[1]*z[328] + z[5]*z[13]*z[328] + z[6]*z[13]*z[327] - z[371]*z[106] - 
  z[372]*z[34] - par__Rhip_y*z[2]*z[327] - z[6]*z[12]*z[328] - z[12]*z[99]*
  z[327];
  z[866] = par__Rfoot_M*(par__Rfoot_CMx*(z[581]+z[393]*z[582])-par__Rfoot_CMy*(
  z[578]-z[393]*z[579]));
  z[867] = z[369]*z[108] + z[370]*z[35] + par__Rankle_y*z[18]*z[328] + par__Rankle_y*
  z[102]*z[327] + z[5]*z[13]*z[328] + z[6]*z[13]*z[327] - z[371]*z[106] - 
  z[372]*z[34] - z[6]*z[12]*z[328] - z[12]*z[99]*z[327];
  z[868] = par__Rfoot_M*(par__Rfoot_CMx*(z[606]+z[393]*z[607])-par__Rfoot_CMy*(
  z[603]-z[393]*z[604]));
  z[869] = z[369]*z[109] + z[370]*z[112] + par__Rankle_y*z[19]*z[327] + 
  par__Rankle_y*z[21]*z[328] + z[5]*z[101]*z[327] + z[7]*z[101]*z[328] - 
  z[371]*z[107] - z[372]*z[111] - z[5]*z[100]*z[328] - z[6]*z[100]*z[327];
  z[870] = par__Rfoot_M*(par__Rfoot_CMx*(z[643]+z[393]*z[645])-par__Rfoot_CMy*(
  z[636]-z[393]*z[638]));
  z[871] = z[369]*z[110] + z[370]*z[113] - z[371]*z[35] - z[372]*z[37];
  z[872] = par__Rfoot_M*(par__Rfoot_CMx*(z[666]+z[393]*z[667])-par__Rfoot_CMy*(
  z[659]-z[393]*z[661]));
  z[873] = z[374]*z[114] + z[375]*z[53] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + par__Lhip_x*z[1]*z[344] + par__Lhip_x*z[2]*z[345] + par__Lhip_y*
  z[1]*z[345] + z[295]*z[91]*z[121] + z[295]*z[94]*z[123] + z[296]*z[81]*
  z[91] + z[296]*z[82]*z[94] + z[297]*z[73]*z[117] + z[297]*z[75]*z[118] + 
  z[298]*z[65]*z[73] + z[298]*z[66]*z[75] + z[52]*z[60]*z[345] + z[53]*z[60]*
  z[344] + z[81]*z[338]*z[345] + z[82]*z[341]*z[345] + z[121]*z[338]*z[344] + 
  z[123]*z[341]*z[344] - z[376]*z[53] - z[377]*z[52] - par__Lhip_y*z[2]*
  z[344] - z[53]*z[59]*z[345] - z[59]*z[114]*z[344];
  z[874] = par__Lfoot_M*(z[94]*(z[564]+z[390]*z[565])+z[91]*(z[561]-z[390]*
  z[562])) + par__Ltibia_M*(z[75]*(z[568]+z[391]*z[569])+z[73]*(z[566]-z[391]*
  z[567])) + par__Lfemur_M*(par__Lfemur_CMx*(z[556]+z[389]*z[557])-par__Lfemur_CMy*(
  z[553]-z[389]*z[554]));
  z[875] = z[374]*z[114] + z[375]*z[53] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + z[295]*z[91]*z[121] + z[295]*z[94]*z[123] + z[296]*z[81]*
  z[91] + z[296]*z[82]*z[94] + z[297]*z[73]*z[117] + z[297]*z[75]*z[118] + 
  z[298]*z[65]*z[73] + z[298]*z[66]*z[75] + z[52]*z[60]*z[345] + z[53]*z[60]*
  z[344] + z[81]*z[338]*z[345] + z[82]*z[341]*z[345] + z[121]*z[338]*z[344] + 
  z[123]*z[341]*z[344] - z[376]*z[53] - z[377]*z[52] - z[53]*z[59]*z[345] - 
  z[59]*z[114]*z[344];
  z[876] = par__Lfoot_M*(z[94]*(z[690]+z[390]*z[691])+z[91]*(z[687]-z[390]*
  z[688])) + par__Ltibia_M*(z[75]*(z[694]+z[391]*z[695])+z[73]*(z[692]-z[391]*
  z[693])) + par__Lfemur_M*(par__Lfemur_CMx*(z[682]+z[389]*z[683])-par__Lfemur_CMy*(
  z[403]-z[389]*z[681]));
  z[877] = par__Lankle_y*z[66]*z[344] + par__Lankle_y*z[68]*z[345] + z[295]*
  z[81]*z[728] + z[295]*z[82]*z[721] + z[295]*z[91]*z[122] + z[295]*z[94]*
  z[124] + z[296]*z[83]*z[728] + z[296]*z[84]*z[721] + z[296]*z[91]*z[126] + 
  z[296]*z[94]*z[127] + z[297]*z[65]*z[716] + z[297]*z[66]*z[73] + z[297]*
  z[66]*z[715] + z[297]*z[75]*z[119] + z[298]*z[67]*z[716] + z[298]*z[68]*
  z[73] + z[298]*z[68]*z[715] + z[298]*z[75]*z[120] + z[52]*z[116]*z[344] + 
  z[54]*z[116]*z[345] + z[81]*z[344]*z[728] + z[82]*z[344]*z[721] + z[83]*
  z[345]*z[728] + z[84]*z[345]*z[721] + z[122]*z[338]*z[344] + z[124]*z[341]*
  z[344] + z[126]*z[338]*z[345] + z[127]*z[341]*z[345] - z[52]*z[115]*z[345] - 
  z[53]*z[115]*z[344];
  z[878] = par__Lfoot_M*(z[721]*z[431]+z[728]*z[430]+z[94]*(z[727]+z[390]*
  z[729])+z[91]*(z[720]-z[390]*z[722])) + par__Ltibia_M*(z[715]*z[437]+z[716]*
  z[436]+z[75]*(z[732]+z[391]*z[733])+z[73]*(z[730]-z[391]*z[731]));
  z[879] = 2*par__Lfoot_M*(z[91]*z[728]+z[94]*z[721]) + 2*par__Ltibia_M*(
  z[73]*z[716]+z[75]*z[715]);
  z[880] = par__Lfoot_M*(z[93]*z[728]+z[96]*z[721]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[715]-par__Ltibia_CMy*z[716]);
  z[881] = par__Lfoot_M*(par__Lfoot_CMx*z[721]-par__Lfoot_CMy*z[728]);
  z[882] = z[295]*z[81]*z[88] + z[295]*z[82]*z[91] + z[295]*z[82]*z[744] + 
  z[295]*z[94]*z[125] + z[296]*z[83]*z[88] + z[296]*z[84]*z[91] + z[296]*
  z[84]*z[744] + z[296]*z[94]*z[128] + z[81]*z[88]*z[344] + z[82]*z[338]*
  z[344] + z[82]*z[344]*z[744] + z[83]*z[88]*z[345] + z[84]*z[338]*z[345] + 
  z[84]*z[345]*z[744] + z[125]*z[341]*z[344] + z[128]*z[341]*z[345];
  z[883] = par__Lfoot_M*(z[88]*z[430]+z[744]*z[431]+z[94]*(z[750]+z[390]*
  z[751])+z[91]*(z[743]-z[390]*z[745]));
  z[884] = par__Lfoot_M*(z[88]*z[91]+z[94]*z[744]);
  z[885] = par__Lfoot_M*(z[88]*z[93]+z[96]*z[744]+par__Lankle_y*z[79]*z[94]+
  par__Lankle_y*z[80]*z[91]);
  z[886] = par__Lfoot_M*(par__Lfoot_CMx*z[744]-par__Lfoot_CMy*z[88]);
  z[887] = z[379]*z[118] + z[380]*z[66] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + par__Lhip_x*z[1]*z[344] + par__Lhip_x*z[2]*z[345] + par__Lhip_y*
  z[1]*z[345] + z[295]*z[93]*z[121] + z[295]*z[96]*z[123] + z[296]*z[81]*
  z[93] + z[296]*z[82]*z[96] + z[52]*z[60]*z[345] + z[53]*z[60]*z[344] + 
  z[81]*z[340]*z[345] + z[82]*z[343]*z[345] + z[121]*z[340]*z[344] + z[123]*
  z[343]*z[344] - z[381]*z[117] - z[382]*z[65] - par__Lhip_y*z[2]*z[344] - 
  z[53]*z[59]*z[345] - z[59]*z[114]*z[344];
  z[888] = par__Lfoot_M*(z[96]*(z[564]+z[390]*z[565])+z[93]*(z[561]-z[390]*
  z[562])) + par__Ltibia_M*(par__Ltibia_CMx*(z[568]+z[391]*z[569])-par__Ltibia_CMy*(
  z[566]-z[391]*z[567]));
  z[889] = z[379]*z[118] + z[380]*z[66] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + z[295]*z[93]*z[121] + z[295]*z[96]*z[123] + z[296]*z[81]*
  z[93] + z[296]*z[82]*z[96] + z[52]*z[60]*z[345] + z[53]*z[60]*z[344] + 
  z[81]*z[340]*z[345] + z[82]*z[343]*z[345] + z[121]*z[340]*z[344] + z[123]*
  z[343]*z[344] - z[381]*z[117] - z[382]*z[65] - z[53]*z[59]*z[345] - z[59]*
  z[114]*z[344];
  z[890] = par__Lfoot_M*(z[96]*(z[690]+z[390]*z[691])+z[93]*(z[687]-z[390]*
  z[688])) + par__Ltibia_M*(par__Ltibia_CMx*(z[694]+z[391]*z[695])-par__Ltibia_CMy*(
  z[692]-z[391]*z[693]));
  z[891] = z[379]*z[119] + z[380]*z[120] + par__Lankle_y*z[66]*z[344] + 
  par__Lankle_y*z[68]*z[345] + z[295]*z[93]*z[122] + z[295]*z[96]*z[124] + 
  z[296]*z[93]*z[126] + z[296]*z[96]*z[127] + z[52]*z[116]*z[344] + z[54]*
  z[116]*z[345] + z[122]*z[340]*z[344] + z[124]*z[343]*z[344] + z[126]*z[340]*
  z[345] + z[127]*z[343]*z[345] - z[381]*z[66] - z[382]*z[68] - z[52]*z[115]*
  z[345] - z[53]*z[115]*z[344];
  z[892] = par__Lfoot_M*(z[96]*(z[727]+z[390]*z[729])+z[93]*(z[720]-z[390]*
  z[722])) + par__Ltibia_M*(par__Ltibia_CMx*(z[732]+z[391]*z[733])-par__Ltibia_CMy*(
  z[730]-z[391]*z[731]));
  z[893] = par__Lankle_y*z[295];
  z[894] = par__Lankle_y*z[296];
  z[895] = z[295]*z[82]*z[93] + z[295]*z[96]*z[125] + z[296]*z[84]*z[93] + 
  z[296]*z[96]*z[128] + z[893]*z[79]*z[82] + z[893]*z[80]*z[81] + z[894]*
  z[79]*z[84] + z[894]*z[80]*z[83] + z[82]*z[340]*z[344] + z[84]*z[340]*
  z[345] + z[125]*z[343]*z[344] + z[128]*z[343]*z[345] + par__Lankle_y*z[79]*
  z[82]*z[344] + par__Lankle_y*z[79]*z[84]*z[345] + par__Lankle_y*z[80]*z[81]*
  z[344] + par__Lankle_y*z[80]*z[83]*z[345];
  z[896] = par__Lfoot_M*(par__Lankle_y*z[79]*z[431]+par__Lankle_y*z[80]*
  z[430]+z[96]*(z[750]+z[390]*z[751])+z[93]*(z[743]-z[390]*z[745]));
  z[897] = par__Lankle_y*par__Lfoot_M;
  z[898] = z[897]*(z[79]*z[96]+z[80]*z[93]);
  z[899] = z[897]*(par__Lfoot_CMx*z[79]-par__Lfoot_CMy*z[80]);
  z[900] = z[384]*z[123] + z[385]*z[82] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + par__Lhip_x*z[1]*z[344] + par__Lhip_x*z[2]*z[345] + par__Lhip_y*
  z[1]*z[345] + z[52]*z[60]*z[345] + z[53]*z[60]*z[344] - z[386]*z[121] - 
  z[387]*z[81] - par__Lhip_y*z[2]*z[344] - z[53]*z[59]*z[345] - z[59]*z[114]*
  z[344];
  z[901] = par__Lfoot_M*(par__Lfoot_CMx*(z[564]+z[390]*z[565])-par__Lfoot_CMy*(
  z[561]-z[390]*z[562]));
  z[902] = z[384]*z[123] + z[385]*z[82] + par__Lankle_y*z[65]*z[345] + par__Lankle_y*
  z[117]*z[344] + z[52]*z[60]*z[345] + z[53]*z[60]*z[344] - z[386]*z[121] - 
  z[387]*z[81] - z[53]*z[59]*z[345] - z[59]*z[114]*z[344];
  z[903] = par__Lfoot_M*(par__Lfoot_CMx*(z[690]+z[390]*z[691])-par__Lfoot_CMy*(
  z[687]-z[390]*z[688]));
  z[904] = z[384]*z[124] + z[385]*z[127] + par__Lankle_y*z[66]*z[344] + 
  par__Lankle_y*z[68]*z[345] + z[52]*z[116]*z[344] + z[54]*z[116]*z[345] - 
  z[386]*z[122] - z[387]*z[126] - z[52]*z[115]*z[345] - z[53]*z[115]*z[344];
  z[905] = par__Lfoot_M*(par__Lfoot_CMx*(z[727]+z[390]*z[729])-par__Lfoot_CMy*(
  z[720]-z[390]*z[722]));
  z[906] = z[384]*z[125] + z[385]*z[128] - z[386]*z[82] - z[387]*z[84];
  z[907] = par__Lfoot_M*(par__Lfoot_CMx*(z[750]+z[390]*z[751])-par__Lfoot_CMy*(
  z[743]-z[390]*z[745]));
  z[908] = (z[552]-z[589]-z[504]*q3pp-z[590]*q2pp)/par__bodyweight;
  z[909] = (z[328]+z[345])/par__bodyweight;
  z[910] = (z[327]+z[344])/par__bodyweight;
  z[911] = (z[596]-z[612]-z[613]*q1pp-z[614]*q4pp-z[615]*q5pp-z[616]*q6pp-
  z[617]*q3pp-z[618]*q2pp)/par__bodyweight;
  z[912] = (z[623]-z[650]-z[651]*q1pp-z[652]*q4pp-z[653]*q5pp-z[654]*q6pp-
  z[655]*q3pp-z[656]*q2pp)/par__bodyweight;
  z[913] = (z[658]-z[668]-2*z[669]*q1pp-z[670]*q4pp-z[671]*q5pp-z[672]*q6pp-
  z[673]*q3pp-z[674]*q2pp)/par__bodyweight;
  z[914] = (z[680]-z[696]-z[697]*q1pp-z[698]*q7pp-z[699]*q8pp-z[700]*q9pp-
  z[701]*q3pp-z[702]*q2pp)/par__bodyweight;
  z[915] = (z[707]-z[734]-z[735]*q1pp-z[736]*q7pp-z[737]*q8pp-z[738]*q9pp-
  z[739]*q3pp-z[740]*q2pp)/par__bodyweight;
  z[916] = (z[742]-z[752]-2*z[753]*q1pp-z[754]*q7pp-z[755]*q8pp-z[756]*q9pp-
  z[757]*q3pp-z[758]*q2pp)/par__bodyweight;
  z[917] = (z[759]-z[760]-z[590]*q1pp-z[761]*q7pp-z[762]*q3pp-z[763]*q2pp-
  z[764]*q9pp-z[765]*q8pp-z[766]*q4pp-z[767]*q6pp-z[768]*q5pp)/par__bodyweight;
  z[918] = (z[769]-z[770]-z[618]*q1pp-z[766]*q4pp-z[767]*q6pp-z[768]*q5pp-
  z[771]*q3pp-z[772]*q2pp)/par__bodyweight;
  z[919] = (z[773]-z[774]-z[656]*q1pp-z[775]*q3pp-z[776]*q2pp-z[777]*q4pp-
  z[778]*q6pp-z[779]*q5pp)/par__bodyweight;
  z[920] = (z[780]-z[781]-2*z[783]*q2pp-z[674]*q1pp-z[782]*q3pp-z[784]*q4pp-
  z[785]*q6pp-z[786]*q5pp)/par__bodyweight;
  z[921] = (z[787]-z[788]-z[702]*q1pp-z[761]*q7pp-z[764]*q9pp-z[765]*q8pp-
  z[789]*q3pp-z[790]*q2pp)/par__bodyweight;
  z[922] = (z[791]-z[792]-z[740]*q1pp-z[793]*q7pp-z[794]*q3pp-z[795]*q2pp-
  z[796]*q9pp-z[797]*q8pp)/par__bodyweight;
  z[923] = (z[798]-z[799]-2*z[802]*q2pp-z[758]*q1pp-z[800]*q7pp-z[801]*q3pp-
  z[803]*q9pp-z[804]*q8pp)/par__bodyweight;
  z[924] = (z[805]-z[806]-z[504]*q1pp-z[506]*q7pp-z[509]*q9pp-z[510]*q8pp-
  z[511]*q4pp-z[512]*q6pp-z[513]*q5pp-z[762]*q2pp-z[807]*q3pp)/par__bodyweight;
  z[925] = (z[808]-z[809]-z[511]*q4pp-z[512]*q6pp-z[513]*q5pp-z[617]*q1pp-
  z[771]*q2pp-z[810]*q3pp)/par__bodyweight;
  z[926] = (z[811]-z[812]-z[655]*q1pp-z[775]*q2pp-z[813]*q3pp-z[814]*q4pp-
  z[815]*q6pp-z[816]*q5pp)/par__bodyweight;
  z[927] = (z[817]-z[818]-2*z[819]*q3pp-z[673]*q1pp-z[782]*q2pp-z[820]*q4pp-
  z[821]*q6pp-z[822]*q5pp)/par__bodyweight;
  z[928] = (z[823]-z[824]-z[506]*q7pp-z[509]*q9pp-z[510]*q8pp-z[701]*q1pp-
  z[789]*q2pp-z[825]*q3pp)/par__bodyweight;
  z[929] = (z[826]-z[827]-z[739]*q1pp-z[794]*q2pp-z[828]*q7pp-z[829]*q3pp-
  z[830]*q9pp-z[831]*q8pp)/par__bodyweight;
  z[930] = (z[832]-z[833]-2*z[835]*q3pp-z[757]*q1pp-z[801]*q2pp-z[834]*q7pp-
  z[836]*q9pp-z[837]*q8pp)/par__bodyweight;
  z[931] = (z[838]-z[839]-z[511]*q3pp-z[766]*q2pp)/par__bodyweight;
  z[932] = z[328]/par__bodyweight;
  z[933] = z[327]/par__bodyweight;
  z[934] = (z[840]-z[841]-z[511]*q3pp-z[614]*q1pp-z[766]*q2pp)/par__bodyweight;
  z[935] = (z[842]-z[843]-z[652]*q1pp-z[777]*q2pp-z[814]*q3pp-z[844]*q4pp-
  z[845]*q5pp-z[846]*q6pp)/par__bodyweight;
  z[936] = (z[847]-z[848]-2*z[849]*q4pp-z[670]*q1pp-z[784]*q2pp-z[820]*q3pp-
  z[850]*q5pp-z[851]*q6pp)/par__bodyweight;
  z[937] = (z[852]-z[853]-z[513]*q3pp-z[768]*q2pp)/par__bodyweight;
  z[938] = (z[854]-z[855]-z[513]*q3pp-z[615]*q1pp-z[768]*q2pp)/par__bodyweight;
  z[939] = (z[856]-z[857]-z[653]*q1pp-z[779]*q2pp-z[816]*q3pp-z[845]*q4pp)/
  par__bodyweight;
  z[940] = (z[860]-z[861]-2*z[863]*q5pp-z[671]*q1pp-z[786]*q2pp-z[822]*q3pp-
  z[850]*q4pp-z[864]*q6pp)/par__bodyweight;
  z[941] = (z[865]-z[866]-z[512]*q3pp-z[767]*q2pp)/par__bodyweight;
  z[942] = (z[867]-z[868]-z[512]*q3pp-z[616]*q1pp-z[767]*q2pp)/par__bodyweight;
  z[943] = (z[869]-z[870]-z[654]*q1pp-z[778]*q2pp-z[815]*q3pp-z[846]*q4pp)/
  par__bodyweight;
  z[944] = (z[871]-z[872]-z[672]*q1pp-z[785]*q2pp-z[821]*q3pp-z[851]*q4pp-
  z[864]*q5pp)/par__bodyweight;
  z[945] = (z[873]-z[874]-z[506]*q3pp-z[761]*q2pp)/par__bodyweight;
  z[946] = z[345]/par__bodyweight;
  z[947] = z[344]/par__bodyweight;
  z[948] = (z[875]-z[876]-z[506]*q3pp-z[698]*q1pp-z[761]*q2pp)/par__bodyweight;
  z[949] = (z[877]-z[878]-z[736]*q1pp-z[793]*q2pp-z[828]*q3pp-z[879]*q7pp-
  z[880]*q8pp-z[881]*q9pp)/par__bodyweight;
  z[950] = (z[882]-z[883]-2*z[884]*q7pp-z[754]*q1pp-z[800]*q2pp-z[834]*q3pp-
  z[885]*q8pp-z[886]*q9pp)/par__bodyweight;
  z[951] = (z[887]-z[888]-z[510]*q3pp-z[765]*q2pp)/par__bodyweight;
  z[952] = (z[889]-z[890]-z[510]*q3pp-z[699]*q1pp-z[765]*q2pp)/par__bodyweight;
  z[953] = (z[891]-z[892]-z[737]*q1pp-z[797]*q2pp-z[831]*q3pp-z[880]*q7pp)/
  par__bodyweight;
  z[954] = (z[895]-z[896]-2*z[898]*q8pp-z[755]*q1pp-z[804]*q2pp-z[837]*q3pp-
  z[885]*q7pp-z[899]*q9pp)/par__bodyweight;
  z[955] = (z[900]-z[901]-z[509]*q3pp-z[764]*q2pp)/par__bodyweight;
  z[956] = (z[902]-z[903]-z[509]*q3pp-z[700]*q1pp-z[764]*q2pp)/par__bodyweight;
  z[957] = (z[904]-z[905]-z[738]*q1pp-z[796]*q2pp-z[830]*q3pp-z[881]*q7pp)/
  par__bodyweight;
  z[958] = (z[906]-z[907]-z[756]*q1pp-z[803]*q2pp-z[836]*q3pp-z[886]*q7pp-
  z[899]*q8pp)/par__bodyweight;
  z[959] = -par__pelvis_CMx*z[2] - par__pelvis_CMy*z[1];
  z[960] = par__pelvis_CMx*z[1] - par__pelvis_CMy*z[2];
  z[961] = (sx*z[959]+sy*z[960])/pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[962] = par__airdrag*(par__pelvis_CMx*s+z[961]*(par__pelvis_CMx*q1p+z[1]*
  q3p-z[2]*q2p));
  z[963] = par__airdrag*(par__pelvis_CMy*s+z[961]*(par__pelvis_CMy*q1p-z[1]*
  q2p-z[2]*q3p));
  z[964] = -z[550]*z[2]*z[961] - z[551]*z[1]*z[961] - par__pelvis_CMx*z[962] - 
  par__pelvis_CMy*z[963];
  z[965] = z[52]*q3p + z[53]*q2p + z[399];
  z[966] = z[965] - z[396] - z[58]*z[389];
  z[967] = z[53]*q3p + z[114]*q2p - z[402];
  z[968] = z[395] + z[967] + z[57]*z[389];
  z[969] = z[81]*q3p + z[121]*q2p + z[421];
  z[970] = z[969] - z[407] - z[95]*z[390];
  z[971] = z[82]*q3p + z[123]*q2p + z[427];
  z[972] = z[406] + z[971] + z[92]*z[390];
  z[973] = z[65]*q3p + z[117]*q2p + z[420];
  z[974] = z[973] - z[433] - z[76]*z[391];
  z[975] = z[66]*q3p + z[118]*q2p + z[419];
  z[976] = z[432] + z[975] + z[74]*z[391];
  z[977] = z[5]*q3p + z[6]*q2p + z[448];
  z[978] = z[977] - z[445] - z[11]*z[392];
  z[979] = z[6]*q3p + z[99]*q2p - z[451];
  z[980] = z[444] + z[979] + z[10]*z[392];
  z[981] = z[34]*q3p + z[106]*q2p + z[470];
  z[982] = z[981] - z[456] - z[48]*z[393];
  z[983] = z[35]*q3p + z[108]*q2p + z[476];
  z[984] = z[455] + z[983] + z[45]*z[393];
  z[985] = z[18]*q3p + z[102]*q2p + z[469];
  z[986] = z[985] - z[482] - z[29]*z[394];
  z[987] = z[19]*q3p + z[103]*q2p + z[468];
  z[988] = z[481] + z[987] + z[27]*z[394];
  z[989] = z[588] - z[98]*q1p - z[488];
  z[990] = z[1]*q2p + z[2]*q3p;
  z[991] = z[97]*q1p + z[487] - z[990];
  z[992] = z[438] - par__pelvis_CMy*q1p - z[990];
  z[993] = z[588] - par__pelvis_CMx*q1p - z[439];
  z[994] = par__Lfemur_M*(z[57]*z[966]+z[58]*z[968]) + par__Lfoot_M*(z[92]*
  z[970]+z[95]*z[972]) + par__Ltibia_M*(z[74]*z[974]+z[76]*z[976]) + par__Rfemur_M*(
  z[10]*z[978]+z[11]*z[980]) + par__Rfoot_M*(z[45]*z[982]+z[48]*z[984]) + 
  par__Rtibia_M*(z[27]*z[986]+z[29]*z[988]) + par__torso_M*(z[97]*z[989]+
  z[98]*z[991]) + par__pelvis_M*(par__pelvis_CMx*z[992]-par__pelvis_CMy*
  z[993]);
  z[995] = sx*(pow(z[1],2)+pow(z[2],2))/pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[996] = par__airdrag*(s*z[2]-z[995]*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*q2p));
  z[997] = par__airdrag*(s*z[1]-z[995]*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*q3p));
  z[998] = par__pelvis_CMx*z[996] + par__pelvis_CMy*z[997] - z[550]*z[2]*
  z[995] - z[551]*z[1]*z[995];
  z[999] = par__Lfemur_M*(z[58]*(z[401]+z[52]*z[389])+z[57]*(z[397]-z[53]*
  z[389])) + par__Lfoot_M*(z[95]*(z[424]+z[81]*z[390])+z[92]*(z[410]-z[82]*
  z[390])) + par__Ltibia_M*(z[76]*(z[409]+z[65]*z[391])+z[74]*(z[408]-z[66]*
  z[391])) + par__Rfemur_M*(z[11]*(z[450]+z[5]*z[392])+z[10]*(z[446]-z[6]*
  z[392])) + par__Rfoot_M*(z[48]*(z[473]+z[34]*z[393])+z[45]*(z[459]-z[35]*
  z[393])) + par__Rtibia_M*(z[29]*(z[458]+z[18]*z[394])+z[27]*(z[457]-z[19]*
  z[394]));
  z[1000] = sy*(pow(z[1],2)+pow(z[2],2))/pow((1.0E-06+pow(sx,2)+pow(sy,2)),
  0.5);
  z[1001] = par__airdrag*(s*z[1]+z[1000]*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*
  q2p));
  z[1002] = par__airdrag*(s*z[2]-z[1000]*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*
  q3p));
  z[1003] = par__pelvis_CMy*z[1002] - z[550]*z[2]*z[1000] - z[551]*z[1]*
  z[1000] - par__pelvis_CMx*z[1001];
  z[1004] = par__Lfemur_M*(z[58]*(z[397]+z[54]*z[389])+z[57]*(z[398]-z[52]*
  z[389])) + par__Lfoot_M*(z[95]*(z[425]+z[83]*z[390])+z[92]*(z[413]-z[84]*
  z[390])) + par__Ltibia_M*(z[76]*(z[412]+z[67]*z[391])+z[74]*(z[411]-z[68]*
  z[391])) + par__Rfemur_M*(z[11]*(z[446]+z[7]*z[392])+z[10]*(z[447]-z[5]*
  z[392])) + par__Rfoot_M*(z[48]*(z[474]+z[36]*z[393])+z[45]*(z[462]-z[37]*
  z[393])) + par__Rtibia_M*(z[29]*(z[461]+z[20]*z[394])+z[27]*(z[460]-z[21]*
  z[394]));
  z[1005] = z[599] - z[445] - par__Rfemur_CMx*z[392];
  z[1006] = z[6]*q3p + z[99]*q2p - z[8]*q1p;
  z[1007] = z[444] + z[1006] - par__Rfemur_CMy*z[392];
  z[1008] = z[9]*z[16] - z[8]*z[17];
  z[1009] = -z[8]*z[16] - z[9]*z[17];
  z[1010] = z[32]*z[1008] + z[33]*z[1009];
  z[1011] = z[34]*q3p + z[106]*q2p + z[1010]*q1p + z[465];
  z[1012] = z[1011] - z[456] - z[47]*z[393];
  z[1013] = z[32]*z[1009] - z[33]*z[1008];
  z[1014] = z[35]*q3p + z[108]*q2p + z[1013]*q1p + z[475];
  z[1015] = z[455] + z[1014] + z[44]*z[393];
  z[1016] = z[18]*q3p + z[102]*q2p + z[1008]*q1p - z[464];
  z[1017] = z[1016] - z[482] - z[28]*z[394];
  z[1018] = z[19]*q3p + z[103]*q2p + z[1009]*q1p - z[463];
  z[1019] = z[481] + z[1018] + z[26]*z[394];
  z[1020] = par__Rfemur_M*(z[10]*z[1005]+z[11]*z[1007]) + par__Rfoot_M*(z[45]*
  z[1012]+z[48]*z[1015]) + par__Rtibia_M*(z[27]*z[1017]+z[29]*z[1019]);
  z[1021] = z[16]*z[101] - z[12]*z[16] - z[13]*z[17] - z[17]*z[100];
  z[1022] = z[12]*z[17] - z[13]*z[16] - z[16]*z[100] - z[17]*z[101];
  z[1023] = -z[32]*z[1021] - z[33]*z[1022];
  z[1024] = z[107]*q2p + z[111]*q3p + z[621]*q1p + z[1023]*q4p + z[471];
  z[1025] = z[1024] - z[456] - z[49]*z[393];
  z[1026] = z[33]*z[1021] - z[32]*z[1022];
  z[1027] = z[109]*q2p + z[112]*q3p + z[622]*q1p + z[1026]*q4p + z[477];
  z[1028] = z[455] + z[1027] + z[46]*z[393];
  z[1029] = z[19]*q2p + z[21]*q3p + z[620]*q1p - z[1021]*q4p;
  z[1030] = z[1029] - z[482] - par__Rtibia_CMx*z[394];
  z[1031] = z[104]*q2p + z[105]*q3p + z[619]*q1p - z[1022]*q4p;
  z[1032] = z[481] + z[1031] - par__Rtibia_CMy*z[394];
  z[1033] = par__Rfoot_M*(z[45]*z[1025]+z[48]*z[1028]) + par__Rtibia_M*(z[27]*
  z[1030]+z[29]*z[1032]);
  z[1034] = z[667] - z[456] - par__Rfoot_CMx*z[393];
  z[1035] = z[455] + z[661] - par__Rfoot_CMy*z[393];
  z[1036] = par__Rfoot_M*(z[45]*z[1034]+z[48]*z[1035]);
  z[1037] = z[683] - z[396] - par__Lfemur_CMx*z[389];
  z[1038] = z[53]*q3p + z[114]*q2p - z[55]*q1p;
  z[1039] = z[395] + z[1038] - par__Lfemur_CMy*z[389];
  z[1040] = z[56]*z[63] - z[55]*z[64];
  z[1041] = -z[55]*z[63] - z[56]*z[64];
  z[1042] = z[79]*z[1040] + z[80]*z[1041];
  z[1043] = z[81]*q3p + z[121]*q2p + z[1042]*q1p + z[416];
  z[1044] = z[1043] - z[407] - z[94]*z[390];
  z[1045] = z[79]*z[1041] - z[80]*z[1040];
  z[1046] = z[82]*q3p + z[123]*q2p + z[1045]*q1p + z[426];
  z[1047] = z[406] + z[1046] + z[91]*z[390];
  z[1048] = z[65]*q3p + z[117]*q2p + z[1040]*q1p - z[415];
  z[1049] = z[1048] - z[433] - z[75]*z[391];
  z[1050] = z[66]*q3p + z[118]*q2p + z[1041]*q1p - z[414];
  z[1051] = z[432] + z[1050] + z[73]*z[391];
  z[1052] = par__Lfemur_M*(z[57]*z[1037]+z[58]*z[1039]) + par__Lfoot_M*(z[92]*
  z[1044]+z[95]*z[1047]) + par__Ltibia_M*(z[74]*z[1049]+z[76]*z[1051]);
  z[1053] = z[63]*z[116] - z[59]*z[63] - z[60]*z[64] - z[64]*z[115];
  z[1054] = z[59]*z[64] - z[60]*z[63] - z[63]*z[115] - z[64]*z[116];
  z[1055] = -z[79]*z[1053] - z[80]*z[1054];
  z[1056] = z[122]*q2p + z[126]*q3p + z[705]*q1p + z[1055]*q7p + z[422];
  z[1057] = z[1056] - z[407] - z[96]*z[390];
  z[1058] = z[80]*z[1053] - z[79]*z[1054];
  z[1059] = z[124]*q2p + z[127]*q3p + z[706]*q1p + z[1058]*q7p + z[428];
  z[1060] = z[406] + z[1059] + z[93]*z[390];
  z[1061] = z[66]*q2p + z[68]*q3p + z[704]*q1p - z[1053]*q7p;
  z[1062] = z[1061] - z[433] - par__Ltibia_CMx*z[391];
  z[1063] = z[119]*q2p + z[120]*q3p + z[703]*q1p - z[1054]*q7p;
  z[1064] = z[432] + z[1063] - par__Ltibia_CMy*z[391];
  z[1065] = par__Lfoot_M*(z[92]*z[1057]+z[95]*z[1060]) + par__Ltibia_M*(z[74]*
  z[1062]+z[76]*z[1064]);
  z[1066] = z[751] - z[407] - par__Lfoot_CMx*z[390];
  z[1067] = z[406] + z[745] - par__Lfoot_CMy*z[390];
  z[1068] = par__Lfoot_M*(z[92]*z[1066]+z[95]*z[1067]);
  z[1069] = z[309]*z[961]*pow(z[1],2) + z[309]*z[961]*pow(z[2],2) + z[1]*
  z[963] + z[2]*z[962];
  z[1070] = par__Lfemur_M*(z[52]*z[966]+z[53]*z[968]) + par__Lfoot_M*(z[81]*
  z[970]+z[82]*z[972]) + par__Ltibia_M*(z[65]*z[974]+z[66]*z[976]) + par__Rfemur_M*(
  z[5]*z[978]+z[6]*z[980]) + par__Rfoot_M*(z[34]*z[982]+z[35]*z[984]) + 
  par__Rtibia_M*(z[18]*z[986]+z[19]*z[988]) + par__pelvis_M*(z[1]*z[993]-z[2]*
  z[992]) + par__torso_M*(z[1]*z[989]-z[2]*z[991]);
  z[1071] = z[309]*z[995]*pow(z[1],2) + z[309]*z[995]*pow(z[2],2) - z[1]*
  z[997] - z[2]*z[996];
  z[1072] = par__Lfemur_M*(z[52]*z[397]+z[53]*z[401]) + par__Lfoot_M*(z[81]*
  z[410]+z[82]*z[424]) + par__Ltibia_M*(z[65]*z[408]+z[66]*z[409]) + par__Rfemur_M*(
  z[5]*z[446]+z[6]*z[450]) + par__Rfoot_M*(z[34]*z[459]+z[35]*z[473]) + 
  par__Rtibia_M*(z[18]*z[457]+z[19]*z[458]);
  z[1073] = z[309]*z[1000]*pow(z[1],2) + z[309]*z[1000]*pow(z[2],2) + z[2]*
  z[1001] - z[1]*z[1002];
  z[1074] = par__Lfemur_M*(z[53]*(z[397]+z[54]*z[389])+z[52]*(z[398]-z[52]*
  z[389])) + par__Lfoot_M*(z[82]*(z[425]+z[83]*z[390])+z[81]*(z[413]-z[84]*
  z[390])) + par__Ltibia_M*(z[66]*(z[412]+z[67]*z[391])+z[65]*(z[411]-z[68]*
  z[391])) + par__Rfemur_M*(z[6]*(z[446]+z[7]*z[392])+z[5]*(z[447]-z[5]*
  z[392])) + par__Rfoot_M*(z[35]*(z[474]+z[36]*z[393])+z[34]*(z[462]-z[37]*
  z[393])) + par__Rtibia_M*(z[19]*(z[461]+z[20]*z[394])+z[18]*(z[460]-z[21]*
  z[394]));
  z[1075] = par__Rfemur_M*(z[5]*z[1005]+z[6]*z[1007]) + par__Rfoot_M*(z[34]*
  z[1012]+z[35]*z[1015]) + par__Rtibia_M*(z[18]*z[1017]+z[19]*z[1019]);
  z[1076] = par__Rfoot_M*(z[34]*z[1025]+z[35]*z[1028]) + par__Rtibia_M*(z[18]*
  z[1030]+z[19]*z[1032]);
  z[1077] = par__Rfoot_M*(z[34]*z[1034]+z[35]*z[1035]);
  z[1078] = par__Lfemur_M*(z[52]*z[1037]+z[53]*z[1039]) + par__Lfoot_M*(z[81]*
  z[1044]+z[82]*z[1047]) + par__Ltibia_M*(z[65]*z[1049]+z[66]*z[1051]);
  z[1079] = par__Lfoot_M*(z[81]*z[1057]+z[82]*z[1060]) + par__Ltibia_M*(z[65]*
  z[1062]+z[66]*z[1064]);
  z[1080] = par__Lfoot_M*(z[81]*z[1066]+z[82]*z[1067]);
  z[1081] = z[2]*z[963] - z[1]*z[962];
  z[1082] = par__Lfemur_M*(z[52]*z[968]+z[54]*z[966]) + par__Lfoot_M*(z[83]*
  z[970]+z[84]*z[972]) + par__Ltibia_M*(z[67]*z[974]+z[68]*z[976]) + par__pelvis_M*(
  z[1]*z[992]+z[2]*z[993]) + par__Rfemur_M*(z[5]*z[980]+z[7]*z[978]) + par__Rfoot_M*(
  z[36]*z[982]+z[37]*z[984]) + par__Rtibia_M*(z[20]*z[986]+z[21]*z[988]) + 
  par__torso_M*(z[1]*z[991]+z[2]*z[989]);
  z[1083] = z[1]*z[996] - z[2]*z[997];
  z[1084] = par__Lfemur_M*(z[52]*(z[401]+z[52]*z[389])+z[54]*(z[397]-z[53]*
  z[389])) + par__Lfoot_M*(z[84]*(z[424]+z[81]*z[390])+z[83]*(z[410]-z[82]*
  z[390])) + par__Ltibia_M*(z[68]*(z[409]+z[65]*z[391])+z[67]*(z[408]-z[66]*
  z[391])) + par__Rfemur_M*(z[5]*(z[450]+z[5]*z[392])+z[7]*(z[446]-z[6]*
  z[392])) + par__Rfoot_M*(z[37]*(z[473]+z[34]*z[393])+z[36]*(z[459]-z[35]*
  z[393])) + par__Rtibia_M*(z[21]*(z[458]+z[18]*z[394])+z[20]*(z[457]-z[19]*
  z[394]));
  z[1085] = -z[1]*z[1001] - z[2]*z[1002];
  z[1086] = par__Lfemur_M*(z[52]*z[397]+z[54]*z[398]) + par__Lfoot_M*(z[83]*
  z[413]+z[84]*z[425]) + par__Ltibia_M*(z[67]*z[411]+z[68]*z[412]) + par__Rfemur_M*(
  z[5]*z[446]+z[7]*z[447]) + par__Rfoot_M*(z[36]*z[462]+z[37]*z[474]) + 
  par__Rtibia_M*(z[20]*z[460]+z[21]*z[461]);
  z[1087] = par__Rfemur_M*(z[5]*z[1007]+z[7]*z[1005]) + par__Rfoot_M*(z[36]*
  z[1012]+z[37]*z[1015]) + par__Rtibia_M*(z[20]*z[1017]+z[21]*z[1019]);
  z[1088] = par__Rfoot_M*(z[36]*z[1025]+z[37]*z[1028]) + par__Rtibia_M*(z[20]*
  z[1030]+z[21]*z[1032]);
  z[1089] = par__Rfoot_M*(z[36]*z[1034]+z[37]*z[1035]);
  z[1090] = par__Lfemur_M*(z[52]*z[1039]+z[54]*z[1037]) + par__Lfoot_M*(z[83]*
  z[1044]+z[84]*z[1047]) + par__Ltibia_M*(z[67]*z[1049]+z[68]*z[1051]);
  z[1091] = par__Lfoot_M*(z[83]*z[1057]+z[84]*z[1060]) + par__Ltibia_M*(z[67]*
  z[1062]+z[68]*z[1064]);
  z[1092] = par__Lfoot_M*(z[83]*z[1066]+z[84]*z[1067]);
  z[1093] = par__Rfoot_M*(z[44]*z[982]+z[47]*z[984]) + par__Rtibia_M*(z[26]*
  z[986]+z[28]*z[988]) + par__Rfemur_M*(par__Rfemur_CMx*z[980]-par__Rfemur_CMy*
  z[978]);
  z[1094] = par__Rfoot_M*(z[47]*(z[473]+z[34]*z[393])+z[44]*(z[459]-z[35]*
  z[393])) + par__Rtibia_M*(z[28]*(z[458]+z[18]*z[394])+z[26]*(z[457]-z[19]*
  z[394])) + par__Rfemur_M*(par__Rfemur_CMx*(z[450]+z[5]*z[392])-par__Rfemur_CMy*(
  z[446]-z[6]*z[392]));
  z[1095] = par__Rfoot_M*(z[47]*(z[474]+z[36]*z[393])+z[44]*(z[462]-z[37]*
  z[393])) + par__Rtibia_M*(z[28]*(z[461]+z[20]*z[394])+z[26]*(z[460]-z[21]*
  z[394])) + par__Rfemur_M*(par__Rfemur_CMx*(z[446]+z[7]*z[392])-par__Rfemur_CMy*(
  z[447]-z[5]*z[392]));
  z[1096] = par__Rfoot_M*(z[44]*z[1012]+z[47]*z[1015]) + par__Rtibia_M*(z[26]*
  z[1017]+z[28]*z[1019]) + par__Rfemur_M*(par__Rfemur_CMx*z[1007]-par__Rfemur_CMy*
  z[1005]);
  z[1097] = par__Rfoot_M*(z[44]*z[1025]+z[47]*z[1028]) + par__Rtibia_M*(z[26]*
  z[1030]+z[28]*z[1032]);
  z[1098] = par__Rfoot_M*(z[44]*z[1034]+z[47]*z[1035]);
  z[1099] = par__Rfoot_M*(z[46]*z[982]+z[49]*z[984]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[988]-par__Rtibia_CMy*z[986]);
  z[1100] = par__Rfoot_M*(z[49]*(z[473]+z[34]*z[393])+z[46]*(z[459]-z[35]*
  z[393])) + par__Rtibia_M*(par__Rtibia_CMx*(z[458]+z[18]*z[394])-par__Rtibia_CMy*(
  z[457]-z[19]*z[394]));
  z[1101] = par__Rfoot_M*(z[49]*(z[474]+z[36]*z[393])+z[46]*(z[462]-z[37]*
  z[393])) + par__Rtibia_M*(par__Rtibia_CMx*(z[461]+z[20]*z[394])-par__Rtibia_CMy*(
  z[460]-z[21]*z[394]));
  z[1102] = par__Rfoot_M*(z[46]*z[1012]+z[49]*z[1015]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[1019]-par__Rtibia_CMy*z[1017]);
  z[1103] = par__Rfoot_M*(z[46]*z[1025]+z[49]*z[1028]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[1032]-par__Rtibia_CMy*z[1030]);
  z[1104] = par__Rfoot_M*(z[46]*z[1034]+z[49]*z[1035]);
  z[1105] = par__Rfoot_M*(par__Rfoot_CMx*z[984]-par__Rfoot_CMy*z[982]);
  z[1106] = par__Rfoot_M*(par__Rfoot_CMx*(z[473]+z[34]*z[393])-par__Rfoot_CMy*(
  z[459]-z[35]*z[393]));
  z[1107] = par__Rfoot_M*(par__Rfoot_CMx*(z[474]+z[36]*z[393])-par__Rfoot_CMy*(
  z[462]-z[37]*z[393]));
  z[1108] = par__Rfoot_M*(par__Rfoot_CMx*z[1015]-par__Rfoot_CMy*z[1012]);
  z[1109] = par__Rfoot_M*(par__Rfoot_CMx*z[1028]-par__Rfoot_CMy*z[1025]);
  z[1110] = par__Rfoot_M*(par__Rfoot_CMx*z[1035]-par__Rfoot_CMy*z[1034]);
  z[1111] = par__Lfoot_M*(z[91]*z[970]+z[94]*z[972]) + par__Ltibia_M*(z[73]*
  z[974]+z[75]*z[976]) + par__Lfemur_M*(par__Lfemur_CMx*z[968]-par__Lfemur_CMy*
  z[966]);
  z[1112] = par__Lfoot_M*(z[94]*(z[424]+z[81]*z[390])+z[91]*(z[410]-z[82]*
  z[390])) + par__Ltibia_M*(z[75]*(z[409]+z[65]*z[391])+z[73]*(z[408]-z[66]*
  z[391])) + par__Lfemur_M*(par__Lfemur_CMx*(z[401]+z[52]*z[389])-par__Lfemur_CMy*(
  z[397]-z[53]*z[389]));
  z[1113] = par__Lfoot_M*(z[94]*(z[425]+z[83]*z[390])+z[91]*(z[413]-z[84]*
  z[390])) + par__Ltibia_M*(z[75]*(z[412]+z[67]*z[391])+z[73]*(z[411]-z[68]*
  z[391])) + par__Lfemur_M*(par__Lfemur_CMx*(z[397]+z[54]*z[389])-par__Lfemur_CMy*(
  z[398]-z[52]*z[389]));
  z[1114] = par__Lfoot_M*(z[91]*z[1044]+z[94]*z[1047]) + par__Ltibia_M*(z[73]*
  z[1049]+z[75]*z[1051]) + par__Lfemur_M*(par__Lfemur_CMx*z[1039]-par__Lfemur_CMy*
  z[1037]);
  z[1115] = par__Lfoot_M*(z[91]*z[1057]+z[94]*z[1060]) + par__Ltibia_M*(z[73]*
  z[1062]+z[75]*z[1064]);
  z[1116] = par__Lfoot_M*(z[91]*z[1066]+z[94]*z[1067]);
  z[1117] = par__Lfoot_M*(z[93]*z[970]+z[96]*z[972]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[976]-par__Ltibia_CMy*z[974]);
  z[1118] = par__Lfoot_M*(z[96]*(z[424]+z[81]*z[390])+z[93]*(z[410]-z[82]*
  z[390])) + par__Ltibia_M*(par__Ltibia_CMx*(z[409]+z[65]*z[391])-par__Ltibia_CMy*(
  z[408]-z[66]*z[391]));
  z[1119] = par__Lfoot_M*(z[96]*(z[425]+z[83]*z[390])+z[93]*(z[413]-z[84]*
  z[390])) + par__Ltibia_M*(par__Ltibia_CMx*(z[412]+z[67]*z[391])-par__Ltibia_CMy*(
  z[411]-z[68]*z[391]));
  z[1120] = par__Lfoot_M*(z[93]*z[1044]+z[96]*z[1047]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[1051]-par__Ltibia_CMy*z[1049]);
  z[1121] = par__Lfoot_M*(z[93]*z[1057]+z[96]*z[1060]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[1064]-par__Ltibia_CMy*z[1062]);
  z[1122] = par__Lfoot_M*(z[93]*z[1066]+z[96]*z[1067]);
  z[1123] = par__Lfoot_M*(par__Lfoot_CMx*z[972]-par__Lfoot_CMy*z[970]);
  z[1124] = par__Lfoot_M*(par__Lfoot_CMx*(z[424]+z[81]*z[390])-par__Lfoot_CMy*(
  z[410]-z[82]*z[390]));
  z[1125] = par__Lfoot_M*(par__Lfoot_CMx*(z[425]+z[83]*z[390])-par__Lfoot_CMy*(
  z[413]-z[84]*z[390]));
  z[1126] = par__Lfoot_M*(par__Lfoot_CMx*z[1047]-par__Lfoot_CMy*z[1044]);
  z[1127] = par__Lfoot_M*(par__Lfoot_CMx*z[1060]-par__Lfoot_CMy*z[1057]);
  z[1128] = par__Lfoot_M*(par__Lfoot_CMx*z[1067]-par__Lfoot_CMy*z[1066]);
  z[1129] = (z[964]-z[994])/par__bodyweight;
  z[1130] = (z[998]-z[999])/par__bodyweight;
  z[1131] = (z[1003]-z[1004])/par__bodyweight;
  z[1132] = z[1020]/par__bodyweight;
  z[1133] = z[1033]/par__bodyweight;
  z[1134] = z[1036]/par__bodyweight;
  z[1135] = z[1052]/par__bodyweight;
  z[1136] = z[1065]/par__bodyweight;
  z[1137] = z[1068]/par__bodyweight;
  z[1138] = (z[1069]-z[1070])/par__bodyweight;
  z[1139] = (z[1071]-z[1072])/par__bodyweight;
  z[1140] = (z[1073]-z[1074])/par__bodyweight;
  z[1141] = z[1075]/par__bodyweight;
  z[1142] = z[1076]/par__bodyweight;
  z[1143] = z[1077]/par__bodyweight;
  z[1144] = z[1078]/par__bodyweight;
  z[1145] = z[1079]/par__bodyweight;
  z[1146] = z[1080]/par__bodyweight;
  z[1147] = (z[1081]-z[1082])/par__bodyweight;
  z[1148] = (z[1083]-z[1084])/par__bodyweight;
  z[1149] = (z[1085]-z[1086])/par__bodyweight;
  z[1150] = z[1087]/par__bodyweight;
  z[1151] = z[1088]/par__bodyweight;
  z[1152] = z[1089]/par__bodyweight;
  z[1153] = z[1090]/par__bodyweight;
  z[1154] = z[1091]/par__bodyweight;
  z[1155] = z[1092]/par__bodyweight;
  z[1156] = z[1093]/par__bodyweight;
  z[1157] = z[1094]/par__bodyweight;
  z[1158] = z[1095]/par__bodyweight;
  z[1159] = z[1096]/par__bodyweight;
  z[1160] = z[1097]/par__bodyweight;
  z[1161] = z[1098]/par__bodyweight;
  z[1162] = z[1099]/par__bodyweight;
  z[1163] = z[1100]/par__bodyweight;
  z[1164] = z[1101]/par__bodyweight;
  z[1165] = z[1102]/par__bodyweight;
  z[1166] = z[1103]/par__bodyweight;
  z[1167] = z[1104]/par__bodyweight;
  z[1168] = z[1105]/par__bodyweight;
  z[1169] = z[1106]/par__bodyweight;
  z[1170] = z[1107]/par__bodyweight;
  z[1171] = z[1108]/par__bodyweight;
  z[1172] = z[1109]/par__bodyweight;
  z[1173] = z[1110]/par__bodyweight;
  z[1174] = z[1111]/par__bodyweight;
  z[1175] = z[1112]/par__bodyweight;
  z[1176] = z[1113]/par__bodyweight;
  z[1177] = z[1114]/par__bodyweight;
  z[1178] = z[1115]/par__bodyweight;
  z[1179] = z[1116]/par__bodyweight;
  z[1180] = z[1117]/par__bodyweight;
  z[1181] = z[1118]/par__bodyweight;
  z[1182] = z[1119]/par__bodyweight;
  z[1183] = z[1120]/par__bodyweight;
  z[1184] = z[1121]/par__bodyweight;
  z[1185] = z[1122]/par__bodyweight;
  z[1186] = z[1123]/par__bodyweight;
  z[1187] = z[1124]/par__bodyweight;
  z[1188] = z[1125]/par__bodyweight;
  z[1189] = z[1126]/par__bodyweight;
  z[1190] = z[1127]/par__bodyweight;
  z[1191] = z[1128]/par__bodyweight;
  z[1192] = z[492]/par__bodyweight;
  z[1193] = z[504]/par__bodyweight;
  z[1194] = z[503]/par__bodyweight;
  z[1195] = z[499]/par__bodyweight;
  z[1196] = z[501]/par__bodyweight;
  z[1197] = z[502]/par__bodyweight;
  z[1198] = z[494]/par__bodyweight;
  z[1199] = z[496]/par__bodyweight;
  z[1200] = z[497]/par__bodyweight;
  z[1201] = z[508]/par__bodyweight;
  z[1202] = z[507]/par__bodyweight;
  z[1203] = z[511]/par__bodyweight;
  z[1204] = z[513]/par__bodyweight;
  z[1205] = z[512]/par__bodyweight;
  z[1206] = z[506]/par__bodyweight;
  z[1207] = z[510]/par__bodyweight;
  z[1208] = z[509]/par__bodyweight;
  z[1209] = z[516]/par__bodyweight;
  z[1210] = z[519]/par__bodyweight;
  z[1211] = z[521]/par__bodyweight;
  z[1212] = z[520]/par__bodyweight;
  z[1213] = z[515]/par__bodyweight;
  z[1214] = z[518]/par__bodyweight;
  z[1215] = z[517]/par__bodyweight;
  z[1216] = z[524]/par__bodyweight;
  z[1217] = z[525]/par__bodyweight;
  z[1218] = z[526]/par__bodyweight;
  z[1219] = z[529]/par__bodyweight;
  z[1220] = z[530]/par__bodyweight;
  z[1221] = z[532]/par__bodyweight;
  z[1222] = z[535]/par__bodyweight;
  z[1223] = z[536]/par__bodyweight;
  z[1224] = z[537]/par__bodyweight;
  z[1225] = z[540]/par__bodyweight;
  z[1226] = z[541]/par__bodyweight;
  z[1227] = z[543]/par__bodyweight;
  z[1228] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[322]-z[35]*z[325]);
  z[1229] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[322]-z[37]*z[325]);
  z[1230] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[339]-z[82]*z[342]);
  z[1231] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[339]-z[84]*z[342]);
  z[1232] = par__bodyweight*(z[34]*z[36]+z[35]*z[37]);
  z[1233] = par__bodyweight*(z[81]*z[83]+z[82]*z[84]);
  z[1234] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[321]-z[35]*z[324]);
  z[1235] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[321]-z[37]*z[324]);
  z[1236] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[323]-z[35]*z[326]);
  z[1237] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[323]-z[37]*z[326]);
  z[1238] = par__bodyweight*(q3+par__Rankle_y*z[18]+par__Rhip_x*z[2]+par__Rhip_y*
  z[1]+z[5]*z[13]-z[6]*z[12]);
  z[1239] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]);
  z[1240] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[338]-z[82]*z[341]);
  z[1241] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[338]-z[84]*z[341]);
  z[1242] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[340]-z[82]*z[343]);
  z[1243] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[340]-z[84]*z[343]);
  z[1244] = par__bodyweight*(q3+par__Lankle_y*z[65]+par__Lhip_x*z[2]+par__Lhip_y*
  z[1]+z[52]*z[60]-z[53]*z[59]);
  z[1245] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]);
  z[1246] = z[1228]/par__bodyweight;
  z[1247] = z[1229]/par__bodyweight;
  z[1248] = z[1230]/par__bodyweight;
  z[1249] = z[1231]/par__bodyweight;
  z[1250] = pow(z[34],2) + pow(z[35],2);
  z[1251] = z[1232]/par__bodyweight;
  z[1252] = pow(z[81],2) + pow(z[82],2);
  z[1253] = z[1233]/par__bodyweight;
  z[1254] = pow(z[36],2) + pow(z[37],2);
  z[1255] = pow(z[83],2) + pow(z[84],2);
  z[1256] = z[1234]/par__bodyweight;
  z[1257] = z[1235]/par__bodyweight;
  z[1258] = z[1236]/par__bodyweight;
  z[1259] = z[1237]/par__bodyweight;
  z[1260] = z[1238]/par__bodyweight;
  z[1261] = z[1239]/par__bodyweight;
  z[1262] = z[1240]/par__bodyweight;
  z[1263] = z[1241]/par__bodyweight;
  z[1264] = z[1242]/par__bodyweight;
  z[1265] = z[1243]/par__bodyweight;
  z[1266] = z[1244]/par__bodyweight;
  z[1267] = z[1245]/par__bodyweight;


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
  df_dq[0][0] = z[908];
  df_dq[0][1] = -z[909];
  df_dq[0][2] = z[910];
  df_dq[0][3] = z[911];
  df_dq[0][4] = z[912];
  df_dq[0][5] = z[913];
  df_dq[0][6] = z[914];
  df_dq[0][7] = z[915];
  df_dq[0][8] = z[916];
  df_dq[1][0] = z[917];
  df_dq[1][1] = 0;
  df_dq[1][2] = 0;
  df_dq[1][3] = z[918];
  df_dq[1][4] = z[919];
  df_dq[1][5] = z[920];
  df_dq[1][6] = z[921];
  df_dq[1][7] = z[922];
  df_dq[1][8] = z[923];
  df_dq[2][0] = z[924];
  df_dq[2][1] = 0;
  df_dq[2][2] = 0;
  df_dq[2][3] = z[925];
  df_dq[2][4] = z[926];
  df_dq[2][5] = z[927];
  df_dq[2][6] = z[928];
  df_dq[2][7] = z[929];
  df_dq[2][8] = z[930];
  df_dq[3][0] = z[931];
  df_dq[3][1] = -z[932];
  df_dq[3][2] = z[933];
  df_dq[3][3] = z[934];
  df_dq[3][4] = z[935];
  df_dq[3][5] = z[936];
  df_dq[3][6] = 0;
  df_dq[3][7] = 0;
  df_dq[3][8] = 0;
  df_dq[4][0] = z[937];
  df_dq[4][1] = -z[932];
  df_dq[4][2] = z[933];
  df_dq[4][3] = z[938];
  df_dq[4][4] = z[939];
  df_dq[4][5] = z[940];
  df_dq[4][6] = 0;
  df_dq[4][7] = 0;
  df_dq[4][8] = 0;
  df_dq[5][0] = z[941];
  df_dq[5][1] = -z[932];
  df_dq[5][2] = z[933];
  df_dq[5][3] = z[942];
  df_dq[5][4] = z[943];
  df_dq[5][5] = z[944];
  df_dq[5][6] = 0;
  df_dq[5][7] = 0;
  df_dq[5][8] = 0;
  df_dq[6][0] = z[945];
  df_dq[6][1] = -z[946];
  df_dq[6][2] = z[947];
  df_dq[6][3] = 0;
  df_dq[6][4] = 0;
  df_dq[6][5] = 0;
  df_dq[6][6] = z[948];
  df_dq[6][7] = z[949];
  df_dq[6][8] = z[950];
  df_dq[7][0] = z[951];
  df_dq[7][1] = -z[946];
  df_dq[7][2] = z[947];
  df_dq[7][3] = 0;
  df_dq[7][4] = 0;
  df_dq[7][5] = 0;
  df_dq[7][6] = z[952];
  df_dq[7][7] = z[953];
  df_dq[7][8] = z[954];
  df_dq[8][0] = z[955];
  df_dq[8][1] = -z[946];
  df_dq[8][2] = z[947];
  df_dq[8][3] = 0;
  df_dq[8][4] = 0;
  df_dq[8][5] = 0;
  df_dq[8][6] = z[956];
  df_dq[8][7] = z[957];
  df_dq[8][8] = z[958];
  df_dqd[0][0] = z[1129];
  df_dqd[0][1] = z[1130];
  df_dqd[0][2] = z[1131];
  df_dqd[0][3] = -z[1132];
  df_dqd[0][4] = -z[1133];
  df_dqd[0][5] = -z[1134];
  df_dqd[0][6] = -z[1135];
  df_dqd[0][7] = -z[1136];
  df_dqd[0][8] = -z[1137];
  df_dqd[1][0] = z[1138];
  df_dqd[1][1] = z[1139];
  df_dqd[1][2] = z[1140];
  df_dqd[1][3] = -z[1141];
  df_dqd[1][4] = -z[1142];
  df_dqd[1][5] = -z[1143];
  df_dqd[1][6] = -z[1144];
  df_dqd[1][7] = -z[1145];
  df_dqd[1][8] = -z[1146];
  df_dqd[2][0] = z[1147];
  df_dqd[2][1] = z[1148];
  df_dqd[2][2] = z[1149];
  df_dqd[2][3] = -z[1150];
  df_dqd[2][4] = -z[1151];
  df_dqd[2][5] = -z[1152];
  df_dqd[2][6] = -z[1153];
  df_dqd[2][7] = -z[1154];
  df_dqd[2][8] = -z[1155];
  df_dqd[3][0] = -z[1156];
  df_dqd[3][1] = -z[1157];
  df_dqd[3][2] = -z[1158];
  df_dqd[3][3] = -z[1159];
  df_dqd[3][4] = -z[1160];
  df_dqd[3][5] = -z[1161];
  df_dqd[3][6] = 0;
  df_dqd[3][7] = 0;
  df_dqd[3][8] = 0;
  df_dqd[4][0] = -z[1162];
  df_dqd[4][1] = -z[1163];
  df_dqd[4][2] = -z[1164];
  df_dqd[4][3] = -z[1165];
  df_dqd[4][4] = -z[1166];
  df_dqd[4][5] = -z[1167];
  df_dqd[4][6] = 0;
  df_dqd[4][7] = 0;
  df_dqd[4][8] = 0;
  df_dqd[5][0] = -z[1168];
  df_dqd[5][1] = -z[1169];
  df_dqd[5][2] = -z[1170];
  df_dqd[5][3] = -z[1171];
  df_dqd[5][4] = -z[1172];
  df_dqd[5][5] = -z[1173];
  df_dqd[5][6] = 0;
  df_dqd[5][7] = 0;
  df_dqd[5][8] = 0;
  df_dqd[6][0] = -z[1174];
  df_dqd[6][1] = -z[1175];
  df_dqd[6][2] = -z[1176];
  df_dqd[6][3] = 0;
  df_dqd[6][4] = 0;
  df_dqd[6][5] = 0;
  df_dqd[6][6] = -z[1177];
  df_dqd[6][7] = -z[1178];
  df_dqd[6][8] = -z[1179];
  df_dqd[7][0] = -z[1180];
  df_dqd[7][1] = -z[1181];
  df_dqd[7][2] = -z[1182];
  df_dqd[7][3] = 0;
  df_dqd[7][4] = 0;
  df_dqd[7][5] = 0;
  df_dqd[7][6] = -z[1183];
  df_dqd[7][7] = -z[1184];
  df_dqd[7][8] = -z[1185];
  df_dqd[8][0] = -z[1186];
  df_dqd[8][1] = -z[1187];
  df_dqd[8][2] = -z[1188];
  df_dqd[8][3] = 0;
  df_dqd[8][4] = 0;
  df_dqd[8][5] = 0;
  df_dqd[8][6] = -z[1189];
  df_dqd[8][7] = -z[1190];
  df_dqd[8][8] = -z[1191];
  df_dqdd[0][0] = -z[1192];
  df_dqdd[0][1] = -z[1193];
  df_dqdd[0][2] = -z[1194];
  df_dqdd[0][3] = -z[1195];
  df_dqdd[0][4] = -z[1196];
  df_dqdd[0][5] = -z[1197];
  df_dqdd[0][6] = -z[1198];
  df_dqdd[0][7] = -z[1199];
  df_dqdd[0][8] = -z[1200];
  df_dqdd[1][0] = -z[1193];
  df_dqdd[1][1] = -z[1201];
  df_dqdd[1][2] = -z[1202];
  df_dqdd[1][3] = -z[1203];
  df_dqdd[1][4] = -z[1204];
  df_dqdd[1][5] = -z[1205];
  df_dqdd[1][6] = -z[1206];
  df_dqdd[1][7] = -z[1207];
  df_dqdd[1][8] = -z[1208];
  df_dqdd[2][0] = -z[1194];
  df_dqdd[2][1] = -z[1202];
  df_dqdd[2][2] = -z[1209];
  df_dqdd[2][3] = -z[1210];
  df_dqdd[2][4] = -z[1211];
  df_dqdd[2][5] = -z[1212];
  df_dqdd[2][6] = -z[1213];
  df_dqdd[2][7] = -z[1214];
  df_dqdd[2][8] = -z[1215];
  df_dqdd[3][0] = -z[1195];
  df_dqdd[3][1] = -z[1203];
  df_dqdd[3][2] = -z[1210];
  df_dqdd[3][3] = -z[1216];
  df_dqdd[3][4] = -z[1217];
  df_dqdd[3][5] = -z[1218];
  df_dqdd[3][6] = 0;
  df_dqdd[3][7] = 0;
  df_dqdd[3][8] = 0;
  df_dqdd[4][0] = -z[1196];
  df_dqdd[4][1] = -z[1204];
  df_dqdd[4][2] = -z[1211];
  df_dqdd[4][3] = -z[1217];
  df_dqdd[4][4] = -z[1219];
  df_dqdd[4][5] = -z[1220];
  df_dqdd[4][6] = 0;
  df_dqdd[4][7] = 0;
  df_dqdd[4][8] = 0;
  df_dqdd[5][0] = -z[1197];
  df_dqdd[5][1] = -z[1205];
  df_dqdd[5][2] = -z[1212];
  df_dqdd[5][3] = -z[1218];
  df_dqdd[5][4] = -z[1220];
  df_dqdd[5][5] = -z[1221];
  df_dqdd[5][6] = 0;
  df_dqdd[5][7] = 0;
  df_dqdd[5][8] = 0;
  df_dqdd[6][0] = -z[1198];
  df_dqdd[6][1] = -z[1206];
  df_dqdd[6][2] = -z[1213];
  df_dqdd[6][3] = 0;
  df_dqdd[6][4] = 0;
  df_dqdd[6][5] = 0;
  df_dqdd[6][6] = -z[1222];
  df_dqdd[6][7] = -z[1223];
  df_dqdd[6][8] = -z[1224];
  df_dqdd[7][0] = -z[1199];
  df_dqdd[7][1] = -z[1207];
  df_dqdd[7][2] = -z[1214];
  df_dqdd[7][3] = 0;
  df_dqdd[7][4] = 0;
  df_dqdd[7][5] = 0;
  df_dqdd[7][6] = -z[1223];
  df_dqdd[7][7] = -z[1225];
  df_dqdd[7][8] = -z[1226];
  df_dqdd[8][0] = -z[1200];
  df_dqdd[8][1] = -z[1208];
  df_dqdd[8][2] = -z[1215];
  df_dqdd[8][3] = 0;
  df_dqdd[8][4] = 0;
  df_dqdd[8][5] = 0;
  df_dqdd[8][6] = -z[1224];
  df_dqdd[8][7] = -z[1226];
  df_dqdd[8][8] = -z[1227];
  df_dG[0][0] = -z[1246];
  df_dG[0][1] = -z[1247];
  df_dG[0][2] = 1;
  df_dG[0][3] = -z[1248];
  df_dG[0][4] = -z[1249];
  df_dG[0][5] = 1;
  df_dG[1][0] = z[1250];
  df_dG[1][1] = z[1251];
  df_dG[1][2] = 0;
  df_dG[1][3] = z[1252];
  df_dG[1][4] = z[1253];
  df_dG[1][5] = 0;
  df_dG[2][0] = z[1251];
  df_dG[2][1] = z[1254];
  df_dG[2][2] = 0;
  df_dG[2][3] = z[1253];
  df_dG[2][4] = z[1255];
  df_dG[2][5] = 0;
  df_dG[3][0] = -z[1256];
  df_dG[3][1] = -z[1257];
  df_dG[3][2] = 1;
  df_dG[3][3] = 0;
  df_dG[3][4] = 0;
  df_dG[3][5] = 0;
  df_dG[4][0] = -z[1258];
  df_dG[4][1] = -z[1259];
  df_dG[4][2] = 1;
  df_dG[4][3] = 0;
  df_dG[4][4] = 0;
  df_dG[4][5] = 0;
  df_dG[5][0] = z[1260];
  df_dG[5][1] = -z[1261];
  df_dG[5][2] = 1;
  df_dG[5][3] = 0;
  df_dG[5][4] = 0;
  df_dG[5][5] = 0;
  df_dG[6][0] = 0;
  df_dG[6][1] = 0;
  df_dG[6][2] = 0;
  df_dG[6][3] = -z[1262];
  df_dG[6][4] = -z[1263];
  df_dG[6][5] = 1;
  df_dG[7][0] = 0;
  df_dG[7][1] = 0;
  df_dG[7][2] = 0;
  df_dG[7][3] = -z[1264];
  df_dG[7][4] = -z[1265];
  df_dG[7][5] = 1;
  df_dG[8][0] = 0;
  df_dG[8][1] = 0;
  df_dG[8][2] = 0;
  df_dG[8][3] = z[1266];
  df_dG[8][4] = -z[1267];
  df_dG[8][5] = 1;

  Encode[0] = 0.0;
  Encode[1] = 0.0;
  Encode[2] = 0.0;
  Encode[3] = 0.0;

  Encode[4] = 0.0;

  Encode[5] = 0.0;
  Encode[6] = 0.0;
  Encode[7] = 0.0;
  Encode[8] = 0.0;

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


