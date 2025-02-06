/* The name of this program is gait10dof18musc/gait10dof18musc_raw.c */
/* Created by Autolev 4.1 on Thu Nov 28 14:33:18 2024 */

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
  par__torso_Izz,par__torso_M,par__wind,G1,G2,G3,G4,G5,G6,q1,q10,q2,q3,q4,q5,q6,
  q7,q8,q9,q1p,q10p,q2p,q3p,q4p,q5p,q6p,q7p,q8p,q9p,q1pp,q10pp,q2pp,q3pp,q4pp,
  q5pp,q6pp,q7pp,q8pp,q9pp;
double   s,sx,sy;
double   Pi,DEGtoRAD,RADtoDEG,z[1348],fk[48],dfk_dq[48][10],fkdot[48],dfkdot_dq[
  48][10],f[10],df_dq[10][10],df_dqd[10][10],df_dqdd[10][10],df_dG[10][6],
  Encode[9];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait10dof18musc/gait10dof18musc_raw.in");
  else sprintf(fileName, "gait10dof18musc/gait10dof18musc_raw.%d", iloop);
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
  par__wind,&G1,&G2,&G3,&G4,&G5,&G6,&q1,&q10,&q2,&q3,&q4,&q5,&q6,&q7,&q8,&q9,&
  q1p,&q10p,&q2p,&q3p,&q4p,&q5p,&q6p,&q7p,&q8p,&q9p,&q1pp,&q10pp,&q2pp,&q3pp,&
  q4pp,&q5pp,&q6pp,&q7pp,&q8pp,&q9pp,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file gait10dof18musc/gait10dof18musc_raw.in" );
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
  z[97] = cos(q10);
  z[98] = sin(q10);
  z[99] = z[1]*z[97] - z[2]*z[98];
  z[100] = -z[1]*z[98] - z[2]*z[97];
  z[101] = z[1]*z[98] + z[2]*z[97];
  z[102] = par__back_x*z[98] - par__back_y*z[97];
  z[103] = par__back_x*z[97] + par__back_y*z[98];
  z[104] = z[102] - par__torso_CMy;
  z[105] = par__torso_CMx + z[103];
  z[106] = z[2]*z[4] - z[1]*z[3];
  z[107] = par__Rknee_x4 + 2*par__Rknee_x3*q5 + 3*par__Rknee_x2*pow(q5,2) + 4*
  par__Rknee_x1*pow(q5,3);
  z[108] = par__Rknee_y4 + 2*par__Rknee_y3*q5 + 3*par__Rknee_y2*pow(q5,2) + 4*
  par__Rknee_y1*pow(q5,3);
  z[109] = z[6]*z[16] + z[17]*z[106];
  z[110] = z[16]*z[106] - z[6]*z[17];
  z[111] = -z[5]*z[16] - z[6]*z[17];
  z[112] = -z[5]*z[17] - z[7]*z[16];
  z[113] = z[32]*z[109] + z[33]*z[110];
  z[114] = z[19]*z[32] + z[33]*z[111];
  z[115] = z[32]*z[110] - z[33]*z[109];
  z[116] = z[32]*z[111] - z[19]*z[33];
  z[117] = -z[18]*z[32] - z[19]*z[33];
  z[118] = z[21]*z[32] + z[33]*z[112];
  z[119] = z[32]*z[112] - z[21]*z[33];
  z[120] = -z[20]*z[32] - z[21]*z[33];
  z[121] = z[2]*z[51] - z[1]*z[50];
  z[122] = par__Lknee_x4 + 2*par__Lknee_x3*q8 + 3*par__Lknee_x2*pow(q8,2) + 4*
  par__Lknee_x1*pow(q8,3);
  z[123] = par__Lknee_y4 + 2*par__Lknee_y3*q8 + 3*par__Lknee_y2*pow(q8,2) + 4*
  par__Lknee_y1*pow(q8,3);
  z[124] = z[53]*z[63] + z[64]*z[121];
  z[125] = z[63]*z[121] - z[53]*z[64];
  z[126] = -z[52]*z[63] - z[53]*z[64];
  z[127] = -z[52]*z[64] - z[54]*z[63];
  z[128] = z[79]*z[124] + z[80]*z[125];
  z[129] = z[66]*z[79] + z[80]*z[126];
  z[130] = z[79]*z[125] - z[80]*z[124];
  z[131] = z[79]*z[126] - z[66]*z[80];
  z[132] = -z[65]*z[79] - z[66]*z[80];
  z[133] = z[68]*z[79] + z[80]*z[127];
  z[134] = z[79]*z[127] - z[68]*z[80];
  z[135] = -z[67]*z[79] - z[68]*z[80];
  z[136] = z[2]*z[98] - z[1]*z[97];
  z[137] = -par__Rhip_x*z[2] - par__Rhip_y*z[1];
  z[138] = par__Rhip_x*z[1] - par__Rhip_y*z[2];
  z[139] = z[6]*z[12] + z[13]*z[106] - par__Rhip_x*z[2] - par__Rhip_y*z[1];
  z[140] = z[6]*z[12] + z[13]*z[106];
  z[141] = z[5]*z[107] + z[6]*z[108];
  z[142] = par__Rhip_x*z[1] + z[5]*z[12] + z[6]*z[13] - par__Rhip_y*z[2];
  z[143] = z[5]*z[12] + z[6]*z[13];
  z[144] = z[5]*z[108] + z[7]*z[107];
  z[145] = par__Rankle_y*z[110] + z[6]*z[12] + z[13]*z[106] - par__Rhip_x*
  z[2] - par__Rhip_y*z[1];
  z[146] = par__Rankle_y*z[110] + z[6]*z[12] + z[13]*z[106];
  z[147] = par__Rankle_y*z[111] + z[5]*z[107] + z[6]*z[108];
  z[148] = par__Rankle_y*z[19] + par__Rhip_x*z[1] + z[5]*z[12] + z[6]*z[13] - 
  par__Rhip_y*z[2];
  z[149] = par__Rankle_y*z[19] + z[5]*z[12] + z[6]*z[13];
  z[150] = par__Rankle_y*z[112] + z[5]*z[108] + z[7]*z[107];
  z[151] = -par__Lhip_x*z[2] - par__Lhip_y*z[1];
  z[152] = par__Lhip_x*z[1] - par__Lhip_y*z[2];
  z[153] = z[53]*z[59] + z[60]*z[121] - par__Lhip_x*z[2] - par__Lhip_y*z[1];
  z[154] = z[53]*z[59] + z[60]*z[121];
  z[155] = z[52]*z[122] + z[53]*z[123];
  z[156] = par__Lhip_x*z[1] + z[52]*z[59] + z[53]*z[60] - par__Lhip_y*z[2];
  z[157] = z[52]*z[59] + z[53]*z[60];
  z[158] = z[52]*z[123] + z[54]*z[122];
  z[159] = par__Lankle_y*z[125] + z[53]*z[59] + z[60]*z[121] - par__Lhip_x*
  z[2] - par__Lhip_y*z[1];
  z[160] = par__Lankle_y*z[125] + z[53]*z[59] + z[60]*z[121];
  z[161] = par__Lankle_y*z[126] + z[52]*z[122] + z[53]*z[123];
  z[162] = par__Lankle_y*z[66] + par__Lhip_x*z[1] + z[52]*z[59] + z[53]*z[60] - 
  par__Lhip_y*z[2];
  z[163] = par__Lankle_y*z[66] + z[52]*z[59] + z[53]*z[60];
  z[164] = par__Lankle_y*z[127] + z[52]*z[123] + z[54]*z[122];
  z[165] = -par__back_x*z[2] - par__back_y*z[1];
  z[166] = par__back_x*z[1] - par__back_y*z[2];
  z[167] = z[2]*q1p;
  z[168] = z[1]*q1p;
  z[169] = q2p - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[170] = q3p + par__Rhip_x*z[1]*q1p - par__Rhip_y*z[2]*q1p;
  z[171] = z[6]*(q1p+q4p);
  z[172] = z[106]*(q1p+q4p);
  z[173] = z[5]*(q1p+q4p);
  z[174] = q2p + z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + 
  z[13]*z[106]*(q1p+q4p) - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[175] = q3p + par__Rhip_x*z[1]*q1p + z[5]*z[108]*q5p + z[7]*z[107]*q5p + 
  z[5]*z[12]*(q1p+q4p) + z[6]*z[13]*(q1p+q4p) - par__Rhip_y*z[2]*q1p;
  z[176] = z[19]*q5p + z[109]*q1p + z[109]*q4p;
  z[177] = z[110]*q1p + z[110]*q4p + z[111]*q5p;
  z[178] = z[18]*q1p + z[18]*q4p + z[21]*q5p;
  z[179] = z[19]*q1p + z[19]*q4p + z[112]*q5p;
  z[180] = q2p + z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + 
  z[13]*z[106]*(q1p+q4p) + par__Rankle_y*(z[110]*q1p+z[110]*q4p+z[111]*q5p) - 
  par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[181] = q3p + par__Rhip_x*z[1]*q1p + z[5]*z[108]*q5p + z[7]*z[107]*q5p + 
  z[5]*z[12]*(q1p+q4p) + z[6]*z[13]*(q1p+q4p) + par__Rankle_y*(z[19]*q1p+
  z[19]*q4p+z[112]*q5p) - par__Rhip_y*z[2]*q1p;
  z[182] = z[35]*q6p + z[113]*q1p + z[113]*q4p + z[114]*q5p;
  z[183] = z[115]*q1p + z[115]*q4p + z[116]*q5p + z[117]*q6p;
  z[184] = z[34]*q1p + z[34]*q4p + z[37]*q6p + z[118]*q5p;
  z[185] = z[35]*q1p + z[35]*q4p + z[119]*q5p + z[120]*q6p;
  z[186] = q2p - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[187] = q3p + par__Lhip_x*z[1]*q1p - par__Lhip_y*z[2]*q1p;
  z[188] = z[53]*(q1p+q7p);
  z[189] = z[121]*(q1p+q7p);
  z[190] = z[52]*(q1p+q7p);
  z[191] = q2p + z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p) - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[192] = q3p + par__Lhip_x*z[1]*q1p + z[52]*z[123]*q8p + z[54]*z[122]*q8p + 
  z[52]*z[59]*(q1p+q7p) + z[53]*z[60]*(q1p+q7p) - par__Lhip_y*z[2]*q1p;
  z[193] = z[66]*q8p + z[124]*q1p + z[124]*q7p;
  z[194] = z[125]*q1p + z[125]*q7p + z[126]*q8p;
  z[195] = z[65]*q1p + z[65]*q7p + z[68]*q8p;
  z[196] = z[66]*q1p + z[66]*q7p + z[127]*q8p;
  z[197] = q2p + z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p) + par__Lankle_y*(z[125]*q1p+z[125]*q7p+z[126]*q8p) - 
  par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[198] = q3p + par__Lhip_x*z[1]*q1p + z[52]*z[123]*q8p + z[54]*z[122]*q8p + 
  z[52]*z[59]*(q1p+q7p) + z[53]*z[60]*(q1p+q7p) + par__Lankle_y*(z[66]*q1p+
  z[66]*q7p+z[127]*q8p) - par__Lhip_y*z[2]*q1p;
  z[199] = z[82]*q9p + z[128]*q1p + z[128]*q7p + z[129]*q8p;
  z[200] = z[130]*q1p + z[130]*q7p + z[131]*q8p + z[132]*q9p;
  z[201] = z[81]*q1p + z[81]*q7p + z[84]*q9p + z[133]*q8p;
  z[202] = z[82]*q1p + z[82]*q7p + z[134]*q8p + z[135]*q9p;
  z[203] = q2p - par__back_x*z[2]*q1p - par__back_y*z[1]*q1p;
  z[204] = q3p + par__back_x*z[1]*q1p - par__back_y*z[2]*q1p;
  z[205] = z[100]*(q1p+q10p);
  z[206] = z[136]*(q1p+q10p);
  z[207] = z[99]*(q1p+q10p);
  z[208] = (par__Rhip_x*z[1]-par__Rhip_y*z[2])*q1p;
  z[209] = (par__Rhip_x*z[2]+par__Rhip_y*z[1])*q1p;
  z[210] = par__Rhip_y*z[2]*q1p + z[6]*z[107]*q5p + z[106]*z[108]*q5p + z[7]*
  z[13]*(q1p+q4p) + z[12]*z[106]*(q1p+q4p) - par__Rhip_x*z[1]*q1p;
  z[211] = z[6]*z[107]*q5p + z[106]*z[108]*q5p + z[7]*z[13]*(q1p+q4p) + z[12]*
  z[106]*(q1p+q4p);
  z[212] = 2*par__Rknee_x3 + 6*par__Rknee_x2*q5 + 12*par__Rknee_x1*pow(q5,2);
  z[213] = 2*par__Rknee_y3 + 6*par__Rknee_y2*q5 + 12*par__Rknee_y1*pow(q5,2);
  z[214] = z[5]*z[212]*q5p + z[6]*z[213]*q5p + z[6]*z[107]*(q1p+q4p) + z[106]*
  z[108]*(q1p+q4p);
  z[215] = z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[106]*(q1p+q4p) - par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[216] = z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[106]*(q1p+q4p);
  z[217] = z[5]*z[213]*q5p + z[7]*z[212]*q5p + z[5]*z[107]*(q1p+q4p) + z[6]*
  z[108]*(q1p+q4p);
  z[218] = z[7]*z[17] + z[16]*z[106];
  z[219] = z[110]*q5p + z[218]*q1p + z[218]*q4p;
  z[220] = z[7]*z[16] - z[17]*z[106];
  z[221] = -z[6]*z[16] - z[17]*z[106];
  z[222] = z[220]*q1p + z[220]*q4p + z[221]*q5p;
  z[223] = z[5]*z[17] - z[6]*z[16];
  z[224] = z[221]*q1p + z[221]*q4p + z[223]*q5p;
  z[225] = z[7]*z[17] - z[5]*z[16];
  z[226] = z[111]*q1p + z[111]*q4p + z[225]*q5p;
  z[227] = par__Rhip_y*z[2]*q1p + z[6]*z[107]*q5p + z[106]*z[108]*q5p + z[7]*
  z[13]*(q1p+q4p) + z[12]*z[106]*(q1p+q4p) + par__Rankle_y*(z[220]*q1p+z[220]*
  q4p+z[221]*q5p) - par__Rhip_x*z[1]*q1p;
  z[228] = z[6]*z[107]*q5p + z[106]*z[108]*q5p + z[7]*z[13]*(q1p+q4p) + z[12]*
  z[106]*(q1p+q4p) + par__Rankle_y*(z[220]*q1p+z[220]*q4p+z[221]*q5p);
  z[229] = z[5]*z[212]*q5p + z[6]*z[213]*q5p + z[6]*z[107]*(q1p+q4p) + z[106]*
  z[108]*(q1p+q4p) + par__Rankle_y*(z[221]*q1p+z[221]*q4p+z[223]*q5p);
  z[230] = z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[106]*(q1p+q4p) + par__Rankle_y*(z[110]*q1p+z[110]*q4p+z[111]*q5p) - 
  par__Rhip_x*z[2]*q1p - par__Rhip_y*z[1]*q1p;
  z[231] = z[5]*z[107]*q5p + z[6]*z[108]*q5p + z[6]*z[12]*(q1p+q4p) + z[13]*
  z[106]*(q1p+q4p) + par__Rankle_y*(z[110]*q1p+z[110]*q4p+z[111]*q5p);
  z[232] = z[5]*z[213]*q5p + z[7]*z[212]*q5p + z[5]*z[107]*(q1p+q4p) + z[6]*
  z[108]*(q1p+q4p) + par__Rankle_y*(z[111]*q1p+z[111]*q4p+z[225]*q5p);
  z[233] = z[32]*z[218] + z[33]*z[220];
  z[234] = z[32]*z[110] + z[33]*z[221];
  z[235] = z[115]*q6p + z[233]*q1p + z[233]*q4p + z[234]*q5p;
  z[236] = z[32]*z[111] + z[33]*z[223];
  z[237] = z[116]*q6p + z[234]*q1p + z[234]*q4p + z[236]*q5p;
  z[238] = z[32]*z[220] - z[33]*z[218];
  z[239] = z[32]*z[221] - z[33]*z[110];
  z[240] = -z[32]*z[109] - z[33]*z[110];
  z[241] = z[238]*q1p + z[238]*q4p + z[239]*q5p + z[240]*q6p;
  z[242] = z[32]*z[223] - z[33]*z[111];
  z[243] = -z[19]*z[32] - z[33]*z[111];
  z[244] = z[239]*q1p + z[239]*q4p + z[242]*q5p + z[243]*q6p;
  z[245] = z[18]*z[33] - z[19]*z[32];
  z[246] = z[240]*q1p + z[240]*q4p + z[243]*q5p + z[245]*q6p;
  z[247] = z[32]*z[112] + z[33]*z[225];
  z[248] = z[114]*q1p + z[114]*q4p + z[119]*q6p + z[247]*q5p;
  z[249] = z[32]*z[225] - z[33]*z[112];
  z[250] = -z[21]*z[32] - z[33]*z[112];
  z[251] = z[116]*q1p + z[116]*q4p + z[249]*q5p + z[250]*q6p;
  z[252] = z[20]*z[33] - z[21]*z[32];
  z[253] = z[117]*q1p + z[117]*q4p + z[250]*q5p + z[252]*q6p;
  z[254] = (par__Lhip_x*z[1]-par__Lhip_y*z[2])*q1p;
  z[255] = (par__Lhip_x*z[2]+par__Lhip_y*z[1])*q1p;
  z[256] = par__Lhip_y*z[2]*q1p + z[53]*z[122]*q8p + z[121]*z[123]*q8p + 
  z[54]*z[60]*(q1p+q7p) + z[59]*z[121]*(q1p+q7p) - par__Lhip_x*z[1]*q1p;
  z[257] = z[53]*z[122]*q8p + z[121]*z[123]*q8p + z[54]*z[60]*(q1p+q7p) + 
  z[59]*z[121]*(q1p+q7p);
  z[258] = 2*par__Lknee_x3 + 6*par__Lknee_x2*q8 + 12*par__Lknee_x1*pow(q8,2);
  z[259] = 2*par__Lknee_y3 + 6*par__Lknee_y2*q8 + 12*par__Lknee_y1*pow(q8,2);
  z[260] = z[52]*z[258]*q8p + z[53]*z[259]*q8p + z[53]*z[122]*(q1p+q7p) + 
  z[121]*z[123]*(q1p+q7p);
  z[261] = z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p) - par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[262] = z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p);
  z[263] = z[52]*z[259]*q8p + z[54]*z[258]*q8p + z[52]*z[122]*(q1p+q7p) + 
  z[53]*z[123]*(q1p+q7p);
  z[264] = z[54]*z[64] + z[63]*z[121];
  z[265] = z[125]*q8p + z[264]*q1p + z[264]*q7p;
  z[266] = z[54]*z[63] - z[64]*z[121];
  z[267] = -z[53]*z[63] - z[64]*z[121];
  z[268] = z[266]*q1p + z[266]*q7p + z[267]*q8p;
  z[269] = z[52]*z[64] - z[53]*z[63];
  z[270] = z[267]*q1p + z[267]*q7p + z[269]*q8p;
  z[271] = z[54]*z[64] - z[52]*z[63];
  z[272] = z[126]*q1p + z[126]*q7p + z[271]*q8p;
  z[273] = par__Lhip_y*z[2]*q1p + z[53]*z[122]*q8p + z[121]*z[123]*q8p + 
  z[54]*z[60]*(q1p+q7p) + z[59]*z[121]*(q1p+q7p) + par__Lankle_y*(z[266]*q1p+
  z[266]*q7p+z[267]*q8p) - par__Lhip_x*z[1]*q1p;
  z[274] = z[53]*z[122]*q8p + z[121]*z[123]*q8p + z[54]*z[60]*(q1p+q7p) + 
  z[59]*z[121]*(q1p+q7p) + par__Lankle_y*(z[266]*q1p+z[266]*q7p+z[267]*q8p);
  z[275] = z[52]*z[258]*q8p + z[53]*z[259]*q8p + z[53]*z[122]*(q1p+q7p) + 
  z[121]*z[123]*(q1p+q7p) + par__Lankle_y*(z[267]*q1p+z[267]*q7p+z[269]*q8p);
  z[276] = z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p) + par__Lankle_y*(z[125]*q1p+z[125]*q7p+z[126]*q8p) - 
  par__Lhip_x*z[2]*q1p - par__Lhip_y*z[1]*q1p;
  z[277] = z[52]*z[122]*q8p + z[53]*z[123]*q8p + z[53]*z[59]*(q1p+q7p) + 
  z[60]*z[121]*(q1p+q7p) + par__Lankle_y*(z[125]*q1p+z[125]*q7p+z[126]*q8p);
  z[278] = z[52]*z[259]*q8p + z[54]*z[258]*q8p + z[52]*z[122]*(q1p+q7p) + 
  z[53]*z[123]*(q1p+q7p) + par__Lankle_y*(z[126]*q1p+z[126]*q7p+z[271]*q8p);
  z[279] = z[79]*z[264] + z[80]*z[266];
  z[280] = z[79]*z[125] + z[80]*z[267];
  z[281] = z[130]*q9p + z[279]*q1p + z[279]*q7p + z[280]*q8p;
  z[282] = z[79]*z[126] + z[80]*z[269];
  z[283] = z[131]*q9p + z[280]*q1p + z[280]*q7p + z[282]*q8p;
  z[284] = z[79]*z[266] - z[80]*z[264];
  z[285] = z[79]*z[267] - z[80]*z[125];
  z[286] = -z[79]*z[124] - z[80]*z[125];
  z[287] = z[284]*q1p + z[284]*q7p + z[285]*q8p + z[286]*q9p;
  z[288] = z[79]*z[269] - z[80]*z[126];
  z[289] = -z[66]*z[79] - z[80]*z[126];
  z[290] = z[285]*q1p + z[285]*q7p + z[288]*q8p + z[289]*q9p;
  z[291] = z[65]*z[80] - z[66]*z[79];
  z[292] = z[286]*q1p + z[286]*q7p + z[289]*q8p + z[291]*q9p;
  z[293] = z[79]*z[127] + z[80]*z[271];
  z[294] = z[129]*q1p + z[129]*q7p + z[134]*q9p + z[293]*q8p;
  z[295] = z[79]*z[271] - z[80]*z[127];
  z[296] = -z[68]*z[79] - z[80]*z[127];
  z[297] = z[131]*q1p + z[131]*q7p + z[295]*q8p + z[296]*q9p;
  z[298] = z[67]*z[80] - z[68]*z[79];
  z[299] = z[132]*q1p + z[132]*q7p + z[296]*q8p + z[298]*q9p;
  z[300] = (par__back_x*z[1]-par__back_y*z[2])*q1p;
  z[301] = (par__back_x*z[2]+par__back_y*z[1])*q1p;
  z[302] = z[7]*(q1p+q4p);
  z[303] = z[54]*(q1p+q7p);
  z[304] = z[101]*(q1p+q10p);
  z[305] = par__gravity_x*par__Lfemur_M;
  z[306] = par__gravity_y*par__Lfemur_M;
  z[307] = par__gravity_x*par__Lfoot_M;
  z[308] = par__gravity_y*par__Lfoot_M;
  z[309] = par__gravity_x*par__Ltibia_M;
  z[310] = par__gravity_y*par__Ltibia_M;
  z[311] = par__gravity_x*par__pelvis_M;
  z[312] = par__gravity_y*par__pelvis_M;
  z[313] = par__gravity_x*par__Rfemur_M;
  z[314] = par__gravity_y*par__Rfemur_M;
  z[315] = par__gravity_x*par__Rfoot_M;
  z[316] = par__gravity_y*par__Rfoot_M;
  z[317] = par__gravity_x*par__Rtibia_M;
  z[318] = par__gravity_y*par__Rtibia_M;
  z[319] = par__gravity_x*par__torso_M;
  z[320] = par__gravity_y*par__torso_M;
  sx = -par__wind - z[1]*(par__pelvis_CMy*q1p-z[1]*q2p) - z[2]*(par__pelvis_CMx*
  q1p-z[2]*q2p);
  sy = z[1]*(par__pelvis_CMx*q1p+z[1]*q3p) - z[2]*(par__pelvis_CMy*q1p-z[2]*
  q3p);
  s = pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[321] = par__airdrag*par__wind;
  z[322] = z[321]*s;
  z[323] = par__airdrag*s*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*q3p);
  z[324] = par__airdrag*s*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*q2p);
  z[325] = z[311] + z[322];
  z[333] = par__Rfoot_CMy + z[44];
  z[334] = par__Rfoot_CMy + z[45];
  z[335] = par__Rfoot_CMy + z[46];
  z[336] = z[47] - par__Rfoot_CMx;
  z[337] = z[48] - par__Rfoot_CMx;
  z[338] = z[49] - par__Rfoot_CMx;
  z[339] = par__bodyweight*G1;
  z[340] = par__bodyweight*G2;
  z[342] = par__bodyweight*G3;
  z[350] = par__Lfoot_CMy + z[91];
  z[351] = par__Lfoot_CMy + z[92];
  z[352] = par__Lfoot_CMy + z[93];
  z[353] = z[94] - par__Lfoot_CMx;
  z[354] = z[95] - par__Lfoot_CMx;
  z[355] = z[96] - par__Lfoot_CMx;
  z[356] = par__bodyweight*G4;
  z[357] = par__bodyweight*G5;
  z[360] = par__bodyweight*G6;
  z[362] = par__pelvis_CMx*z[312];
  z[363] = par__pelvis_CMy*z[312];
  z[364] = z[342] + z[360] + z[362]*z[1] + q3*z[339] + q3*z[356] + par__Lankle_y*
  z[65]*z[356] + par__Lankle_y*z[67]*z[357] + par__Lhip_x*z[2]*z[356] + 
  par__Lhip_y*z[1]*z[356] + par__Lhip_y*z[2]*z[357] + par__Rankle_y*z[18]*
  z[339] + par__Rankle_y*z[20]*z[340] + par__Rhip_x*z[2]*z[339] + par__Rhip_y*
  z[1]*z[339] + par__Rhip_y*z[2]*z[340] + z[305]*z[52]*z[57] + z[305]*z[53]*
  z[58] + z[306]*z[52]*z[58] + z[306]*z[54]*z[57] + z[307]*z[81]*z[92] + 
  z[307]*z[82]*z[95] + z[308]*z[83]*z[92] + z[308]*z[84]*z[95] + z[309]*z[65]*
  z[74] + z[309]*z[66]*z[76] + z[310]*z[67]*z[74] + z[310]*z[68]*z[76] + 
  z[313]*z[5]*z[10] + z[313]*z[6]*z[11] + z[314]*z[5]*z[11] + z[314]*z[7]*
  z[10] + z[315]*z[34]*z[45] + z[315]*z[35]*z[48] + z[316]*z[36]*z[45] + 
  z[316]*z[37]*z[48] + z[317]*z[18]*z[27] + z[317]*z[19]*z[29] + z[318]*z[20]*
  z[27] + z[318]*z[21]*z[29] + z[319]*z[99]*z[104] + z[319]*z[100]*z[105] + 
  z[320]*z[99]*z[105] + z[320]*z[101]*z[104] + z[5]*z[13]*z[339] + z[7]*z[13]*
  z[340] + z[34]*z[334]*z[339] + z[35]*z[337]*z[339] + z[36]*z[334]*z[340] + 
  z[37]*z[337]*z[340] + z[52]*z[60]*z[356] + z[54]*z[60]*z[357] + z[81]*
  z[351]*z[356] + z[82]*z[354]*z[356] + z[83]*z[351]*z[357] + z[84]*z[354]*
  z[357] - z[363]*z[2] - q2*z[340] - q2*z[357] - par__Lhip_x*z[1]*z[357] - 
  par__pelvis_CMx*z[2]*z[325] - par__pelvis_CMy*z[1]*z[325] - par__Rhip_x*
  z[1]*z[340] - z[5]*z[12]*z[340] - z[6]*z[12]*z[339] - z[52]*z[59]*z[357] - 
  z[53]*z[59]*z[356] - par__pelvis_CMx*z[324] - par__pelvis_CMy*z[323];
  z[365] = z[305]*pow(z[52],2) + z[305]*pow(z[53],2) + z[306]*z[52]*z[53] + 
  z[306]*z[52]*z[54] + z[307]*pow(z[81],2) + z[307]*pow(z[82],2) + z[308]*
  z[81]*z[83] + z[308]*z[82]*z[84] + z[309]*pow(z[65],2) + z[309]*pow(z[66],2) + 
  z[310]*z[65]*z[67] + z[310]*z[66]*z[68] + z[313]*pow(z[5],2) + z[313]*pow(
  z[6],2) + z[314]*z[5]*z[6] + z[314]*z[5]*z[7] + z[315]*pow(z[34],2) + 
  z[315]*pow(z[35],2) + z[316]*z[34]*z[36] + z[316]*z[35]*z[37] + z[317]*pow(
  z[18],2) + z[317]*pow(z[19],2) + z[318]*z[18]*z[20] + z[318]*z[19]*z[21] + 
  z[319]*pow(z[99],2) + z[319]*pow(z[100],2) + z[320]*z[99]*z[100] + z[320]*
  z[99]*z[101] + z[34]*z[36]*z[340] + z[35]*z[37]*z[340] + z[81]*z[83]*z[357] + 
  z[82]*z[84]*z[357] + z[325]*pow(z[1],2) + z[325]*pow(z[2],2) + z[339]*pow(
  z[34],2) + z[339]*pow(z[35],2) + z[356]*pow(z[81],2) + z[356]*pow(z[82],2) + 
  z[1]*z[323] + z[2]*z[324];
  z[366] = z[305]*z[52]*z[53] + z[305]*z[52]*z[54] + z[306]*pow(z[52],2) + 
  z[306]*pow(z[54],2) + z[307]*z[81]*z[83] + z[307]*z[82]*z[84] + z[308]*pow(
  z[83],2) + z[308]*pow(z[84],2) + z[309]*z[65]*z[67] + z[309]*z[66]*z[68] + 
  z[310]*pow(z[67],2) + z[310]*pow(z[68],2) + z[312]*pow(z[1],2) + z[312]*
  pow(z[2],2) + z[313]*z[5]*z[6] + z[313]*z[5]*z[7] + z[314]*pow(z[5],2) + 
  z[314]*pow(z[7],2) + z[315]*z[34]*z[36] + z[315]*z[35]*z[37] + z[316]*pow(
  z[36],2) + z[316]*pow(z[37],2) + z[317]*z[18]*z[20] + z[317]*z[19]*z[21] + 
  z[318]*pow(z[20],2) + z[318]*pow(z[21],2) + z[319]*z[99]*z[100] + z[319]*
  z[99]*z[101] + z[320]*pow(z[99],2) + z[320]*pow(z[101],2) + z[34]*z[36]*
  z[339] + z[35]*z[37]*z[339] + z[81]*z[83]*z[356] + z[82]*z[84]*z[356] + 
  z[340]*pow(z[36],2) + z[340]*pow(z[37],2) + z[357]*pow(z[83],2) + z[357]*
  pow(z[84],2) + z[2]*z[323] - z[1]*z[324];
  z[367] = par__Rfemur_CMx*z[313];
  z[368] = par__Rfemur_CMx*z[314];
  z[369] = par__Rfemur_CMy*z[313];
  z[370] = par__Rfemur_CMy*z[314];
  z[371] = z[342] + z[367]*z[6] + z[368]*z[5] + q3*z[339] + par__Rankle_y*
  z[18]*z[339] + par__Rankle_y*z[20]*z[340] + par__Rhip_x*z[2]*z[339] + 
  par__Rhip_y*z[1]*z[339] + par__Rhip_y*z[2]*z[340] + z[315]*z[34]*z[44] + 
  z[315]*z[35]*z[47] + z[316]*z[36]*z[44] + z[316]*z[37]*z[47] + z[317]*z[18]*
  z[26] + z[317]*z[19]*z[28] + z[318]*z[20]*z[26] + z[318]*z[21]*z[28] + z[5]*
  z[13]*z[339] + z[7]*z[13]*z[340] + z[34]*z[333]*z[339] + z[35]*z[336]*
  z[339] + z[36]*z[333]*z[340] + z[37]*z[336]*z[340] - z[369]*z[5] - z[370]*
  z[7] - q2*z[340] - par__Rhip_x*z[1]*z[340] - z[5]*z[12]*z[340] - z[6]*z[12]*
  z[339];
  z[372] = par__Rtibia_CMx*z[317];
  z[373] = par__Rtibia_CMx*z[318];
  z[374] = par__Rtibia_CMy*z[317];
  z[375] = par__Rtibia_CMy*z[318];
  z[376] = z[342] + z[372]*z[19] + z[373]*z[21] + q3*z[339] + par__Rankle_y*
  z[18]*z[339] + par__Rankle_y*z[20]*z[340] + par__Rhip_x*z[2]*z[339] + 
  par__Rhip_y*z[1]*z[339] + par__Rhip_y*z[2]*z[340] + z[315]*z[34]*z[46] + 
  z[315]*z[35]*z[49] + z[316]*z[36]*z[46] + z[316]*z[37]*z[49] + z[5]*z[13]*
  z[339] + z[7]*z[13]*z[340] + z[34]*z[335]*z[339] + z[35]*z[338]*z[339] + 
  z[36]*z[335]*z[340] + z[37]*z[338]*z[340] - z[374]*z[18] - z[375]*z[20] - 
  q2*z[340] - par__Rhip_x*z[1]*z[340] - z[5]*z[12]*z[340] - z[6]*z[12]*z[339];
  z[377] = par__Rfoot_CMx*z[315];
  z[378] = par__Rfoot_CMx*z[316];
  z[379] = par__Rfoot_CMy*z[315];
  z[380] = par__Rfoot_CMy*z[316];
  z[381] = z[342] + z[377]*z[35] + z[378]*z[37] + q3*z[339] + par__Rankle_y*
  z[18]*z[339] + par__Rankle_y*z[20]*z[340] + par__Rhip_x*z[2]*z[339] + 
  par__Rhip_y*z[1]*z[339] + par__Rhip_y*z[2]*z[340] + z[5]*z[13]*z[339] + 
  z[7]*z[13]*z[340] - z[379]*z[34] - z[380]*z[36] - q2*z[340] - par__Rhip_x*
  z[1]*z[340] - z[5]*z[12]*z[340] - z[6]*z[12]*z[339];
  z[382] = par__Lfemur_CMx*z[305];
  z[383] = par__Lfemur_CMx*z[306];
  z[384] = par__Lfemur_CMy*z[305];
  z[385] = par__Lfemur_CMy*z[306];
  z[386] = z[360] + z[382]*z[53] + z[383]*z[52] + q3*z[356] + par__Lankle_y*
  z[65]*z[356] + par__Lankle_y*z[67]*z[357] + par__Lhip_x*z[2]*z[356] + 
  par__Lhip_y*z[1]*z[356] + par__Lhip_y*z[2]*z[357] + z[307]*z[81]*z[91] + 
  z[307]*z[82]*z[94] + z[308]*z[83]*z[91] + z[308]*z[84]*z[94] + z[309]*z[65]*
  z[73] + z[309]*z[66]*z[75] + z[310]*z[67]*z[73] + z[310]*z[68]*z[75] + 
  z[52]*z[60]*z[356] + z[54]*z[60]*z[357] + z[81]*z[350]*z[356] + z[82]*
  z[353]*z[356] + z[83]*z[350]*z[357] + z[84]*z[353]*z[357] - z[384]*z[52] - 
  z[385]*z[54] - q2*z[357] - par__Lhip_x*z[1]*z[357] - z[52]*z[59]*z[357] - 
  z[53]*z[59]*z[356];
  z[387] = par__Ltibia_CMx*z[309];
  z[388] = par__Ltibia_CMx*z[310];
  z[389] = par__Ltibia_CMy*z[309];
  z[390] = par__Ltibia_CMy*z[310];
  z[391] = z[360] + z[387]*z[66] + z[388]*z[68] + q3*z[356] + par__Lankle_y*
  z[65]*z[356] + par__Lankle_y*z[67]*z[357] + par__Lhip_x*z[2]*z[356] + 
  par__Lhip_y*z[1]*z[356] + par__Lhip_y*z[2]*z[357] + z[307]*z[81]*z[93] + 
  z[307]*z[82]*z[96] + z[308]*z[83]*z[93] + z[308]*z[84]*z[96] + z[52]*z[60]*
  z[356] + z[54]*z[60]*z[357] + z[81]*z[352]*z[356] + z[82]*z[355]*z[356] + 
  z[83]*z[352]*z[357] + z[84]*z[355]*z[357] - z[389]*z[65] - z[390]*z[67] - 
  q2*z[357] - par__Lhip_x*z[1]*z[357] - z[52]*z[59]*z[357] - z[53]*z[59]*
  z[356];
  z[392] = par__Lfoot_CMx*z[307];
  z[393] = par__Lfoot_CMx*z[308];
  z[394] = par__Lfoot_CMy*z[307];
  z[395] = par__Lfoot_CMy*z[308];
  z[396] = z[360] + z[392]*z[82] + z[393]*z[84] + q3*z[356] + par__Lankle_y*
  z[65]*z[356] + par__Lankle_y*z[67]*z[357] + par__Lhip_x*z[2]*z[356] + 
  par__Lhip_y*z[1]*z[356] + par__Lhip_y*z[2]*z[357] + z[52]*z[60]*z[356] + 
  z[54]*z[60]*z[357] - z[394]*z[81] - z[395]*z[83] - q2*z[357] - par__Lhip_x*
  z[1]*z[357] - z[52]*z[59]*z[357] - z[53]*z[59]*z[356];
  z[397] = par__torso_CMx*z[319];
  z[398] = par__torso_CMx*z[320];
  z[399] = par__torso_CMy*z[319];
  z[400] = par__torso_CMy*z[320];
  z[401] = z[397]*z[100] + z[398]*z[99] - z[399]*z[99] - z[400]*z[101];
  z[402] = q1p + q7p;
  z[403] = q1p + q7p + q8p + q9p;
  z[404] = q1p + q7p + q8p;
  z[405] = q1p + q4p;
  z[406] = q1p + q4p + q5p + q6p;
  z[407] = q1p + q4p + q5p;
  z[408] = q1p + q10p;
  z[409] = z[52]*q2p + z[54]*q3p + z[57]*q1p - par__Lfemur_CMy*q7p;
  z[410] = par__Lfemur_CMx*q7p + z[52]*q3p + z[53]*q2p + z[58]*q1p;
  z[411] = -z[1]*z[51]*q1p - z[1]*z[51]*q7p - z[2]*z[50]*q1p - z[2]*z[50]*q7p;
  z[412] = z[1]*z[50]*q1p + z[1]*z[50]*q7p - z[2]*z[51]*q1p - z[2]*z[51]*q7p;
  z[413] = (par__Lhip_x*z[50]+par__Lhip_y*z[51])*q7p;
  z[414] = q1p*z[413] + q2p*z[411] + q3p*z[412];
  z[415] = z[2]*z[51]*q1p + z[2]*z[51]*q7p - z[1]*z[50]*q1p - z[1]*z[50]*q7p;
  z[416] = (par__Lhip_x*z[51]-par__Lhip_y*z[50])*q7p;
  z[417] = q2p*z[415] + q3p*z[411] - q1p*z[416];
  z[418] = z[414] - z[402]*z[410];
  z[419] = z[417] + z[402]*z[409];
  z[420] = z[81]*q2p + z[83]*q3p + z[91]*q7p + z[92]*q1p + z[93]*q8p - par__Lfoot_CMy*
  q9p;
  z[421] = par__Lfoot_CMx*q9p + z[82]*q2p + z[84]*q3p + z[94]*q7p + z[95]*q1p + 
  z[96]*q8p;
  z[422] = z[53]*z[63]*q8p + z[63]*z[411] + z[64]*z[415] - z[52]*z[64]*q8p;
  z[423] = z[63]*z[415] - z[52]*z[63]*q8p - z[53]*z[64]*q8p - z[64]*z[411];
  z[424] = z[66]*z[79]*q9p + z[79]*z[422] + z[80]*z[423] - z[65]*z[80]*q9p;
  z[425] = z[52]*z[63]*q8p + z[63]*z[412] + z[64]*z[411] - z[54]*z[64]*q8p;
  z[426] = z[63]*z[411] - z[52]*z[64]*q8p - z[54]*z[63]*q8p - z[64]*z[412];
  z[427] = z[68]*z[79]*q9p + z[79]*z[425] + z[80]*z[426] - z[67]*z[80]*q9p;
  z[428] = (z[59]*z[64]-z[60]*z[63]-z[63]*z[122]-z[64]*z[123])*q8p;
  z[429] = (z[63]*z[123]-z[59]*z[63]-z[60]*z[64]-z[64]*z[122])*q8p;
  z[430] = z[71]*z[79]*q9p - z[77]*z[80]*q9p - z[79]*z[429] - z[80]*z[428];
  z[431] = z[122]*q8p - z[416];
  z[432] = z[413] - z[123]*q8p;
  z[433] = z[63]*z[431] - z[61]*z[63]*q8p - z[62]*z[64]*q8p - z[64]*z[432];
  z[434] = z[62]*z[63]*q8p + z[63]*z[432] + z[64]*z[431] - z[61]*z[64]*q8p;
  z[435] = z[72]*z[79]*q9p + z[79]*z[434] + z[80]*z[433] - z[78]*z[80]*q9p;
  z[436] = par__Lankle_y*z[80]*q9p;
  z[437] = q1p*z[435] + q2p*z[424] + q3p*z[427] + q7p*z[430] + q8p*z[436];
  z[438] = z[79]*z[423] - z[65]*z[79]*q9p - z[66]*z[80]*q9p - z[80]*z[422];
  z[439] = z[79]*z[426] - z[67]*z[79]*q9p - z[68]*z[80]*q9p - z[80]*z[425];
  z[440] = z[80]*z[429] - z[71]*z[80]*q9p - z[77]*z[79]*q9p - z[79]*z[428];
  z[441] = z[79]*z[433] - z[72]*z[80]*q9p - z[78]*z[79]*q9p - z[80]*z[434];
  z[442] = par__Lankle_y*z[79]*q9p;
  z[443] = q1p*z[441] + q2p*z[438] + q3p*z[439] + q7p*z[440] + q8p*z[442];
  z[444] = z[437] - z[403]*z[421];
  z[445] = z[443] + z[403]*z[420];
  z[446] = z[65]*q2p + z[67]*q3p + z[73]*q7p + z[74]*q1p - par__Ltibia_CMy*
  q8p;
  z[447] = par__Ltibia_CMx*q8p + z[66]*q2p + z[68]*q3p + z[75]*q7p + z[76]*
  q1p;
  z[448] = q1p*z[434] + q2p*z[422] + q3p*z[425] - q7p*z[429];
  z[449] = q1p*z[433] + q2p*z[423] + q3p*z[426] - q7p*z[428];
  z[450] = z[448] - z[404]*z[447];
  z[451] = z[449] + z[404]*z[446];
  z[452] = z[1]*q2p + z[2]*q3p - par__pelvis_CMy*q1p;
  z[453] = par__pelvis_CMx*q1p + z[1]*q3p - z[2]*q2p;
  z[454] = q1p*(z[1]*q3p-z[2]*q2p);
  z[455] = q1p*(z[1]*q2p+z[2]*q3p);
  z[456] = z[454] - q1p*z[453];
  z[457] = q1p*z[452] - z[455];
  z[458] = z[5]*q2p + z[7]*q3p + z[10]*q1p - par__Rfemur_CMy*q4p;
  z[459] = par__Rfemur_CMx*q4p + z[5]*q3p + z[6]*q2p + z[11]*q1p;
  z[460] = -z[1]*z[4]*q1p - z[1]*z[4]*q4p - z[2]*z[3]*q1p - z[2]*z[3]*q4p;
  z[461] = z[1]*z[3]*q1p + z[1]*z[3]*q4p - z[2]*z[4]*q1p - z[2]*z[4]*q4p;
  z[462] = (par__Rhip_x*z[3]+par__Rhip_y*z[4])*q4p;
  z[463] = q1p*z[462] + q2p*z[460] + q3p*z[461];
  z[464] = z[2]*z[4]*q1p + z[2]*z[4]*q4p - z[1]*z[3]*q1p - z[1]*z[3]*q4p;
  z[465] = (par__Rhip_x*z[4]-par__Rhip_y*z[3])*q4p;
  z[466] = q2p*z[464] + q3p*z[460] - q1p*z[465];
  z[467] = z[463] - z[405]*z[459];
  z[468] = z[466] + z[405]*z[458];
  z[469] = z[34]*q2p + z[36]*q3p + z[44]*q4p + z[45]*q1p + z[46]*q5p - par__Rfoot_CMy*
  q6p;
  z[470] = par__Rfoot_CMx*q6p + z[35]*q2p + z[37]*q3p + z[47]*q4p + z[48]*q1p + 
  z[49]*q5p;
  z[471] = z[6]*z[16]*q5p + z[16]*z[460] + z[17]*z[464] - z[5]*z[17]*q5p;
  z[472] = z[16]*z[464] - z[5]*z[16]*q5p - z[6]*z[17]*q5p - z[17]*z[460];
  z[473] = z[19]*z[32]*q6p + z[32]*z[471] + z[33]*z[472] - z[18]*z[33]*q6p;
  z[474] = z[5]*z[16]*q5p + z[16]*z[461] + z[17]*z[460] - z[7]*z[17]*q5p;
  z[475] = z[16]*z[460] - z[5]*z[17]*q5p - z[7]*z[16]*q5p - z[17]*z[461];
  z[476] = z[21]*z[32]*q6p + z[32]*z[474] + z[33]*z[475] - z[20]*z[33]*q6p;
  z[477] = (z[12]*z[17]-z[13]*z[16]-z[16]*z[107]-z[17]*z[108])*q5p;
  z[478] = (z[16]*z[108]-z[12]*z[16]-z[13]*z[17]-z[17]*z[107])*q5p;
  z[479] = z[24]*z[32]*q6p - z[30]*z[33]*q6p - z[32]*z[478] - z[33]*z[477];
  z[480] = z[107]*q5p - z[465];
  z[481] = z[462] - z[108]*q5p;
  z[482] = z[16]*z[480] - z[14]*z[16]*q5p - z[15]*z[17]*q5p - z[17]*z[481];
  z[483] = z[15]*z[16]*q5p + z[16]*z[481] + z[17]*z[480] - z[14]*z[17]*q5p;
  z[484] = z[25]*z[32]*q6p + z[32]*z[483] + z[33]*z[482] - z[31]*z[33]*q6p;
  z[485] = par__Rankle_y*z[33]*q6p;
  z[486] = q1p*z[484] + q2p*z[473] + q3p*z[476] + q4p*z[479] + q5p*z[485];
  z[487] = z[32]*z[472] - z[18]*z[32]*q6p - z[19]*z[33]*q6p - z[33]*z[471];
  z[488] = z[32]*z[475] - z[20]*z[32]*q6p - z[21]*z[33]*q6p - z[33]*z[474];
  z[489] = z[33]*z[478] - z[24]*z[33]*q6p - z[30]*z[32]*q6p - z[32]*z[477];
  z[490] = z[32]*z[482] - z[25]*z[33]*q6p - z[31]*z[32]*q6p - z[33]*z[483];
  z[491] = par__Rankle_y*z[32]*q6p;
  z[492] = q1p*z[490] + q2p*z[487] + q3p*z[488] + q4p*z[489] + q5p*z[491];
  z[493] = z[486] - z[406]*z[470];
  z[494] = z[492] + z[406]*z[469];
  z[495] = z[18]*q2p + z[20]*q3p + z[26]*q4p + z[27]*q1p - par__Rtibia_CMy*
  q5p;
  z[496] = par__Rtibia_CMx*q5p + z[19]*q2p + z[21]*q3p + z[28]*q4p + z[29]*
  q1p;
  z[497] = q1p*z[483] + q2p*z[471] + q3p*z[474] - q4p*z[478];
  z[498] = q1p*z[482] + q2p*z[472] + q3p*z[475] - q4p*z[477];
  z[499] = z[497] - z[407]*z[496];
  z[500] = z[498] + z[407]*z[495];
  z[501] = z[99]*q2p + z[101]*q3p + z[104]*q1p - par__torso_CMy*q10p;
  z[502] = par__torso_CMx*q10p + z[99]*q3p + z[100]*q2p + z[105]*q1p;
  z[503] = -z[1]*z[98]*q1p - z[1]*z[98]*q10p - z[2]*z[97]*q1p - z[2]*z[97]*
  q10p;
  z[504] = z[1]*z[97]*q1p + z[1]*z[97]*q10p - z[2]*z[98]*q1p - z[2]*z[98]*
  q10p;
  z[505] = (par__back_x*z[97]+par__back_y*z[98])*q10p;
  z[506] = q1p*z[505] + q2p*z[503] + q3p*z[504];
  z[507] = z[2]*z[98]*q1p + z[2]*z[98]*q10p - z[1]*z[97]*q1p - z[1]*z[97]*
  q10p;
  z[508] = (par__back_x*z[98]-par__back_y*z[97])*q10p;
  z[509] = q2p*z[507] + q3p*z[503] - q1p*z[508];
  z[510] = z[506] - z[408]*z[502];
  z[511] = z[509] + z[408]*z[501];
  z[512] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz + par__pelvis_Izz + 
  par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz + par__torso_Izz + 
  par__pelvis_M*(pow(par__pelvis_CMx,2)+pow(par__pelvis_CMy,2));
  z[513] = z[512] + par__Lfemur_M*(pow(z[57],2)+pow(z[58],2)) + par__Lfoot_M*(
  pow(z[92],2)+pow(z[95],2)) + par__Ltibia_M*(pow(z[74],2)+pow(z[76],2)) + 
  par__Rfemur_M*(pow(z[10],2)+pow(z[11],2)) + par__Rfoot_M*(pow(z[45],2)+pow(
  z[48],2)) + par__Rtibia_M*(pow(z[27],2)+pow(z[29],2)) + par__torso_M*(pow(
  z[104],2)+pow(z[105],2));
  z[514] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz;
  z[515] = z[514] + par__Lfoot_M*(z[91]*z[92]+z[94]*z[95]) + par__Ltibia_M*(
  z[73]*z[74]+z[75]*z[76]) + par__Lfemur_M*(par__Lfemur_CMx*z[58]-par__Lfemur_CMy*
  z[57]);
  z[516] = par__Lfoot_Izz + par__Ltibia_Izz;
  z[517] = z[516] + par__Lfoot_M*(z[92]*z[93]+z[95]*z[96]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[76]-par__Ltibia_CMy*z[74]);
  z[518] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[95]-par__Lfoot_CMy*
  z[92]);
  z[519] = par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz;
  z[520] = z[519] + par__Rfoot_M*(z[44]*z[45]+z[47]*z[48]) + par__Rtibia_M*(
  z[26]*z[27]+z[28]*z[29]) + par__Rfemur_M*(par__Rfemur_CMx*z[11]-par__Rfemur_CMy*
  z[10]);
  z[521] = par__Rfoot_Izz + par__Rtibia_Izz;
  z[522] = z[521] + par__Rfoot_M*(z[45]*z[46]+z[48]*z[49]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[29]-par__Rtibia_CMy*z[27]);
  z[523] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[48]-par__Rfoot_CMy*
  z[45]);
  z[524] = par__torso_Izz + par__torso_M*(par__torso_CMx*z[105]-par__torso_CMy*
  z[104]);
  z[525] = par__Lfemur_M*(z[52]*z[58]+z[54]*z[57]) + par__Lfoot_M*(z[83]*
  z[92]+z[84]*z[95]) + par__Ltibia_M*(z[67]*z[74]+z[68]*z[76]) + par__Rfemur_M*(
  z[5]*z[11]+z[7]*z[10]) + par__Rfoot_M*(z[36]*z[45]+z[37]*z[48]) + par__Rtibia_M*(
  z[20]*z[27]+z[21]*z[29]) + par__torso_M*(z[99]*z[105]+z[101]*z[104]) + 
  par__pelvis_M*(par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2]);
  z[526] = par__Lfemur_M*(z[52]*z[57]+z[53]*z[58]) + par__Lfoot_M*(z[81]*
  z[92]+z[82]*z[95]) + par__Ltibia_M*(z[65]*z[74]+z[66]*z[76]) + par__Rfemur_M*(
  z[5]*z[10]+z[6]*z[11]) + par__Rfoot_M*(z[34]*z[45]+z[35]*z[48]) + par__Rtibia_M*(
  z[18]*z[27]+z[19]*z[29]) + par__torso_M*(z[99]*z[104]+z[100]*z[105]) - 
  par__pelvis_M*(par__pelvis_CMx*z[2]+par__pelvis_CMy*z[1]);
  z[527] = par__Lfemur_M*(z[57]*z[418]+z[58]*z[419]) + par__Lfoot_M*(z[92]*
  z[444]+z[95]*z[445]) + par__Ltibia_M*(z[74]*z[450]+z[76]*z[451]) + par__Rfemur_M*(
  z[10]*z[467]+z[11]*z[468]) + par__Rfoot_M*(z[45]*z[493]+z[48]*z[494]) + 
  par__Rtibia_M*(z[27]*z[499]+z[29]*z[500]) + par__torso_M*(z[104]*z[510]+
  z[105]*z[511]) + par__pelvis_M*(par__pelvis_CMx*z[457]-par__pelvis_CMy*
  z[456]);
  z[528] = par__Lfoot_M*(z[81]*z[91]+z[82]*z[94]) + par__Ltibia_M*(z[65]*
  z[73]+z[66]*z[75]) + par__Lfemur_M*(par__Lfemur_CMx*z[53]-par__Lfemur_CMy*
  z[52]);
  z[529] = par__Lfemur_M*z[52]*(z[53]+z[54]) + par__Rfemur_M*z[5]*(z[6]+z[7]) + 
  par__torso_M*z[99]*(z[100]+z[101]) + par__Lfoot_M*(z[81]*z[83]+z[82]*z[84]) + 
  par__Ltibia_M*(z[65]*z[67]+z[66]*z[68]) + par__Rfoot_M*(z[34]*z[36]+z[35]*
  z[37]) + par__Rtibia_M*(z[18]*z[20]+z[19]*z[21]);
  z[530] = par__Lfemur_M*(pow(z[52],2)+pow(z[53],2)) + par__Lfoot_M*(pow(
  z[81],2)+pow(z[82],2)) + par__Ltibia_M*(pow(z[65],2)+pow(z[66],2)) + par__pelvis_M*(
  pow(z[1],2)+pow(z[2],2)) + par__Rfemur_M*(pow(z[5],2)+pow(z[6],2)) + par__Rfoot_M*(
  pow(z[34],2)+pow(z[35],2)) + par__Rtibia_M*(pow(z[18],2)+pow(z[19],2)) + 
  par__torso_M*(pow(z[99],2)+pow(z[100],2));
  z[531] = par__Lfoot_M*(par__Lfoot_CMx*z[82]-par__Lfoot_CMy*z[81]);
  z[532] = par__Lfoot_M*(z[81]*z[93]+z[82]*z[96]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[66]-par__Ltibia_CMy*z[65]);
  z[533] = par__Rfoot_M*(z[34]*z[44]+z[35]*z[47]) + par__Rtibia_M*(z[18]*
  z[26]+z[19]*z[28]) + par__Rfemur_M*(par__Rfemur_CMx*z[6]-par__Rfemur_CMy*
  z[5]);
  z[534] = par__Rfoot_M*(par__Rfoot_CMx*z[35]-par__Rfoot_CMy*z[34]);
  z[535] = par__Rfoot_M*(z[34]*z[46]+z[35]*z[49]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[19]-par__Rtibia_CMy*z[18]);
  z[536] = par__torso_M*(par__torso_CMx*z[100]-par__torso_CMy*z[99]);
  z[537] = par__Lfemur_M*(z[52]*z[418]+z[53]*z[419]) + par__Lfoot_M*(z[81]*
  z[444]+z[82]*z[445]) + par__Ltibia_M*(z[65]*z[450]+z[66]*z[451]) + par__Rfemur_M*(
  z[5]*z[467]+z[6]*z[468]) + par__Rfoot_M*(z[34]*z[493]+z[35]*z[494]) + 
  par__Rtibia_M*(z[18]*z[499]+z[19]*z[500]) + par__torso_M*(z[99]*z[510]+
  z[100]*z[511]) + par__pelvis_M*(z[1]*z[456]-z[2]*z[457]);
  z[538] = par__Lfoot_M*(z[83]*z[91]+z[84]*z[94]) + par__Ltibia_M*(z[67]*
  z[73]+z[68]*z[75]) + par__Lfemur_M*(par__Lfemur_CMx*z[52]-par__Lfemur_CMy*
  z[54]);
  z[539] = par__Lfemur_M*(pow(z[52],2)+pow(z[54],2)) + par__Lfoot_M*(pow(
  z[83],2)+pow(z[84],2)) + par__Ltibia_M*(pow(z[67],2)+pow(z[68],2)) + par__pelvis_M*(
  pow(z[1],2)+pow(z[2],2)) + par__Rfemur_M*(pow(z[5],2)+pow(z[7],2)) + par__Rfoot_M*(
  pow(z[36],2)+pow(z[37],2)) + par__Rtibia_M*(pow(z[20],2)+pow(z[21],2)) + 
  par__torso_M*(pow(z[99],2)+pow(z[101],2));
  z[540] = par__Lfoot_M*(par__Lfoot_CMx*z[84]-par__Lfoot_CMy*z[83]);
  z[541] = par__Lfoot_M*(z[83]*z[93]+z[84]*z[96]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[68]-par__Ltibia_CMy*z[67]);
  z[542] = par__Rfoot_M*(z[36]*z[44]+z[37]*z[47]) + par__Rtibia_M*(z[20]*
  z[26]+z[21]*z[28]) + par__Rfemur_M*(par__Rfemur_CMx*z[5]-par__Rfemur_CMy*
  z[7]);
  z[543] = par__Rfoot_M*(par__Rfoot_CMx*z[37]-par__Rfoot_CMy*z[36]);
  z[544] = par__Rfoot_M*(z[36]*z[46]+z[37]*z[49]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[21]-par__Rtibia_CMy*z[20]);
  z[545] = par__torso_M*(par__torso_CMx*z[99]-par__torso_CMy*z[101]);
  z[546] = par__Lfemur_M*(z[52]*z[419]+z[54]*z[418]) + par__Lfoot_M*(z[83]*
  z[444]+z[84]*z[445]) + par__Ltibia_M*(z[67]*z[450]+z[68]*z[451]) + par__pelvis_M*(
  z[1]*z[457]+z[2]*z[456]) + par__Rfemur_M*(z[5]*z[468]+z[7]*z[467]) + par__Rfoot_M*(
  z[36]*z[493]+z[37]*z[494]) + par__Rtibia_M*(z[20]*z[499]+z[21]*z[500]) + 
  par__torso_M*(z[99]*z[511]+z[101]*z[510]);
  z[547] = par__Rfemur_Izz + par__Rfoot_Izz + par__Rtibia_Izz + par__Rfemur_M*(
  pow(par__Rfemur_CMx,2)+pow(par__Rfemur_CMy,2));
  z[548] = z[547] + par__Rfoot_M*(pow(z[44],2)+pow(z[47],2)) + par__Rtibia_M*(
  pow(z[26],2)+pow(z[28],2));
  z[549] = z[521] + par__Rfoot_M*(z[44]*z[46]+z[47]*z[49]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[28]-par__Rtibia_CMy*z[26]);
  z[550] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[47]-par__Rfoot_CMy*
  z[44]);
  z[551] = par__Rfoot_M*(z[44]*z[493]+z[47]*z[494]) + par__Rtibia_M*(z[26]*
  z[499]+z[28]*z[500]) + par__Rfemur_M*(par__Rfemur_CMx*z[468]-par__Rfemur_CMy*
  z[467]);
  z[552] = par__Rfoot_Izz + par__Rtibia_Izz + par__Rtibia_M*(pow(par__Rtibia_CMx,
  2)+pow(par__Rtibia_CMy,2));
  z[553] = z[552] + par__Rfoot_M*(pow(z[46],2)+pow(z[49],2));
  z[554] = par__Rfoot_Izz + par__Rfoot_M*(par__Rfoot_CMx*z[49]-par__Rfoot_CMy*
  z[46]);
  z[555] = par__Rfoot_M*(z[46]*z[493]+z[49]*z[494]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[500]-par__Rtibia_CMy*z[499]);
  z[556] = par__Rfoot_Izz + par__Rfoot_M*(pow(par__Rfoot_CMx,2)+pow(par__Rfoot_CMy,
  2));
  z[557] = par__Rfoot_M*(par__Rfoot_CMx*z[494]-par__Rfoot_CMy*z[493]);
  z[558] = par__Lfemur_Izz + par__Lfoot_Izz + par__Ltibia_Izz + par__Lfemur_M*(
  pow(par__Lfemur_CMx,2)+pow(par__Lfemur_CMy,2));
  z[559] = z[558] + par__Lfoot_M*(pow(z[91],2)+pow(z[94],2)) + par__Ltibia_M*(
  pow(z[73],2)+pow(z[75],2));
  z[560] = z[516] + par__Lfoot_M*(z[91]*z[93]+z[94]*z[96]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[75]-par__Ltibia_CMy*z[73]);
  z[561] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[94]-par__Lfoot_CMy*
  z[91]);
  z[562] = par__Lfoot_M*(z[91]*z[444]+z[94]*z[445]) + par__Ltibia_M*(z[73]*
  z[450]+z[75]*z[451]) + par__Lfemur_M*(par__Lfemur_CMx*z[419]-par__Lfemur_CMy*
  z[418]);
  z[563] = par__Lfoot_Izz + par__Ltibia_Izz + par__Ltibia_M*(pow(par__Ltibia_CMx,
  2)+pow(par__Ltibia_CMy,2));
  z[564] = z[563] + par__Lfoot_M*(pow(z[93],2)+pow(z[96],2));
  z[565] = par__Lfoot_Izz + par__Lfoot_M*(par__Lfoot_CMx*z[96]-par__Lfoot_CMy*
  z[93]);
  z[566] = par__Lfoot_M*(z[93]*z[444]+z[96]*z[445]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[451]-par__Ltibia_CMy*z[450]);
  z[567] = par__Lfoot_Izz + par__Lfoot_M*(pow(par__Lfoot_CMx,2)+pow(par__Lfoot_CMy,
  2));
  z[568] = par__Lfoot_M*(par__Lfoot_CMx*z[445]-par__Lfoot_CMy*z[444]);
  z[569] = par__torso_Izz + par__torso_M*(pow(par__torso_CMx,2)+pow(par__torso_CMy,
  2));
  z[570] = par__torso_M*(par__torso_CMx*z[511]-par__torso_CMy*z[510]);
  z[571] = (par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2])*q1p;
  z[572] = (par__pelvis_CMx*z[2]+par__pelvis_CMy*z[1])*q1p;
  z[573] = (sx*z[571]+sy*z[572])/pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[574] = par__airdrag*(s*(z[1]*q2p+z[2]*q3p)+(par__pelvis_CMx*q1p+z[1]*q3p-
  z[2]*q2p)*z[573]);
  z[575] = par__airdrag*(s*(z[1]*q3p-z[2]*q2p)+(par__pelvis_CMy*q1p-z[1]*q2p-
  z[2]*q3p)*z[573]);
  z[576] = par__pelvis_CMx*z[321];
  z[577] = par__pelvis_CMy*z[321];
  z[578] = par__Lankle_y*z[65]*z[357] + par__Lankle_y*z[124]*z[356] + par__Lhip_x*
  z[1]*z[356] + par__Lhip_x*z[2]*z[357] + par__Lhip_y*z[1]*z[357] + par__pelvis_CMy*
  z[2]*z[325] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*z[109]*z[339] + 
  par__Rhip_x*z[1]*z[339] + par__Rhip_x*z[2]*z[340] + par__Rhip_y*z[1]*z[340] + 
  z[305]*z[53]*z[57] + z[305]*z[58]*z[121] + z[306]*z[52]*z[57] + z[306]*
  z[53]*z[58] + z[307]*z[92]*z[128] + z[307]*z[95]*z[130] + z[308]*z[81]*
  z[92] + z[308]*z[82]*z[95] + z[309]*z[74]*z[124] + z[309]*z[76]*z[125] + 
  z[310]*z[65]*z[74] + z[310]*z[66]*z[76] + z[313]*z[6]*z[10] + z[313]*z[11]*
  z[106] + z[314]*z[5]*z[10] + z[314]*z[6]*z[11] + z[315]*z[45]*z[113] + 
  z[315]*z[48]*z[115] + z[316]*z[34]*z[45] + z[316]*z[35]*z[48] + z[317]*
  z[27]*z[109] + z[317]*z[29]*z[110] + z[318]*z[18]*z[27] + z[318]*z[19]*
  z[29] + z[319]*z[100]*z[104] + z[319]*z[105]*z[136] + z[320]*z[99]*z[104] + 
  z[320]*z[100]*z[105] + z[5]*z[13]*z[340] + z[6]*z[13]*z[339] + z[34]*z[334]*
  z[340] + z[35]*z[337]*z[340] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] + 
  z[81]*z[351]*z[357] + z[82]*z[354]*z[357] + z[113]*z[334]*z[339] + z[115]*
  z[337]*z[339] + z[128]*z[351]*z[356] + z[130]*z[354]*z[356] + par__pelvis_CMx*
  z[574] + par__pelvis_CMy*z[575] + z[576]*z[2]*z[573] + z[577]*z[1]*z[573] - 
  z[362]*z[2] - z[363]*z[1] - par__Lhip_y*z[2]*z[356] - par__pelvis_CMx*z[1]*
  z[325] - par__Rhip_y*z[2]*z[339] - z[6]*z[12]*z[340] - z[12]*z[106]*z[339] - 
  z[53]*z[59]*z[357] - z[59]*z[121]*z[356];
  z[579] = q2p*z[415] + q3p*z[411];
  z[580] = z[53]*q3p + z[121]*q2p;
  z[581] = z[1]*z[51]*q1p + z[1]*z[51]*q7p + z[2]*z[50]*q1p + z[2]*z[50]*q7p;
  z[582] = q2p*z[581] + q3p*z[415];
  z[583] = z[52]*q3p + z[53]*q2p;
  z[584] = z[63]*z[121]*q8p + z[63]*z[415] + z[64]*z[581] - z[53]*z[64]*q8p;
  z[585] = z[63]*z[581] - z[53]*z[63]*q8p - z[64]*z[121]*q8p - z[64]*z[415];
  z[586] = z[79]*z[125]*q9p + z[79]*z[584] + z[80]*z[585] - z[80]*z[124]*q9p;
  z[587] = q2p*z[586] + q3p*z[424];
  z[588] = z[82]*q3p + z[130]*q2p;
  z[589] = z[79]*z[585] - z[79]*z[124]*q9p - z[80]*z[125]*q9p - z[80]*z[584];
  z[590] = q2p*z[589] + q3p*z[438];
  z[591] = z[81]*q3p + z[128]*q2p;
  z[592] = q2p*z[584] + q3p*z[422];
  z[593] = z[66]*q3p + z[125]*q2p;
  z[594] = q2p*z[585] + q3p*z[423];
  z[595] = z[65]*q3p + z[124]*q2p;
  z[596] = q2p*z[464] + q3p*z[460];
  z[597] = z[6]*q3p + z[106]*q2p;
  z[598] = z[1]*z[4]*q1p + z[1]*z[4]*q4p + z[2]*z[3]*q1p + z[2]*z[3]*q4p;
  z[599] = q2p*z[598] + q3p*z[464];
  z[600] = z[5]*q3p + z[6]*q2p;
  z[601] = z[16]*z[106]*q5p + z[16]*z[464] + z[17]*z[598] - z[6]*z[17]*q5p;
  z[602] = z[16]*z[598] - z[6]*z[16]*q5p - z[17]*z[106]*q5p - z[17]*z[464];
  z[603] = z[32]*z[110]*q6p + z[32]*z[601] + z[33]*z[602] - z[33]*z[109]*q6p;
  z[604] = q2p*z[603] + q3p*z[473];
  z[605] = z[35]*q3p + z[115]*q2p;
  z[606] = z[32]*z[602] - z[32]*z[109]*q6p - z[33]*z[110]*q6p - z[33]*z[601];
  z[607] = q2p*z[606] + q3p*z[487];
  z[608] = z[34]*q3p + z[113]*q2p;
  z[609] = q2p*z[601] + q3p*z[471];
  z[610] = z[19]*q3p + z[110]*q2p;
  z[611] = q2p*z[602] + q3p*z[472];
  z[612] = z[18]*q3p + z[109]*q2p;
  z[613] = q2p*z[507] + q3p*z[503];
  z[614] = z[100]*q3p + z[136]*q2p;
  z[615] = z[1]*z[98]*q1p + z[1]*z[98]*q10p + z[2]*z[97]*q1p + z[2]*z[97]*
  q10p;
  z[616] = q2p*z[615] + q3p*z[507];
  z[617] = z[99]*q3p + z[100]*q2p;
  z[618] = z[1]*q3p - z[2]*q2p;
  z[619] = -z[1]*q2p - z[2]*q3p;
  z[620] = par__Lfemur_M*(z[58]*(z[582]+z[402]*z[583])+z[57]*(z[579]-z[402]*
  z[580])) + par__Lfoot_M*(z[95]*(z[590]+z[403]*z[591])+z[92]*(z[587]-z[403]*
  z[588])) + par__Ltibia_M*(z[76]*(z[594]+z[404]*z[595])+z[74]*(z[592]-z[404]*
  z[593])) + par__Rfemur_M*(z[11]*(z[599]+z[405]*z[600])+z[10]*(z[596]-z[405]*
  z[597])) + par__Rfoot_M*(z[48]*(z[607]+z[406]*z[608])+z[45]*(z[604]-z[406]*
  z[605])) + par__Rtibia_M*(z[29]*(z[611]+z[407]*z[612])+z[27]*(z[609]-z[407]*
  z[610])) + par__torso_M*(z[105]*(z[616]+z[408]*z[617])+z[104]*(z[613]-
  z[408]*z[614])) + par__pelvis_M*(par__pelvis_CMy*(z[455]+q1p*z[619])-par__pelvis_CMx*(
  z[454]-q1p*z[618]));
  z[621] = par__Lfemur_M*(z[53]*z[57]+z[58]*z[121]) + par__Lfoot_M*(z[92]*
  z[128]+z[95]*z[130]) + par__Ltibia_M*(z[74]*z[124]+z[76]*z[125]) + par__Rfemur_M*(
  z[6]*z[10]+z[11]*z[106]) + par__Rfoot_M*(z[45]*z[113]+z[48]*z[115]) + 
  par__Rtibia_M*(z[27]*z[109]+z[29]*z[110]) + par__torso_M*(z[100]*z[104]+
  z[105]*z[136]) - par__pelvis_M*(par__pelvis_CMx*z[1]-par__pelvis_CMy*z[2]);
  z[622] = par__Rhip_y*z[3] - par__Rhip_x*z[4];
  z[623] = z[16]*z[622] - z[9]*z[17];
  z[624] = z[9]*z[16] + z[17]*z[622];
  z[625] = z[32]*z[624] + z[33]*z[623];
  z[626] = z[32]*z[623] - z[33]*z[624];
  z[627] = par__Rankle_y*z[18]*z[340] + par__Rankle_y*z[109]*z[339] + z[313]*
  z[5]*z[9] + z[313]*z[6]*z[10] + z[313]*z[6]*z[622] + z[313]*z[11]*z[106] + 
  z[314]*z[5]*z[10] + z[314]*z[5]*z[622] + z[314]*z[6]*z[11] + z[314]*z[7]*
  z[9] + z[315]*z[34]*z[625] + z[315]*z[35]*z[626] + z[315]*z[45]*z[113] + 
  z[315]*z[48]*z[115] + z[316]*z[34]*z[45] + z[316]*z[35]*z[48] + z[316]*
  z[36]*z[625] + z[316]*z[37]*z[626] + z[317]*z[18]*z[624] + z[317]*z[19]*
  z[623] + z[317]*z[27]*z[109] + z[317]*z[29]*z[110] + z[318]*z[18]*z[27] + 
  z[318]*z[19]*z[29] + z[318]*z[20]*z[624] + z[318]*z[21]*z[623] + z[5]*z[13]*
  z[340] + z[6]*z[13]*z[339] + z[34]*z[334]*z[340] + z[34]*z[339]*z[625] + 
  z[35]*z[337]*z[340] + z[35]*z[339]*z[626] + z[36]*z[340]*z[625] + z[37]*
  z[340]*z[626] + z[113]*z[334]*z[339] + z[115]*z[337]*z[339] - z[6]*z[12]*
  z[340] - z[12]*z[106]*z[339];
  z[628] = z[6]*q3p + z[106]*q2p + z[622]*q1p;
  z[629] = q2p*z[598] + q3p*z[464] - q1p*z[462];
  z[630] = z[5]*q3p + z[6]*q2p + z[9]*q1p;
  z[631] = z[16]*z[622]*q5p - z[9]*z[17]*q5p - z[16]*z[465] - z[17]*z[462];
  z[632] = z[17]*z[465] - z[9]*z[16]*q5p - z[17]*z[622]*q5p - z[16]*z[462];
  z[633] = z[32]*z[623]*q6p + z[32]*z[631] + z[33]*z[632] - z[33]*z[624]*q6p;
  z[634] = q1p*z[633] + q2p*z[603] + q3p*z[473];
  z[635] = z[35]*q3p + z[115]*q2p + z[626]*q1p;
  z[636] = z[32]*z[632] - z[32]*z[624]*q6p - z[33]*z[623]*q6p - z[33]*z[631];
  z[637] = q1p*z[636] + q2p*z[606] + q3p*z[487];
  z[638] = z[34]*q3p + z[113]*q2p + z[625]*q1p;
  z[639] = q1p*z[631] + q2p*z[601] + q3p*z[471];
  z[640] = z[19]*q3p + z[110]*q2p + z[623]*q1p;
  z[641] = q1p*z[632] + q2p*z[602] + q3p*z[472];
  z[642] = z[18]*q3p + z[109]*q2p + z[624]*q1p;
  z[643] = par__Rfemur_M*(z[9]*z[467]+z[622]*z[468]+z[11]*(z[629]+z[405]*
  z[630])+z[10]*(z[466]-z[405]*z[628])) + par__Rfoot_M*(z[625]*z[493]+z[626]*
  z[494]+z[48]*(z[637]+z[406]*z[638])+z[45]*(z[634]-z[406]*z[635])) + par__Rtibia_M*(
  z[623]*z[500]+z[624]*z[499]+z[29]*(z[641]+z[407]*z[642])+z[27]*(z[639]-
  z[407]*z[640]));
  z[644] = 2*par__Rfemur_M*(z[9]*z[10]+z[11]*z[622]) + 2*par__Rfoot_M*(z[45]*
  z[625]+z[48]*z[626]) + 2*par__Rtibia_M*(z[27]*z[624]+z[29]*z[623]);
  z[645] = par__Rfoot_M*(z[44]*z[625]+z[47]*z[626]) + par__Rtibia_M*(z[26]*
  z[624]+z[28]*z[623]) + par__Rfemur_M*(par__Rfemur_CMx*z[622]-par__Rfemur_CMy*
  z[9]);
  z[646] = par__Rfoot_M*(z[46]*z[625]+z[49]*z[626]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[623]-par__Rtibia_CMy*z[624]);
  z[647] = par__Rfoot_M*(par__Rfoot_CMx*z[626]-par__Rfoot_CMy*z[625]);
  z[648] = par__Rfemur_M*(z[5]*z[10]+z[5]*z[622]+z[6]*z[11]+z[7]*z[9]) + 
  par__Rfoot_M*(z[34]*z[45]+z[35]*z[48]+z[36]*z[625]+z[37]*z[626]) + par__Rtibia_M*(
  z[18]*z[27]+z[19]*z[29]+z[20]*z[624]+z[21]*z[623]);
  z[649] = par__Rfemur_M*(z[5]*z[9]+z[6]*z[10]+z[6]*z[622]+z[11]*z[106]) + 
  par__Rfoot_M*(z[34]*z[625]+z[35]*z[626]+z[45]*z[113]+z[48]*z[115]) + par__Rtibia_M*(
  z[18]*z[624]+z[19]*z[623]+z[27]*z[109]+z[29]*z[110]);
  z[650] = z[16]*z[107] + z[17]*z[108] - z[14]*z[16] - z[15]*z[17];
  z[651] = z[15]*z[16] + z[17]*z[107] - z[14]*z[17] - z[16]*z[108];
  z[652] = z[32]*z[651] + z[33]*z[650];
  z[653] = z[32]*z[650] - z[33]*z[651];
  z[654] = par__Rankle_y*z[19]*z[339] + par__Rankle_y*z[21]*z[340] + z[315]*
  z[34]*z[652] + z[315]*z[35]*z[653] + z[315]*z[45]*z[114] + z[315]*z[48]*
  z[116] + z[316]*z[36]*z[652] + z[316]*z[37]*z[653] + z[316]*z[45]*z[118] + 
  z[316]*z[48]*z[119] + z[317]*z[18]*z[651] + z[317]*z[19]*z[27] + z[317]*
  z[19]*z[650] + z[317]*z[29]*z[111] + z[318]*z[20]*z[651] + z[318]*z[21]*
  z[27] + z[318]*z[21]*z[650] + z[318]*z[29]*z[112] + z[5]*z[108]*z[339] + 
  z[7]*z[108]*z[340] + z[34]*z[339]*z[652] + z[35]*z[339]*z[653] + z[36]*
  z[340]*z[652] + z[37]*z[340]*z[653] + z[114]*z[334]*z[339] + z[116]*z[337]*
  z[339] + z[118]*z[334]*z[340] + z[119]*z[337]*z[340] - z[5]*z[107]*z[340] - 
  z[6]*z[107]*z[339];
  z[655] = z[16]*z[107]*q5p + z[17]*z[108]*q5p + z[17]*z[212]*q5p + z[16]*
  z[480] - z[14]*z[16]*q5p - z[15]*z[17]*q5p - z[16]*z[213]*q5p - z[17]*
  z[481];
  z[656] = z[14]*z[17]*q5p + z[16]*z[108]*q5p + z[16]*z[212]*q5p + z[17]*
  z[213]*q5p - z[15]*z[16]*q5p - z[17]*z[107]*q5p - z[16]*z[481] - z[17]*
  z[480];
  z[657] = z[32]*z[650]*q6p + z[32]*z[655] + z[33]*z[656] - z[33]*z[651]*q6p;
  z[658] = z[5]*z[17]*q5p - z[6]*z[16]*q5p - z[16]*z[460] - z[17]*z[464];
  z[659] = z[32]*z[111]*q6p + z[32]*z[472] + z[33]*z[658] - z[19]*z[33]*q6p;
  z[660] = z[7]*z[17]*q5p - z[5]*z[16]*q5p - z[16]*z[461] - z[17]*z[460];
  z[661] = z[32]*z[112]*q6p + z[32]*z[475] + z[33]*z[660] - z[21]*z[33]*q6p;
  z[662] = z[13]*z[16] + z[16]*z[107] + z[17]*z[108] - z[12]*z[17];
  z[663] = z[12]*z[16] + z[13]*z[17] + z[17]*z[107] - z[16]*z[108];
  z[664] = (z[12]*z[17]+z[16]*z[213]-2*z[16]*z[107]-2*z[17]*z[108]-z[13]*
  z[16]-z[17]*z[212])*q5p;
  z[665] = (z[12]*z[16]+z[13]*z[17]+2*z[17]*z[107]-2*z[16]*z[108]-z[16]*
  z[212]-z[17]*z[213])*q5p;
  z[666] = z[32]*z[662]*q6p - z[33]*z[663]*q6p - z[32]*z[664] - z[33]*z[665];
  z[667] = q1p*z[657] + q2p*z[659] + q3p*z[661] + q4p*z[666];
  z[668] = z[32]*z[662] - z[33]*z[663];
  z[669] = z[116]*q2p + z[119]*q3p + z[653]*q1p + z[668]*q4p;
  z[670] = z[32]*z[656] - z[32]*z[651]*q6p - z[33]*z[650]*q6p - z[33]*z[655];
  z[671] = z[32]*z[658] - z[19]*z[32]*q6p - z[33]*z[111]*q6p - z[33]*z[472];
  z[672] = z[32]*z[660] - z[21]*z[32]*q6p - z[33]*z[112]*q6p - z[33]*z[475];
  z[673] = z[33]*z[664] - z[32]*z[663]*q6p - z[33]*z[662]*q6p - z[32]*z[665];
  z[674] = q1p*z[670] + q2p*z[671] + q3p*z[672] + q4p*z[673];
  z[675] = z[32]*z[663] + z[33]*z[662];
  z[676] = z[114]*q2p + z[118]*q3p + z[652]*q1p + z[675]*q4p;
  z[677] = q1p*z[655] + q2p*z[472] + q3p*z[475] - q4p*z[664];
  z[678] = z[111]*q2p + z[112]*q3p + z[650]*q1p + z[662]*q4p;
  z[679] = q1p*z[656] + q2p*z[658] + q3p*z[660] - q4p*z[665];
  z[680] = z[19]*q2p + z[21]*q3p + z[651]*q1p + z[663]*q4p;
  z[681] = par__Rfoot_M*(z[652]*z[493]+z[653]*z[494]+z[48]*(z[674]+z[406]*
  z[676])+z[45]*(z[667]-z[406]*z[669])) + par__Rtibia_M*(z[650]*z[500]+z[651]*
  z[499]+z[29]*(z[679]+z[407]*z[680])+z[27]*(z[677]-z[407]*z[678]));
  z[682] = 2*par__Rfoot_M*(z[45]*z[652]+z[48]*z[653]) + 2*par__Rtibia_M*(
  z[27]*z[651]+z[29]*z[650]);
  z[683] = par__Rfoot_M*(z[44]*z[652]+z[45]*z[675]+z[47]*z[653]+z[48]*z[668]) + 
  par__Rtibia_M*(z[26]*z[651]+z[27]*z[663]+z[28]*z[650]+z[29]*z[662]);
  z[684] = par__Rfoot_M*(z[46]*z[652]+z[49]*z[653]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[650]-par__Rtibia_CMy*z[651]);
  z[685] = par__Rfoot_M*(par__Rfoot_CMx*z[653]-par__Rfoot_CMy*z[652]);
  z[686] = par__Rfoot_M*(z[36]*z[652]+z[37]*z[653]+z[45]*z[118]+z[48]*z[119]) + 
  par__Rtibia_M*(z[20]*z[651]+z[21]*z[27]+z[21]*z[650]+z[29]*z[112]);
  z[687] = par__Rfoot_M*(z[34]*z[652]+z[35]*z[653]+z[45]*z[114]+z[48]*z[116]) + 
  par__Rtibia_M*(z[18]*z[651]+z[19]*z[27]+z[19]*z[650]+z[29]*z[111]);
  z[688] = -z[25]*z[33] - z[31]*z[32];
  z[689] = z[315]*z[34]*z[42] + z[315]*z[35]*z[45] + z[315]*z[35]*z[688] + 
  z[315]*z[48]*z[117] + z[316]*z[36]*z[42] + z[316]*z[37]*z[45] + z[316]*
  z[37]*z[688] + z[316]*z[48]*z[120] + z[34]*z[42]*z[339] + z[35]*z[334]*
  z[339] + z[35]*z[339]*z[688] + z[36]*z[42]*z[340] + z[37]*z[334]*z[340] + 
  z[37]*z[340]*z[688] + z[117]*z[337]*z[339] + z[120]*z[337]*z[340];
  z[690] = par__Rankle_y*z[32]*q5p*q6p + q1p*z[490] + q2p*z[487] + q3p*z[488] + 
  q4p*z[489];
  z[691] = -z[24]*z[33] - z[30]*z[32];
  z[692] = z[117]*q2p + z[120]*q3p + z[688]*q1p + z[691]*q4p + par__Rankle_y*
  z[32]*q5p;
  z[693] = z[31]*z[33]*q6p - z[25]*z[32]*q6p - z[32]*z[483] - z[33]*z[482];
  z[694] = z[18]*z[33]*q6p - z[19]*z[32]*q6p - z[32]*z[471] - z[33]*z[472];
  z[695] = z[20]*z[33]*q6p - z[21]*z[32]*q6p - z[32]*z[474] - z[33]*z[475];
  z[696] = z[30]*z[33]*q6p + z[32]*z[478] + z[33]*z[477] - z[24]*z[32]*q6p;
  z[697] = q1p*z[693] + q2p*z[694] + q3p*z[695] + q4p*z[696] - par__Rankle_y*
  z[33]*q5p*q6p;
  z[698] = z[35]*q2p + z[37]*q3p + z[41]*q4p + z[42]*q1p + par__Rankle_y*
  z[33]*q5p;
  z[699] = par__Rfoot_M*(z[42]*z[493]+z[688]*z[494]+z[48]*(z[697]+z[406]*
  z[698])+z[45]*(z[690]-z[406]*z[692]));
  z[700] = par__Rfoot_M*(z[42]*z[45]+z[48]*z[688]);
  z[701] = par__Rfoot_M*(z[41]*z[45]+z[42]*z[44]+z[47]*z[688]+z[48]*z[691]);
  z[702] = par__Rfoot_M*(z[42]*z[46]+z[49]*z[688]+par__Rankle_y*z[32]*z[48]+
  par__Rankle_y*z[33]*z[45]);
  z[703] = par__Rfoot_M*(par__Rfoot_CMx*z[688]-par__Rfoot_CMy*z[42]);
  z[704] = par__Rfoot_M*(z[36]*z[42]+z[37]*z[45]+z[37]*z[688]+z[48]*z[120]);
  z[705] = par__Rfoot_M*(z[34]*z[42]+z[35]*z[45]+z[35]*z[688]+z[48]*z[117]);
  z[706] = par__Lhip_y*z[50] - par__Lhip_x*z[51];
  z[707] = z[63]*z[706] - z[56]*z[64];
  z[708] = z[56]*z[63] + z[64]*z[706];
  z[709] = z[79]*z[708] + z[80]*z[707];
  z[710] = z[79]*z[707] - z[80]*z[708];
  z[711] = par__Lankle_y*z[65]*z[357] + par__Lankle_y*z[124]*z[356] + z[305]*
  z[52]*z[56] + z[305]*z[53]*z[57] + z[305]*z[53]*z[706] + z[305]*z[58]*
  z[121] + z[306]*z[52]*z[57] + z[306]*z[52]*z[706] + z[306]*z[53]*z[58] + 
  z[306]*z[54]*z[56] + z[307]*z[81]*z[709] + z[307]*z[82]*z[710] + z[307]*
  z[92]*z[128] + z[307]*z[95]*z[130] + z[308]*z[81]*z[92] + z[308]*z[82]*
  z[95] + z[308]*z[83]*z[709] + z[308]*z[84]*z[710] + z[309]*z[65]*z[708] + 
  z[309]*z[66]*z[707] + z[309]*z[74]*z[124] + z[309]*z[76]*z[125] + z[310]*
  z[65]*z[74] + z[310]*z[66]*z[76] + z[310]*z[67]*z[708] + z[310]*z[68]*
  z[707] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] + z[81]*z[351]*z[357] + 
  z[81]*z[356]*z[709] + z[82]*z[354]*z[357] + z[82]*z[356]*z[710] + z[83]*
  z[357]*z[709] + z[84]*z[357]*z[710] + z[128]*z[351]*z[356] + z[130]*z[354]*
  z[356] - z[53]*z[59]*z[357] - z[59]*z[121]*z[356];
  z[712] = z[53]*q3p + z[121]*q2p + z[706]*q1p;
  z[713] = q2p*z[581] + q3p*z[415] - q1p*z[413];
  z[714] = z[52]*q3p + z[53]*q2p + z[56]*q1p;
  z[715] = z[63]*z[706]*q8p - z[56]*z[64]*q8p - z[63]*z[416] - z[64]*z[413];
  z[716] = z[64]*z[416] - z[56]*z[63]*q8p - z[64]*z[706]*q8p - z[63]*z[413];
  z[717] = z[79]*z[707]*q9p + z[79]*z[715] + z[80]*z[716] - z[80]*z[708]*q9p;
  z[718] = q1p*z[717] + q2p*z[586] + q3p*z[424];
  z[719] = z[82]*q3p + z[130]*q2p + z[710]*q1p;
  z[720] = z[79]*z[716] - z[79]*z[708]*q9p - z[80]*z[707]*q9p - z[80]*z[715];
  z[721] = q1p*z[720] + q2p*z[589] + q3p*z[438];
  z[722] = z[81]*q3p + z[128]*q2p + z[709]*q1p;
  z[723] = q1p*z[715] + q2p*z[584] + q3p*z[422];
  z[724] = z[66]*q3p + z[125]*q2p + z[707]*q1p;
  z[725] = q1p*z[716] + q2p*z[585] + q3p*z[423];
  z[726] = z[65]*q3p + z[124]*q2p + z[708]*q1p;
  z[727] = par__Lfemur_M*(z[56]*z[418]+z[706]*z[419]+z[58]*(z[713]+z[402]*
  z[714])+z[57]*(z[417]-z[402]*z[712])) + par__Lfoot_M*(z[709]*z[444]+z[710]*
  z[445]+z[95]*(z[721]+z[403]*z[722])+z[92]*(z[718]-z[403]*z[719])) + par__Ltibia_M*(
  z[707]*z[451]+z[708]*z[450]+z[76]*(z[725]+z[404]*z[726])+z[74]*(z[723]-
  z[404]*z[724]));
  z[728] = 2*par__Lfemur_M*(z[56]*z[57]+z[58]*z[706]) + 2*par__Lfoot_M*(z[92]*
  z[709]+z[95]*z[710]) + 2*par__Ltibia_M*(z[74]*z[708]+z[76]*z[707]);
  z[729] = par__Lfoot_M*(z[91]*z[709]+z[94]*z[710]) + par__Ltibia_M*(z[73]*
  z[708]+z[75]*z[707]) + par__Lfemur_M*(par__Lfemur_CMx*z[706]-par__Lfemur_CMy*
  z[56]);
  z[730] = par__Lfoot_M*(z[93]*z[709]+z[96]*z[710]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[707]-par__Ltibia_CMy*z[708]);
  z[731] = par__Lfoot_M*(par__Lfoot_CMx*z[710]-par__Lfoot_CMy*z[709]);
  z[732] = par__Lfemur_M*(z[52]*z[57]+z[52]*z[706]+z[53]*z[58]+z[54]*z[56]) + 
  par__Lfoot_M*(z[81]*z[92]+z[82]*z[95]+z[83]*z[709]+z[84]*z[710]) + par__Ltibia_M*(
  z[65]*z[74]+z[66]*z[76]+z[67]*z[708]+z[68]*z[707]);
  z[733] = par__Lfemur_M*(z[52]*z[56]+z[53]*z[57]+z[53]*z[706]+z[58]*z[121]) + 
  par__Lfoot_M*(z[81]*z[709]+z[82]*z[710]+z[92]*z[128]+z[95]*z[130]) + par__Ltibia_M*(
  z[65]*z[708]+z[66]*z[707]+z[74]*z[124]+z[76]*z[125]);
  z[734] = z[63]*z[122] + z[64]*z[123] - z[61]*z[63] - z[62]*z[64];
  z[735] = z[62]*z[63] + z[64]*z[122] - z[61]*z[64] - z[63]*z[123];
  z[736] = z[79]*z[735] + z[80]*z[734];
  z[737] = z[79]*z[734] - z[80]*z[735];
  z[738] = par__Lankle_y*z[66]*z[356] + par__Lankle_y*z[68]*z[357] + z[307]*
  z[81]*z[736] + z[307]*z[82]*z[737] + z[307]*z[92]*z[129] + z[307]*z[95]*
  z[131] + z[308]*z[83]*z[736] + z[308]*z[84]*z[737] + z[308]*z[92]*z[133] + 
  z[308]*z[95]*z[134] + z[309]*z[65]*z[735] + z[309]*z[66]*z[74] + z[309]*
  z[66]*z[734] + z[309]*z[76]*z[126] + z[310]*z[67]*z[735] + z[310]*z[68]*
  z[74] + z[310]*z[68]*z[734] + z[310]*z[76]*z[127] + z[52]*z[123]*z[356] + 
  z[54]*z[123]*z[357] + z[81]*z[356]*z[736] + z[82]*z[356]*z[737] + z[83]*
  z[357]*z[736] + z[84]*z[357]*z[737] + z[129]*z[351]*z[356] + z[131]*z[354]*
  z[356] + z[133]*z[351]*z[357] + z[134]*z[354]*z[357] - z[52]*z[122]*z[357] - 
  z[53]*z[122]*z[356];
  z[739] = z[63]*z[122]*q8p + z[64]*z[123]*q8p + z[64]*z[258]*q8p + z[63]*
  z[431] - z[61]*z[63]*q8p - z[62]*z[64]*q8p - z[63]*z[259]*q8p - z[64]*
  z[432];
  z[740] = z[61]*z[64]*q8p + z[63]*z[123]*q8p + z[63]*z[258]*q8p + z[64]*
  z[259]*q8p - z[62]*z[63]*q8p - z[64]*z[122]*q8p - z[63]*z[432] - z[64]*
  z[431];
  z[741] = z[79]*z[734]*q9p + z[79]*z[739] + z[80]*z[740] - z[80]*z[735]*q9p;
  z[742] = z[52]*z[64]*q8p - z[53]*z[63]*q8p - z[63]*z[411] - z[64]*z[415];
  z[743] = z[79]*z[126]*q9p + z[79]*z[423] + z[80]*z[742] - z[66]*z[80]*q9p;
  z[744] = z[54]*z[64]*q8p - z[52]*z[63]*q8p - z[63]*z[412] - z[64]*z[411];
  z[745] = z[79]*z[127]*q9p + z[79]*z[426] + z[80]*z[744] - z[68]*z[80]*q9p;
  z[746] = z[60]*z[63] + z[63]*z[122] + z[64]*z[123] - z[59]*z[64];
  z[747] = z[59]*z[63] + z[60]*z[64] + z[64]*z[122] - z[63]*z[123];
  z[748] = (z[59]*z[64]+z[63]*z[259]-2*z[63]*z[122]-2*z[64]*z[123]-z[60]*
  z[63]-z[64]*z[258])*q8p;
  z[749] = (z[59]*z[63]+z[60]*z[64]+2*z[64]*z[122]-2*z[63]*z[123]-z[63]*
  z[258]-z[64]*z[259])*q8p;
  z[750] = z[79]*z[746]*q9p - z[80]*z[747]*q9p - z[79]*z[748] - z[80]*z[749];
  z[751] = q1p*z[741] + q2p*z[743] + q3p*z[745] + q7p*z[750];
  z[752] = z[79]*z[746] - z[80]*z[747];
  z[753] = z[131]*q2p + z[134]*q3p + z[737]*q1p + z[752]*q7p;
  z[754] = z[79]*z[740] - z[79]*z[735]*q9p - z[80]*z[734]*q9p - z[80]*z[739];
  z[755] = z[79]*z[742] - z[66]*z[79]*q9p - z[80]*z[126]*q9p - z[80]*z[423];
  z[756] = z[79]*z[744] - z[68]*z[79]*q9p - z[80]*z[127]*q9p - z[80]*z[426];
  z[757] = z[80]*z[748] - z[79]*z[747]*q9p - z[80]*z[746]*q9p - z[79]*z[749];
  z[758] = q1p*z[754] + q2p*z[755] + q3p*z[756] + q7p*z[757];
  z[759] = z[79]*z[747] + z[80]*z[746];
  z[760] = z[129]*q2p + z[133]*q3p + z[736]*q1p + z[759]*q7p;
  z[761] = q1p*z[739] + q2p*z[423] + q3p*z[426] - q7p*z[748];
  z[762] = z[126]*q2p + z[127]*q3p + z[734]*q1p + z[746]*q7p;
  z[763] = q1p*z[740] + q2p*z[742] + q3p*z[744] - q7p*z[749];
  z[764] = z[66]*q2p + z[68]*q3p + z[735]*q1p + z[747]*q7p;
  z[765] = par__Lfoot_M*(z[736]*z[444]+z[737]*z[445]+z[95]*(z[758]+z[403]*
  z[760])+z[92]*(z[751]-z[403]*z[753])) + par__Ltibia_M*(z[734]*z[451]+z[735]*
  z[450]+z[76]*(z[763]+z[404]*z[764])+z[74]*(z[761]-z[404]*z[762]));
  z[766] = 2*par__Lfoot_M*(z[92]*z[736]+z[95]*z[737]) + 2*par__Ltibia_M*(
  z[74]*z[735]+z[76]*z[734]);
  z[767] = par__Lfoot_M*(z[91]*z[736]+z[92]*z[759]+z[94]*z[737]+z[95]*z[752]) + 
  par__Ltibia_M*(z[73]*z[735]+z[74]*z[747]+z[75]*z[734]+z[76]*z[746]);
  z[768] = par__Lfoot_M*(z[93]*z[736]+z[96]*z[737]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[734]-par__Ltibia_CMy*z[735]);
  z[769] = par__Lfoot_M*(par__Lfoot_CMx*z[737]-par__Lfoot_CMy*z[736]);
  z[770] = par__Lfoot_M*(z[83]*z[736]+z[84]*z[737]+z[92]*z[133]+z[95]*z[134]) + 
  par__Ltibia_M*(z[67]*z[735]+z[68]*z[74]+z[68]*z[734]+z[76]*z[127]);
  z[771] = par__Lfoot_M*(z[81]*z[736]+z[82]*z[737]+z[92]*z[129]+z[95]*z[131]) + 
  par__Ltibia_M*(z[65]*z[735]+z[66]*z[74]+z[66]*z[734]+z[76]*z[126]);
  z[772] = -z[72]*z[80] - z[78]*z[79];
  z[773] = z[307]*z[81]*z[89] + z[307]*z[82]*z[92] + z[307]*z[82]*z[772] + 
  z[307]*z[95]*z[132] + z[308]*z[83]*z[89] + z[308]*z[84]*z[92] + z[308]*
  z[84]*z[772] + z[308]*z[95]*z[135] + z[81]*z[89]*z[356] + z[82]*z[351]*
  z[356] + z[82]*z[356]*z[772] + z[83]*z[89]*z[357] + z[84]*z[351]*z[357] + 
  z[84]*z[357]*z[772] + z[132]*z[354]*z[356] + z[135]*z[354]*z[357];
  z[774] = par__Lankle_y*z[79]*q8p*q9p + q1p*z[441] + q2p*z[438] + q3p*z[439] + 
  q7p*z[440];
  z[775] = -z[71]*z[80] - z[77]*z[79];
  z[776] = z[132]*q2p + z[135]*q3p + z[772]*q1p + z[775]*q7p + par__Lankle_y*
  z[79]*q8p;
  z[777] = z[78]*z[80]*q9p - z[72]*z[79]*q9p - z[79]*z[434] - z[80]*z[433];
  z[778] = z[65]*z[80]*q9p - z[66]*z[79]*q9p - z[79]*z[422] - z[80]*z[423];
  z[779] = z[67]*z[80]*q9p - z[68]*z[79]*q9p - z[79]*z[425] - z[80]*z[426];
  z[780] = z[77]*z[80]*q9p + z[79]*z[429] + z[80]*z[428] - z[71]*z[79]*q9p;
  z[781] = q1p*z[777] + q2p*z[778] + q3p*z[779] + q7p*z[780] - par__Lankle_y*
  z[80]*q8p*q9p;
  z[782] = z[82]*q2p + z[84]*q3p + z[88]*q7p + z[89]*q1p + par__Lankle_y*
  z[80]*q8p;
  z[783] = par__Lfoot_M*(z[89]*z[444]+z[772]*z[445]+z[95]*(z[781]+z[403]*
  z[782])+z[92]*(z[774]-z[403]*z[776]));
  z[784] = par__Lfoot_M*(z[89]*z[92]+z[95]*z[772]);
  z[785] = par__Lfoot_M*(z[88]*z[92]+z[89]*z[91]+z[94]*z[772]+z[95]*z[775]);
  z[786] = par__Lfoot_M*(z[89]*z[93]+z[96]*z[772]+par__Lankle_y*z[79]*z[95]+
  par__Lankle_y*z[80]*z[92]);
  z[787] = par__Lfoot_M*(par__Lfoot_CMx*z[772]-par__Lfoot_CMy*z[89]);
  z[788] = par__Lfoot_M*(z[83]*z[89]+z[84]*z[92]+z[84]*z[772]+z[95]*z[135]);
  z[789] = par__Lfoot_M*(z[81]*z[89]+z[82]*z[92]+z[82]*z[772]+z[95]*z[132]);
  z[790] = par__back_y*z[97] - par__back_x*z[98];
  z[791] = z[319]*z[99]*z[103] + z[319]*z[100]*z[104] + z[319]*z[100]*z[790] + 
  z[319]*z[105]*z[136] + z[320]*z[99]*z[104] + z[320]*z[99]*z[790] + z[320]*
  z[100]*z[105] + z[320]*z[101]*z[103];
  z[792] = z[100]*q3p + z[136]*q2p + z[790]*q1p;
  z[793] = q2p*z[615] + q3p*z[507] - q1p*z[505];
  z[794] = z[99]*q3p + z[100]*q2p + z[103]*q1p;
  z[795] = par__torso_M*(z[103]*z[510]+z[790]*z[511]+z[105]*(z[793]+z[408]*
  z[794])+z[104]*(z[509]-z[408]*z[792]));
  z[796] = par__torso_M*(z[103]*z[104]+z[105]*z[790]);
  z[797] = par__torso_M*(par__torso_CMx*z[790]-par__torso_CMy*z[103]);
  z[798] = par__torso_M*(z[99]*z[104]+z[99]*z[790]+z[100]*z[105]+z[101]*
  z[103]);
  z[799] = par__torso_M*(z[99]*z[103]+z[100]*z[104]+z[100]*z[790]+z[105]*
  z[136]);
  z[800] = z[306]*z[52]*z[121] + z[306]*pow(z[52],2) + z[306]*z[53]*z[54] + 
  z[306]*pow(z[53],2) + z[308]*pow(z[81],2) + z[308]*pow(z[82],2) + z[308]*
  z[83]*z[128] + z[308]*z[84]*z[130] + z[310]*pow(z[65],2) + z[310]*pow(z[66],
  2) + z[310]*z[67]*z[124] + z[310]*z[68]*z[125] + z[314]*z[5]*z[106] + 
  z[314]*pow(z[5],2) + z[314]*z[6]*z[7] + z[314]*pow(z[6],2) + z[316]*pow(
  z[34],2) + z[316]*pow(z[35],2) + z[316]*z[36]*z[113] + z[316]*z[37]*z[115] + 
  z[318]*pow(z[18],2) + z[318]*pow(z[19],2) + z[318]*z[20]*z[109] + z[318]*
  z[21]*z[110] + z[320]*z[99]*z[136] + z[320]*pow(z[99],2) + z[320]*z[100]*
  z[101] + z[320]*pow(z[100],2) + z[36]*z[113]*z[340] + z[37]*z[115]*z[340] + 
  z[83]*z[128]*z[357] + z[84]*z[130]*z[357] + z[340]*pow(z[34],2) + z[340]*
  pow(z[35],2) + z[357]*pow(z[81],2) + z[357]*pow(z[82],2) + 2*z[305]*z[52]*
  z[53] + 2*z[305]*z[53]*z[121] + 2*z[307]*z[81]*z[128] + 2*z[307]*z[82]*
  z[130] + 2*z[309]*z[65]*z[124] + 2*z[309]*z[66]*z[125] + 2*z[313]*z[5]*z[6] + 
  2*z[313]*z[6]*z[106] + 2*z[315]*z[34]*z[113] + 2*z[315]*z[35]*z[115] + 2*
  z[317]*z[18]*z[109] + 2*z[317]*z[19]*z[110] + 2*z[319]*z[99]*z[100] + 2*
  z[319]*z[100]*z[136] + 2*z[34]*z[113]*z[339] + 2*z[35]*z[115]*z[339] + 2*
  z[81]*z[128]*z[356] + 2*z[82]*z[130]*z[356] + z[1]*z[324] - z[1]*z[575] - 
  z[2]*z[323] - z[2]*z[574] - z[321]*pow(z[1],2)*z[573] - z[321]*pow(z[2],2)*
  z[573];
  z[801] = par__Lfemur_M*(z[53]*z[418]+z[121]*z[419]+z[53]*(z[582]+z[402]*
  z[583])+z[52]*(z[579]-z[402]*z[580])) + par__Lfoot_M*(z[128]*z[444]+z[130]*
  z[445]+z[82]*(z[590]+z[403]*z[591])+z[81]*(z[587]-z[403]*z[588])) + par__Ltibia_M*(
  z[124]*z[450]+z[125]*z[451]+z[66]*(z[594]+z[404]*z[595])+z[65]*(z[592]-
  z[404]*z[593])) + par__Rfemur_M*(z[6]*z[467]+z[106]*z[468]+z[6]*(z[599]+
  z[405]*z[600])+z[5]*(z[596]-z[405]*z[597])) + par__Rfoot_M*(z[113]*z[493]+
  z[115]*z[494]+z[35]*(z[607]+z[406]*z[608])+z[34]*(z[604]-z[406]*z[605])) + 
  par__Rtibia_M*(z[109]*z[499]+z[110]*z[500]+z[19]*(z[611]+z[407]*z[612])+
  z[18]*(z[609]-z[407]*z[610])) + par__torso_M*(z[100]*z[510]+z[136]*z[511]+
  z[100]*(z[616]+z[408]*z[617])+z[99]*(z[613]-z[408]*z[614])) - par__pelvis_M*(
  z[1]*z[457]+z[2]*z[456]+z[1]*(z[455]+q1p*z[619])-z[2]*(z[454]-q1p*z[618]));
  z[802] = par__Lfoot_M*(z[91]*z[128]+z[94]*z[130]) + par__Ltibia_M*(z[73]*
  z[124]+z[75]*z[125]) + par__Lfemur_M*(par__Lfemur_CMx*z[121]-par__Lfemur_CMy*
  z[53]);
  z[803] = par__Lfemur_M*z[52]*(z[52]+z[121]) + par__Lfemur_M*z[53]*(z[53]+
  z[54]) + par__Rfemur_M*z[5]*(z[5]+z[106]) + par__Rfemur_M*z[6]*(z[6]+z[7]) + 
  par__torso_M*z[99]*(z[99]+z[136]) + par__torso_M*z[100]*(z[100]+z[101]) + 
  par__Lfoot_M*(pow(z[81],2)+pow(z[82],2)+z[83]*z[128]+z[84]*z[130]) + par__Ltibia_M*(
  pow(z[65],2)+pow(z[66],2)+z[67]*z[124]+z[68]*z[125]) + par__Rfoot_M*(pow(
  z[34],2)+pow(z[35],2)+z[36]*z[113]+z[37]*z[115]) + par__Rtibia_M*(pow(z[18],
  2)+pow(z[19],2)+z[20]*z[109]+z[21]*z[110]);
  z[804] = 2*par__Lfemur_M*z[53]*(z[52]+z[121]) + 2*par__Rfemur_M*z[6]*(z[5]+
  z[106]) + 2*par__torso_M*z[100]*(z[99]+z[136]) + 2*par__Lfoot_M*(z[81]*
  z[128]+z[82]*z[130]) + 2*par__Ltibia_M*(z[65]*z[124]+z[66]*z[125]) + 2*
  par__Rfoot_M*(z[34]*z[113]+z[35]*z[115]) + 2*par__Rtibia_M*(z[18]*z[109]+
  z[19]*z[110]);
  z[805] = par__Lfoot_M*(par__Lfoot_CMx*z[130]-par__Lfoot_CMy*z[128]);
  z[806] = par__Lfoot_M*(z[93]*z[128]+z[96]*z[130]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[125]-par__Ltibia_CMy*z[124]);
  z[807] = par__Rfoot_M*(z[44]*z[113]+z[47]*z[115]) + par__Rtibia_M*(z[26]*
  z[109]+z[28]*z[110]) + par__Rfemur_M*(par__Rfemur_CMx*z[106]-par__Rfemur_CMy*
  z[6]);
  z[808] = par__Rfoot_M*(par__Rfoot_CMx*z[115]-par__Rfoot_CMy*z[113]);
  z[809] = par__Rfoot_M*(z[46]*z[113]+z[49]*z[115]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[110]-par__Rtibia_CMy*z[109]);
  z[810] = par__torso_M*(par__torso_CMx*z[136]-par__torso_CMy*z[100]);
  z[811] = z[314]*z[5]*z[106] + z[314]*pow(z[5],2) + z[314]*z[6]*z[7] + 
  z[314]*pow(z[6],2) + z[316]*pow(z[34],2) + z[316]*pow(z[35],2) + z[316]*
  z[36]*z[113] + z[316]*z[37]*z[115] + z[318]*pow(z[18],2) + z[318]*pow(z[19],
  2) + z[318]*z[20]*z[109] + z[318]*z[21]*z[110] + z[36]*z[113]*z[340] + 
  z[37]*z[115]*z[340] + z[340]*pow(z[34],2) + z[340]*pow(z[35],2) + 2*z[313]*
  z[5]*z[6] + 2*z[313]*z[6]*z[106] + 2*z[315]*z[34]*z[113] + 2*z[315]*z[35]*
  z[115] + 2*z[317]*z[18]*z[109] + 2*z[317]*z[19]*z[110] + 2*z[34]*z[113]*
  z[339] + 2*z[35]*z[115]*z[339];
  z[812] = par__Rfemur_M*(z[6]*z[467]+z[106]*z[468]+z[6]*(z[629]+z[405]*
  z[630])+z[5]*(z[466]-z[405]*z[628])) + par__Rfoot_M*(z[113]*z[493]+z[115]*
  z[494]+z[35]*(z[637]+z[406]*z[638])+z[34]*(z[634]-z[406]*z[635])) + par__Rtibia_M*(
  z[109]*z[499]+z[110]*z[500]+z[19]*(z[641]+z[407]*z[642])+z[18]*(z[639]-
  z[407]*z[640]));
  z[813] = par__Rfemur_M*z[5]*(z[5]+z[106]) + par__Rfemur_M*z[6]*(z[6]+z[7]) + 
  par__Rfoot_M*(pow(z[34],2)+pow(z[35],2)+z[36]*z[113]+z[37]*z[115]) + par__Rtibia_M*(
  pow(z[18],2)+pow(z[19],2)+z[20]*z[109]+z[21]*z[110]);
  z[814] = 2*par__Rfemur_M*z[6]*(z[5]+z[106]) + 2*par__Rfoot_M*(z[34]*z[113]+
  z[35]*z[115]) + 2*par__Rtibia_M*(z[18]*z[109]+z[19]*z[110]);
  z[815] = z[316]*z[34]*z[118] + z[316]*z[35]*z[119] + z[316]*z[36]*z[114] + 
  z[316]*z[37]*z[116] + z[318]*z[18]*z[21] + z[318]*z[19]*z[20] + z[318]*
  z[19]*z[112] + z[318]*z[21]*z[111] + z[34]*z[118]*z[340] + z[35]*z[119]*
  z[340] + z[36]*z[114]*z[340] + z[37]*z[116]*z[340] + 2*z[315]*z[34]*z[114] + 
  2*z[315]*z[35]*z[116] + 2*z[317]*z[18]*z[19] + 2*z[317]*z[19]*z[111] + 2*
  z[34]*z[114]*z[339] + 2*z[35]*z[116]*z[339];
  z[816] = par__Rfoot_M*(z[114]*z[493]+z[116]*z[494]+z[35]*(z[674]+z[406]*
  z[676])+z[34]*(z[667]-z[406]*z[669])) + par__Rtibia_M*(z[19]*z[499]+z[111]*
  z[500]+z[19]*(z[679]+z[407]*z[680])+z[18]*(z[677]-z[407]*z[678]));
  z[817] = par__Rfoot_M*(z[34]*z[118]+z[35]*z[119]+z[36]*z[114]+z[37]*z[116]) + 
  par__Rtibia_M*(z[18]*z[21]+z[19]*z[20]+z[19]*z[112]+z[21]*z[111]);
  z[818] = 2*par__Rtibia_M*z[19]*(z[18]+z[111]) + 2*par__Rfoot_M*(z[34]*
  z[114]+z[35]*z[116]);
  z[819] = par__Rfoot_M*(z[34]*z[675]+z[35]*z[668]+z[44]*z[114]+z[47]*z[116]) + 
  par__Rtibia_M*(z[18]*z[663]+z[19]*z[26]+z[19]*z[662]+z[28]*z[111]);
  z[820] = par__Rfoot_M*(par__Rfoot_CMx*z[116]-par__Rfoot_CMy*z[114]);
  z[821] = par__Rfoot_M*(z[46]*z[114]+z[49]*z[116]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[111]-par__Rtibia_CMy*z[19]);
  z[822] = z[316]*z[34]*z[37] + z[316]*z[35]*z[36] + z[316]*z[35]*z[120] + 
  z[316]*z[37]*z[117] + z[34]*z[37]*z[340] + z[35]*z[36]*z[340] + z[35]*
  z[120]*z[340] + z[37]*z[117]*z[340] + 2*z[315]*z[34]*z[35] + 2*z[315]*z[35]*
  z[117] + 2*z[34]*z[35]*z[339] + 2*z[35]*z[117]*z[339];
  z[823] = par__Rfoot_M*(z[35]*z[493]+z[117]*z[494]+z[35]*(z[697]+z[406]*
  z[698])+z[34]*(z[690]-z[406]*z[692]));
  z[824] = par__Rfoot_M*(z[34]*z[37]+z[35]*z[36]+z[35]*z[120]+z[37]*z[117]);
  z[825] = par__Rfoot_M*z[35]*(z[34]+z[117]);
  z[826] = par__Rfoot_M*(z[34]*z[41]+z[35]*z[44]+z[35]*z[691]+z[47]*z[117]);
  z[827] = par__Rfoot_M*(par__Rfoot_CMx*z[117]-par__Rfoot_CMy*z[35]);
  z[828] = par__Rfoot_M*(z[35]*z[46]+z[49]*z[117]+par__Rankle_y*z[32]*z[35]+
  par__Rankle_y*z[33]*z[34]);
  z[829] = z[306]*z[52]*z[121] + z[306]*pow(z[52],2) + z[306]*z[53]*z[54] + 
  z[306]*pow(z[53],2) + z[308]*pow(z[81],2) + z[308]*pow(z[82],2) + z[308]*
  z[83]*z[128] + z[308]*z[84]*z[130] + z[310]*pow(z[65],2) + z[310]*pow(z[66],
  2) + z[310]*z[67]*z[124] + z[310]*z[68]*z[125] + z[83]*z[128]*z[357] + 
  z[84]*z[130]*z[357] + z[357]*pow(z[81],2) + z[357]*pow(z[82],2) + 2*z[305]*
  z[52]*z[53] + 2*z[305]*z[53]*z[121] + 2*z[307]*z[81]*z[128] + 2*z[307]*
  z[82]*z[130] + 2*z[309]*z[65]*z[124] + 2*z[309]*z[66]*z[125] + 2*z[81]*
  z[128]*z[356] + 2*z[82]*z[130]*z[356];
  z[830] = par__Lfemur_M*(z[53]*z[418]+z[121]*z[419]+z[53]*(z[713]+z[402]*
  z[714])+z[52]*(z[417]-z[402]*z[712])) + par__Lfoot_M*(z[128]*z[444]+z[130]*
  z[445]+z[82]*(z[721]+z[403]*z[722])+z[81]*(z[718]-z[403]*z[719])) + par__Ltibia_M*(
  z[124]*z[450]+z[125]*z[451]+z[66]*(z[725]+z[404]*z[726])+z[65]*(z[723]-
  z[404]*z[724]));
  z[831] = par__Lfemur_M*z[52]*(z[52]+z[121]) + par__Lfemur_M*z[53]*(z[53]+
  z[54]) + par__Lfoot_M*(pow(z[81],2)+pow(z[82],2)+z[83]*z[128]+z[84]*z[130]) + 
  par__Ltibia_M*(pow(z[65],2)+pow(z[66],2)+z[67]*z[124]+z[68]*z[125]);
  z[832] = 2*par__Lfemur_M*z[53]*(z[52]+z[121]) + 2*par__Lfoot_M*(z[81]*
  z[128]+z[82]*z[130]) + 2*par__Ltibia_M*(z[65]*z[124]+z[66]*z[125]);
  z[833] = z[308]*z[81]*z[133] + z[308]*z[82]*z[134] + z[308]*z[83]*z[129] + 
  z[308]*z[84]*z[131] + z[310]*z[65]*z[68] + z[310]*z[66]*z[67] + z[310]*
  z[66]*z[127] + z[310]*z[68]*z[126] + z[81]*z[133]*z[357] + z[82]*z[134]*
  z[357] + z[83]*z[129]*z[357] + z[84]*z[131]*z[357] + 2*z[307]*z[81]*z[129] + 
  2*z[307]*z[82]*z[131] + 2*z[309]*z[65]*z[66] + 2*z[309]*z[66]*z[126] + 2*
  z[81]*z[129]*z[356] + 2*z[82]*z[131]*z[356];
  z[834] = par__Lfoot_M*(z[129]*z[444]+z[131]*z[445]+z[82]*(z[758]+z[403]*
  z[760])+z[81]*(z[751]-z[403]*z[753])) + par__Ltibia_M*(z[66]*z[450]+z[126]*
  z[451]+z[66]*(z[763]+z[404]*z[764])+z[65]*(z[761]-z[404]*z[762]));
  z[835] = par__Lfoot_M*(z[81]*z[759]+z[82]*z[752]+z[91]*z[129]+z[94]*z[131]) + 
  par__Ltibia_M*(z[65]*z[747]+z[66]*z[73]+z[66]*z[746]+z[75]*z[126]);
  z[836] = par__Lfoot_M*(z[81]*z[133]+z[82]*z[134]+z[83]*z[129]+z[84]*z[131]) + 
  par__Ltibia_M*(z[65]*z[68]+z[66]*z[67]+z[66]*z[127]+z[68]*z[126]);
  z[837] = 2*par__Ltibia_M*z[66]*(z[65]+z[126]) + 2*par__Lfoot_M*(z[81]*
  z[129]+z[82]*z[131]);
  z[838] = par__Lfoot_M*(par__Lfoot_CMx*z[131]-par__Lfoot_CMy*z[129]);
  z[839] = par__Lfoot_M*(z[93]*z[129]+z[96]*z[131]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[126]-par__Ltibia_CMy*z[66]);
  z[840] = z[308]*z[81]*z[84] + z[308]*z[82]*z[83] + z[308]*z[82]*z[135] + 
  z[308]*z[84]*z[132] + z[81]*z[84]*z[357] + z[82]*z[83]*z[357] + z[82]*
  z[135]*z[357] + z[84]*z[132]*z[357] + 2*z[307]*z[81]*z[82] + 2*z[307]*z[82]*
  z[132] + 2*z[81]*z[82]*z[356] + 2*z[82]*z[132]*z[356];
  z[841] = par__Lfoot_M*(z[82]*z[444]+z[132]*z[445]+z[82]*(z[781]+z[403]*
  z[782])+z[81]*(z[774]-z[403]*z[776]));
  z[842] = par__Lfoot_M*(z[81]*z[88]+z[82]*z[91]+z[82]*z[775]+z[94]*z[132]);
  z[843] = par__Lfoot_M*(z[81]*z[84]+z[82]*z[83]+z[82]*z[135]+z[84]*z[132]);
  z[844] = par__Lfoot_M*z[82]*(z[81]+z[132]);
  z[845] = par__Lfoot_M*(par__Lfoot_CMx*z[132]-par__Lfoot_CMy*z[82]);
  z[846] = par__Lfoot_M*(z[82]*z[93]+z[96]*z[132]+par__Lankle_y*z[79]*z[82]+
  par__Lankle_y*z[80]*z[81]);
  z[847] = z[320]*z[99]*z[136] + z[320]*pow(z[99],2) + z[320]*z[100]*z[101] + 
  z[320]*pow(z[100],2) + 2*z[319]*z[99]*z[100] + 2*z[319]*z[100]*z[136];
  z[848] = par__torso_M*(z[100]*z[510]+z[136]*z[511]+z[100]*(z[793]+z[408]*
  z[794])+z[99]*(z[509]-z[408]*z[792]));
  z[849] = par__torso_M*(z[99]*(z[99]+z[136])+z[100]*(z[100]+z[101]));
  z[850] = par__torso_M*z[100]*(z[99]+z[136]);
  z[851] = z[305]*z[52]*z[121] + z[305]*pow(z[52],2) + z[305]*z[53]*z[54] + 
  z[305]*pow(z[53],2) + z[307]*pow(z[81],2) + z[307]*pow(z[82],2) + z[307]*
  z[83]*z[128] + z[307]*z[84]*z[130] + z[309]*pow(z[65],2) + z[309]*pow(z[66],
  2) + z[309]*z[67]*z[124] + z[309]*z[68]*z[125] + z[313]*z[5]*z[106] + 
  z[313]*pow(z[5],2) + z[313]*z[6]*z[7] + z[313]*pow(z[6],2) + z[315]*pow(
  z[34],2) + z[315]*pow(z[35],2) + z[315]*z[36]*z[113] + z[315]*z[37]*z[115] + 
  z[317]*pow(z[18],2) + z[317]*pow(z[19],2) + z[317]*z[20]*z[109] + z[317]*
  z[21]*z[110] + z[319]*z[99]*z[136] + z[319]*pow(z[99],2) + z[319]*z[100]*
  z[101] + z[319]*pow(z[100],2) + z[36]*z[113]*z[339] + z[37]*z[115]*z[339] + 
  z[83]*z[128]*z[356] + z[84]*z[130]*z[356] + z[339]*pow(z[34],2) + z[339]*
  pow(z[35],2) + z[356]*pow(z[81],2) + z[356]*pow(z[82],2) + 2*z[306]*z[52]*
  z[53] + 2*z[306]*z[52]*z[54] + 2*z[308]*z[81]*z[83] + 2*z[308]*z[82]*z[84] + 
  2*z[310]*z[65]*z[67] + 2*z[310]*z[66]*z[68] + 2*z[314]*z[5]*z[6] + 2*z[314]*
  z[5]*z[7] + 2*z[316]*z[34]*z[36] + 2*z[316]*z[35]*z[37] + 2*z[318]*z[18]*
  z[20] + 2*z[318]*z[19]*z[21] + 2*z[320]*z[99]*z[100] + 2*z[320]*z[99]*
  z[101] + 2*z[34]*z[36]*z[340] + 2*z[35]*z[37]*z[340] + 2*z[81]*z[83]*z[357] + 
  2*z[82]*z[84]*z[357] + z[1]*z[323] + z[1]*z[574] + z[2]*z[324] - z[2]*
  z[575];
  z[852] = par__Lfemur_M*(z[52]*z[418]+z[53]*z[419]+z[52]*(z[582]+z[402]*
  z[583])+z[54]*(z[579]-z[402]*z[580])) + par__Lfoot_M*(z[81]*z[444]+z[82]*
  z[445]+z[84]*(z[590]+z[403]*z[591])+z[83]*(z[587]-z[403]*z[588])) + par__Ltibia_M*(
  z[65]*z[450]+z[66]*z[451]+z[68]*(z[594]+z[404]*z[595])+z[67]*(z[592]-z[404]*
  z[593])) + par__Rfemur_M*(z[5]*z[467]+z[6]*z[468]+z[5]*(z[599]+z[405]*
  z[600])+z[7]*(z[596]-z[405]*z[597])) + par__Rfoot_M*(z[34]*z[493]+z[35]*
  z[494]+z[37]*(z[607]+z[406]*z[608])+z[36]*(z[604]-z[406]*z[605])) + par__Rtibia_M*(
  z[18]*z[499]+z[19]*z[500]+z[21]*(z[611]+z[407]*z[612])+z[20]*(z[609]-z[407]*
  z[610])) + par__torso_M*(z[99]*z[510]+z[100]*z[511]+z[99]*(z[616]+z[408]*
  z[617])+z[101]*(z[613]-z[408]*z[614])) + par__pelvis_M*(z[1]*z[456]-z[2]*
  z[457]-z[2]*(z[455]+q1p*z[619])-z[1]*(z[454]-q1p*z[618]));
  z[853] = 2*par__Lfemur_M*z[52]*(z[53]+z[54]) + 2*par__Rfemur_M*z[5]*(z[6]+
  z[7]) + 2*par__torso_M*z[99]*(z[100]+z[101]) + 2*par__Lfoot_M*(z[81]*z[83]+
  z[82]*z[84]) + 2*par__Ltibia_M*(z[65]*z[67]+z[66]*z[68]) + 2*par__Rfoot_M*(
  z[34]*z[36]+z[35]*z[37]) + 2*par__Rtibia_M*(z[18]*z[20]+z[19]*z[21]);
  z[854] = z[313]*z[5]*z[106] + z[313]*pow(z[5],2) + z[313]*z[6]*z[7] + 
  z[313]*pow(z[6],2) + z[315]*pow(z[34],2) + z[315]*pow(z[35],2) + z[315]*
  z[36]*z[113] + z[315]*z[37]*z[115] + z[317]*pow(z[18],2) + z[317]*pow(z[19],
  2) + z[317]*z[20]*z[109] + z[317]*z[21]*z[110] + z[36]*z[113]*z[339] + 
  z[37]*z[115]*z[339] + z[339]*pow(z[34],2) + z[339]*pow(z[35],2) + 2*z[314]*
  z[5]*z[6] + 2*z[314]*z[5]*z[7] + 2*z[316]*z[34]*z[36] + 2*z[316]*z[35]*
  z[37] + 2*z[318]*z[18]*z[20] + 2*z[318]*z[19]*z[21] + 2*z[34]*z[36]*z[340] + 
  2*z[35]*z[37]*z[340];
  z[855] = par__Rfemur_M*(z[5]*z[467]+z[6]*z[468]+z[5]*(z[629]+z[405]*z[630])+
  z[7]*(z[466]-z[405]*z[628])) + par__Rfoot_M*(z[34]*z[493]+z[35]*z[494]+
  z[37]*(z[637]+z[406]*z[638])+z[36]*(z[634]-z[406]*z[635])) + par__Rtibia_M*(
  z[18]*z[499]+z[19]*z[500]+z[21]*(z[641]+z[407]*z[642])+z[20]*(z[639]-z[407]*
  z[640]));
  z[856] = 2*par__Rfemur_M*z[5]*(z[6]+z[7]) + 2*par__Rfoot_M*(z[34]*z[36]+
  z[35]*z[37]) + 2*par__Rtibia_M*(z[18]*z[20]+z[19]*z[21]);
  z[857] = z[315]*z[34]*z[118] + z[315]*z[35]*z[119] + z[315]*z[36]*z[114] + 
  z[315]*z[37]*z[116] + z[317]*z[18]*z[21] + z[317]*z[19]*z[20] + z[317]*
  z[19]*z[112] + z[317]*z[21]*z[111] + z[34]*z[118]*z[339] + z[35]*z[119]*
  z[339] + z[36]*z[114]*z[339] + z[37]*z[116]*z[339] + 2*z[316]*z[36]*z[118] + 
  2*z[316]*z[37]*z[119] + 2*z[318]*z[20]*z[21] + 2*z[318]*z[21]*z[112] + 2*
  z[36]*z[118]*z[340] + 2*z[37]*z[119]*z[340];
  z[858] = par__Rfoot_M*(z[118]*z[493]+z[119]*z[494]+z[37]*(z[674]+z[406]*
  z[676])+z[36]*(z[667]-z[406]*z[669])) + par__Rtibia_M*(z[21]*z[499]+z[112]*
  z[500]+z[21]*(z[679]+z[407]*z[680])+z[20]*(z[677]-z[407]*z[678]));
  z[859] = 2*par__Rtibia_M*z[21]*(z[20]+z[112]) + 2*par__Rfoot_M*(z[36]*
  z[118]+z[37]*z[119]);
  z[860] = par__Rfoot_M*(z[36]*z[675]+z[37]*z[668]+z[44]*z[118]+z[47]*z[119]) + 
  par__Rtibia_M*(z[20]*z[663]+z[21]*z[26]+z[21]*z[662]+z[28]*z[112]);
  z[861] = par__Rfoot_M*(par__Rfoot_CMx*z[119]-par__Rfoot_CMy*z[118]);
  z[862] = par__Rfoot_M*(z[46]*z[118]+z[49]*z[119]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[112]-par__Rtibia_CMy*z[21]);
  z[863] = z[315]*z[34]*z[37] + z[315]*z[35]*z[36] + z[315]*z[35]*z[120] + 
  z[315]*z[37]*z[117] + z[34]*z[37]*z[339] + z[35]*z[36]*z[339] + z[35]*
  z[120]*z[339] + z[37]*z[117]*z[339] + 2*z[316]*z[36]*z[37] + 2*z[316]*z[37]*
  z[120] + 2*z[36]*z[37]*z[340] + 2*z[37]*z[120]*z[340];
  z[864] = par__Rfoot_M*(z[37]*z[493]+z[120]*z[494]+z[37]*(z[697]+z[406]*
  z[698])+z[36]*(z[690]-z[406]*z[692]));
  z[865] = par__Rfoot_M*z[37]*(z[36]+z[120]);
  z[866] = par__Rfoot_M*(z[36]*z[41]+z[37]*z[44]+z[37]*z[691]+z[47]*z[120]);
  z[867] = par__Rfoot_M*(par__Rfoot_CMx*z[120]-par__Rfoot_CMy*z[37]);
  z[868] = par__Rfoot_M*(z[37]*z[46]+z[49]*z[120]+par__Rankle_y*z[32]*z[37]+
  par__Rankle_y*z[33]*z[36]);
  z[869] = z[305]*z[52]*z[121] + z[305]*pow(z[52],2) + z[305]*z[53]*z[54] + 
  z[305]*pow(z[53],2) + z[307]*pow(z[81],2) + z[307]*pow(z[82],2) + z[307]*
  z[83]*z[128] + z[307]*z[84]*z[130] + z[309]*pow(z[65],2) + z[309]*pow(z[66],
  2) + z[309]*z[67]*z[124] + z[309]*z[68]*z[125] + z[83]*z[128]*z[356] + 
  z[84]*z[130]*z[356] + z[356]*pow(z[81],2) + z[356]*pow(z[82],2) + 2*z[306]*
  z[52]*z[53] + 2*z[306]*z[52]*z[54] + 2*z[308]*z[81]*z[83] + 2*z[308]*z[82]*
  z[84] + 2*z[310]*z[65]*z[67] + 2*z[310]*z[66]*z[68] + 2*z[81]*z[83]*z[357] + 
  2*z[82]*z[84]*z[357];
  z[870] = par__Lfemur_M*(z[52]*z[418]+z[53]*z[419]+z[52]*(z[713]+z[402]*
  z[714])+z[54]*(z[417]-z[402]*z[712])) + par__Lfoot_M*(z[81]*z[444]+z[82]*
  z[445]+z[84]*(z[721]+z[403]*z[722])+z[83]*(z[718]-z[403]*z[719])) + par__Ltibia_M*(
  z[65]*z[450]+z[66]*z[451]+z[68]*(z[725]+z[404]*z[726])+z[67]*(z[723]-z[404]*
  z[724]));
  z[871] = 2*par__Lfemur_M*z[52]*(z[53]+z[54]) + 2*par__Lfoot_M*(z[81]*z[83]+
  z[82]*z[84]) + 2*par__Ltibia_M*(z[65]*z[67]+z[66]*z[68]);
  z[872] = z[307]*z[81]*z[133] + z[307]*z[82]*z[134] + z[307]*z[83]*z[129] + 
  z[307]*z[84]*z[131] + z[309]*z[65]*z[68] + z[309]*z[66]*z[67] + z[309]*
  z[66]*z[127] + z[309]*z[68]*z[126] + z[81]*z[133]*z[356] + z[82]*z[134]*
  z[356] + z[83]*z[129]*z[356] + z[84]*z[131]*z[356] + 2*z[308]*z[83]*z[133] + 
  2*z[308]*z[84]*z[134] + 2*z[310]*z[67]*z[68] + 2*z[310]*z[68]*z[127] + 2*
  z[83]*z[133]*z[357] + 2*z[84]*z[134]*z[357];
  z[873] = par__Lfoot_M*(z[133]*z[444]+z[134]*z[445]+z[84]*(z[758]+z[403]*
  z[760])+z[83]*(z[751]-z[403]*z[753])) + par__Ltibia_M*(z[68]*z[450]+z[127]*
  z[451]+z[68]*(z[763]+z[404]*z[764])+z[67]*(z[761]-z[404]*z[762]));
  z[874] = par__Lfoot_M*(z[83]*z[759]+z[84]*z[752]+z[91]*z[133]+z[94]*z[134]) + 
  par__Ltibia_M*(z[67]*z[747]+z[68]*z[73]+z[68]*z[746]+z[75]*z[127]);
  z[875] = 2*par__Ltibia_M*z[68]*(z[67]+z[127]) + 2*par__Lfoot_M*(z[83]*
  z[133]+z[84]*z[134]);
  z[876] = par__Lfoot_M*(par__Lfoot_CMx*z[134]-par__Lfoot_CMy*z[133]);
  z[877] = par__Lfoot_M*(z[93]*z[133]+z[96]*z[134]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[127]-par__Ltibia_CMy*z[68]);
  z[878] = z[307]*z[81]*z[84] + z[307]*z[82]*z[83] + z[307]*z[82]*z[135] + 
  z[307]*z[84]*z[132] + z[81]*z[84]*z[356] + z[82]*z[83]*z[356] + z[82]*
  z[135]*z[356] + z[84]*z[132]*z[356] + 2*z[308]*z[83]*z[84] + 2*z[308]*z[84]*
  z[135] + 2*z[83]*z[84]*z[357] + 2*z[84]*z[135]*z[357];
  z[879] = par__Lfoot_M*(z[84]*z[444]+z[135]*z[445]+z[84]*(z[781]+z[403]*
  z[782])+z[83]*(z[774]-z[403]*z[776]));
  z[880] = par__Lfoot_M*(z[83]*z[88]+z[84]*z[91]+z[84]*z[775]+z[94]*z[135]);
  z[881] = par__Lfoot_M*z[84]*(z[83]+z[135]);
  z[882] = par__Lfoot_M*(par__Lfoot_CMx*z[135]-par__Lfoot_CMy*z[84]);
  z[883] = par__Lfoot_M*(z[84]*z[93]+z[96]*z[135]+par__Lankle_y*z[79]*z[84]+
  par__Lankle_y*z[80]*z[83]);
  z[884] = z[319]*z[99]*z[136] + z[319]*pow(z[99],2) + z[319]*z[100]*z[101] + 
  z[319]*pow(z[100],2) + 2*z[320]*z[99]*z[100] + 2*z[320]*z[99]*z[101];
  z[885] = par__torso_M*(z[99]*z[510]+z[100]*z[511]+z[99]*(z[793]+z[408]*
  z[794])+z[101]*(z[509]-z[408]*z[792]));
  z[886] = par__torso_M*z[99]*(z[100]+z[101]);
  z[887] = z[367]*z[106] + z[368]*z[6] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + par__Rhip_x*z[1]*z[339] + par__Rhip_x*z[2]*z[340] + par__Rhip_y*
  z[1]*z[340] + z[315]*z[44]*z[113] + z[315]*z[47]*z[115] + z[316]*z[34]*
  z[44] + z[316]*z[35]*z[47] + z[317]*z[26]*z[109] + z[317]*z[28]*z[110] + 
  z[318]*z[18]*z[26] + z[318]*z[19]*z[28] + z[5]*z[13]*z[340] + z[6]*z[13]*
  z[339] + z[34]*z[333]*z[340] + z[35]*z[336]*z[340] + z[113]*z[333]*z[339] + 
  z[115]*z[336]*z[339] - z[369]*z[6] - z[370]*z[5] - par__Rhip_y*z[2]*z[339] - 
  z[6]*z[12]*z[340] - z[12]*z[106]*z[339];
  z[888] = par__Rfoot_M*(z[47]*(z[607]+z[406]*z[608])+z[44]*(z[604]-z[406]*
  z[605])) + par__Rtibia_M*(z[28]*(z[611]+z[407]*z[612])+z[26]*(z[609]-z[407]*
  z[610])) + par__Rfemur_M*(par__Rfemur_CMx*(z[599]+z[405]*z[600])-par__Rfemur_CMy*(
  z[596]-z[405]*z[597]));
  z[889] = z[367]*z[106] + z[368]*z[6] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + z[315]*z[44]*z[113] + z[315]*z[47]*z[115] + z[316]*z[34]*
  z[44] + z[316]*z[35]*z[47] + z[317]*z[26]*z[109] + z[317]*z[28]*z[110] + 
  z[318]*z[18]*z[26] + z[318]*z[19]*z[28] + z[5]*z[13]*z[340] + z[6]*z[13]*
  z[339] + z[34]*z[333]*z[340] + z[35]*z[336]*z[340] + z[113]*z[333]*z[339] + 
  z[115]*z[336]*z[339] - z[369]*z[6] - z[370]*z[5] - z[6]*z[12]*z[340] - 
  z[12]*z[106]*z[339];
  z[890] = par__Rfoot_M*(z[47]*(z[637]+z[406]*z[638])+z[44]*(z[634]-z[406]*
  z[635])) + par__Rtibia_M*(z[28]*(z[641]+z[407]*z[642])+z[26]*(z[639]-z[407]*
  z[640])) + par__Rfemur_M*(par__Rfemur_CMx*(z[629]+z[405]*z[630])-par__Rfemur_CMy*(
  z[466]-z[405]*z[628]));
  z[891] = par__Rankle_y*z[19]*z[339] + par__Rankle_y*z[21]*z[340] + z[315]*
  z[34]*z[675] + z[315]*z[35]*z[668] + z[315]*z[44]*z[114] + z[315]*z[47]*
  z[116] + z[316]*z[36]*z[675] + z[316]*z[37]*z[668] + z[316]*z[44]*z[118] + 
  z[316]*z[47]*z[119] + z[317]*z[18]*z[663] + z[317]*z[19]*z[26] + z[317]*
  z[19]*z[662] + z[317]*z[28]*z[111] + z[318]*z[20]*z[663] + z[318]*z[21]*
  z[26] + z[318]*z[21]*z[662] + z[318]*z[28]*z[112] + z[5]*z[108]*z[339] + 
  z[7]*z[108]*z[340] + z[34]*z[339]*z[675] + z[35]*z[339]*z[668] + z[36]*
  z[340]*z[675] + z[37]*z[340]*z[668] + z[114]*z[333]*z[339] + z[116]*z[336]*
  z[339] + z[118]*z[333]*z[340] + z[119]*z[336]*z[340] - z[5]*z[107]*z[340] - 
  z[6]*z[107]*z[339];
  z[892] = par__Rfoot_M*(z[668]*z[494]+z[675]*z[493]+z[47]*(z[674]+z[406]*
  z[676])+z[44]*(z[667]-z[406]*z[669])) + par__Rtibia_M*(z[662]*z[500]+z[663]*
  z[499]+z[28]*(z[679]+z[407]*z[680])+z[26]*(z[677]-z[407]*z[678]));
  z[893] = 2*par__Rfoot_M*(z[44]*z[675]+z[47]*z[668]) + 2*par__Rtibia_M*(
  z[26]*z[663]+z[28]*z[662]);
  z[894] = par__Rfoot_M*(z[46]*z[675]+z[49]*z[668]) + par__Rtibia_M*(par__Rtibia_CMx*
  z[662]-par__Rtibia_CMy*z[663]);
  z[895] = par__Rfoot_M*(par__Rfoot_CMx*z[668]-par__Rfoot_CMy*z[675]);
  z[896] = z[315]*z[34]*z[41] + z[315]*z[35]*z[44] + z[315]*z[35]*z[691] + 
  z[315]*z[47]*z[117] + z[316]*z[36]*z[41] + z[316]*z[37]*z[44] + z[316]*
  z[37]*z[691] + z[316]*z[47]*z[120] + z[34]*z[41]*z[339] + z[35]*z[333]*
  z[339] + z[35]*z[339]*z[691] + z[36]*z[41]*z[340] + z[37]*z[333]*z[340] + 
  z[37]*z[340]*z[691] + z[117]*z[336]*z[339] + z[120]*z[336]*z[340];
  z[897] = par__Rfoot_M*(z[41]*z[493]+z[691]*z[494]+z[47]*(z[697]+z[406]*
  z[698])+z[44]*(z[690]-z[406]*z[692]));
  z[898] = par__Rfoot_M*(z[41]*z[44]+z[47]*z[691]);
  z[899] = par__Rfoot_M*(z[41]*z[46]+z[49]*z[691]+par__Rankle_y*z[32]*z[47]+
  par__Rankle_y*z[33]*z[44]);
  z[900] = par__Rfoot_M*(par__Rfoot_CMx*z[691]-par__Rfoot_CMy*z[41]);
  z[901] = z[372]*z[110] + z[373]*z[19] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + par__Rhip_x*z[1]*z[339] + par__Rhip_x*z[2]*z[340] + par__Rhip_y*
  z[1]*z[340] + z[315]*z[46]*z[113] + z[315]*z[49]*z[115] + z[316]*z[34]*
  z[46] + z[316]*z[35]*z[49] + z[5]*z[13]*z[340] + z[6]*z[13]*z[339] + z[34]*
  z[335]*z[340] + z[35]*z[338]*z[340] + z[113]*z[335]*z[339] + z[115]*z[338]*
  z[339] - z[374]*z[109] - z[375]*z[18] - par__Rhip_y*z[2]*z[339] - z[6]*
  z[12]*z[340] - z[12]*z[106]*z[339];
  z[902] = par__Rfoot_M*(z[49]*(z[607]+z[406]*z[608])+z[46]*(z[604]-z[406]*
  z[605])) + par__Rtibia_M*(par__Rtibia_CMx*(z[611]+z[407]*z[612])-par__Rtibia_CMy*(
  z[609]-z[407]*z[610]));
  z[903] = z[372]*z[110] + z[373]*z[19] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + z[315]*z[46]*z[113] + z[315]*z[49]*z[115] + z[316]*z[34]*
  z[46] + z[316]*z[35]*z[49] + z[5]*z[13]*z[340] + z[6]*z[13]*z[339] + z[34]*
  z[335]*z[340] + z[35]*z[338]*z[340] + z[113]*z[335]*z[339] + z[115]*z[338]*
  z[339] - z[374]*z[109] - z[375]*z[18] - z[6]*z[12]*z[340] - z[12]*z[106]*
  z[339];
  z[904] = par__Rfoot_M*(z[49]*(z[637]+z[406]*z[638])+z[46]*(z[634]-z[406]*
  z[635])) + par__Rtibia_M*(par__Rtibia_CMx*(z[641]+z[407]*z[642])-par__Rtibia_CMy*(
  z[639]-z[407]*z[640]));
  z[905] = z[372]*z[111] + z[373]*z[112] + par__Rankle_y*z[19]*z[339] + 
  par__Rankle_y*z[21]*z[340] + z[315]*z[46]*z[114] + z[315]*z[49]*z[116] + 
  z[316]*z[46]*z[118] + z[316]*z[49]*z[119] + z[5]*z[108]*z[339] + z[7]*
  z[108]*z[340] + z[114]*z[335]*z[339] + z[116]*z[338]*z[339] + z[118]*z[335]*
  z[340] + z[119]*z[338]*z[340] - z[374]*z[19] - z[375]*z[21] - z[5]*z[107]*
  z[340] - z[6]*z[107]*z[339];
  z[906] = par__Rfoot_M*(z[49]*(z[674]+z[406]*z[676])+z[46]*(z[667]-z[406]*
  z[669])) + par__Rtibia_M*(par__Rtibia_CMx*(z[679]+z[407]*z[680])-par__Rtibia_CMy*(
  z[677]-z[407]*z[678]));
  z[907] = par__Rankle_y*z[315];
  z[908] = par__Rankle_y*z[316];
  z[909] = z[315]*z[35]*z[46] + z[315]*z[49]*z[117] + z[316]*z[37]*z[46] + 
  z[316]*z[49]*z[120] + z[907]*z[32]*z[35] + z[907]*z[33]*z[34] + z[908]*
  z[32]*z[37] + z[908]*z[33]*z[36] + z[35]*z[335]*z[339] + z[37]*z[335]*
  z[340] + z[117]*z[338]*z[339] + z[120]*z[338]*z[340] + par__Rankle_y*z[32]*
  z[35]*z[339] + par__Rankle_y*z[32]*z[37]*z[340] + par__Rankle_y*z[33]*z[34]*
  z[339] + par__Rankle_y*z[33]*z[36]*z[340];
  z[910] = par__Rfoot_M*(par__Rankle_y*z[32]*z[494]+par__Rankle_y*z[33]*
  z[493]+z[49]*(z[697]+z[406]*z[698])+z[46]*(z[690]-z[406]*z[692]));
  z[911] = par__Rankle_y*par__Rfoot_M;
  z[912] = z[911]*(z[32]*z[49]+z[33]*z[46]);
  z[913] = z[911]*(par__Rfoot_CMx*z[32]-par__Rfoot_CMy*z[33]);
  z[914] = z[377]*z[115] + z[378]*z[35] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + par__Rhip_x*z[1]*z[339] + par__Rhip_x*z[2]*z[340] + par__Rhip_y*
  z[1]*z[340] + z[5]*z[13]*z[340] + z[6]*z[13]*z[339] - z[379]*z[113] - 
  z[380]*z[34] - par__Rhip_y*z[2]*z[339] - z[6]*z[12]*z[340] - z[12]*z[106]*
  z[339];
  z[915] = par__Rfoot_M*(par__Rfoot_CMx*(z[607]+z[406]*z[608])-par__Rfoot_CMy*(
  z[604]-z[406]*z[605]));
  z[916] = z[377]*z[115] + z[378]*z[35] + par__Rankle_y*z[18]*z[340] + par__Rankle_y*
  z[109]*z[339] + z[5]*z[13]*z[340] + z[6]*z[13]*z[339] - z[379]*z[113] - 
  z[380]*z[34] - z[6]*z[12]*z[340] - z[12]*z[106]*z[339];
  z[917] = par__Rfoot_M*(par__Rfoot_CMx*(z[637]+z[406]*z[638])-par__Rfoot_CMy*(
  z[634]-z[406]*z[635]));
  z[918] = z[377]*z[116] + z[378]*z[119] + par__Rankle_y*z[19]*z[339] + 
  par__Rankle_y*z[21]*z[340] + z[5]*z[108]*z[339] + z[7]*z[108]*z[340] - 
  z[379]*z[114] - z[380]*z[118] - z[5]*z[107]*z[340] - z[6]*z[107]*z[339];
  z[919] = par__Rfoot_M*(par__Rfoot_CMx*(z[674]+z[406]*z[676])-par__Rfoot_CMy*(
  z[667]-z[406]*z[669]));
  z[920] = z[377]*z[117] + z[378]*z[120] - z[379]*z[35] - z[380]*z[37];
  z[921] = par__Rfoot_M*(par__Rfoot_CMx*(z[697]+z[406]*z[698])-par__Rfoot_CMy*(
  z[690]-z[406]*z[692]));
  z[922] = z[382]*z[121] + z[383]*z[53] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + par__Lhip_x*z[1]*z[356] + par__Lhip_x*z[2]*z[357] + par__Lhip_y*
  z[1]*z[357] + z[307]*z[91]*z[128] + z[307]*z[94]*z[130] + z[308]*z[81]*
  z[91] + z[308]*z[82]*z[94] + z[309]*z[73]*z[124] + z[309]*z[75]*z[125] + 
  z[310]*z[65]*z[73] + z[310]*z[66]*z[75] + z[52]*z[60]*z[357] + z[53]*z[60]*
  z[356] + z[81]*z[350]*z[357] + z[82]*z[353]*z[357] + z[128]*z[350]*z[356] + 
  z[130]*z[353]*z[356] - z[384]*z[53] - z[385]*z[52] - par__Lhip_y*z[2]*
  z[356] - z[53]*z[59]*z[357] - z[59]*z[121]*z[356];
  z[923] = par__Lfoot_M*(z[94]*(z[590]+z[403]*z[591])+z[91]*(z[587]-z[403]*
  z[588])) + par__Ltibia_M*(z[75]*(z[594]+z[404]*z[595])+z[73]*(z[592]-z[404]*
  z[593])) + par__Lfemur_M*(par__Lfemur_CMx*(z[582]+z[402]*z[583])-par__Lfemur_CMy*(
  z[579]-z[402]*z[580]));
  z[924] = z[382]*z[121] + z[383]*z[53] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + z[307]*z[91]*z[128] + z[307]*z[94]*z[130] + z[308]*z[81]*
  z[91] + z[308]*z[82]*z[94] + z[309]*z[73]*z[124] + z[309]*z[75]*z[125] + 
  z[310]*z[65]*z[73] + z[310]*z[66]*z[75] + z[52]*z[60]*z[357] + z[53]*z[60]*
  z[356] + z[81]*z[350]*z[357] + z[82]*z[353]*z[357] + z[128]*z[350]*z[356] + 
  z[130]*z[353]*z[356] - z[384]*z[53] - z[385]*z[52] - z[53]*z[59]*z[357] - 
  z[59]*z[121]*z[356];
  z[925] = par__Lfoot_M*(z[94]*(z[721]+z[403]*z[722])+z[91]*(z[718]-z[403]*
  z[719])) + par__Ltibia_M*(z[75]*(z[725]+z[404]*z[726])+z[73]*(z[723]-z[404]*
  z[724])) + par__Lfemur_M*(par__Lfemur_CMx*(z[713]+z[402]*z[714])-par__Lfemur_CMy*(
  z[417]-z[402]*z[712]));
  z[926] = par__Lankle_y*z[66]*z[356] + par__Lankle_y*z[68]*z[357] + z[307]*
  z[81]*z[759] + z[307]*z[82]*z[752] + z[307]*z[91]*z[129] + z[307]*z[94]*
  z[131] + z[308]*z[83]*z[759] + z[308]*z[84]*z[752] + z[308]*z[91]*z[133] + 
  z[308]*z[94]*z[134] + z[309]*z[65]*z[747] + z[309]*z[66]*z[73] + z[309]*
  z[66]*z[746] + z[309]*z[75]*z[126] + z[310]*z[67]*z[747] + z[310]*z[68]*
  z[73] + z[310]*z[68]*z[746] + z[310]*z[75]*z[127] + z[52]*z[123]*z[356] + 
  z[54]*z[123]*z[357] + z[81]*z[356]*z[759] + z[82]*z[356]*z[752] + z[83]*
  z[357]*z[759] + z[84]*z[357]*z[752] + z[129]*z[350]*z[356] + z[131]*z[353]*
  z[356] + z[133]*z[350]*z[357] + z[134]*z[353]*z[357] - z[52]*z[122]*z[357] - 
  z[53]*z[122]*z[356];
  z[927] = par__Lfoot_M*(z[752]*z[445]+z[759]*z[444]+z[94]*(z[758]+z[403]*
  z[760])+z[91]*(z[751]-z[403]*z[753])) + par__Ltibia_M*(z[746]*z[451]+z[747]*
  z[450]+z[75]*(z[763]+z[404]*z[764])+z[73]*(z[761]-z[404]*z[762]));
  z[928] = 2*par__Lfoot_M*(z[91]*z[759]+z[94]*z[752]) + 2*par__Ltibia_M*(
  z[73]*z[747]+z[75]*z[746]);
  z[929] = par__Lfoot_M*(z[93]*z[759]+z[96]*z[752]) + par__Ltibia_M*(par__Ltibia_CMx*
  z[746]-par__Ltibia_CMy*z[747]);
  z[930] = par__Lfoot_M*(par__Lfoot_CMx*z[752]-par__Lfoot_CMy*z[759]);
  z[931] = z[307]*z[81]*z[88] + z[307]*z[82]*z[91] + z[307]*z[82]*z[775] + 
  z[307]*z[94]*z[132] + z[308]*z[83]*z[88] + z[308]*z[84]*z[91] + z[308]*
  z[84]*z[775] + z[308]*z[94]*z[135] + z[81]*z[88]*z[356] + z[82]*z[350]*
  z[356] + z[82]*z[356]*z[775] + z[83]*z[88]*z[357] + z[84]*z[350]*z[357] + 
  z[84]*z[357]*z[775] + z[132]*z[353]*z[356] + z[135]*z[353]*z[357];
  z[932] = par__Lfoot_M*(z[88]*z[444]+z[775]*z[445]+z[94]*(z[781]+z[403]*
  z[782])+z[91]*(z[774]-z[403]*z[776]));
  z[933] = par__Lfoot_M*(z[88]*z[91]+z[94]*z[775]);
  z[934] = par__Lfoot_M*(z[88]*z[93]+z[96]*z[775]+par__Lankle_y*z[79]*z[94]+
  par__Lankle_y*z[80]*z[91]);
  z[935] = par__Lfoot_M*(par__Lfoot_CMx*z[775]-par__Lfoot_CMy*z[88]);
  z[936] = z[387]*z[125] + z[388]*z[66] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + par__Lhip_x*z[1]*z[356] + par__Lhip_x*z[2]*z[357] + par__Lhip_y*
  z[1]*z[357] + z[307]*z[93]*z[128] + z[307]*z[96]*z[130] + z[308]*z[81]*
  z[93] + z[308]*z[82]*z[96] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] + 
  z[81]*z[352]*z[357] + z[82]*z[355]*z[357] + z[128]*z[352]*z[356] + z[130]*
  z[355]*z[356] - z[389]*z[124] - z[390]*z[65] - par__Lhip_y*z[2]*z[356] - 
  z[53]*z[59]*z[357] - z[59]*z[121]*z[356];
  z[937] = par__Lfoot_M*(z[96]*(z[590]+z[403]*z[591])+z[93]*(z[587]-z[403]*
  z[588])) + par__Ltibia_M*(par__Ltibia_CMx*(z[594]+z[404]*z[595])-par__Ltibia_CMy*(
  z[592]-z[404]*z[593]));
  z[938] = z[387]*z[125] + z[388]*z[66] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + z[307]*z[93]*z[128] + z[307]*z[96]*z[130] + z[308]*z[81]*
  z[93] + z[308]*z[82]*z[96] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] + 
  z[81]*z[352]*z[357] + z[82]*z[355]*z[357] + z[128]*z[352]*z[356] + z[130]*
  z[355]*z[356] - z[389]*z[124] - z[390]*z[65] - z[53]*z[59]*z[357] - z[59]*
  z[121]*z[356];
  z[939] = par__Lfoot_M*(z[96]*(z[721]+z[403]*z[722])+z[93]*(z[718]-z[403]*
  z[719])) + par__Ltibia_M*(par__Ltibia_CMx*(z[725]+z[404]*z[726])-par__Ltibia_CMy*(
  z[723]-z[404]*z[724]));
  z[940] = z[387]*z[126] + z[388]*z[127] + par__Lankle_y*z[66]*z[356] + 
  par__Lankle_y*z[68]*z[357] + z[307]*z[93]*z[129] + z[307]*z[96]*z[131] + 
  z[308]*z[93]*z[133] + z[308]*z[96]*z[134] + z[52]*z[123]*z[356] + z[54]*
  z[123]*z[357] + z[129]*z[352]*z[356] + z[131]*z[355]*z[356] + z[133]*z[352]*
  z[357] + z[134]*z[355]*z[357] - z[389]*z[66] - z[390]*z[68] - z[52]*z[122]*
  z[357] - z[53]*z[122]*z[356];
  z[941] = par__Lfoot_M*(z[96]*(z[758]+z[403]*z[760])+z[93]*(z[751]-z[403]*
  z[753])) + par__Ltibia_M*(par__Ltibia_CMx*(z[763]+z[404]*z[764])-par__Ltibia_CMy*(
  z[761]-z[404]*z[762]));
  z[942] = par__Lankle_y*z[307];
  z[943] = par__Lankle_y*z[308];
  z[944] = z[307]*z[82]*z[93] + z[307]*z[96]*z[132] + z[308]*z[84]*z[93] + 
  z[308]*z[96]*z[135] + z[942]*z[79]*z[82] + z[942]*z[80]*z[81] + z[943]*
  z[79]*z[84] + z[943]*z[80]*z[83] + z[82]*z[352]*z[356] + z[84]*z[352]*
  z[357] + z[132]*z[355]*z[356] + z[135]*z[355]*z[357] + par__Lankle_y*z[79]*
  z[82]*z[356] + par__Lankle_y*z[79]*z[84]*z[357] + par__Lankle_y*z[80]*z[81]*
  z[356] + par__Lankle_y*z[80]*z[83]*z[357];
  z[945] = par__Lfoot_M*(par__Lankle_y*z[79]*z[445]+par__Lankle_y*z[80]*
  z[444]+z[96]*(z[781]+z[403]*z[782])+z[93]*(z[774]-z[403]*z[776]));
  z[946] = par__Lankle_y*par__Lfoot_M;
  z[947] = z[946]*(z[79]*z[96]+z[80]*z[93]);
  z[948] = z[946]*(par__Lfoot_CMx*z[79]-par__Lfoot_CMy*z[80]);
  z[949] = z[392]*z[130] + z[393]*z[82] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + par__Lhip_x*z[1]*z[356] + par__Lhip_x*z[2]*z[357] + par__Lhip_y*
  z[1]*z[357] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] - z[394]*z[128] - 
  z[395]*z[81] - par__Lhip_y*z[2]*z[356] - z[53]*z[59]*z[357] - z[59]*z[121]*
  z[356];
  z[950] = par__Lfoot_M*(par__Lfoot_CMx*(z[590]+z[403]*z[591])-par__Lfoot_CMy*(
  z[587]-z[403]*z[588]));
  z[951] = z[392]*z[130] + z[393]*z[82] + par__Lankle_y*z[65]*z[357] + par__Lankle_y*
  z[124]*z[356] + z[52]*z[60]*z[357] + z[53]*z[60]*z[356] - z[394]*z[128] - 
  z[395]*z[81] - z[53]*z[59]*z[357] - z[59]*z[121]*z[356];
  z[952] = par__Lfoot_M*(par__Lfoot_CMx*(z[721]+z[403]*z[722])-par__Lfoot_CMy*(
  z[718]-z[403]*z[719]));
  z[953] = z[392]*z[131] + z[393]*z[134] + par__Lankle_y*z[66]*z[356] + 
  par__Lankle_y*z[68]*z[357] + z[52]*z[123]*z[356] + z[54]*z[123]*z[357] - 
  z[394]*z[129] - z[395]*z[133] - z[52]*z[122]*z[357] - z[53]*z[122]*z[356];
  z[954] = par__Lfoot_M*(par__Lfoot_CMx*(z[758]+z[403]*z[760])-par__Lfoot_CMy*(
  z[751]-z[403]*z[753]));
  z[955] = z[392]*z[132] + z[393]*z[135] - z[394]*z[82] - z[395]*z[84];
  z[956] = par__Lfoot_M*(par__Lfoot_CMx*(z[781]+z[403]*z[782])-par__Lfoot_CMy*(
  z[774]-z[403]*z[776]));
  z[957] = z[397]*z[136] + z[398]*z[100] - z[399]*z[100] - z[400]*z[99];
  z[958] = par__torso_M*(par__torso_CMx*(z[616]+z[408]*z[617])-par__torso_CMy*(
  z[613]-z[408]*z[614]));
  z[959] = par__torso_M*(par__torso_CMx*(z[793]+z[408]*z[794])-par__torso_CMy*(
  z[509]-z[408]*z[792]));
  z[960] = (z[578]-z[620]-z[526]*q3pp-z[621]*q2pp)/par__bodyweight;
  z[961] = (z[340]+z[357])/par__bodyweight;
  z[962] = (z[339]+z[356])/par__bodyweight;
  z[963] = (z[627]-z[643]-z[644]*q1pp-z[645]*q4pp-z[646]*q5pp-z[647]*q6pp-
  z[648]*q3pp-z[649]*q2pp)/par__bodyweight;
  z[964] = (z[654]-z[681]-z[682]*q1pp-z[683]*q4pp-z[684]*q5pp-z[685]*q6pp-
  z[686]*q3pp-z[687]*q2pp)/par__bodyweight;
  z[965] = (z[689]-z[699]-2*z[700]*q1pp-z[701]*q4pp-z[702]*q5pp-z[703]*q6pp-
  z[704]*q3pp-z[705]*q2pp)/par__bodyweight;
  z[966] = (z[711]-z[727]-z[728]*q1pp-z[729]*q7pp-z[730]*q8pp-z[731]*q9pp-
  z[732]*q3pp-z[733]*q2pp)/par__bodyweight;
  z[967] = (z[738]-z[765]-z[766]*q1pp-z[767]*q7pp-z[768]*q8pp-z[769]*q9pp-
  z[770]*q3pp-z[771]*q2pp)/par__bodyweight;
  z[968] = (z[773]-z[783]-2*z[784]*q1pp-z[785]*q7pp-z[786]*q8pp-z[787]*q9pp-
  z[788]*q3pp-z[789]*q2pp)/par__bodyweight;
  z[969] = (z[791]-z[795]-2*z[796]*q1pp-z[797]*q10pp-z[798]*q3pp-z[799]*q2pp)/
  par__bodyweight;
  z[970] = (z[800]-z[801]-z[621]*q1pp-z[802]*q7pp-z[803]*q3pp-z[804]*q2pp-
  z[805]*q9pp-z[806]*q8pp-z[807]*q4pp-z[808]*q6pp-z[809]*q5pp-z[810]*q10pp)/
  par__bodyweight;
  z[971] = (z[811]-z[812]-z[649]*q1pp-z[807]*q4pp-z[808]*q6pp-z[809]*q5pp-
  z[813]*q3pp-z[814]*q2pp)/par__bodyweight;
  z[972] = (z[815]-z[816]-z[687]*q1pp-z[817]*q3pp-z[818]*q2pp-z[819]*q4pp-
  z[820]*q6pp-z[821]*q5pp)/par__bodyweight;
  z[973] = (z[822]-z[823]-2*z[825]*q2pp-z[705]*q1pp-z[824]*q3pp-z[826]*q4pp-
  z[827]*q6pp-z[828]*q5pp)/par__bodyweight;
  z[974] = (z[829]-z[830]-z[733]*q1pp-z[802]*q7pp-z[805]*q9pp-z[806]*q8pp-
  z[831]*q3pp-z[832]*q2pp)/par__bodyweight;
  z[975] = (z[833]-z[834]-z[771]*q1pp-z[835]*q7pp-z[836]*q3pp-z[837]*q2pp-
  z[838]*q9pp-z[839]*q8pp)/par__bodyweight;
  z[976] = (z[840]-z[841]-2*z[844]*q2pp-z[789]*q1pp-z[842]*q7pp-z[843]*q3pp-
  z[845]*q9pp-z[846]*q8pp)/par__bodyweight;
  z[977] = (z[847]-z[848]-2*z[850]*q2pp-z[799]*q1pp-z[810]*q10pp-z[849]*q3pp)/
  par__bodyweight;
  z[978] = (z[851]-z[852]-z[526]*q1pp-z[528]*q7pp-z[531]*q9pp-z[532]*q8pp-
  z[533]*q4pp-z[534]*q6pp-z[535]*q5pp-z[536]*q10pp-z[803]*q2pp-z[853]*q3pp)/
  par__bodyweight;
  z[979] = (z[854]-z[855]-z[533]*q4pp-z[534]*q6pp-z[535]*q5pp-z[648]*q1pp-
  z[813]*q2pp-z[856]*q3pp)/par__bodyweight;
  z[980] = (z[857]-z[858]-z[686]*q1pp-z[817]*q2pp-z[859]*q3pp-z[860]*q4pp-
  z[861]*q6pp-z[862]*q5pp)/par__bodyweight;
  z[981] = (z[863]-z[864]-2*z[865]*q3pp-z[704]*q1pp-z[824]*q2pp-z[866]*q4pp-
  z[867]*q6pp-z[868]*q5pp)/par__bodyweight;
  z[982] = (z[869]-z[870]-z[528]*q7pp-z[531]*q9pp-z[532]*q8pp-z[732]*q1pp-
  z[831]*q2pp-z[871]*q3pp)/par__bodyweight;
  z[983] = (z[872]-z[873]-z[770]*q1pp-z[836]*q2pp-z[874]*q7pp-z[875]*q3pp-
  z[876]*q9pp-z[877]*q8pp)/par__bodyweight;
  z[984] = (z[878]-z[879]-2*z[881]*q3pp-z[788]*q1pp-z[843]*q2pp-z[880]*q7pp-
  z[882]*q9pp-z[883]*q8pp)/par__bodyweight;
  z[985] = (z[884]-z[885]-2*z[886]*q3pp-z[536]*q10pp-z[798]*q1pp-z[849]*q2pp)/
  par__bodyweight;
  z[986] = (z[887]-z[888]-z[533]*q3pp-z[807]*q2pp)/par__bodyweight;
  z[987] = z[340]/par__bodyweight;
  z[988] = z[339]/par__bodyweight;
  z[989] = (z[889]-z[890]-z[533]*q3pp-z[645]*q1pp-z[807]*q2pp)/par__bodyweight;
  z[990] = (z[891]-z[892]-z[683]*q1pp-z[819]*q2pp-z[860]*q3pp-z[893]*q4pp-
  z[894]*q5pp-z[895]*q6pp)/par__bodyweight;
  z[991] = (z[896]-z[897]-2*z[898]*q4pp-z[701]*q1pp-z[826]*q2pp-z[866]*q3pp-
  z[899]*q5pp-z[900]*q6pp)/par__bodyweight;
  z[992] = (z[901]-z[902]-z[535]*q3pp-z[809]*q2pp)/par__bodyweight;
  z[993] = (z[903]-z[904]-z[535]*q3pp-z[646]*q1pp-z[809]*q2pp)/par__bodyweight;
  z[994] = (z[905]-z[906]-z[684]*q1pp-z[821]*q2pp-z[862]*q3pp-z[894]*q4pp)/
  par__bodyweight;
  z[995] = (z[909]-z[910]-2*z[912]*q5pp-z[702]*q1pp-z[828]*q2pp-z[868]*q3pp-
  z[899]*q4pp-z[913]*q6pp)/par__bodyweight;
  z[996] = (z[914]-z[915]-z[534]*q3pp-z[808]*q2pp)/par__bodyweight;
  z[997] = (z[916]-z[917]-z[534]*q3pp-z[647]*q1pp-z[808]*q2pp)/par__bodyweight;
  z[998] = (z[918]-z[919]-z[685]*q1pp-z[820]*q2pp-z[861]*q3pp-z[895]*q4pp)/
  par__bodyweight;
  z[999] = (z[920]-z[921]-z[703]*q1pp-z[827]*q2pp-z[867]*q3pp-z[900]*q4pp-
  z[913]*q5pp)/par__bodyweight;
  z[1000] = (z[922]-z[923]-z[528]*q3pp-z[802]*q2pp)/par__bodyweight;
  z[1001] = z[357]/par__bodyweight;
  z[1002] = z[356]/par__bodyweight;
  z[1003] = (z[924]-z[925]-z[528]*q3pp-z[729]*q1pp-z[802]*q2pp)/par__bodyweight;
  z[1004] = (z[926]-z[927]-z[767]*q1pp-z[835]*q2pp-z[874]*q3pp-z[928]*q7pp-
  z[929]*q8pp-z[930]*q9pp)/par__bodyweight;
  z[1005] = (z[931]-z[932]-2*z[933]*q7pp-z[785]*q1pp-z[842]*q2pp-z[880]*q3pp-
  z[934]*q8pp-z[935]*q9pp)/par__bodyweight;
  z[1006] = (z[936]-z[937]-z[532]*q3pp-z[806]*q2pp)/par__bodyweight;
  z[1007] = (z[938]-z[939]-z[532]*q3pp-z[730]*q1pp-z[806]*q2pp)/par__bodyweight;
  z[1008] = (z[940]-z[941]-z[768]*q1pp-z[839]*q2pp-z[877]*q3pp-z[929]*q7pp)/
  par__bodyweight;
  z[1009] = (z[944]-z[945]-2*z[947]*q8pp-z[786]*q1pp-z[846]*q2pp-z[883]*q3pp-
  z[934]*q7pp-z[948]*q9pp)/par__bodyweight;
  z[1010] = (z[949]-z[950]-z[531]*q3pp-z[805]*q2pp)/par__bodyweight;
  z[1011] = (z[951]-z[952]-z[531]*q3pp-z[731]*q1pp-z[805]*q2pp)/par__bodyweight;
  z[1012] = (z[953]-z[954]-z[769]*q1pp-z[838]*q2pp-z[876]*q3pp-z[930]*q7pp)/
  par__bodyweight;
  z[1013] = (z[955]-z[956]-z[787]*q1pp-z[845]*q2pp-z[882]*q3pp-z[935]*q7pp-
  z[948]*q8pp)/par__bodyweight;
  z[1014] = (z[957]-z[958]-z[536]*q3pp-z[810]*q2pp)/par__bodyweight;
  z[1015] = (z[957]-z[959]-z[536]*q3pp-z[797]*q1pp-z[810]*q2pp)/par__bodyweight;
  z[1016] = -par__pelvis_CMx*z[2] - par__pelvis_CMy*z[1];
  z[1017] = par__pelvis_CMx*z[1] - par__pelvis_CMy*z[2];
  z[1018] = (sx*z[1016]+sy*z[1017])/pow((1.0E-06+pow(sx,2)+pow(sy,2)),0.5);
  z[1019] = par__airdrag*(par__pelvis_CMx*s+z[1018]*(par__pelvis_CMx*q1p+z[1]*
  q3p-z[2]*q2p));
  z[1020] = par__airdrag*(par__pelvis_CMy*s+z[1018]*(par__pelvis_CMy*q1p-z[1]*
  q2p-z[2]*q3p));
  z[1021] = -z[576]*z[2]*z[1018] - z[577]*z[1]*z[1018] - par__pelvis_CMx*
  z[1019] - par__pelvis_CMy*z[1020];
  z[1022] = z[52]*q3p + z[53]*q2p + z[413];
  z[1023] = z[1022] - z[410] - z[58]*z[402];
  z[1024] = z[53]*q3p + z[121]*q2p - z[416];
  z[1025] = z[409] + z[1024] + z[57]*z[402];
  z[1026] = z[81]*q3p + z[128]*q2p + z[435];
  z[1027] = z[1026] - z[421] - z[95]*z[403];
  z[1028] = z[82]*q3p + z[130]*q2p + z[441];
  z[1029] = z[420] + z[1028] + z[92]*z[403];
  z[1030] = z[65]*q3p + z[124]*q2p + z[434];
  z[1031] = z[1030] - z[447] - z[76]*z[404];
  z[1032] = z[66]*q3p + z[125]*q2p + z[433];
  z[1033] = z[446] + z[1032] + z[74]*z[404];
  z[1034] = z[5]*q3p + z[6]*q2p + z[462];
  z[1035] = z[1034] - z[459] - z[11]*z[405];
  z[1036] = z[6]*q3p + z[106]*q2p - z[465];
  z[1037] = z[458] + z[1036] + z[10]*z[405];
  z[1038] = z[34]*q3p + z[113]*q2p + z[484];
  z[1039] = z[1038] - z[470] - z[48]*z[406];
  z[1040] = z[35]*q3p + z[115]*q2p + z[490];
  z[1041] = z[469] + z[1040] + z[45]*z[406];
  z[1042] = z[18]*q3p + z[109]*q2p + z[483];
  z[1043] = z[1042] - z[496] - z[29]*z[407];
  z[1044] = z[19]*q3p + z[110]*q2p + z[482];
  z[1045] = z[495] + z[1044] + z[27]*z[407];
  z[1046] = z[99]*q3p + z[100]*q2p + z[505];
  z[1047] = z[1046] - z[502] - z[105]*z[408];
  z[1048] = z[100]*q3p + z[136]*q2p - z[508];
  z[1049] = z[501] + z[1048] + z[104]*z[408];
  z[1050] = z[1]*q2p + z[2]*q3p;
  z[1051] = z[452] - par__pelvis_CMy*q1p - z[1050];
  z[1052] = z[618] - par__pelvis_CMx*q1p - z[453];
  z[1053] = par__Lfemur_M*(z[57]*z[1023]+z[58]*z[1025]) + par__Lfoot_M*(z[92]*
  z[1027]+z[95]*z[1029]) + par__Ltibia_M*(z[74]*z[1031]+z[76]*z[1033]) + 
  par__Rfemur_M*(z[10]*z[1035]+z[11]*z[1037]) + par__Rfoot_M*(z[45]*z[1039]+
  z[48]*z[1041]) + par__Rtibia_M*(z[27]*z[1043]+z[29]*z[1045]) + par__torso_M*(
  z[104]*z[1047]+z[105]*z[1049]) + par__pelvis_M*(par__pelvis_CMx*z[1051]-
  par__pelvis_CMy*z[1052]);
  z[1054] = sx*(pow(z[1],2)+pow(z[2],2))/pow((1.0E-06+pow(sx,2)+pow(sy,2)),
  0.5);
  z[1055] = par__airdrag*(s*z[2]-z[1054]*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*
  q2p));
  z[1056] = par__airdrag*(s*z[1]-z[1054]*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*
  q3p));
  z[1057] = par__pelvis_CMx*z[1055] + par__pelvis_CMy*z[1056] - z[576]*z[2]*
  z[1054] - z[577]*z[1]*z[1054];
  z[1058] = par__Lfemur_M*(z[58]*(z[415]+z[52]*z[402])+z[57]*(z[411]-z[53]*
  z[402])) + par__Lfoot_M*(z[95]*(z[438]+z[81]*z[403])+z[92]*(z[424]-z[82]*
  z[403])) + par__Ltibia_M*(z[76]*(z[423]+z[65]*z[404])+z[74]*(z[422]-z[66]*
  z[404])) + par__Rfemur_M*(z[11]*(z[464]+z[5]*z[405])+z[10]*(z[460]-z[6]*
  z[405])) + par__Rfoot_M*(z[48]*(z[487]+z[34]*z[406])+z[45]*(z[473]-z[35]*
  z[406])) + par__Rtibia_M*(z[29]*(z[472]+z[18]*z[407])+z[27]*(z[471]-z[19]*
  z[407])) + par__torso_M*(z[105]*(z[507]+z[99]*z[408])+z[104]*(z[503]-z[100]*
  z[408]));
  z[1059] = sy*(pow(z[1],2)+pow(z[2],2))/pow((1.0E-06+pow(sx,2)+pow(sy,2)),
  0.5);
  z[1060] = par__airdrag*(s*z[1]+z[1059]*(par__pelvis_CMx*q1p+z[1]*q3p-z[2]*
  q2p));
  z[1061] = par__airdrag*(s*z[2]-z[1059]*(par__pelvis_CMy*q1p-z[1]*q2p-z[2]*
  q3p));
  z[1062] = par__pelvis_CMy*z[1061] - z[576]*z[2]*z[1059] - z[577]*z[1]*
  z[1059] - par__pelvis_CMx*z[1060];
  z[1063] = par__Lfemur_M*(z[58]*(z[411]+z[54]*z[402])+z[57]*(z[412]-z[52]*
  z[402])) + par__Lfoot_M*(z[95]*(z[439]+z[83]*z[403])+z[92]*(z[427]-z[84]*
  z[403])) + par__Ltibia_M*(z[76]*(z[426]+z[67]*z[404])+z[74]*(z[425]-z[68]*
  z[404])) + par__Rfemur_M*(z[11]*(z[460]+z[7]*z[405])+z[10]*(z[461]-z[5]*
  z[405])) + par__Rfoot_M*(z[48]*(z[488]+z[36]*z[406])+z[45]*(z[476]-z[37]*
  z[406])) + par__Rtibia_M*(z[29]*(z[475]+z[20]*z[407])+z[27]*(z[474]-z[21]*
  z[407])) + par__torso_M*(z[105]*(z[503]+z[101]*z[408])+z[104]*(z[504]-z[99]*
  z[408]));
  z[1064] = z[630] - z[459] - par__Rfemur_CMx*z[405];
  z[1065] = z[6]*q3p + z[106]*q2p - z[8]*q1p;
  z[1066] = z[458] + z[1065] - par__Rfemur_CMy*z[405];
  z[1067] = z[9]*z[16] - z[8]*z[17];
  z[1068] = -z[8]*z[16] - z[9]*z[17];
  z[1069] = z[32]*z[1067] + z[33]*z[1068];
  z[1070] = z[34]*q3p + z[113]*q2p + z[1069]*q1p + z[479];
  z[1071] = z[1070] - z[470] - z[47]*z[406];
  z[1072] = z[32]*z[1068] - z[33]*z[1067];
  z[1073] = z[35]*q3p + z[115]*q2p + z[1072]*q1p + z[489];
  z[1074] = z[469] + z[1073] + z[44]*z[406];
  z[1075] = z[18]*q3p + z[109]*q2p + z[1067]*q1p - z[478];
  z[1076] = z[1075] - z[496] - z[28]*z[407];
  z[1077] = z[19]*q3p + z[110]*q2p + z[1068]*q1p - z[477];
  z[1078] = z[495] + z[1077] + z[26]*z[407];
  z[1079] = par__Rfemur_M*(z[10]*z[1064]+z[11]*z[1066]) + par__Rfoot_M*(z[45]*
  z[1071]+z[48]*z[1074]) + par__Rtibia_M*(z[27]*z[1076]+z[29]*z[1078]);
  z[1080] = z[16]*z[108] - z[12]*z[16] - z[13]*z[17] - z[17]*z[107];
  z[1081] = z[12]*z[17] - z[13]*z[16] - z[16]*z[107] - z[17]*z[108];
  z[1082] = -z[32]*z[1080] - z[33]*z[1081];
  z[1083] = z[114]*q2p + z[118]*q3p + z[652]*q1p + z[1082]*q4p + z[485];
  z[1084] = z[1083] - z[470] - z[49]*z[406];
  z[1085] = z[33]*z[1080] - z[32]*z[1081];
  z[1086] = z[116]*q2p + z[119]*q3p + z[653]*q1p + z[1085]*q4p + z[491];
  z[1087] = z[469] + z[1086] + z[46]*z[406];
  z[1088] = z[19]*q2p + z[21]*q3p + z[651]*q1p - z[1080]*q4p;
  z[1089] = z[1088] - z[496] - par__Rtibia_CMx*z[407];
  z[1090] = z[111]*q2p + z[112]*q3p + z[650]*q1p - z[1081]*q4p;
  z[1091] = z[495] + z[1090] - par__Rtibia_CMy*z[407];
  z[1092] = par__Rfoot_M*(z[45]*z[1084]+z[48]*z[1087]) + par__Rtibia_M*(z[27]*
  z[1089]+z[29]*z[1091]);
  z[1093] = z[698] - z[470] - par__Rfoot_CMx*z[406];
  z[1094] = z[469] + z[692] - par__Rfoot_CMy*z[406];
  z[1095] = par__Rfoot_M*(z[45]*z[1093]+z[48]*z[1094]);
  z[1096] = z[714] - z[410] - par__Lfemur_CMx*z[402];
  z[1097] = z[53]*q3p + z[121]*q2p - z[55]*q1p;
  z[1098] = z[409] + z[1097] - par__Lfemur_CMy*z[402];
  z[1099] = z[56]*z[63] - z[55]*z[64];
  z[1100] = -z[55]*z[63] - z[56]*z[64];
  z[1101] = z[79]*z[1099] + z[80]*z[1100];
  z[1102] = z[81]*q3p + z[128]*q2p + z[1101]*q1p + z[430];
  z[1103] = z[1102] - z[421] - z[94]*z[403];
  z[1104] = z[79]*z[1100] - z[80]*z[1099];
  z[1105] = z[82]*q3p + z[130]*q2p + z[1104]*q1p + z[440];
  z[1106] = z[420] + z[1105] + z[91]*z[403];
  z[1107] = z[65]*q3p + z[124]*q2p + z[1099]*q1p - z[429];
  z[1108] = z[1107] - z[447] - z[75]*z[404];
  z[1109] = z[66]*q3p + z[125]*q2p + z[1100]*q1p - z[428];
  z[1110] = z[446] + z[1109] + z[73]*z[404];
  z[1111] = par__Lfemur_M*(z[57]*z[1096]+z[58]*z[1098]) + par__Lfoot_M*(z[92]*
  z[1103]+z[95]*z[1106]) + par__Ltibia_M*(z[74]*z[1108]+z[76]*z[1110]);
  z[1112] = z[63]*z[123] - z[59]*z[63] - z[60]*z[64] - z[64]*z[122];
  z[1113] = z[59]*z[64] - z[60]*z[63] - z[63]*z[122] - z[64]*z[123];
  z[1114] = -z[79]*z[1112] - z[80]*z[1113];
  z[1115] = z[129]*q2p + z[133]*q3p + z[736]*q1p + z[1114]*q7p + z[436];
  z[1116] = z[1115] - z[421] - z[96]*z[403];
  z[1117] = z[80]*z[1112] - z[79]*z[1113];
  z[1118] = z[131]*q2p + z[134]*q3p + z[737]*q1p + z[1117]*q7p + z[442];
  z[1119] = z[420] + z[1118] + z[93]*z[403];
  z[1120] = z[66]*q2p + z[68]*q3p + z[735]*q1p - z[1112]*q7p;
  z[1121] = z[1120] - z[447] - par__Ltibia_CMx*z[404];
  z[1122] = z[126]*q2p + z[127]*q3p + z[734]*q1p - z[1113]*q7p;
  z[1123] = z[446] + z[1122] - par__Ltibia_CMy*z[404];
  z[1124] = par__Lfoot_M*(z[92]*z[1116]+z[95]*z[1119]) + par__Ltibia_M*(z[74]*
  z[1121]+z[76]*z[1123]);
  z[1125] = z[782] - z[421] - par__Lfoot_CMx*z[403];
  z[1126] = z[420] + z[776] - par__Lfoot_CMy*z[403];
  z[1127] = par__Lfoot_M*(z[92]*z[1125]+z[95]*z[1126]);
  z[1128] = z[794] - z[502] - par__torso_CMx*z[408];
  z[1129] = z[100]*q3p + z[136]*q2p - z[102]*q1p;
  z[1130] = z[501] + z[1129] - par__torso_CMy*z[408];
  z[1131] = par__torso_M*(z[104]*z[1128]+z[105]*z[1130]);
  z[1132] = z[321]*z[1018]*pow(z[1],2) + z[321]*z[1018]*pow(z[2],2) + z[1]*
  z[1020] + z[2]*z[1019];
  z[1133] = par__Lfemur_M*(z[52]*z[1023]+z[53]*z[1025]) + par__Lfoot_M*(z[81]*
  z[1027]+z[82]*z[1029]) + par__Ltibia_M*(z[65]*z[1031]+z[66]*z[1033]) + 
  par__Rfemur_M*(z[5]*z[1035]+z[6]*z[1037]) + par__Rfoot_M*(z[34]*z[1039]+
  z[35]*z[1041]) + par__Rtibia_M*(z[18]*z[1043]+z[19]*z[1045]) + par__torso_M*(
  z[99]*z[1047]+z[100]*z[1049]) + par__pelvis_M*(z[1]*z[1052]-z[2]*z[1051]);
  z[1134] = z[321]*z[1054]*pow(z[1],2) + z[321]*z[1054]*pow(z[2],2) - z[1]*
  z[1056] - z[2]*z[1055];
  z[1135] = par__Lfemur_M*(z[52]*z[411]+z[53]*z[415]) + par__Lfoot_M*(z[81]*
  z[424]+z[82]*z[438]) + par__Ltibia_M*(z[65]*z[422]+z[66]*z[423]) + par__Rfemur_M*(
  z[5]*z[460]+z[6]*z[464]) + par__Rfoot_M*(z[34]*z[473]+z[35]*z[487]) + 
  par__Rtibia_M*(z[18]*z[471]+z[19]*z[472]) + par__torso_M*(z[99]*z[503]+
  z[100]*z[507]);
  z[1136] = z[321]*z[1059]*pow(z[1],2) + z[321]*z[1059]*pow(z[2],2) + z[2]*
  z[1060] - z[1]*z[1061];
  z[1137] = par__Lfemur_M*(z[53]*(z[411]+z[54]*z[402])+z[52]*(z[412]-z[52]*
  z[402])) + par__Lfoot_M*(z[82]*(z[439]+z[83]*z[403])+z[81]*(z[427]-z[84]*
  z[403])) + par__Ltibia_M*(z[66]*(z[426]+z[67]*z[404])+z[65]*(z[425]-z[68]*
  z[404])) + par__Rfemur_M*(z[6]*(z[460]+z[7]*z[405])+z[5]*(z[461]-z[5]*
  z[405])) + par__Rfoot_M*(z[35]*(z[488]+z[36]*z[406])+z[34]*(z[476]-z[37]*
  z[406])) + par__Rtibia_M*(z[19]*(z[475]+z[20]*z[407])+z[18]*(z[474]-z[21]*
  z[407])) + par__torso_M*(z[100]*(z[503]+z[101]*z[408])+z[99]*(z[504]-z[99]*
  z[408]));
  z[1138] = par__Rfemur_M*(z[5]*z[1064]+z[6]*z[1066]) + par__Rfoot_M*(z[34]*
  z[1071]+z[35]*z[1074]) + par__Rtibia_M*(z[18]*z[1076]+z[19]*z[1078]);
  z[1139] = par__Rfoot_M*(z[34]*z[1084]+z[35]*z[1087]) + par__Rtibia_M*(z[18]*
  z[1089]+z[19]*z[1091]);
  z[1140] = par__Rfoot_M*(z[34]*z[1093]+z[35]*z[1094]);
  z[1141] = par__Lfemur_M*(z[52]*z[1096]+z[53]*z[1098]) + par__Lfoot_M*(z[81]*
  z[1103]+z[82]*z[1106]) + par__Ltibia_M*(z[65]*z[1108]+z[66]*z[1110]);
  z[1142] = par__Lfoot_M*(z[81]*z[1116]+z[82]*z[1119]) + par__Ltibia_M*(z[65]*
  z[1121]+z[66]*z[1123]);
  z[1143] = par__Lfoot_M*(z[81]*z[1125]+z[82]*z[1126]);
  z[1144] = par__torso_M*(z[99]*z[1128]+z[100]*z[1130]);
  z[1145] = z[2]*z[1020] - z[1]*z[1019];
  z[1146] = par__Lfemur_M*(z[52]*z[1025]+z[54]*z[1023]) + par__Lfoot_M*(z[83]*
  z[1027]+z[84]*z[1029]) + par__Ltibia_M*(z[67]*z[1031]+z[68]*z[1033]) + 
  par__pelvis_M*(z[1]*z[1051]+z[2]*z[1052]) + par__Rfemur_M*(z[5]*z[1037]+
  z[7]*z[1035]) + par__Rfoot_M*(z[36]*z[1039]+z[37]*z[1041]) + par__Rtibia_M*(
  z[20]*z[1043]+z[21]*z[1045]) + par__torso_M*(z[99]*z[1049]+z[101]*z[1047]);
  z[1147] = z[1]*z[1055] - z[2]*z[1056];
  z[1148] = par__Lfemur_M*(z[52]*(z[415]+z[52]*z[402])+z[54]*(z[411]-z[53]*
  z[402])) + par__Lfoot_M*(z[84]*(z[438]+z[81]*z[403])+z[83]*(z[424]-z[82]*
  z[403])) + par__Ltibia_M*(z[68]*(z[423]+z[65]*z[404])+z[67]*(z[422]-z[66]*
  z[404])) + par__Rfemur_M*(z[5]*(z[464]+z[5]*z[405])+z[7]*(z[460]-z[6]*
  z[405])) + par__Rfoot_M*(z[37]*(z[487]+z[34]*z[406])+z[36]*(z[473]-z[35]*
  z[406])) + par__Rtibia_M*(z[21]*(z[472]+z[18]*z[407])+z[20]*(z[471]-z[19]*
  z[407])) + par__torso_M*(z[99]*(z[507]+z[99]*z[408])+z[101]*(z[503]-z[100]*
  z[408]));
  z[1149] = -z[1]*z[1060] - z[2]*z[1061];
  z[1150] = par__Lfemur_M*(z[52]*z[411]+z[54]*z[412]) + par__Lfoot_M*(z[83]*
  z[427]+z[84]*z[439]) + par__Ltibia_M*(z[67]*z[425]+z[68]*z[426]) + par__Rfemur_M*(
  z[5]*z[460]+z[7]*z[461]) + par__Rfoot_M*(z[36]*z[476]+z[37]*z[488]) + 
  par__Rtibia_M*(z[20]*z[474]+z[21]*z[475]) + par__torso_M*(z[99]*z[503]+
  z[101]*z[504]);
  z[1151] = par__Rfemur_M*(z[5]*z[1066]+z[7]*z[1064]) + par__Rfoot_M*(z[36]*
  z[1071]+z[37]*z[1074]) + par__Rtibia_M*(z[20]*z[1076]+z[21]*z[1078]);
  z[1152] = par__Rfoot_M*(z[36]*z[1084]+z[37]*z[1087]) + par__Rtibia_M*(z[20]*
  z[1089]+z[21]*z[1091]);
  z[1153] = par__Rfoot_M*(z[36]*z[1093]+z[37]*z[1094]);
  z[1154] = par__Lfemur_M*(z[52]*z[1098]+z[54]*z[1096]) + par__Lfoot_M*(z[83]*
  z[1103]+z[84]*z[1106]) + par__Ltibia_M*(z[67]*z[1108]+z[68]*z[1110]);
  z[1155] = par__Lfoot_M*(z[83]*z[1116]+z[84]*z[1119]) + par__Ltibia_M*(z[67]*
  z[1121]+z[68]*z[1123]);
  z[1156] = par__Lfoot_M*(z[83]*z[1125]+z[84]*z[1126]);
  z[1157] = par__torso_M*(z[99]*z[1130]+z[101]*z[1128]);
  z[1158] = par__Rfoot_M*(z[44]*z[1039]+z[47]*z[1041]) + par__Rtibia_M*(z[26]*
  z[1043]+z[28]*z[1045]) + par__Rfemur_M*(par__Rfemur_CMx*z[1037]-par__Rfemur_CMy*
  z[1035]);
  z[1159] = par__Rfoot_M*(z[47]*(z[487]+z[34]*z[406])+z[44]*(z[473]-z[35]*
  z[406])) + par__Rtibia_M*(z[28]*(z[472]+z[18]*z[407])+z[26]*(z[471]-z[19]*
  z[407])) + par__Rfemur_M*(par__Rfemur_CMx*(z[464]+z[5]*z[405])-par__Rfemur_CMy*(
  z[460]-z[6]*z[405]));
  z[1160] = par__Rfoot_M*(z[47]*(z[488]+z[36]*z[406])+z[44]*(z[476]-z[37]*
  z[406])) + par__Rtibia_M*(z[28]*(z[475]+z[20]*z[407])+z[26]*(z[474]-z[21]*
  z[407])) + par__Rfemur_M*(par__Rfemur_CMx*(z[460]+z[7]*z[405])-par__Rfemur_CMy*(
  z[461]-z[5]*z[405]));
  z[1161] = par__Rfoot_M*(z[44]*z[1071]+z[47]*z[1074]) + par__Rtibia_M*(z[26]*
  z[1076]+z[28]*z[1078]) + par__Rfemur_M*(par__Rfemur_CMx*z[1066]-par__Rfemur_CMy*
  z[1064]);
  z[1162] = par__Rfoot_M*(z[44]*z[1084]+z[47]*z[1087]) + par__Rtibia_M*(z[26]*
  z[1089]+z[28]*z[1091]);
  z[1163] = par__Rfoot_M*(z[44]*z[1093]+z[47]*z[1094]);
  z[1164] = par__Rfoot_M*(z[46]*z[1039]+z[49]*z[1041]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[1045]-par__Rtibia_CMy*z[1043]);
  z[1165] = par__Rfoot_M*(z[49]*(z[487]+z[34]*z[406])+z[46]*(z[473]-z[35]*
  z[406])) + par__Rtibia_M*(par__Rtibia_CMx*(z[472]+z[18]*z[407])-par__Rtibia_CMy*(
  z[471]-z[19]*z[407]));
  z[1166] = par__Rfoot_M*(z[49]*(z[488]+z[36]*z[406])+z[46]*(z[476]-z[37]*
  z[406])) + par__Rtibia_M*(par__Rtibia_CMx*(z[475]+z[20]*z[407])-par__Rtibia_CMy*(
  z[474]-z[21]*z[407]));
  z[1167] = par__Rfoot_M*(z[46]*z[1071]+z[49]*z[1074]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[1078]-par__Rtibia_CMy*z[1076]);
  z[1168] = par__Rfoot_M*(z[46]*z[1084]+z[49]*z[1087]) + par__Rtibia_M*(
  par__Rtibia_CMx*z[1091]-par__Rtibia_CMy*z[1089]);
  z[1169] = par__Rfoot_M*(z[46]*z[1093]+z[49]*z[1094]);
  z[1170] = par__Rfoot_M*(par__Rfoot_CMx*z[1041]-par__Rfoot_CMy*z[1039]);
  z[1171] = par__Rfoot_M*(par__Rfoot_CMx*(z[487]+z[34]*z[406])-par__Rfoot_CMy*(
  z[473]-z[35]*z[406]));
  z[1172] = par__Rfoot_M*(par__Rfoot_CMx*(z[488]+z[36]*z[406])-par__Rfoot_CMy*(
  z[476]-z[37]*z[406]));
  z[1173] = par__Rfoot_M*(par__Rfoot_CMx*z[1074]-par__Rfoot_CMy*z[1071]);
  z[1174] = par__Rfoot_M*(par__Rfoot_CMx*z[1087]-par__Rfoot_CMy*z[1084]);
  z[1175] = par__Rfoot_M*(par__Rfoot_CMx*z[1094]-par__Rfoot_CMy*z[1093]);
  z[1176] = par__Lfoot_M*(z[91]*z[1027]+z[94]*z[1029]) + par__Ltibia_M*(z[73]*
  z[1031]+z[75]*z[1033]) + par__Lfemur_M*(par__Lfemur_CMx*z[1025]-par__Lfemur_CMy*
  z[1023]);
  z[1177] = par__Lfoot_M*(z[94]*(z[438]+z[81]*z[403])+z[91]*(z[424]-z[82]*
  z[403])) + par__Ltibia_M*(z[75]*(z[423]+z[65]*z[404])+z[73]*(z[422]-z[66]*
  z[404])) + par__Lfemur_M*(par__Lfemur_CMx*(z[415]+z[52]*z[402])-par__Lfemur_CMy*(
  z[411]-z[53]*z[402]));
  z[1178] = par__Lfoot_M*(z[94]*(z[439]+z[83]*z[403])+z[91]*(z[427]-z[84]*
  z[403])) + par__Ltibia_M*(z[75]*(z[426]+z[67]*z[404])+z[73]*(z[425]-z[68]*
  z[404])) + par__Lfemur_M*(par__Lfemur_CMx*(z[411]+z[54]*z[402])-par__Lfemur_CMy*(
  z[412]-z[52]*z[402]));
  z[1179] = par__Lfoot_M*(z[91]*z[1103]+z[94]*z[1106]) + par__Ltibia_M*(z[73]*
  z[1108]+z[75]*z[1110]) + par__Lfemur_M*(par__Lfemur_CMx*z[1098]-par__Lfemur_CMy*
  z[1096]);
  z[1180] = par__Lfoot_M*(z[91]*z[1116]+z[94]*z[1119]) + par__Ltibia_M*(z[73]*
  z[1121]+z[75]*z[1123]);
  z[1181] = par__Lfoot_M*(z[91]*z[1125]+z[94]*z[1126]);
  z[1182] = par__Lfoot_M*(z[93]*z[1027]+z[96]*z[1029]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[1033]-par__Ltibia_CMy*z[1031]);
  z[1183] = par__Lfoot_M*(z[96]*(z[438]+z[81]*z[403])+z[93]*(z[424]-z[82]*
  z[403])) + par__Ltibia_M*(par__Ltibia_CMx*(z[423]+z[65]*z[404])-par__Ltibia_CMy*(
  z[422]-z[66]*z[404]));
  z[1184] = par__Lfoot_M*(z[96]*(z[439]+z[83]*z[403])+z[93]*(z[427]-z[84]*
  z[403])) + par__Ltibia_M*(par__Ltibia_CMx*(z[426]+z[67]*z[404])-par__Ltibia_CMy*(
  z[425]-z[68]*z[404]));
  z[1185] = par__Lfoot_M*(z[93]*z[1103]+z[96]*z[1106]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[1110]-par__Ltibia_CMy*z[1108]);
  z[1186] = par__Lfoot_M*(z[93]*z[1116]+z[96]*z[1119]) + par__Ltibia_M*(
  par__Ltibia_CMx*z[1123]-par__Ltibia_CMy*z[1121]);
  z[1187] = par__Lfoot_M*(z[93]*z[1125]+z[96]*z[1126]);
  z[1188] = par__Lfoot_M*(par__Lfoot_CMx*z[1029]-par__Lfoot_CMy*z[1027]);
  z[1189] = par__Lfoot_M*(par__Lfoot_CMx*(z[438]+z[81]*z[403])-par__Lfoot_CMy*(
  z[424]-z[82]*z[403]));
  z[1190] = par__Lfoot_M*(par__Lfoot_CMx*(z[439]+z[83]*z[403])-par__Lfoot_CMy*(
  z[427]-z[84]*z[403]));
  z[1191] = par__Lfoot_M*(par__Lfoot_CMx*z[1106]-par__Lfoot_CMy*z[1103]);
  z[1192] = par__Lfoot_M*(par__Lfoot_CMx*z[1119]-par__Lfoot_CMy*z[1116]);
  z[1193] = par__Lfoot_M*(par__Lfoot_CMx*z[1126]-par__Lfoot_CMy*z[1125]);
  z[1194] = par__torso_M*(par__torso_CMx*z[1049]-par__torso_CMy*z[1047]);
  z[1195] = par__torso_M*(par__torso_CMx*(z[507]+z[99]*z[408])-par__torso_CMy*(
  z[503]-z[100]*z[408]));
  z[1196] = par__torso_M*(par__torso_CMx*(z[503]+z[101]*z[408])-par__torso_CMy*(
  z[504]-z[99]*z[408]));
  z[1197] = par__torso_M*(par__torso_CMx*z[1130]-par__torso_CMy*z[1128]);
  z[1198] = (z[1021]-z[1053])/par__bodyweight;
  z[1199] = (z[1057]-z[1058])/par__bodyweight;
  z[1200] = (z[1062]-z[1063])/par__bodyweight;
  z[1201] = z[1079]/par__bodyweight;
  z[1202] = z[1092]/par__bodyweight;
  z[1203] = z[1095]/par__bodyweight;
  z[1204] = z[1111]/par__bodyweight;
  z[1205] = z[1124]/par__bodyweight;
  z[1206] = z[1127]/par__bodyweight;
  z[1207] = z[1131]/par__bodyweight;
  z[1208] = (z[1132]-z[1133])/par__bodyweight;
  z[1209] = (z[1134]-z[1135])/par__bodyweight;
  z[1210] = (z[1136]-z[1137])/par__bodyweight;
  z[1211] = z[1138]/par__bodyweight;
  z[1212] = z[1139]/par__bodyweight;
  z[1213] = z[1140]/par__bodyweight;
  z[1214] = z[1141]/par__bodyweight;
  z[1215] = z[1142]/par__bodyweight;
  z[1216] = z[1143]/par__bodyweight;
  z[1217] = z[1144]/par__bodyweight;
  z[1218] = (z[1145]-z[1146])/par__bodyweight;
  z[1219] = (z[1147]-z[1148])/par__bodyweight;
  z[1220] = (z[1149]-z[1150])/par__bodyweight;
  z[1221] = z[1151]/par__bodyweight;
  z[1222] = z[1152]/par__bodyweight;
  z[1223] = z[1153]/par__bodyweight;
  z[1224] = z[1154]/par__bodyweight;
  z[1225] = z[1155]/par__bodyweight;
  z[1226] = z[1156]/par__bodyweight;
  z[1227] = z[1157]/par__bodyweight;
  z[1228] = z[1158]/par__bodyweight;
  z[1229] = z[1159]/par__bodyweight;
  z[1230] = z[1160]/par__bodyweight;
  z[1231] = z[1161]/par__bodyweight;
  z[1232] = z[1162]/par__bodyweight;
  z[1233] = z[1163]/par__bodyweight;
  z[1234] = z[1164]/par__bodyweight;
  z[1235] = z[1165]/par__bodyweight;
  z[1236] = z[1166]/par__bodyweight;
  z[1237] = z[1167]/par__bodyweight;
  z[1238] = z[1168]/par__bodyweight;
  z[1239] = z[1169]/par__bodyweight;
  z[1240] = z[1170]/par__bodyweight;
  z[1241] = z[1171]/par__bodyweight;
  z[1242] = z[1172]/par__bodyweight;
  z[1243] = z[1173]/par__bodyweight;
  z[1244] = z[1174]/par__bodyweight;
  z[1245] = z[1175]/par__bodyweight;
  z[1246] = z[1176]/par__bodyweight;
  z[1247] = z[1177]/par__bodyweight;
  z[1248] = z[1178]/par__bodyweight;
  z[1249] = z[1179]/par__bodyweight;
  z[1250] = z[1180]/par__bodyweight;
  z[1251] = z[1181]/par__bodyweight;
  z[1252] = z[1182]/par__bodyweight;
  z[1253] = z[1183]/par__bodyweight;
  z[1254] = z[1184]/par__bodyweight;
  z[1255] = z[1185]/par__bodyweight;
  z[1256] = z[1186]/par__bodyweight;
  z[1257] = z[1187]/par__bodyweight;
  z[1258] = z[1188]/par__bodyweight;
  z[1259] = z[1189]/par__bodyweight;
  z[1260] = z[1190]/par__bodyweight;
  z[1261] = z[1191]/par__bodyweight;
  z[1262] = z[1192]/par__bodyweight;
  z[1263] = z[1193]/par__bodyweight;
  z[1264] = z[1194]/par__bodyweight;
  z[1265] = z[1195]/par__bodyweight;
  z[1266] = z[1196]/par__bodyweight;
  z[1267] = z[1197]/par__bodyweight;
  z[1268] = z[513]/par__bodyweight;
  z[1269] = z[526]/par__bodyweight;
  z[1270] = z[525]/par__bodyweight;
  z[1271] = z[520]/par__bodyweight;
  z[1272] = z[522]/par__bodyweight;
  z[1273] = z[523]/par__bodyweight;
  z[1274] = z[515]/par__bodyweight;
  z[1275] = z[517]/par__bodyweight;
  z[1276] = z[518]/par__bodyweight;
  z[1277] = z[524]/par__bodyweight;
  z[1278] = z[530]/par__bodyweight;
  z[1279] = z[529]/par__bodyweight;
  z[1280] = z[533]/par__bodyweight;
  z[1281] = z[535]/par__bodyweight;
  z[1282] = z[534]/par__bodyweight;
  z[1283] = z[528]/par__bodyweight;
  z[1284] = z[532]/par__bodyweight;
  z[1285] = z[531]/par__bodyweight;
  z[1286] = z[536]/par__bodyweight;
  z[1287] = z[539]/par__bodyweight;
  z[1288] = z[542]/par__bodyweight;
  z[1289] = z[544]/par__bodyweight;
  z[1290] = z[543]/par__bodyweight;
  z[1291] = z[538]/par__bodyweight;
  z[1292] = z[541]/par__bodyweight;
  z[1293] = z[540]/par__bodyweight;
  z[1294] = z[545]/par__bodyweight;
  z[1295] = z[548]/par__bodyweight;
  z[1296] = z[549]/par__bodyweight;
  z[1297] = z[550]/par__bodyweight;
  z[1298] = z[553]/par__bodyweight;
  z[1299] = z[554]/par__bodyweight;
  z[1300] = z[556]/par__bodyweight;
  z[1301] = z[559]/par__bodyweight;
  z[1302] = z[560]/par__bodyweight;
  z[1303] = z[561]/par__bodyweight;
  z[1304] = z[564]/par__bodyweight;
  z[1305] = z[565]/par__bodyweight;
  z[1306] = z[567]/par__bodyweight;
  z[1307] = z[569]/par__bodyweight;
  z[1308] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[334]-z[35]*z[337]);
  z[1309] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[334]-z[37]*z[337]);
  z[1310] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[351]-z[82]*z[354]);
  z[1311] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[351]-z[84]*z[354]);
  z[1312] = par__bodyweight*(z[34]*z[36]+z[35]*z[37]);
  z[1313] = par__bodyweight*(z[81]*z[83]+z[82]*z[84]);
  z[1314] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[333]-z[35]*z[336]);
  z[1315] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[333]-z[37]*z[336]);
  z[1316] = par__bodyweight*(z[6]*z[12]-q3-par__Rankle_y*z[18]-par__Rhip_x*
  z[2]-par__Rhip_y*z[1]-z[5]*z[13]-z[34]*z[335]-z[35]*z[338]);
  z[1317] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]-z[36]*z[335]-z[37]*z[338]);
  z[1318] = par__bodyweight*(q3+par__Rankle_y*z[18]+par__Rhip_x*z[2]+par__Rhip_y*
  z[1]+z[5]*z[13]-z[6]*z[12]);
  z[1319] = par__bodyweight*(q2+par__Rhip_x*z[1]+z[5]*z[12]-par__Rankle_y*
  z[20]-par__Rhip_y*z[2]-z[7]*z[13]);
  z[1320] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[350]-z[82]*z[353]);
  z[1321] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[350]-z[84]*z[353]);
  z[1322] = par__bodyweight*(z[53]*z[59]-q3-par__Lankle_y*z[65]-par__Lhip_x*
  z[2]-par__Lhip_y*z[1]-z[52]*z[60]-z[81]*z[352]-z[82]*z[355]);
  z[1323] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]-z[83]*z[352]-z[84]*z[355]);
  z[1324] = par__bodyweight*(q3+par__Lankle_y*z[65]+par__Lhip_x*z[2]+par__Lhip_y*
  z[1]+z[52]*z[60]-z[53]*z[59]);
  z[1325] = par__bodyweight*(q2+par__Lhip_x*z[1]+z[52]*z[59]-par__Lankle_y*
  z[67]-par__Lhip_y*z[2]-z[54]*z[60]);
  z[1326] = z[1308]/par__bodyweight;
  z[1327] = z[1309]/par__bodyweight;
  z[1328] = z[1310]/par__bodyweight;
  z[1329] = z[1311]/par__bodyweight;
  z[1330] = pow(z[34],2) + pow(z[35],2);
  z[1331] = z[1312]/par__bodyweight;
  z[1332] = pow(z[81],2) + pow(z[82],2);
  z[1333] = z[1313]/par__bodyweight;
  z[1334] = pow(z[36],2) + pow(z[37],2);
  z[1335] = pow(z[83],2) + pow(z[84],2);
  z[1336] = z[1314]/par__bodyweight;
  z[1337] = z[1315]/par__bodyweight;
  z[1338] = z[1316]/par__bodyweight;
  z[1339] = z[1317]/par__bodyweight;
  z[1340] = z[1318]/par__bodyweight;
  z[1341] = z[1319]/par__bodyweight;
  z[1342] = z[1320]/par__bodyweight;
  z[1343] = z[1321]/par__bodyweight;
  z[1344] = z[1322]/par__bodyweight;
  z[1345] = z[1323]/par__bodyweight;
  z[1346] = z[1324]/par__bodyweight;
  z[1347] = z[1325]/par__bodyweight;


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
  fk[44] = z[99];
  fk[45] = z[100];
  fk[46] = z[101];
  fk[47] = z[99];
  dfk_dq[0][0] = 0;
  dfk_dq[0][1] = 1;
  dfk_dq[0][2] = 0;
  dfk_dq[0][3] = 0;
  dfk_dq[0][4] = 0;
  dfk_dq[0][5] = 0;
  dfk_dq[0][6] = 0;
  dfk_dq[0][7] = 0;
  dfk_dq[0][8] = 0;
  dfk_dq[0][9] = 0;
  dfk_dq[1][0] = 0;
  dfk_dq[1][1] = 0;
  dfk_dq[1][2] = 1;
  dfk_dq[1][3] = 0;
  dfk_dq[1][4] = 0;
  dfk_dq[1][5] = 0;
  dfk_dq[1][6] = 0;
  dfk_dq[1][7] = 0;
  dfk_dq[1][8] = 0;
  dfk_dq[1][9] = 0;
  dfk_dq[2][0] = -z[2];
  dfk_dq[2][1] = 0;
  dfk_dq[2][2] = 0;
  dfk_dq[2][3] = 0;
  dfk_dq[2][4] = 0;
  dfk_dq[2][5] = 0;
  dfk_dq[2][6] = 0;
  dfk_dq[2][7] = 0;
  dfk_dq[2][8] = 0;
  dfk_dq[2][9] = 0;
  dfk_dq[3][0] = -z[1];
  dfk_dq[3][1] = 0;
  dfk_dq[3][2] = 0;
  dfk_dq[3][3] = 0;
  dfk_dq[3][4] = 0;
  dfk_dq[3][5] = 0;
  dfk_dq[3][6] = 0;
  dfk_dq[3][7] = 0;
  dfk_dq[3][8] = 0;
  dfk_dq[3][9] = 0;
  dfk_dq[4][0] = z[1];
  dfk_dq[4][1] = 0;
  dfk_dq[4][2] = 0;
  dfk_dq[4][3] = 0;
  dfk_dq[4][4] = 0;
  dfk_dq[4][5] = 0;
  dfk_dq[4][6] = 0;
  dfk_dq[4][7] = 0;
  dfk_dq[4][8] = 0;
  dfk_dq[4][9] = 0;
  dfk_dq[5][0] = -z[2];
  dfk_dq[5][1] = 0;
  dfk_dq[5][2] = 0;
  dfk_dq[5][3] = 0;
  dfk_dq[5][4] = 0;
  dfk_dq[5][5] = 0;
  dfk_dq[5][6] = 0;
  dfk_dq[5][7] = 0;
  dfk_dq[5][8] = 0;
  dfk_dq[5][9] = 0;
  dfk_dq[6][0] = z[137];
  dfk_dq[6][1] = 1;
  dfk_dq[6][2] = 0;
  dfk_dq[6][3] = 0;
  dfk_dq[6][4] = 0;
  dfk_dq[6][5] = 0;
  dfk_dq[6][6] = 0;
  dfk_dq[6][7] = 0;
  dfk_dq[6][8] = 0;
  dfk_dq[6][9] = 0;
  dfk_dq[7][0] = z[138];
  dfk_dq[7][1] = 0;
  dfk_dq[7][2] = 1;
  dfk_dq[7][3] = 0;
  dfk_dq[7][4] = 0;
  dfk_dq[7][5] = 0;
  dfk_dq[7][6] = 0;
  dfk_dq[7][7] = 0;
  dfk_dq[7][8] = 0;
  dfk_dq[7][9] = 0;
  dfk_dq[8][0] = z[6];
  dfk_dq[8][1] = 0;
  dfk_dq[8][2] = 0;
  dfk_dq[8][3] = z[6];
  dfk_dq[8][4] = 0;
  dfk_dq[8][5] = 0;
  dfk_dq[8][6] = 0;
  dfk_dq[8][7] = 0;
  dfk_dq[8][8] = 0;
  dfk_dq[8][9] = 0;
  dfk_dq[9][0] = z[106];
  dfk_dq[9][1] = 0;
  dfk_dq[9][2] = 0;
  dfk_dq[9][3] = z[106];
  dfk_dq[9][4] = 0;
  dfk_dq[9][5] = 0;
  dfk_dq[9][6] = 0;
  dfk_dq[9][7] = 0;
  dfk_dq[9][8] = 0;
  dfk_dq[9][9] = 0;
  dfk_dq[10][0] = z[5];
  dfk_dq[10][1] = 0;
  dfk_dq[10][2] = 0;
  dfk_dq[10][3] = z[5];
  dfk_dq[10][4] = 0;
  dfk_dq[10][5] = 0;
  dfk_dq[10][6] = 0;
  dfk_dq[10][7] = 0;
  dfk_dq[10][8] = 0;
  dfk_dq[10][9] = 0;
  dfk_dq[11][0] = z[6];
  dfk_dq[11][1] = 0;
  dfk_dq[11][2] = 0;
  dfk_dq[11][3] = z[6];
  dfk_dq[11][4] = 0;
  dfk_dq[11][5] = 0;
  dfk_dq[11][6] = 0;
  dfk_dq[11][7] = 0;
  dfk_dq[11][8] = 0;
  dfk_dq[11][9] = 0;
  dfk_dq[12][0] = z[139];
  dfk_dq[12][1] = 1;
  dfk_dq[12][2] = 0;
  dfk_dq[12][3] = z[140];
  dfk_dq[12][4] = z[141];
  dfk_dq[12][5] = 0;
  dfk_dq[12][6] = 0;
  dfk_dq[12][7] = 0;
  dfk_dq[12][8] = 0;
  dfk_dq[12][9] = 0;
  dfk_dq[13][0] = z[142];
  dfk_dq[13][1] = 0;
  dfk_dq[13][2] = 1;
  dfk_dq[13][3] = z[143];
  dfk_dq[13][4] = z[144];
  dfk_dq[13][5] = 0;
  dfk_dq[13][6] = 0;
  dfk_dq[13][7] = 0;
  dfk_dq[13][8] = 0;
  dfk_dq[13][9] = 0;
  dfk_dq[14][0] = z[109];
  dfk_dq[14][1] = 0;
  dfk_dq[14][2] = 0;
  dfk_dq[14][3] = z[109];
  dfk_dq[14][4] = z[19];
  dfk_dq[14][5] = 0;
  dfk_dq[14][6] = 0;
  dfk_dq[14][7] = 0;
  dfk_dq[14][8] = 0;
  dfk_dq[14][9] = 0;
  dfk_dq[15][0] = z[110];
  dfk_dq[15][1] = 0;
  dfk_dq[15][2] = 0;
  dfk_dq[15][3] = z[110];
  dfk_dq[15][4] = z[111];
  dfk_dq[15][5] = 0;
  dfk_dq[15][6] = 0;
  dfk_dq[15][7] = 0;
  dfk_dq[15][8] = 0;
  dfk_dq[15][9] = 0;
  dfk_dq[16][0] = z[18];
  dfk_dq[16][1] = 0;
  dfk_dq[16][2] = 0;
  dfk_dq[16][3] = z[18];
  dfk_dq[16][4] = z[21];
  dfk_dq[16][5] = 0;
  dfk_dq[16][6] = 0;
  dfk_dq[16][7] = 0;
  dfk_dq[16][8] = 0;
  dfk_dq[16][9] = 0;
  dfk_dq[17][0] = z[19];
  dfk_dq[17][1] = 0;
  dfk_dq[17][2] = 0;
  dfk_dq[17][3] = z[19];
  dfk_dq[17][4] = z[112];
  dfk_dq[17][5] = 0;
  dfk_dq[17][6] = 0;
  dfk_dq[17][7] = 0;
  dfk_dq[17][8] = 0;
  dfk_dq[17][9] = 0;
  dfk_dq[18][0] = z[145];
  dfk_dq[18][1] = 1;
  dfk_dq[18][2] = 0;
  dfk_dq[18][3] = z[146];
  dfk_dq[18][4] = z[147];
  dfk_dq[18][5] = 0;
  dfk_dq[18][6] = 0;
  dfk_dq[18][7] = 0;
  dfk_dq[18][8] = 0;
  dfk_dq[18][9] = 0;
  dfk_dq[19][0] = z[148];
  dfk_dq[19][1] = 0;
  dfk_dq[19][2] = 1;
  dfk_dq[19][3] = z[149];
  dfk_dq[19][4] = z[150];
  dfk_dq[19][5] = 0;
  dfk_dq[19][6] = 0;
  dfk_dq[19][7] = 0;
  dfk_dq[19][8] = 0;
  dfk_dq[19][9] = 0;
  dfk_dq[20][0] = z[113];
  dfk_dq[20][1] = 0;
  dfk_dq[20][2] = 0;
  dfk_dq[20][3] = z[113];
  dfk_dq[20][4] = z[114];
  dfk_dq[20][5] = z[35];
  dfk_dq[20][6] = 0;
  dfk_dq[20][7] = 0;
  dfk_dq[20][8] = 0;
  dfk_dq[20][9] = 0;
  dfk_dq[21][0] = z[115];
  dfk_dq[21][1] = 0;
  dfk_dq[21][2] = 0;
  dfk_dq[21][3] = z[115];
  dfk_dq[21][4] = z[116];
  dfk_dq[21][5] = z[117];
  dfk_dq[21][6] = 0;
  dfk_dq[21][7] = 0;
  dfk_dq[21][8] = 0;
  dfk_dq[21][9] = 0;
  dfk_dq[22][0] = z[34];
  dfk_dq[22][1] = 0;
  dfk_dq[22][2] = 0;
  dfk_dq[22][3] = z[34];
  dfk_dq[22][4] = z[118];
  dfk_dq[22][5] = z[37];
  dfk_dq[22][6] = 0;
  dfk_dq[22][7] = 0;
  dfk_dq[22][8] = 0;
  dfk_dq[22][9] = 0;
  dfk_dq[23][0] = z[35];
  dfk_dq[23][1] = 0;
  dfk_dq[23][2] = 0;
  dfk_dq[23][3] = z[35];
  dfk_dq[23][4] = z[119];
  dfk_dq[23][5] = z[120];
  dfk_dq[23][6] = 0;
  dfk_dq[23][7] = 0;
  dfk_dq[23][8] = 0;
  dfk_dq[23][9] = 0;
  dfk_dq[24][0] = z[151];
  dfk_dq[24][1] = 1;
  dfk_dq[24][2] = 0;
  dfk_dq[24][3] = 0;
  dfk_dq[24][4] = 0;
  dfk_dq[24][5] = 0;
  dfk_dq[24][6] = 0;
  dfk_dq[24][7] = 0;
  dfk_dq[24][8] = 0;
  dfk_dq[24][9] = 0;
  dfk_dq[25][0] = z[152];
  dfk_dq[25][1] = 0;
  dfk_dq[25][2] = 1;
  dfk_dq[25][3] = 0;
  dfk_dq[25][4] = 0;
  dfk_dq[25][5] = 0;
  dfk_dq[25][6] = 0;
  dfk_dq[25][7] = 0;
  dfk_dq[25][8] = 0;
  dfk_dq[25][9] = 0;
  dfk_dq[26][0] = z[53];
  dfk_dq[26][1] = 0;
  dfk_dq[26][2] = 0;
  dfk_dq[26][3] = 0;
  dfk_dq[26][4] = 0;
  dfk_dq[26][5] = 0;
  dfk_dq[26][6] = z[53];
  dfk_dq[26][7] = 0;
  dfk_dq[26][8] = 0;
  dfk_dq[26][9] = 0;
  dfk_dq[27][0] = z[121];
  dfk_dq[27][1] = 0;
  dfk_dq[27][2] = 0;
  dfk_dq[27][3] = 0;
  dfk_dq[27][4] = 0;
  dfk_dq[27][5] = 0;
  dfk_dq[27][6] = z[121];
  dfk_dq[27][7] = 0;
  dfk_dq[27][8] = 0;
  dfk_dq[27][9] = 0;
  dfk_dq[28][0] = z[52];
  dfk_dq[28][1] = 0;
  dfk_dq[28][2] = 0;
  dfk_dq[28][3] = 0;
  dfk_dq[28][4] = 0;
  dfk_dq[28][5] = 0;
  dfk_dq[28][6] = z[52];
  dfk_dq[28][7] = 0;
  dfk_dq[28][8] = 0;
  dfk_dq[28][9] = 0;
  dfk_dq[29][0] = z[53];
  dfk_dq[29][1] = 0;
  dfk_dq[29][2] = 0;
  dfk_dq[29][3] = 0;
  dfk_dq[29][4] = 0;
  dfk_dq[29][5] = 0;
  dfk_dq[29][6] = z[53];
  dfk_dq[29][7] = 0;
  dfk_dq[29][8] = 0;
  dfk_dq[29][9] = 0;
  dfk_dq[30][0] = z[153];
  dfk_dq[30][1] = 1;
  dfk_dq[30][2] = 0;
  dfk_dq[30][3] = 0;
  dfk_dq[30][4] = 0;
  dfk_dq[30][5] = 0;
  dfk_dq[30][6] = z[154];
  dfk_dq[30][7] = z[155];
  dfk_dq[30][8] = 0;
  dfk_dq[30][9] = 0;
  dfk_dq[31][0] = z[156];
  dfk_dq[31][1] = 0;
  dfk_dq[31][2] = 1;
  dfk_dq[31][3] = 0;
  dfk_dq[31][4] = 0;
  dfk_dq[31][5] = 0;
  dfk_dq[31][6] = z[157];
  dfk_dq[31][7] = z[158];
  dfk_dq[31][8] = 0;
  dfk_dq[31][9] = 0;
  dfk_dq[32][0] = z[124];
  dfk_dq[32][1] = 0;
  dfk_dq[32][2] = 0;
  dfk_dq[32][3] = 0;
  dfk_dq[32][4] = 0;
  dfk_dq[32][5] = 0;
  dfk_dq[32][6] = z[124];
  dfk_dq[32][7] = z[66];
  dfk_dq[32][8] = 0;
  dfk_dq[32][9] = 0;
  dfk_dq[33][0] = z[125];
  dfk_dq[33][1] = 0;
  dfk_dq[33][2] = 0;
  dfk_dq[33][3] = 0;
  dfk_dq[33][4] = 0;
  dfk_dq[33][5] = 0;
  dfk_dq[33][6] = z[125];
  dfk_dq[33][7] = z[126];
  dfk_dq[33][8] = 0;
  dfk_dq[33][9] = 0;
  dfk_dq[34][0] = z[65];
  dfk_dq[34][1] = 0;
  dfk_dq[34][2] = 0;
  dfk_dq[34][3] = 0;
  dfk_dq[34][4] = 0;
  dfk_dq[34][5] = 0;
  dfk_dq[34][6] = z[65];
  dfk_dq[34][7] = z[68];
  dfk_dq[34][8] = 0;
  dfk_dq[34][9] = 0;
  dfk_dq[35][0] = z[66];
  dfk_dq[35][1] = 0;
  dfk_dq[35][2] = 0;
  dfk_dq[35][3] = 0;
  dfk_dq[35][4] = 0;
  dfk_dq[35][5] = 0;
  dfk_dq[35][6] = z[66];
  dfk_dq[35][7] = z[127];
  dfk_dq[35][8] = 0;
  dfk_dq[35][9] = 0;
  dfk_dq[36][0] = z[159];
  dfk_dq[36][1] = 1;
  dfk_dq[36][2] = 0;
  dfk_dq[36][3] = 0;
  dfk_dq[36][4] = 0;
  dfk_dq[36][5] = 0;
  dfk_dq[36][6] = z[160];
  dfk_dq[36][7] = z[161];
  dfk_dq[36][8] = 0;
  dfk_dq[36][9] = 0;
  dfk_dq[37][0] = z[162];
  dfk_dq[37][1] = 0;
  dfk_dq[37][2] = 1;
  dfk_dq[37][3] = 0;
  dfk_dq[37][4] = 0;
  dfk_dq[37][5] = 0;
  dfk_dq[37][6] = z[163];
  dfk_dq[37][7] = z[164];
  dfk_dq[37][8] = 0;
  dfk_dq[37][9] = 0;
  dfk_dq[38][0] = z[128];
  dfk_dq[38][1] = 0;
  dfk_dq[38][2] = 0;
  dfk_dq[38][3] = 0;
  dfk_dq[38][4] = 0;
  dfk_dq[38][5] = 0;
  dfk_dq[38][6] = z[128];
  dfk_dq[38][7] = z[129];
  dfk_dq[38][8] = z[82];
  dfk_dq[38][9] = 0;
  dfk_dq[39][0] = z[130];
  dfk_dq[39][1] = 0;
  dfk_dq[39][2] = 0;
  dfk_dq[39][3] = 0;
  dfk_dq[39][4] = 0;
  dfk_dq[39][5] = 0;
  dfk_dq[39][6] = z[130];
  dfk_dq[39][7] = z[131];
  dfk_dq[39][8] = z[132];
  dfk_dq[39][9] = 0;
  dfk_dq[40][0] = z[81];
  dfk_dq[40][1] = 0;
  dfk_dq[40][2] = 0;
  dfk_dq[40][3] = 0;
  dfk_dq[40][4] = 0;
  dfk_dq[40][5] = 0;
  dfk_dq[40][6] = z[81];
  dfk_dq[40][7] = z[133];
  dfk_dq[40][8] = z[84];
  dfk_dq[40][9] = 0;
  dfk_dq[41][0] = z[82];
  dfk_dq[41][1] = 0;
  dfk_dq[41][2] = 0;
  dfk_dq[41][3] = 0;
  dfk_dq[41][4] = 0;
  dfk_dq[41][5] = 0;
  dfk_dq[41][6] = z[82];
  dfk_dq[41][7] = z[134];
  dfk_dq[41][8] = z[135];
  dfk_dq[41][9] = 0;
  dfk_dq[42][0] = z[165];
  dfk_dq[42][1] = 1;
  dfk_dq[42][2] = 0;
  dfk_dq[42][3] = 0;
  dfk_dq[42][4] = 0;
  dfk_dq[42][5] = 0;
  dfk_dq[42][6] = 0;
  dfk_dq[42][7] = 0;
  dfk_dq[42][8] = 0;
  dfk_dq[42][9] = 0;
  dfk_dq[43][0] = z[166];
  dfk_dq[43][1] = 0;
  dfk_dq[43][2] = 1;
  dfk_dq[43][3] = 0;
  dfk_dq[43][4] = 0;
  dfk_dq[43][5] = 0;
  dfk_dq[43][6] = 0;
  dfk_dq[43][7] = 0;
  dfk_dq[43][8] = 0;
  dfk_dq[43][9] = 0;
  dfk_dq[44][0] = z[100];
  dfk_dq[44][1] = 0;
  dfk_dq[44][2] = 0;
  dfk_dq[44][3] = 0;
  dfk_dq[44][4] = 0;
  dfk_dq[44][5] = 0;
  dfk_dq[44][6] = 0;
  dfk_dq[44][7] = 0;
  dfk_dq[44][8] = 0;
  dfk_dq[44][9] = z[100];
  dfk_dq[45][0] = z[136];
  dfk_dq[45][1] = 0;
  dfk_dq[45][2] = 0;
  dfk_dq[45][3] = 0;
  dfk_dq[45][4] = 0;
  dfk_dq[45][5] = 0;
  dfk_dq[45][6] = 0;
  dfk_dq[45][7] = 0;
  dfk_dq[45][8] = 0;
  dfk_dq[45][9] = z[136];
  dfk_dq[46][0] = z[99];
  dfk_dq[46][1] = 0;
  dfk_dq[46][2] = 0;
  dfk_dq[46][3] = 0;
  dfk_dq[46][4] = 0;
  dfk_dq[46][5] = 0;
  dfk_dq[46][6] = 0;
  dfk_dq[46][7] = 0;
  dfk_dq[46][8] = 0;
  dfk_dq[46][9] = z[99];
  dfk_dq[47][0] = z[100];
  dfk_dq[47][1] = 0;
  dfk_dq[47][2] = 0;
  dfk_dq[47][3] = 0;
  dfk_dq[47][4] = 0;
  dfk_dq[47][5] = 0;
  dfk_dq[47][6] = 0;
  dfk_dq[47][7] = 0;
  dfk_dq[47][8] = 0;
  dfk_dq[47][9] = z[100];
  fkdot[0] = q2p;
  fkdot[1] = q3p;
  fkdot[2] = -z[167];
  fkdot[3] = -z[168];
  fkdot[4] = z[168];
  fkdot[5] = -z[167];
  fkdot[6] = z[169];
  fkdot[7] = z[170];
  fkdot[8] = z[171];
  fkdot[9] = z[172];
  fkdot[10] = z[173];
  fkdot[11] = z[171];
  fkdot[12] = z[174];
  fkdot[13] = z[175];
  fkdot[14] = z[176];
  fkdot[15] = z[177];
  fkdot[16] = z[178];
  fkdot[17] = z[179];
  fkdot[18] = z[180];
  fkdot[19] = z[181];
  fkdot[20] = z[182];
  fkdot[21] = z[183];
  fkdot[22] = z[184];
  fkdot[23] = z[185];
  fkdot[24] = z[186];
  fkdot[25] = z[187];
  fkdot[26] = z[188];
  fkdot[27] = z[189];
  fkdot[28] = z[190];
  fkdot[29] = z[188];
  fkdot[30] = z[191];
  fkdot[31] = z[192];
  fkdot[32] = z[193];
  fkdot[33] = z[194];
  fkdot[34] = z[195];
  fkdot[35] = z[196];
  fkdot[36] = z[197];
  fkdot[37] = z[198];
  fkdot[38] = z[199];
  fkdot[39] = z[200];
  fkdot[40] = z[201];
  fkdot[41] = z[202];
  fkdot[42] = z[203];
  fkdot[43] = z[204];
  fkdot[44] = z[205];
  fkdot[45] = z[206];
  fkdot[46] = z[207];
  fkdot[47] = z[205];
  dfkdot_dq[0][0] = 0;
  dfkdot_dq[0][1] = 0;
  dfkdot_dq[0][2] = 0;
  dfkdot_dq[0][3] = 0;
  dfkdot_dq[0][4] = 0;
  dfkdot_dq[0][5] = 0;
  dfkdot_dq[0][6] = 0;
  dfkdot_dq[0][7] = 0;
  dfkdot_dq[0][8] = 0;
  dfkdot_dq[0][9] = 0;
  dfkdot_dq[1][0] = 0;
  dfkdot_dq[1][1] = 0;
  dfkdot_dq[1][2] = 0;
  dfkdot_dq[1][3] = 0;
  dfkdot_dq[1][4] = 0;
  dfkdot_dq[1][5] = 0;
  dfkdot_dq[1][6] = 0;
  dfkdot_dq[1][7] = 0;
  dfkdot_dq[1][8] = 0;
  dfkdot_dq[1][9] = 0;
  dfkdot_dq[2][0] = -z[168];
  dfkdot_dq[2][1] = 0;
  dfkdot_dq[2][2] = 0;
  dfkdot_dq[2][3] = 0;
  dfkdot_dq[2][4] = 0;
  dfkdot_dq[2][5] = 0;
  dfkdot_dq[2][6] = 0;
  dfkdot_dq[2][7] = 0;
  dfkdot_dq[2][8] = 0;
  dfkdot_dq[2][9] = 0;
  dfkdot_dq[3][0] = z[167];
  dfkdot_dq[3][1] = 0;
  dfkdot_dq[3][2] = 0;
  dfkdot_dq[3][3] = 0;
  dfkdot_dq[3][4] = 0;
  dfkdot_dq[3][5] = 0;
  dfkdot_dq[3][6] = 0;
  dfkdot_dq[3][7] = 0;
  dfkdot_dq[3][8] = 0;
  dfkdot_dq[3][9] = 0;
  dfkdot_dq[4][0] = -z[167];
  dfkdot_dq[4][1] = 0;
  dfkdot_dq[4][2] = 0;
  dfkdot_dq[4][3] = 0;
  dfkdot_dq[4][4] = 0;
  dfkdot_dq[4][5] = 0;
  dfkdot_dq[4][6] = 0;
  dfkdot_dq[4][7] = 0;
  dfkdot_dq[4][8] = 0;
  dfkdot_dq[4][9] = 0;
  dfkdot_dq[5][0] = -z[168];
  dfkdot_dq[5][1] = 0;
  dfkdot_dq[5][2] = 0;
  dfkdot_dq[5][3] = 0;
  dfkdot_dq[5][4] = 0;
  dfkdot_dq[5][5] = 0;
  dfkdot_dq[5][6] = 0;
  dfkdot_dq[5][7] = 0;
  dfkdot_dq[5][8] = 0;
  dfkdot_dq[5][9] = 0;
  dfkdot_dq[6][0] = -z[208];
  dfkdot_dq[6][1] = 0;
  dfkdot_dq[6][2] = 0;
  dfkdot_dq[6][3] = 0;
  dfkdot_dq[6][4] = 0;
  dfkdot_dq[6][5] = 0;
  dfkdot_dq[6][6] = 0;
  dfkdot_dq[6][7] = 0;
  dfkdot_dq[6][8] = 0;
  dfkdot_dq[6][9] = 0;
  dfkdot_dq[7][0] = -z[209];
  dfkdot_dq[7][1] = 0;
  dfkdot_dq[7][2] = 0;
  dfkdot_dq[7][3] = 0;
  dfkdot_dq[7][4] = 0;
  dfkdot_dq[7][5] = 0;
  dfkdot_dq[7][6] = 0;
  dfkdot_dq[7][7] = 0;
  dfkdot_dq[7][8] = 0;
  dfkdot_dq[7][9] = 0;
  dfkdot_dq[8][0] = z[172];
  dfkdot_dq[8][1] = 0;
  dfkdot_dq[8][2] = 0;
  dfkdot_dq[8][3] = z[172];
  dfkdot_dq[8][4] = 0;
  dfkdot_dq[8][5] = 0;
  dfkdot_dq[8][6] = 0;
  dfkdot_dq[8][7] = 0;
  dfkdot_dq[8][8] = 0;
  dfkdot_dq[8][9] = 0;
  dfkdot_dq[9][0] = z[302];
  dfkdot_dq[9][1] = 0;
  dfkdot_dq[9][2] = 0;
  dfkdot_dq[9][3] = z[302];
  dfkdot_dq[9][4] = 0;
  dfkdot_dq[9][5] = 0;
  dfkdot_dq[9][6] = 0;
  dfkdot_dq[9][7] = 0;
  dfkdot_dq[9][8] = 0;
  dfkdot_dq[9][9] = 0;
  dfkdot_dq[10][0] = z[171];
  dfkdot_dq[10][1] = 0;
  dfkdot_dq[10][2] = 0;
  dfkdot_dq[10][3] = z[171];
  dfkdot_dq[10][4] = 0;
  dfkdot_dq[10][5] = 0;
  dfkdot_dq[10][6] = 0;
  dfkdot_dq[10][7] = 0;
  dfkdot_dq[10][8] = 0;
  dfkdot_dq[10][9] = 0;
  dfkdot_dq[11][0] = z[172];
  dfkdot_dq[11][1] = 0;
  dfkdot_dq[11][2] = 0;
  dfkdot_dq[11][3] = z[172];
  dfkdot_dq[11][4] = 0;
  dfkdot_dq[11][5] = 0;
  dfkdot_dq[11][6] = 0;
  dfkdot_dq[11][7] = 0;
  dfkdot_dq[11][8] = 0;
  dfkdot_dq[11][9] = 0;
  dfkdot_dq[12][0] = z[210];
  dfkdot_dq[12][1] = 0;
  dfkdot_dq[12][2] = 0;
  dfkdot_dq[12][3] = z[211];
  dfkdot_dq[12][4] = z[214];
  dfkdot_dq[12][5] = 0;
  dfkdot_dq[12][6] = 0;
  dfkdot_dq[12][7] = 0;
  dfkdot_dq[12][8] = 0;
  dfkdot_dq[12][9] = 0;
  dfkdot_dq[13][0] = z[215];
  dfkdot_dq[13][1] = 0;
  dfkdot_dq[13][2] = 0;
  dfkdot_dq[13][3] = z[216];
  dfkdot_dq[13][4] = z[217];
  dfkdot_dq[13][5] = 0;
  dfkdot_dq[13][6] = 0;
  dfkdot_dq[13][7] = 0;
  dfkdot_dq[13][8] = 0;
  dfkdot_dq[13][9] = 0;
  dfkdot_dq[14][0] = z[219];
  dfkdot_dq[14][1] = 0;
  dfkdot_dq[14][2] = 0;
  dfkdot_dq[14][3] = z[219];
  dfkdot_dq[14][4] = z[177];
  dfkdot_dq[14][5] = 0;
  dfkdot_dq[14][6] = 0;
  dfkdot_dq[14][7] = 0;
  dfkdot_dq[14][8] = 0;
  dfkdot_dq[14][9] = 0;
  dfkdot_dq[15][0] = z[222];
  dfkdot_dq[15][1] = 0;
  dfkdot_dq[15][2] = 0;
  dfkdot_dq[15][3] = z[222];
  dfkdot_dq[15][4] = z[224];
  dfkdot_dq[15][5] = 0;
  dfkdot_dq[15][6] = 0;
  dfkdot_dq[15][7] = 0;
  dfkdot_dq[15][8] = 0;
  dfkdot_dq[15][9] = 0;
  dfkdot_dq[16][0] = z[176];
  dfkdot_dq[16][1] = 0;
  dfkdot_dq[16][2] = 0;
  dfkdot_dq[16][3] = z[176];
  dfkdot_dq[16][4] = z[179];
  dfkdot_dq[16][5] = 0;
  dfkdot_dq[16][6] = 0;
  dfkdot_dq[16][7] = 0;
  dfkdot_dq[16][8] = 0;
  dfkdot_dq[16][9] = 0;
  dfkdot_dq[17][0] = z[177];
  dfkdot_dq[17][1] = 0;
  dfkdot_dq[17][2] = 0;
  dfkdot_dq[17][3] = z[177];
  dfkdot_dq[17][4] = z[226];
  dfkdot_dq[17][5] = 0;
  dfkdot_dq[17][6] = 0;
  dfkdot_dq[17][7] = 0;
  dfkdot_dq[17][8] = 0;
  dfkdot_dq[17][9] = 0;
  dfkdot_dq[18][0] = z[227];
  dfkdot_dq[18][1] = 0;
  dfkdot_dq[18][2] = 0;
  dfkdot_dq[18][3] = z[228];
  dfkdot_dq[18][4] = z[229];
  dfkdot_dq[18][5] = 0;
  dfkdot_dq[18][6] = 0;
  dfkdot_dq[18][7] = 0;
  dfkdot_dq[18][8] = 0;
  dfkdot_dq[18][9] = 0;
  dfkdot_dq[19][0] = z[230];
  dfkdot_dq[19][1] = 0;
  dfkdot_dq[19][2] = 0;
  dfkdot_dq[19][3] = z[231];
  dfkdot_dq[19][4] = z[232];
  dfkdot_dq[19][5] = 0;
  dfkdot_dq[19][6] = 0;
  dfkdot_dq[19][7] = 0;
  dfkdot_dq[19][8] = 0;
  dfkdot_dq[19][9] = 0;
  dfkdot_dq[20][0] = z[235];
  dfkdot_dq[20][1] = 0;
  dfkdot_dq[20][2] = 0;
  dfkdot_dq[20][3] = z[235];
  dfkdot_dq[20][4] = z[237];
  dfkdot_dq[20][5] = z[183];
  dfkdot_dq[20][6] = 0;
  dfkdot_dq[20][7] = 0;
  dfkdot_dq[20][8] = 0;
  dfkdot_dq[20][9] = 0;
  dfkdot_dq[21][0] = z[241];
  dfkdot_dq[21][1] = 0;
  dfkdot_dq[21][2] = 0;
  dfkdot_dq[21][3] = z[241];
  dfkdot_dq[21][4] = z[244];
  dfkdot_dq[21][5] = z[246];
  dfkdot_dq[21][6] = 0;
  dfkdot_dq[21][7] = 0;
  dfkdot_dq[21][8] = 0;
  dfkdot_dq[21][9] = 0;
  dfkdot_dq[22][0] = z[182];
  dfkdot_dq[22][1] = 0;
  dfkdot_dq[22][2] = 0;
  dfkdot_dq[22][3] = z[182];
  dfkdot_dq[22][4] = z[248];
  dfkdot_dq[22][5] = z[185];
  dfkdot_dq[22][6] = 0;
  dfkdot_dq[22][7] = 0;
  dfkdot_dq[22][8] = 0;
  dfkdot_dq[22][9] = 0;
  dfkdot_dq[23][0] = z[183];
  dfkdot_dq[23][1] = 0;
  dfkdot_dq[23][2] = 0;
  dfkdot_dq[23][3] = z[183];
  dfkdot_dq[23][4] = z[251];
  dfkdot_dq[23][5] = z[253];
  dfkdot_dq[23][6] = 0;
  dfkdot_dq[23][7] = 0;
  dfkdot_dq[23][8] = 0;
  dfkdot_dq[23][9] = 0;
  dfkdot_dq[24][0] = -z[254];
  dfkdot_dq[24][1] = 0;
  dfkdot_dq[24][2] = 0;
  dfkdot_dq[24][3] = 0;
  dfkdot_dq[24][4] = 0;
  dfkdot_dq[24][5] = 0;
  dfkdot_dq[24][6] = 0;
  dfkdot_dq[24][7] = 0;
  dfkdot_dq[24][8] = 0;
  dfkdot_dq[24][9] = 0;
  dfkdot_dq[25][0] = -z[255];
  dfkdot_dq[25][1] = 0;
  dfkdot_dq[25][2] = 0;
  dfkdot_dq[25][3] = 0;
  dfkdot_dq[25][4] = 0;
  dfkdot_dq[25][5] = 0;
  dfkdot_dq[25][6] = 0;
  dfkdot_dq[25][7] = 0;
  dfkdot_dq[25][8] = 0;
  dfkdot_dq[25][9] = 0;
  dfkdot_dq[26][0] = z[189];
  dfkdot_dq[26][1] = 0;
  dfkdot_dq[26][2] = 0;
  dfkdot_dq[26][3] = 0;
  dfkdot_dq[26][4] = 0;
  dfkdot_dq[26][5] = 0;
  dfkdot_dq[26][6] = z[189];
  dfkdot_dq[26][7] = 0;
  dfkdot_dq[26][8] = 0;
  dfkdot_dq[26][9] = 0;
  dfkdot_dq[27][0] = z[303];
  dfkdot_dq[27][1] = 0;
  dfkdot_dq[27][2] = 0;
  dfkdot_dq[27][3] = 0;
  dfkdot_dq[27][4] = 0;
  dfkdot_dq[27][5] = 0;
  dfkdot_dq[27][6] = z[303];
  dfkdot_dq[27][7] = 0;
  dfkdot_dq[27][8] = 0;
  dfkdot_dq[27][9] = 0;
  dfkdot_dq[28][0] = z[188];
  dfkdot_dq[28][1] = 0;
  dfkdot_dq[28][2] = 0;
  dfkdot_dq[28][3] = 0;
  dfkdot_dq[28][4] = 0;
  dfkdot_dq[28][5] = 0;
  dfkdot_dq[28][6] = z[188];
  dfkdot_dq[28][7] = 0;
  dfkdot_dq[28][8] = 0;
  dfkdot_dq[28][9] = 0;
  dfkdot_dq[29][0] = z[189];
  dfkdot_dq[29][1] = 0;
  dfkdot_dq[29][2] = 0;
  dfkdot_dq[29][3] = 0;
  dfkdot_dq[29][4] = 0;
  dfkdot_dq[29][5] = 0;
  dfkdot_dq[29][6] = z[189];
  dfkdot_dq[29][7] = 0;
  dfkdot_dq[29][8] = 0;
  dfkdot_dq[29][9] = 0;
  dfkdot_dq[30][0] = z[256];
  dfkdot_dq[30][1] = 0;
  dfkdot_dq[30][2] = 0;
  dfkdot_dq[30][3] = 0;
  dfkdot_dq[30][4] = 0;
  dfkdot_dq[30][5] = 0;
  dfkdot_dq[30][6] = z[257];
  dfkdot_dq[30][7] = z[260];
  dfkdot_dq[30][8] = 0;
  dfkdot_dq[30][9] = 0;
  dfkdot_dq[31][0] = z[261];
  dfkdot_dq[31][1] = 0;
  dfkdot_dq[31][2] = 0;
  dfkdot_dq[31][3] = 0;
  dfkdot_dq[31][4] = 0;
  dfkdot_dq[31][5] = 0;
  dfkdot_dq[31][6] = z[262];
  dfkdot_dq[31][7] = z[263];
  dfkdot_dq[31][8] = 0;
  dfkdot_dq[31][9] = 0;
  dfkdot_dq[32][0] = z[265];
  dfkdot_dq[32][1] = 0;
  dfkdot_dq[32][2] = 0;
  dfkdot_dq[32][3] = 0;
  dfkdot_dq[32][4] = 0;
  dfkdot_dq[32][5] = 0;
  dfkdot_dq[32][6] = z[265];
  dfkdot_dq[32][7] = z[194];
  dfkdot_dq[32][8] = 0;
  dfkdot_dq[32][9] = 0;
  dfkdot_dq[33][0] = z[268];
  dfkdot_dq[33][1] = 0;
  dfkdot_dq[33][2] = 0;
  dfkdot_dq[33][3] = 0;
  dfkdot_dq[33][4] = 0;
  dfkdot_dq[33][5] = 0;
  dfkdot_dq[33][6] = z[268];
  dfkdot_dq[33][7] = z[270];
  dfkdot_dq[33][8] = 0;
  dfkdot_dq[33][9] = 0;
  dfkdot_dq[34][0] = z[193];
  dfkdot_dq[34][1] = 0;
  dfkdot_dq[34][2] = 0;
  dfkdot_dq[34][3] = 0;
  dfkdot_dq[34][4] = 0;
  dfkdot_dq[34][5] = 0;
  dfkdot_dq[34][6] = z[193];
  dfkdot_dq[34][7] = z[196];
  dfkdot_dq[34][8] = 0;
  dfkdot_dq[34][9] = 0;
  dfkdot_dq[35][0] = z[194];
  dfkdot_dq[35][1] = 0;
  dfkdot_dq[35][2] = 0;
  dfkdot_dq[35][3] = 0;
  dfkdot_dq[35][4] = 0;
  dfkdot_dq[35][5] = 0;
  dfkdot_dq[35][6] = z[194];
  dfkdot_dq[35][7] = z[272];
  dfkdot_dq[35][8] = 0;
  dfkdot_dq[35][9] = 0;
  dfkdot_dq[36][0] = z[273];
  dfkdot_dq[36][1] = 0;
  dfkdot_dq[36][2] = 0;
  dfkdot_dq[36][3] = 0;
  dfkdot_dq[36][4] = 0;
  dfkdot_dq[36][5] = 0;
  dfkdot_dq[36][6] = z[274];
  dfkdot_dq[36][7] = z[275];
  dfkdot_dq[36][8] = 0;
  dfkdot_dq[36][9] = 0;
  dfkdot_dq[37][0] = z[276];
  dfkdot_dq[37][1] = 0;
  dfkdot_dq[37][2] = 0;
  dfkdot_dq[37][3] = 0;
  dfkdot_dq[37][4] = 0;
  dfkdot_dq[37][5] = 0;
  dfkdot_dq[37][6] = z[277];
  dfkdot_dq[37][7] = z[278];
  dfkdot_dq[37][8] = 0;
  dfkdot_dq[37][9] = 0;
  dfkdot_dq[38][0] = z[281];
  dfkdot_dq[38][1] = 0;
  dfkdot_dq[38][2] = 0;
  dfkdot_dq[38][3] = 0;
  dfkdot_dq[38][4] = 0;
  dfkdot_dq[38][5] = 0;
  dfkdot_dq[38][6] = z[281];
  dfkdot_dq[38][7] = z[283];
  dfkdot_dq[38][8] = z[200];
  dfkdot_dq[38][9] = 0;
  dfkdot_dq[39][0] = z[287];
  dfkdot_dq[39][1] = 0;
  dfkdot_dq[39][2] = 0;
  dfkdot_dq[39][3] = 0;
  dfkdot_dq[39][4] = 0;
  dfkdot_dq[39][5] = 0;
  dfkdot_dq[39][6] = z[287];
  dfkdot_dq[39][7] = z[290];
  dfkdot_dq[39][8] = z[292];
  dfkdot_dq[39][9] = 0;
  dfkdot_dq[40][0] = z[199];
  dfkdot_dq[40][1] = 0;
  dfkdot_dq[40][2] = 0;
  dfkdot_dq[40][3] = 0;
  dfkdot_dq[40][4] = 0;
  dfkdot_dq[40][5] = 0;
  dfkdot_dq[40][6] = z[199];
  dfkdot_dq[40][7] = z[294];
  dfkdot_dq[40][8] = z[202];
  dfkdot_dq[40][9] = 0;
  dfkdot_dq[41][0] = z[200];
  dfkdot_dq[41][1] = 0;
  dfkdot_dq[41][2] = 0;
  dfkdot_dq[41][3] = 0;
  dfkdot_dq[41][4] = 0;
  dfkdot_dq[41][5] = 0;
  dfkdot_dq[41][6] = z[200];
  dfkdot_dq[41][7] = z[297];
  dfkdot_dq[41][8] = z[299];
  dfkdot_dq[41][9] = 0;
  dfkdot_dq[42][0] = -z[300];
  dfkdot_dq[42][1] = 0;
  dfkdot_dq[42][2] = 0;
  dfkdot_dq[42][3] = 0;
  dfkdot_dq[42][4] = 0;
  dfkdot_dq[42][5] = 0;
  dfkdot_dq[42][6] = 0;
  dfkdot_dq[42][7] = 0;
  dfkdot_dq[42][8] = 0;
  dfkdot_dq[42][9] = 0;
  dfkdot_dq[43][0] = -z[301];
  dfkdot_dq[43][1] = 0;
  dfkdot_dq[43][2] = 0;
  dfkdot_dq[43][3] = 0;
  dfkdot_dq[43][4] = 0;
  dfkdot_dq[43][5] = 0;
  dfkdot_dq[43][6] = 0;
  dfkdot_dq[43][7] = 0;
  dfkdot_dq[43][8] = 0;
  dfkdot_dq[43][9] = 0;
  dfkdot_dq[44][0] = z[206];
  dfkdot_dq[44][1] = 0;
  dfkdot_dq[44][2] = 0;
  dfkdot_dq[44][3] = 0;
  dfkdot_dq[44][4] = 0;
  dfkdot_dq[44][5] = 0;
  dfkdot_dq[44][6] = 0;
  dfkdot_dq[44][7] = 0;
  dfkdot_dq[44][8] = 0;
  dfkdot_dq[44][9] = z[206];
  dfkdot_dq[45][0] = z[304];
  dfkdot_dq[45][1] = 0;
  dfkdot_dq[45][2] = 0;
  dfkdot_dq[45][3] = 0;
  dfkdot_dq[45][4] = 0;
  dfkdot_dq[45][5] = 0;
  dfkdot_dq[45][6] = 0;
  dfkdot_dq[45][7] = 0;
  dfkdot_dq[45][8] = 0;
  dfkdot_dq[45][9] = z[304];
  dfkdot_dq[46][0] = z[205];
  dfkdot_dq[46][1] = 0;
  dfkdot_dq[46][2] = 0;
  dfkdot_dq[46][3] = 0;
  dfkdot_dq[46][4] = 0;
  dfkdot_dq[46][5] = 0;
  dfkdot_dq[46][6] = 0;
  dfkdot_dq[46][7] = 0;
  dfkdot_dq[46][8] = 0;
  dfkdot_dq[46][9] = z[205];
  dfkdot_dq[47][0] = z[206];
  dfkdot_dq[47][1] = 0;
  dfkdot_dq[47][2] = 0;
  dfkdot_dq[47][3] = 0;
  dfkdot_dq[47][4] = 0;
  dfkdot_dq[47][5] = 0;
  dfkdot_dq[47][6] = 0;
  dfkdot_dq[47][7] = 0;
  dfkdot_dq[47][8] = 0;
  dfkdot_dq[47][9] = z[206];
  f[0] = (z[364]-z[527]-z[513]*q1pp-z[515]*q7pp-z[517]*q8pp-z[518]*q9pp-
  z[520]*q4pp-z[522]*q5pp-z[523]*q6pp-z[524]*q10pp-z[525]*q3pp-z[526]*q2pp)/
  par__bodyweight;
  f[1] = (z[365]-z[537]-z[526]*q1pp-z[528]*q7pp-z[529]*q3pp-z[530]*q2pp-
  z[531]*q9pp-z[532]*q8pp-z[533]*q4pp-z[534]*q6pp-z[535]*q5pp-z[536]*q10pp)/
  par__bodyweight;
  f[2] = (z[366]-z[546]-z[525]*q1pp-z[529]*q2pp-z[538]*q7pp-z[539]*q3pp-
  z[540]*q9pp-z[541]*q8pp-z[542]*q4pp-z[543]*q6pp-z[544]*q5pp-z[545]*q10pp)/
  par__bodyweight;
  f[3] = (z[371]-z[551]-z[520]*q1pp-z[533]*q2pp-z[542]*q3pp-z[548]*q4pp-
  z[549]*q5pp-z[550]*q6pp)/par__bodyweight;
  f[4] = (z[376]-z[555]-z[522]*q1pp-z[535]*q2pp-z[544]*q3pp-z[549]*q4pp-
  z[553]*q5pp-z[554]*q6pp)/par__bodyweight;
  f[5] = (z[381]-z[557]-z[556]*q6pp-z[523]*q1pp-z[534]*q2pp-z[543]*q3pp-
  z[550]*q4pp-z[554]*q5pp)/par__bodyweight;
  f[6] = (z[386]-z[562]-z[515]*q1pp-z[528]*q2pp-z[538]*q3pp-z[559]*q7pp-
  z[560]*q8pp-z[561]*q9pp)/par__bodyweight;
  f[7] = (z[391]-z[566]-z[517]*q1pp-z[532]*q2pp-z[541]*q3pp-z[560]*q7pp-
  z[564]*q8pp-z[565]*q9pp)/par__bodyweight;
  f[8] = (z[396]-z[568]-z[567]*q9pp-z[518]*q1pp-z[531]*q2pp-z[540]*q3pp-
  z[561]*q7pp-z[565]*q8pp)/par__bodyweight;
  f[9] = (z[401]-z[570]-z[569]*q10pp-z[524]*q1pp-z[536]*q2pp-z[545]*q3pp)/
  par__bodyweight;
  df_dq[0][0] = z[960];
  df_dq[0][1] = -z[961];
  df_dq[0][2] = z[962];
  df_dq[0][3] = z[963];
  df_dq[0][4] = z[964];
  df_dq[0][5] = z[965];
  df_dq[0][6] = z[966];
  df_dq[0][7] = z[967];
  df_dq[0][8] = z[968];
  df_dq[0][9] = z[969];
  df_dq[1][0] = z[970];
  df_dq[1][1] = 0;
  df_dq[1][2] = 0;
  df_dq[1][3] = z[971];
  df_dq[1][4] = z[972];
  df_dq[1][5] = z[973];
  df_dq[1][6] = z[974];
  df_dq[1][7] = z[975];
  df_dq[1][8] = z[976];
  df_dq[1][9] = z[977];
  df_dq[2][0] = z[978];
  df_dq[2][1] = 0;
  df_dq[2][2] = 0;
  df_dq[2][3] = z[979];
  df_dq[2][4] = z[980];
  df_dq[2][5] = z[981];
  df_dq[2][6] = z[982];
  df_dq[2][7] = z[983];
  df_dq[2][8] = z[984];
  df_dq[2][9] = z[985];
  df_dq[3][0] = z[986];
  df_dq[3][1] = -z[987];
  df_dq[3][2] = z[988];
  df_dq[3][3] = z[989];
  df_dq[3][4] = z[990];
  df_dq[3][5] = z[991];
  df_dq[3][6] = 0;
  df_dq[3][7] = 0;
  df_dq[3][8] = 0;
  df_dq[3][9] = 0;
  df_dq[4][0] = z[992];
  df_dq[4][1] = -z[987];
  df_dq[4][2] = z[988];
  df_dq[4][3] = z[993];
  df_dq[4][4] = z[994];
  df_dq[4][5] = z[995];
  df_dq[4][6] = 0;
  df_dq[4][7] = 0;
  df_dq[4][8] = 0;
  df_dq[4][9] = 0;
  df_dq[5][0] = z[996];
  df_dq[5][1] = -z[987];
  df_dq[5][2] = z[988];
  df_dq[5][3] = z[997];
  df_dq[5][4] = z[998];
  df_dq[5][5] = z[999];
  df_dq[5][6] = 0;
  df_dq[5][7] = 0;
  df_dq[5][8] = 0;
  df_dq[5][9] = 0;
  df_dq[6][0] = z[1000];
  df_dq[6][1] = -z[1001];
  df_dq[6][2] = z[1002];
  df_dq[6][3] = 0;
  df_dq[6][4] = 0;
  df_dq[6][5] = 0;
  df_dq[6][6] = z[1003];
  df_dq[6][7] = z[1004];
  df_dq[6][8] = z[1005];
  df_dq[6][9] = 0;
  df_dq[7][0] = z[1006];
  df_dq[7][1] = -z[1001];
  df_dq[7][2] = z[1002];
  df_dq[7][3] = 0;
  df_dq[7][4] = 0;
  df_dq[7][5] = 0;
  df_dq[7][6] = z[1007];
  df_dq[7][7] = z[1008];
  df_dq[7][8] = z[1009];
  df_dq[7][9] = 0;
  df_dq[8][0] = z[1010];
  df_dq[8][1] = -z[1001];
  df_dq[8][2] = z[1002];
  df_dq[8][3] = 0;
  df_dq[8][4] = 0;
  df_dq[8][5] = 0;
  df_dq[8][6] = z[1011];
  df_dq[8][7] = z[1012];
  df_dq[8][8] = z[1013];
  df_dq[8][9] = 0;
  df_dq[9][0] = z[1014];
  df_dq[9][1] = 0;
  df_dq[9][2] = 0;
  df_dq[9][3] = 0;
  df_dq[9][4] = 0;
  df_dq[9][5] = 0;
  df_dq[9][6] = 0;
  df_dq[9][7] = 0;
  df_dq[9][8] = 0;
  df_dq[9][9] = z[1015];
  df_dqd[0][0] = z[1198];
  df_dqd[0][1] = z[1199];
  df_dqd[0][2] = z[1200];
  df_dqd[0][3] = -z[1201];
  df_dqd[0][4] = -z[1202];
  df_dqd[0][5] = -z[1203];
  df_dqd[0][6] = -z[1204];
  df_dqd[0][7] = -z[1205];
  df_dqd[0][8] = -z[1206];
  df_dqd[0][9] = -z[1207];
  df_dqd[1][0] = z[1208];
  df_dqd[1][1] = z[1209];
  df_dqd[1][2] = z[1210];
  df_dqd[1][3] = -z[1211];
  df_dqd[1][4] = -z[1212];
  df_dqd[1][5] = -z[1213];
  df_dqd[1][6] = -z[1214];
  df_dqd[1][7] = -z[1215];
  df_dqd[1][8] = -z[1216];
  df_dqd[1][9] = -z[1217];
  df_dqd[2][0] = z[1218];
  df_dqd[2][1] = z[1219];
  df_dqd[2][2] = z[1220];
  df_dqd[2][3] = -z[1221];
  df_dqd[2][4] = -z[1222];
  df_dqd[2][5] = -z[1223];
  df_dqd[2][6] = -z[1224];
  df_dqd[2][7] = -z[1225];
  df_dqd[2][8] = -z[1226];
  df_dqd[2][9] = -z[1227];
  df_dqd[3][0] = -z[1228];
  df_dqd[3][1] = -z[1229];
  df_dqd[3][2] = -z[1230];
  df_dqd[3][3] = -z[1231];
  df_dqd[3][4] = -z[1232];
  df_dqd[3][5] = -z[1233];
  df_dqd[3][6] = 0;
  df_dqd[3][7] = 0;
  df_dqd[3][8] = 0;
  df_dqd[3][9] = 0;
  df_dqd[4][0] = -z[1234];
  df_dqd[4][1] = -z[1235];
  df_dqd[4][2] = -z[1236];
  df_dqd[4][3] = -z[1237];
  df_dqd[4][4] = -z[1238];
  df_dqd[4][5] = -z[1239];
  df_dqd[4][6] = 0;
  df_dqd[4][7] = 0;
  df_dqd[4][8] = 0;
  df_dqd[4][9] = 0;
  df_dqd[5][0] = -z[1240];
  df_dqd[5][1] = -z[1241];
  df_dqd[5][2] = -z[1242];
  df_dqd[5][3] = -z[1243];
  df_dqd[5][4] = -z[1244];
  df_dqd[5][5] = -z[1245];
  df_dqd[5][6] = 0;
  df_dqd[5][7] = 0;
  df_dqd[5][8] = 0;
  df_dqd[5][9] = 0;
  df_dqd[6][0] = -z[1246];
  df_dqd[6][1] = -z[1247];
  df_dqd[6][2] = -z[1248];
  df_dqd[6][3] = 0;
  df_dqd[6][4] = 0;
  df_dqd[6][5] = 0;
  df_dqd[6][6] = -z[1249];
  df_dqd[6][7] = -z[1250];
  df_dqd[6][8] = -z[1251];
  df_dqd[6][9] = 0;
  df_dqd[7][0] = -z[1252];
  df_dqd[7][1] = -z[1253];
  df_dqd[7][2] = -z[1254];
  df_dqd[7][3] = 0;
  df_dqd[7][4] = 0;
  df_dqd[7][5] = 0;
  df_dqd[7][6] = -z[1255];
  df_dqd[7][7] = -z[1256];
  df_dqd[7][8] = -z[1257];
  df_dqd[7][9] = 0;
  df_dqd[8][0] = -z[1258];
  df_dqd[8][1] = -z[1259];
  df_dqd[8][2] = -z[1260];
  df_dqd[8][3] = 0;
  df_dqd[8][4] = 0;
  df_dqd[8][5] = 0;
  df_dqd[8][6] = -z[1261];
  df_dqd[8][7] = -z[1262];
  df_dqd[8][8] = -z[1263];
  df_dqd[8][9] = 0;
  df_dqd[9][0] = -z[1264];
  df_dqd[9][1] = -z[1265];
  df_dqd[9][2] = -z[1266];
  df_dqd[9][3] = 0;
  df_dqd[9][4] = 0;
  df_dqd[9][5] = 0;
  df_dqd[9][6] = 0;
  df_dqd[9][7] = 0;
  df_dqd[9][8] = 0;
  df_dqd[9][9] = -z[1267];
  df_dqdd[0][0] = -z[1268];
  df_dqdd[0][1] = -z[1269];
  df_dqdd[0][2] = -z[1270];
  df_dqdd[0][3] = -z[1271];
  df_dqdd[0][4] = -z[1272];
  df_dqdd[0][5] = -z[1273];
  df_dqdd[0][6] = -z[1274];
  df_dqdd[0][7] = -z[1275];
  df_dqdd[0][8] = -z[1276];
  df_dqdd[0][9] = -z[1277];
  df_dqdd[1][0] = -z[1269];
  df_dqdd[1][1] = -z[1278];
  df_dqdd[1][2] = -z[1279];
  df_dqdd[1][3] = -z[1280];
  df_dqdd[1][4] = -z[1281];
  df_dqdd[1][5] = -z[1282];
  df_dqdd[1][6] = -z[1283];
  df_dqdd[1][7] = -z[1284];
  df_dqdd[1][8] = -z[1285];
  df_dqdd[1][9] = -z[1286];
  df_dqdd[2][0] = -z[1270];
  df_dqdd[2][1] = -z[1279];
  df_dqdd[2][2] = -z[1287];
  df_dqdd[2][3] = -z[1288];
  df_dqdd[2][4] = -z[1289];
  df_dqdd[2][5] = -z[1290];
  df_dqdd[2][6] = -z[1291];
  df_dqdd[2][7] = -z[1292];
  df_dqdd[2][8] = -z[1293];
  df_dqdd[2][9] = -z[1294];
  df_dqdd[3][0] = -z[1271];
  df_dqdd[3][1] = -z[1280];
  df_dqdd[3][2] = -z[1288];
  df_dqdd[3][3] = -z[1295];
  df_dqdd[3][4] = -z[1296];
  df_dqdd[3][5] = -z[1297];
  df_dqdd[3][6] = 0;
  df_dqdd[3][7] = 0;
  df_dqdd[3][8] = 0;
  df_dqdd[3][9] = 0;
  df_dqdd[4][0] = -z[1272];
  df_dqdd[4][1] = -z[1281];
  df_dqdd[4][2] = -z[1289];
  df_dqdd[4][3] = -z[1296];
  df_dqdd[4][4] = -z[1298];
  df_dqdd[4][5] = -z[1299];
  df_dqdd[4][6] = 0;
  df_dqdd[4][7] = 0;
  df_dqdd[4][8] = 0;
  df_dqdd[4][9] = 0;
  df_dqdd[5][0] = -z[1273];
  df_dqdd[5][1] = -z[1282];
  df_dqdd[5][2] = -z[1290];
  df_dqdd[5][3] = -z[1297];
  df_dqdd[5][4] = -z[1299];
  df_dqdd[5][5] = -z[1300];
  df_dqdd[5][6] = 0;
  df_dqdd[5][7] = 0;
  df_dqdd[5][8] = 0;
  df_dqdd[5][9] = 0;
  df_dqdd[6][0] = -z[1274];
  df_dqdd[6][1] = -z[1283];
  df_dqdd[6][2] = -z[1291];
  df_dqdd[6][3] = 0;
  df_dqdd[6][4] = 0;
  df_dqdd[6][5] = 0;
  df_dqdd[6][6] = -z[1301];
  df_dqdd[6][7] = -z[1302];
  df_dqdd[6][8] = -z[1303];
  df_dqdd[6][9] = 0;
  df_dqdd[7][0] = -z[1275];
  df_dqdd[7][1] = -z[1284];
  df_dqdd[7][2] = -z[1292];
  df_dqdd[7][3] = 0;
  df_dqdd[7][4] = 0;
  df_dqdd[7][5] = 0;
  df_dqdd[7][6] = -z[1302];
  df_dqdd[7][7] = -z[1304];
  df_dqdd[7][8] = -z[1305];
  df_dqdd[7][9] = 0;
  df_dqdd[8][0] = -z[1276];
  df_dqdd[8][1] = -z[1285];
  df_dqdd[8][2] = -z[1293];
  df_dqdd[8][3] = 0;
  df_dqdd[8][4] = 0;
  df_dqdd[8][5] = 0;
  df_dqdd[8][6] = -z[1303];
  df_dqdd[8][7] = -z[1305];
  df_dqdd[8][8] = -z[1306];
  df_dqdd[8][9] = 0;
  df_dqdd[9][0] = -z[1277];
  df_dqdd[9][1] = -z[1286];
  df_dqdd[9][2] = -z[1294];
  df_dqdd[9][3] = 0;
  df_dqdd[9][4] = 0;
  df_dqdd[9][5] = 0;
  df_dqdd[9][6] = 0;
  df_dqdd[9][7] = 0;
  df_dqdd[9][8] = 0;
  df_dqdd[9][9] = -z[1307];
  df_dG[0][0] = -z[1326];
  df_dG[0][1] = -z[1327];
  df_dG[0][2] = 1;
  df_dG[0][3] = -z[1328];
  df_dG[0][4] = -z[1329];
  df_dG[0][5] = 1;
  df_dG[1][0] = z[1330];
  df_dG[1][1] = z[1331];
  df_dG[1][2] = 0;
  df_dG[1][3] = z[1332];
  df_dG[1][4] = z[1333];
  df_dG[1][5] = 0;
  df_dG[2][0] = z[1331];
  df_dG[2][1] = z[1334];
  df_dG[2][2] = 0;
  df_dG[2][3] = z[1333];
  df_dG[2][4] = z[1335];
  df_dG[2][5] = 0;
  df_dG[3][0] = -z[1336];
  df_dG[3][1] = -z[1337];
  df_dG[3][2] = 1;
  df_dG[3][3] = 0;
  df_dG[3][4] = 0;
  df_dG[3][5] = 0;
  df_dG[4][0] = -z[1338];
  df_dG[4][1] = -z[1339];
  df_dG[4][2] = 1;
  df_dG[4][3] = 0;
  df_dG[4][4] = 0;
  df_dG[4][5] = 0;
  df_dG[5][0] = z[1340];
  df_dG[5][1] = -z[1341];
  df_dG[5][2] = 1;
  df_dG[5][3] = 0;
  df_dG[5][4] = 0;
  df_dG[5][5] = 0;
  df_dG[6][0] = 0;
  df_dG[6][1] = 0;
  df_dG[6][2] = 0;
  df_dG[6][3] = -z[1342];
  df_dG[6][4] = -z[1343];
  df_dG[6][5] = 1;
  df_dG[7][0] = 0;
  df_dG[7][1] = 0;
  df_dG[7][2] = 0;
  df_dG[7][3] = -z[1344];
  df_dG[7][4] = -z[1345];
  df_dG[7][5] = 1;
  df_dG[8][0] = 0;
  df_dG[8][1] = 0;
  df_dG[8][2] = 0;
  df_dG[8][3] = z[1346];
  df_dG[8][4] = -z[1347];
  df_dG[8][5] = 1;
  df_dG[9][0] = 0;
  df_dG[9][1] = 0;
  df_dG[9][2] = 0;
  df_dG[9][3] = 0;
  df_dG[9][4] = 0;
  df_dG[9][5] = 0;

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


