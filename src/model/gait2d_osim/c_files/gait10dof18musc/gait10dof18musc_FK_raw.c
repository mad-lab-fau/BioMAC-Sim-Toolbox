/* The name of this program is gait10dof18musc/gait10dof18musc_FK_raw.c */
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

double   par__back_x,par__back_y,par__Lankle_y,par__Lhip_x,par__Lhip_y,par__Lknee_x1,
  par__Lknee_x2,par__Lknee_x3,par__Lknee_x4,par__Lknee_x5,par__Lknee_y1,par__Lknee_y2,
  par__Lknee_y3,par__Lknee_y4,par__Lknee_y5,par__Rankle_y,par__Rhip_x,par__Rhip_y,
  par__Rknee_x1,par__Rknee_x2,par__Rknee_x3,par__Rknee_x4,par__Rknee_x5,par__Rknee_y1,
  par__Rknee_y2,par__Rknee_y3,par__Rknee_y4,par__Rknee_y5,q1,q10,q2,q3,q4,q5,q6,
  q7,q8,q9,q1p,q10p,q2p,q3p,q4p,q5p,q6p,q7p,q8p,q9p;
double   Pi,DEGtoRAD,RADtoDEG,z[305],fk[48],dfk_dq[48][10],fkdot[48],dfkdot_dq[
  48][10],Encode[4];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait10dof18musc/gait10dof18musc_FK_raw.in");
  else sprintf(fileName, "gait10dof18musc/gait10dof18musc_FK_raw.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&par__back_x,&par__back_y,&par__Lankle_y,&par__Lhip_x,&par__Lhip_y,&
  par__Lknee_x1,&par__Lknee_x2,&par__Lknee_x3,&par__Lknee_x4,&par__Lknee_x5,&par__Lknee_y1,&
  par__Lknee_y2,&par__Lknee_y3,&par__Lknee_y4,&par__Lknee_y5,&par__Rankle_y,&par__Rhip_x,&
  par__Rhip_y,&par__Rknee_x1,&par__Rknee_x2,&par__Rknee_x3,&par__Rknee_x4,&par__Rknee_x5,&
  par__Rknee_y1,&par__Rknee_y2,&par__Rknee_y3,&par__Rknee_y4,&par__Rknee_y5,&q1,&
  q10,&q2,&q3,&q4,&q5,&q6,&q7,&q8,&q9,&q1p,&q10p,&q2p,&q3p,&q4p,&q5p,&q6p,&q7p,&
  q8p,&q9p,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file gait10dof18musc/gait10dof18musc_FK_raw.in" );
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
  z[12] = par__Rknee_x5 + par__Rknee_x4*q5 + par__Rknee_x1*pow(q5,4) + par__Rknee_x2*
  pow(q5,3) + par__Rknee_x3*pow(q5,2);
  z[13] = par__Rknee_y5 + par__Rknee_y4*q5 + par__Rknee_y1*pow(q5,4) + par__Rknee_y2*
  pow(q5,3) + par__Rknee_y3*pow(q5,2);
  z[16] = cos(q5);
  z[17] = sin(q5);
  z[18] = z[5]*z[16] + z[6]*z[17];
  z[19] = z[6]*z[16] - z[5]*z[17];
  z[20] = z[5]*z[17] + z[7]*z[16];
  z[21] = z[5]*z[16] - z[7]*z[17];
  z[32] = cos(q6);
  z[33] = sin(q6);
  z[34] = z[18]*z[32] + z[19]*z[33];
  z[35] = z[19]*z[32] - z[18]*z[33];
  z[36] = z[20]*z[32] + z[21]*z[33];
  z[37] = z[21]*z[32] - z[20]*z[33];
  z[50] = cos(q7);
  z[51] = sin(q7);
  z[52] = z[1]*z[50] - z[2]*z[51];
  z[53] = -z[1]*z[51] - z[2]*z[50];
  z[54] = z[1]*z[51] + z[2]*z[50];
  z[59] = par__Lknee_x5 + par__Lknee_x4*q8 + par__Lknee_x1*pow(q8,4) + par__Lknee_x2*
  pow(q8,3) + par__Lknee_x3*pow(q8,2);
  z[60] = par__Lknee_y5 + par__Lknee_y4*q8 + par__Lknee_y1*pow(q8,4) + par__Lknee_y2*
  pow(q8,3) + par__Lknee_y3*pow(q8,2);
  z[63] = cos(q8);
  z[64] = sin(q8);
  z[65] = z[52]*z[63] + z[53]*z[64];
  z[66] = z[53]*z[63] - z[52]*z[64];
  z[67] = z[52]*z[64] + z[54]*z[63];
  z[68] = z[52]*z[63] - z[54]*z[64];
  z[79] = cos(q9);
  z[80] = sin(q9);
  z[81] = z[65]*z[79] + z[66]*z[80];
  z[82] = z[66]*z[79] - z[65]*z[80];
  z[83] = z[67]*z[79] + z[68]*z[80];
  z[84] = z[68]*z[79] - z[67]*z[80];
  z[97] = cos(q10);
  z[98] = sin(q10);
  z[99] = z[1]*z[97] - z[2]*z[98];
  z[100] = -z[1]*z[98] - z[2]*z[97];
  z[101] = z[1]*z[98] + z[2]*z[97];
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

  Encode[0] = 0.0;
  Encode[1] = 0.0;
  Encode[2] = 0.0;
  Encode[3] = 0.0;

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


