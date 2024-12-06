/* The name of this program is gait2d_osim/gait2d_osim_FK_raw.c */
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

double   par__back_x,par__back_y,par__Lankle_y,par__Lhip_x,par__Lhip_y,par__Lknee_x1,
  par__Lknee_x2,par__Lknee_x3,par__Lknee_x4,par__Lknee_x5,par__Lknee_y1,par__Lknee_y2,
  par__Lknee_y3,par__Lknee_y4,par__Lknee_y5,par__Rankle_y,par__Rhip_x,par__Rhip_y,
  par__Rknee_x1,par__Rknee_x2,par__Rknee_x3,par__Rknee_x4,par__Rknee_x5,par__Rknee_y1,
  par__Rknee_y2,par__Rknee_y3,par__Rknee_y4,par__Rknee_y5,q1,q2,q3,q4,q5,q6,q7,
  q8,q9,q1p,q2p,q3p,q4p,q5p,q6p,q7p,q8p,q9p;
double   Pi,DEGtoRAD,RADtoDEG,z[293],fk[48],dfk_dq[48][9],fkdot[48],dfkdot_dq[
  48][9],Encode[4];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait2d_osim/gait2d_osim_FK_raw.in");
  else sprintf(fileName, "gait2d_osim/gait2d_osim_FK_raw.%d", iloop);
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
  q2,&q3,&q4,&q5,&q6,&q7,&q8,&q9,&q1p,&q2p,&q3p,&q4p,&q5p,&q6p,&q7p,&q8p,&q9p,
  NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file gait2d_osim/gait2d_osim_FK_raw.in" );
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


