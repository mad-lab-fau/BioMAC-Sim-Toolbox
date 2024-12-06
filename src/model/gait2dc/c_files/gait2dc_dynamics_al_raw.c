/* The name of this program is gait2dc_dynamics_al_raw.c */
/* Created by Autolev 4.1 on Thu Oct 24 13:28:07 2024 */

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

double   par__airdrag,par__bodyweight,par__FootCMx,par__FootCMy,par__FootInertia,
  par__FootMass,par__gravity,par__ShankCMy,par__ShankInertia,par__ShankLen,par__ShankMass,
  par__slope,par__ThighCMy,par__ThighInertia,par__ThighLen,par__ThighMass,par__TrunkCMy,
  par__TrunkInertia,par__TrunkMass,par__wind,LFootFx,LFootFy,LFootMz,mom1,mom2,
  mom3,mom4,mom5,mom6,q1,q2,q3,q4,q5,q6,q7,q8,q9,RFootFx,RFootFy,RFootMz,q1p,q2p,
  q3p,q4p,q5p,q6p,q7p,q8p,q9p,q1pp,q2pp,q3pp,q4pp,q5pp,q6pp,q7pp,q8pp,q9pp;
double   Pi,DEGtoRAD,RADtoDEG,z[634],fk[42],dfk_dq[42][9],fkdot[42],dfkdot_dq[
  42][9],Zero[9],dz_dq[9][9],dz_dqd[9][9],dz_dqdd[9][9],dz_dmom[9][6],dz_dgrf[
  9][6],Encode[10];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "gait2dc_dynamics_al_raw.in");
  else sprintf(fileName, "gait2dc_dynamics_al_raw.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&par__airdrag,&par__bodyweight,&par__FootCMx,&par__FootCMy,&par__FootInertia,&
  par__FootMass,&par__gravity,&par__ShankCMy,&par__ShankInertia,&par__ShankLen,&
  par__ShankMass,&par__slope,&par__ThighCMy,&par__ThighInertia,&par__ThighLen,&
  par__ThighMass,&par__TrunkCMy,&par__TrunkInertia,&par__TrunkMass,&par__wind,&
  LFootFx,&LFootFy,&LFootMz,&mom1,&mom2,&mom3,&mom4,&mom5,&mom6,&q1,&q2,&q3,&q4,&
  q5,&q6,&q7,&q8,&q9,&RFootFx,&RFootFy,&RFootMz,&q1p,&q2p,&q3p,&q4p,&q5p,&q6p,&
  q7p,&q8p,&q9p,&q1pp,&q2pp,&q3pp,&q4pp,&q5pp,&q6pp,&q7pp,&q8pp,&q9pp,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file gait2dc_dynamics_al_raw.in" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (void)
{

/* Evaluate constants */
  z[1] = cos(q3);
  z[2] = sin(q3);
  z[3] = par__TrunkCMy*q3p;
  z[4] = q3p*z[3];
  z[5] = cos(q4);
  z[6] = sin(q4);
  z[7] = z[1]*z[5] - z[2]*z[6];
  z[8] = -z[1]*z[6] - z[2]*z[5];
  z[9] = z[1]*z[6] + z[2]*z[5];
  z[10] = par__ThighCMy*(q3p+q4p);
  z[11] = q3p + q4p;
  z[12] = z[10]*z[11];
  z[13] = cos(q5);
  z[14] = sin(q5);
  z[15] = par__ThighLen*(q3p+q4p);
  z[16] = z[11]*z[15];
  z[17] = z[7]*z[13] + z[8]*z[14];
  z[18] = z[8]*z[13] - z[7]*z[14];
  z[19] = z[7]*z[14] + z[9]*z[13];
  z[20] = z[7]*z[13] - z[9]*z[14];
  z[21] = par__ThighLen*z[13];
  z[22] = par__ThighLen*z[14];
  z[23] = z[21] - par__ShankCMy;
  z[24] = z[14]*z[16];
  z[25] = z[13]*z[16];
  z[26] = par__ShankCMy*(q3p+q4p+q5p);
  z[27] = q3p + q4p + q5p;
  z[28] = z[25] - z[26]*z[27];
  z[29] = cos(q6);
  z[30] = sin(q6);
  z[31] = par__ShankLen + z[21];
  z[32] = par__ShankLen*(q3p+q4p+q5p);
  z[33] = z[25] + z[27]*z[32];
  z[34] = z[17]*z[29] + z[18]*z[30];
  z[35] = z[18]*z[29] - z[17]*z[30];
  z[36] = z[19]*z[29] + z[20]*z[30];
  z[37] = z[20]*z[29] - z[19]*z[30];
  z[38] = par__ShankLen*z[29];
  z[39] = z[29]*z[31] - z[22]*z[30];
  z[40] = -z[22]*z[29] - z[30]*z[31];
  z[41] = par__ShankLen*z[30];
  z[42] = z[38] - par__FootCMy;
  z[43] = z[39] - par__FootCMy;
  z[44] = par__FootCMx + z[40];
  z[45] = par__FootCMx - z[41];
  z[46] = z[29]*z[24] + z[30]*z[33];
  z[47] = z[29]*z[33] - z[30]*z[24];
  z[48] = par__FootCMx*(q3p+q4p+q5p+q6p);
  z[49] = par__FootCMy*(q3p+q4p+q5p+q6p);
  z[50] = q3p + q4p + q5p + q6p;
  z[51] = z[46] - z[48]*z[50];
  z[52] = z[47] - z[49]*z[50];
  z[53] = cos(q7);
  z[54] = sin(q7);
  z[55] = z[1]*z[53] - z[2]*z[54];
  z[56] = -z[1]*z[54] - z[2]*z[53];
  z[57] = z[1]*z[54] + z[2]*z[53];
  z[58] = par__ThighCMy*(q3p+q7p);
  z[59] = q3p + q7p;
  z[60] = z[58]*z[59];
  z[61] = cos(q8);
  z[62] = sin(q8);
  z[63] = par__ThighLen*(q3p+q7p);
  z[64] = z[59]*z[63];
  z[65] = z[55]*z[61] + z[56]*z[62];
  z[66] = z[56]*z[61] - z[55]*z[62];
  z[67] = z[55]*z[62] + z[57]*z[61];
  z[68] = z[55]*z[61] - z[57]*z[62];
  z[69] = par__ThighLen*z[61];
  z[70] = par__ThighLen*z[62];
  z[71] = z[69] - par__ShankCMy;
  z[72] = z[62]*z[64];
  z[73] = z[61]*z[64];
  z[74] = par__ShankCMy*(q3p+q7p+q8p);
  z[75] = q3p + q7p + q8p;
  z[76] = z[73] - z[74]*z[75];
  z[77] = cos(q9);
  z[78] = sin(q9);
  z[79] = par__ShankLen + z[69];
  z[80] = par__ShankLen*(q3p+q7p+q8p);
  z[81] = z[73] + z[75]*z[80];
  z[82] = z[65]*z[77] + z[66]*z[78];
  z[83] = z[66]*z[77] - z[65]*z[78];
  z[84] = z[67]*z[77] + z[68]*z[78];
  z[85] = z[68]*z[77] - z[67]*z[78];
  z[86] = par__ShankLen*z[77];
  z[87] = z[77]*z[79] - z[70]*z[78];
  z[88] = -z[70]*z[77] - z[78]*z[79];
  z[89] = par__ShankLen*z[78];
  z[90] = z[86] - par__FootCMy;
  z[91] = z[87] - par__FootCMy;
  z[92] = par__FootCMx + z[88];
  z[93] = par__FootCMx - z[89];
  z[94] = z[77]*z[72] + z[78]*z[81];
  z[95] = z[77]*z[81] - z[78]*z[72];
  z[96] = par__FootCMx*(q3p+q7p+q8p+q9p);
  z[97] = par__FootCMy*(q3p+q7p+q8p+q9p);
  z[98] = q3p + q7p + q8p + q9p;
  z[99] = z[94] - z[96]*z[98];
  z[100] = z[95] - z[97]*z[98];
  z[101] = z[2]*z[6] - z[1]*z[5];
  z[102] = z[8]*z[13] + z[14]*z[101];
  z[103] = z[13]*z[101] - z[8]*z[14];
  z[104] = -z[7]*z[13] - z[8]*z[14];
  z[105] = -z[7]*z[14] - z[9]*z[13];
  z[106] = z[29]*z[102] + z[30]*z[103];
  z[107] = z[18]*z[29] + z[30]*z[104];
  z[108] = z[29]*z[103] - z[30]*z[102];
  z[109] = z[29]*z[104] - z[18]*z[30];
  z[110] = -z[17]*z[29] - z[18]*z[30];
  z[111] = z[20]*z[29] + z[30]*z[105];
  z[112] = z[29]*z[105] - z[20]*z[30];
  z[113] = -z[19]*z[29] - z[20]*z[30];
  z[114] = z[2]*z[54] - z[1]*z[53];
  z[115] = z[56]*z[61] + z[62]*z[114];
  z[116] = z[61]*z[114] - z[56]*z[62];
  z[117] = -z[55]*z[61] - z[56]*z[62];
  z[118] = -z[55]*z[62] - z[57]*z[61];
  z[119] = z[77]*z[115] + z[78]*z[116];
  z[120] = z[66]*z[77] + z[78]*z[117];
  z[121] = z[77]*z[116] - z[78]*z[115];
  z[122] = z[77]*z[117] - z[66]*z[78];
  z[123] = -z[65]*z[77] - z[66]*z[78];
  z[124] = z[68]*z[77] + z[78]*z[118];
  z[125] = z[77]*z[118] - z[68]*z[78];
  z[126] = -z[67]*z[77] - z[68]*z[78];
  z[127] = par__ThighLen*z[101];
  z[128] = par__ThighLen*z[8];
  z[129] = -par__ShankLen*z[103] - par__ThighLen*z[101];
  z[130] = par__ShankLen*z[104];
  z[131] = -par__ShankLen*z[18] - par__ThighLen*z[8];
  z[132] = par__ShankLen*z[105];
  z[133] = par__ThighLen*z[114];
  z[134] = par__ThighLen*z[56];
  z[135] = -par__ShankLen*z[116] - par__ThighLen*z[114];
  z[136] = par__ShankLen*z[117];
  z[137] = -par__ShankLen*z[66] - par__ThighLen*z[56];
  z[138] = par__ShankLen*z[118];
  z[139] = z[2]*q3p;
  z[140] = z[1]*q3p;
  z[141] = z[8]*(q3p+q4p);
  z[142] = z[101]*(q3p+q4p);
  z[143] = z[7]*(q3p+q4p);
  z[144] = q1p - par__ThighLen*z[101]*(q3p+q4p);
  z[145] = q2p - par__ThighLen*z[8]*(q3p+q4p);
  z[146] = z[18]*q5p + z[102]*q3p + z[102]*q4p;
  z[147] = z[103]*q3p + z[103]*q4p + z[104]*q5p;
  z[148] = z[17]*q3p + z[17]*q4p + z[20]*q5p;
  z[149] = z[18]*q3p + z[18]*q4p + z[105]*q5p;
  z[150] = q1p - par__ThighLen*z[101]*(q3p+q4p) - par__ShankLen*(z[103]*q3p+
  z[103]*q4p+z[104]*q5p);
  z[151] = q2p - par__ThighLen*z[8]*(q3p+q4p) - par__ShankLen*(z[18]*q3p+
  z[18]*q4p+z[105]*q5p);
  z[152] = z[35]*q6p + z[106]*q3p + z[106]*q4p + z[107]*q5p;
  z[153] = z[108]*q3p + z[108]*q4p + z[109]*q5p + z[110]*q6p;
  z[154] = z[34]*q3p + z[34]*q4p + z[37]*q6p + z[111]*q5p;
  z[155] = z[35]*q3p + z[35]*q4p + z[112]*q5p + z[113]*q6p;
  z[156] = z[56]*(q3p+q7p);
  z[157] = z[114]*(q3p+q7p);
  z[158] = z[55]*(q3p+q7p);
  z[159] = q1p - par__ThighLen*z[114]*(q3p+q7p);
  z[160] = q2p - par__ThighLen*z[56]*(q3p+q7p);
  z[161] = z[66]*q8p + z[115]*q3p + z[115]*q7p;
  z[162] = z[116]*q3p + z[116]*q7p + z[117]*q8p;
  z[163] = z[65]*q3p + z[65]*q7p + z[68]*q8p;
  z[164] = z[66]*q3p + z[66]*q7p + z[118]*q8p;
  z[165] = q1p - par__ThighLen*z[114]*(q3p+q7p) - par__ShankLen*(z[116]*q3p+
  z[116]*q7p+z[117]*q8p);
  z[166] = q2p - par__ThighLen*z[56]*(q3p+q7p) - par__ShankLen*(z[66]*q3p+
  z[66]*q7p+z[118]*q8p);
  z[167] = z[83]*q9p + z[119]*q3p + z[119]*q7p + z[120]*q8p;
  z[168] = z[121]*q3p + z[121]*q7p + z[122]*q8p + z[123]*q9p;
  z[169] = z[82]*q3p + z[82]*q7p + z[85]*q9p + z[124]*q8p;
  z[170] = z[83]*q3p + z[83]*q7p + z[125]*q8p + z[126]*q9p;
  z[171] = par__ThighLen*z[9]*(q3p+q4p);
  z[172] = par__ThighLen*z[101]*(q3p+q4p);
  z[173] = z[9]*z[14] + z[13]*z[101];
  z[174] = z[103]*q5p + z[173]*q3p + z[173]*q4p;
  z[175] = z[9]*z[13] - z[14]*z[101];
  z[176] = -z[8]*z[13] - z[14]*z[101];
  z[177] = z[175]*q3p + z[175]*q4p + z[176]*q5p;
  z[178] = z[7]*z[14] - z[8]*z[13];
  z[179] = z[176]*q3p + z[176]*q4p + z[178]*q5p;
  z[180] = z[9]*z[14] - z[7]*z[13];
  z[181] = z[104]*q3p + z[104]*q4p + z[180]*q5p;
  z[182] = -par__ThighLen*z[9]*(q3p+q4p) - par__ShankLen*(z[175]*q3p+z[175]*
  q4p+z[176]*q5p);
  z[183] = par__ShankLen*(z[176]*q3p+z[176]*q4p+z[178]*q5p);
  z[184] = -par__ThighLen*z[101]*(q3p+q4p) - par__ShankLen*(z[103]*q3p+z[103]*
  q4p+z[104]*q5p);
  z[185] = par__ShankLen*(z[104]*q3p+z[104]*q4p+z[180]*q5p);
  z[186] = z[29]*z[173] + z[30]*z[175];
  z[187] = z[29]*z[103] + z[30]*z[176];
  z[188] = z[108]*q6p + z[186]*q3p + z[186]*q4p + z[187]*q5p;
  z[189] = z[29]*z[104] + z[30]*z[178];
  z[190] = z[109]*q6p + z[187]*q3p + z[187]*q4p + z[189]*q5p;
  z[191] = z[29]*z[175] - z[30]*z[173];
  z[192] = z[29]*z[176] - z[30]*z[103];
  z[193] = -z[29]*z[102] - z[30]*z[103];
  z[194] = z[191]*q3p + z[191]*q4p + z[192]*q5p + z[193]*q6p;
  z[195] = z[29]*z[178] - z[30]*z[104];
  z[196] = -z[18]*z[29] - z[30]*z[104];
  z[197] = z[192]*q3p + z[192]*q4p + z[195]*q5p + z[196]*q6p;
  z[198] = z[17]*z[30] - z[18]*z[29];
  z[199] = z[193]*q3p + z[193]*q4p + z[196]*q5p + z[198]*q6p;
  z[200] = z[29]*z[105] + z[30]*z[180];
  z[201] = z[107]*q3p + z[107]*q4p + z[112]*q6p + z[200]*q5p;
  z[202] = z[29]*z[180] - z[30]*z[105];
  z[203] = -z[20]*z[29] - z[30]*z[105];
  z[204] = z[109]*q3p + z[109]*q4p + z[202]*q5p + z[203]*q6p;
  z[205] = z[19]*z[30] - z[20]*z[29];
  z[206] = z[110]*q3p + z[110]*q4p + z[203]*q5p + z[205]*q6p;
  z[207] = par__ThighLen*z[57]*(q3p+q7p);
  z[208] = par__ThighLen*z[114]*(q3p+q7p);
  z[209] = z[57]*z[62] + z[61]*z[114];
  z[210] = z[116]*q8p + z[209]*q3p + z[209]*q7p;
  z[211] = z[57]*z[61] - z[62]*z[114];
  z[212] = -z[56]*z[61] - z[62]*z[114];
  z[213] = z[211]*q3p + z[211]*q7p + z[212]*q8p;
  z[214] = z[55]*z[62] - z[56]*z[61];
  z[215] = z[212]*q3p + z[212]*q7p + z[214]*q8p;
  z[216] = z[57]*z[62] - z[55]*z[61];
  z[217] = z[117]*q3p + z[117]*q7p + z[216]*q8p;
  z[218] = -par__ThighLen*z[57]*(q3p+q7p) - par__ShankLen*(z[211]*q3p+z[211]*
  q7p+z[212]*q8p);
  z[219] = par__ShankLen*(z[212]*q3p+z[212]*q7p+z[214]*q8p);
  z[220] = -par__ThighLen*z[114]*(q3p+q7p) - par__ShankLen*(z[116]*q3p+z[116]*
  q7p+z[117]*q8p);
  z[221] = par__ShankLen*(z[117]*q3p+z[117]*q7p+z[216]*q8p);
  z[222] = z[77]*z[209] + z[78]*z[211];
  z[223] = z[77]*z[116] + z[78]*z[212];
  z[224] = z[121]*q9p + z[222]*q3p + z[222]*q7p + z[223]*q8p;
  z[225] = z[77]*z[117] + z[78]*z[214];
  z[226] = z[122]*q9p + z[223]*q3p + z[223]*q7p + z[225]*q8p;
  z[227] = z[77]*z[211] - z[78]*z[209];
  z[228] = z[77]*z[212] - z[78]*z[116];
  z[229] = -z[77]*z[115] - z[78]*z[116];
  z[230] = z[227]*q3p + z[227]*q7p + z[228]*q8p + z[229]*q9p;
  z[231] = z[77]*z[214] - z[78]*z[117];
  z[232] = -z[66]*z[77] - z[78]*z[117];
  z[233] = z[228]*q3p + z[228]*q7p + z[231]*q8p + z[232]*q9p;
  z[234] = z[65]*z[78] - z[66]*z[77];
  z[235] = z[229]*q3p + z[229]*q7p + z[232]*q8p + z[234]*q9p;
  z[236] = z[77]*z[118] + z[78]*z[216];
  z[237] = z[120]*q3p + z[120]*q7p + z[125]*q9p + z[236]*q8p;
  z[238] = z[77]*z[216] - z[78]*z[118];
  z[239] = -z[68]*z[77] - z[78]*z[118];
  z[240] = z[122]*q3p + z[122]*q7p + z[238]*q8p + z[239]*q9p;
  z[241] = z[67]*z[78] - z[68]*z[77];
  z[242] = z[123]*q3p + z[123]*q7p + z[239]*q8p + z[241]*q9p;
  z[243] = z[9]*(q3p+q4p);
  z[244] = z[57]*(q3p+q7p);
  z[245] = sin(par__slope);
  z[246] = par__FootMass*par__gravity*z[245];
  z[247] = cos(par__slope);
  z[248] = par__FootMass*par__gravity*z[247];
  z[249] = par__gravity*par__ShankMass*z[245];
  z[250] = par__gravity*par__ShankMass*z[247];
  z[251] = par__gravity*par__ThighMass*z[245];
  z[252] = par__gravity*par__ThighMass*z[247];
  z[253] = par__gravity*par__TrunkMass*z[245];
  z[254] = par__gravity*par__TrunkMass*z[247];
  z[255] = z[1]*q1p + z[2]*q2p - par__TrunkCMy*q3p;
  z[256] = z[1]*q2p - z[2]*q1p;
  z[257] = par__airdrag*par__wind;
  z[258] = pow(par__wind,2);
  z[259] = z[257]*pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*
  z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[260] = par__airdrag*(par__TrunkCMy*q3p-z[1]*q1p-z[2]*q2p)*pow((z[258]+
  pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),
  0.5);
  z[261] = par__airdrag*(z[1]*q2p-z[2]*q1p)*pow((z[258]+pow(z[255],2)+pow(
  z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[262] = z[259] - z[253];
  z[266] = par__FootCMy + z[42];
  z[267] = par__FootCMy + z[43];
  z[268] = z[44] - par__FootCMx;
  z[269] = z[45] - par__FootCMx;
  z[270] = par__bodyweight*RFootFx;
  z[271] = par__bodyweight*RFootFy;
  z[273] = par__bodyweight*RFootMz;
  z[277] = par__FootCMy + z[90];
  z[278] = par__FootCMy + z[91];
  z[279] = z[92] - par__FootCMx;
  z[280] = z[93] - par__FootCMx;
  z[281] = par__bodyweight*LFootFx;
  z[282] = par__bodyweight*LFootFy;
  z[285] = par__bodyweight*LFootMz;
  z[287] = mom1 - mom2;
  z[288] = mom3 + z[273];
  z[289] = mom2 - mom3;
  z[290] = -mom1 - mom4;
  z[291] = mom4 - mom5;
  z[292] = mom6 + z[285];
  z[293] = mom5 - mom6;
  z[294] = z[34]*z[36]*z[271] + z[35]*z[37]*z[271] + z[82]*z[84]*z[282] + 
  z[83]*z[85]*z[282] + z[270]*pow(z[34],2) + z[270]*pow(z[35],2) + z[281]*
  pow(z[82],2) + z[281]*pow(z[83],2) + z[1]*z[260] + z[2]*z[261] + pow(z[1],2)*
  z[262] + pow(z[2],2)*z[262] - z[246]*pow(z[34],2) - z[246]*pow(z[35],2) - 
  z[246]*pow(z[82],2) - z[246]*pow(z[83],2) - z[248]*z[34]*z[36] - z[248]*
  z[35]*z[37] - z[248]*z[82]*z[84] - z[248]*z[83]*z[85] - z[249]*pow(z[17],2) - 
  z[249]*pow(z[18],2) - z[249]*pow(z[65],2) - z[249]*pow(z[66],2) - z[250]*
  z[17]*z[19] - z[250]*z[18]*z[20] - z[250]*z[65]*z[67] - z[250]*z[66]*z[68] - 
  z[251]*pow(z[7],2) - z[251]*pow(z[8],2) - z[251]*pow(z[55],2) - z[251]*pow(
  z[56],2) - z[252]*z[7]*z[8] - z[252]*z[7]*z[9] - z[252]*z[55]*z[56] - 
  z[252]*z[55]*z[57];
  z[295] = z[34]*z[36]*z[270] + z[35]*z[37]*z[270] + z[82]*z[84]*z[281] + 
  z[83]*z[85]*z[281] + z[271]*pow(z[36],2) + z[271]*pow(z[37],2) + z[282]*
  pow(z[84],2) + z[282]*pow(z[85],2) + z[2]*z[260] - z[246]*z[34]*z[36] - 
  z[246]*z[35]*z[37] - z[246]*z[82]*z[84] - z[246]*z[83]*z[85] - z[248]*pow(
  z[36],2) - z[248]*pow(z[37],2) - z[248]*pow(z[84],2) - z[248]*pow(z[85],2) - 
  z[249]*z[17]*z[19] - z[249]*z[18]*z[20] - z[249]*z[65]*z[67] - z[249]*z[66]*
  z[68] - z[250]*pow(z[19],2) - z[250]*pow(z[20],2) - z[250]*pow(z[67],2) - 
  z[250]*pow(z[68],2) - z[251]*z[7]*z[8] - z[251]*z[7]*z[9] - z[251]*z[55]*
  z[56] - z[251]*z[55]*z[57] - z[252]*pow(z[7],2) - z[252]*pow(z[9],2) - 
  z[252]*pow(z[55],2) - z[252]*pow(z[57],2) - z[254]*pow(z[1],2) - z[254]*
  pow(z[2],2) - z[1]*z[261];
  z[296] = z[287] + z[288] + z[289] + z[290] + z[291] + z[292] + z[293] + q2*
  z[270] + q2*z[281] + z[249]*z[18]*z[22] + z[249]*z[66]*z[70] + z[250]*z[20]*
  z[22] + z[250]*z[68]*z[70] + z[34]*z[267]*z[270] + z[35]*z[268]*z[270] + 
  z[36]*z[267]*z[271] + z[37]*z[268]*z[271] + z[82]*z[278]*z[281] + z[83]*
  z[279]*z[281] + z[84]*z[278]*z[282] + z[85]*z[279]*z[282] + par__ThighCMy*(
  z[251]*z[7]+z[252]*z[9]) + par__ThighCMy*(z[251]*z[55]+z[252]*z[57]) + 
  par__TrunkCMy*(z[254]*z[2]-z[260]-z[1]*z[262]) - q1*z[271] - q1*z[282] - 
  par__ShankLen*z[17]*z[270] - par__ShankLen*z[19]*z[271] - par__ShankLen*
  z[65]*z[281] - par__ShankLen*z[67]*z[282] - par__ThighLen*z[7]*z[270] - 
  par__ThighLen*z[9]*z[271] - par__ThighLen*z[55]*z[281] - par__ThighLen*
  z[57]*z[282] - z[246]*z[34]*z[43] - z[246]*z[35]*z[44] - z[246]*z[82]*z[91] - 
  z[246]*z[83]*z[92] - z[248]*z[36]*z[43] - z[248]*z[37]*z[44] - z[248]*z[84]*
  z[91] - z[248]*z[85]*z[92] - z[249]*z[17]*z[23] - z[249]*z[65]*z[71] - 
  z[250]*z[19]*z[23] - z[250]*z[67]*z[71];
  z[297] = z[287] + z[288] + z[289] + q2*z[270] + z[249]*z[18]*z[22] + z[250]*
  z[20]*z[22] + z[34]*z[267]*z[270] + z[35]*z[268]*z[270] + z[36]*z[267]*
  z[271] + z[37]*z[268]*z[271] + par__ThighCMy*(z[251]*z[7]+z[252]*z[9]) - q1*
  z[271] - par__ShankLen*z[17]*z[270] - par__ShankLen*z[19]*z[271] - par__ThighLen*
  z[7]*z[270] - par__ThighLen*z[9]*z[271] - z[246]*z[34]*z[43] - z[246]*z[35]*
  z[44] - z[248]*z[36]*z[43] - z[248]*z[37]*z[44] - z[249]*z[17]*z[23] - 
  z[250]*z[19]*z[23];
  z[298] = z[288] + z[289] + q2*z[270] + z[34]*z[266]*z[270] + z[35]*z[269]*
  z[270] + z[36]*z[266]*z[271] + z[37]*z[269]*z[271] + par__ShankCMy*(z[249]*
  z[17]+z[250]*z[19]) - q1*z[271] - par__ShankLen*z[17]*z[270] - par__ShankLen*
  z[19]*z[271] - par__ThighLen*z[7]*z[270] - par__ThighLen*z[9]*z[271] - 
  z[246]*z[34]*z[42] - z[246]*z[35]*z[45] - z[248]*z[36]*z[42] - z[248]*z[37]*
  z[45];
  z[299] = par__FootCMy*z[246];
  z[300] = par__FootCMy*z[248];
  z[301] = par__FootCMx*z[246];
  z[302] = par__FootCMx*z[248];
  z[303] = z[288] + z[299]*z[34] + z[300]*z[36] + q2*z[270] - z[301]*z[35] - 
  z[302]*z[37] - q1*z[271] - par__ShankLen*z[17]*z[270] - par__ShankLen*z[19]*
  z[271] - par__ThighLen*z[7]*z[270] - par__ThighLen*z[9]*z[271];
  z[304] = z[291] + z[292] + z[293] + q2*z[281] + z[249]*z[66]*z[70] + z[250]*
  z[68]*z[70] + z[82]*z[278]*z[281] + z[83]*z[279]*z[281] + z[84]*z[278]*
  z[282] + z[85]*z[279]*z[282] + par__ThighCMy*(z[251]*z[55]+z[252]*z[57]) - 
  q1*z[282] - par__ShankLen*z[65]*z[281] - par__ShankLen*z[67]*z[282] - 
  par__ThighLen*z[55]*z[281] - par__ThighLen*z[57]*z[282] - z[246]*z[82]*
  z[91] - z[246]*z[83]*z[92] - z[248]*z[84]*z[91] - z[248]*z[85]*z[92] - 
  z[249]*z[65]*z[71] - z[250]*z[67]*z[71];
  z[305] = z[292] + z[293] + q2*z[281] + z[82]*z[277]*z[281] + z[83]*z[280]*
  z[281] + z[84]*z[277]*z[282] + z[85]*z[280]*z[282] + par__ShankCMy*(z[249]*
  z[65]+z[250]*z[67]) - q1*z[282] - par__ShankLen*z[65]*z[281] - par__ShankLen*
  z[67]*z[282] - par__ThighLen*z[55]*z[281] - par__ThighLen*z[57]*z[282] - 
  z[246]*z[82]*z[90] - z[246]*z[83]*z[93] - z[248]*z[84]*z[90] - z[248]*z[85]*
  z[93];
  z[306] = z[292] + z[299]*z[82] + z[300]*z[84] + q2*z[281] - z[301]*z[83] - 
  z[302]*z[85] - q1*z[282] - par__ShankLen*z[65]*z[281] - par__ShankLen*z[67]*
  z[282] - par__ThighLen*z[55]*z[281] - par__ThighLen*z[57]*z[282];
  z[307] = par__FootMass*(par__FootCMx*z[83]-par__FootCMy*z[82]);
  z[308] = par__FootMass*(pow(z[34],2)+pow(z[35],2)) + par__FootMass*(pow(
  z[82],2)+pow(z[83],2)) + par__ShankMass*(pow(z[17],2)+pow(z[18],2)) + 
  par__ShankMass*(pow(z[65],2)+pow(z[66],2)) + par__ThighMass*(pow(z[7],2)+
  pow(z[8],2)) + par__ThighMass*(pow(z[55],2)+pow(z[56],2)) + par__TrunkMass*(
  pow(z[1],2)+pow(z[2],2));
  z[309] = par__ThighMass*z[7]*(z[8]+z[9]) + par__ThighMass*z[55]*(z[56]+
  z[57]) + par__FootMass*(z[34]*z[36]+z[35]*z[37]) + par__FootMass*(z[82]*
  z[84]+z[83]*z[85]) + par__ShankMass*(z[17]*z[19]+z[18]*z[20]) + par__ShankMass*(
  z[65]*z[67]+z[66]*z[68]);
  z[310] = par__ThighCMy*par__ThighMass;
  z[311] = par__TrunkCMy*par__TrunkMass;
  z[312] = par__FootMass*(z[34]*z[43]+z[35]*z[44]) + par__FootMass*(z[82]*
  z[91]+z[83]*z[92]) + par__ShankMass*(z[17]*z[23]-z[18]*z[22]) + par__ShankMass*(
  z[65]*z[71]-z[66]*z[70]) - z[310]*z[7] - z[310]*z[55] - z[311]*z[1];
  z[313] = par__FootMass*(z[82]*z[91]+z[83]*z[92]) + par__ShankMass*(z[65]*
  z[71]-z[66]*z[70]) - z[310]*z[55];
  z[314] = par__ShankCMy*par__ShankMass;
  z[315] = par__FootMass*(z[82]*z[90]+z[83]*z[93]) - z[314]*z[65];
  z[316] = par__FootMass*(par__FootCMx*z[35]-par__FootCMy*z[34]);
  z[317] = par__FootMass*(z[34]*z[43]+z[35]*z[44]) + par__ShankMass*(z[17]*
  z[23]-z[18]*z[22]) - z[310]*z[7];
  z[318] = par__FootMass*(z[34]*z[42]+z[35]*z[45]) - z[314]*z[17];
  z[319] = par__TrunkMass*z[2]*z[4] + par__FootMass*(z[34]*z[51]+z[35]*z[52]) + 
  par__FootMass*(z[82]*z[99]+z[83]*z[100]) + par__ShankMass*(z[17]*z[24]+
  z[18]*z[28]) + par__ShankMass*(z[65]*z[72]+z[66]*z[76]) - par__ThighMass*
  z[8]*z[12] - par__ThighMass*z[56]*z[60];
  z[320] = par__FootMass*(par__FootCMx*z[85]-par__FootCMy*z[84]);
  z[321] = par__FootMass*(pow(z[36],2)+pow(z[37],2)) + par__FootMass*(pow(
  z[84],2)+pow(z[85],2)) + par__ShankMass*(pow(z[19],2)+pow(z[20],2)) + 
  par__ShankMass*(pow(z[67],2)+pow(z[68],2)) + par__ThighMass*(pow(z[7],2)+
  pow(z[9],2)) + par__ThighMass*(pow(z[55],2)+pow(z[57],2)) + par__TrunkMass*(
  pow(z[1],2)+pow(z[2],2));
  z[322] = par__FootMass*(z[36]*z[43]+z[37]*z[44]) + par__FootMass*(z[84]*
  z[91]+z[85]*z[92]) + par__ShankMass*(z[19]*z[23]-z[20]*z[22]) + par__ShankMass*(
  z[67]*z[71]-z[68]*z[70]) - z[310]*z[9] - z[310]*z[57] - z[311]*z[2];
  z[323] = par__FootMass*(z[84]*z[91]+z[85]*z[92]) + par__ShankMass*(z[67]*
  z[71]-z[68]*z[70]) - z[310]*z[57];
  z[324] = par__FootMass*(z[84]*z[90]+z[85]*z[93]) - z[314]*z[67];
  z[325] = par__FootMass*(par__FootCMx*z[37]-par__FootCMy*z[36]);
  z[326] = par__FootMass*(z[36]*z[43]+z[37]*z[44]) + par__ShankMass*(z[19]*
  z[23]-z[20]*z[22]) - z[310]*z[9];
  z[327] = par__FootMass*(z[36]*z[42]+z[37]*z[45]) - z[314]*z[19];
  z[328] = par__FootMass*(z[36]*z[51]+z[37]*z[52]) + par__FootMass*(z[84]*
  z[99]+z[85]*z[100]) + par__ShankMass*(z[19]*z[24]+z[20]*z[28]) + par__ShankMass*(
  z[67]*z[72]+z[68]*z[76]) - par__ThighMass*z[7]*z[12] - par__ThighMass*z[55]*
  z[60] - par__TrunkMass*z[1]*z[4];
  z[329] = par__TrunkInertia + 2*par__FootInertia + 2*par__ShankInertia + 2*
  par__ThighInertia + par__TrunkMass*pow(par__TrunkCMy,2) + 2*par__ThighMass*
  pow(par__ThighCMy,2);
  z[330] = z[329] + par__FootMass*(pow(z[43],2)+pow(z[44],2)) + par__FootMass*(
  pow(z[91],2)+pow(z[92],2)) + par__ShankMass*(pow(z[22],2)+pow(z[23],2)) + 
  par__ShankMass*(pow(z[70],2)+pow(z[71],2));
  z[331] = par__FootInertia + par__ShankInertia + par__ThighInertia + par__ThighMass*
  pow(par__ThighCMy,2);
  z[332] = z[331] + par__FootMass*(pow(z[91],2)+pow(z[92],2)) + par__ShankMass*(
  pow(z[70],2)+pow(z[71],2));
  z[333] = par__FootInertia + par__ShankInertia;
  z[334] = z[333] + par__FootMass*(z[90]*z[91]+z[92]*z[93]) - z[314]*z[71];
  z[335] = par__FootInertia + par__FootMass*(par__FootCMx*z[92]-par__FootCMy*
  z[91]);
  z[336] = z[331] + par__FootMass*(pow(z[43],2)+pow(z[44],2)) + par__ShankMass*(
  pow(z[22],2)+pow(z[23],2));
  z[337] = z[333] + par__FootMass*(z[42]*z[43]+z[44]*z[45]) - z[314]*z[23];
  z[338] = par__FootInertia + par__FootMass*(par__FootCMx*z[44]-par__FootCMy*
  z[43]);
  z[339] = par__FootMass*(z[43]*z[51]+z[44]*z[52]) + par__FootMass*(z[91]*
  z[99]+z[92]*z[100]) - par__ShankMass*(z[22]*z[28]-z[23]*z[24]) - par__ShankMass*(
  z[70]*z[76]-z[71]*z[72]);
  z[340] = par__FootMass*(z[43]*z[51]+z[44]*z[52]) - par__ShankMass*(z[22]*
  z[28]-z[23]*z[24]);
  z[341] = par__FootInertia + par__ShankInertia + par__ShankMass*pow(par__ShankCMy,
  2);
  z[342] = z[341] + par__FootMass*(pow(z[42],2)+pow(z[45],2));
  z[343] = par__FootInertia + par__FootMass*(par__FootCMx*z[45]-par__FootCMy*
  z[42]);
  z[344] = par__FootMass*(z[42]*z[51]+z[45]*z[52]) - z[314]*z[24];
  z[345] = par__FootInertia + par__FootMass*(pow(par__FootCMx,2)+pow(par__FootCMy,
  2));
  z[346] = par__FootMass*(par__FootCMx*z[52]-par__FootCMy*z[51]);
  z[347] = par__FootMass*(z[91]*z[99]+z[92]*z[100]) - par__ShankMass*(z[70]*
  z[76]-z[71]*z[72]);
  z[348] = z[341] + par__FootMass*(pow(z[90],2)+pow(z[93],2));
  z[349] = par__FootInertia + par__FootMass*(par__FootCMx*z[93]-par__FootCMy*
  z[90]);
  z[350] = par__FootMass*(z[90]*z[99]+z[93]*z[100]) - z[314]*z[72];
  z[351] = par__FootMass*(par__FootCMx*z[100]-par__FootCMy*z[99]);
  z[352] = -z[1]*q1p - z[2]*q2p;
  z[353] = par__airdrag*((par__TrunkCMy*q3p-z[1]*q1p-z[2]*q2p)*(z[255]*z[256]+
  z[256]*z[352]+par__wind*z[2]*z[255]+par__wind*z[2]*z[352])+(z[1]*q2p-z[2]*
  q1p)*(2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-pow(z[256],2)-2*par__wind*
  z[2]*z[256]))/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*
  z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[354] = par__airdrag*((z[1]*q2p-z[2]*q1p)*(z[255]*z[256]+z[256]*z[352]+
  par__wind*z[2]*z[255]+par__wind*z[2]*z[352])+(z[1]*q1p+z[2]*q2p)*(2*par__wind*
  z[1]*z[255]-z[258]-pow(z[255],2)-pow(z[256],2)-2*par__wind*z[2]*z[256]))/
  pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*
  z[1]*z[255]),0.5);
  z[355] = z[257]*(z[255]*z[256]+z[256]*z[352]+par__wind*z[2]*z[255]+par__wind*
  z[2]*z[352])/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*
  z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[356] = z[36]*z[106]*z[271] + z[37]*z[108]*z[271] + z[84]*z[119]*z[282] + 
  z[85]*z[121]*z[282] + z[271]*pow(z[34],2) + z[271]*pow(z[35],2) + z[282]*
  pow(z[82],2) + z[282]*pow(z[83],2) + 2*z[34]*z[106]*z[270] + 2*z[35]*z[108]*
  z[270] + 2*z[82]*z[119]*z[281] + 2*z[83]*z[121]*z[281] + z[1]*z[261] + z[1]*
  z[353] + z[2]*z[354] + pow(z[1],2)*z[355] + pow(z[2],2)*z[355] - 2*z[246]*
  z[34]*z[106] - 2*z[246]*z[35]*z[108] - 2*z[246]*z[82]*z[119] - 2*z[246]*
  z[83]*z[121] - 2*z[249]*z[17]*z[102] - 2*z[249]*z[18]*z[103] - 2*z[249]*
  z[65]*z[115] - 2*z[249]*z[66]*z[116] - 2*z[251]*z[7]*z[8] - 2*z[251]*z[8]*
  z[101] - 2*z[251]*z[55]*z[56] - 2*z[251]*z[56]*z[114] - z[248]*pow(z[34],2) - 
  z[248]*pow(z[35],2) - z[248]*z[36]*z[106] - z[248]*z[37]*z[108] - z[248]*
  pow(z[82],2) - z[248]*pow(z[83],2) - z[248]*z[84]*z[119] - z[248]*z[85]*
  z[121] - z[250]*pow(z[17],2) - z[250]*pow(z[18],2) - z[250]*z[19]*z[102] - 
  z[250]*z[20]*z[103] - z[250]*pow(z[65],2) - z[250]*pow(z[66],2) - z[250]*
  z[67]*z[115] - z[250]*z[68]*z[116] - z[252]*z[7]*z[101] - z[252]*pow(z[7],2) - 
  z[252]*z[8]*z[9] - z[252]*pow(z[8],2) - z[252]*z[55]*z[114] - z[252]*pow(
  z[55],2) - z[252]*z[56]*z[57] - z[252]*pow(z[56],2) - z[2]*z[260];
  z[357] = par__TrunkMass*z[1]*z[4] + par__FootMass*(z[106]*z[51]+z[108]*
  z[52]) + par__FootMass*(z[119]*z[99]+z[121]*z[100]) + par__ShankMass*(
  z[102]*z[24]+z[103]*z[28]) + par__ShankMass*(z[115]*z[72]+z[116]*z[76]) - 
  par__ThighMass*z[101]*z[12] - par__ThighMass*z[114]*z[60];
  z[358] = par__FootMass*(par__FootCMx*z[121]-par__FootCMy*z[119]);
  z[359] = 2*par__ThighMass*z[8]*(z[7]+z[101]) + 2*par__ThighMass*z[56]*(
  z[55]+z[114]) + 2*par__FootMass*(z[34]*z[106]+z[35]*z[108]) + 2*par__FootMass*(
  z[82]*z[119]+z[83]*z[121]) + 2*par__ShankMass*(z[17]*z[102]+z[18]*z[103]) + 
  2*par__ShankMass*(z[65]*z[115]+z[66]*z[116]);
  z[360] = par__ThighMass*z[7]*(z[7]+z[101]) + par__ThighMass*z[8]*(z[8]+z[9]) + 
  par__ThighMass*z[55]*(z[55]+z[114]) + par__ThighMass*z[56]*(z[56]+z[57]) + 
  par__FootMass*(pow(z[34],2)+pow(z[35],2)+z[36]*z[106]+z[37]*z[108]) + 
  par__FootMass*(pow(z[82],2)+pow(z[83],2)+z[84]*z[119]+z[85]*z[121]) + 
  par__ShankMass*(pow(z[17],2)+pow(z[18],2)+z[19]*z[102]+z[20]*z[103]) + 
  par__ShankMass*(pow(z[65],2)+pow(z[66],2)+z[67]*z[115]+z[68]*z[116]);
  z[361] = z[311]*z[2] + par__FootMass*(z[43]*z[106]+z[44]*z[108]) + par__FootMass*(
  z[91]*z[119]+z[92]*z[121]) - z[310]*z[8] - z[310]*z[56] - par__ShankMass*(
  z[22]*z[103]-z[23]*z[102]) - par__ShankMass*(z[70]*z[116]-z[71]*z[115]);
  z[362] = par__FootMass*(z[91]*z[119]+z[92]*z[121]) - z[310]*z[56] - par__ShankMass*(
  z[70]*z[116]-z[71]*z[115]);
  z[363] = par__FootMass*(z[90]*z[119]+z[93]*z[121]) - z[314]*z[115];
  z[364] = par__FootMass*(par__FootCMx*z[108]-par__FootCMy*z[106]);
  z[365] = par__FootMass*(z[43]*z[106]+z[44]*z[108]) - z[310]*z[8] - par__ShankMass*(
  z[22]*z[103]-z[23]*z[102]);
  z[366] = par__FootMass*(z[42]*z[106]+z[45]*z[108]) - z[314]*z[102];
  z[367] = z[36]*z[106]*z[271] + z[37]*z[108]*z[271] + z[271]*pow(z[34],2) + 
  z[271]*pow(z[35],2) + 2*z[34]*z[106]*z[270] + 2*z[35]*z[108]*z[270] - 2*
  z[246]*z[34]*z[106] - 2*z[246]*z[35]*z[108] - 2*z[249]*z[17]*z[102] - 2*
  z[249]*z[18]*z[103] - 2*z[251]*z[7]*z[8] - 2*z[251]*z[8]*z[101] - z[248]*
  pow(z[34],2) - z[248]*pow(z[35],2) - z[248]*z[36]*z[106] - z[248]*z[37]*
  z[108] - z[250]*pow(z[17],2) - z[250]*pow(z[18],2) - z[250]*z[19]*z[102] - 
  z[250]*z[20]*z[103] - z[252]*z[7]*z[101] - z[252]*pow(z[7],2) - z[252]*z[8]*
  z[9] - z[252]*pow(z[8],2);
  z[368] = par__FootMass*(z[106]*z[51]+z[108]*z[52]) + par__ShankMass*(z[102]*
  z[24]+z[103]*z[28]) - par__ThighMass*z[101]*z[12];
  z[369] = 2*par__ThighMass*z[8]*(z[7]+z[101]) + 2*par__FootMass*(z[34]*
  z[106]+z[35]*z[108]) + 2*par__ShankMass*(z[17]*z[102]+z[18]*z[103]);
  z[370] = par__ThighMass*z[7]*(z[7]+z[101]) + par__ThighMass*z[8]*(z[8]+z[9]) + 
  par__FootMass*(pow(z[34],2)+pow(z[35],2)+z[36]*z[106]+z[37]*z[108]) + 
  par__ShankMass*(pow(z[17],2)+pow(z[18],2)+z[19]*z[102]+z[20]*z[103]);
  z[371] = z[34]*z[111]*z[271] + z[35]*z[112]*z[271] + z[36]*z[107]*z[271] + 
  z[37]*z[109]*z[271] + 2*z[34]*z[107]*z[270] + 2*z[35]*z[109]*z[270] - 2*
  z[246]*z[34]*z[107] - 2*z[246]*z[35]*z[109] - 2*z[249]*z[17]*z[18] - 2*
  z[249]*z[18]*z[104] - z[248]*z[34]*z[111] - z[248]*z[35]*z[112] - z[248]*
  z[36]*z[107] - z[248]*z[37]*z[109] - z[250]*z[17]*z[20] - z[250]*z[18]*
  z[19] - z[250]*z[18]*z[105] - z[250]*z[20]*z[104];
  z[372] = (z[13]*z[29]-z[14]*z[30])*z[16];
  z[373] = (z[13]*z[30]+z[14]*z[29])*z[16];
  z[374] = par__ShankMass*(z[18]*z[24]+z[104]*z[28]+z[13]*z[17]*z[16]-z[14]*
  z[18]*z[16]) - par__FootMass*(z[35]*z[373]-z[34]*z[372]-z[107]*z[51]-z[109]*
  z[52]);
  z[375] = 2*par__ShankMass*z[18]*(z[17]+z[104]) + 2*par__FootMass*(z[34]*
  z[107]+z[35]*z[109]);
  z[376] = par__FootMass*(z[34]*z[111]+z[35]*z[112]+z[36]*z[107]+z[37]*z[109]) + 
  par__ShankMass*(z[17]*z[20]+z[18]*z[19]+z[18]*z[105]+z[20]*z[104]);
  z[377] = par__ThighLen*(z[13]*z[30]+z[14]*z[29]);
  z[378] = par__ThighLen*(z[13]*z[29]-z[14]*z[30]);
  z[379] = par__ShankMass*(z[18]*z[23]-z[22]*z[104]-par__ThighLen*z[13]*z[18]-
  par__ThighLen*z[14]*z[17]) - par__FootMass*(z[34]*z[377]+z[35]*z[378]-z[43]*
  z[107]-z[44]*z[109]);
  z[380] = par__FootMass*(par__FootCMx*z[109]-par__FootCMy*z[107]);
  z[381] = par__FootMass*(z[42]*z[107]+z[45]*z[109]) - z[314]*z[18];
  z[382] = z[34]*z[37]*z[271] + z[35]*z[36]*z[271] + z[35]*z[113]*z[271] + 
  z[37]*z[110]*z[271] + 2*z[34]*z[35]*z[270] + 2*z[35]*z[110]*z[270] - 2*
  z[246]*z[34]*z[35] - 2*z[246]*z[35]*z[110] - z[248]*z[34]*z[37] - z[248]*
  z[35]*z[36] - z[248]*z[35]*z[113] - z[248]*z[37]*z[110];
  z[383] = -z[29]*z[24] - z[30]*z[33];
  z[384] = par__FootMass*(z[34]*z[47]+z[35]*z[51]+z[35]*z[383]+z[110]*z[52]);
  z[385] = par__FootMass*z[35]*(z[34]+z[110]);
  z[386] = par__FootMass*(z[34]*z[37]+z[35]*z[36]+z[35]*z[113]+z[37]*z[110]);
  z[387] = z[22]*z[30] - z[29]*z[31];
  z[388] = par__FootMass*(z[34]*z[40]+z[35]*z[43]+z[35]*z[387]+z[44]*z[110]);
  z[389] = par__FootMass*(par__FootCMx*z[110]-par__FootCMy*z[35]);
  z[390] = par__FootMass*(z[35]*z[42]+z[45]*z[110]-par__ShankLen*z[29]*z[35]-
  par__ShankLen*z[30]*z[34]);
  z[391] = z[84]*z[119]*z[282] + z[85]*z[121]*z[282] + z[282]*pow(z[82],2) + 
  z[282]*pow(z[83],2) + 2*z[82]*z[119]*z[281] + 2*z[83]*z[121]*z[281] - 2*
  z[246]*z[82]*z[119] - 2*z[246]*z[83]*z[121] - 2*z[249]*z[65]*z[115] - 2*
  z[249]*z[66]*z[116] - 2*z[251]*z[55]*z[56] - 2*z[251]*z[56]*z[114] - z[248]*
  pow(z[82],2) - z[248]*pow(z[83],2) - z[248]*z[84]*z[119] - z[248]*z[85]*
  z[121] - z[250]*pow(z[65],2) - z[250]*pow(z[66],2) - z[250]*z[67]*z[115] - 
  z[250]*z[68]*z[116] - z[252]*z[55]*z[114] - z[252]*pow(z[55],2) - z[252]*
  z[56]*z[57] - z[252]*pow(z[56],2);
  z[392] = par__FootMass*(z[119]*z[99]+z[121]*z[100]) + par__ShankMass*(
  z[115]*z[72]+z[116]*z[76]) - par__ThighMass*z[114]*z[60];
  z[393] = 2*par__ThighMass*z[56]*(z[55]+z[114]) + 2*par__FootMass*(z[82]*
  z[119]+z[83]*z[121]) + 2*par__ShankMass*(z[65]*z[115]+z[66]*z[116]);
  z[394] = par__ThighMass*z[55]*(z[55]+z[114]) + par__ThighMass*z[56]*(z[56]+
  z[57]) + par__FootMass*(pow(z[82],2)+pow(z[83],2)+z[84]*z[119]+z[85]*z[121]) + 
  par__ShankMass*(pow(z[65],2)+pow(z[66],2)+z[67]*z[115]+z[68]*z[116]);
  z[395] = z[82]*z[124]*z[282] + z[83]*z[125]*z[282] + z[84]*z[120]*z[282] + 
  z[85]*z[122]*z[282] + 2*z[82]*z[120]*z[281] + 2*z[83]*z[122]*z[281] - 2*
  z[246]*z[82]*z[120] - 2*z[246]*z[83]*z[122] - 2*z[249]*z[65]*z[66] - 2*
  z[249]*z[66]*z[117] - z[248]*z[82]*z[124] - z[248]*z[83]*z[125] - z[248]*
  z[84]*z[120] - z[248]*z[85]*z[122] - z[250]*z[65]*z[68] - z[250]*z[66]*
  z[67] - z[250]*z[66]*z[118] - z[250]*z[68]*z[117];
  z[396] = (z[61]*z[77]-z[62]*z[78])*z[64];
  z[397] = (z[61]*z[78]+z[62]*z[77])*z[64];
  z[398] = par__ShankMass*(z[66]*z[72]+z[117]*z[76]+z[61]*z[65]*z[64]-z[62]*
  z[66]*z[64]) - par__FootMass*(z[83]*z[397]-z[82]*z[396]-z[120]*z[99]-z[122]*
  z[100]);
  z[399] = par__FootMass*(par__FootCMx*z[122]-par__FootCMy*z[120]);
  z[400] = 2*par__ShankMass*z[66]*(z[65]+z[117]) + 2*par__FootMass*(z[82]*
  z[120]+z[83]*z[122]);
  z[401] = par__FootMass*(z[82]*z[124]+z[83]*z[125]+z[84]*z[120]+z[85]*z[122]) + 
  par__ShankMass*(z[65]*z[68]+z[66]*z[67]+z[66]*z[118]+z[68]*z[117]);
  z[402] = par__ThighLen*(z[61]*z[78]+z[62]*z[77]);
  z[403] = par__ThighLen*(z[61]*z[77]-z[62]*z[78]);
  z[404] = par__ShankMass*(z[66]*z[71]-z[70]*z[117]-par__ThighLen*z[61]*z[66]-
  par__ThighLen*z[62]*z[65]) - par__FootMass*(z[82]*z[402]+z[83]*z[403]-z[91]*
  z[120]-z[92]*z[122]);
  z[405] = par__FootMass*(z[90]*z[120]+z[93]*z[122]) - z[314]*z[66];
  z[406] = z[82]*z[85]*z[282] + z[83]*z[84]*z[282] + z[83]*z[126]*z[282] + 
  z[85]*z[123]*z[282] + 2*z[82]*z[83]*z[281] + 2*z[83]*z[123]*z[281] - 2*
  z[246]*z[82]*z[83] - 2*z[246]*z[83]*z[123] - z[248]*z[82]*z[85] - z[248]*
  z[83]*z[84] - z[248]*z[83]*z[126] - z[248]*z[85]*z[123];
  z[407] = -z[77]*z[72] - z[78]*z[81];
  z[408] = par__FootMass*(z[82]*z[95]+z[83]*z[99]+z[83]*z[407]+z[123]*z[100]);
  z[409] = par__FootMass*(par__FootCMx*z[123]-par__FootCMy*z[83]);
  z[410] = par__FootMass*z[83]*(z[82]+z[123]);
  z[411] = par__FootMass*(z[82]*z[85]+z[83]*z[84]+z[83]*z[126]+z[85]*z[123]);
  z[412] = z[70]*z[78] - z[77]*z[79];
  z[413] = par__FootMass*(z[82]*z[88]+z[83]*z[91]+z[83]*z[412]+z[92]*z[123]);
  z[414] = par__FootMass*(z[83]*z[90]+z[93]*z[123]-par__ShankLen*z[77]*z[83]-
  par__ShankLen*z[78]*z[82]);
  z[415] = z[36]*z[106]*z[270] + z[37]*z[108]*z[270] + z[84]*z[119]*z[281] + 
  z[85]*z[121]*z[281] + z[270]*pow(z[34],2) + z[270]*pow(z[35],2) + z[281]*
  pow(z[82],2) + z[281]*pow(z[83],2) + 2*z[34]*z[36]*z[271] + 2*z[35]*z[37]*
  z[271] + 2*z[82]*z[84]*z[282] + 2*z[83]*z[85]*z[282] + z[1]*z[260] + z[2]*
  z[261] + z[2]*z[353] - 2*z[248]*z[34]*z[36] - 2*z[248]*z[35]*z[37] - 2*
  z[248]*z[82]*z[84] - 2*z[248]*z[83]*z[85] - 2*z[250]*z[17]*z[19] - 2*z[250]*
  z[18]*z[20] - 2*z[250]*z[65]*z[67] - 2*z[250]*z[66]*z[68] - 2*z[252]*z[7]*
  z[8] - 2*z[252]*z[7]*z[9] - 2*z[252]*z[55]*z[56] - 2*z[252]*z[55]*z[57] - 
  z[246]*pow(z[34],2) - z[246]*pow(z[35],2) - z[246]*z[36]*z[106] - z[246]*
  z[37]*z[108] - z[246]*pow(z[82],2) - z[246]*pow(z[83],2) - z[246]*z[84]*
  z[119] - z[246]*z[85]*z[121] - z[249]*pow(z[17],2) - z[249]*pow(z[18],2) - 
  z[249]*z[19]*z[102] - z[249]*z[20]*z[103] - z[249]*pow(z[65],2) - z[249]*
  pow(z[66],2) - z[249]*z[67]*z[115] - z[249]*z[68]*z[116] - z[251]*z[7]*
  z[101] - z[251]*pow(z[7],2) - z[251]*z[8]*z[9] - z[251]*pow(z[8],2) - 
  z[251]*z[55]*z[114] - z[251]*pow(z[55],2) - z[251]*z[56]*z[57] - z[251]*
  pow(z[56],2) - z[1]*z[354];
  z[416] = 2*par__ThighMass*z[7]*(z[8]+z[9]) + 2*par__ThighMass*z[55]*(z[56]+
  z[57]) + 2*par__FootMass*(z[34]*z[36]+z[35]*z[37]) + 2*par__FootMass*(z[82]*
  z[84]+z[83]*z[85]) + 2*par__ShankMass*(z[17]*z[19]+z[18]*z[20]) + 2*par__ShankMass*(
  z[65]*z[67]+z[66]*z[68]);
  z[417] = z[36]*z[106]*z[270] + z[37]*z[108]*z[270] + z[270]*pow(z[34],2) + 
  z[270]*pow(z[35],2) + 2*z[34]*z[36]*z[271] + 2*z[35]*z[37]*z[271] - 2*
  z[248]*z[34]*z[36] - 2*z[248]*z[35]*z[37] - 2*z[250]*z[17]*z[19] - 2*z[250]*
  z[18]*z[20] - 2*z[252]*z[7]*z[8] - 2*z[252]*z[7]*z[9] - z[246]*pow(z[34],2) - 
  z[246]*pow(z[35],2) - z[246]*z[36]*z[106] - z[246]*z[37]*z[108] - z[249]*
  pow(z[17],2) - z[249]*pow(z[18],2) - z[249]*z[19]*z[102] - z[249]*z[20]*
  z[103] - z[251]*z[7]*z[101] - z[251]*pow(z[7],2) - z[251]*z[8]*z[9] - 
  z[251]*pow(z[8],2);
  z[418] = par__FootMass*(z[34]*z[51]+z[35]*z[52]) + par__ShankMass*(z[17]*
  z[24]+z[18]*z[28]) - par__ThighMass*z[8]*z[12];
  z[419] = 2*par__ThighMass*z[7]*(z[8]+z[9]) + 2*par__FootMass*(z[34]*z[36]+
  z[35]*z[37]) + 2*par__ShankMass*(z[17]*z[19]+z[18]*z[20]);
  z[420] = z[34]*z[111]*z[270] + z[35]*z[112]*z[270] + z[36]*z[107]*z[270] + 
  z[37]*z[109]*z[270] + 2*z[36]*z[111]*z[271] + 2*z[37]*z[112]*z[271] - 2*
  z[248]*z[36]*z[111] - 2*z[248]*z[37]*z[112] - 2*z[250]*z[19]*z[20] - 2*
  z[250]*z[20]*z[105] - z[246]*z[34]*z[111] - z[246]*z[35]*z[112] - z[246]*
  z[36]*z[107] - z[246]*z[37]*z[109] - z[249]*z[17]*z[20] - z[249]*z[18]*
  z[19] - z[249]*z[18]*z[105] - z[249]*z[20]*z[104];
  z[421] = par__ShankMass*(z[20]*z[24]+z[105]*z[28]+z[13]*z[19]*z[16]-z[14]*
  z[20]*z[16]) - par__FootMass*(z[37]*z[373]-z[36]*z[372]-z[111]*z[51]-z[112]*
  z[52]);
  z[422] = 2*par__ShankMass*z[20]*(z[19]+z[105]) + 2*par__FootMass*(z[36]*
  z[111]+z[37]*z[112]);
  z[423] = par__ShankMass*(z[20]*z[23]-z[22]*z[105]-par__ThighLen*z[13]*z[20]-
  par__ThighLen*z[14]*z[19]) - par__FootMass*(z[36]*z[377]+z[37]*z[378]-z[43]*
  z[111]-z[44]*z[112]);
  z[424] = par__FootMass*(par__FootCMx*z[112]-par__FootCMy*z[111]);
  z[425] = par__FootMass*(z[42]*z[111]+z[45]*z[112]) - z[314]*z[20];
  z[426] = z[34]*z[37]*z[270] + z[35]*z[36]*z[270] + z[35]*z[113]*z[270] + 
  z[37]*z[110]*z[270] + 2*z[36]*z[37]*z[271] + 2*z[37]*z[113]*z[271] - 2*
  z[248]*z[36]*z[37] - 2*z[248]*z[37]*z[113] - z[246]*z[34]*z[37] - z[246]*
  z[35]*z[36] - z[246]*z[35]*z[113] - z[246]*z[37]*z[110];
  z[427] = par__FootMass*(z[36]*z[47]+z[37]*z[51]+z[37]*z[383]+z[113]*z[52]);
  z[428] = par__FootMass*z[37]*(z[36]+z[113]);
  z[429] = par__FootMass*(z[36]*z[40]+z[37]*z[43]+z[37]*z[387]+z[44]*z[113]);
  z[430] = par__FootMass*(par__FootCMx*z[113]-par__FootCMy*z[37]);
  z[431] = par__FootMass*(z[37]*z[42]+z[45]*z[113]-par__ShankLen*z[29]*z[37]-
  par__ShankLen*z[30]*z[36]);
  z[432] = z[84]*z[119]*z[281] + z[85]*z[121]*z[281] + z[281]*pow(z[82],2) + 
  z[281]*pow(z[83],2) + 2*z[82]*z[84]*z[282] + 2*z[83]*z[85]*z[282] - 2*
  z[248]*z[82]*z[84] - 2*z[248]*z[83]*z[85] - 2*z[250]*z[65]*z[67] - 2*z[250]*
  z[66]*z[68] - 2*z[252]*z[55]*z[56] - 2*z[252]*z[55]*z[57] - z[246]*pow(
  z[82],2) - z[246]*pow(z[83],2) - z[246]*z[84]*z[119] - z[246]*z[85]*z[121] - 
  z[249]*pow(z[65],2) - z[249]*pow(z[66],2) - z[249]*z[67]*z[115] - z[249]*
  z[68]*z[116] - z[251]*z[55]*z[114] - z[251]*pow(z[55],2) - z[251]*z[56]*
  z[57] - z[251]*pow(z[56],2);
  z[433] = par__FootMass*(z[82]*z[99]+z[83]*z[100]) + par__ShankMass*(z[65]*
  z[72]+z[66]*z[76]) - par__ThighMass*z[56]*z[60];
  z[434] = 2*par__ThighMass*z[55]*(z[56]+z[57]) + 2*par__FootMass*(z[82]*
  z[84]+z[83]*z[85]) + 2*par__ShankMass*(z[65]*z[67]+z[66]*z[68]);
  z[435] = z[82]*z[124]*z[281] + z[83]*z[125]*z[281] + z[84]*z[120]*z[281] + 
  z[85]*z[122]*z[281] + 2*z[84]*z[124]*z[282] + 2*z[85]*z[125]*z[282] - 2*
  z[248]*z[84]*z[124] - 2*z[248]*z[85]*z[125] - 2*z[250]*z[67]*z[68] - 2*
  z[250]*z[68]*z[118] - z[246]*z[82]*z[124] - z[246]*z[83]*z[125] - z[246]*
  z[84]*z[120] - z[246]*z[85]*z[122] - z[249]*z[65]*z[68] - z[249]*z[66]*
  z[67] - z[249]*z[66]*z[118] - z[249]*z[68]*z[117];
  z[436] = par__ShankMass*(z[68]*z[72]+z[118]*z[76]+z[61]*z[67]*z[64]-z[62]*
  z[68]*z[64]) - par__FootMass*(z[85]*z[397]-z[84]*z[396]-z[124]*z[99]-z[125]*
  z[100]);
  z[437] = par__FootMass*(par__FootCMx*z[125]-par__FootCMy*z[124]);
  z[438] = 2*par__ShankMass*z[68]*(z[67]+z[118]) + 2*par__FootMass*(z[84]*
  z[124]+z[85]*z[125]);
  z[439] = par__ShankMass*(z[68]*z[71]-z[70]*z[118]-par__ThighLen*z[61]*z[68]-
  par__ThighLen*z[62]*z[67]) - par__FootMass*(z[84]*z[402]+z[85]*z[403]-z[91]*
  z[124]-z[92]*z[125]);
  z[440] = par__FootMass*(z[90]*z[124]+z[93]*z[125]) - z[314]*z[68];
  z[441] = z[82]*z[85]*z[281] + z[83]*z[84]*z[281] + z[83]*z[126]*z[281] + 
  z[85]*z[123]*z[281] + 2*z[84]*z[85]*z[282] + 2*z[85]*z[126]*z[282] - 2*
  z[248]*z[84]*z[85] - 2*z[248]*z[85]*z[126] - z[246]*z[82]*z[85] - z[246]*
  z[83]*z[84] - z[246]*z[83]*z[126] - z[246]*z[85]*z[123];
  z[442] = par__FootMass*(z[84]*z[95]+z[85]*z[99]+z[85]*z[407]+z[126]*z[100]);
  z[443] = par__FootMass*(par__FootCMx*z[126]-par__FootCMy*z[85]);
  z[444] = par__FootMass*z[85]*(z[84]+z[126]);
  z[445] = par__FootMass*(z[84]*z[88]+z[85]*z[91]+z[85]*z[412]+z[92]*z[126]);
  z[446] = par__FootMass*(z[85]*z[90]+z[93]*z[126]-par__ShankLen*z[77]*z[85]-
  par__ShankLen*z[78]*z[84]);
  z[447] = z[249]*z[22]*z[103] + z[249]*z[70]*z[116] + z[250]*z[18]*z[22] + 
  z[250]*z[66]*z[70] + z[34]*z[267]*z[271] + z[35]*z[268]*z[271] + z[82]*
  z[278]*z[282] + z[83]*z[279]*z[282] + z[106]*z[267]*z[270] + z[108]*z[268]*
  z[270] + z[119]*z[278]*z[281] + z[121]*z[279]*z[281] + par__ThighCMy*(
  z[251]*z[8]+z[252]*z[7]) + par__ThighCMy*(z[251]*z[56]+z[252]*z[55]) - 
  par__ShankLen*z[17]*z[271] - par__ShankLen*z[65]*z[282] - par__ShankLen*
  z[102]*z[270] - par__ShankLen*z[115]*z[281] - par__ThighLen*z[7]*z[271] - 
  par__ThighLen*z[8]*z[270] - par__ThighLen*z[55]*z[282] - par__ThighLen*
  z[56]*z[281] - z[246]*z[43]*z[106] - z[246]*z[44]*z[108] - z[246]*z[91]*
  z[119] - z[246]*z[92]*z[121] - z[248]*z[34]*z[43] - z[248]*z[35]*z[44] - 
  z[248]*z[82]*z[91] - z[248]*z[83]*z[92] - z[249]*z[23]*z[102] - z[249]*
  z[71]*z[115] - z[250]*z[17]*z[23] - z[250]*z[65]*z[71] - par__TrunkCMy*(
  z[353]+z[1]*z[355]-z[254]*z[1]-z[2]*z[262]);
  z[448] = z[249]*z[22]*z[103] + z[250]*z[18]*z[22] + z[34]*z[267]*z[271] + 
  z[35]*z[268]*z[271] + z[106]*z[267]*z[270] + z[108]*z[268]*z[270] + par__ThighCMy*(
  z[251]*z[8]+z[252]*z[7]) - par__ShankLen*z[17]*z[271] - par__ShankLen*
  z[102]*z[270] - par__ThighLen*z[7]*z[271] - par__ThighLen*z[8]*z[270] - 
  z[246]*z[43]*z[106] - z[246]*z[44]*z[108] - z[248]*z[34]*z[43] - z[248]*
  z[35]*z[44] - z[249]*z[23]*z[102] - z[250]*z[17]*z[23];
  z[449] = par__ThighLen*z[249];
  z[450] = par__ThighLen*z[250];
  z[451] = z[246]*z[34]*z[377] + z[246]*z[35]*z[378] + z[248]*z[36]*z[377] + 
  z[248]*z[37]*z[378] + z[249]*z[22]*z[104] + z[250]*z[22]*z[105] + z[449]*
  z[13]*z[18] + z[449]*z[14]*z[17] + z[450]*z[13]*z[20] + z[450]*z[14]*z[19] + 
  z[107]*z[267]*z[270] + z[109]*z[268]*z[270] + z[111]*z[267]*z[271] + z[112]*
  z[268]*z[271] - par__ShankLen*z[18]*z[270] - par__ShankLen*z[20]*z[271] - 
  z[246]*z[43]*z[107] - z[246]*z[44]*z[109] - z[248]*z[43]*z[111] - z[248]*
  z[44]*z[112] - z[249]*z[18]*z[23] - z[250]*z[20]*z[23] - z[34]*z[270]*
  z[377] - z[35]*z[270]*z[378] - z[36]*z[271]*z[377] - z[37]*z[271]*z[378];
  z[452] = par__FootMass*(z[43]*z[372]-z[44]*z[373]-z[377]*z[51]-z[378]*z[52]) - 
  par__ShankMass*(par__ThighLen*z[13]*z[28]+par__ThighLen*z[14]*z[24]-z[13]*
  z[23]*z[16]-z[14]*z[22]*z[16]);
  z[453] = par__ShankMass*par__ThighLen;
  z[454] = 2*z[453]*(z[13]*z[22]-z[14]*z[23]) - 2*par__FootMass*(z[43]*z[377]+
  z[44]*z[378]);
  z[455] = par__ThighLen*z[314];
  z[456] = z[455]*z[14] - par__FootMass*(z[42]*z[377]+z[45]*z[378]);
  z[457] = par__FootMass*(par__FootCMx*z[378]-par__FootCMy*z[377]);
  z[458] = z[34]*z[40]*z[270] + z[35]*z[267]*z[270] + z[35]*z[270]*z[387] + 
  z[36]*z[40]*z[271] + z[37]*z[267]*z[271] + z[37]*z[271]*z[387] + z[110]*
  z[268]*z[270] + z[113]*z[268]*z[271] - z[246]*z[34]*z[40] - z[246]*z[35]*
  z[43] - z[246]*z[35]*z[387] - z[246]*z[44]*z[110] - z[248]*z[36]*z[40] - 
  z[248]*z[37]*z[43] - z[248]*z[37]*z[387] - z[248]*z[44]*z[113];
  z[459] = par__FootMass*(z[40]*z[51]+z[43]*z[47]+z[44]*z[383]+z[387]*z[52]);
  z[460] = par__FootMass*(z[40]*z[43]+z[44]*z[387]);
  z[461] = par__FootMass*(z[40]*z[42]+z[45]*z[387]-par__ShankLen*z[29]*z[44]-
  par__ShankLen*z[30]*z[43]);
  z[462] = par__FootMass*(par__FootCMx*z[387]-par__FootCMy*z[40]);
  z[463] = z[249]*z[70]*z[116] + z[250]*z[66]*z[70] + z[82]*z[278]*z[282] + 
  z[83]*z[279]*z[282] + z[119]*z[278]*z[281] + z[121]*z[279]*z[281] + par__ThighCMy*(
  z[251]*z[56]+z[252]*z[55]) - par__ShankLen*z[65]*z[282] - par__ShankLen*
  z[115]*z[281] - par__ThighLen*z[55]*z[282] - par__ThighLen*z[56]*z[281] - 
  z[246]*z[91]*z[119] - z[246]*z[92]*z[121] - z[248]*z[82]*z[91] - z[248]*
  z[83]*z[92] - z[249]*z[71]*z[115] - z[250]*z[65]*z[71];
  z[464] = z[246]*z[82]*z[402] + z[246]*z[83]*z[403] + z[248]*z[84]*z[402] + 
  z[248]*z[85]*z[403] + z[249]*z[70]*z[117] + z[250]*z[70]*z[118] + z[449]*
  z[61]*z[66] + z[449]*z[62]*z[65] + z[450]*z[61]*z[68] + z[450]*z[62]*z[67] + 
  z[120]*z[278]*z[281] + z[122]*z[279]*z[281] + z[124]*z[278]*z[282] + z[125]*
  z[279]*z[282] - par__ShankLen*z[66]*z[281] - par__ShankLen*z[68]*z[282] - 
  z[246]*z[91]*z[120] - z[246]*z[92]*z[122] - z[248]*z[91]*z[124] - z[248]*
  z[92]*z[125] - z[249]*z[66]*z[71] - z[250]*z[68]*z[71] - z[82]*z[281]*
  z[402] - z[83]*z[281]*z[403] - z[84]*z[282]*z[402] - z[85]*z[282]*z[403];
  z[465] = par__FootMass*(z[91]*z[396]-z[92]*z[397]-z[402]*z[99]-z[403]*
  z[100]) - par__ShankMass*(par__ThighLen*z[61]*z[76]+par__ThighLen*z[62]*
  z[72]-z[61]*z[71]*z[64]-z[62]*z[70]*z[64]);
  z[466] = 2*z[453]*(z[61]*z[70]-z[62]*z[71]) - 2*par__FootMass*(z[91]*z[402]+
  z[92]*z[403]);
  z[467] = z[455]*z[62] - par__FootMass*(z[90]*z[402]+z[93]*z[403]);
  z[468] = par__FootMass*(par__FootCMx*z[403]-par__FootCMy*z[402]);
  z[469] = z[82]*z[88]*z[281] + z[83]*z[278]*z[281] + z[83]*z[281]*z[412] + 
  z[84]*z[88]*z[282] + z[85]*z[278]*z[282] + z[85]*z[282]*z[412] + z[123]*
  z[279]*z[281] + z[126]*z[279]*z[282] - z[246]*z[82]*z[88] - z[246]*z[83]*
  z[91] - z[246]*z[83]*z[412] - z[246]*z[92]*z[123] - z[248]*z[84]*z[88] - 
  z[248]*z[85]*z[91] - z[248]*z[85]*z[412] - z[248]*z[92]*z[126];
  z[470] = par__FootMass*(z[88]*z[99]+z[91]*z[95]+z[92]*z[407]+z[412]*z[100]);
  z[471] = par__FootMass*(z[88]*z[91]+z[92]*z[412]);
  z[472] = par__FootMass*(z[88]*z[90]+z[93]*z[412]-par__ShankLen*z[77]*z[92]-
  par__ShankLen*z[78]*z[91]);
  z[473] = par__FootMass*(par__FootCMx*z[412]-par__FootCMy*z[88]);
  z[474] = z[34]*z[266]*z[271] + z[35]*z[269]*z[271] + z[106]*z[266]*z[270] + 
  z[108]*z[269]*z[270] + par__ShankCMy*(z[249]*z[102]+z[250]*z[17]) - par__ShankLen*
  z[17]*z[271] - par__ShankLen*z[102]*z[270] - par__ThighLen*z[7]*z[271] - 
  par__ThighLen*z[8]*z[270] - z[246]*z[42]*z[106] - z[246]*z[45]*z[108] - 
  z[248]*z[34]*z[42] - z[248]*z[35]*z[45];
  z[475] = z[107]*z[266]*z[270] + z[109]*z[269]*z[270] + z[111]*z[266]*z[271] + 
  z[112]*z[269]*z[271] + par__ShankCMy*(z[249]*z[18]+z[250]*z[20]) - par__ShankLen*
  z[18]*z[270] - par__ShankLen*z[20]*z[271] - z[246]*z[42]*z[107] - z[246]*
  z[45]*z[109] - z[248]*z[42]*z[111] - z[248]*z[45]*z[112];
  z[476] = par__FootMass*(z[42]*z[372]-z[45]*z[373]) - z[314]*z[13]*z[16];
  z[477] = par__ShankLen*z[246];
  z[478] = par__ShankLen*z[248];
  z[479] = z[477]*z[29]*z[35] + z[477]*z[30]*z[34] + z[478]*z[29]*z[37] + 
  z[478]*z[30]*z[36] + z[35]*z[266]*z[270] + z[37]*z[266]*z[271] + z[110]*
  z[269]*z[270] + z[113]*z[269]*z[271] - z[246]*z[35]*z[42] - z[246]*z[45]*
  z[110] - z[248]*z[37]*z[42] - z[248]*z[45]*z[113] - par__ShankLen*z[29]*
  z[35]*z[270] - par__ShankLen*z[29]*z[37]*z[271] - par__ShankLen*z[30]*z[34]*
  z[270] - par__ShankLen*z[30]*z[36]*z[271];
  z[480] = par__FootMass*(z[42]*z[47]+z[45]*z[383]-par__ShankLen*z[29]*z[52]-
  par__ShankLen*z[30]*z[51]);
  z[481] = par__FootMass*par__ShankLen;
  z[482] = z[481]*(z[29]*z[45]+z[30]*z[42]);
  z[483] = z[481]*(par__FootCMx*z[29]-par__FootCMy*z[30]);
  z[484] = z[299]*z[106] + z[300]*z[34] - z[301]*z[108] - z[302]*z[35] - 
  par__ShankLen*z[17]*z[271] - par__ShankLen*z[102]*z[270] - par__ThighLen*
  z[7]*z[271] - par__ThighLen*z[8]*z[270];
  z[485] = z[299]*z[107] + z[300]*z[111] - z[301]*z[109] - z[302]*z[112] - 
  par__ShankLen*z[18]*z[270] - par__ShankLen*z[20]*z[271];
  z[486] = par__FootMass*(par__FootCMx*z[373]+par__FootCMy*z[372]);
  z[487] = z[299]*z[35] + z[300]*z[37] - z[301]*z[110] - z[302]*z[113];
  z[488] = par__FootMass*(par__FootCMx*z[383]-par__FootCMy*z[47]);
  z[489] = z[82]*z[277]*z[282] + z[83]*z[280]*z[282] + z[119]*z[277]*z[281] + 
  z[121]*z[280]*z[281] + par__ShankCMy*(z[249]*z[115]+z[250]*z[65]) - par__ShankLen*
  z[65]*z[282] - par__ShankLen*z[115]*z[281] - par__ThighLen*z[55]*z[282] - 
  par__ThighLen*z[56]*z[281] - z[246]*z[90]*z[119] - z[246]*z[93]*z[121] - 
  z[248]*z[82]*z[90] - z[248]*z[83]*z[93];
  z[490] = z[120]*z[277]*z[281] + z[122]*z[280]*z[281] + z[124]*z[277]*z[282] + 
  z[125]*z[280]*z[282] + par__ShankCMy*(z[249]*z[66]+z[250]*z[68]) - par__ShankLen*
  z[66]*z[281] - par__ShankLen*z[68]*z[282] - z[246]*z[90]*z[120] - z[246]*
  z[93]*z[122] - z[248]*z[90]*z[124] - z[248]*z[93]*z[125];
  z[491] = par__FootMass*(z[90]*z[396]-z[93]*z[397]) - z[314]*z[61]*z[64];
  z[492] = z[477]*z[77]*z[83] + z[477]*z[78]*z[82] + z[478]*z[77]*z[85] + 
  z[478]*z[78]*z[84] + z[83]*z[277]*z[281] + z[85]*z[277]*z[282] + z[123]*
  z[280]*z[281] + z[126]*z[280]*z[282] - z[246]*z[83]*z[90] - z[246]*z[93]*
  z[123] - z[248]*z[85]*z[90] - z[248]*z[93]*z[126] - par__ShankLen*z[77]*
  z[83]*z[281] - par__ShankLen*z[77]*z[85]*z[282] - par__ShankLen*z[78]*z[82]*
  z[281] - par__ShankLen*z[78]*z[84]*z[282];
  z[493] = par__FootMass*(z[90]*z[95]+z[93]*z[407]-par__ShankLen*z[77]*z[100]-
  par__ShankLen*z[78]*z[99]);
  z[494] = z[481]*(z[77]*z[93]+z[78]*z[90]);
  z[495] = z[481]*(par__FootCMx*z[77]-par__FootCMy*z[78]);
  z[496] = z[299]*z[119] + z[300]*z[82] - z[301]*z[121] - z[302]*z[83] - 
  par__ShankLen*z[65]*z[282] - par__ShankLen*z[115]*z[281] - par__ThighLen*
  z[55]*z[282] - par__ThighLen*z[56]*z[281];
  z[497] = z[299]*z[120] + z[300]*z[124] - z[301]*z[122] - z[302]*z[125] - 
  par__ShankLen*z[66]*z[281] - par__ShankLen*z[68]*z[282];
  z[498] = par__FootMass*(par__FootCMx*z[397]+par__FootCMy*z[396]);
  z[499] = z[299]*z[83] + z[300]*z[85] - z[301]*z[123] - z[302]*z[126];
  z[500] = par__FootMass*(par__FootCMx*z[407]-par__FootCMy*z[95]);
  z[501] = 0.01*z[356] - 0.01*z[357] - 0.01*z[358]*q9pp - 0.01*z[359]*q1pp - 
  0.01*z[360]*q2pp - 0.01*z[361]*q3pp - 0.01*z[362]*q7pp - 0.01*z[363]*q8pp - 
  0.01*z[364]*q6pp - 0.01*z[365]*q4pp - 0.01*z[366]*q5pp;
  z[502] = 0.01*z[367] - 0.01*z[368] - 0.01*z[364]*q6pp - 0.01*z[365]*q3pp - 
  0.01*z[365]*q4pp - 0.01*z[366]*q5pp - 0.01*z[369]*q1pp - 0.01*z[370]*q2pp;
  z[503] = 0.01*z[371] - 0.01*z[374] - 0.01*z[375]*q1pp - 0.01*z[376]*q2pp - 
  0.01*z[379]*q3pp - 0.01*z[379]*q4pp - 0.01*z[380]*q6pp - 0.01*z[381]*q5pp;
  z[504] = 0.01*z[382] - 0.01*z[384] - 0.02*z[385]*q1pp - 0.01*z[386]*q2pp - 
  0.01*z[388]*q3pp - 0.01*z[388]*q4pp - 0.01*z[389]*q6pp - 0.01*z[390]*q5pp;
  z[505] = 0.01*z[391] - 0.01*z[392] - 0.01*z[358]*q9pp - 0.01*z[362]*q3pp - 
  0.01*z[362]*q7pp - 0.01*z[363]*q8pp - 0.01*z[393]*q1pp - 0.01*z[394]*q2pp;
  z[506] = 0.01*z[395] - 0.01*z[398] - 0.01*z[399]*q9pp - 0.01*z[400]*q1pp - 
  0.01*z[401]*q2pp - 0.01*z[404]*q3pp - 0.01*z[404]*q7pp - 0.01*z[405]*q8pp;
  z[507] = 0.01*z[406] - 0.01*z[408] - 0.02*z[410]*q1pp - 0.01*z[409]*q9pp - 
  0.01*z[411]*q2pp - 0.01*z[413]*q3pp - 0.01*z[413]*q7pp - 0.01*z[414]*q8pp;
  z[508] = 0.01*z[415] - 0.01*z[319] - 0.01*z[307]*q9pp - 0.01*z[312]*q3pp - 
  0.01*z[313]*q7pp - 0.01*z[315]*q8pp - 0.01*z[316]*q6pp - 0.01*z[317]*q4pp - 
  0.01*z[318]*q5pp - 0.01*z[360]*q1pp - 0.01*z[416]*q2pp;
  z[509] = 0.01*z[417] - 0.01*z[418] - 0.01*z[316]*q6pp - 0.01*z[317]*q3pp - 
  0.01*z[317]*q4pp - 0.01*z[318]*q5pp - 0.01*z[370]*q1pp - 0.01*z[419]*q2pp;
  z[510] = 0.01*z[420] - 0.01*z[421] - 0.01*z[376]*q1pp - 0.01*z[422]*q2pp - 
  0.01*z[423]*q3pp - 0.01*z[423]*q4pp - 0.01*z[424]*q6pp - 0.01*z[425]*q5pp;
  z[511] = 0.01*z[426] - 0.01*z[427] - 0.02*z[428]*q2pp - 0.01*z[386]*q1pp - 
  0.01*z[429]*q3pp - 0.01*z[429]*q4pp - 0.01*z[430]*q6pp - 0.01*z[431]*q5pp;
  z[512] = 0.01*z[432] - 0.01*z[433] - 0.01*z[307]*q9pp - 0.01*z[313]*q3pp - 
  0.01*z[313]*q7pp - 0.01*z[315]*q8pp - 0.01*z[394]*q1pp - 0.01*z[434]*q2pp;
  z[513] = 0.01*z[435] - 0.01*z[436] - 0.01*z[401]*q1pp - 0.01*z[437]*q9pp - 
  0.01*z[438]*q2pp - 0.01*z[439]*q3pp - 0.01*z[439]*q7pp - 0.01*z[440]*q8pp;
  z[514] = 0.01*z[441] - 0.01*z[442] - 0.02*z[444]*q2pp - 0.01*z[411]*q1pp - 
  0.01*z[443]*q9pp - 0.01*z[445]*q3pp - 0.01*z[445]*q7pp - 0.01*z[446]*q8pp;
  z[515] = -0.01*z[271] - 0.01*z[282];
  z[516] = 0.01*z[270] + 0.01*z[281];
  z[517] = 0.01*z[447] - 0.01*z[312]*q2pp - 0.01*z[361]*q1pp;
  z[518] = 0.01*z[448] - 0.01*z[317]*q2pp - 0.01*z[365]*q1pp;
  z[519] = 0.01*z[451] + 0.01*z[457]*q6pp - 0.01*z[452] - 0.01*z[379]*q1pp - 
  0.01*z[423]*q2pp - 0.01*z[454]*q3pp - 0.01*z[454]*q4pp - 0.01*z[456]*q5pp;
  z[520] = 0.01*z[458] - 0.01*z[459] - 0.02*z[460]*q3pp - 0.02*z[460]*q4pp - 
  0.01*z[388]*q1pp - 0.01*z[429]*q2pp - 0.01*z[461]*q5pp - 0.01*z[462]*q6pp;
  z[521] = 0.01*z[463] - 0.01*z[313]*q2pp - 0.01*z[362]*q1pp;
  z[522] = 0.01*z[464] + 0.01*z[468]*q9pp - 0.01*z[465] - 0.01*z[404]*q1pp - 
  0.01*z[439]*q2pp - 0.01*z[466]*q3pp - 0.01*z[466]*q7pp - 0.01*z[467]*q8pp;
  z[523] = 0.01*z[469] - 0.01*z[470] - 0.02*z[471]*q3pp - 0.02*z[471]*q7pp - 
  0.01*z[413]*q1pp - 0.01*z[445]*q2pp - 0.01*z[472]*q8pp - 0.01*z[473]*q9pp;
  z[524] = 0.01*z[474] - 0.01*z[318]*q2pp - 0.01*z[366]*q1pp;
  z[525] = 0.01*z[475] - 0.01*z[476] - 0.01*z[381]*q1pp - 0.01*z[425]*q2pp - 
  0.01*z[456]*q3pp - 0.01*z[456]*q4pp;
  z[526] = 0.01*z[479] + 0.01*z[483]*q6pp + 0.02*z[482]*q5pp - 0.01*z[480] - 
  0.01*z[390]*q1pp - 0.01*z[431]*q2pp - 0.01*z[461]*q3pp - 0.01*z[461]*q4pp;
  z[527] = 0.01*z[484] - 0.01*z[316]*q2pp - 0.01*z[364]*q1pp;
  z[528] = 0.01*z[485] + 0.01*z[486] + 0.01*z[457]*q3pp + 0.01*z[457]*q4pp - 
  0.01*z[380]*q1pp - 0.01*z[424]*q2pp;
  z[529] = 0.01*z[487] + 0.01*z[483]*q5pp - 0.01*z[488] - 0.01*z[389]*q1pp - 
  0.01*z[430]*q2pp - 0.01*z[462]*q3pp - 0.01*z[462]*q4pp;
  z[530] = 0.01*z[489] - 0.01*z[315]*q2pp - 0.01*z[363]*q1pp;
  z[531] = 0.01*z[490] - 0.01*z[491] - 0.01*z[405]*q1pp - 0.01*z[440]*q2pp - 
  0.01*z[467]*q3pp - 0.01*z[467]*q7pp;
  z[532] = 0.01*z[492] + 0.01*z[495]*q9pp + 0.02*z[494]*q8pp - 0.01*z[493] - 
  0.01*z[414]*q1pp - 0.01*z[446]*q2pp - 0.01*z[472]*q3pp - 0.01*z[472]*q7pp;
  z[533] = 0.01*z[496] - 0.01*z[307]*q2pp - 0.01*z[358]*q1pp;
  z[534] = 0.01*z[497] + 0.01*z[498] + 0.01*z[468]*q3pp + 0.01*z[468]*q7pp - 
  0.01*z[399]*q1pp - 0.01*z[437]*q2pp;
  z[535] = 0.01*z[499] + 0.01*z[495]*q8pp - 0.01*z[500] - 0.01*z[409]*q1pp - 
  0.01*z[443]*q2pp - 0.01*z[473]*q3pp - 0.01*z[473]*q7pp;
  z[536] = par__airdrag*(z[1]*(2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-
  pow(z[256],2)-2*par__wind*z[2]*z[256])+(par__TrunkCMy*q3p-z[1]*q1p-z[2]*q2p)*(
  z[1]*z[255]-par__wind*pow(z[1],2)-par__wind*pow(z[2],2)-z[2]*z[256]))/pow((
  z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*
  z[255]),0.5);
  z[537] = par__airdrag*(z[2]*(2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-
  pow(z[256],2)-2*par__wind*z[2]*z[256])+(z[1]*q2p-z[2]*q1p)*(z[1]*z[255]-
  par__wind*pow(z[1],2)-par__wind*pow(z[2],2)-z[2]*z[256]))/pow((z[258]+pow(
  z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),
  0.5);
  z[538] = z[257]*(z[1]*z[255]-par__wind*pow(z[1],2)-par__wind*pow(z[2],2)-
  z[2]*z[256])/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*
  z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[539] = z[1]*z[536] + z[2]*z[537] + pow(z[1],2)*z[538] + pow(z[2],2)*
  z[538];
  z[540] = par__airdrag*((par__TrunkCMy*q3p-z[1]*q1p-z[2]*q2p)*(z[1]*z[256]+
  z[2]*z[255])+z[2]*(2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-pow(z[256],
  2)-2*par__wind*z[2]*z[256]))/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*
  par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[541] = par__airdrag*((z[1]*q2p-z[2]*q1p)*(z[1]*z[256]+z[2]*z[255])-z[1]*(
  2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-pow(z[256],2)-2*par__wind*z[2]*
  z[256]))/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*
  par__wind*z[1]*z[255]),0.5);
  z[542] = z[257]*(z[1]*z[256]+z[2]*z[255])/pow((z[258]+pow(z[255],2)+pow(
  z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[543] = z[1]*z[540] + z[2]*z[541] + pow(z[1],2)*z[542] + pow(z[2],2)*
  z[542];
  z[544] = par__airdrag*par__TrunkCMy;
  z[545] = z[544]*(2*par__wind*z[1]*z[255]-z[258]-pow(z[255],2)-pow(z[256],2)-
  2*par__wind*z[2]*z[256]-(par__TrunkCMy*q3p-z[1]*q1p-z[2]*q2p)*(par__wind*
  z[1]-z[255]))/pow((z[258]+pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*
  z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[546] = z[544]*(z[1]*q2p-z[2]*q1p)*(par__wind*z[1]-z[255])/pow((z[258]+
  pow(z[255],2)+pow(z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),
  0.5);
  z[547] = par__TrunkCMy*z[257];
  z[548] = z[547]*(par__wind*z[1]-z[255])/pow((z[258]+pow(z[255],2)+pow(
  z[256],2)+2*par__wind*z[2]*z[256]-2*par__wind*z[1]*z[255]),0.5);
  z[549] = z[2]*z[546] + pow(z[1],2)*z[548] + pow(z[2],2)*z[548] - z[1]*
  z[545];
  z[550] = z[14]*(z[15]+par__ThighLen*z[11]);
  z[551] = z[13]*(z[15]+par__ThighLen*z[11]);
  z[552] = z[32] + z[551] + par__ShankLen*z[27];
  z[553] = z[29]*z[550] + z[30]*z[552];
  z[554] = z[553] - z[48] - par__FootCMx*z[50];
  z[555] = z[29]*z[552] - z[30]*z[550];
  z[556] = z[555] - z[49] - par__FootCMy*z[50];
  z[557] = z[62]*(z[63]+par__ThighLen*z[59]);
  z[558] = z[61]*(z[63]+par__ThighLen*z[59]);
  z[559] = z[80] + z[558] + par__ShankLen*z[75];
  z[560] = z[77]*z[557] + z[78]*z[559];
  z[561] = z[560] - z[96] - par__FootCMx*z[98];
  z[562] = z[77]*z[559] - z[78]*z[557];
  z[563] = z[562] - z[97] - par__FootCMy*z[98];
  z[564] = z[551] - z[26] - par__ShankCMy*z[27];
  z[565] = z[558] - z[74] - par__ShankCMy*z[75];
  z[566] = par__FootMass*(z[34]*z[554]+z[35]*z[556]) + par__FootMass*(z[82]*
  z[561]+z[83]*z[563]) + par__ShankMass*(z[17]*z[550]+z[18]*z[564]) + par__ShankMass*(
  z[65]*z[557]+z[66]*z[565]) + par__TrunkMass*z[2]*(par__TrunkCMy*q3p+z[3]) - 
  par__ThighMass*z[8]*(z[10]+par__ThighCMy*z[11]) - par__ThighMass*z[56]*(
  z[58]+par__ThighCMy*z[59]);
  z[567] = par__FootMass*(z[34]*z[554]+z[35]*z[556]) + par__ShankMass*(z[17]*
  z[550]+z[18]*z[564]) - par__ThighMass*z[8]*(z[10]+par__ThighCMy*z[11]);
  z[568] = z[30]*(z[32]+par__ShankLen*z[27]);
  z[569] = z[568] - z[48] - par__FootCMx*z[50];
  z[570] = z[29]*(z[32]+par__ShankLen*z[27]);
  z[571] = z[570] - z[49] - par__FootCMy*z[50];
  z[572] = par__FootMass*(z[34]*z[569]+z[35]*z[571]) - par__ShankMass*z[18]*(
  z[26]+par__ShankCMy*z[27]);
  z[573] = par__FootMass*(z[34]*(z[48]+par__FootCMx*z[50])+z[35]*(z[49]+
  par__FootCMy*z[50]));
  z[574] = par__FootMass*(z[82]*z[561]+z[83]*z[563]) + par__ShankMass*(z[65]*
  z[557]+z[66]*z[565]) - par__ThighMass*z[56]*(z[58]+par__ThighCMy*z[59]);
  z[575] = z[78]*(z[80]+par__ShankLen*z[75]);
  z[576] = z[575] - z[96] - par__FootCMx*z[98];
  z[577] = z[77]*(z[80]+par__ShankLen*z[75]);
  z[578] = z[577] - z[97] - par__FootCMy*z[98];
  z[579] = par__FootMass*(z[82]*z[576]+z[83]*z[578]) - par__ShankMass*z[66]*(
  z[74]+par__ShankCMy*z[75]);
  z[580] = par__FootMass*(z[82]*(z[96]+par__FootCMx*z[98])+z[83]*(z[97]+
  par__FootCMy*z[98]));
  z[581] = z[2]*z[536] - z[1]*z[537];
  z[582] = z[2]*z[540] - z[1]*z[541];
  z[583] = -z[1]*z[546] - z[2]*z[545];
  z[584] = par__FootMass*(z[36]*z[554]+z[37]*z[556]) + par__FootMass*(z[84]*
  z[561]+z[85]*z[563]) + par__ShankMass*(z[19]*z[550]+z[20]*z[564]) + par__ShankMass*(
  z[67]*z[557]+z[68]*z[565]) - par__ThighMass*z[7]*(z[10]+par__ThighCMy*z[11]) - 
  par__ThighMass*z[55]*(z[58]+par__ThighCMy*z[59]) - par__TrunkMass*z[1]*(
  par__TrunkCMy*q3p+z[3]);
  z[585] = par__FootMass*(z[36]*z[554]+z[37]*z[556]) + par__ShankMass*(z[19]*
  z[550]+z[20]*z[564]) - par__ThighMass*z[7]*(z[10]+par__ThighCMy*z[11]);
  z[586] = par__FootMass*(z[36]*z[569]+z[37]*z[571]) - par__ShankMass*z[20]*(
  z[26]+par__ShankCMy*z[27]);
  z[587] = par__FootMass*(z[36]*(z[48]+par__FootCMx*z[50])+z[37]*(z[49]+
  par__FootCMy*z[50]));
  z[588] = par__FootMass*(z[84]*z[561]+z[85]*z[563]) + par__ShankMass*(z[67]*
  z[557]+z[68]*z[565]) - par__ThighMass*z[55]*(z[58]+par__ThighCMy*z[59]);
  z[589] = par__FootMass*(z[84]*z[576]+z[85]*z[578]) - par__ShankMass*z[68]*(
  z[74]+par__ShankCMy*z[75]);
  z[590] = par__FootMass*(z[84]*(z[96]+par__FootCMx*z[98])+z[85]*(z[97]+
  par__FootCMy*z[98]));
  z[591] = par__TrunkCMy*(z[536]+z[1]*z[538]);
  z[592] = par__TrunkCMy*(z[540]+z[1]*z[542]);
  z[593] = par__TrunkCMy*(z[545]-z[1]*z[548]);
  z[594] = par__FootMass*(z[43]*z[554]+z[44]*z[556]) + par__FootMass*(z[91]*
  z[561]+z[92]*z[563]) - par__ShankMass*(z[22]*z[564]-z[23]*z[550]) - par__ShankMass*(
  z[70]*z[565]-z[71]*z[557]);
  z[595] = par__FootMass*(z[43]*z[554]+z[44]*z[556]) - par__ShankMass*(z[22]*
  z[564]-z[23]*z[550]);
  z[596] = par__FootMass*(z[43]*z[569]+z[44]*z[571]) + par__ShankMass*z[22]*(
  z[26]+par__ShankCMy*z[27]);
  z[597] = par__FootMass*(z[43]*(z[48]+par__FootCMx*z[50])+z[44]*(z[49]+
  par__FootCMy*z[50]));
  z[598] = par__FootMass*(z[91]*z[561]+z[92]*z[563]) - par__ShankMass*(z[70]*
  z[565]-z[71]*z[557]);
  z[599] = par__FootMass*(z[91]*z[576]+z[92]*z[578]) + par__ShankMass*z[70]*(
  z[74]+par__ShankCMy*z[75]);
  z[600] = par__FootMass*(z[91]*(z[96]+par__FootCMx*z[98])+z[92]*(z[97]+
  par__FootCMy*z[98]));
  z[601] = par__FootMass*(z[42]*z[554]+z[45]*z[556]) - z[314]*z[550];
  z[602] = par__FootMass*(z[42]*z[569]+z[45]*z[571]);
  z[603] = par__FootMass*(z[42]*(z[48]+par__FootCMx*z[50])+z[45]*(z[49]+
  par__FootCMy*z[50]));
  z[604] = par__FootMass*(par__FootCMx*z[556]-par__FootCMy*z[554]);
  z[605] = par__FootMass*(par__FootCMx*z[571]-par__FootCMy*z[569]);
  z[606] = par__FootMass*(par__FootCMx*z[49]-par__FootCMy*z[48]);
  z[607] = par__FootMass*(z[90]*z[561]+z[93]*z[563]) - z[314]*z[557];
  z[608] = par__FootMass*(z[90]*z[576]+z[93]*z[578]);
  z[609] = par__FootMass*(z[90]*(z[96]+par__FootCMx*z[98])+z[93]*(z[97]+
  par__FootCMy*z[98]));
  z[610] = par__FootMass*(par__FootCMx*z[563]-par__FootCMy*z[561]);
  z[611] = par__FootMass*(par__FootCMx*z[578]-par__FootCMy*z[576]);
  z[612] = par__FootMass*(par__FootCMx*z[97]-par__FootCMy*z[96]);
  z[613] = 0.01*z[549] - 0.01*z[566];
  z[614] = 0.01*z[583] - 0.01*z[584];
  z[615] = 0.01*z[593] - 0.01*z[594];
  z[616] = par__bodyweight*(z[34]*z[36]+z[35]*z[37]);
  z[617] = par__bodyweight*(z[82]*z[84]+z[83]*z[85]);
  z[618] = par__bodyweight*(par__ShankLen*z[17]+par__ThighLen*z[7]-q2-z[34]*
  z[267]-z[35]*z[268]);
  z[619] = par__bodyweight*(q1+par__ShankLen*z[19]+par__ThighLen*z[9]-z[36]*
  z[267]-z[37]*z[268]);
  z[620] = par__bodyweight*(par__ShankLen*z[65]+par__ThighLen*z[55]-q2-z[82]*
  z[278]-z[83]*z[279]);
  z[621] = par__bodyweight*(q1+par__ShankLen*z[67]+par__ThighLen*z[57]-z[84]*
  z[278]-z[85]*z[279]);
  z[622] = par__bodyweight*(par__ShankLen*z[17]+par__ThighLen*z[7]-q2-z[34]*
  z[266]-z[35]*z[269]);
  z[623] = par__bodyweight*(q1+par__ShankLen*z[19]+par__ThighLen*z[9]-z[36]*
  z[266]-z[37]*z[269]);
  z[624] = par__bodyweight*(q2-par__ShankLen*z[17]-par__ThighLen*z[7]);
  z[625] = par__bodyweight*(q1+par__ShankLen*z[19]+par__ThighLen*z[9]);
  z[626] = par__bodyweight*(par__ShankLen*z[65]+par__ThighLen*z[55]-q2-z[82]*
  z[277]-z[83]*z[280]);
  z[627] = par__bodyweight*(q1+par__ShankLen*z[67]+par__ThighLen*z[57]-z[84]*
  z[277]-z[85]*z[280]);
  z[628] = par__bodyweight*(q2-par__ShankLen*z[65]-par__ThighLen*z[55]);
  z[629] = par__bodyweight*(q1+par__ShankLen*z[67]+par__ThighLen*z[57]);
  z[630] = par__bodyweight*(pow(z[34],2)+pow(z[35],2));
  z[631] = par__bodyweight*(pow(z[82],2)+pow(z[83],2));
  z[632] = par__bodyweight*(pow(z[36],2)+pow(z[37],2));
  z[633] = par__bodyweight*(pow(z[84],2)+pow(z[85],2));


/* Evaluate output quantities */
}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[] )
{
int      i1;

/* Write output to screen and to output file(s) */
  fk[0] = q1;
  fk[1] = q2;
  fk[2] = z[1];
  fk[3] = -z[2];
  fk[4] = z[2];
  fk[5] = z[1];
  fk[6] = q1;
  fk[7] = q2;
  fk[8] = z[7];
  fk[9] = z[8];
  fk[10] = z[9];
  fk[11] = z[7];
  fk[12] = q1 - par__ThighLen*z[8];
  fk[13] = q2 - par__ThighLen*z[7];
  fk[14] = z[17];
  fk[15] = z[18];
  fk[16] = z[19];
  fk[17] = z[20];
  fk[18] = q1 - par__ShankLen*z[18] - par__ThighLen*z[8];
  fk[19] = q2 - par__ShankLen*z[20] - par__ThighLen*z[7];
  fk[20] = z[34];
  fk[21] = z[35];
  fk[22] = z[36];
  fk[23] = z[37];
  fk[24] = q1;
  fk[25] = q2;
  fk[26] = z[55];
  fk[27] = z[56];
  fk[28] = z[57];
  fk[29] = z[55];
  fk[30] = q1 - par__ThighLen*z[56];
  fk[31] = q2 - par__ThighLen*z[55];
  fk[32] = z[65];
  fk[33] = z[66];
  fk[34] = z[67];
  fk[35] = z[68];
  fk[36] = q1 - par__ShankLen*z[66] - par__ThighLen*z[56];
  fk[37] = q2 - par__ShankLen*z[68] - par__ThighLen*z[55];
  fk[38] = z[82];
  fk[39] = z[83];
  fk[40] = z[84];
  fk[41] = z[85];
  dfk_dq[0][0] = 1;
  dfk_dq[0][1] = 0;
  dfk_dq[0][2] = 0;
  dfk_dq[0][3] = 0;
  dfk_dq[0][4] = 0;
  dfk_dq[0][5] = 0;
  dfk_dq[0][6] = 0;
  dfk_dq[0][7] = 0;
  dfk_dq[0][8] = 0;
  dfk_dq[1][0] = 0;
  dfk_dq[1][1] = 1;
  dfk_dq[1][2] = 0;
  dfk_dq[1][3] = 0;
  dfk_dq[1][4] = 0;
  dfk_dq[1][5] = 0;
  dfk_dq[1][6] = 0;
  dfk_dq[1][7] = 0;
  dfk_dq[1][8] = 0;
  dfk_dq[2][0] = 0;
  dfk_dq[2][1] = 0;
  dfk_dq[2][2] = -z[2];
  dfk_dq[2][3] = 0;
  dfk_dq[2][4] = 0;
  dfk_dq[2][5] = 0;
  dfk_dq[2][6] = 0;
  dfk_dq[2][7] = 0;
  dfk_dq[2][8] = 0;
  dfk_dq[3][0] = 0;
  dfk_dq[3][1] = 0;
  dfk_dq[3][2] = -z[1];
  dfk_dq[3][3] = 0;
  dfk_dq[3][4] = 0;
  dfk_dq[3][5] = 0;
  dfk_dq[3][6] = 0;
  dfk_dq[3][7] = 0;
  dfk_dq[3][8] = 0;
  dfk_dq[4][0] = 0;
  dfk_dq[4][1] = 0;
  dfk_dq[4][2] = z[1];
  dfk_dq[4][3] = 0;
  dfk_dq[4][4] = 0;
  dfk_dq[4][5] = 0;
  dfk_dq[4][6] = 0;
  dfk_dq[4][7] = 0;
  dfk_dq[4][8] = 0;
  dfk_dq[5][0] = 0;
  dfk_dq[5][1] = 0;
  dfk_dq[5][2] = -z[2];
  dfk_dq[5][3] = 0;
  dfk_dq[5][4] = 0;
  dfk_dq[5][5] = 0;
  dfk_dq[5][6] = 0;
  dfk_dq[5][7] = 0;
  dfk_dq[5][8] = 0;
  dfk_dq[6][0] = 1;
  dfk_dq[6][1] = 0;
  dfk_dq[6][2] = 0;
  dfk_dq[6][3] = 0;
  dfk_dq[6][4] = 0;
  dfk_dq[6][5] = 0;
  dfk_dq[6][6] = 0;
  dfk_dq[6][7] = 0;
  dfk_dq[6][8] = 0;
  dfk_dq[7][0] = 0;
  dfk_dq[7][1] = 1;
  dfk_dq[7][2] = 0;
  dfk_dq[7][3] = 0;
  dfk_dq[7][4] = 0;
  dfk_dq[7][5] = 0;
  dfk_dq[7][6] = 0;
  dfk_dq[7][7] = 0;
  dfk_dq[7][8] = 0;
  dfk_dq[8][0] = 0;
  dfk_dq[8][1] = 0;
  dfk_dq[8][2] = z[8];
  dfk_dq[8][3] = z[8];
  dfk_dq[8][4] = 0;
  dfk_dq[8][5] = 0;
  dfk_dq[8][6] = 0;
  dfk_dq[8][7] = 0;
  dfk_dq[8][8] = 0;
  dfk_dq[9][0] = 0;
  dfk_dq[9][1] = 0;
  dfk_dq[9][2] = z[101];
  dfk_dq[9][3] = z[101];
  dfk_dq[9][4] = 0;
  dfk_dq[9][5] = 0;
  dfk_dq[9][6] = 0;
  dfk_dq[9][7] = 0;
  dfk_dq[9][8] = 0;
  dfk_dq[10][0] = 0;
  dfk_dq[10][1] = 0;
  dfk_dq[10][2] = z[7];
  dfk_dq[10][3] = z[7];
  dfk_dq[10][4] = 0;
  dfk_dq[10][5] = 0;
  dfk_dq[10][6] = 0;
  dfk_dq[10][7] = 0;
  dfk_dq[10][8] = 0;
  dfk_dq[11][0] = 0;
  dfk_dq[11][1] = 0;
  dfk_dq[11][2] = z[8];
  dfk_dq[11][3] = z[8];
  dfk_dq[11][4] = 0;
  dfk_dq[11][5] = 0;
  dfk_dq[11][6] = 0;
  dfk_dq[11][7] = 0;
  dfk_dq[11][8] = 0;
  dfk_dq[12][0] = 1;
  dfk_dq[12][1] = 0;
  dfk_dq[12][2] = -z[127];
  dfk_dq[12][3] = -z[127];
  dfk_dq[12][4] = 0;
  dfk_dq[12][5] = 0;
  dfk_dq[12][6] = 0;
  dfk_dq[12][7] = 0;
  dfk_dq[12][8] = 0;
  dfk_dq[13][0] = 0;
  dfk_dq[13][1] = 1;
  dfk_dq[13][2] = -z[128];
  dfk_dq[13][3] = -z[128];
  dfk_dq[13][4] = 0;
  dfk_dq[13][5] = 0;
  dfk_dq[13][6] = 0;
  dfk_dq[13][7] = 0;
  dfk_dq[13][8] = 0;
  dfk_dq[14][0] = 0;
  dfk_dq[14][1] = 0;
  dfk_dq[14][2] = z[102];
  dfk_dq[14][3] = z[102];
  dfk_dq[14][4] = z[18];
  dfk_dq[14][5] = 0;
  dfk_dq[14][6] = 0;
  dfk_dq[14][7] = 0;
  dfk_dq[14][8] = 0;
  dfk_dq[15][0] = 0;
  dfk_dq[15][1] = 0;
  dfk_dq[15][2] = z[103];
  dfk_dq[15][3] = z[103];
  dfk_dq[15][4] = z[104];
  dfk_dq[15][5] = 0;
  dfk_dq[15][6] = 0;
  dfk_dq[15][7] = 0;
  dfk_dq[15][8] = 0;
  dfk_dq[16][0] = 0;
  dfk_dq[16][1] = 0;
  dfk_dq[16][2] = z[17];
  dfk_dq[16][3] = z[17];
  dfk_dq[16][4] = z[20];
  dfk_dq[16][5] = 0;
  dfk_dq[16][6] = 0;
  dfk_dq[16][7] = 0;
  dfk_dq[16][8] = 0;
  dfk_dq[17][0] = 0;
  dfk_dq[17][1] = 0;
  dfk_dq[17][2] = z[18];
  dfk_dq[17][3] = z[18];
  dfk_dq[17][4] = z[105];
  dfk_dq[17][5] = 0;
  dfk_dq[17][6] = 0;
  dfk_dq[17][7] = 0;
  dfk_dq[17][8] = 0;
  dfk_dq[18][0] = 1;
  dfk_dq[18][1] = 0;
  dfk_dq[18][2] = z[129];
  dfk_dq[18][3] = z[129];
  dfk_dq[18][4] = -z[130];
  dfk_dq[18][5] = 0;
  dfk_dq[18][6] = 0;
  dfk_dq[18][7] = 0;
  dfk_dq[18][8] = 0;
  dfk_dq[19][0] = 0;
  dfk_dq[19][1] = 1;
  dfk_dq[19][2] = z[131];
  dfk_dq[19][3] = z[131];
  dfk_dq[19][4] = -z[132];
  dfk_dq[19][5] = 0;
  dfk_dq[19][6] = 0;
  dfk_dq[19][7] = 0;
  dfk_dq[19][8] = 0;
  dfk_dq[20][0] = 0;
  dfk_dq[20][1] = 0;
  dfk_dq[20][2] = z[106];
  dfk_dq[20][3] = z[106];
  dfk_dq[20][4] = z[107];
  dfk_dq[20][5] = z[35];
  dfk_dq[20][6] = 0;
  dfk_dq[20][7] = 0;
  dfk_dq[20][8] = 0;
  dfk_dq[21][0] = 0;
  dfk_dq[21][1] = 0;
  dfk_dq[21][2] = z[108];
  dfk_dq[21][3] = z[108];
  dfk_dq[21][4] = z[109];
  dfk_dq[21][5] = z[110];
  dfk_dq[21][6] = 0;
  dfk_dq[21][7] = 0;
  dfk_dq[21][8] = 0;
  dfk_dq[22][0] = 0;
  dfk_dq[22][1] = 0;
  dfk_dq[22][2] = z[34];
  dfk_dq[22][3] = z[34];
  dfk_dq[22][4] = z[111];
  dfk_dq[22][5] = z[37];
  dfk_dq[22][6] = 0;
  dfk_dq[22][7] = 0;
  dfk_dq[22][8] = 0;
  dfk_dq[23][0] = 0;
  dfk_dq[23][1] = 0;
  dfk_dq[23][2] = z[35];
  dfk_dq[23][3] = z[35];
  dfk_dq[23][4] = z[112];
  dfk_dq[23][5] = z[113];
  dfk_dq[23][6] = 0;
  dfk_dq[23][7] = 0;
  dfk_dq[23][8] = 0;
  dfk_dq[24][0] = 1;
  dfk_dq[24][1] = 0;
  dfk_dq[24][2] = 0;
  dfk_dq[24][3] = 0;
  dfk_dq[24][4] = 0;
  dfk_dq[24][5] = 0;
  dfk_dq[24][6] = 0;
  dfk_dq[24][7] = 0;
  dfk_dq[24][8] = 0;
  dfk_dq[25][0] = 0;
  dfk_dq[25][1] = 1;
  dfk_dq[25][2] = 0;
  dfk_dq[25][3] = 0;
  dfk_dq[25][4] = 0;
  dfk_dq[25][5] = 0;
  dfk_dq[25][6] = 0;
  dfk_dq[25][7] = 0;
  dfk_dq[25][8] = 0;
  dfk_dq[26][0] = 0;
  dfk_dq[26][1] = 0;
  dfk_dq[26][2] = z[56];
  dfk_dq[26][3] = 0;
  dfk_dq[26][4] = 0;
  dfk_dq[26][5] = 0;
  dfk_dq[26][6] = z[56];
  dfk_dq[26][7] = 0;
  dfk_dq[26][8] = 0;
  dfk_dq[27][0] = 0;
  dfk_dq[27][1] = 0;
  dfk_dq[27][2] = z[114];
  dfk_dq[27][3] = 0;
  dfk_dq[27][4] = 0;
  dfk_dq[27][5] = 0;
  dfk_dq[27][6] = z[114];
  dfk_dq[27][7] = 0;
  dfk_dq[27][8] = 0;
  dfk_dq[28][0] = 0;
  dfk_dq[28][1] = 0;
  dfk_dq[28][2] = z[55];
  dfk_dq[28][3] = 0;
  dfk_dq[28][4] = 0;
  dfk_dq[28][5] = 0;
  dfk_dq[28][6] = z[55];
  dfk_dq[28][7] = 0;
  dfk_dq[28][8] = 0;
  dfk_dq[29][0] = 0;
  dfk_dq[29][1] = 0;
  dfk_dq[29][2] = z[56];
  dfk_dq[29][3] = 0;
  dfk_dq[29][4] = 0;
  dfk_dq[29][5] = 0;
  dfk_dq[29][6] = z[56];
  dfk_dq[29][7] = 0;
  dfk_dq[29][8] = 0;
  dfk_dq[30][0] = 1;
  dfk_dq[30][1] = 0;
  dfk_dq[30][2] = -z[133];
  dfk_dq[30][3] = 0;
  dfk_dq[30][4] = 0;
  dfk_dq[30][5] = 0;
  dfk_dq[30][6] = -z[133];
  dfk_dq[30][7] = 0;
  dfk_dq[30][8] = 0;
  dfk_dq[31][0] = 0;
  dfk_dq[31][1] = 1;
  dfk_dq[31][2] = -z[134];
  dfk_dq[31][3] = 0;
  dfk_dq[31][4] = 0;
  dfk_dq[31][5] = 0;
  dfk_dq[31][6] = -z[134];
  dfk_dq[31][7] = 0;
  dfk_dq[31][8] = 0;
  dfk_dq[32][0] = 0;
  dfk_dq[32][1] = 0;
  dfk_dq[32][2] = z[115];
  dfk_dq[32][3] = 0;
  dfk_dq[32][4] = 0;
  dfk_dq[32][5] = 0;
  dfk_dq[32][6] = z[115];
  dfk_dq[32][7] = z[66];
  dfk_dq[32][8] = 0;
  dfk_dq[33][0] = 0;
  dfk_dq[33][1] = 0;
  dfk_dq[33][2] = z[116];
  dfk_dq[33][3] = 0;
  dfk_dq[33][4] = 0;
  dfk_dq[33][5] = 0;
  dfk_dq[33][6] = z[116];
  dfk_dq[33][7] = z[117];
  dfk_dq[33][8] = 0;
  dfk_dq[34][0] = 0;
  dfk_dq[34][1] = 0;
  dfk_dq[34][2] = z[65];
  dfk_dq[34][3] = 0;
  dfk_dq[34][4] = 0;
  dfk_dq[34][5] = 0;
  dfk_dq[34][6] = z[65];
  dfk_dq[34][7] = z[68];
  dfk_dq[34][8] = 0;
  dfk_dq[35][0] = 0;
  dfk_dq[35][1] = 0;
  dfk_dq[35][2] = z[66];
  dfk_dq[35][3] = 0;
  dfk_dq[35][4] = 0;
  dfk_dq[35][5] = 0;
  dfk_dq[35][6] = z[66];
  dfk_dq[35][7] = z[118];
  dfk_dq[35][8] = 0;
  dfk_dq[36][0] = 1;
  dfk_dq[36][1] = 0;
  dfk_dq[36][2] = z[135];
  dfk_dq[36][3] = 0;
  dfk_dq[36][4] = 0;
  dfk_dq[36][5] = 0;
  dfk_dq[36][6] = z[135];
  dfk_dq[36][7] = -z[136];
  dfk_dq[36][8] = 0;
  dfk_dq[37][0] = 0;
  dfk_dq[37][1] = 1;
  dfk_dq[37][2] = z[137];
  dfk_dq[37][3] = 0;
  dfk_dq[37][4] = 0;
  dfk_dq[37][5] = 0;
  dfk_dq[37][6] = z[137];
  dfk_dq[37][7] = -z[138];
  dfk_dq[37][8] = 0;
  dfk_dq[38][0] = 0;
  dfk_dq[38][1] = 0;
  dfk_dq[38][2] = z[119];
  dfk_dq[38][3] = 0;
  dfk_dq[38][4] = 0;
  dfk_dq[38][5] = 0;
  dfk_dq[38][6] = z[119];
  dfk_dq[38][7] = z[120];
  dfk_dq[38][8] = z[83];
  dfk_dq[39][0] = 0;
  dfk_dq[39][1] = 0;
  dfk_dq[39][2] = z[121];
  dfk_dq[39][3] = 0;
  dfk_dq[39][4] = 0;
  dfk_dq[39][5] = 0;
  dfk_dq[39][6] = z[121];
  dfk_dq[39][7] = z[122];
  dfk_dq[39][8] = z[123];
  dfk_dq[40][0] = 0;
  dfk_dq[40][1] = 0;
  dfk_dq[40][2] = z[82];
  dfk_dq[40][3] = 0;
  dfk_dq[40][4] = 0;
  dfk_dq[40][5] = 0;
  dfk_dq[40][6] = z[82];
  dfk_dq[40][7] = z[124];
  dfk_dq[40][8] = z[85];
  dfk_dq[41][0] = 0;
  dfk_dq[41][1] = 0;
  dfk_dq[41][2] = z[83];
  dfk_dq[41][3] = 0;
  dfk_dq[41][4] = 0;
  dfk_dq[41][5] = 0;
  dfk_dq[41][6] = z[83];
  dfk_dq[41][7] = z[125];
  dfk_dq[41][8] = z[126];
  fkdot[0] = q1p;
  fkdot[1] = q2p;
  fkdot[2] = -z[139];
  fkdot[3] = -z[140];
  fkdot[4] = z[140];
  fkdot[5] = -z[139];
  fkdot[6] = q1p;
  fkdot[7] = q2p;
  fkdot[8] = z[141];
  fkdot[9] = z[142];
  fkdot[10] = z[143];
  fkdot[11] = z[141];
  fkdot[12] = z[144];
  fkdot[13] = z[145];
  fkdot[14] = z[146];
  fkdot[15] = z[147];
  fkdot[16] = z[148];
  fkdot[17] = z[149];
  fkdot[18] = z[150];
  fkdot[19] = z[151];
  fkdot[20] = z[152];
  fkdot[21] = z[153];
  fkdot[22] = z[154];
  fkdot[23] = z[155];
  fkdot[24] = q1p;
  fkdot[25] = q2p;
  fkdot[26] = z[156];
  fkdot[27] = z[157];
  fkdot[28] = z[158];
  fkdot[29] = z[156];
  fkdot[30] = z[159];
  fkdot[31] = z[160];
  fkdot[32] = z[161];
  fkdot[33] = z[162];
  fkdot[34] = z[163];
  fkdot[35] = z[164];
  fkdot[36] = z[165];
  fkdot[37] = z[166];
  fkdot[38] = z[167];
  fkdot[39] = z[168];
  fkdot[40] = z[169];
  fkdot[41] = z[170];
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
  dfkdot_dq[2][0] = 0;
  dfkdot_dq[2][1] = 0;
  dfkdot_dq[2][2] = -z[140];
  dfkdot_dq[2][3] = 0;
  dfkdot_dq[2][4] = 0;
  dfkdot_dq[2][5] = 0;
  dfkdot_dq[2][6] = 0;
  dfkdot_dq[2][7] = 0;
  dfkdot_dq[2][8] = 0;
  dfkdot_dq[3][0] = 0;
  dfkdot_dq[3][1] = 0;
  dfkdot_dq[3][2] = z[139];
  dfkdot_dq[3][3] = 0;
  dfkdot_dq[3][4] = 0;
  dfkdot_dq[3][5] = 0;
  dfkdot_dq[3][6] = 0;
  dfkdot_dq[3][7] = 0;
  dfkdot_dq[3][8] = 0;
  dfkdot_dq[4][0] = 0;
  dfkdot_dq[4][1] = 0;
  dfkdot_dq[4][2] = -z[139];
  dfkdot_dq[4][3] = 0;
  dfkdot_dq[4][4] = 0;
  dfkdot_dq[4][5] = 0;
  dfkdot_dq[4][6] = 0;
  dfkdot_dq[4][7] = 0;
  dfkdot_dq[4][8] = 0;
  dfkdot_dq[5][0] = 0;
  dfkdot_dq[5][1] = 0;
  dfkdot_dq[5][2] = -z[140];
  dfkdot_dq[5][3] = 0;
  dfkdot_dq[5][4] = 0;
  dfkdot_dq[5][5] = 0;
  dfkdot_dq[5][6] = 0;
  dfkdot_dq[5][7] = 0;
  dfkdot_dq[5][8] = 0;
  dfkdot_dq[6][0] = 0;
  dfkdot_dq[6][1] = 0;
  dfkdot_dq[6][2] = 0;
  dfkdot_dq[6][3] = 0;
  dfkdot_dq[6][4] = 0;
  dfkdot_dq[6][5] = 0;
  dfkdot_dq[6][6] = 0;
  dfkdot_dq[6][7] = 0;
  dfkdot_dq[6][8] = 0;
  dfkdot_dq[7][0] = 0;
  dfkdot_dq[7][1] = 0;
  dfkdot_dq[7][2] = 0;
  dfkdot_dq[7][3] = 0;
  dfkdot_dq[7][4] = 0;
  dfkdot_dq[7][5] = 0;
  dfkdot_dq[7][6] = 0;
  dfkdot_dq[7][7] = 0;
  dfkdot_dq[7][8] = 0;
  dfkdot_dq[8][0] = 0;
  dfkdot_dq[8][1] = 0;
  dfkdot_dq[8][2] = z[142];
  dfkdot_dq[8][3] = z[142];
  dfkdot_dq[8][4] = 0;
  dfkdot_dq[8][5] = 0;
  dfkdot_dq[8][6] = 0;
  dfkdot_dq[8][7] = 0;
  dfkdot_dq[8][8] = 0;
  dfkdot_dq[9][0] = 0;
  dfkdot_dq[9][1] = 0;
  dfkdot_dq[9][2] = z[243];
  dfkdot_dq[9][3] = z[243];
  dfkdot_dq[9][4] = 0;
  dfkdot_dq[9][5] = 0;
  dfkdot_dq[9][6] = 0;
  dfkdot_dq[9][7] = 0;
  dfkdot_dq[9][8] = 0;
  dfkdot_dq[10][0] = 0;
  dfkdot_dq[10][1] = 0;
  dfkdot_dq[10][2] = z[141];
  dfkdot_dq[10][3] = z[141];
  dfkdot_dq[10][4] = 0;
  dfkdot_dq[10][5] = 0;
  dfkdot_dq[10][6] = 0;
  dfkdot_dq[10][7] = 0;
  dfkdot_dq[10][8] = 0;
  dfkdot_dq[11][0] = 0;
  dfkdot_dq[11][1] = 0;
  dfkdot_dq[11][2] = z[142];
  dfkdot_dq[11][3] = z[142];
  dfkdot_dq[11][4] = 0;
  dfkdot_dq[11][5] = 0;
  dfkdot_dq[11][6] = 0;
  dfkdot_dq[11][7] = 0;
  dfkdot_dq[11][8] = 0;
  dfkdot_dq[12][0] = 0;
  dfkdot_dq[12][1] = 0;
  dfkdot_dq[12][2] = -z[171];
  dfkdot_dq[12][3] = -z[171];
  dfkdot_dq[12][4] = 0;
  dfkdot_dq[12][5] = 0;
  dfkdot_dq[12][6] = 0;
  dfkdot_dq[12][7] = 0;
  dfkdot_dq[12][8] = 0;
  dfkdot_dq[13][0] = 0;
  dfkdot_dq[13][1] = 0;
  dfkdot_dq[13][2] = -z[172];
  dfkdot_dq[13][3] = -z[172];
  dfkdot_dq[13][4] = 0;
  dfkdot_dq[13][5] = 0;
  dfkdot_dq[13][6] = 0;
  dfkdot_dq[13][7] = 0;
  dfkdot_dq[13][8] = 0;
  dfkdot_dq[14][0] = 0;
  dfkdot_dq[14][1] = 0;
  dfkdot_dq[14][2] = z[174];
  dfkdot_dq[14][3] = z[174];
  dfkdot_dq[14][4] = z[147];
  dfkdot_dq[14][5] = 0;
  dfkdot_dq[14][6] = 0;
  dfkdot_dq[14][7] = 0;
  dfkdot_dq[14][8] = 0;
  dfkdot_dq[15][0] = 0;
  dfkdot_dq[15][1] = 0;
  dfkdot_dq[15][2] = z[177];
  dfkdot_dq[15][3] = z[177];
  dfkdot_dq[15][4] = z[179];
  dfkdot_dq[15][5] = 0;
  dfkdot_dq[15][6] = 0;
  dfkdot_dq[15][7] = 0;
  dfkdot_dq[15][8] = 0;
  dfkdot_dq[16][0] = 0;
  dfkdot_dq[16][1] = 0;
  dfkdot_dq[16][2] = z[146];
  dfkdot_dq[16][3] = z[146];
  dfkdot_dq[16][4] = z[149];
  dfkdot_dq[16][5] = 0;
  dfkdot_dq[16][6] = 0;
  dfkdot_dq[16][7] = 0;
  dfkdot_dq[16][8] = 0;
  dfkdot_dq[17][0] = 0;
  dfkdot_dq[17][1] = 0;
  dfkdot_dq[17][2] = z[147];
  dfkdot_dq[17][3] = z[147];
  dfkdot_dq[17][4] = z[181];
  dfkdot_dq[17][5] = 0;
  dfkdot_dq[17][6] = 0;
  dfkdot_dq[17][7] = 0;
  dfkdot_dq[17][8] = 0;
  dfkdot_dq[18][0] = 0;
  dfkdot_dq[18][1] = 0;
  dfkdot_dq[18][2] = z[182];
  dfkdot_dq[18][3] = z[182];
  dfkdot_dq[18][4] = -z[183];
  dfkdot_dq[18][5] = 0;
  dfkdot_dq[18][6] = 0;
  dfkdot_dq[18][7] = 0;
  dfkdot_dq[18][8] = 0;
  dfkdot_dq[19][0] = 0;
  dfkdot_dq[19][1] = 0;
  dfkdot_dq[19][2] = z[184];
  dfkdot_dq[19][3] = z[184];
  dfkdot_dq[19][4] = -z[185];
  dfkdot_dq[19][5] = 0;
  dfkdot_dq[19][6] = 0;
  dfkdot_dq[19][7] = 0;
  dfkdot_dq[19][8] = 0;
  dfkdot_dq[20][0] = 0;
  dfkdot_dq[20][1] = 0;
  dfkdot_dq[20][2] = z[188];
  dfkdot_dq[20][3] = z[188];
  dfkdot_dq[20][4] = z[190];
  dfkdot_dq[20][5] = z[153];
  dfkdot_dq[20][6] = 0;
  dfkdot_dq[20][7] = 0;
  dfkdot_dq[20][8] = 0;
  dfkdot_dq[21][0] = 0;
  dfkdot_dq[21][1] = 0;
  dfkdot_dq[21][2] = z[194];
  dfkdot_dq[21][3] = z[194];
  dfkdot_dq[21][4] = z[197];
  dfkdot_dq[21][5] = z[199];
  dfkdot_dq[21][6] = 0;
  dfkdot_dq[21][7] = 0;
  dfkdot_dq[21][8] = 0;
  dfkdot_dq[22][0] = 0;
  dfkdot_dq[22][1] = 0;
  dfkdot_dq[22][2] = z[152];
  dfkdot_dq[22][3] = z[152];
  dfkdot_dq[22][4] = z[201];
  dfkdot_dq[22][5] = z[155];
  dfkdot_dq[22][6] = 0;
  dfkdot_dq[22][7] = 0;
  dfkdot_dq[22][8] = 0;
  dfkdot_dq[23][0] = 0;
  dfkdot_dq[23][1] = 0;
  dfkdot_dq[23][2] = z[153];
  dfkdot_dq[23][3] = z[153];
  dfkdot_dq[23][4] = z[204];
  dfkdot_dq[23][5] = z[206];
  dfkdot_dq[23][6] = 0;
  dfkdot_dq[23][7] = 0;
  dfkdot_dq[23][8] = 0;
  dfkdot_dq[24][0] = 0;
  dfkdot_dq[24][1] = 0;
  dfkdot_dq[24][2] = 0;
  dfkdot_dq[24][3] = 0;
  dfkdot_dq[24][4] = 0;
  dfkdot_dq[24][5] = 0;
  dfkdot_dq[24][6] = 0;
  dfkdot_dq[24][7] = 0;
  dfkdot_dq[24][8] = 0;
  dfkdot_dq[25][0] = 0;
  dfkdot_dq[25][1] = 0;
  dfkdot_dq[25][2] = 0;
  dfkdot_dq[25][3] = 0;
  dfkdot_dq[25][4] = 0;
  dfkdot_dq[25][5] = 0;
  dfkdot_dq[25][6] = 0;
  dfkdot_dq[25][7] = 0;
  dfkdot_dq[25][8] = 0;
  dfkdot_dq[26][0] = 0;
  dfkdot_dq[26][1] = 0;
  dfkdot_dq[26][2] = z[157];
  dfkdot_dq[26][3] = 0;
  dfkdot_dq[26][4] = 0;
  dfkdot_dq[26][5] = 0;
  dfkdot_dq[26][6] = z[157];
  dfkdot_dq[26][7] = 0;
  dfkdot_dq[26][8] = 0;
  dfkdot_dq[27][0] = 0;
  dfkdot_dq[27][1] = 0;
  dfkdot_dq[27][2] = z[244];
  dfkdot_dq[27][3] = 0;
  dfkdot_dq[27][4] = 0;
  dfkdot_dq[27][5] = 0;
  dfkdot_dq[27][6] = z[244];
  dfkdot_dq[27][7] = 0;
  dfkdot_dq[27][8] = 0;
  dfkdot_dq[28][0] = 0;
  dfkdot_dq[28][1] = 0;
  dfkdot_dq[28][2] = z[156];
  dfkdot_dq[28][3] = 0;
  dfkdot_dq[28][4] = 0;
  dfkdot_dq[28][5] = 0;
  dfkdot_dq[28][6] = z[156];
  dfkdot_dq[28][7] = 0;
  dfkdot_dq[28][8] = 0;
  dfkdot_dq[29][0] = 0;
  dfkdot_dq[29][1] = 0;
  dfkdot_dq[29][2] = z[157];
  dfkdot_dq[29][3] = 0;
  dfkdot_dq[29][4] = 0;
  dfkdot_dq[29][5] = 0;
  dfkdot_dq[29][6] = z[157];
  dfkdot_dq[29][7] = 0;
  dfkdot_dq[29][8] = 0;
  dfkdot_dq[30][0] = 0;
  dfkdot_dq[30][1] = 0;
  dfkdot_dq[30][2] = -z[207];
  dfkdot_dq[30][3] = 0;
  dfkdot_dq[30][4] = 0;
  dfkdot_dq[30][5] = 0;
  dfkdot_dq[30][6] = -z[207];
  dfkdot_dq[30][7] = 0;
  dfkdot_dq[30][8] = 0;
  dfkdot_dq[31][0] = 0;
  dfkdot_dq[31][1] = 0;
  dfkdot_dq[31][2] = -z[208];
  dfkdot_dq[31][3] = 0;
  dfkdot_dq[31][4] = 0;
  dfkdot_dq[31][5] = 0;
  dfkdot_dq[31][6] = -z[208];
  dfkdot_dq[31][7] = 0;
  dfkdot_dq[31][8] = 0;
  dfkdot_dq[32][0] = 0;
  dfkdot_dq[32][1] = 0;
  dfkdot_dq[32][2] = z[210];
  dfkdot_dq[32][3] = 0;
  dfkdot_dq[32][4] = 0;
  dfkdot_dq[32][5] = 0;
  dfkdot_dq[32][6] = z[210];
  dfkdot_dq[32][7] = z[162];
  dfkdot_dq[32][8] = 0;
  dfkdot_dq[33][0] = 0;
  dfkdot_dq[33][1] = 0;
  dfkdot_dq[33][2] = z[213];
  dfkdot_dq[33][3] = 0;
  dfkdot_dq[33][4] = 0;
  dfkdot_dq[33][5] = 0;
  dfkdot_dq[33][6] = z[213];
  dfkdot_dq[33][7] = z[215];
  dfkdot_dq[33][8] = 0;
  dfkdot_dq[34][0] = 0;
  dfkdot_dq[34][1] = 0;
  dfkdot_dq[34][2] = z[161];
  dfkdot_dq[34][3] = 0;
  dfkdot_dq[34][4] = 0;
  dfkdot_dq[34][5] = 0;
  dfkdot_dq[34][6] = z[161];
  dfkdot_dq[34][7] = z[164];
  dfkdot_dq[34][8] = 0;
  dfkdot_dq[35][0] = 0;
  dfkdot_dq[35][1] = 0;
  dfkdot_dq[35][2] = z[162];
  dfkdot_dq[35][3] = 0;
  dfkdot_dq[35][4] = 0;
  dfkdot_dq[35][5] = 0;
  dfkdot_dq[35][6] = z[162];
  dfkdot_dq[35][7] = z[217];
  dfkdot_dq[35][8] = 0;
  dfkdot_dq[36][0] = 0;
  dfkdot_dq[36][1] = 0;
  dfkdot_dq[36][2] = z[218];
  dfkdot_dq[36][3] = 0;
  dfkdot_dq[36][4] = 0;
  dfkdot_dq[36][5] = 0;
  dfkdot_dq[36][6] = z[218];
  dfkdot_dq[36][7] = -z[219];
  dfkdot_dq[36][8] = 0;
  dfkdot_dq[37][0] = 0;
  dfkdot_dq[37][1] = 0;
  dfkdot_dq[37][2] = z[220];
  dfkdot_dq[37][3] = 0;
  dfkdot_dq[37][4] = 0;
  dfkdot_dq[37][5] = 0;
  dfkdot_dq[37][6] = z[220];
  dfkdot_dq[37][7] = -z[221];
  dfkdot_dq[37][8] = 0;
  dfkdot_dq[38][0] = 0;
  dfkdot_dq[38][1] = 0;
  dfkdot_dq[38][2] = z[224];
  dfkdot_dq[38][3] = 0;
  dfkdot_dq[38][4] = 0;
  dfkdot_dq[38][5] = 0;
  dfkdot_dq[38][6] = z[224];
  dfkdot_dq[38][7] = z[226];
  dfkdot_dq[38][8] = z[168];
  dfkdot_dq[39][0] = 0;
  dfkdot_dq[39][1] = 0;
  dfkdot_dq[39][2] = z[230];
  dfkdot_dq[39][3] = 0;
  dfkdot_dq[39][4] = 0;
  dfkdot_dq[39][5] = 0;
  dfkdot_dq[39][6] = z[230];
  dfkdot_dq[39][7] = z[233];
  dfkdot_dq[39][8] = z[235];
  dfkdot_dq[40][0] = 0;
  dfkdot_dq[40][1] = 0;
  dfkdot_dq[40][2] = z[167];
  dfkdot_dq[40][3] = 0;
  dfkdot_dq[40][4] = 0;
  dfkdot_dq[40][5] = 0;
  dfkdot_dq[40][6] = z[167];
  dfkdot_dq[40][7] = z[237];
  dfkdot_dq[40][8] = z[170];
  dfkdot_dq[41][0] = 0;
  dfkdot_dq[41][1] = 0;
  dfkdot_dq[41][2] = z[168];
  dfkdot_dq[41][3] = 0;
  dfkdot_dq[41][4] = 0;
  dfkdot_dq[41][5] = 0;
  dfkdot_dq[41][6] = z[168];
  dfkdot_dq[41][7] = z[240];
  dfkdot_dq[41][8] = z[242];
  Zero[0] = 0.01*z[294] - 0.01*z[319] - 0.01*z[307]*q9pp - 0.01*z[308]*q1pp - 
  0.01*z[309]*q2pp - 0.01*z[312]*q3pp - 0.01*z[313]*q7pp - 0.01*z[315]*q8pp - 
  0.01*z[316]*q6pp - 0.01*z[317]*q4pp - 0.01*z[318]*q5pp;
  Zero[1] = 0.01*z[295] - 0.01*z[328] - 0.01*z[309]*q1pp - 0.01*z[320]*q9pp - 
  0.01*z[321]*q2pp - 0.01*z[322]*q3pp - 0.01*z[323]*q7pp - 0.01*z[324]*q8pp - 
  0.01*z[325]*q6pp - 0.01*z[326]*q4pp - 0.01*z[327]*q5pp;
  Zero[2] = 0.01*z[296] - 0.01*z[339] - 0.01*z[312]*q1pp - 0.01*z[322]*q2pp - 
  0.01*z[330]*q3pp - 0.01*z[332]*q7pp - 0.01*z[334]*q8pp - 0.01*z[335]*q9pp - 
  0.01*z[336]*q4pp - 0.01*z[337]*q5pp - 0.01*z[338]*q6pp;
  Zero[3] = 0.01*z[297] - 0.01*z[340] - 0.01*z[317]*q1pp - 0.01*z[326]*q2pp - 
  0.01*z[336]*q3pp - 0.01*z[336]*q4pp - 0.01*z[337]*q5pp - 0.01*z[338]*q6pp;
  Zero[4] = 0.01*z[298] - 0.01*z[344] - 0.01*z[318]*q1pp - 0.01*z[327]*q2pp - 
  0.01*z[337]*q3pp - 0.01*z[337]*q4pp - 0.01*z[342]*q5pp - 0.01*z[343]*q6pp;
  Zero[5] = 0.01*z[303] - 0.01*z[346] - 0.01*z[345]*q6pp - 0.01*z[316]*q1pp - 
  0.01*z[325]*q2pp - 0.01*z[338]*q3pp - 0.01*z[338]*q4pp - 0.01*z[343]*q5pp;
  Zero[6] = 0.01*z[304] - 0.01*z[347] - 0.01*z[313]*q1pp - 0.01*z[323]*q2pp - 
  0.01*z[332]*q3pp - 0.01*z[332]*q7pp - 0.01*z[334]*q8pp - 0.01*z[335]*q9pp;
  Zero[7] = 0.01*z[305] - 0.01*z[350] - 0.01*z[315]*q1pp - 0.01*z[324]*q2pp - 
  0.01*z[334]*q3pp - 0.01*z[334]*q7pp - 0.01*z[348]*q8pp - 0.01*z[349]*q9pp;
  Zero[8] = 0.01*z[306] - 0.01*z[351] - 0.01*z[345]*q9pp - 0.01*z[307]*q1pp - 
  0.01*z[320]*q2pp - 0.01*z[335]*q3pp - 0.01*z[335]*q7pp - 0.01*z[349]*q8pp;
  dz_dq[0][0] = 0;
  dz_dq[0][1] = 0;
  dz_dq[0][2] = z[501];
  dz_dq[0][3] = z[502];
  dz_dq[0][4] = z[503];
  dz_dq[0][5] = z[504];
  dz_dq[0][6] = z[505];
  dz_dq[0][7] = z[506];
  dz_dq[0][8] = z[507];
  dz_dq[1][0] = 0;
  dz_dq[1][1] = 0;
  dz_dq[1][2] = z[508];
  dz_dq[1][3] = z[509];
  dz_dq[1][4] = z[510];
  dz_dq[1][5] = z[511];
  dz_dq[1][6] = z[512];
  dz_dq[1][7] = z[513];
  dz_dq[1][8] = z[514];
  dz_dq[2][0] = z[515];
  dz_dq[2][1] = z[516];
  dz_dq[2][2] = z[517];
  dz_dq[2][3] = z[518];
  dz_dq[2][4] = z[519];
  dz_dq[2][5] = z[520];
  dz_dq[2][6] = z[521];
  dz_dq[2][7] = z[522];
  dz_dq[2][8] = z[523];
  dz_dq[3][0] = -0.01*z[271];
  dz_dq[3][1] = 0.01*z[270];
  dz_dq[3][2] = z[518];
  dz_dq[3][3] = z[518];
  dz_dq[3][4] = z[519];
  dz_dq[3][5] = z[520];
  dz_dq[3][6] = 0;
  dz_dq[3][7] = 0;
  dz_dq[3][8] = 0;
  dz_dq[4][0] = -0.01*z[271];
  dz_dq[4][1] = 0.01*z[270];
  dz_dq[4][2] = z[524];
  dz_dq[4][3] = z[524];
  dz_dq[4][4] = z[525];
  dz_dq[4][5] = z[526];
  dz_dq[4][6] = 0;
  dz_dq[4][7] = 0;
  dz_dq[4][8] = 0;
  dz_dq[5][0] = -0.01*z[271];
  dz_dq[5][1] = 0.01*z[270];
  dz_dq[5][2] = z[527];
  dz_dq[5][3] = z[527];
  dz_dq[5][4] = z[528];
  dz_dq[5][5] = z[529];
  dz_dq[5][6] = 0;
  dz_dq[5][7] = 0;
  dz_dq[5][8] = 0;
  dz_dq[6][0] = -0.01*z[282];
  dz_dq[6][1] = 0.01*z[281];
  dz_dq[6][2] = z[521];
  dz_dq[6][3] = 0;
  dz_dq[6][4] = 0;
  dz_dq[6][5] = 0;
  dz_dq[6][6] = z[521];
  dz_dq[6][7] = z[522];
  dz_dq[6][8] = z[523];
  dz_dq[7][0] = -0.01*z[282];
  dz_dq[7][1] = 0.01*z[281];
  dz_dq[7][2] = z[530];
  dz_dq[7][3] = 0;
  dz_dq[7][4] = 0;
  dz_dq[7][5] = 0;
  dz_dq[7][6] = z[530];
  dz_dq[7][7] = z[531];
  dz_dq[7][8] = z[532];
  dz_dq[8][0] = -0.01*z[282];
  dz_dq[8][1] = 0.01*z[281];
  dz_dq[8][2] = z[533];
  dz_dq[8][3] = 0;
  dz_dq[8][4] = 0;
  dz_dq[8][5] = 0;
  dz_dq[8][6] = z[533];
  dz_dq[8][7] = z[534];
  dz_dq[8][8] = z[535];
  dz_dqd[0][0] = 0.01*z[539];
  dz_dqd[0][1] = 0.01*z[543];
  dz_dqd[0][2] = z[613];
  dz_dqd[0][3] = -0.01*z[567];
  dz_dqd[0][4] = -0.01*z[572];
  dz_dqd[0][5] = 0.01*z[573];
  dz_dqd[0][6] = -0.01*z[574];
  dz_dqd[0][7] = -0.01*z[579];
  dz_dqd[0][8] = 0.01*z[580];
  dz_dqd[1][0] = 0.01*z[581];
  dz_dqd[1][1] = 0.01*z[582];
  dz_dqd[1][2] = z[614];
  dz_dqd[1][3] = -0.01*z[585];
  dz_dqd[1][4] = -0.01*z[586];
  dz_dqd[1][5] = 0.01*z[587];
  dz_dqd[1][6] = -0.01*z[588];
  dz_dqd[1][7] = -0.01*z[589];
  dz_dqd[1][8] = 0.01*z[590];
  dz_dqd[2][0] = -0.01*z[591];
  dz_dqd[2][1] = -0.01*z[592];
  dz_dqd[2][2] = z[615];
  dz_dqd[2][3] = -0.01*z[595];
  dz_dqd[2][4] = -0.01*z[596];
  dz_dqd[2][5] = 0.01*z[597];
  dz_dqd[2][6] = -0.01*z[598];
  dz_dqd[2][7] = -0.01*z[599];
  dz_dqd[2][8] = 0.01*z[600];
  dz_dqd[3][0] = 0;
  dz_dqd[3][1] = 0;
  dz_dqd[3][2] = -0.01*z[595];
  dz_dqd[3][3] = -0.01*z[595];
  dz_dqd[3][4] = -0.01*z[596];
  dz_dqd[3][5] = 0.01*z[597];
  dz_dqd[3][6] = 0;
  dz_dqd[3][7] = 0;
  dz_dqd[3][8] = 0;
  dz_dqd[4][0] = 0;
  dz_dqd[4][1] = 0;
  dz_dqd[4][2] = -0.01*z[601];
  dz_dqd[4][3] = -0.01*z[601];
  dz_dqd[4][4] = -0.01*z[602];
  dz_dqd[4][5] = 0.01*z[603];
  dz_dqd[4][6] = 0;
  dz_dqd[4][7] = 0;
  dz_dqd[4][8] = 0;
  dz_dqd[5][0] = 0;
  dz_dqd[5][1] = 0;
  dz_dqd[5][2] = -0.01*z[604];
  dz_dqd[5][3] = -0.01*z[604];
  dz_dqd[5][4] = -0.01*z[605];
  dz_dqd[5][5] = 0.01*z[606];
  dz_dqd[5][6] = 0;
  dz_dqd[5][7] = 0;
  dz_dqd[5][8] = 0;
  dz_dqd[6][0] = 0;
  dz_dqd[6][1] = 0;
  dz_dqd[6][2] = -0.01*z[598];
  dz_dqd[6][3] = 0;
  dz_dqd[6][4] = 0;
  dz_dqd[6][5] = 0;
  dz_dqd[6][6] = -0.01*z[598];
  dz_dqd[6][7] = -0.01*z[599];
  dz_dqd[6][8] = 0.01*z[600];
  dz_dqd[7][0] = 0;
  dz_dqd[7][1] = 0;
  dz_dqd[7][2] = -0.01*z[607];
  dz_dqd[7][3] = 0;
  dz_dqd[7][4] = 0;
  dz_dqd[7][5] = 0;
  dz_dqd[7][6] = -0.01*z[607];
  dz_dqd[7][7] = -0.01*z[608];
  dz_dqd[7][8] = 0.01*z[609];
  dz_dqd[8][0] = 0;
  dz_dqd[8][1] = 0;
  dz_dqd[8][2] = -0.01*z[610];
  dz_dqd[8][3] = 0;
  dz_dqd[8][4] = 0;
  dz_dqd[8][5] = 0;
  dz_dqd[8][6] = -0.01*z[610];
  dz_dqd[8][7] = -0.01*z[611];
  dz_dqd[8][8] = 0.01*z[612];
  dz_dqdd[0][0] = -0.01*z[308];
  dz_dqdd[0][1] = -0.01*z[309];
  dz_dqdd[0][2] = -0.01*z[312];
  dz_dqdd[0][3] = -0.01*z[317];
  dz_dqdd[0][4] = -0.01*z[318];
  dz_dqdd[0][5] = -0.01*z[316];
  dz_dqdd[0][6] = -0.01*z[313];
  dz_dqdd[0][7] = -0.01*z[315];
  dz_dqdd[0][8] = -0.01*z[307];
  dz_dqdd[1][0] = -0.01*z[309];
  dz_dqdd[1][1] = -0.01*z[321];
  dz_dqdd[1][2] = -0.01*z[322];
  dz_dqdd[1][3] = -0.01*z[326];
  dz_dqdd[1][4] = -0.01*z[327];
  dz_dqdd[1][5] = -0.01*z[325];
  dz_dqdd[1][6] = -0.01*z[323];
  dz_dqdd[1][7] = -0.01*z[324];
  dz_dqdd[1][8] = -0.01*z[320];
  dz_dqdd[2][0] = -0.01*z[312];
  dz_dqdd[2][1] = -0.01*z[322];
  dz_dqdd[2][2] = -0.01*z[330];
  dz_dqdd[2][3] = -0.01*z[336];
  dz_dqdd[2][4] = -0.01*z[337];
  dz_dqdd[2][5] = -0.01*z[338];
  dz_dqdd[2][6] = -0.01*z[332];
  dz_dqdd[2][7] = -0.01*z[334];
  dz_dqdd[2][8] = -0.01*z[335];
  dz_dqdd[3][0] = -0.01*z[317];
  dz_dqdd[3][1] = -0.01*z[326];
  dz_dqdd[3][2] = -0.01*z[336];
  dz_dqdd[3][3] = -0.01*z[336];
  dz_dqdd[3][4] = -0.01*z[337];
  dz_dqdd[3][5] = -0.01*z[338];
  dz_dqdd[3][6] = 0;
  dz_dqdd[3][7] = 0;
  dz_dqdd[3][8] = 0;
  dz_dqdd[4][0] = -0.01*z[318];
  dz_dqdd[4][1] = -0.01*z[327];
  dz_dqdd[4][2] = -0.01*z[337];
  dz_dqdd[4][3] = -0.01*z[337];
  dz_dqdd[4][4] = -0.01*z[342];
  dz_dqdd[4][5] = -0.01*z[343];
  dz_dqdd[4][6] = 0;
  dz_dqdd[4][7] = 0;
  dz_dqdd[4][8] = 0;
  dz_dqdd[5][0] = -0.01*z[316];
  dz_dqdd[5][1] = -0.01*z[325];
  dz_dqdd[5][2] = -0.01*z[338];
  dz_dqdd[5][3] = -0.01*z[338];
  dz_dqdd[5][4] = -0.01*z[343];
  dz_dqdd[5][5] = -0.01*z[345];
  dz_dqdd[5][6] = 0;
  dz_dqdd[5][7] = 0;
  dz_dqdd[5][8] = 0;
  dz_dqdd[6][0] = -0.01*z[313];
  dz_dqdd[6][1] = -0.01*z[323];
  dz_dqdd[6][2] = -0.01*z[332];
  dz_dqdd[6][3] = 0;
  dz_dqdd[6][4] = 0;
  dz_dqdd[6][5] = 0;
  dz_dqdd[6][6] = -0.01*z[332];
  dz_dqdd[6][7] = -0.01*z[334];
  dz_dqdd[6][8] = -0.01*z[335];
  dz_dqdd[7][0] = -0.01*z[315];
  dz_dqdd[7][1] = -0.01*z[324];
  dz_dqdd[7][2] = -0.01*z[334];
  dz_dqdd[7][3] = 0;
  dz_dqdd[7][4] = 0;
  dz_dqdd[7][5] = 0;
  dz_dqdd[7][6] = -0.01*z[334];
  dz_dqdd[7][7] = -0.01*z[348];
  dz_dqdd[7][8] = -0.01*z[349];
  dz_dqdd[8][0] = -0.01*z[307];
  dz_dqdd[8][1] = -0.01*z[320];
  dz_dqdd[8][2] = -0.01*z[335];
  dz_dqdd[8][3] = 0;
  dz_dqdd[8][4] = 0;
  dz_dqdd[8][5] = 0;
  dz_dqdd[8][6] = -0.01*z[335];
  dz_dqdd[8][7] = -0.01*z[349];
  dz_dqdd[8][8] = -0.01*z[345];
  dz_dmom[0][0] = 0;
  dz_dmom[0][1] = 0;
  dz_dmom[0][2] = 0;
  dz_dmom[0][3] = 0;
  dz_dmom[0][4] = 0;
  dz_dmom[0][5] = 0;
  dz_dmom[1][0] = 0;
  dz_dmom[1][1] = 0;
  dz_dmom[1][2] = 0;
  dz_dmom[1][3] = 0;
  dz_dmom[1][4] = 0;
  dz_dmom[1][5] = 0;
  dz_dmom[2][0] = 0;
  dz_dmom[2][1] = 0;
  dz_dmom[2][2] = 0;
  dz_dmom[2][3] = 0;
  dz_dmom[2][4] = 0;
  dz_dmom[2][5] = 0;
  dz_dmom[3][0] = 0.01;
  dz_dmom[3][1] = 0;
  dz_dmom[3][2] = 0;
  dz_dmom[3][3] = 0;
  dz_dmom[3][4] = 0;
  dz_dmom[3][5] = 0;
  dz_dmom[4][0] = 0;
  dz_dmom[4][1] = 0.01;
  dz_dmom[4][2] = 0;
  dz_dmom[4][3] = 0;
  dz_dmom[4][4] = 0;
  dz_dmom[4][5] = 0;
  dz_dmom[5][0] = 0;
  dz_dmom[5][1] = 0;
  dz_dmom[5][2] = 0.01;
  dz_dmom[5][3] = 0;
  dz_dmom[5][4] = 0;
  dz_dmom[5][5] = 0;
  dz_dmom[6][0] = 0;
  dz_dmom[6][1] = 0;
  dz_dmom[6][2] = 0;
  dz_dmom[6][3] = 0.01;
  dz_dmom[6][4] = 0;
  dz_dmom[6][5] = 0;
  dz_dmom[7][0] = 0;
  dz_dmom[7][1] = 0;
  dz_dmom[7][2] = 0;
  dz_dmom[7][3] = 0;
  dz_dmom[7][4] = 0.01;
  dz_dmom[7][5] = 0;
  dz_dmom[8][0] = 0;
  dz_dmom[8][1] = 0;
  dz_dmom[8][2] = 0;
  dz_dmom[8][3] = 0;
  dz_dmom[8][4] = 0;
  dz_dmom[8][5] = 0.01;
  dz_dgrf[0][0] = 0.01*z[630];
  dz_dgrf[0][1] = 0.01*z[616];
  dz_dgrf[0][2] = 0;
  dz_dgrf[0][3] = 0.01*z[631];
  dz_dgrf[0][4] = 0.01*z[617];
  dz_dgrf[0][5] = 0;
  dz_dgrf[1][0] = 0.01*z[616];
  dz_dgrf[1][1] = 0.01*z[632];
  dz_dgrf[1][2] = 0;
  dz_dgrf[1][3] = 0.01*z[617];
  dz_dgrf[1][4] = 0.01*z[633];
  dz_dgrf[1][5] = 0;
  dz_dgrf[2][0] = -0.01*z[618];
  dz_dgrf[2][1] = -0.01*z[619];
  dz_dgrf[2][2] = 0.01*par__bodyweight;
  dz_dgrf[2][3] = -0.01*z[620];
  dz_dgrf[2][4] = -0.01*z[621];
  dz_dgrf[2][5] = 0.01*par__bodyweight;
  dz_dgrf[3][0] = -0.01*z[618];
  dz_dgrf[3][1] = -0.01*z[619];
  dz_dgrf[3][2] = 0.01*par__bodyweight;
  dz_dgrf[3][3] = 0;
  dz_dgrf[3][4] = 0;
  dz_dgrf[3][5] = 0;
  dz_dgrf[4][0] = -0.01*z[622];
  dz_dgrf[4][1] = -0.01*z[623];
  dz_dgrf[4][2] = 0.01*par__bodyweight;
  dz_dgrf[4][3] = 0;
  dz_dgrf[4][4] = 0;
  dz_dgrf[4][5] = 0;
  dz_dgrf[5][0] = 0.01*z[624];
  dz_dgrf[5][1] = -0.01*z[625];
  dz_dgrf[5][2] = 0.01*par__bodyweight;
  dz_dgrf[5][3] = 0;
  dz_dgrf[5][4] = 0;
  dz_dgrf[5][5] = 0;
  dz_dgrf[6][0] = 0;
  dz_dgrf[6][1] = 0;
  dz_dgrf[6][2] = 0;
  dz_dgrf[6][3] = -0.01*z[620];
  dz_dgrf[6][4] = -0.01*z[621];
  dz_dgrf[6][5] = 0.01*par__bodyweight;
  dz_dgrf[7][0] = 0;
  dz_dgrf[7][1] = 0;
  dz_dgrf[7][2] = 0;
  dz_dgrf[7][3] = -0.01*z[626];
  dz_dgrf[7][4] = -0.01*z[627];
  dz_dgrf[7][5] = 0.01*par__bodyweight;
  dz_dgrf[8][0] = 0;
  dz_dgrf[8][1] = 0;
  dz_dgrf[8][2] = 0;
  dz_dgrf[8][3] = 0.01*z[628];
  dz_dgrf[8][4] = -0.01*z[629];
  dz_dgrf[8][5] = 0.01*par__bodyweight;

  Encode[0] = 0.0;
  Encode[1] = 0.0;
  Encode[2] = 0.0;
  Encode[3] = 0.0;

  Encode[4] = 0.0;
  Encode[5] = 0.0;
  Encode[6] = 0.0;
  Encode[7] = 0.0;
  Encode[8] = 0.0;
  Encode[9] = 0.0;

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


