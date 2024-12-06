/* The name of this program is contact_al_raw.c */
/* Created by Autolev 4.1 on Thu Oct 24 13:52:04 2024 */

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

double   contact__b,contact__c,contact__kxx,contact__kxy,contact__kyy,contact__LambdaX,
  contact__LambdaY,contact__x,contact__y,fk1,fk2,fk3,fk4,fk5,fk6,fkdot1,fkdot2,
  fkdot3,fkdot4,fkdot5,fkdot6,xc1,xc2,xc3,xc4,xcdot1,xcdot2;
double   d,d1,d2,sd1,sd2,Ffriction,Fnormal;
double   Pi,DEGtoRAD,RADtoDEG,z[34],f[4],df_dfk[4][6],df_dfkdot[4][6],df_dxc[4][
  4],df_dxcdot[4][2],Encode[5];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "contact_al_raw.in");
  else sprintf(fileName, "contact_al_raw.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&contact__b,&contact__c,&contact__kxx,&contact__kxy,&contact__kyy,&
  contact__LambdaX,&contact__LambdaY,&contact__x,&contact__y,&fk1,&fk2,&fk3,&fk4,&
  fk5,&fk6,&fkdot1,&fkdot2,&fkdot3,&fkdot4,&fkdot5,&fkdot6,&xc1,&xc2,&xc3,&xc4,&
  xcdot1,&xcdot2,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file contact_al_raw.in" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (void)
{

/* Evaluate constants */
  d1 = contact__c*xc4*(-1+10*xcdot1) - xc3;
  d2 = contact__c*xc4*(1+10*xcdot1) - xc3;
  sd1 = pow((contact__LambdaX+pow(d1,2)),0.5);
  sd2 = pow((contact__LambdaX+pow(d2,2)),0.5);
  Ffriction = sd1 - 2*xc3 - sd2;
  d = pow((contact__LambdaY+pow(xc4,2)+40000*pow(xc2,2)),0.5);
  Fnormal = xc4 + 200*xc2 - d;
  z[1] = fk3*xc3 + fk5*xc4 + contact__kxx*(contact__x+fk3*(fk1-xc1)+fk5*(fk2-
  xc2)) + contact__kxy*(contact__y+fk4*(fk1-xc1)+fk6*(fk2-xc2)) + contact__b*(
  fk3*(fkdot1-xcdot1)+fk5*(fkdot2-xcdot2)+fkdot3*(fk1-xc1)+fkdot5*(fk2-xc2));
  z[2] = fk4*xc3 + fk6*xc4 + contact__kxy*(contact__x+fk3*(fk1-xc1)+fk5*(fk2-
  xc2)) + contact__kyy*(contact__y+fk4*(fk1-xc1)+fk6*(fk2-xc2)) + contact__b*(
  fk4*(fkdot1-xcdot1)+fk6*(fkdot2-xcdot2)+fkdot4*(fk1-xc1)+fkdot6*(fk2-xc2));
  z[3] = contact__b*fkdot3 + contact__kxx*fk3 + contact__kxy*fk4;
  z[4] = contact__b*fkdot5 + contact__kxx*fk5 + contact__kxy*fk6;
  z[5] = xc3 + contact__b*(fkdot1-xcdot1) + contact__kxx*(fk1-xc1);
  z[6] = xc4 + contact__b*(fkdot2-xcdot2) + contact__kxx*(fk2-xc2);
  z[7] = contact__b*fkdot4 + contact__kxy*fk3 + contact__kyy*fk4;
  z[8] = contact__b*fkdot6 + contact__kxy*fk5 + contact__kyy*fk6;
  z[9] = xc3 + contact__b*(fkdot1-xcdot1) + contact__kyy*(fk1-xc1);
  z[10] = xc4 + contact__b*(fkdot2-xcdot2) + contact__kyy*(fk2-xc2);
  z[11] = contact__kxy*(fk1-xc1);
  z[12] = contact__kxy*(fk2-xc2);
  z[13] = contact__b*fk3;
  z[14] = contact__b*fk5;
  z[15] = contact__b*(fk1-xc1);
  z[16] = contact__b*(fk2-xc2);
  z[17] = contact__b*fk4;
  z[18] = contact__b*fk6;
  z[19] = -2 + d2/pow((contact__LambdaX+pow(d2,2)),0.5) - d1/pow((contact__LambdaX+
  pow(d1,2)),0.5);
  z[20] = contact__c*d1*(-1+10*xcdot1)/pow((contact__LambdaX+pow(d1,2)),0.5);
  z[21] = contact__c*d2*(1+10*xcdot1)/pow((contact__LambdaX+pow(d2,2)),0.5);
  z[22] = xc2/pow((contact__LambdaY+pow(xc4,2)+40000*pow(xc2,2)),0.5);
  z[23] = xc4/pow((contact__LambdaY+pow(xc4,2)+40000*pow(xc2,2)),0.5);
  z[24] = -contact__b*fkdot3 - contact__kxx*fk3 - contact__kxy*fk4;
  z[25] = -contact__b*fkdot5 - contact__kxx*fk5 - contact__kxy*fk6;
  z[26] = -contact__b*fkdot4 - contact__kxy*fk3 - contact__kyy*fk4;
  z[27] = -contact__b*fkdot6 - contact__kxy*fk5 - contact__kyy*fk6;
  z[28] = z[20] - z[21];
  z[29] = 200 - 40000*z[22];
  z[30] = 1 - z[23];
  z[31] = contact__c*d1*xc4/pow((contact__LambdaX+pow(d1,2)),0.5);
  z[32] = contact__c*d2*xc4/pow((contact__LambdaX+pow(d2,2)),0.5);
  z[33] = 10*z[31] - 10*z[32];


/* Evaluate output quantities */
}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[] )
{
int      i1;

/* Write output to screen and to output file(s) */
  f[0] = Ffriction;
  f[1] = Fnormal;
  f[2] = z[1];
  f[3] = z[2];
  df_dfk[0][0] = 0;
  df_dfk[0][1] = 0;
  df_dfk[0][2] = 0;
  df_dfk[0][3] = 0;
  df_dfk[0][4] = 0;
  df_dfk[0][5] = 0;
  df_dfk[1][0] = 0;
  df_dfk[1][1] = 0;
  df_dfk[1][2] = 0;
  df_dfk[1][3] = 0;
  df_dfk[1][4] = 0;
  df_dfk[1][5] = 0;
  df_dfk[2][0] = z[3];
  df_dfk[2][1] = z[4];
  df_dfk[2][2] = z[5];
  df_dfk[2][3] = z[11];
  df_dfk[2][4] = z[6];
  df_dfk[2][5] = z[12];
  df_dfk[3][0] = z[7];
  df_dfk[3][1] = z[8];
  df_dfk[3][2] = z[11];
  df_dfk[3][3] = z[9];
  df_dfk[3][4] = z[12];
  df_dfk[3][5] = z[10];
  df_dfkdot[0][0] = 0;
  df_dfkdot[0][1] = 0;
  df_dfkdot[0][2] = 0;
  df_dfkdot[0][3] = 0;
  df_dfkdot[0][4] = 0;
  df_dfkdot[0][5] = 0;
  df_dfkdot[1][0] = 0;
  df_dfkdot[1][1] = 0;
  df_dfkdot[1][2] = 0;
  df_dfkdot[1][3] = 0;
  df_dfkdot[1][4] = 0;
  df_dfkdot[1][5] = 0;
  df_dfkdot[2][0] = z[13];
  df_dfkdot[2][1] = z[14];
  df_dfkdot[2][2] = z[15];
  df_dfkdot[2][3] = 0;
  df_dfkdot[2][4] = z[16];
  df_dfkdot[2][5] = 0;
  df_dfkdot[3][0] = z[17];
  df_dfkdot[3][1] = z[18];
  df_dfkdot[3][2] = 0;
  df_dfkdot[3][3] = z[15];
  df_dfkdot[3][4] = 0;
  df_dfkdot[3][5] = z[16];
  df_dxc[0][0] = 0;
  df_dxc[0][1] = 0;
  df_dxc[0][2] = z[19];
  df_dxc[0][3] = z[28];
  df_dxc[1][0] = 0;
  df_dxc[1][1] = z[29];
  df_dxc[1][2] = 0;
  df_dxc[1][3] = z[30];
  df_dxc[2][0] = z[24];
  df_dxc[2][1] = z[25];
  df_dxc[2][2] = fk3;
  df_dxc[2][3] = fk5;
  df_dxc[3][0] = z[26];
  df_dxc[3][1] = z[27];
  df_dxc[3][2] = fk4;
  df_dxc[3][3] = fk6;
  df_dxcdot[0][0] = z[33];
  df_dxcdot[0][1] = 0;
  df_dxcdot[1][0] = 0;
  df_dxcdot[1][1] = 0;
  df_dxcdot[2][0] = -z[13];
  df_dxcdot[2][1] = -z[14];
  df_dxcdot[3][0] = -z[17];
  df_dxcdot[3][1] = -z[18];

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


