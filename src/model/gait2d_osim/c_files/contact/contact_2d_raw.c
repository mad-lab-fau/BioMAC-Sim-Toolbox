/* The name of this program is contact_raw.c */
/* Created by Autolev 4.1 on Thu Nov 28 14:40:12 2024 */

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

double   contact__xp,contact__yp,fk1,fk2,fk3,fk4,fk5,fk6,x1,x2,x3,x4,xdot3,xdot4;
double   Fx,Fy,xc,xcloc,yc,ycloc,f1,f2,f3,f4;
double   Pi,DEGtoRAD,RADtoDEG,z[7],f[4],df_dfk[4][6],df_dfkdot[4][6],df_dx[4][
  4],df_dxdot[4][4],Encode[5];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[1];
int      iloop;

/* Open input and output files */
for(iloop=0;  iloop<=0;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "contact_raw.in");
  else sprintf(fileName, "contact_raw.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&contact__xp,&contact__yp,&fk1,&fk2,&fk3,&fk4,&fk5,&fk6,&x1,&x2,&
  x3,&x4,&xdot3,&xdot4,NULL);

/* Write heading(s) to output file(s) */

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate output quantities */
eqns1();
output( Fptr );

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file contact_raw.in" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (void)
{

/* Evaluate constants */
  Fx = x1;
  Fy = x2;
  xc = x3;
  yc = x4;
  xcloc = -fk3*(fk1-xc) - fk5*(fk2-yc);
  ycloc = -fk4*(fk1-xc) - fk6*(fk2-yc);
  f1 = xcloc - contact__xp;
  f2 = ycloc - contact__yp;
  f3 = Fy - 37.5*(-1.333333333333333+xdot4)*(yc-pow((1.0E-06+pow(yc,2)),0.5));
  f4 = Fx + 100*Fy*xdot3/pow((1+10000*pow(xdot3,2)),0.5);
  z[1] = xc - fk1;
  z[2] = yc - fk2;
  z[3] = (-1.333333333333333+xdot4)*(1-yc/pow((1.0E-06+pow(yc,2)),0.5));
  z[4] = xdot3/pow((1+10000*pow(xdot3,2)),0.5);
  z[5] = 37.5*pow((1.0E-06+pow(yc,2)),0.5) - 37.5*yc;
  z[6] = Fy/pow((1+10000*pow(xdot3,2)),1.5);


/* Evaluate output quantities */
}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[] )
{
int      i1;

/* Write output to screen and to output file(s) */
  f[0] = f1;
  f[1] = f2;
  f[2] = f3;
  f[3] = f4;
  df_dfk[0][0] = -fk3;
  df_dfk[0][1] = -fk5;
  df_dfk[0][2] = z[1];
  df_dfk[0][3] = 0;
  df_dfk[0][4] = z[2];
  df_dfk[0][5] = 0;
  df_dfk[1][0] = -fk4;
  df_dfk[1][1] = -fk6;
  df_dfk[1][2] = 0;
  df_dfk[1][3] = z[1];
  df_dfk[1][4] = 0;
  df_dfk[1][5] = z[2];
  df_dfk[2][0] = 0;
  df_dfk[2][1] = 0;
  df_dfk[2][2] = 0;
  df_dfk[2][3] = 0;
  df_dfk[2][4] = 0;
  df_dfk[2][5] = 0;
  df_dfk[3][0] = 0;
  df_dfk[3][1] = 0;
  df_dfk[3][2] = 0;
  df_dfk[3][3] = 0;
  df_dfk[3][4] = 0;
  df_dfk[3][5] = 0;
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
  df_dfkdot[2][0] = 0;
  df_dfkdot[2][1] = 0;
  df_dfkdot[2][2] = 0;
  df_dfkdot[2][3] = 0;
  df_dfkdot[2][4] = 0;
  df_dfkdot[2][5] = 0;
  df_dfkdot[3][0] = 0;
  df_dfkdot[3][1] = 0;
  df_dfkdot[3][2] = 0;
  df_dfkdot[3][3] = 0;
  df_dfkdot[3][4] = 0;
  df_dfkdot[3][5] = 0;
  df_dx[0][0] = 0;
  df_dx[0][1] = 0;
  df_dx[0][2] = fk3;
  df_dx[0][3] = fk5;
  df_dx[1][0] = 0;
  df_dx[1][1] = 0;
  df_dx[1][2] = fk4;
  df_dx[1][3] = fk6;
  df_dx[2][0] = 0;
  df_dx[2][1] = 1;
  df_dx[2][2] = 0;
  df_dx[2][3] = -37.5*z[3];
  df_dx[3][0] = 1;
  df_dx[3][1] = 100*z[4];
  df_dx[3][2] = 0;
  df_dx[3][3] = 0;
  df_dxdot[0][0] = 0;
  df_dxdot[0][1] = 0;
  df_dxdot[0][2] = 0;
  df_dxdot[0][3] = 0;
  df_dxdot[1][0] = 0;
  df_dxdot[1][1] = 0;
  df_dxdot[1][2] = 0;
  df_dxdot[1][3] = 0;
  df_dxdot[2][0] = 0;
  df_dxdot[2][1] = 0;
  df_dxdot[2][2] = 0;
  df_dxdot[2][3] = z[5];
  df_dxdot[3][0] = 0;
  df_dxdot[3][1] = 0;
  df_dxdot[3][2] = 100*z[6];
  df_dxdot[3][3] = 0;

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


