/* The name of this program is contact_raw.c */
/* Created by Autolev 4.1 on Fri Aug 23 14:38:37 2024 */

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

double   contact__xp,contact__yp,contact__zp,fk1,fk10,fk11,fk12,fk2,fk3,fk4,fk5,
  fk6,fk7,fk8,fk9,x1,x2,x3,x4,x5,x6,xdot4,xdot5,xdot6;
double   Fx,Fy,Fz,xc,xcloc,yc,ycloc,zc,zcloc,f1,f2,f3,f4,f5,f6;
double   Pi,DEGtoRAD,RADtoDEG,z[10],f[6],df_dfk[6][12],df_dfkdot[6][12],
  df_dx[6][6],df_dxdot[6][6],Encode[5];

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
readf(Fptr[0],&contact__xp,&contact__yp,&contact__zp,&fk1,&fk10,&fk11,&fk12,&fk2,&
  fk3,&fk4,&fk5,&fk6,&fk7,&fk8,&fk9,&x1,&x2,&x3,&x4,&x5,&x6,&xdot4,&xdot5,&xdot6,
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
puts( "\n Input is in the file contact_raw.in" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (void)
{

/* Evaluate constants */
  Fx = x1;
  Fy = x2;
  Fz = x3;
  xc = x4;
  yc = x5;
  zc = x6;
  xcloc = -fk10*(fk3-zc) - fk4*(fk1-xc) - fk7*(fk2-yc);
  ycloc = -fk11*(fk3-zc) - fk5*(fk1-xc) - fk8*(fk2-yc);
  zcloc = -fk12*(fk3-zc) - fk6*(fk1-xc) - fk9*(fk2-yc);
  f1 = xcloc - contact__xp;
  f2 = ycloc - contact__yp;
  f3 = zcloc - contact__zp;
  f4 = Fy - 37.5*(-1.333333333333333+xdot5)*(yc-pow((1.0E-06+pow(yc,2)),0.5));
  f5 = Fx + 100*Fy*xdot4/pow((1+10000*pow(xdot4,2)),0.5);
  f6 = Fz + 100*Fy*xdot6/pow((1+10000*pow(xdot6,2)),0.5);
  z[1] = xc - fk1;
  z[2] = yc - fk2;
  z[3] = zc - fk3;
  z[4] = (-1.333333333333333+xdot5)*(1-yc/pow((1.0E-06+pow(yc,2)),0.5));
  z[5] = xdot4/pow((1+10000*pow(xdot4,2)),0.5);
  z[6] = xdot6/pow((1+10000*pow(xdot6,2)),0.5);
  z[7] = 37.5*pow((1.0E-06+pow(yc,2)),0.5) - 37.5*yc;
  z[8] = Fy/pow((1+10000*pow(xdot4,2)),1.5);
  z[9] = Fy/pow((1+10000*pow(xdot6,2)),1.5);


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
  f[4] = f5;
  f[5] = f6;
  df_dfk[0][0] = -fk4;
  df_dfk[0][1] = -fk7;
  df_dfk[0][2] = -fk10;
  df_dfk[0][3] = z[1];
  df_dfk[0][4] = 0;
  df_dfk[0][5] = 0;
  df_dfk[0][6] = z[2];
  df_dfk[0][7] = 0;
  df_dfk[0][8] = 0;
  df_dfk[0][9] = z[3];
  df_dfk[0][10] = 0;
  df_dfk[0][11] = 0;
  df_dfk[1][0] = -fk5;
  df_dfk[1][1] = -fk8;
  df_dfk[1][2] = -fk11;
  df_dfk[1][3] = 0;
  df_dfk[1][4] = z[1];
  df_dfk[1][5] = 0;
  df_dfk[1][6] = 0;
  df_dfk[1][7] = z[2];
  df_dfk[1][8] = 0;
  df_dfk[1][9] = 0;
  df_dfk[1][10] = z[3];
  df_dfk[1][11] = 0;
  df_dfk[2][0] = -fk6;
  df_dfk[2][1] = -fk9;
  df_dfk[2][2] = -fk12;
  df_dfk[2][3] = 0;
  df_dfk[2][4] = 0;
  df_dfk[2][5] = z[1];
  df_dfk[2][6] = 0;
  df_dfk[2][7] = 0;
  df_dfk[2][8] = z[2];
  df_dfk[2][9] = 0;
  df_dfk[2][10] = 0;
  df_dfk[2][11] = z[3];
  df_dfk[3][0] = 0;
  df_dfk[3][1] = 0;
  df_dfk[3][2] = 0;
  df_dfk[3][3] = 0;
  df_dfk[3][4] = 0;
  df_dfk[3][5] = 0;
  df_dfk[3][6] = 0;
  df_dfk[3][7] = 0;
  df_dfk[3][8] = 0;
  df_dfk[3][9] = 0;
  df_dfk[3][10] = 0;
  df_dfk[3][11] = 0;
  df_dfk[4][0] = 0;
  df_dfk[4][1] = 0;
  df_dfk[4][2] = 0;
  df_dfk[4][3] = 0;
  df_dfk[4][4] = 0;
  df_dfk[4][5] = 0;
  df_dfk[4][6] = 0;
  df_dfk[4][7] = 0;
  df_dfk[4][8] = 0;
  df_dfk[4][9] = 0;
  df_dfk[4][10] = 0;
  df_dfk[4][11] = 0;
  df_dfk[5][0] = 0;
  df_dfk[5][1] = 0;
  df_dfk[5][2] = 0;
  df_dfk[5][3] = 0;
  df_dfk[5][4] = 0;
  df_dfk[5][5] = 0;
  df_dfk[5][6] = 0;
  df_dfk[5][7] = 0;
  df_dfk[5][8] = 0;
  df_dfk[5][9] = 0;
  df_dfk[5][10] = 0;
  df_dfk[5][11] = 0;
  df_dfkdot[0][0] = 0;
  df_dfkdot[0][1] = 0;
  df_dfkdot[0][2] = 0;
  df_dfkdot[0][3] = 0;
  df_dfkdot[0][4] = 0;
  df_dfkdot[0][5] = 0;
  df_dfkdot[0][6] = 0;
  df_dfkdot[0][7] = 0;
  df_dfkdot[0][8] = 0;
  df_dfkdot[0][9] = 0;
  df_dfkdot[0][10] = 0;
  df_dfkdot[0][11] = 0;
  df_dfkdot[1][0] = 0;
  df_dfkdot[1][1] = 0;
  df_dfkdot[1][2] = 0;
  df_dfkdot[1][3] = 0;
  df_dfkdot[1][4] = 0;
  df_dfkdot[1][5] = 0;
  df_dfkdot[1][6] = 0;
  df_dfkdot[1][7] = 0;
  df_dfkdot[1][8] = 0;
  df_dfkdot[1][9] = 0;
  df_dfkdot[1][10] = 0;
  df_dfkdot[1][11] = 0;
  df_dfkdot[2][0] = 0;
  df_dfkdot[2][1] = 0;
  df_dfkdot[2][2] = 0;
  df_dfkdot[2][3] = 0;
  df_dfkdot[2][4] = 0;
  df_dfkdot[2][5] = 0;
  df_dfkdot[2][6] = 0;
  df_dfkdot[2][7] = 0;
  df_dfkdot[2][8] = 0;
  df_dfkdot[2][9] = 0;
  df_dfkdot[2][10] = 0;
  df_dfkdot[2][11] = 0;
  df_dfkdot[3][0] = 0;
  df_dfkdot[3][1] = 0;
  df_dfkdot[3][2] = 0;
  df_dfkdot[3][3] = 0;
  df_dfkdot[3][4] = 0;
  df_dfkdot[3][5] = 0;
  df_dfkdot[3][6] = 0;
  df_dfkdot[3][7] = 0;
  df_dfkdot[3][8] = 0;
  df_dfkdot[3][9] = 0;
  df_dfkdot[3][10] = 0;
  df_dfkdot[3][11] = 0;
  df_dfkdot[4][0] = 0;
  df_dfkdot[4][1] = 0;
  df_dfkdot[4][2] = 0;
  df_dfkdot[4][3] = 0;
  df_dfkdot[4][4] = 0;
  df_dfkdot[4][5] = 0;
  df_dfkdot[4][6] = 0;
  df_dfkdot[4][7] = 0;
  df_dfkdot[4][8] = 0;
  df_dfkdot[4][9] = 0;
  df_dfkdot[4][10] = 0;
  df_dfkdot[4][11] = 0;
  df_dfkdot[5][0] = 0;
  df_dfkdot[5][1] = 0;
  df_dfkdot[5][2] = 0;
  df_dfkdot[5][3] = 0;
  df_dfkdot[5][4] = 0;
  df_dfkdot[5][5] = 0;
  df_dfkdot[5][6] = 0;
  df_dfkdot[5][7] = 0;
  df_dfkdot[5][8] = 0;
  df_dfkdot[5][9] = 0;
  df_dfkdot[5][10] = 0;
  df_dfkdot[5][11] = 0;
  df_dx[0][0] = 0;
  df_dx[0][1] = 0;
  df_dx[0][2] = 0;
  df_dx[0][3] = fk4;
  df_dx[0][4] = fk7;
  df_dx[0][5] = fk10;
  df_dx[1][0] = 0;
  df_dx[1][1] = 0;
  df_dx[1][2] = 0;
  df_dx[1][3] = fk5;
  df_dx[1][4] = fk8;
  df_dx[1][5] = fk11;
  df_dx[2][0] = 0;
  df_dx[2][1] = 0;
  df_dx[2][2] = 0;
  df_dx[2][3] = fk6;
  df_dx[2][4] = fk9;
  df_dx[2][5] = fk12;
  df_dx[3][0] = 0;
  df_dx[3][1] = 1;
  df_dx[3][2] = 0;
  df_dx[3][3] = 0;
  df_dx[3][4] = -37.5*z[4];
  df_dx[3][5] = 0;
  df_dx[4][0] = 1;
  df_dx[4][1] = 100*z[5];
  df_dx[4][2] = 0;
  df_dx[4][3] = 0;
  df_dx[4][4] = 0;
  df_dx[4][5] = 0;
  df_dx[5][0] = 0;
  df_dx[5][1] = 100*z[6];
  df_dx[5][2] = 1;
  df_dx[5][3] = 0;
  df_dx[5][4] = 0;
  df_dx[5][5] = 0;
  df_dxdot[0][0] = 0;
  df_dxdot[0][1] = 0;
  df_dxdot[0][2] = 0;
  df_dxdot[0][3] = 0;
  df_dxdot[0][4] = 0;
  df_dxdot[0][5] = 0;
  df_dxdot[1][0] = 0;
  df_dxdot[1][1] = 0;
  df_dxdot[1][2] = 0;
  df_dxdot[1][3] = 0;
  df_dxdot[1][4] = 0;
  df_dxdot[1][5] = 0;
  df_dxdot[2][0] = 0;
  df_dxdot[2][1] = 0;
  df_dxdot[2][2] = 0;
  df_dxdot[2][3] = 0;
  df_dxdot[2][4] = 0;
  df_dxdot[2][5] = 0;
  df_dxdot[3][0] = 0;
  df_dxdot[3][1] = 0;
  df_dxdot[3][2] = 0;
  df_dxdot[3][3] = 0;
  df_dxdot[3][4] = z[7];
  df_dxdot[3][5] = 0;
  df_dxdot[4][0] = 0;
  df_dxdot[4][1] = 0;
  df_dxdot[4][2] = 0;
  df_dxdot[4][3] = 100*z[8];
  df_dxdot[4][4] = 0;
  df_dxdot[4][5] = 0;
  df_dxdot[5][0] = 0;
  df_dxdot[5][1] = 0;
  df_dxdot[5][2] = 0;
  df_dxdot[5][3] = 0;
  df_dxdot[5][4] = 0;
  df_dxdot[5][5] = 100*z[9];

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


