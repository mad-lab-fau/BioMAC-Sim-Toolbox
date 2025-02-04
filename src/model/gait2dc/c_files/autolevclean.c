/** 
 *
 * @file autolevclean.c
 * @brief MEX function to extract clean C source code from C files generated by Autolevs
 *
 * @author Antonie J. (Ton) van den Bogert
 * @date January 25, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mex.h"
#include "gait2dc.h"
#define STRLEN 100 				// source lines are never longer than 100 characters

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	FILE *fid1, *fid2;
	char line[STRLEN];
	int i, copying;
	char *ptr;

    //==================================================================
	printf("Cleaning Autolev output from gait2dc_dyn.al...\n");

	// open the file that came from Autolev
	if ((fid1 = fopen("gait2dc_dyn_al_raw.c","r")) == NULL) {
		printf("Could not open gait2dc_dyn_al_raw.c\n");
		exit(1);
	}
	
// write the clean C file
	if ((fid2 = fopen("gait2dc_dyn_al.c","w")) == NULL) {
		printf("Could not write gait2dc_dyn_al.c\n");
		exit(1);
	}
	fprintf(fid2,"// This file was generated by autolevclean.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait2dc.h\"\n");
	fprintf(fid2,"void gait2dc_dynamics_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n");
	fprintf(fid2,"   double mom[NMOM], double grf[NGRF], double Zero[NDOF], double dz_dq[NDOF][NDOF],\n");
	fprintf(fid2,"   double dz_dqd[NDOF][NDOF], double dz_dqdd[NDOF][NDOF],\n");
	fprintf(fid2,"   double dz_dmom[NDOF][NMOM], double dz_dgrf[NDOF][NGRF],\n");
	fprintf(fid2,"   double fk[NFK], double dfk_dq[NFK][NDOF],\n");
	fprintf(fid2,"   double fkdot[NFK], double dfkdot_dq[NFK][NDOF]) {\n");;		
	
	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}
	
	// generate C code to copy joint moments into scalar variables
	for (i=0; i<NMOM; i++) {
		fprintf(fid2,"\tdouble mom%1d = mom[%1d];\n", i+1,i);			
	}

	// generate C code to copy grf into scalar variables
	fprintf(fid2,"\tdouble RFootFx = grf[0];\n");			
	fprintf(fid2,"\tdouble RFootFy = grf[1];\n");			
	fprintf(fid2,"\tdouble RFootMz = grf[2];\n");			
	fprintf(fid2,"\tdouble LFootFx = grf[3];\n");			
	fprintf(fid2,"\tdouble LFootFy = grf[4];\n");			
	fprintf(fid2,"\tdouble LFootMz = grf[5];\n");			

	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0))				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
	}
	
	// close the input file
	fclose(fid1);
	
	// close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);
	

    //==================================================================
    printf("Cleaning Autolev output from gait2dc_FK.al...\n");

	// open the file that came from Autolev
	if ((fid1 = fopen("gait2dc_FK_al_raw.c","r")) == NULL) {
		printf("Could not open gait2dc_FK_al_raw.c\n");
		exit(1);
	}
	
    // write the clean C file
	if ((fid2 = fopen("gait2dc_FK_al.c","w")) == NULL) {
		printf("Could not write gait2dc_FK_al.c\n");
		exit(1);
	}
	fprintf(fid2,"// This file was generated by autolevclean.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait2dc.h\"\n");
	fprintf(fid2,"void gait2dc_FK_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n");
	fprintf(fid2,"   double fk[NFK], double dfk_dq[NFK][NDOF],\n");
	fprintf(fid2,"   double fkdot[NFK], double dfkdot_dq[NFK][NDOF]) {\n");		
	
	// generate C code to copy qinto scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);
        fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);
	}

	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0))				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
	}

    // close the input file
	fclose(fid1);

	// close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);

    //==================================================================
    printf("Cleaning Autolev output from gait2dc_stick.al...\n");

	// open the file that came from Autolev
	if ((fid1 = fopen("gait2dc_stick_al_raw.c","r")) == NULL) {
		printf("Could not open gait2dc_stick_al_raw.c\n");
		exit(1);
	}
	
    // write the clean C file
	if ((fid2 = fopen("gait2dc_stick_al.c","w")) == NULL) {
		printf("Could not write gait2dc_stick_al.c\n");
		exit(1);
	}
	fprintf(fid2,"// This file was generated by autolevclean.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait2dc.h\"\n");
	fprintf(fid2,"void gait2dc_stick_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n");
	fprintf(fid2,"   double fk[NFK], double dfk_dq[NFK][NDOF],\n");
	fprintf(fid2,"   double fkdot[NFK], double dfkdot_dq[NFK][NDOF],\n");	
	fprintf(fid2,"   double Stick[NSTICK][2]) {\n");		
	
	// generate C code to copy qinto scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);
	}

	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0))				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
	}
	
	// close the input file
	fclose(fid1);
	
	// close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);


    //==================================================================
	printf("Cleaning Autolev output from contact.al...\n");

	// open the file that came from Autolev
	if ((fid1 = fopen("contact_al_raw.c","r")) == NULL) {
		printf("Could not open deform_al_raw.c\n");
		exit(1);
	}
	
	// write the clean C file
	if ((fid2 = fopen("contact_al.c","w")) == NULL) {
		printf("Could not write deform_al.c\n");
		exit(1);
	}
	fprintf(fid2,"// This file was generated by autolevclean2.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait2dc.h\"\n");
	fprintf(fid2,"void contact_al(contactprop* contact,\n");										
	fprintf(fid2,"\tdouble fk[6], double fkdot[6], double xc[4], double xcdot[2],\n");			
	fprintf(fid2,"\tdouble f[4], double df_dfk[4][6], double df_dfkdot[4][6],\n");
	fprintf(fid2,"\tdouble df_dxc[4][4], double df_dxcdot[4][2]) {\n\n");		
	
	// generate C code to copy input arrays into scalar variables
	for (i=0; i<6; i++) {
		fprintf(fid2,"\tdouble fk%1d = fk[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble fkdot%1d = fkdot[%1d];\n", i+1,i);	
	}
	for (i=0; i<4; i++) {
		fprintf(fid2,"\tdouble xc%1d = xc[%1d];\n", i+1,i);		
	}
	for (i=0; i<2; i++) {
		fprintf(fid2,"\tdouble xcdot%1d = xcdot[%1d];\n", i+1,i);		
	}
	
	// declare the temporary variables
	fprintf(fid2,"\tdouble d,d1,d2,sd1,sd2,Ffriction,Fnormal;\n");
	
	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0))				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// convert contact__ into contact->
			ptr = strstr(line, "contact__");
			while (ptr != NULL) {
				strncpy(ptr, "contact->", 9);
				ptr = strstr(line, "contact__");
			}
			fputs(line, fid2);
		}
	}
	
	// close the input file
	fclose(fid1);
	
	// close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);
    
    //==================================================================
	printf("Cleaning Autolev output from gait2dc_acc.al...\n");
		
    // open the file that came from Autolev
	if ((fid1 = fopen("gait2dc_acc_al_raw.c","r")) == NULL) {
		printf("Could not open gait2dc_acc_al_raw.c\n");
		exit(1);
	}
	
    // write the clean C file
	if ((fid2 = fopen("gait2dc_acc_al.c","w")) == NULL) {
		printf("Could not write gait2dc_acc_al.c\n");
		exit(1);
	}
	fprintf(fid2,"// This file was generated by autolevclean.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait2dc.h\"\n");
	fprintf(fid2,"void acc_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n");
	fprintf(fid2,"   double acc[42], double dacc_dq[42][NDOF],\n");
	fprintf(fid2,"   double dacc_dqd[42][NDOF], double dacc_dqdd[42][NDOF]) {\n");

	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}
	
	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// make sure there is enough room for all Zs
		}
		if (strcmp(line, "/* Evaluate constants */\n") == 0) 				// Z[] code starts here
			copying = 1;
		else if (strcmp(line, "/* Evaluate output quantities */\n") == 0)	// Z[] code ends here
			copying = 0;
		else if (strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0)		// kinematics code starts here
			copying = 1;
		else if (strcmp(line, "  Encode[0] = 0.0;\n") == 0)	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
	}
	
    // close the input file
	fclose(fid1);
	
    // close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);
			
}			// END OF MAIN PROGRAM
