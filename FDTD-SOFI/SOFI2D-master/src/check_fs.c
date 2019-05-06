/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*-------------------------------------------------------------
 *  Check file system of SOFI2D (Part I)
 *
 *  Check
 *  - command line parameters
 *  - existence of parameter file(s) and if they are readable
 *  - existence of directories and if they are read-/writable
 *  
 *  Part II:
 *    depending on the content of the input file(s) further
 *    checks are done in read_par() resp. read_par_auto()
 *
 *  -----------------------------------------------------------*/


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "fd.h"


void check_fs(FILE *fp, int argc, char *fileinp)
{

	FILE *fpinp;
	char modestr[40], infostr[40], fdflt[40]="", errmsg[80];
	int runmd=0, finplen, fserr=0;


	/********************************************/
	/* Check number of commandline parameters   */
	/********************************************/
	
	if (argc < 2)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: Missing commandline parameter for SOFI2D!\n");
		fprintf(fp, "         Parameter: name of input file");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (argc > 2)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: Too many commandline parameters for SOFI2D (%d)!\n", argc-1);
		fprintf(fp, "         SOFI2D only accepts one parameter (name of input file)!");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	
	
	/********************************************/
	/* Check parameter file(s)                  */
	/********************************************/
	
	/* >>> input file for standard usage of SOFI2D */
	if (access(fileinp,0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The stated input file does not exist!\n");
		fprintf(fp, "         File name: <%s>", fileinp);
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access(fileinp,4) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The stated input file does not have read access!\n");
		fprintf(fp, "         File name: <%s>", fileinp);
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	
	
	/* >>> input file for auto mode of SOFI2D */
	/* Open parameter-file to check if auto mode or not*/
	fpinp = fopen(fileinp,"r");
	fscanf(fpinp, "%s %s = %i", infostr, modestr, &runmd);
	fclose(fpinp);
	if (runmd == 1)
	{
		/* construct file name for default input file, for example sofi2D_auto.inp -> sofi2D_auto.default */
		finplen = strlen(fileinp);
		strncpy(fdflt, fileinp, finplen-3);
		fdflt[finplen] = '\0';
		strcat(fdflt, "default");
		
		/* check access */
		if (access(fdflt,0) != 0)
		{
			fprintf(fp, "\n==================================================================\n");
			fprintf(fp, "  ERROR: The following parameter file does not exist!\n");
			fprintf(fp, "         File name: <%s>\n", fdflt);
			fprintf(fp, "         Auto mode of SOFI2D requires this file!");
			fprintf(fp, "\n==================================================================\n");
			fserr = 1;
		}
		else if (access(fdflt,4) != 0)
		{
			fprintf(fp, "\n==================================================================\n");
			fprintf(fp, "  ERROR: The following parameter file does not have read access!\n");
			fprintf(fp, "         File name: <%s>\n", fdflt);
			fprintf(fp, "         Auto mode of SOFI2D requires this file!");
			fprintf(fp, "\n==================================================================\n");
			fserr = 1;
		}
	}	
	
	
	/********************************************/
	/* Check subdirectories                     */
	/********************************************/
	if (access("su",0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/su/> for seismogram storage does\n");
		fprintf(fp, "         not exist!");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access("su",2) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/su/> does not have write access!\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
		
	if (access("snap",0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/snap/> for snapshot storage does\n");
		fprintf(fp, "         not exist!");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access("snap",2) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/snap/> does not have write access!\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	
	if (access("model",0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/model/> for model storage does\n");
		fprintf(fp, "         not exist!");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access("model",6) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/model/> does not have read and/or\n");
		fprintf(fp, "         write access!");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	
	if (access("log",0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The log directory <../par/log/> does not exist!\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access("log",2) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The subdirectory <../par/log/> does not have write access!\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}

	/********************************************/
	/* Check directory ../fsofi2D/par/       */
	/********************************************/
	if (access("../par",0) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The directory <../par/> does not exist\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}
	else if (access("../par",2) != 0)
	{
		fprintf(fp, "\n==================================================================\n");
		fprintf(fp, "  ERROR: The directory <../par/> does not have write access!\n");
		fprintf(fp, "\n==================================================================\n");
		fserr = 1;
	}



	/********************************************/
	/* ERROR                                    */
	/********************************************/
	if (fserr)
	{
		fprintf(fp, "\n");
		sprintf(errmsg, "\n  in: <check_fs.c> \n");
		declare_error(errmsg);
	}	

}
