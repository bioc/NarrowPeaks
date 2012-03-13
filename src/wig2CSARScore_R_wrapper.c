#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <stdio.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

int wig2CSARScore(char** file_Name,int* nbChr, int* chrL,char** filenames,char** chr);

void wig2CSARScore_R_wrapper(char** file_Name,int* nbChr, int* chrL,char** filenames,char** chr, int *digits) {
	*digits = wig2CSARScore(file_Name,nbChr,chrL,filenames,chr);
}
