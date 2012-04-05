#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <stdio.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Memory.h>

int wig2CSARScore(const char* file_Name,int* nbChr, int* chrL,char** filenames, char** chr);

SEXP wig2CSARScore_R_wrapper(SEXP file_Name, SEXP nbChr, SEXP chrL)
{
    int c_nbChr = INTEGER(nbChr)[0];
    char /* R_alloc'd memory is automatically freed when we return to R */
        **c_filenames = (char **) R_alloc(c_nbChr, sizeof(char *)),
        **c_chr = (char **) R_alloc(c_nbChr, sizeof(char *));
    SEXP result, r_filenames, r_chr, r_digits;

    int c_digits = wig2CSARScore(CHAR(STRING_ELT(file_Name, 0)),
                                 INTEGER(nbChr), INTEGER(chrL),
                                 c_filenames, c_chr);

    /* 'SEXP's need to be PROTECTed from garbage collection */
    result = PROTECT(allocVector(VECSXP, 4)); /* explicit PROTECTion */
    r_chr = allocVector(STRSXP, c_nbChr);
    SET_VECTOR_ELT(result, 0, r_chr); /* 1st element -- chr; implicit PROTECTion */
    SET_VECTOR_ELT(result, 1, chrL);  /* 2nd element -- chrL */
    r_filenames = allocVector(STRSXP, c_nbChr),
    SET_VECTOR_ELT(result, 2, r_filenames); /* 3rd element -- filenames */
    r_digits = allocVector(REALSXP, 1);
    SET_VECTOR_ELT(result, 3, r_digits); /* 4th element -- digits */

    /* copy from C -> R */
    for (int i = 0; i < c_nbChr; ++i) {
        SET_STRING_ELT(r_chr, i, mkChar(c_chr[i]));
        SET_STRING_ELT(r_filenames, i, mkChar(c_filenames[i]));
    }
    REAL(r_digits)[0] = c_digits;

    /* clean up */
    UNPROTECT(1);
    return result;
}

