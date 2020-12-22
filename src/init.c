#include <R_ext/Rdynload.h>
#include "RBF.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_NativePrimitiveArgType kernel_huber_pos_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_tukey_pos_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_huber_lin_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType kernel_tukey_lin_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType huber_pos_t[] = {
    INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType tukey_pos_t[] = {
    INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP
};


static const R_CMethodDef CEntries[]  = {
    CDEF(kernel_huber_pos),
    CDEF(kernel_tukey_pos),
    CDEF(kernel_huber_lin),
    CDEF(kernel_tukey_lin),
    CDEF(huber_pos),
    CDEF(tukey_pos),
    {NULL, NULL, 0}
};


void R_init_RBF(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
