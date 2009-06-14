//**********************************************************************************/
// file name	mmcm.resamp.h
// purpose		MMCM wrap function for R shared library (header file)
// license    GPL-3
// Copyright (c) 2009, Kengo NAGASHIMA and Yasunori SATO.
//**********************************************************************************/
typedef unsigned uint32_t;

struct mmcmdat {
	double param;
	uint32_t clsrnd;
};

void init_gen_rand(uint32_t);
uint32_t gen_rand32(void);
 
void mmcm_rwrap(double *, int *, double *, int *, int *, int *, int *, int *, int *);

int stat_resamp(struct mmcmdat *, double *, long, long, int, int, int, int *);
int mean_get(struct mmcmdat *, int, int, int, long *, double, double *);
int rmean_get(struct mmcmdat *, int, int, int, long *, double, double *);
int gsize_chk(struct mmcmdat *, int, int, int *, long *, double *);
int stat_calc(double *, double *, double *, int, int, double *);
int stat_denom(double *, int, int, double *);
int clsrnd_cmp(const void *_p0, const void *_p1);
