/*
    Copyright (C) 2009, Kengo NAGASHIMA and Yasunori SATO.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

//**********************************************************************************/
// file name    mmcm.resamp.c
// purpose      MMCM wrap function for R shared library
// license      GPL-3
// Copyright (c) 2009, Kengo NAGASHIMA and Yasunori SATO.
//**********************************************************************************/
#define MYMAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MYMIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./mmcm.resamp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

//**********************************************************************************/
// function name  mmcm_rwrap
// purpose        MMCM wrap function for R
// argument       *rdat        - for rdat.clsrnd (class information; before resampling)
//                             - clsrnd[sample_size]
//                *param       - for rdat.param (measurement)
//                             - param[sample_size]
//                *ctr_mat     - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                *resamp_size - No. of resampling repetition
//                *seed        - seed for psude-random number
//                *class_dim   - No. of groups
//                *contr_dim   - No. of contrast
//                *sample_size - whole sample size
// return value   *pval        - resampling P-Value (count)
//                             - p[contr_dim]
//**********************************************************************************/
void mmcm_rwrap(double *param, int *clsrnd, double *ctr_mat, long *resamp_size,
    long *seed, int *class_dim, int *contr_dim, int *sample_size,
    long *pval) {

  int i;
  struct mmcmdat *rdat;

  for (i = 0; i < *contr_dim; i++) {
    pval[i] = 0;
  }

  if ((rdat = (struct mmcmdat *) malloc(sizeof(struct mmcmdat) * *sample_size)) == NULL) {
    REprintf("Error malloc for rdat\n");
    exit(1);
  }

  for (i = 0; i < *sample_size; i++) {
    rdat[i].param = param[i];
    rdat[i].clsrnd = clsrnd[i];
  }

  // !!!!!!!!!!
  // calculate resampling P-Value (count) for MMCM
  //   argument
  //     dataset, coefficient matrix, resampling size, seed
  //     No. of class, No. of contrast, sample size
  stat_resamp(rdat, ctr_mat, *resamp_size, *seed, *class_dim, *contr_dim,
  *sample_size, pval);

  free(rdat);

}

//**********************************************************************************/
// function name   stat_resamp
// purpose         calculate resampling P-Value (count) for MMCM
// argument        *rdat       - dataset
//                             - rdat.clsrnd (class information; before resampling)
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 *ctr_mat    - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                 resamp_size - No. of resampling repetition
//                 seed        - seed for psude-random number
//                 class_dim   - No. of groups
//                 contr_dim   - No. of contrast
//                 sample_size - whole sample size
// return value    *pval       - resampling P-Value (count)
//                             - p[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int stat_resamp(struct mmcmdat *rdat, double *ctr_mat, long resamp_size,
    long seed, int class_dim, int contr_dim, int sample_size, long *pval) {

  int i, class_max = -1;
  long j, *class_size;
  double all_sum = 0.0, *ctr_denom, *class_mean, *t_d, *rt_d;

  if ((class_size = (long *) malloc(sizeof(long) * class_dim)) == NULL) {
    REprintf("Error malloc for class_size\n");
    exit(EXIT_FAILURE);
  }
  if ((class_mean = (double *) malloc(sizeof(double) * class_dim)) == NULL) {
    REprintf("Error malloc for class_mean\n");
    exit(EXIT_FAILURE);
  }
  if ((ctr_denom = (double *) malloc(sizeof(double) * contr_dim)) == NULL) {
    REprintf("Error malloc for ctr_denom\n");
    exit(EXIT_FAILURE);
  }
  if ((t_d = (double *) malloc(sizeof(double) * contr_dim)) == NULL) {
    REprintf("Error malloc for t_d\n");
    exit(EXIT_FAILURE);
  }
  if ((rt_d = (double *) malloc(sizeof(double) * contr_dim)) == NULL) {
    REprintf("Error malloc for rt_d\n");
    exit(EXIT_FAILURE);
  }

  // calculate denominator of modified contrast statistics
  stat_denom(ctr_mat, class_dim, contr_dim, ctr_denom);

  // calculate
  //   class of the maximum sample size
  //   sample_size of each class
  //   sum of all measurement
  qsort(rdat, sample_size, sizeof(struct mmcmdat), clsrnd_cmp);
  gsize_chk(rdat, class_dim, sample_size, &class_max, class_size, &all_sum);

  // calculate modified contrast statistics (sample)
  //   \bar{Y} = class_mean
  //   T'_k = frac{c^t_k \bar{Y}}{\sqrt{c^t_k c_k}}
  mean_get(rdat, class_dim, sample_size, class_max, class_size, all_sum,
      class_mean);
  stat_calc(ctr_mat, ctr_denom, class_mean, class_dim, contr_dim, t_d);

  // if openmp parallelization
#ifdef _OPENMP
#pragma omp parallel shared(pval)
  {
#pragma omp for firstprivate(rdat,t_d,rt_d,i,class_mean)
#endif
  // resampling (resamp_size repetition)
  init_gen_rand(seed);
  for (j = 0; j < resamp_size; j++) {

    // generate psude-random number
    for (i = 0; i < sample_size; i++) {
      rdat[i].clsrnd = gen_rand32();
    }
    
    // calculate modified contrast statistics (resampling)
    //   approach 1 (this program)
    //     rmean_get -> stat_calc
    //   approach 2 (common method)
    //     sort by rdat[i].clsrnd -> mean_get -> stat_calc
    rmean_get(rdat, class_dim, sample_size, class_max, class_size, all_sum,
        class_mean);
    stat_calc(ctr_mat, ctr_denom, class_mean, class_dim, contr_dim, rt_d);

    // calculate modified maximum contrast statistics (resampling)
    //   T'_{max} = \max \{ T'_k \}
    for (i = 1; i < contr_dim; i++) {
      rt_d[0] = MYMAX(rt_d[0], rt_d[i]);
    }

    // compare "resampling" T'_{max} > "sample" T'_i
    //   P-Value = pval[i] / resamp_size
    //   return value are counts of "resampling" T'_{max} > "sample" T'_i
    // TODO 080113
        //   return value are counts of "resampling" T'_{max} > "sample" T'_{max}
    for (i = 0; i < contr_dim; i++) {
      if (rt_d[0] >= t_d[i])
        pval[i]++;
    }
  }
#ifdef _OPENMP
}
#endif

  free(class_size);
  free(ctr_denom);
  free(class_mean);
  free(t_d);
  free(rt_d);

  return 0;

}

//**********************************************************************************/
// function name   stat_calc
// purpose         calculate "modified contrast statistics"
// argument        *ctr_mat    - coefficient matrix of contrast statistics
//                             - ctr_mat[class_dim * contr_dim]
//                 *ctr_denom  - denominator of modified contrast statistics
//                             - ctr_denom[contr_dim]
//                 *class_mean - "sample" or "resampling" mean of each class
//                             - class_mean[class_dim]
//                 class_dim   - No. of groups
//                 contr_dim   - No. of contrast
// return value    *t_d        - modified contrast statistics
//                             - t_d[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int stat_calc(double *ctr_mat, double *ctr_denom, double *class_mean,
    int class_dim, int contr_dim, double *t_d) {

  int i, j;
  double sum;

  // modified contrast statistics
  // T'_k = frac{c^t_k \bar{Y}}{\sqrt{c^t_k c_k}}
  for (i = 0; i < contr_dim; i++) {
    sum = 0;
    for (j = 0; j < class_dim; j++) {
      sum += class_mean[j] * ctr_mat[class_dim * i + j];
    }

    t_d[i] = fabs(sum) / ctr_denom[i];

  }

  return 0;

}

//**********************************************************************************/
// function name   stat_denom
// purpose         calculate denominator of "modified contrast statistics"
// argument        *ctr_mat   - coefficient matrix of contrast statistics
//                            - ctr_mat[class_dim * contr_dim]
//                 class_dim  - No. of groups
//                 contr_dim  - No. of contrast
// return value    *ctr_denom - denominator of modified contrast statistics
//                            - ctr_denom[contr_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int stat_denom(double *ctr_mat, int class_dim, int contr_dim, double *ctr_denom) {

  int i, j;
  double sum;

  // denominator of modified contrast statistics
  // \sqrt{c^t_k c_k} = \sum_{i=1}^{contr_dim} c_{i}^2
  for (i = 0; i < contr_dim; i++) {
    sum = 0;
    for (j = 0; j < class_dim; j ++) {
      sum += pow(ctr_mat[class_dim * i + j], 2);
    }
    ctr_denom[i] = sqrt(sum);
  }

  return 0;

}

//**********************************************************************************/
// function name   gsize_chk
// purpose         calculate sample size of each class
//                 calculate class has the maximum sample size
//                 calculate sum of all measurements
// argument        *rdat      - dataset
//                            - rdat.clsrnd (class information; before resampling)
//                            - rdat.param (measurement)
//                            - rdat[sample_size]
//                class_dim   - No. of groups
//                sample_size - whole sample size
// return value   *class_max  - class of the maximum sample size
//                *class_size - sample_size of each class
//                            - class_size[class_dim]
//                *all_sum    - sum of all measurement
// error code     always 0 yet, this version
//**********************************************************************************/
int gsize_chk(struct mmcmdat *rdat, int class_dim, int sample_size,
    int *class_max, long *class_size, double *all_sum) {

  long i, max = 0, tmp = 1;
  double sum = 0;

  // calculate class_size[i] & all_sum
  for (i = 0; i < sample_size; i++) {
    sum += rdat[i].param;
    if (tmp != rdat[i].clsrnd) {
      class_size[tmp-1] = i;
      tmp++;
    }
  }
  class_size[class_dim-1] = sample_size;

  *all_sum = sum;

  // calculate class_max & cnvert class_size[i]
  for (i = class_dim - 1; i > 0; i--) {
    class_size[i] = class_size[i] - class_size[i-1];
    if (max < class_size[i]) {
      max = class_size[i];
      *class_max = (int)i;
    }
  }
  if (max < class_size[0]) {
    *class_max = 0;
  }

  return 0;

}

//**********************************************************************************/
// function name   mean_get
// purpose         calculate "sample" mean of each class (class_mean[i])
// argument        *rdat       - dataset
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 class_dim   - No. of groups
//                 sample_size - whole sample size
//                 class_max   - class of the maximum sample size
//                 *class_size - sample_size of each class
//                             - class_size[class_dim]
//                 all_sum     - sum of all measurement
// return value    *class_mean - "sample" mean of each class
//                             - class_mean[class_dim]
// error code    always 0 yet, this version
//**********************************************************************************/
int mean_get(struct mmcmdat *rdat, int class_dim, int sample_size,
    int class_max, long *class_size, double all_sum, double *class_mean) {

  int i;
  long j, nsum = sample_size, start = 0, end = 0;
  double sum = 0, dsum = all_sum;

  // calculate [i] class mean by class_size[i]
  for (i = 0; i < class_dim; i++) {
    end += class_size[i];
    sum = 0;
    if (i != class_max) {
      for (j = start; j < end; j++) {
        sum += rdat[j].param;
      }
      class_mean[i] = sum / class_size[i];
      dsum -= sum;
      nsum -= class_size[i];
    }
    start += class_size[i];
  }

  class_mean[class_max] = dsum / nsum;

  return 0;

}

//**********************************************************************************/
// function name   rmean_get
// purpose         calculate "resampling" mean of each class (class_mean[i])
// argument        *rdat       - dataset
//                             - rdat.clsrnd (pseudo-random number for resampling)
//                             - rdat.param (measurement)
//                             - rdat[sample_size]
//                 class_dim   - No. of groups
//                 sample_size - whole sample size
//                 class_max   - class of the maximum sample size
//                 *class_size - sample_size of each class
//                             - class_size[class_dim]
//                 all_sum     - sum of all measurement
// return value    *class_mean - "resampling" mean of each class
//                             - class_mean[class_dim]
// error code      always 0 yet, this version
//**********************************************************************************/
int rmean_get(struct mmcmdat *rdat, int class_dim, int sample_size,
    int class_max, long *class_size, double all_sum, double *class_mean) {

  int i;
  long j, k, pre_param = -1, now_param = 0, nsum = sample_size;
  double sum, dsum = all_sum;

  for (i = 0; i < class_dim; i++) {

    sum = 0;

    if (i != class_max) {

      for (j = 0; j < class_size[i]; j++) {

        // only first sample
        // select maximum pseudo-random number
        if (pre_param == -1) {

          for (k = 0; k < sample_size; k++) {
            if (rdat[now_param].clsrnd < rdat[k].clsrnd)
              now_param = k;
          }
        }
        // other sample
        // select next largest pseudo-random number
        else {
          for (k = 0; k < sample_size; k++) {
            if (rdat[pre_param].clsrnd > rdat[k].clsrnd) {
              now_param = k;
              break;
            }
          }
          for (k = 0; k < sample_size; k++) {
            if (rdat[now_param].clsrnd < rdat[k].clsrnd
                && rdat[pre_param].clsrnd > rdat[k].clsrnd)
              now_param = k;
          }
        }
        sum += rdat[now_param].param;
        pre_param = now_param;
      }

      // calculate [i] class mean
      class_mean[i] = sum / class_size[i];

      dsum -= sum;
      nsum -= class_size[i];

    }
  }

  // calculat last class mean (class has the maximum sample size)
  // reduce calculation amount
  class_mean[class_max] = dsum / nsum;

  return 0;

}

//**********************************************************************************/
// function name  clsrnd_cmp
// purpose        comparison function in qsort()
//                sort by mmcmdat.clsrnd value in ascending order
//**********************************************************************************/
int clsrnd_cmp(const void *_p0, const void *_p1) {
  struct mmcmdat *p0 = (struct mmcmdat *)_p0;
  struct mmcmdat *p1 = (struct mmcmdat *)_p1;
  if (p0->clsrnd < p1->clsrnd)
    return -1; // in ascending order
  else if (p0->clsrnd > p1->clsrnd)
    return 1; // in ascending order
  else
    return 0;
}
