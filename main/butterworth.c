#include "butterworth.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Funkcia na výpočet Butterworth koeficientov
void butterworth_lowpass(int order, double fc, double fs, double **b, double **a) {
    double wc = 2.0 * M_PI * fc / fs;
    double wc_tan = tan(wc / 2.0);
    
    *b = (double *)malloc((order + 1) * sizeof(double));
    *a = (double *)malloc((order + 1) * sizeof(double));

    (*b)[0] = wc_tan / (1.0 + wc_tan);
    (*b)[1] = (*b)[0];
    (*a)[0] = 1.0;
    (*a)[1] = (wc_tan - 1.0) / (1.0 + wc_tan);

    for (int i = 2; i <= order; i++) {
        (*b)[i] = 0;
        (*a)[i] = 0;
    }
}

// Funkcia na aplikáciu filtra (dopredný a spätný prechod - filtfilt)
void filtfilt(double *b, double *a, int order, double *input, double *output, int n) {
    double *y = (double *)malloc(n * sizeof(double));
    memset(y, 0, n * sizeof(double));

    for (int i = 0; i < n; i++) {
        y[i] = b[0] * input[i];
        for (int j = 1; j <= order; j++) {
            if (i >= j) {
                y[i] += b[j] * input[i - j] - a[j] * y[i - j];
            }
        }
    }

    for (int i = 0; i < n; i++) {
        output[i] = y[n - i - 1]; // Obrátenie výsledku pre filtfilt
    }

    free(y);
}
