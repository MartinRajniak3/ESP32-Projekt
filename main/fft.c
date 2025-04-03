#include "fft.h"
#include <math.h>
#include <stdlib.h>

void fft(double *x_in, double *y_real, double *y_imag, int N) {
    if (N <= 1) return;

    // Rozdelenie na párne a nepárne indexy
    double *even_real = (double *)malloc(N / 2 * sizeof(double));
    double *even_imag = (double *)malloc(N / 2 * sizeof(double));
    double *odd_real = (double *)malloc(N / 2 * sizeof(double));
    double *odd_imag = (double *)malloc(N / 2 * sizeof(double));

    for (int i = 0; i < N / 2; i++) {
        even_real[i] = x_in[2 * i];
        even_imag[i] = 0.0;
        odd_real[i] = x_in[2 * i + 1];
        odd_imag[i] = 0.0;
    }

    fft(even_real, even_real, even_imag, N / 2);
    fft(odd_real, odd_real, odd_imag, N / 2);

    // Vypočítanie FFT
    for (int k = 0; k < N / 2; k++) {
        double t_real = cos(-2 * M_PI * k / N) * odd_real[k] - sin(-2 * M_PI * k / N) * odd_imag[k];
        double t_imag = sin(-2 * M_PI * k / N) * odd_real[k] + cos(-2 * M_PI * k / N) * odd_imag[k];

        y_real[k] = even_real[k] + t_real;
        y_imag[k] = even_imag[k] + t_imag;
        y_real[k + N / 2] = even_real[k] - t_real;
        y_imag[k + N / 2] = even_imag[k] - t_imag;
    }

    free(even_real);
    free(even_imag);
    free(odd_real);
    free(odd_imag);
}
