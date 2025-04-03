#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

void butterworth_lowpass(int order, double fc, double fs, double **b, double **a);
void filtfilt(double *b, double *a, int order, double *input, double *output, int n);

#endif
