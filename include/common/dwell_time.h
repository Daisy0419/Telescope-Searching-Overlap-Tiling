#pragma once

// AIRIS Dwelltime
// float dwell_times_airis(float theta, float flux = 3.0, int snr = 5, float sig_attenuation=1.0);
float dwell_times_airis(float theta, float flux = 4.0, int snr = 5.0, float sig_attenuation=1.0);


// AIRMASS-ONLY Dwelltime
double dwell_time_airmass(double zenithAngle, double dwelltime_zenith = 1.0);