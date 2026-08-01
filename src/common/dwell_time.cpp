#include "dwell_time.h"
#include <cmath>
#include <iostream>
#include <vector>

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)


// AIRIS Dwell Models
// Returns the flux in eV the 645-675 nm band given Band function parameters
// Assumes amplitude A is in ph/cm^2/s/keV
float airis_flux_from_band(float A, float E_peak, float alpha, float beta) {
    constexpr float lower_wavelength = 645.0; // nm
    constexpr float upper_wavelength = 675.0; // nm
    constexpr float hc = 1239.84193; // eV*nm
    constexpr float e_lower = hc / upper_wavelength / 1000; // keV
    constexpr float e_upper = hc / lower_wavelength / 1000; // keV
    constexpr int num_steps = 100;
    constexpr float energy_step = (e_upper - e_lower) / num_steps;

    float total_flux = 0.0;
    for (float e = e_lower + energy_step / 2;  e < e_upper; e += energy_step)
    {
        //From https://iopscience.iop.org/article/10.3847/1538-4357/abf24d/pdf Eqn. 3
        // float band_value = A *
        //     std::pow(e / 100, beta) *
        //     std::exp(beta - alpha) *
        //     std::pow((alpha - beta) * E_peak / (100 * (2 + alpha)), alpha - beta);

        //From https://articles.adsabs.harvard.edu/pdf/1993ApJ...413..281B Eqn. 1
        float E_0 = E_peak / (2 + alpha);
        float band_value = A *
            std::pow(e / 100, alpha) * std::exp(-e / E_0);

        float photon_flux = band_value * energy_step;
        float energy_flux = photon_flux * e * 1000; // convert keV to eV
        total_flux += energy_flux;
    }
    
    // std::cout << "Total flux in band: " << total_flux << " eV/cm^2/s\n";
    return total_flux;

}

// Apply atmospheric attenuation based on polar angle theta in degrees
float airis_attenuate_flux(float flux, float theta, float sig_attenuation=1.0) {
    constexpr float altitude = 40000.0; // in meters
    constexpr float scale_height = 8435.0; // in meters
    float relative_pressure = std::exp(-altitude / scale_height);
    float air_mass = relative_pressure / (
                     std::cos(theta * M_PI / 180.0) +
                     0.50572 * std::pow(96.07995 - theta, -1.6364) );
    float attenuation_factor = std::exp(-0.04558 * air_mass);
    float attenuated_flux = flux * attenuation_factor * sig_attenuation;
    
    // std::cout << "Total attenuated flux in band: " << attenuated_flux << " eV/cm^2/s\n";

    return attenuated_flux;
}

// Convert flux in eV/cm^2/s to magnitude
float airis_flux_to_magnitude(float flux) {
    constexpr float bandpass_width = 300.0; // in Angstroms
    constexpr float zero_point_magnitude = 1198.36974; // in eV/cm^2/s/Angstrom
    constexpr float zero_point_flux = zero_point_magnitude * bandpass_width; // in eV/cm^2/s
    float magnitude = -2.5 * std::log10(flux / zero_point_flux);
    // std::cout << "Computed magnitude: " << magnitude << "\n";
    return magnitude;
}

float dwell_times_airis(float theta, float flux, int snr, float sig_attenuation) {
    // float flux = airis_flux_from_band(A, E_peak, alpha, beta);
    float attenuated_flux = airis_attenuate_flux(flux, theta, sig_attenuation);
    float magnitude = airis_flux_to_magnitude(attenuated_flux);

    return std::pow(
        snr / 897603.0 *
        std::pow(10.0, 0.4 * magnitude * 1.06596),
        1.60636
    );
}



// AIR-MASS-Only Dwell Models
// Computes air mass using the Kasten & Young model
double computeAirMass(double zenithAngleDegrees) {
    if (zenithAngleDegrees >= 90.0) {
        std::cerr << "Zenith angle must be less than 90 degrees.";
    }
    double theta = DEG_TO_RAD * zenithAngleDegrees;
    double numerator = 1.0;
    double denominator = std::cos(theta) + 0.50572 * std::pow(96.07995 - zenithAngleDegrees, -1.6364);
    return numerator / denominator;
}

// Computes the light attenuation scaling factor for a given air mass
double computeScalingFactor(double airMass) {
    return 1.1129 * std::exp(-0.107 * airMass);
}

// Computes the scaled dwell time needed to achieve same SNR as reference dwell time t_ref at s=1
double computeDwellTime(double t_ref, double scalingFactor) {
    return t_ref / (scalingFactor * scalingFactor);
}


double dwell_time_airmass(double zenithAngle, double dwelltime_zenith) {
    if (zenithAngle >= 90.0) {
        std::cerr << "Zenith angle must be less than 90 degrees.";
    }
    double airMass = computeAirMass(zenithAngle);
    double scalingFactor = computeScalingFactor(airMass);
    double dwell_time = computeDwellTime(dwelltime_zenith, scalingFactor);
    return dwell_time;
}
