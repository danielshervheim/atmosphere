// Daniel Shervheim, 2019
// danielshervheim.com

#include "coefficients.h"

#include <cmath>

// These have to be re-declared in here for some reason...
constexpr double Coefficients::kSpectralIrradianceTable[];
constexpr double Coefficients::kOzoneAbsorptionCrossSections[];
constexpr double Coefficients::kCie2DegColorMatchingTable[];

double Coefficients::GetRayleighAbsorptionCoefficient(double lambda)
{
    // Molecules do not absorb light, so the coefficient is zero.
    return 0.0;
}

double Coefficients::GetRayleighScatteringCoefficient(double lambda)
{
    double refractive_index = 1.0003;
    double molecular_density = 2.545e25;
    double lambda_meters = lambda * 1e-9;

    // Source:
    // http://publications.lib.chalmers.se/records/fulltext/203057/203057.pdf
    double coefficient = 8.0 * pow(M_PI, 3.0);
    coefficient *= pow(pow(refractive_index, 2.0) - 1.0, 2.0);
    coefficient /= (3.0 * molecular_density * pow(lambda_meters, 4.0));
    return coefficient;
}

double Coefficients::GetRayleighExtinctionCoefficient(double lambda)
{
    double absorption = GetRayleighAbsorptionCoefficient(lambda);
    double scattering = GetRayleighScatteringCoefficient(lambda);
    return absorption + scattering;
}


double Coefficients::GetMieAbsorptionCoefficient(double lambda)
{
    // TODO: justify this choice.
    return GetMieScatteringCoefficient(lambda) / 9.0;
}

double Coefficients::GetMieScatteringCoefficient(double lambda)
{
    // Particles scatter light equally, with no dependance on wavelength.
    // Hence the constant coefficient.
    return 2e-6;
}

double Coefficients::GetMieExtinctionCoefficient(double lambda)
{
    double absorption = GetMieAbsorptionCoefficient(lambda);
    double scattering = GetMieScatteringCoefficient(lambda);
    return absorption + scattering;
}


double Coefficients::GetOzoneAbsorptionCoefficient(double lambda)
{
    // Convert the wavelength into an index to sample the table.
    int index = (int)lambda - 215;
    if (index < 0 || index >= 785)
    {
        // TODO: throw error?
        return 0.0;
    }

    return kOzoneAbsorptionCrossSections[index];
}

double Coefficients::GetOzoneScatteringCoefficient(double lambda)
{
    // Ozone does not scatter light, so the coefficient is zero.
    return 0.0;
}

double Coefficients::GetOzoneExtinctionCoefficient(double lambda)
{
    double absorption = GetOzoneAbsorptionCoefficient(lambda);
    double scattering = GetOzoneScatteringCoefficient(lambda);
    return absorption + scattering;
}

double Coefficients::SampleCie2DegColorMatchingTable(double lambda, int column)
{
    if (lambda < 360.0 || lambda >= 830.0)
    {
        // TODO: error?
        return 0.0;
    }

    double u = (lambda - 360) / 5.0;
    int row = floor(u);
    u -= row;
    double a = kCie2DegColorMatchingTable[4 * row + column];
    double b = kCie2DegColorMatchingTable[4 * (row + 1) + column];
    return a + (b-a)*u;
}

double Coefficients::GetSpectralIrradiance(double lambda)
{
    int index = (int)lambda - 280;
    if (index < 0 || index >= 721)
    {
        // TODO: error?
        return 0.0;
    }
    return kSpectralIrradianceTable[index];
}
