// Daniel Shervheim, 2019
// danielshervheim.com

#include "atmosphere.h"
#include "coefficients.h"
#include "utilities.h"

#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
    Atmosphere atmosphere;
    std::string output_path = "";

    // TODO: parse the command line input for additional parameters.
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-o") == 0)
        {
            if (i+1 < argc)
            {
                output_path = argv[i+1];
            }
            else
            {
                // TODO: Error.
            }
        }
    }

    atmosphere.PrecomputeTable();

    // Convert the Rayleigh table to a float array and write it to a file.
    float* rayleigh = Utilities::ConvertDoubleArrayToFloatArray(
        atmosphere.GetPrecomputedRayleighTable(),
        atmosphere.GetTableLength()
    );
    if (!Utilities::WriteFloatArrayToFile(output_path + "rayleigh.bin", rayleigh, atmosphere.GetTableLength()))
    {
        std::cout << "ERROR: Could not write rayleigh table to file." << std::endl;
    }

    // Convert the Mie table to a float array and write it to a file.
    float* mie = Utilities::ConvertDoubleArrayToFloatArray(
        atmosphere.GetPrecomputedMieTable(),
        atmosphere.GetTableLength()
    );
    if (!Utilities::WriteFloatArrayToFile(output_path + "mie.bin", mie, atmosphere.GetTableLength()))
    {
        std::cout << "ERROR: Could not write mie table to file." << std::endl;
    }

    delete [] rayleigh;
    delete [] mie;

    double minR, maxR, minM, maxM;
    atmosphere.GetNormalizationFactorsRayleigh(minR, maxR);
    atmosphere.GetNormalizationFactorsMie(minM, maxM);

    double r, g, b;
    atmosphere.GetSpectralToRGBConversionConstants(r, g, b);

    // Create a markdown file to write information about the precomputation into.
    std::ofstream info(output_path + "results.txt", std::ios::out);
    if(!info)
    {
        std::cout << "ERROR: Could not write information to file." << std::endl;
    }

    // If normalization was enabled, write the formulas needed to reconstruct
    // the correct values from the table into the file.
    if (atmosphere.normalize_precomputation_results_)
    {
        info << "NORMALIZATION" << std::endl;
        info << std::scientific << std::endl;

        info << "rayleigh = rayleigh_from_table*(" << maxR << " - " << minR << ") + " << minR << std::endl;
        info << "mie = mie_from_table*(" << maxM << " - " << minM << ") + " << minM << std::endl;
        info << std::endl;
    }

    // Write the spectral to RGB conversion constants into the table.
    info << "SPECTRAL TO RGB" << std::endl;
    info << std::endl;

    info << "R = " << r << std::endl;
    info << "G = " << g << std::endl;
    info << "B = " << b << std::endl;
    info << std::endl;

    // Write the spectral irradiance values into the table.
    info << "SPECTRAL IRRADIANCE" << std::endl;
    info << std::endl;

    info << "R = " << Coefficients::GetSpectralIrradiance(atmosphere.lambda_r_) << std::endl;
    info << "G = " << Coefficients::GetSpectralIrradiance(atmosphere.lambda_g_) << std::endl;
    info << "B = " << Coefficients::GetSpectralIrradiance(atmosphere.lambda_b_) << std::endl;

    // Close the file.
    info.close();
}
