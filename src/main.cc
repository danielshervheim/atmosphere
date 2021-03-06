// Daniel Shervheim, 2019
// danielshervheim.com

#include "atmosphere.h"
#include "coefficients.h"
#include "utilities.h"

#include <iostream>
#include <fstream>

#include <ImfRgbaFile.h>
#include <ImfRgba.h>
#include <half.h>

int main(int argc, char* argv[])
{
    Atmosphere atmosphere;
    std::string output_path = "";

    bool writeToExr = false;

    // TODO: parse the command line input for additional parameters.
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-o") == 0)
        {
            if (i+1 < argc)
            {
                output_path = argv[i+1];

                // Verify output path ends in /
                if (argv[i+1][output_path.size()-1] != '/')
                {
                    std::cout << "ERROR: provided ouput path must end in \"/\"." << std::endl;
                    return -1;
                }

                // Verify output path exists.
                if (!Utilities::DirectoryExists(argv[i+1]))
                {
                    std::cout << "ERROR: provided ouput path does not exist." << std::endl;
                    return -1;
                }
            }
            else
            {
                std::cout << "ERROR: -o must be followed by a space and a valid path." << std::endl;
                return -1;
            }
        }

        if (strcmp(argv[i], "-exr") == 0)
        {
            writeToExr = true;
        }

        if (strcmp(argv[i], "-n") == 0)
        {
            atmosphere.normalize_precomputation_results_ = true;
        }
    }

    atmosphere.PrecomputeTable();

    if (writeToExr)
    {
        const int DIM = atmosphere.table_dimension_;
        Imf::Rgba *pixels = new Imf::Rgba[DIM * DIM];

        // Write the rayleigh table to an OpenEXR file.
        double* rayleigh = atmosphere.GetPrecomputedRayleighTable();
        char rayleigh_path[output_path.size() + 12 + 1];
        strcpy(rayleigh_path, (output_path + "rayleigh.exr").c_str());
        for (int i = 0; i < DIM*DIM*3; i += 3)
        {
            int j = i / 3;
            pixels[j] = Imf::Rgba(rayleigh[i+0], rayleigh[i+1], rayleigh[i+2], 1.0);
        }

        try
        {
            Imf::RgbaOutputFile rayleigh_file(rayleigh_path, DIM, DIM, Imf::WRITE_RGBA);
            rayleigh_file.setFrameBuffer(pixels, 1, DIM);
            rayleigh_file.writePixels(DIM);
        }
        catch (const std::exception &exc)
        {
            (void)exc;
            std::cout << "ERROR: failed to write rayleigh table to .exr file." << std::endl;
        }

        // Write the mie table to an OpenEXR file.
        double* mie = atmosphere.GetPrecomputedMieTable();
        char mie_path[output_path.size() + 7 + 1];
        strcpy(mie_path, (output_path + "mie.exr").c_str());
        for (int i = 0; i < DIM*DIM*3; i += 3)
        {
            int j = i / 3;
            pixels[j] = Imf::Rgba(mie[i+0], mie[i+1], mie[i+2], 1.0);
        }
        try
        {
            Imf::RgbaOutputFile mie_file(mie_path, DIM, DIM, Imf::WRITE_RGBA);
            mie_file.setFrameBuffer(pixels, 1, DIM);
            mie_file.writePixels(DIM);
        }
        catch (const std::exception &exc)
        {
            (void)exc;
            std::cout << "ERROR: failed to write rayleigh table to .exr file." << std::endl;
        }

        delete [] pixels;

        // Write transmittance to file.
        pixels = new Imf::Rgba[DIM];
        double* transmittance = atmosphere.GetPrecomputedTransmittanceTable();
        char transmittance_path[output_path.size() + 17 + 1];
        strcpy(transmittance_path, (output_path + "transmittance.exr").c_str());
        for (int i = 0; i < DIM*3; i += 3)
        {
            int j = i / 3;
            pixels[j] = Imf::Rgba(transmittance[i+0], transmittance[i+1], transmittance[i+2], 1.0);
        }
        try
        {
            Imf::RgbaOutputFile transmittance_file(transmittance_path, DIM, 1, Imf::WRITE_RGBA);
            transmittance_file.setFrameBuffer(pixels, 1, DIM);
            transmittance_file.writePixels(1);
        }
        catch (const std::exception &exc)
        {
            (void)exc;
            std::cout << "ERROR: failed to write transmittance table to .exr file." << std::endl;
        }
        delete [] pixels;
    }
    else
    {
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
    }

    double minR, maxR, minM, maxM;
    atmosphere.GetNormalizationFactorsRayleigh(minR, maxR);
    atmosphere.GetNormalizationFactorsMie(minM, maxM);

    double r, g, b;
    atmosphere.GetSpectralToRGBConversionConstants(r, g, b);

    // Create a markdown file to write information about the precomputation into.
    std::ofstream info(output_path + "results.txt", std::ios::out);
    if(!info)
    {
        std::cout << "ERROR: Could not write precomputed parameter results to file." << std::endl;
        return -1;
    }

    // If normalization was enabled, write the formulas needed to reconstruct
    // the correct values from the table into the file.
    if (atmosphere.normalize_precomputation_results_)
    {
        info << "NORMALIZATION" << std::endl;
        info << std::scientific << std::endl;

        info << "Rayleigh min, max = (" << minR << ", " << maxR << ")" << std::endl;
        info << "Mie min, max = (" << minM << ", " << maxM << ")" << std::endl;

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

    return 0;
}
