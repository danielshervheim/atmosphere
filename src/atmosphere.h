#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "vec2.h"

class Atmosphere
{
public:
    // Precomputation parameters.
    int table_dimension_ = 64;
    double lambda_r_ = 650.0;
    double lambda_g_ = 510.0;
    double lambda_b_ = 475.0;
    bool normalize_precomputation_results_ = true;

    // Integration parameters.
    int view_ray_integration_steps_ = 64;
    int light_ray_integration_steps_ = 64;
    int spectral_to_rgb_integration_steps_ = 32;

    // Spectral to RGB conversion parameters.
    double alpha_ = 2.0;  // TODO: explain what this is.

    // Planetary properties.
    double radius_planet_ = 6371000.0;
    double radius_atmosphere_ = 6471000.0;

    double scale_height_rayleigh_ = 8000.0;  // TODO: explain what these are.
    double scale_height_mie_ = 1200.0;

    Atmosphere();
    ~Atmosphere();

    void PrecomputeTable();

    int GetTableLength();
    double* GetPrecomputedRayleighTable();
    double* GetPrecomputedMieTable();

    void GetNormalizationFactorsRayleigh(double& min, double& max);
    void GetNormalizationFactorsMie(double& min, double& max);

    void GetSpectralToRGBConversionConstants(double& k_r, double& k_g, double& k_b);
private:
    double* precomputed_rayleigh_table_ = nullptr;
    double* precomputed_mie_table_ = nullptr;

    double normalization_min_rayleigh_;
    double normalization_max_rayleigh_;

    double normalization_min_mie_;
    double normalization_max_mie_;

    void PrecomputeTableCell(vec2 view_dir, vec2 light_dir, double lambda, double& rayleigh, double& mie);

    void NormalizeTable();
};

#endif  // ATMOSPHERE_H_
