// Daniel Shervheim, 2019
// danielshervheim.com

#include "atmosphere.h"

#include "coefficients.h"
#include "utilities.h"

#include <cmath>  // acos
#include <iostream>  // cout
#include <stdexcept>  // runtime_error

Atmosphere::Atmosphere()
{

}

Atmosphere::~Atmosphere()
{
    if (precomputed_rayleigh_table_ != nullptr)
    {
        delete [] precomputed_rayleigh_table_;
    }

    if (precomputed_mie_table_ != nullptr)
    {
        delete [] precomputed_mie_table_;
    }
}

void Atmosphere::PrecomputeTable()
{
    const int kTableLength = GetTableLength();

    if (precomputed_rayleigh_table_ != nullptr)
    {
        delete [] precomputed_rayleigh_table_;
    }
    precomputed_rayleigh_table_ = new double[kTableLength];

    if (precomputed_mie_table_ != nullptr)
    {
        delete [] precomputed_mie_table_;
    }
    precomputed_mie_table_ = new double[kTableLength];

    for (int i = 0; i < kTableLength; i += 3)
    {
        int j = i / 3;

        // Convert the 1D index to 2D indices.
        int x = j % table_dimension_;
        int y = j / table_dimension_;

        // For both, we use sin() for x and cos() for y, since we are
        // measuring the angle from the zenith, rather than the horizon.

        double u = (x+0.5) / (double)(table_dimension_);
        double theta_view_dir = acos(TextureCoordinateToCosTheta(u));
        vec2 view_dir = vec2(sin(theta_view_dir), cos(theta_view_dir));

        double v = (y+0.5) / (double)(table_dimension_);
        double theta_light_dir = acos(TextureCoordinateToCosTheta(v));
        vec2 light_dir = vec2(sin(theta_light_dir), cos(theta_light_dir));

        double rayleigh = 0.0;
        double mie = 0.0;

        PrecomputeTableCell(view_dir, light_dir, lambda_r_, rayleigh, mie);
        precomputed_rayleigh_table_[i+0] = rayleigh;
        precomputed_mie_table_[i+0] = mie;

        PrecomputeTableCell(view_dir, light_dir, lambda_g_, rayleigh, mie);
        precomputed_rayleigh_table_[i+1] = rayleigh;
        precomputed_mie_table_[i+1] = mie;

        PrecomputeTableCell(view_dir, light_dir, lambda_b_, rayleigh, mie);
        precomputed_rayleigh_table_[i+2] = rayleigh;
        precomputed_mie_table_[i+2] = mie;
    }

    if (normalize_precomputation_results_)
    {
        NormalizeTable();
    }

    // Compute transmittance
    if (precomputed_transmittance_table_ != nullptr)
    {
        delete [] precomputed_transmittance_table_;
    }
    precomputed_transmittance_table_ = new double[table_dimension_*3];

    for (int i = 0; i < table_dimension_*3; i += 3)
    {
        int j = i/3;
        double u = (j+0.5) / (double)(table_dimension_);
        double theta_view_dir = acos(TextureCoordinateToCosTheta(u));
        vec2 view_dir = vec2(sin(theta_view_dir), cos(theta_view_dir));

        double transmittance = 0.0;

        PrecomputeTransmittance(view_dir, lambda_r_, transmittance);
        precomputed_transmittance_table_[i+0] = transmittance;

        PrecomputeTransmittance(view_dir, lambda_g_, transmittance);
        precomputed_transmittance_table_[i+1] = transmittance;

        PrecomputeTransmittance(view_dir, lambda_b_, transmittance);
        precomputed_transmittance_table_[i+2] = transmittance;
    }
}

int Atmosphere::GetTableLength()
{
    return table_dimension_ * table_dimension_ * 3;
}

double* Atmosphere::GetPrecomputedRayleighTable()
{
    return precomputed_rayleigh_table_;
}

double* Atmosphere::GetPrecomputedMieTable()
{
    return precomputed_mie_table_;
}

double* Atmosphere::GetPrecomputedTransmittanceTable()
{
    return precomputed_transmittance_table_;
}

void Atmosphere::GetNormalizationFactorsRayleigh(double& min, double& max)
{
    min = normalization_min_rayleigh_;
    max = normalization_max_rayleigh_;
}

void Atmosphere::GetNormalizationFactorsMie(double& min, double& max)
{
    min = normalization_min_mie_;
    max = normalization_max_mie_;
}

void Atmosphere::GetSpectralToRGBConversionConstants(double& k_r, double& k_g, double& k_b)
{
    k_r = 0.0;
    k_g = 0.0;
    k_b = 0.0;

    double lambda = 360.0;
    const double kDeltaLambda = (830.0-360.0) / (double)spectral_to_rgb_integration_steps_;

    for (int i = 0; i < spectral_to_rgb_integration_steps_; i++)
    {
        double x_bar = Coefficients::SampleCie2DegColorMatchingTable(lambda, 1);
        double y_bar = Coefficients::SampleCie2DegColorMatchingTable(lambda, 2);
        double z_bar = Coefficients::SampleCie2DegColorMatchingTable(lambda, 3);

        double r_bar, g_bar, b_bar;
        Utilities::ConvertXyzToRgb(x_bar, y_bar, z_bar, r_bar, g_bar, b_bar);

        k_r += r_bar * (pow(lambda, -alpha_) / pow(lambda_r_, -alpha_)) * kDeltaLambda;
        k_g += g_bar * (pow(lambda, -alpha_) / pow(lambda_g_, -alpha_)) * kDeltaLambda;
        k_b += b_bar * (pow(lambda, -alpha_) / pow(lambda_b_, -alpha_)) * kDeltaLambda;

        lambda += kDeltaLambda;
    }
}

void Atmosphere::PrecomputeTableCell(vec2 view_dir, vec2 light_dir, double lambda, double& rayleigh, double& mie)
{
    rayleigh = 0.0;
    mie = 0.0;

    // `pa`, the viewing position. We offset the planet radius by the height
    // of an average human to get a typical viewing position.
    vec2 pa = vec2(0.0, radius_planet_ + 1.8);

    double dist_to_pb;
    double t_min, t_max;

    // Intersect the viewing ray with the atmosphere.
    if (!Utilities::RayCircleIntersection(pa, view_dir, vec2(0.0, 0.0), radius_atmosphere_, t_min, t_max))
    {
        std::cerr << "pa outside the atmosphere. ensure that the atmosphere radius is greater than the planet radius." << std::endl;
    }

    dist_to_pb = t_max;

    // Intersect the viewing ray with the ground. This ensures that we do not
    // integrate through the planet to the atmosphere on the opposite side.
    if (Utilities::RayCircleIntersection(pa, view_dir, vec2(0.0, 0.0), radius_planet_, t_min, t_max))
    {
        dist_to_pb = fmin(dist_to_pb, t_min);
    }

    // Compute `pb`, the intersection of the `view_dir` ray with the atmosphere
    // from `pa`.
    vec2 pb = pa + (view_dir * dist_to_pb);

    const double kViewDelta = vec2::magnitude(pb - pa) / (double)view_ray_integration_steps_;

    double optical_depth_rayleigh = 0.0;
    double optical_depth_mie = 0.0;
    double optical_depth_ozone = 0.0;

    // Integrate from `pa` to `pb`.
    for (int i = 0; i < view_ray_integration_steps_; i++)
    {
        double t = (double)i / (double)(view_ray_integration_steps_ - 1.0);
        vec2 p = vec2::lerp(pa, pb, t);

        // Compute the height of the point `p` from the surface of the planet.
        double height_p = vec2::magnitude(p) - radius_planet_;

        // NOTE: the ozone density function is an approximation. In reality,
        // ozone does not decrease exponentially with height.
        optical_depth_rayleigh += exp(-height_p / scale_height_rayleigh_) * kViewDelta;
        optical_depth_mie += exp(-height_p / scale_height_mie_) * kViewDelta;
        optical_depth_ozone += exp(-height_p / scale_height_rayleigh_) * 6e-7 * kViewDelta;

        // Compute the optical depth and transmittance from the `pa` to `p`.
        double optical_depth_pa_to_p =
            optical_depth_rayleigh * Coefficients::GetRayleighExtinctionCoefficient(lambda) +
            optical_depth_mie * Coefficients::GetMieExtinctionCoefficient(lambda) +
            optical_depth_ozone * Coefficients::GetOzoneExtinctionCoefficient(lambda);
        double transmittance_pa_to_p = exp(-optical_depth_pa_to_p);

        // Intersect the light ray with the atmosphere. We add 1 meter to the
        // atmosphere radius, otherwise the intersection test fails when p = pb).
        if (!Utilities::RayCircleIntersection(p, light_dir, vec2(0.0, 0.0), radius_atmosphere_+1.0, t_min, t_max))
        {
            std::cerr << "p is outside the atmosphere." << std::endl;
        }

        double dist_to_pc = t_max;
        double transmittance_p_to_pc = 0.0;

        // Intersect the light ray with the planet. We only want to integrate
        // from `p` to `pc` if the light ray is not blocked by the planet.
        if (!Utilities::RayCircleIntersection(p, light_dir, vec2(0.0, 0.0), radius_planet_, t_min, t_max))
        {
            // Compute `pc`, the intersection of the `light_dir` ray with the atmosphere
            // from `p`.
            vec2 pc = p + (light_dir * dist_to_pc);

            // Compute the transmittance between `p` and `pc`.
            const double kDeltaLight = vec2::magnitude(pc - p) / (double)light_ray_integration_steps_;
            double odr_p_to_pc = 0.0;
            double odm_p_to_pc = 0.0;
            double odo_p_to_pc = 0.0;

            for (int j = 0; j < light_ray_integration_steps_; j++)
            {
                double t_2 = (double)j / (double)(light_ray_integration_steps_ - 1.0);
                vec2 p_2 = vec2::lerp(p, pc, t_2);
                double height_p_2 = vec2::magnitude(p_2) - radius_planet_;

                odr_p_to_pc += exp(-height_p_2 / scale_height_rayleigh_) * kDeltaLight;
                odm_p_to_pc += exp(-height_p_2 / scale_height_mie_) * kDeltaLight;
                odo_p_to_pc += exp(-height_p_2 / scale_height_rayleigh_) * 6e-7 * kDeltaLight;
            }

            double optical_depth_p_to_pc =
                odr_p_to_pc * Coefficients::GetRayleighExtinctionCoefficient(lambda) +
                odm_p_to_pc * Coefficients::GetMieExtinctionCoefficient(lambda) +
                odo_p_to_pc * Coefficients::GetOzoneExtinctionCoefficient(lambda);
            transmittance_p_to_pc = exp(-optical_depth_p_to_pc);
        }

        double transmittance_pa_to_p_to_pc = transmittance_pa_to_p * transmittance_p_to_pc;

        rayleigh += exp(-height_p / scale_height_rayleigh_) * transmittance_pa_to_p_to_pc * kViewDelta;
        mie += exp(-height_p / scale_height_mie_) * transmittance_pa_to_p_to_pc * kViewDelta;;
    }

    rayleigh *= Coefficients::GetRayleighScatteringCoefficient(lambda) / (4.0 * M_PI);
    mie *= Coefficients::GetMieScatteringCoefficient(lambda) / (4.0 * M_PI);

    // We defer the inclusion of the phase and spectral intensity terms to
    // the pixel shader, because they are constant anyways.
}

void Atmosphere::PrecomputeTransmittance(vec2 view_dir, double lambda, double& transmittance)
{
    transmittance = 0.0;

    // `pa`, the viewing position. We offset the planet radius by the height
    // of an average human to get a typical viewing position.
    vec2 pa = vec2(0.0, radius_planet_ + 1.8);

    double dist_to_pb;
    double t_min, t_max;

    // Intersect the viewing ray with the atmosphere.
    if (!Utilities::RayCircleIntersection(pa, view_dir, vec2(0.0, 0.0), radius_atmosphere_, t_min, t_max))
    {
        std::cerr << "pa outside the atmosphere. ensure that the atmosphere radius is greater than the planet radius." << std::endl;
    }

    dist_to_pb = t_max;

    // Intersect the viewing ray with the ground. This ensures that we do not
    // integrate through the planet to the atmosphere on the opposite side.
    if (Utilities::RayCircleIntersection(pa, view_dir, vec2(0.0, 0.0), radius_planet_, t_min, t_max))
    {
        transmittance = 0.0;
        return;
    }

    // Compute `pb`, the intersection of the `view_dir` ray with the atmosphere
    // from `pa`.
    vec2 pb = pa + (view_dir * dist_to_pb);

    const double kViewDelta = vec2::magnitude(pb - pa) / (double)view_ray_integration_steps_;

    double optical_depth_rayleigh = 0.0;
    double optical_depth_mie = 0.0;
    double optical_depth_ozone = 0.0;

    // Integrate from `pa` to `pb`.
    for (int i = 0; i < view_ray_integration_steps_; i++)
    {
        double t = (double)i / (double)(view_ray_integration_steps_ - 1.0);
        vec2 p = vec2::lerp(pa, pb, t);

        // Compute the height of the point `p` from the surface of the planet.
        double height_p = vec2::magnitude(p) - radius_planet_;

        optical_depth_rayleigh += exp(-height_p / scale_height_rayleigh_) * kViewDelta;
        optical_depth_mie += exp (-height_p / scale_height_mie_) * kViewDelta;
        optical_depth_ozone += exp(-height_p / scale_height_rayleigh_) * 6e-7 * kViewDelta;
    }

    optical_depth_rayleigh *= Coefficients::GetRayleighExtinctionCoefficient(lambda);
    optical_depth_mie *= Coefficients::GetMieExtinctionCoefficient(lambda);
    optical_depth_ozone *= Coefficients::GetOzoneExtinctionCoefficient(lambda);

    double optical_depth = optical_depth_rayleigh + optical_depth_mie + optical_depth_ozone;

    transmittance = exp(-optical_depth);
}

void Atmosphere::NormalizeTable()
{
    const int kTableLength = GetTableLength();

    // Set the min and max to the first element in the array.
    normalization_min_rayleigh_ = precomputed_rayleigh_table_[0];
    normalization_max_rayleigh_ = precomputed_rayleigh_table_[0];
    normalization_min_mie_ = precomputed_mie_table_[0];
    normalization_max_mie_ = precomputed_mie_table_[0];

    // Find the minimum and maximum values in each array.
    for (int i = 0; i < kTableLength; i++)
    {
        normalization_min_rayleigh_ = fmin(normalization_min_rayleigh_, precomputed_rayleigh_table_[i]);
        normalization_max_rayleigh_ = fmax(normalization_max_rayleigh_, precomputed_rayleigh_table_[i]);
        normalization_min_mie_ = fmin(normalization_min_mie_, precomputed_mie_table_[i]);
        normalization_max_mie_ = fmax(normalization_max_mie_, precomputed_mie_table_[i]);
    }

    // Compute the range in values for each array.
    double range_rayleigh = normalization_max_rayleigh_ - normalization_min_rayleigh_;
    double range_mie = normalization_max_mie_ - normalization_min_mie_;

    // Remap each value from min...max to 0...1.
    for (int i = 0; i < kTableLength; i++)
    {
        precomputed_rayleigh_table_[i] = (precomputed_rayleigh_table_[i] +
            normalization_min_rayleigh_) / range_rayleigh;

        precomputed_mie_table_[i] = (precomputed_mie_table_[i] +
            normalization_min_mie_) / range_mie;
    }
}

double Atmosphere::CosThetaToTextureCoordinate(double cosTheta)
{
    // Favors precision near horizon
    double sign = cosTheta < 0.0 ? -1.0 : 1.0;
    return (sign*pow(abs(cosTheta), 1.0/texture_exponent_) + 1.0)/2.0;

    // Linear remapping
    // return cosTheta*0.5 + 0.5;
}

double Atmosphere::TextureCoordinateToCosTheta(double textureCoordinate)
{
    // Favors precision near horizon
    double remapped = textureCoordinate*2.0 - 1.0;
    double sign = remapped < 0.0 ? -1.0 : 1.0;
    return sign * pow(abs(remapped), texture_exponent_);

    // Linear remapping
    // return textureCoordinate*2.0 - 1.0;
}
