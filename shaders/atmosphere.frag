const vec2 RAYLEIGH_NORM = vec2(0.0, 0.05588319);
const vec2 MIE_NORM = vec2(0.0, .02527083);
const vec3 SPECTRAL_IRRADIANCE = vec3(1.526, 1.91, 2.08) / 10.0;
const vec3 SPECTRAL_TO_RGB = vec3(133.3209, 88.51855, 112.7552);
const float SUN_ANGULAR_RADIUS = 0.004675034;

uniform vec3 lightDir;

uniform sampler2D rayleighTexture;
uniform sampler2D mieTexture;

uniform float exposure;
uniform bool rayleigh;
uniform bool mie;
uniform float mieG;

varying vec3 vWorldPosition;

float RayleighPhaseFunction(float cosTheta)
{
    // Original rayleigh phase function.
    // return 0.75 * (1 + pow(cosTheta, 2));

    // Modified to better account for sun-view azimuth as described in Section 4.1 of:
    // http://publications.lib.chalmers.se/records/fulltext/203057/203057.pdf
    return 0.8 * (1.4 + 0.5*cosTheta);
}

float MiePhaseFunction(float cosTheta, float g)
{
    float g2 = g * g;
    float t2 = cosTheta * cosTheta;
    float result = 3.0 / 2.0;
    result *= (1.0 - g2) / (2.0 + g2);
    result *= (1.0 + t2) / pow(1.0 + g2 - 2.0*g*t2, 3.0/2.0);
    return result;
}

void main()
{
    vec3 viewDir = normalize(vWorldPosition - cameraPosition);

    // Calculate the view-zenith and sun-zenith angles.
    float cosV = dot(viewDir, vec3(0, -1, 0));
    float cosL = dot(lightDir, vec3(0, 1, 0));

    // Convert the angles to texture coordinates using the parameterization function.
    // Note: we use abs+sign to avoid negative roots!
    float u = 0.5 * (1.0 + sign(cosV)*pow(abs(cosV), 1.0/3.0));
    float v = 0.5 * (1.0 + sign(cosL)*pow(abs(cosL), 1.0/3.0));

    // Sample the textures.
    vec3 rayleighVal = texture2D(rayleighTexture, vec2(u, v)).rgb;
    vec3 mieVal = texture2D(mieTexture, vec2(u, v)).rgb;

    // Remap the values.
    rayleighVal = rayleighVal*(RAYLEIGH_NORM.y-RAYLEIGH_NORM.x) + RAYLEIGH_NORM.x;
    mieVal = mieVal*(MIE_NORM.y-MIE_NORM.x) + MIE_NORM.x;

    // Calculate the view-sun angle for the phase function.
    // Note: we clamp it between [0, 1] or else we would get the sun
    // on both sides of the light direction.
    float cosTheta = dot(viewDir, lightDir);
    cosTheta = saturate(cosTheta);

    // Apply the phase function.
    rayleighVal *= RayleighPhaseFunction(cosTheta);
    mieVal *= MiePhaseFunction(cosTheta, mieG);

    // Compute the scattering, and apply the spectral intensity to
    // get the spectral radiance for this fragment.
    vec3 radiance = vec3(0.0);
    if (rayleigh)
    {
        radiance += rayleighVal;
    }
    if (mie)
    {
        radiance += mieVal;
    }
    radiance *= SPECTRAL_IRRADIANCE * vec3(exposure);

    // Multiply by the SPECTRAL_TO_RGB conversion constants to convert
    // the spectral radiance to RGB values. These values must be
    // supplied from the precomputation process.
    vec3 rgb = radiance * SPECTRAL_TO_RGB;

    if (acos(cosTheta) < SUN_ANGULAR_RADIUS)
    {
        // TODO: this is not physically correct... Only works for intensity <= 1.
        // Looks nice though!
        rgb /= SPECTRAL_IRRADIANCE * vec3(exposure);
    }

    // Tonemap the resulting RGB samples into valid RGB ranges.
    rgb = pow(vec3(1.0) - exp(-rgb), vec3(1.0/2.2));

    gl_FragColor = vec4(rgb, 1.0);
}
