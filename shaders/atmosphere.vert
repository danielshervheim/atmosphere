varying vec3 vWorldPosition;

void main()
{
    // Center the skysphere over the camera position.
    vec4 posWS = vec4(cameraPosition + position, 1.0);
    vec4 posVS = viewMatrix * posWS;
    vec4 posCS = projectionMatrix * posVS;

    vWorldPosition = posWS.xyz;

    // Set z to w, in other words as far as possible.
    // This makes the skysphere render behind everything else.
    gl_Position = posCS.xyww;
}
