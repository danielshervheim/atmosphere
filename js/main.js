// Daniel Shervheim, 2019
// danielshervheim.com

// Constants.
const TEXTURE_DIMENSION = 64;

// Content to be loaded.
var rayleighTexture;
var mieTexture;
var vertexShader;
var fragmentShader;

// Scene contents.
var renderer;
var scene;

// The main timer.
var clock;
var dt;

// Mouse and touch control variables.
var mousePosition = new THREE.Vector2(0.0, 0.0);
var mouseDown = false;
var rotationTarget = 'light';

// The main camera.
var camera;
var cameraPitch = 0.32;
var cameraYaw = 0.0;

// The main directional light.
var light;
var lightPitch = -0.48;
var lightYaw = -2.9;

// Sky sphere.
var skySphere;

// Setup a loading manager.
var manager = new THREE.LoadingManager();
manager.onProgress = function(url, itemsLoaded, itemsTotal)
{
    console.log( 'Loading file: ' + url + '.\nLoaded ' + itemsLoaded + ' of ' + itemsTotal + ' files.' );
}

manager.onLoad = function()
{
    init();
}

// Setup a loader.
THREE.Cache.enabled = true;
var loader = new THREE.FileLoader(manager);

// Load in all the files we need.
loader.setResponseType('arraybuffer');
loader.load(
    'textures/rayleigh.bin',
    function (response)
    {
        rayleighTexture = new THREE.DataTexture(new Float32Array(response),
            TEXTURE_DIMENSION, TEXTURE_DIMENSION, THREE.RGBFormat
        );
        rayleighTexture.type = THREE.FloatType;
        rayleighTexture.minFilter = THREE.LinearFilter;
        rayleighTexture.magFilter = THREE.LinearFilter;
        rayleighTexture.encoding = THREE.LinearEncoding;
        rayleighTexture.flipY = true;
        rayleighTexture.needsUpdate = true;
    }
);
loader.load(
    'textures/mie.bin',
    function (response)
    {
        mieTexture = new THREE.DataTexture(new Float32Array(response),
            TEXTURE_DIMENSION, TEXTURE_DIMENSION, THREE.RGBFormat
        );
        mieTexture.type = THREE.FloatType;
        mieTexture.minFilter = THREE.LinearFilter;
        mieTexture.magFilter = THREE.LinearFilter;
        mieTexture.encoding = THREE.LinearEncoding;
        mieTexture.flipY = true;
        mieTexture.needsUpdate = true;
    }
);
loader.setResponseType('text');
loader.load(
    'shaders/atmosphere.vert',
    function (response)
    {
        vertexShader = response;
    }
);
loader.load(
    'shaders/atmosphere.frag',
    function (response)
    {
        fragmentShader = response;
    }
);

function init()
{
    // Create the clock.
    clock = new THREE.Clock();

    // Setup renderer.
    var canvas = document.querySelector('#c');
    renderer = new THREE.WebGLRenderer({canvas});
    renderer.getContext().getExtension('OES_texture_float');
    renderer.getContext().getExtension('OES_texture_float_linear');

    // Create the scene.
    scene = new THREE.Scene();

    // Setup camera and light.
    camera = new THREE.PerspectiveCamera(60.0, 2.0, 0.1, 100.0);
    scene.add(camera);
    light = new THREE.DirectionalLight(0xFFFFFF, 1);
    scene.add(light);

    // Setup the skysphere.
    var geo = new THREE.SphereGeometry(1, 32, 24);
    var mat = new THREE.ShaderMaterial({
        uniforms:
        {
            lightDir: { value: new THREE.Vector3(0, 0, 1).applyQuaternion(light.quaternion) },
            rayleighTexture: { value: rayleighTexture },
            mieTexture: { value: mieTexture },
            exposure: { value: 1.0 },
            rayleighEnabled: { value: true },
            mieEnabled: { value: true },
            mieG: { value: 0.85 },
        },
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
    });
    mat.side = THREE.BackSide;
    skySphere = new THREE.Mesh(geo, mat);
    skySphere.frustumCulled = false;
    scene.add(skySphere);

    // Create the GUI.
    var guiOptions = {
        exposure : 1.0,
        rayleigh : true,
        mie : true,
        mieG : 0.85,
        rotate: rotationTarget,
    };
    var gui = new dat.gui.GUI();
    gui.remember(guiOptions);
    gui.add(guiOptions, 'rayleigh').name("Rayleigh Scattering").onChange( function ()
    {
        skySphere.material.uniforms.rayleighEnabled.value = guiOptions.rayleigh;
    });
    gui.add(guiOptions, 'mie').name("Mie Scattering").onChange( function ()
    {
        skySphere.material.uniforms.mieEnabled.value = guiOptions.mie;
    });
    gui.add(guiOptions, 'mieG').min(-1.0).max(1.0).name("Mie Asymmetry Parameter").onChange( function ()
    {
        skySphere.material.uniforms.mieG.value = guiOptions.mieG;
    });
    gui.add(guiOptions, 'exposure').name("Exposure").onChange( function ()
    {
        skySphere.material.uniforms.exposure.value = guiOptions.exposure;
    });
    gui.add(guiOptions, 'rotate', [ 'light', 'camera' ]).name("Rotate").onChange( function ()
    {
        rotationTarget = guiOptions.rotate;
    });

    // Install mouse event handlers.
    canvas.addEventListener('mousedown', function (e)
    {
        mousePosition = new THREE.Vector2(e.offsetX, e.offsetY);
        mouseDown = (e.button == 0);
    });
    canvas.addEventListener('mousemove', function (e)
    {
        if (mouseDown)
        {
            var deltaPosition = new THREE.Vector2(e.offsetX - mousePosition.x, e.offsetY - mousePosition.y);
            mousePosition = new THREE.Vector2(e.offsetX, e.offsetY);

            if (rotationTarget == 'camera')
            {
                rotateCamera(deltaPosition, 5.0);
            }
            if (rotationTarget == 'light')
            {
                rotateLight(deltaPosition, 5.0);
            }
        }
    });
    document.addEventListener('mouseup', function (e)
    {
        mouseDown = false;
    });

    // Call cam and light update functions once to set initial position.
    rotateCamera(null, 0.0);
    rotateLight(null, 0.0);

    // Start the first update call.
    requestAnimationFrame(update);
}

function update()
{
    // Resize the canvas and renderer if the window changed.
    if (resizeRendererToDisplaySize(renderer))
    {
        const canvas = renderer.domElement;
        camera.aspect = canvas.clientWidth / canvas.clientHeight;
        camera.updateProjectionMatrix();
    }

    dt = clock.getDelta();

    renderer.render(scene, camera);
    requestAnimationFrame(update);
}

function resizeRendererToDisplaySize(renderer)
{
    const canvas = renderer.domElement;
    const pixelRatio = window.devicePixelRatio;
    const width  = canvas.clientWidth  * pixelRatio | 0;
    const height = canvas.clientHeight * pixelRatio | 0;
    if (canvas.width != width || canvas.height != height)
    {
        renderer.setSize(width, height, false);
        return true;
    }
    return false;
}

function rotateCamera(deltaPosition, speed)
{
    if (deltaPosition != null)
    {
        cameraPitch += THREE.Math.degToRad(deltaPosition.y) * dt * speed;
        cameraPitch = THREE.Math.clamp(cameraPitch, -Math.PI/2.0 + 0.01, Math.PI/2.0 - 0.01);

        cameraYaw += THREE.Math.degToRad(deltaPosition.x) * dt * speed;
        cameraYaw %= Math.PI * 2.0;
    }

    var lookAtPos = sphericalToCartesian(cameraPitch - Math.PI/2.0, Math.PI/2.0 - cameraYaw, 1.0);
    camera.lookAt(lookAtPos.add(camera.position));
}

function rotateLight(deltaPosition, speed)
{
    if (deltaPosition != null)
    {
        lightPitch += THREE.Math.degToRad(deltaPosition.y) * dt * speed;
        lightPitch = THREE.Math.clamp(lightPitch, -Math.PI/2.0 + 0.01, Math.PI/2.0 - 0.01);

        lightYaw -= THREE.Math.degToRad(deltaPosition.x) * dt * speed;
        lightYaw %= Math.PI*2.0;
    }

    var lookAtPos = sphericalToCartesian(lightPitch - Math.PI/2.0, Math.PI/2.0 - lightYaw, 1.0);
    light.lookAt(lookAtPos.add(light.position));

    var lightForward = new THREE.Vector3(0, 0, 1).applyQuaternion(light.quaternion);
    skySphere.material.uniforms.lightDir.value = lightForward;
}

function sphericalToCartesian(theta, gamma, radius)
{
    return new THREE.Vector3(
        radius*Math.sin(theta) * Math.cos(gamma),
        radius*Math.cos(theta),
        radius*Math.sin(theta) * Math.sin(gamma)
    );
}
