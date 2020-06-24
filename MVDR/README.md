# sap-offline-beamformer-design
A toolbox to design FIR beamformer coefficients based on the array geometry

_N.B. This repository depends on submodules so when cloning it you should do:_

```
git clone --recurse-submodules https://github.com/ImperialCollegeLondon/sap-offline-beamformer-design.git
```

## Usage
The main function to call is `design_beamformer`. It takes a single struct as input which specifies the required and optional parameters. An example of how to do this is provided in `example_script.m`.

Mininally, you mush specify

- sample rate
- look direction (as a cartesian direction vector)
- the array properties

The complexity of the array is wrapped up using the ElobesMicrophoneArray class. For a uniform linear array, just use the example script and change the number of sensors and their spacing.
