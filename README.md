# DOPpy
Read binary files (.bdd) of DOP measurement devices

The DOP measurement device series by [Signal Processing](https://www.signal-processing.com/) implements the ultrasound Doppler velocimetry (UDV) measurement method to record velocity profiles in liquids. The measurement data can be saved in an ascii file (.ADD = ascii DOP data) or a binary file (.BDD = binary DOP data). While the ascii is human readable and relatively easy to import, it is generally larger than the binary format and does only contain the time, position and velocity information of the measurement. The binary format contains much more information (measurement parameters & settings) as well as having a more compact file-size. The downside is of course that the binary structure has to be decoded. The job of DOPpy is to read the binary data of a BDD-file and provide it for further analysis within a Python console or script. 

Currently only files of the DOP2000 and DOP3000/3010 systems are supported. DOP4000 files are similar to the DOP3000 formats, but were not tested.

## Requirements
DOPpy was written in Python 3.6 but should be compatible with Python 2.7. Required modules are `numpy` and `matplotlib`.

## Usage
This section is a quick overview of the features of DOPpy. For more information refer to the doc-strings of the respective functions.

Put the file `DOPpy.py` somwhere Python can find it (e.g. the folder of your script or a folder listed in the PATH veriable). Call the function `DOP(fname)` from the module with the absolute or relative path to the BDD-file `fname`. The function determines which version (DOP2000 or 3000) the file is and returns a `DOP2000` or `DOP3000` instance, respectively.

```python
import DOP from DOPpy
bdd = DOP('path/to/file.bdd')
```

### Getting Parameters & Data
These classes have all available informations from the read file stored. They are subclasses of `DOPBase`. Every parameter of the file can be retrieved either by using the method `bdd.getParam(param)` or by treating the instance as a dictionary `bdd[param]`, where `param` is a string with the name of the parameter. All available parameters can be retrieved with the `bdd.keys()` method. Since there can be a lot of parameters, a list of all parameters that contain a search string can be retrieved by calling ``bdd.keysSearch(searchterm)``. Channel-specific parameters are named with a prefix `'ch1_'` followed by the name (e.g. `'ch2_prf'` gives the pulse repetition frequency of channel 2). The available channel-specific parameter names can be retrieved with `bdd.keysChannel(channel)`.

Additionally there are some predefined methods to retrieve often used data:

- `bdd.getChannels()` returns the list of channels that were used in the measurement.
- `bdd.getProfileType(channel)` returns the profile types that where recorded in a channel.
- `bdd.getTime(channel)` returns the timestamps of the channel in seconds.
- `bdd.getDepth(channel)` returns the gate depths of the channel in millimeter.
- `bdd.getVelocity(channel)` returns the measured velocity of the channel in meter per second as a 2D numpy-array. The first dimension is the time, the second dimension the depth.
- `bdd.getEcho(channel)` returns the echo amplitude of the channel as a 2D numpy-array. The first dimension is the time, the second dimension the depth.
- `bdd.getChannelParam(param, channel)` returns the value of the channel-specific parameter given in `param` of the channel (see `bdd.keysChannel`).

For all these methods `channel` can be an integer or a list of integers. In the latter case a list of values is returned. If `channel` is omitted or set to `None`, it defaults to all recorded channels (see `bdd.getChannels`). These methods raise an error if the data is not available for at least one requested channel.

- `bdd.printSettings(channel)` prints relevant operation parameters that were used during the measurement with the given channel. The mandatory argument `channel` is an integer.

### Displaying Measurements
While DOPpy is primarily designed to import the data into python a couple of functions for quick visualisation of the measured data are available:

- `bdd.contour(profile, channel)` plots a color-coded contour-plot of the profile with type `profile` for the specified channel over time and depth. See `bdd.getChannels` and `bdd.getProfileType` for available channels and their profile types, respectively.
- `bdd.replay(profile, channel)` plots an animation of the profile-snapshots with type `profile` for the specified channel over the depth. See `bdd.getChannels` and `bdd.getProfileType` for available channels and their profile types, respectively.

For both methods `channel` can be an integer or a list of integers. They raise an error if the data is not available for at least one requested channel.


## Notes / Known bugs
- 2D/3D-component measurements are currently not supported.
- Very old DOP2000 file versions might not be read correctly. The script was tested with file version `'BINWDOPV4.06.1'` .

## Changelog
v2.10:
- Optimized code for reading the BDD-files (for `DOP2000` and `DOP3000`), which is now less memory intensive. This allows larger files to be read faster and with less possibility of an `Memory Error`. The raw measurement data from the file are by default not saved any more. This behaviour can be changed by setting the `saveMeas` argument to `True` in the `DOP`-function.
- Transformed certain parameters from their raw value to sensible ones used by the DOP-software. The raw parameter values are stored with an suffix `'_file'` to the parameter name. Examples: emitting power, TGC mode, sampling volume (`DOP3000` only), spatial filter bandwidth, sensitivity. 
- Added the `DOPBase.printSettings` method to display important operation parameters used in the measurement.
- Corrected axis labels in `DOPBase.contour` and `DOPBase.replay` and fixed handling of `timerange` and `depthrange` arguments  in `DOPBase.contour`.
- When decoding the comment string any character that does not fit the used codec is now ignored. This behaviour can be changed with the `decode_errors` keyword-argument for the function `DOP`.
      
v2.04:
- Corrected time-calculation for `DOP2000`.

v2.03:
- Imports true division from `__future__` to avoid incompatibilities for use with Python 2.7 when dividing by integers.
- Corrected a typo in `keysChannel`.
    
v2.02:
- Fixed Doppler angle in velocity calculation of `DOP3000` and added Doppler angle in velocity calculation of `DOP2000`.
- Added `DOPBase.removeAliasing` to remove aliasing effects for smooth velocity data.
- Added the option to plot custom horizontal and vertical lines in the `DOPBase.replay` method (see the Keyword-Argument section in the method's documentation).

v2.01:
- Removed aliasing from echo-profile.
- Added `maxtimes` argument to `DOPBase.contour` to avoid RAM overflow.

v2.00:
- Initial "release"
