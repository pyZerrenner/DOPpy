# DOPpy
Read binary files (.bdd) of DOP measurement devices

The DOP measurement device series by [Signal Processing](https://www.signal-processing.com/) implements the ultrasonic Doppler velocimetry (UDV) measurement method to record velocity profiles in liquids. The measurement data can be saved in an ascii file (.ADD = ascii DOP data) or a binary file (.BDD = binary DOP data). While the ascii is human readable and relatively easy to import, it is generally larger than the binary format and does only contain the time, position and velocity information of the measurement. The binary format contains much more information (measurement parameters & settings) as well as having a more compact file-size. The downside is of course that the binary structure has to be decoded. The Job of DOPpy is to read the binary data of a BDD-file and provide it for further analysis within a Python console or scripts. 

Currently only files of the DOP2000 and DOP3000/3010 systems are supported. DOP4000 files are similar to the DOP3000 formats, but were not tested.

## Requirements
DOPpy was written in Python 3.6 but should be compatible with Python 2.7. Required modules are `numpy` and `matplotlib`.

## Usage
This section is a quick overview of the features of DOPpy. For more information refer to the doc-strings of the respective functions.

Call the function `DOP(fname)` with the absolute or relative path to the BDD-file `fname`. The function determines which version (DOP2000 or 3000) the file is and returns a `DOP2000` or `DOP3000` instance, respectively.

```python
bdd = DOP(fname)
```

### Getting Parameters & Data
These classes have all available informations from the read file stored. They are subclasses of `DOPBase`. Every parameter of the file can be retrieved either by using the method `bdd.getParam(param)` or by treating the instance as a dictionary `bdd[param]`, where `param` is a string with the name of the parameter. All available parameters can be retrieved with the `bdd.keys()` method. Since there can be a lot of parameters, a list of all parameters that contain a search string can be retrieved by calling ``DOPBase.keysSearch(searchterm)``. Channel-specific parameters are named with a prefix `'ch1_'` followed by the name (e.g. `'ch2_prf'` gives the pulse repetition frequency of channel 2). The available channel-specific parameter names can be retrieved with `bdd.keysChannel(channel)`.

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
- `bdd.replay(profile, channel)` plots an animation of the profile-snapshots with type `profile` for the specified channel over the depth. See `DOPBase.getChannels` and `DOPBase.getProfileType` for available channels and their profile types, respectively.

For both methods `channel` can be an integer or a list of integers. They raise an error if the data is not available for at least one requested channel.


## Notes / Known bugs
- 2D/3D-component measurements are currently not supported.
- Very old DOP2000 file versions might not be read correctly. The script was tested with file version `'BINWDOPV4.06.1'` .
