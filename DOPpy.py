# -*- coding: utf-8 -*-
"""Classes for reading BDD-files of DOP 2000 and 3000/3010

| Version: 2.02
| Date: 2016-12-12

Usage
=====
    Call the function ``DOP(fname)`` with the absolute or relative path to the
    BDD-file `fname`. The function determines which version (DOP2000 or 3000)
    the file is and returns a `DOP2000` or `DOP3000` instance,
    respectively.

Getting Parameters
==================
    These classes have all available informations from the read file stored.
    They are subclasses of `DOPBase`. Every parameter of the file can be
    retrieved either by using the method ``DOPBase.getParam(param)`` or by
    treating the instance as a dictionary ``DOPBase[param]``, where `param` is
    the name of the parameter. All available parameters can be retrieved with
    the ``DOPBase.keys()`` method. Since there can be a lot of parameters, a
    list of all parameters that contain a search string can be retrieved by
    calling ``DOPBase.keysSearch(searchterm)``. Channel-specific parameters are
    named with a prefix ``'ch1_'`` followed by the name (e.g. ``'ch2_prf'``).
    The available channel-specific parameters can be retrieved with
    ``DOPBase.keysChannel(channel)``.

    Additionally there are some predefined methods to retrieve often used data:

    * ``DOPBase.getChannels()`` returns the list of channels that were used in
      the measurement.
    * ``DOPBase.getProfileType(channel)`` returns the profile types that where
      recorded in a channel.
    * ``DOPBase.getTime(channel)`` returns the timestamps of the channel in s.
    * ``DOPBase.getDepth(channel)`` returns the gate depths of the channel in
      mm.
    * ``DOPBase.getVelocity(channel)`` returns the measured velocity of the
      channel in m/s as a 2D-array of the form ``velocity[time, depth]``.
    * ``DOPBase.getEcho(channel)`` returns the echo amplitude of the channel as
      a 2D-array of the form ``echo[time, depth]``.
    * ``DOPBase.getChannelParam(param, channel)`` returns the value of the
      channel-specific parameter given in `param` of the channel (see
      `DOPBase.keysChannel`).

    For all methods `channel` can be an integer or a list of integers. In the
    latter case a list of values is returned. If `channel` is omitted or set to
    ``None``, it defaults to all recorded channels (see `DOPBase.getChannels`).
    These methods raise an error if the data is not available for at least one
    requested channel.

Displaying Measurements
=======================
    For a quick visualisation of the measured data the following methods are
    available:

    * ``DOPBase.contour(profile, channel)`` plots a color-coded contour-plot of
      the profile with type `profile` for the specified channel over time and
      depth. See `DOPBase.getChannels` and `DOPBase.getProfileType` for
      available channels and their profile types, respectively.
    * ``DOPBase.replay(profile, channel)`` plots an animation of the
      profile-snapshots with type `profile` for the specified channel over the
      depth. See `DOPBase.getChannels` and `DOPBase.getProfileType` for
      available channels and their profile types, respectively.

    For both methods `channel` can be an integer or a list of integers. They
    raise an error if the data is not available for at least one requested
    channel.


Notes
=====
    * 2D/3D-component measurements can be read but are not supported.
    * Very old DOP 2000 file versions might not be read correctly. The script
      was tested with file version 'BINWDOPV4.06.1'.


Changelog
=========
v2.00:
    * Initial release
v2.01:
    * Removed aliasing from echo-profile.
    * Added `maxtimes` argument to ``DOPBase.contour`` to avoid RAM overflow.
v2.02:
    * Fixed Doppler angle in velocity calculation of `DOP3000` and added
      Doppler angle in velocity calculation of `DOP2000`.
    * Added ``DOPBase.removeAliasing`` to remove aliasing effects for smooth
      velocity data.
    * Added the option to plot custom horizontal and vertical lines in the
      ``DOPBase.replay`` method (see the Keyword-Argument section in the
      method's documentation).
"""



import struct
from warnings import warn
import numpy as np
import bz2
import gzip
import matplotlib.pyplot as plt



class DOPBase(object):
    """ Base class for DOP measurements

    Subclasses need to implement the `_read` and `_refine` methods.

    Custom format keys
    ==================
    For unpacking the binary data of the BDD-file, `_readParam` accepts
    additional fomrat keys:

    v: Verbose (1 byte)
        Values are not unpacked. The format cannot be combined with other
        formats from the struct module.
    m: Machine code (1 bit)
        Read exactly one bit. A preceeding integer is interpreted as the
        buffer size in bytes. A following integer is interpreted as offset
        in bits inside the buffer. Example: ``'4m10'`` reads 4 bytes and
        extracts the 10th bit value. the value is returned as boolean.
        The format cannot be combined with other formats from the struct
        module.
    """
    _codec = 'cp1252'  # file codec

    def __init__(self, fname):
        """ Read a DOP binary file (*.BDD)

        The file may be a *.gz or *.bz2-archive.

        Arguments:
        ==========
        fname: str
            Path to the BDD-file.
        """
        self._fname = fname
        self._file = None
        self._values = {}
        self._chToMeas = [[] for i in range(10)]

        if self._fname.endswith('.bz2'):
            self._file = bz2.BZ2File(self._fname, 'rb')
        elif self._fname.endswith('.gz'):
            self._file = gzip.GzipFile(self._fname, 'rb')
        else:
            self._file = open(self._fname, 'rb')

        self._read()
        self._file.close()

        self._refine()


    def _read(self):
        """ Read the data in the given BDD file
        """
        # This method is implemented by subclasses.
        # Use methods `_readParam` to get a parameter value from the file
        # Use methods `setParam` and `getParam` to modify parameter values.
        pass


    def _refine(self):
        """ Process data read from the BDD file
        """
        # This method is implemented by subclasses.
        # Use methods `setParam` and `getParam` to modify parameter values.

        # This method should define the parameters:
        #   'channelUsed': List of channel numbers that were used

        # This method should define these channel-parameters for each used
        # channel (i.e. these have an additional prefix as defined in
        # `DOPBase._prefixChannel`, e.g. 'ch1_profTypeMax'):
        #   'profTypeName': List of profile names that were recorded
        #   'veloMax': Maximum velocity in m/s, velocity range is +/- veloMax
        #   'time': Time array in seconds
        #   'depth': Gate depth array in millimeter
        #   'velo'/'echo'/...: All recorded profiles as named in 'profTypeName'
        pass


    def keys(self):
        """ Returns list of available parameters
        """
        return self._values.keys()


    def keysSearch(self, searchterm):
        """ Returns list of parameters that contain the search term.

        Arguments:
        ==========
        searchterm: str
            String that is searched for in all available keys.
        """
        res = []
        keys = self.keys()
        for k in keys:
            if searchterm in k:
                res.append(k)
        return res


    def keysChannel(self, channel=None):
        """ Returns list of channel-specific parameters

        The channel-specific parameters can be accessed by the parameter name
        ``'ch<n>_<param>'``, where `<n>` is the channel number and `<param>` is
        the channel-specific parameter name (e.g. ``'ch3_prf'`` gives the pulse
        repetition frequency of channel 3).

        Attention: If multiple channels are requested the returned parameters
        may not be present for all channels.

        Arguments:
        ==========
        channel: int, list or None
            Channel number (1-10) or list of channel numbers whose available
            parameters should be returned. If None `channel` defaults to all
            used channels.
        """
        if channel is None:
            channel = self.getChannels()
        else:
            if isinstance(channel, int):
                channel = [channel]

        res = []
        for ch in channel:
            preCh = self._prefixChannel(ch)
            keys = self.keySearch(preCh)
            for k in keys:
                k = k[len(preCh):]
                if k not in res:
                    res.append(k)

        return res


    def __contains__(self, key):
        """ Returns whether parameter is available
        """
        return self._values.__contains__(key)


    def _prefixChannel(self, channel):
        """ Returns parameter name prefix for a channel """
        return 'ch{:d}_'.format(channel)


    def _prefixMeas(self, meas):
        """ Returns parameter name prefix for a measurement """
        return 'meas{:d}_'.format(meas)


    def _prefixProfile(self, profile):
        """ Returns parameter name prefix for a profile """
        return 'prof{:d}_'.format(profile)


    def _readParam(self, param, offset, fmt):
        """ Read and return a single parameter from file

        Call only if self._file is an opened file.

        Arguments:
        ==========
        param: str
            Parameter name
        offset:
            Absolute offset from the start of the file.
        fmt: str
            Data-format string of the parameter (see help(struct) and
            help(DOP3000))

        Returns:
        ========
        value:
            The value read from the file. Its type is determined by `fmt`.
        """
        self._file.seek(offset)

        if 'm' not in fmt and 'v' not in fmt:
            # normal struct format
            size = struct.calcsize(fmt)
            value = struct.unpack(fmt, self._file.read(size))
            if len(value) == 1:
                value = value[0]
        elif 'v' in fmt:
            # special verbose format (don't unpack)
            size = int(fmt[:-1])
            value = self._file.read(size)
        elif 'm' in fmt:
            # special machine code format (read a single bit)
            ind = fmt.index('m')
            # buffer size in bytes
            if ind != 0:
                size = int(fmt[:ind])
            else:
                size = 1
            # bit position inside buffer
            if ind != len(fmt)-1:
                pos = int(fmt[ind+1:])
            else:
                pos = 0

            value = self._byteToBit(self._file.read(size))[pos]
            value = bool(int(value))

        self.setParam(param, value)

        return value


    def setParam(self, param, value):
        """ Set the value of a parameter
        """
        self._values[param] = value


    def getParam(self, param):
        """ Returns the value of a parameter
        """
        return self._values[param]


    def __getitem__(self, param):
        """ Return the value of a parameter
        """
        return self.getParam(param)


    def removeAliasing(self, jumpSize=.8):
        """ Removes aliasing effects by detecting large jumps in the velocity

        This method assumes that the first gate has a velocity that does not
        show an aliasing effect.

        Do not use this method with very noisy data (e.g. when the DOP's
        sensitivity is set to 'Very high'). The algorithm cannot distinguish
        between noise spikes and aliasing jumps. Compare the results with the
        ``DOPBase.replay`` method to ensure that the data is correct.

        This change in the velocity profile is currently not reversible. The
        BDD-file has to be re-read if the original profile is to be recovered.

        Arguments:
        ==========
        jumpSize: float
            Gives the threshold, at which a jump in velocity is considered an
            aliasing effect. The value is in units of the whole velocity range:
            if the velocity range is e.g. +/- 40 mm/s, then any difference of
            neighbouring velocities above ``(jumpSize*80)`` mm/s is considered
            an aliasing effect. The maximum velocity for all channels can be
            retrieved by calling the method ``getChannelParam('veloMax')``.
        """
        for ch in self.getChannels():
            preCh = self._prefixChannel(ch)
            vmax = self.getParam(preCh+'veloMax')
            jump = 2*vmax*jumpSize

            velo = self.getParam(preCh+'velo').copy()
            veloDiff = np.diff(velo, axis=1)

            posJump = zip(*np.where(veloDiff > jump))
            negJump = zip(*np.where(veloDiff < -jump))

            for ti, di in posJump:
                velo[ti,di+1:] -= 2*vmax
            for ti, di in negJump:
                velo[ti,di+1:] += 2*vmax

            self.setParam(preCh+'velo', velo)

    def getChannels(self):
        """ Returns the channels in use
        """
        return self.getParam('channelUsed')


    def getChannelParam(self, param, channel=None):
        """ Returns the parameter value for given channels

        Arguments:
        ==========
        param: str
            Parameter name. Must be a channel-parameter.
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.
        """
        if channel is None:
            channel = self.getChannels()

        try:
            value = []
            for ch in channel:
                preCh = self._prefixChannel(ch)
                value.append(self.getParam(preCh+param))
        except TypeError:
            preCh = self._prefixChannel(channel)
            value = self.getParam(preCh+param)

        return value


    def getProfileType(self, channel=None):
        """ Returns the available profile types for given channels

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.

        Returns:
        ========
        profileTypes: list
            If `channel` is an integer a list of strings is returned. If
            `channel` is a list of ints, a list is returned with each element
            corresponding to the same element given in `channel`.
        """
        return self.getChannelParam('profTypeName', channel)


    def getTime(self, channel=None):
        """ Returns the timestamps for given channels in seconds

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.

        Returns:
        ========
        time: array or list
            If `channel` is an integer an array is returned. If `channel` is a
            list of ints, a list of arrays is returned with each element
            corresponding to the same element given in `channel`.
        """
        return self.getChannelParam('time', channel)


    def getDepth(self, channel=None):
        """ Returns the gate depths for given channels in mm

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.

        Returns:
        ========
        depth: array or list
            If `channel` is an integer an array is returned. If `channel` is a
            list of ints, a list of arrays is returned with each element
            corresponding to the same element given in `channel`.
        """
        return self.getChannelParam('depth', channel)


    def getVelocity(self, channel=None):
        """ Returns the velocity field for given channels in m/s

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.

        Returns:
        ========
        velo: array or list
            If `channel` is an integer a 2d-array of the form
            ``velo[time, depth]`` is returned. If `channel` is a list of ints,
            a list of 2d-arrays is returned with each element corresponding to
            the same element given in `channel`.
        """
        return self.getChannelParam('velo', channel)


    def getEcho(self, channel=None):
        """ Returns the echo field for given channels

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers. If
            ``None`` is given all available channels are used.

        Returns:
        ========
        echo: array or list
            If `channel` is an integer a 2d-array of the form
            ``echo[time, depth]`` is returned. If `channel` is a list of ints,
            a list of 2d-arrays is returned with each element corresponding to
            the same element given in `channel`.
        """
        return self.getChannelParam('echo', channel)


    def contour(self, profile, channel=None, timerange=slice(None),
                depthrange=slice(None), maxtimes=1000, **kw):
        """ Show the contour plot of a profile over time and depth

        Attention:
            Use this method only as a quick way to visualise the UDV-data.
            To save RAM, a maximum of only 1000 timesteps are selected by
            default in the time range given. This possible exclusion of
            intermediate timesteps may create errors in the visualisation.
            See the explanation of `maxtimes` for more information.

        Arguments:
        ==========
        profile: str
            Profile type to be plotted.
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers to be
            plotted. If ``None``, all available cannels are taken.
        timerange: slice or list
            Slice or indices of the time-array to be plotted. ``slice(None)``
            uses the whole available time range.
        depthrange: slice or list
            Slice or indices of the depth-array to be plotted. ``slice(None)``
            uses the whole available depth range.
        maxtimes: int
            No more than `maxtimes` timesteps are used for plotting. This can
            be used if the length of the measurement is very large, which could
            cause the RAM to be overloaded by the ``pyplot.contourf`` function.
            Any value below 1 plots ALL timesteps (use at own risk). Be aware
            that this option can discard intermediate timesteps from plotting
            completely, which could cause artefacts or errors in the
            visualisation.

        Keyword-Arguments:
        ==================
        levelN: int or None
            Number of levels for the contour plot. Is ignored if set to None.
        All other keyword arguments are passed to the pyplot.contourf function.
        """
        levelN = kw.pop('levelN', None) # number of levels

        if channel is None:
            channel = self.getChannels()

        # number of channels
        try:
            chN = len(channel)
        except:
            chN = 1
            channel = [channel]

        time = self.getTime(channel)
        depth = self.getDepth(channel)
        data = self.getChannelParam(profile, channel)

        fig, ax = plt.subplots(chN, 1, squeeze=False, sharex=True)
        ax[0,-1].set_xlabel('Time [s]')

        for ci, ch in enumerate(channel):
            ax[ci, 0].set_ylabel('Depth [mm]')

            X, Y = np.meshgrid(time[ci][timerange], depth[ci][depthrange])
            Z = data[ci][timerange, depthrange].T

            tN = X.shape[1]
            if 1 <= maxtimes < tN:
                inds = np.round(np.linspace(0, tN-1, maxtimes)).astype(int)
                X = X[:,inds]
                Y = Y[:,inds]
                Z = Z[:,inds]


            if levelN is None:
                cont = ax[ci, 0].contourf(X, Y, Z, **kw)
            else:
                cont = ax[ci, 0].contourf(X, Y, Z, levelN, **kw)
            plt.colorbar(cont, ax=ax[ci, 0])


    def replay(self, profile, channel, start=0, end=-1, fps=None,
               showMean=False, showRunMean=False, **kw):
        """ Show an animated timeline of a profile

        Arguments:
        ==========
        profile: str
            Profile type to be shown. See self.getProfileType for available
            options.
        channel: int or list
            Channel number (1-10) or list of channel numbers to be plotted.
            The channels given must contain the profiles of type `profile`.
        start: int
            Starting time index. Can be negative.
        end: int
            Ending time index. If less than `start` the animation lasts until
            the last timestep.
        fps: float
            Speed of the animation in frames per second. If None the average
            fps of the first given channel is taken (i.e. a somewhat realistic
            speed).
        showMean: bool
            Plot the time-average over the timeline (from `start` to `end`).
        showRunMean: bool
            Plot the running time-average over the timeline (from `start` to
            the current frame)

        Keyword-Arguments:
        ==================
        animStyle: dict
            Plotting style dictionary of the snapshot lines. These are passed
            to pyplot.plot as keyword-arguments.
        meanStyle: dict
            Plotting style dictionary of the mean lines. These are passed
            to pyplot.plot as keyword-arguments.
        runMeanStyle: dict
            Plotting style dictionary of the running mean lines. These are
            passed to pyplot.plot as keyword-arguments.
        hlines: list or tuple
            List of positions at which horizontal lines are plotted. By default
            the lines are plotted over the whole width of the plot. If multiple
            channels are given in `channel`, a tuple of lists for each
            respective channel may be given. If a single list is given, its
            values are plotted for all channels.
        hlinesStyle: dict
            Keyword dictionary passed to ``matplotlib.pyplot.hlines`` alongside
            the positions given in the `hlines` argument. May include `'xmin'`-
            and `'xmax'`-entries to define the lines' start and end. If
            multiple channels are given in `channel`, a tuple of dicts for each
            respective channel may be given. If a single dict is given, its
            entries are used for all channels.
        vlines: list or tuple
            List of positions at which vertical lines are plotted. By default
            the lines are plotted over the whole height of the plot. If
            multiple channels are given in `channel`, a tuple of lists for each
            respective channel may be given. If a single list is given, its
            values are plotted for all channels.
        vlinesStyle: dict
            Keyword dictionary passed to ``matplotlib.pyplot.hlines`` alongside
            the positions given in the `hlines` argument. May include `'ymin'`-
            and `'ymax'`-entries to define the lines' start and end. If
            multiple channels are given in `channel`, a tuple of dicts for each
            respective channel may be given. If a single dict is given, its
            entries are used for all channels.
        """
        def getData(channel, index='all'):
            data = self.getChannelParam(profile, channel)
            if index != 'all':
                data = data[index, :]
            return data

        # number of channels
        try:
            chN = len(channel)
        except:
            chN = 1
            channel = [channel]

        time = self.getTime(channel[0])
        if end <= start:
            end = len(time)

        if fps is None:
            fps = 1./np.mean(np.ediff1d(time))

        animStyle = {'color': 'b', 'linestyle': '-', 'marker': ''}
        animStyle.update(kw.pop('animStyle', {}))
        meanStyle = {'color': 'r', 'linestyle': '-', 'marker': ''}
        meanStyle.update(kw.pop('meanStyle', {}))
        runMeanStyle = {'color': 'g', 'linestyle': '-', 'marker': ''}
        runMeanStyle.update(kw.pop('runMeanStyle', {}))

        hlines = kw.pop('hlines', [])
        vlines = kw.pop('vlines', [])

        fig, ax = plt.subplots(chN, 1, squeeze=False)

        lines = [None]*chN
        meanLines = [None]*chN
        runMeanLines = [None]*chN
        titleText = 'Profile {:d} of ' + '{:d}\n'.format(end-start) + \
                    'Timestamp: {:.2f} s'
        for ci, ch in enumerate(channel):
            a = ax[ci,0]
            a.grid()
            a.set_ylim(np.min(getData(ch, slice(start,end))),
                       np.max(getData(ch, slice(start,end))))

            lines[ci], = a.plot(self.getDepth(ch), getData(ch, start),
                                **animStyle)

            if showRunMean:
                runMeanLines[ci], = a.plot(self.getDepth(ch),
                                           getData(ch, start),
                                           **runMeanStyle)
            if showMean:
                data = np.mean(getData(ch, slice(start,end)), axis=0)
                meanLines[ci], = a.plot(self.getDepth(ch), data,
                                        **meanStyle)

            if len(hlines) > 0:
                if isinstance(hlines[0], (list, tuple)):
                    yvals = hlines[0]
                else:
                    yvals = hlines

                hlinesStyle = dict(zip(['xmin', 'xmax'], a.get_xlim()))
                hlinesStyle.update(kw.pop('hlinesStyle', {}))

                a.hlines(yvals, **hlinesStyle)

            if len(vlines) > 0:
                if isinstance(vlines[0], (list, tuple)):
                    xvals = vlines[0]
                else:
                    xvals = vlines

                vlinesStyle = dict(zip(['ymin', 'ymax'], a.get_ylim()))
                vlinesStyle.update(kw.pop('vlinesStyle', {}))

                a.vlines(xvals, **vlinesStyle)


        for i, ti in enumerate(range(start, end)):
            plt.pause(1./fps)
            for ci, ch in enumerate(channel):
                if not plt.fignum_exists(fig.number):
                    return

                lines[ci].set_ydata(getData(ch, ti))
                ax[0,0].set_title(titleText.format(i+1, time[ti]))
                if showRunMean:
                    data = np.array(runMeanLines[ci].get_ydata())
                    data = (getData(ch, ti) + data*i) / (i+1)
                    runMeanLines[ci].set_ydata(data)



class DOP2000(DOPBase):
    """ Class to read DOP2000 files

    Based on the DOP2000 User's manual
    """

    ### Parameters with fixed offsets
    _fixedParam = [
        # information block at offset 0
        ['version', 0+0, '16s'],
        ['parameter', 0+16, '1024s'],
        ['comment', 0+1040, '476s'],

        # hardware information block at offset 1516
        ['hardware', 1516+0, '20B'], # exact decoding unknown

        # operation parameter block at offset 1536
        ['mainFreq', 1536+0, 'I'],
        ['prf', 1536+4, 'I'],  # in us
        ['emitFreqOff', 1536+8, 'I'],  # emitFreq is at offset 279+emitFreqOff
        ['burstLengh', 1536+12, 'I'],  # in cycles
        ['emitPower', 1536+16, 'I'],  # 0,1,2=low,medium,high
        ['gateN', 1536+20, 'I'],
        ['emitN', 1536+24, 'I'],
        ['_internalUse_2', 1536+28, 'I'],  # exact decoding unknown
        ['_internalUse_3', 1536+32, 'I'],  # exact decoding unknown
        ['soundSpeed', 1536+36, 'I'],  # in m/s
        ['veloScale', 1536+40, 'I'],
        ['timeRes', 1536+44, 'I'],  # in ns
        ['timeGate1', 1536+48, 'I'],  # in ns
        ['dopplerAngle', 1536+52, 'i'],  # in degree, 'I' also possible
        ['unit', 1536+56, 'I'],  # use doppler angle
        ['veloOffset', 1536+60, 'i'],
        ['tgcStart', 1536+64, 'I'],  # in dB
        ['tgcEnd', 1536+68, 'I'],  # in dB
        ['gateAutoN', 1536+72, '?xxx'],
        ['gateAutoNmax', 1536+76, 'I'],
        ['memoryProfileN', 1536+80, 'I'],
        ['skipProfileN', 1536+84, 'I'],  # exact decoding unknown
        ['fftWindow', 1536+88, 'I'],  # exact decoding unknown
        ['gateFirst', 1536+92, 'I'],
        ['gateLast', 1536+96, 'I'],
        ['gateInfoAll', 1536+100, '?xxx'],
        ['profType', 1536+104, 'I'],  # 0=velo, 1=echo,
            # 2=energy, 3=gate FFT, 10=velo+echo, 11=velo+energy,
            # 12=velo+v(t), 13=velo+flow rate, 14=velo+history,
            # 15=echo+history, 16=energy+history
        ['fftPointN', 1536+108, 'I'],
        ['_internalUse_4', 1536+112, 'I'],  # exact decoding unknown
        ['_internalUse_5', 1536+116, 'I'],  # exact decoding unknown
        ['moduleScale', 1536+120, 'I'],
        ['filter', 1536+124, 'Bxxx'],  # 0=none, 1=moveAvg, 2=median
        ['filterMAvgRejectZero', 1536+128, '?xxx'],
        ['filterMAvgProfileN', 1536+132, 'I'],
        ['filterMedianProfileN', 1536+136, 'I'],
        ['fileFormat', 1536+140, 'I'],  # 0=binary, 1=ASCII, 2=both
        ['fileInclComment', 1536+144, '?xxx'],
        ['fileInclDepth', 1536+148, '?xxx'],
        ['statistics', 1536+152, '?xxx'],
        ['triggerMode', 1536+156, 'I'],  # 0=no, 1=button, 2=BNC, 3=BNC low
        ['triggerDisplayWait', 1536+160, '?xxx'],
        ['triggerProfileN', 1536+164, 'I'],  # by sequence
        ['triggerSequenceN', 1536+168, 'I'],
        ['triggerDelay', 1536+172, 'I'],  # in us, exact decoding unknown
        ['triggerAutoRearm', 1536+176, '?xxx'],
        ['triggerAutoRecord', 1536+180, '?xxx'],
        ['multi', 1536+184, '?xxx'],
        ['_internalUse_6', 1536+188, 'I'],  # exact decoding unknown
        ['paramRerecordStart', 1536+192, 'I'],
        ['paramRerecordEnd', 1536+196, 'I'],
        ['multiSequenceN', 1536+200, 'I'],
        ['_internalUse_7', 1536+204, 'I'],  # exact decoding unknown
        ['_internalUse_8', 1536+208, 'I'],  # exact decoding unknown
        ['_internalUse_9', 1536+212, 'I'],  # exact decoding unknown
        ['_internalUse_10', 1536+216, 'I'],  # exact decoding unknown
        ['multiInitInTrigger', 1536+220, '?xxx'],
        ['iqGateN', 1536+224, 'I'],
        ['iqEmitN', 1536+228, 'I'],
        ['graphLCursor1Index', 1536+232, 'I'],
        ['graphLCursor2Index', 1536+236, 'I'],
        ['graphRCursor1Index', 1536+240, 'I'],
        ['graphTimeScaleVelo', 1536+244, 'I'],
        ['graphTimeScaleFlow', 1536+248, 'I'],
        ['flowScale', 1536+252, 'I'],
        ['flowUnit', 1536+256, 'I'],  # 0=ml/min, 1=l/min, 2=dl/s
        ['tgcMode', 1536+260, 'I'],  # 0=slope, 1=custom
        ['tgcCustomMin', 1536+264, 'i'],
        ['tgcCustomMax', 1536+268, 'i'],
        ['historyClassN', 1536+272, 'I'],
        ['histogramBinN', 1536+276, 'I'],
        ['emitFreq1', 1536+280, 'I'],
        ['emitFreq2', 1536+284, 'I'],
        ['emitFreq3', 1536+288, 'I'],
        ['emitFreq4', 1536+292, 'I'],
        ['emitFreq5', 1536+296, 'I'],
        ['tgcCellDist', 1536+300, 'I'],  # in ns
        ['filterLowPass', 1536+304, 'I'],
        ['_internalUse_11', 1536+308, 'I'],  # exact decoding unknown
        ['gainOverall', 1536+312, 'I'],
        ['emitRecieveSameBNC', 1536+316, 'I'],  # 0=yes
        ['_internalUse_12to48', 1536+320, 37*'I'],  # exact decoding unknown
        ['udvfAreaWidth', 1536+468, 'I'],
        ['udvfAreaHeight', 1536+472, 'I'],
        ['_internalUse_49', 1536+476, 'I'],  # exact decoding unknown
        ['udvfVectorDist', 1536+480, 'I'],
        ['udvfVectorSize', 1536+484, 'I'],
        ['_internalUse_50', 1536+488, 'I'],  # exact decoding unknown
        ['udvfAreaWeqH', 1536+492, '?xxx'],
        ['udvf3d', 1536+496, '?xxx'],
        ['udvfVeloOffsetTR1', 1536+500, 'I'],
        ['udvfVeloOffsetTR2', 1536+504, 'I'],
        ['udvfVeloOffsetTR3', 1536+508, 'I'],
        ['udvfVeloScaleTR1', 1536+512, 'I'],
        ['udvfVeloScaleTR2', 1536+516, 'I'],
        ['udvfVeloScaleTR3', 1536+520, 'I'],
        ['udvfCustomAngleEmitReciev', 1536+524, 'i'], # 'I' also possible
        ['udvfCustomDistEmitRef', 1536+528, 'I'],
        ['udvfCustomDistReceivRef', 1536+532, 'I'],
        ['udvfCustomFreq', 1536+536, 'I'],
        ['udvfProbeNumber', 1536+540, 'I'],  # 0=custom
        ['udvfProbeMode', 1536+544, 'I'],
        ['udvf2d', 1536+548, '?xxx'],
        ['_internalUse_51', 1536+552, 'I'],  # exact decoding unknown
        ['udvf2dDistEmitReceiv', 1536+556, 'I'],  # in 10*mm
        ['udvf2dVeloScaleTR1', 1536+560, 'I'],
        ['udvf2dVeloScaleTR2', 1536+564, 'I'],
        ['udvf2dVeloOffsetTR1', 1536+568, 'I'],
        ['udvf2dVeloOffsetTR2', 1536+572, 'I'],
        ['udvf2dGate1Time', 1536+576, 'I'],  # in ns
        ['udvf2dRes', 1536+580, 'I'],  # in ns
        ['udvf2dGateSkipNDisplay', 1536+584, 'I'],
        ['udvf2dUmax', 1536+588, 'i'],
        ['udvf2dUmin', 1536+592, 'i'],
        ['udvf2dVmax', 1536+596, 'i'],
        ['udvf2dVmin', 1536+600, 'i'],
        ['udvf2dEmitRadius', 1536+604, 'I'],  # in 10*mm
        ['udvf2dEmitHalfAngle', 1536+608, 'i'],  # in degree
        ['udvf2dReceivRadius', 1536+612, 'I'],  # in 10*mm
        ['udvf2dReceivHalfAngle', 1536+616, 'i'],  # in degree
        ['udvf2dDopplerAngle', 1536+620, 'i'],  # in degree
        ['udvf2dEnergyMaxConst', 1536+624, 'i'],  # *1e4
        ['udvf2dEnergyMaxXCoeff', 1536+628, 'i'],  # *1e4

        # multiplexer block at offset 2560
        ['multi_channelUsed', 2560+0, 10*'?xxx'],
        ['multi_profN', 2560+40, '10I'],
        ['multi_prf', 2560+80, '10I'],  # in us
        ['multi_gateN', 2560+120, '10I'],
        ['multi_timeRes', 2560+160, '10I'],  # in ns
        ['multi_emitFreq', 2560+200, '10I'],  # in kHz
        ['multi_emitPower', 2560+240, '10I'],  # 0,1,2=low,medium,high
        ['multi_timeGate1', 2560+280, '10I'],  # in ns
        ['multi_tgcSart', 2560+320, '10i'],
        ['multi_tgcEnd', 2560+360, '10i'],
        ['multi_tgcMode', 2560+400, '10I'],  # 0: slope, 1: custom
        ['multi_burstLength', 2560+440, '10I'],
        ['multi_sensitivity', 2560+480, '10I'],
        ['multi_emitN', 2560+520, '10I'],
        ['multi_veloScale', 2560+560, '10I'],
        ['multi_veloOffset', 2560+600, '10i'],
        ['multi_moduleScale', 2560+640, '10I'],
        ['multi_profType', 2560+680, '10I'],
        ['multi_dopplerAngle', 2560+720, '10i'],
        ['multi_unit', 2560+760, 10*'I'],

        # custom TGC block at offset 3360
        ['tgcChannel1', 3360+0*1024, '1024B'],
        ['tgcChannel2', 3360+1*1024, '1024B'],
        ['tgcChannel3', 3360+2*1024, '1024B'],
        ['tgcChannel4', 3360+3*1024, '1024B'],
        ['tgcChannel5', 3360+4*1024, '1024B'],
        ['tgcChannel6', 3360+5*1024, '1024B'],
        ['tgcChannel7', 3360+6*1024, '1024B'],
        ['tgcChannel8', 3360+7*1024, '1024B'],
        ['tgcChannel9', 3360+8*1024, '1024B'],
        ['tgcChannel10', 3360+9*1024, '1024B']
        ]

    ### Measurement parameters
    _measBaseOffset = 13600
    _measLength = ['length', 0, 'H']
    _measParam = [
        # measurement blocks at offset 13600
        # negative offsets are measured from the end of the block
        ['data', 2, '{:d}b'],  # size dependent on operation parameters
        ['timeStamp', -12, 'I'],  # in us
        ['flowRate', -8, 'I'],  # in ml/min
        ['triggerSequence', -4, 'B'],
        ['channel', -3, 'B'],
        ['length2', -2, 'H'],
        ]
    # calculate the length of the fixed part of every measurement block
    _measFixedLen = struct.calcsize(_measLength[2])
    for param, offset, fmt in _measParam:
        if param != 'data':
            _measFixedLen += struct.calcsize(fmt)

    ### Parameter conversion dictionaries
    _profileTypeNames = {
        0: ['velo'],
        1: ['echo'],
        2: ['energy'],
        3: ['gate FFT'],
        10: ['velo', 'echo'],
        11: ['velo', 'energy'],
        12: ['velo', 'v(t)'],
        13: ['velo', 'flow rate'],
        14: ['velo', 'history'],
        15: ['echo', 'history'],
        16: ['energy', 'history']
        }


    def _read(self):
        """ Read the data in the given BDD file
        """
        ### Read parameters at fixed positions
        for param, offset, fmt in self._fixedParam:
            self._readParam(param, offset, fmt)

        ### Read measurement blocks
        measStart = self._measBaseOffset  # first measurement start
        eof = self._file.seek(0,2)  # end of file

        meas = 1  # current measurement number
        exhausted = False  # set True once file end is reached
        while not exhausted:
            preMeas = self._prefixMeas(meas)

            # read measurement length
            param, offset, fmt = self._measLength
            measStart += offset
            measLen = self._readParam(preMeas+param, measStart, fmt)
            measEnd = measStart + measLen

            # check for end of file
            if measEnd >= eof:
                exhausted = True

            for param, offset, fmt in self._measParam:
                if param == 'data':
                    # set size of data
                    fmt = fmt.format(measLen-self._measFixedLen)

                # set offset from start or end of measurement
                if offset >= 0:
                    offset += measStart
                elif offset < 0:
                    offset += measEnd

                value = self._readParam(preMeas+param, offset, fmt)

                if param == 'channel':
                    # store which measurement belongs to which channel
                    self._chToMeas[value-1].append(meas)

            meas += 1
            measStart = measEnd


    def _refine(self):
        """ Process data read from the BDD file
        """
        ### process information block
        for param in ['version', 'parameter', 'comment']:
            val = self.getParam(param).decode(self._codec)
            val = val.strip('\00\r\n')
            self.setParam(param, val)

        ### process parameter block
        chUsed = np.where(self.getParam('multi_channelUsed'))[0]+1
        self.setParam('channelUsed', chUsed)

        ### process multiplexer block
        for k in list(self.keys()):
            if not k.startswith('multi_'):
                continue
            param = k[6:]
            data = self.getParam(k)
            for i, d in enumerate(data):
                preCh = self._prefixChannel(i+1)
                self.setParam(preCh+param, d)

        ### decode the profileType number into the respective profile types
        # front channel
        pT = self.getParam('profType')
        if pT in self._profileTypeNames:
            self.setParam('profTypeName', self._profileTypeNames[pT])
        else:
            warn('Profile type {} is unknown'.format(pT))
            self.setParam('profTypeName', '')

        # multiplexer channels
        for ci in range(10):
            ch = ci+1
            preCh = self._prefixChannel(ch)
            pT = self.getParam(preCh+'profType')

            if pT in self._profileTypeNames:
                self.setParam(preCh+'profTypeName',
                              self._profileTypeNames[pT])
            else:
                warn('Profile type {} of '.format(pT) +
                     'multiplexer-channel {} is unknown'.format(ch))
                self.setParam(preCh+'profTypeName', '')

        ### assemble measurement blocks into channels and process profiles
        if self.getParam('multi'):
            self._mode = 'multi'
        elif self.getParam('udvf2d'):
            self._mode = 'udvf2d'
        elif self.getParam('udvf3d'):
            self._mode = 'udvf3d'
        else:
            self._mode = 'default'

        for chInd, measList in enumerate(self._chToMeas):
            if measList == []:
                # no profiles recorded for this channel
                continue

            ch = chInd+1  # channel number
            preCh = 'ch{:d}_'.format(ch)
            gateN = self._getGateN(ch)

            # number of profiles
            if self._mode == 'multi':
                profileTypeName = self.getParam(preCh+'profTypeName')
                profileN = len(profileTypeName)
            if self._mode == 'udvf2d':
                profileN = 2
                profileTypeName = ['velo']*profileN
            if self._mode == 'udfv3d':
                profileN = 3
                profileTypeName = ['velo']*profileN
            if self._mode == 'default':
                profileTypeName = self.getParam('profTypeName')
                profileN = len(profileTypeName)

            # assemble profiles and informations
            profile = np.zeros((profileN, len(measList), gateN))
            flowRate = np.zeros(len(measList))
            triggerSeq = np.zeros(len(measList))
            time = np.zeros(len(measList))

            # After 2**32-1 us the time-value (int32) overflows (returns to 0).
            timeOverflow = 2**32-1 # in us (about 1.193 h)
            # This has to be corrected by an additional offset,
            # that gets incremented each time an overflow happens.
            timeOffset = 0

            for i, meas in enumerate(measList):
                preMeas = self._prefixMeas(meas)

                time[i] = self.getParam(preMeas+'timeStamp') + timeOffset
                # detect and correct overflow of the time
                if i > 0 and time[i-1]-time[i] < 0:
                    time[i] += timeOverflow
                    timeOffset += timeOverflow

                prof = self.getParam(preMeas+'data')
                profile[:,i,:] = np.reshape(prof, (profileN, gateN))

                flowRate[i] = self.getParam(preMeas + 'flowRate')
                triggerSeq[i] = self.getParam(preMeas + 'triggerSequence')

            self.setParam(preCh+'depth', self._calcGateDepth(ch))
            self.setParam(preCh+'time', time*1e-9)
            self.setParam(preCh+'flowRate', flowRate)
            self.setParam(preCh+'triggerSeq', triggerSeq)
            self.setParam(preCh+'profile', profile)

            for i, pName in enumerate(profileTypeName):
                pArray = profile[i,:,:]

                if pName == 'velo':
                    pArray, vmax = self._calcVelo(pArray, ch)
                    self.setParam(preCh+'veloMax', vmax)
                elif pName == 'echo':
                    pArray = self._calcEcho(pArray, ch)

                self.setParam(preCh+pName, pArray)


    def _getGateN(self, channel):
        """ Returns the number of gates for a given channel
        """
        if self._mode == 'multi':
            return self.getParam('multi_gateN')[channel-1]
        else:
            return self.getParam('gateN')


    def _calcGateDepth(self, channel):
        """ Gate depths in mm for a given channel
        """
        if self._mode == 'multi':
            soundSpeed = self.getParam('soundSpeed')
            timeRes = self.getParam('multi_timeRes')[channel-1]
            timeGate1 = self.getParam('multi_timeGate1')[channel-1]
            gateN = self.getParam('multi_gateN')[channel-1]
        else:
            soundSpeed = self.getParam('soundSpeed')
            timeRes = self.getParam('timeRes')
            timeGate1 = self.getParam('timeGate1')
            gateN = self.getParam('gateN')

        gateDepth = (np.arange(gateN)*timeRes+timeGate1)*1e-6*soundSpeed/2.

        return gateDepth


    def _calcVelo(self, data, ch):
        """ Velocity and maximum velocity in m/s
        """
        data = np.array(data)

        if self._mode == 'default':
            veloOffset = self.getParam('veloOffset')
            prfPeriod = self.getParam('prfPeriod')
            veloScale = self.getParam('veloScale')
            soundSpeed = self.getParam('soundSpeed')
            emitFreqOff = self.getParam('emitFreqOff')
            emitFreqNumber = int((emitFreqOff-1)/4+1)
            emitFreq = self.getParam('emitFreq{:d}'.format(emitFreqNumber))

            if self.getParam('unit'):
                angle = self.getParam('dopplerAngle')*np.pi/180.
            else:
                angle = 0

        elif self._mode == 'multi':
            chi = ch-1
            veloOffset = self.getParam('multi_veloOffset')[chi]
            prfPeriod = self.getParam('multi_prf')[chi]
            veloScale = self.getParam('multi_veloScale')[chi]
            soundSpeed = self.getParam('soundSpeed')
            emitFreq = self.getParam('multi_emitFreq')[chi]

            if self.getParam('multi_unit')[chi]:
                angle = self.getParam('multi_dopplerAngle')[chi]*np.pi/180.
            else:
                angle = 0

        # correct offset
        data[data + veloOffset > 127] -= 256
        data[data + veloOffset < -128] += 256

        def velocity(data):
            # calculate doppler frequency
            fdoppler = data*1e6 / (128*prfPeriod*veloScale)
            # calculate velocity
            velo = fdoppler*soundSpeed / (2e3*np.cos(angle)*emitFreq)
            return velo

        return velocity(data), velocity(128)


    def _calcEcho(self, data, ch):
        """ Echo amplitude from 0 to 255
        """
        data = np.array(data)

        return data+(data<0)*256



class DOP3000(DOPBase):
    """ Class to read DOP3000 or DOP3010 files

    Based on the DOP3000-3010 User's manual version 6.0, revision: 1
    """
    _standardDepth = 'File'

    ### Parameters at fixed offsets
    _fixedParam = [
        ['version', 0, '16s'],
        ['comment', 16, '512s'],
        ['hardware', 528, '20B'], # exact decoding unknown
        ] + \
        [['ch'+str(i+1)+'_tgc', 10788+i*2048, '2048b'] for i in range(10)]

    ### Channel operation parameters
    # Parameters receive the prefix 'ch<n>_' where <n> is the channel number
    # starting at 1 (see also _prefixChannel).
    _operationBaseOffset = 548
    _operationBlockLen = 4*256
    _operationParam = [
        ['emitFreq', 4*0, 'i'],  # in kHz
        ['assistMode', 4*1, 'i'],  # 1 if True
        ['userDepth', 4*2, 'i'],  # in mm
        ['userVelo', 4*3, 'i'],  # in mm/s
        ['qualityFactor', 4*4, 'i'],  # (0 to 1000)
        ['prf', 4*5, 'i'],  # in us
        ['emitEnable', 4*6, 'i'],  # 1 if True
        ['emitPower', 4*7, 'i'],  # 0 = Low, 1 = Medium, 2 = High
        ['burstLength', 4*8, 'i'],  #
        ['gate1', 4*9, 'i'],  #
        ['resolution', 4*10, 'i'],  # (n+1)*0.166 ns
        ['resolutionAuto', 4*11, 'i'],  #
        ['gateNAuto', 4*12, 'i'],  # 1 if True
        ['gateN', 4*13, 'i'],  # 1 if True
        ['emitNprofile', 4*14, 'i'],  # emissions per profile
        ['veloScale', 4*15, 'i'],  # 0 to 3141
        ['wallFilterCoeff', 4*16, 'i'],  # in 1e-3
        ['emitNstabil', 4*17, 'i'],  #
        ['sensitivity', 4*18, 'i'],  #
        ['soundSpeed', 4*19, 'i'],  # in m/s
        ['dopplerAngle', 4*20, 'i'],  # in deg
        ['moduleScale', 4*21, 'i'],  # 256 to 2048
        ['veloOffset', 4*22, 'i'],  #
        ['tgcMode', 4*23, 'i'],  # 0 = uniform, 1 = slope, 2 = auto, 3 = custom
        ['tgcStart', 4*24, 'i'],  # 0 to 256
        ['tgcEnd', 4*25, 'i'],  # 0 to 256
        ['tgcGateSize', 4*26, 'i'],  # 0 = 0.666 ns, 1 = 1.333 ns
        ['bandwidth', 4*27, 'i'],  # 0 = 50 kHz, 1 = 100 kHz, ..., 5 = 300 kHz
        ['gainOverall', 4*28, 'i'],  # 0 = 0 dB, 1 = 6 dB, 2 = 14 dB, 3 = 20 dB
        ['aquisitionRate', 4*29, 'i'],  # 0 = 6 MHz, 1: 12 or 40 MHz
        ['_internal_1', 4*30, 'i'],  #
        ['profileType', 4*31, 'i'],
        ['display', 4*32, 'i'],  # 0 = horizontal, 1 = vertical
        ['profileNmax', 4*33, 'i'],  #
        #['trigger', 4*34, '4v'],  #
        ['triggerExternal', 4*34, '1m0'],  #
        ['triggerLowLevel', 4*34, '1m1'],  #
        ['profileNtrigger', 4*35, 'i'],  #
        ['filter', 4*36, 'i'],  # 0 = None, 1 = Moving average, 2 = Median
        ['movAvgProfileN', 4*37, 'i'],  #
        ['movAvgRejectZero', 4*38, 'i'],  #
        ['profileNrecord', 4*39, 'i'],  #
        #['record', 4*40, '4v'],  #
        ['recordBinary', 4*40, '1m0'],  #
        ['recordAscii', 4*40, '1m1'],  #
        ['recordStatAscii', 4*40, '1m2'],  #
        ['_internal_2', 4*41, 'i'],  #
        ['_internal_3', 4*42, 'i'],  #
        ['useInOutProbe', 4*43, 'i'],  #
        ['_internal_4', 4*44, 'i'],  #
        ['_internal_5', 4*45, 'i'],  #
        ['hardwareDelay', 4*46, 'i'],  # in ns
        ['triggerDelay', 4*47, 'i'],  # in ms
        ['profileNblock', 4*48, 'i'],  #
        ['blockN', 4*49, 'i'],  #
        ['recordAuto', 4*50, 'i'],  #
        ['multi_profileNblock', 4*51, 'i'],  #
        # ['multi_info', 4*52, '4v'],  #
        ['multi', 4*52, '4m0'],  #
        ['udvmd', 4*52, '4m1'],  #
        ['udv3d', 4*52, '4m2'],  # 2d otherwise
        ['useCh1', 4*52, '4m3'],  #
        ['useCh2', 4*52, '4m4'],  #
        ['useCh3', 4*52, '4m5'],  #
        ['useCh4', 4*52, '4m6'],  #
        ['useCh5', 4*52, '4m7'],  #
        ['useCh6', 4*52, '4m8'],  #
        ['useCh7', 4*52, '4m9'],  #
        ['useCh8', 4*52, '4m10'],  #
        ['useCh9', 4*52, '4m11'],  #
        ['useCh10', 4*52, '4m12'],  #
        ['channelFirst1', 4*52, '4m14'],  #
        ['channelFirst2', 4*52, '4m15'],  #
        ['channelFirst3', 4*52, '4m16'],  #
        ['channelFirst4', 4*52, '4m17'],  #
        ['channelFirst5', 4*52, '4m18'],  #
        ['channelFirst6', 4*52, '4m19'],  #
        ['channelFirst7', 4*52, '4m20'],  #
        ['channelFirst8', 4*52, '4m21'],  #
        ['channelFirst9', 4*52, '4m22'],  #
        ['channelFirst10', 4*52, '4m23'],  #
        ['channelFirst11', 4*52, '4m24'],  #
        ['channelFirst12', 4*52, '4m25'],  #
        ['channelFirst13', 4*52, '4m26'],  #
        ['channelFirst14', 4*52, '4m27'],  #
        ['channelFirst15', 4*52, '4m28'],  #
        ['channelFirst16', 4*52, '4m29'],  #
        ['channelFirst17', 4*52, '4m30'],  #s
        ['multi_blockN', 4*53, 'i'],  #
        ['medianProfileN', 4*54, 'i'],  #
        ['echoMode', 4*55, 'i'],  #
        ['-internal_6', 4*56, 'i'],  #
        ['_internal_7', 4*57, 'i'],  #
        ['gateFftPointN', 4*58, 'i'],  #
        ['gateFftHaming', 4*59, 'i'],  #
        ['_internal_8', 4*60, 'i'],  #
        ['vtTimeScale', 4*61, 'i'],  # in ms
        ['soundSpeedDist', 4*62, 'i'],  # in 0.1 m
        ['udvmd_emitProbe', 4*63, 'i'],  #
        ['udvmd_receiveProbe', 4*64, 'i'],  #
        ['udvmd_distProbe', 4*65, 'i'],  # in 0.1 mm
        ['udvmd_dopplerAngle', 4*66, 'i'],  # in deg
        ['udvmd_veloXmax', 4*67, 'i'],  # in 0.1 mm/s
        ['udvmd_veloYmax', 4*68, 'i'],  # in 0.1 mm/s
        ['udvmd_veloZmax', 4*69, 'i'],  # in 0.1 mm/s
        ['udvmd_tr1FreqOff', 4*70, 'i'],  # -127 to 128
        ['udvmd_tr2FreqOff', 4*71, 'i'],  # -127 to 128
        ['udvmd_tr3FreqOff', 4*72, 'i'],  # -127 to 128
        ['udvmd_tr1ScaleAngle', 4*73, 'i'],  # in 1e-3 rad
        ['udvmd_tr2ScaleAngle', 4*74, 'i'],  # in 1e-3 1000
        ['udvmd_tr3ScaleAngle', 4*75, 'i'],  # in 1e-3 1000
        ['udvmd_tr1ModuleScale', 4*76, 'i'],  #
        ['udvmd_tr2ModuleScale', 4*77, 'i'],  #
        ['udvmd_tr3ModuleScale', 4*78, 'i'],  #
        ['udvmd_depthFrom', 4*79, 'i'],  # in mm
        ['udvmd_depthTo', 4*80, 'i'],  # in mm
        ['udvmd_resolution', 4*81, 'i'],  # in 0.1 mm
        ['udvmd_qualityFactor', 4*82, 'i'],  #
        ['channelCurrent', 4*83, 'i'],  #
        ['profileSkip', 4*84, 'i'],  #
        ['simul_veloScaleFactor', 4*85, 'i'],  #
        ['liquidAttenuation', 4*86, 'i'],  # in 10 dB/m
        ['simul_gateIQ', 4*87, 'i'],  #
        ['_internal_9', 4*88, 'i'],  #
        ['udvmd_vectorSkipN', 4*89, 'i'],  #
        ['rawData_prfN', 4*90, 'i'],  #
        ['rawData_gateN', 4*91, 'i'],  #
        ['rawData_gate1', 4*92, 'i'],  #
        ['gainOverall2', 4*93, 'i'],  # in dB
        ['cursor1', 4*94, 'i'],  # on gate
        ['cursor2', 4*95, 'i'],  # on gate
        ['cursor3', 4*96, 'i'],  # on gate
        ['cursor4', 4*97, 'i'],  # on gate
        ['udv2d3dsim_from', 4*98, 'i'],  #
        ['udv2d3dsim_to', 4*99, 'i'],  #
        ['udv2d3dsim_sourceXpos', 4*100, 'i'],  #
        ['udv2d3dsim_sourceYpos', 4*101, 'i'],  #
        ['udv2d3dsim_sourceZpos', 4*102, 'i'],  #
        ['simul_attenuationMax', 4*103, 'i'],  # in dB
        ['udv2d3dsim_partDensity', 4*104, 'i'],  #
        ['udv2d3dsim_veloFieldType', 4*105, 'i'],  #
        ['udv2d3dsim_rotSpeed', 4*106, 'i'],  #
        ['udv2d3dsim_tr1Shift', 4*107, 'i'],  #
        ['udv2d3dsim_tr2Shift', 4*108, 'i'],  #
        ['udv2d3dsim_tr3Shift', 4*109, 'i'],  #
        ['udv2d3dsim_prfAdapt', 4*110, 'i'],  #
        ['flowRateUnit', 4*111, 'i'],  # 0 = ml/min, 1 = ml/s, 2 = dl/s
        ['flowRateScale', 4*112, 'i'],  #
        ['aliasAutorCorr1', 4*113, 'i'],  #
        ['aliasAutorCorr2', 4*114, 'i'],  #
        ['wallSoundSpeed', 4*115, 'i'],  #
        ['wallLength', 4*116, 'i'],  #
        ['couplantSoundSpeed', 4*117, 'i'],  #
        ['couplantLength', 4*118, 'i'],  #
        ['emitFreqDivider', 4*119, 'i'],  #
        ['demodFreqDivider', 4*120, 'i'],  #
        ['dop4000_crtlReg1', 4*121, 'i'],  #
        ['dop4000_crtlReg2', 4*122, 'i'],  #
        ['histClassN', 4*123, 'i'],  #
        ]

    ### Measurement parameter
    # Negative offsets are measured from the end of the measurement.
    # The following parameters receive the prefix 'meas<n>_' where <n> is the
    # measurement number starting at 1 (see also _prefixMeas).
    _measBaseOffset = 31268
    _measLen = ['length', 0, 'H']
    _measInfoParam = [
        ['timeStamp', -12, 'I'],  # in ms/10
        ['block', -8, 'H'],
        ['markByte', -6, 'B'],
        ['triggerState', -5, 'B'],
        ['_internal', -4, 'B'],  # no meaning
        ['channel', -3, 'B'],
        ['length2', -2, 'H']
        ]
    # The following parameters receive the prefix 'meas<n>_prof<m>_' where <n>
    # is the measurement number starting at 1 (see also _prefixMeas) and <m>
    # is the profile number starting at 1 (see also _prefixProfile).
    _measProfParam = [
        ['length', 0, 'H'],  # in bytes
        ['type', 2, 'B'],
        ['data', 3, '{length:.0f}{fmt:}'],
        ]

    ### Parameter conversion dictionaries
    _profileTypeNames = {
         0: 'velo',
         1: 'echo',
         2: 'energy',
         3: 'gateFFT',
         4: 'phase',  # in deg/10
         5: 'vitson',
         6: 'frequency',
         7: 'frequencyTR1',
         8: 'frequencyTR2',
         9: 'frequencyTR3',
         10: 'echoTR1',
         11: 'echoTR2',
         12: 'echoTR3',
         13: 'energyTR1',
         14: 'energyTR2',
         15: 'energyTR3',
         16: 'velo_mm/s/10',
         17: 'velo_mm/s',
         18: 'flow',  # in ml/min
         19: 'diameter',  # in 10
         20: 'aliasingRefer',
         21: 'aliasingReferTR1',
         22: 'aliasingReferTR2',
         23: 'aliasingReferTR3',
         24: 'frequency_10Hz',
         25: 'depth',  # in mm/10
         28: 'tgc',
         29: 'iq',
         30: 'angle',  # in deg/10
         }

    _profileTypeFmt = {
         0 : 'b', 1 : 'B', 2 : 'b', 3 : 'b', 4 : 'h', 5 : 'b', 6 : 'b',
         7 : 'b', 8 : 'b', 9 : 'b', 10: 'b', 11: 'b', 12: 'b', 13: 'b',
         14: 'b', 15: 'b', 16: 'h', 17: 'h', 18: 'i', 19: 'h', 20: 'b',
         21: 'b', 22: 'b', 23: 'b', 24: 'h', 25: 'h', 28: 'b', 29: 'b',
         30: 'h',
         }


    def _byteToBit(self, byteString):
        """ Convert a buffer of bytes to a bit representation
        """
        res = ''
        if len(byteString) != 0:
            for byte in byteString:
                if not isinstance(byte, int):
                    byte = ord(byte)
                for i in range(8):
                    res += str((byte >> i) & 1)
        return res


    def _readMeas(self, meas, measStart):
        """ Read a measurement from file

        Arguments:
        ==========
        meas: int
            Measurement number
        measStart: int
            Offset of the measurement in the file

        Returns:
        ========
        measEnd: int
            Offset of the measurement end (= start of the next measurement)
            in bytes.
        """
        preMeas = self._prefixMeas(meas)

        # read profile length
        param, offset, fmt = self._measLen
        measStart += offset
        measLen = self._readParam(preMeas+param, measStart, fmt)
        measEnd = measStart + measLen

        profile = 1  # current profile
        profStart = measStart + struct.calcsize(fmt)
        profLen = 1
        while profLen != 0:
            preProf = preMeas + self._prefixProfile(profile)
            profEnd = profStart

            # read profile length
            param, offset, fmt = self._measProfParam[0]
            profEnd += struct.calcsize(fmt)
            profLen = self._readParam(preProf+param, profStart+offset, fmt)
            if profLen == 0:
                # no more profiles in this measurement
                break

            # read profile type
            param, offset, fmt = self._measProfParam[1]
            profEnd += struct.calcsize(fmt)
            profType = self._readParam(preProf+param, profStart+offset, fmt)
            profFmt = self._profileTypeFmt[profType]
            self.setParam(preProf+'format', profFmt)

            # read profile data
            param, offset, fmt = self._measProfParam[2]
            fmt = fmt.format(length=profLen/struct.calcsize(profFmt),
                             fmt=profFmt)
            profEnd += struct.calcsize(fmt)
            self._readParam(preProf+param, profStart+offset, fmt)

            # next profile
            profile += 1
            profStart = profEnd

        self.setParam(preMeas+'profN', profile-1) # number of profiles

        # read measurement information
        for param, offset, fmt in self._measInfoParam:
            # set offset from end of measurement
            offset += measEnd
            value = self._readParam(preMeas+param, offset, fmt)

            if param == 'channel':
                # store which measurement belongs to which channel
                self._chToMeas[value-1].append(meas)

        if self.getParam(preMeas+'length') != self.getParam(preMeas+'length2'):
            warn('Lengths in measurement {:d} do not match!'.format(meas))

        return measEnd


    def _read(self):
        """ Read the data in the given BDD file
        """
        self._file.seek(0,2)
        eof = self._file.tell()  # end of file

        ### read parameters at fixed positions
        for param, offset, fmt in self._fixedParam:
            self._readParam(param, offset, fmt)

        ### read channel parameters
        for ch in range(1,11):
            preCh = self._prefixChannel(ch)
            baseOffset = self._operationBaseOffset + \
                         (ch-1)*self._operationBlockLen
            for param, offset, fmt in self._operationParam:
                self._readParam(preCh+param, baseOffset+offset, fmt)

        ### read measured profiles
        measStart = self._measBaseOffset  # first block start

        meas = 1  # current measurement number
        exhausted = False  # set True once file end is reached
        while not exhausted:
            measEnd = self._readMeas(meas, measStart)

            # check for end of file
            if measEnd >= eof:
                exhausted = True

            measStart = measEnd
            meas += 1

        self.setParam('measN', meas-1)


    def _refine(self):
        """ Process data read from the BDD file
        """
        ### process information parameters
        for param in ['version', 'comment']:
            val = self.getParam(param).decode(self._codec)
            val = val.strip('\00\r\n')
            self.setParam(param, val)

        channelUsed = []
        ### assemble measurements for each channel
        for chInd, measList in enumerate(self._chToMeas):
            ch = chInd+1  # channel number
            preCh = self._prefixChannel(ch)#

            if len(measList) == 0:
                # no profiles recorded for this channel
                continue
            else:
                channelUsed.append(ch)

            # assemble profiles and informations
            # After 2**32-1 us the time-value (int32) overflows (returns to 0).
            timeOverflow = 2**32-1 # in us (about 1.193 h)
            # This has to be corrected by an additional offset,
            # that gets incremented each time an overflow occurs.
            timeOffset = 0

            time = []
            triggerState = []
            depth = []
            profName = []
            profData = []
            for i, meas in enumerate(measList):

                preMeas = self._prefixMeas(meas)

                # assemble measurement specific information
                profN = self.getParam(preMeas+'profN')

                for profile in range(1,profN+1):
                    preProf = preMeas + self._prefixProfile(profile)
                    type_ = self.getParam(preProf+'type')

                    name = self._profileTypeNames[type_]
                    self.setParam(preProf+'name', name)

                    data = self.getParam(preProf+'data')

                    if name == 'depth':
                        depth = np.array(data)/10
                    else:
                        if len(profData) < profile:
                            profData.append([data])
                            profName.append(name)
                        else:
                            profData[profile-1].append(data)

                        if profile == 1:
                            triggerState.append(
                                self.getParam(preMeas+'triggerState'))
                            time.append(self.getParam(preMeas+'timeStamp') +
                                        timeOffset)
                            # detect and correct overflow of the time
                            if len(time) > 1:
                                if time[-1]-time[-2] < 0:
                                    time[-1] += timeOverflow
                                    timeOffset += timeOverflow

            self.setParam(preCh+'depthFile', np.array(depth))
            self.setParam(preCh+'depthCalc', self._calcDepth(ch))
            self.setParam(preCh+'time', np.array(time)*1e-4)
            self.setParam(preCh+'triggerState', np.array(triggerState))
            self.setParam(preCh+'profTypeName', profName)

            for ip, data in enumerate(profData):
                name = profName[ip]

                if name == 'velo':
                    data, vmax = self._calcVelo(data, ch)
                    self.setParam(preCh+'veloMax', vmax)

                self.setParam(preCh+name, np.array(data))

        self.setParam('channelUsed', channelUsed)


    def _calcDepth(self, channel):
        """ Calculate gate depth in mm for a given channel
        """
        preCh = self._prefixChannel(channel)

        gateN = self.getParam(preCh+'gateN')
        gate1 = self.getParam(preCh+'gate1')  # Par[9]
        resolution = self.getParam(preCh+'resolution')  # Par[10]
        soundSpeed = self.getParam(preCh+'soundSpeed')  # Par[19]
        aquisitionRate = self.getParam(preCh+'aquisitionRate') % 256 # Par[29]
        hardwareDelay = self.getParam(preCh+'hardwareDelay')  # Par[46]

        if aquisitionRate == 0:
            aquisitionRate = 6e3
        else:
            raise ValueError('I don\'t know what aquisition rate ' +
                             '{} means.'.fromat(aquisitionRate))

        gateNumber = np.arange(1,gateN+1)
        term1 = (gate1+(resolution+1)*(gateNumber-1)) / (2*aquisitionRate)
        term2 = hardwareDelay / 2e6
        depth = soundSpeed * (term1-term2)

        return depth


    def _calcVelo(self, data, channel):
        """ Velocity and maximum velocity in m/s
        """
        preCh = self._prefixChannel(channel)
        data = np.array(data)

        # correct offset
        veloOffset = self.getParam(preCh+'veloOffset')
        data[data + veloOffset > 127] -= 256
        data[data + veloOffset < -128] += 256

        def velocity(data):
            # calculate doppler frequency [Hz]
            prf = self.getParam(preCh+'prf')
            veloScale = self.getParam(preCh+'veloScale')
            fdoppler = data*veloScale*1e3 / (256*np.pi*prf)

            # calculate velocity [mm/s]
            emitFreq = self.getParam(preCh+'emitFreq')
            soundSpeed = self.getParam(preCh+'soundSpeed')
            dopplerAngle = self.getParam(preCh+'dopplerAngle') * np.pi/180.
            velo = fdoppler*soundSpeed / (2e3*np.cos(dopplerAngle)*emitFreq)

            return velo

        return velocity(data), velocity(128)


    def getDepth(self, channel=None, version='File'):
        """ Returns the gate depths for given channels in mm

        Arguments:
        ==========
        channel: int or list
            A channel number (1 to 10) or a list of channel numbers.
        version: {'File', 'Calc'}
            Whether the values stored in the file or values calculated
            separately should be returned.

        Returns:
        ========
        depth: array or list
            If `channel` is an integer an array is returned. If `channel` is a
            list of ints, a list of arrays is returned with each element
            corresponding to the same element given in `channel`.
        """
        return self.getChannelParam('depth'+version, channel)



def DOP(fname):
    """ Reads a binary DOP-file (*.BDD)

    DOP2000 and DOP3000/3010-files are supported. The file may be a *.gz or
    *.bz2-archive.

    Arguments:
    ==========
    fname: str
        Absolute or relative path to the file.
    """
    # open file
    if fname.endswith('.bz2'):
        f = bz2.BZ2File(fname, 'rb')
    elif fname.endswith('.gz'):
        f = gzip.GzipFile(fname, 'rb')
    else:
        f = open(fname, 'rb')

    # read file version
    version = f.readline()
    f.close()

    # process file version
    if version.startswith(b'BINWDOPV'):
        return DOP2000(fname)
    elif version.startswith(b'BINUDOPV'):
        return DOP3000(fname)
    else:
        raise Exception('BDD version {} '.format(repr(version)) +
                        'of file {} '.format(repr(fname)) + 'is unknown.')
