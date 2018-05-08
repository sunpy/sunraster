# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from astropy import units as u
from pkg_resources import resource_filename


class Obsid(object):
    """A class to convert the IRIS OBS ID to human-readable format.

    Parameters
    ----------
    obsid: `int`
        IRIS OBS ID to query. Needs to be a valid OBS ID, with 10 digits.

    Notes
    -----
    Currently, only OBS IDs with version numbers 36, 38, and 40 are
    supported. Calibration OBS IDs (version 42) are not yet supported.

    Examples
    --------
    Quickly show OBS ID parameters:

    >>> obsid.Obsid(3675508564)
    IRIS OBS ID 3675508564
    ----------------------
    Description:            Large dense 96-step raster 31.35x120 96s
    SJI filters:                                    C II   Mg II w s
    SJI field of view:                                       120x120
    Exposure time:                                               8.0 s
    Binning:                               Spatial x 1, Spectral x 1
    FUV binning:                         FUV spectrally rebinned x 4
    SJI cadence:                                 SJI cadence default
    Compression:                                Lossless compression
    Linelist:                                       Flare linelist 1

    Access values:

    >>> data = obsid.Obsid(3675508564)
    >>> data.exptime
    <Quantity 8. s>
    >>> data.linelist
    'Flare linelist 1'
    >>> data.raster_fov
    '31.35x120'
    >>> data.raster_step
    0.33
    """

    def __init__(self, obsid):
        self.obsid = obsid
        self.field_keys = {'Large linelist': 'linelist',
                           'Default compression': 'compression',
                           'Lossy compression': 'compression',
                           'Non-simultaneous readout': 'readout',
                           'SJI cadence default': 'sji_cadence',
                           'SJI cadence 10s': 'sji_cadence',
                           'FUV binned same as NUV': 'fuv_binning',
                           'Spatial x 1, Spectral x 1': 'binning',
                           'Exposure 1s': 'exptime',
                           'C II   Si IV   Mg II h/k   Mg II w   ': 'sjis'}
        self._versions = [36, 38, 40]
        for value in self.field_keys.values():
            setattr(self, value, None)
        self.read_obsid(obsid)

    def __repr__(self):
        quants = {value: getattr(self, value)
                  for value in (list(self.field_keys.values()) +
                  ['obsid', 'spec_cadence', 'raster_fov', 'raster_fulldesc',
                   'sji_fov'])}
        return ("IRIS OBS ID {obsid}\n"
                "----------------------\n"
                "Description:       {raster_fulldesc:>45}\n"
                "SJI filters:       {sjis:>45}\n"
                "SJI field of view: {sji_fov:>45}\n"
                "Exposure time:     {exptime:>45}\n"
                "Binning:           {binning:>45}\n"
                "FUV binning:       {fuv_binning:>45}\n"
                "SJI cadence:       {sji_cadence:>45}\n"
                "Compression:       {compression:>45}\n"
                "Linelist:          {linelist:>45}").format(**quants)

    @staticmethod
    def __exptime_to_quant(exptime):
        """
        Converts an 'exptime' string (used in IRIS tables and OBS_DESC)
        to a Quantity instance in seconds.
        """
        if exptime == "Exposure 1s":
            return 1. * u.s
        else:
            return float(exptime.split(' x ')[1]) * u.s

    def read_obsid(self, obsid):
        """
        Reads different fields from OBS ID number.
        """
        # here choose between tables
        if len(str(obsid)) != 10:
            raise ValueError("Invalid OBS ID: must have 10 digits.")
        version = int(str(obsid)[:2])
        if version not in self._versions:
            raise ValueError("Invalid OBS ID: two first digits must one of"
                             " {0}".format(self._versions))
        obsid = int(str(obsid)[2:])  # version digits are no longer needed
        table1 = pd.read_csv(resource_filename('irispy',
                                               'data/v%i-table10.csv' % version))
        table2 = pd.read_csv(resource_filename('irispy',
                                               'data/v%i-table2000.csv' % version))
        id_raster = int(str(obsid)[-2:])
        try:
            meta = table1.where(table1['OBS-ID'] == id_raster).dropna().iloc[0]
        except IndexError:
            raise ValueError("Invalid OBS ID: last two numbers must be between"
                             " {0} and {1}".format(table1['OBS-ID'].min(),
                                                   table1['OBS-ID'].max()))
        self.raster_step = meta['Raster step']
        self.raster_fov = meta['Raster FOV']
        self.spec_cadence = meta['Spectral cadence']
        self.sji_fov = meta['SJI FOV']
        self.raster_desc = meta['Description']
        self.raster_fulldesc = '%s %s %s' % (self.raster_desc, self.raster_fov,
                                             self.spec_cadence)
        self.raster_meta = meta.copy()
        field_ranges = np.concatenate([  # find all dividers between fields
            table2.where(table2['OBS ID'] == 0).dropna(how='all').index,
            np.array([len(table2)])])
        # field indices, start from largest and subtract
        for ri, rf in zip(field_ranges[-2::-1], field_ranges[:0:-1]):
            table = table2.iloc[ri:rf]
            for i in np.arange(ri, rf)[::-1]:
                index = i
                tmp = table['OBS ID'].loc[i]
                if (obsid - tmp) >= 0:
                    obsid -= tmp
                    break
            desc = table['Size + description']
            # Save values for attributes but also table options as function of OBS ID
            if desc.iloc[0] in self.field_keys:
                attr_name = self.field_keys[desc.iloc[0]]
                if attr_name == 'exptime':
                    options = (self.__exptime_to_quant(a) for a in list(desc.values))
                    options = dict(zip(options, table['OBS ID']))
                    attr_value = self.__exptime_to_quant(desc.loc[index])
                else:
                    options = dict(zip(desc, table['OBS ID']))
                    attr_value = desc.loc[index].strip()
                setattr(self, attr_name + '_options', options)
                setattr(self, attr_name, attr_value)
