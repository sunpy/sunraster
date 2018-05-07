# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from pkg_resources import resource_filename


class Obsid(object):
    """
    A class to convert the IRIS OBS ID to human-readable format.
    """
    def __init__(self, obsid):
        self.obsid = obsid
        self.field_keys = {'Large linelist': 'linelist',
                           'Default compression': 'compression',
                           'Non-simultaneous readout': 'readout',
                           'SJI cadence default': 'sji_cadence',
                           'FUV binned same as NUV': 'fuv_binning',
                           'Spatial x 1, Spectral x 1': 'binning',
                           'Exposure 1s': 'exptime',
                           'C II   Si IV   Mg II h/k   Mg II w   ': 'sjis'}
        for value in self.field_keys.values():
            setattr(self, value, None)
        self.read_obsid(obsid)

    def __repr__(self):
        quants = {value: getattr(self, value)
                  for value in (list(self.field_keys.values()) +
                  ['obsid', 'spec_cadence', 'raster_fov', 'raster_fulldesc',
                   'sji_fov'])}
        return ("IRIS OBSID {obsid}\n"
                "---------------------\n"
                "Description:       {raster_fulldesc:>45}\n"
                "SJI filters:       {sjis:>45}\n"
                "SJI field of view: {sji_fov:>45}\n"
                "Exposure time:     {exptime:>45}\n"
                "Binning:           {binning:>45}\n"
                "FUV binning:       {fuv_binning:>45}\n"
                "SJI cadence:       {sji_cadence:>45}\n"
                "Compression:       {compression:>45}\n"
                "Linelist:          {linelist:>45}").format(**quants)

    def read_obsid(self, obsid):
        """
        Reads different fields from OBSID number.
        """
        # here choose between tables
        version = int(str(obsid)[:2])
        obsid = int(str(obsid)[2:])
        TAB1 = resource_filename('irispy', 'data/v%i-table10.csv' % version)
        TAB2 = resource_filename('irispy', 'data/v%i-table2000.csv' % version)
        table1 = pd.read_csv(TAB1)
        table2 = pd.read_csv(TAB2)
        id_raster = int(str(obsid)[-2:])
        try:
            meta = table1.where(table1['OBS-ID'] == id_raster).dropna().iloc[0]
        except IndexError:
            raise ValueError("Invalid OBSID: last two numbers must be between"
                             "{0} and {1}".format(table1['OBS-ID'].min(),
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
        for r_i, r_f in zip(field_ranges[-2::-1], field_ranges[:0:-1]):
            table = table2.iloc[r_i:r_f]
            for i in np.arange(r_i, r_f)[::-1]:
                index = i
                tmp = table['OBS ID'].loc[i]
                if (obsid - tmp) >= 0:
                    obsid -= tmp
                    break
            description = table['Size + description']
            if description.iloc[0] in self.field_keys:
                setattr(self, self.field_keys[description.iloc[0]],
                        description.loc[index].strip())
