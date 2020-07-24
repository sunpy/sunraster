import textwrap

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.coordinates import HeliographicStonyhurst
from ndcube import NDCollection

from sunraster.meta import Meta, SlitSpectrographMetaABC


__all__ = ["SPICEMeta"]


class SPICEMeta(Meta, metaclass=SlitSpectrographMetaABC):
    # ---------- SPICE-specific convenience methods ----------
    def _get_unit(self, key):
        try:
            return [s.split("]") for s in self.get_comment(key).split("[")[1:]][0][:-1][0]
        except IndexError:
            return None

    def _construct_quantity(self, key):
        val = self.get(key)
        if val:
            val *= u.Unit(self._get_unit(key))
        return val

    def _construct_time(self, key):
        val = self.get(key)
        scale = self._get_unit(key).lower()
        if val:
            val = Time(val, format="fits", scale=scale)
        return val

    def __str__(self):
        return textwrap.dedent(f"""\
                SPICEMeta
                ---------
                Observatory:\t\t{self.observatory}
                Instrument:\t\t{self.instrument}
                Detector:\t\t{self.detector}
                Spectral Window:\t{self.spectral_window}
                Date:\t\t\t{self.date_reference}
                OBS ID:\t\t\t{self.observing_mode_id}
                """)

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    # ---------- Inherited ABC properties ----------
    @property
    def spectral_window(self):
        spectral_window = self.get("EXTNAME")
        redundant_txt = "WINDOW"
        if redundant_txt in spectral_window:
            spectral_window = np.asanyarray(spectral_window.split("_"))
            idx = np.array([redundant_txt not in window_chunk for window_chunk in spectral_window])
            spectral_window = "_".join(spectral_window[idx])
        return spectral_window

    @property
    def detector(self):
        return self.get("DETECTOR")

    @property
    def instrument(self):
        return self.get("INSTRUME")

    @property
    def observatory(self):
        return self.get("OBSRVTRY")

    @property
    def processing_level(self):
        return self.get("LEVEL")

    @property
    def rsun_meters(self):
        return self._construct_quantity("RSUN_REF")

    @property
    def rsun_angular(self):
        return self._construct_quantity("RSUN_ARC")

    @property
    def observing_mode_id(self):
        return self.get("SPIOBSID")

    @property
    def observatory_radial_velocity(self):
        return self.get("OBS_VR")

    @property
    def distance_to_sun(self):
        return self._construct_quantity("DSUN_OBS")

    @property
    def date_reference(self):
        return self._construct_time("DATE-OBS")

    @property
    def date_start(self):
        return self._construct_time("DATE-BEG")

    @property
    def date_end(self):
        return self._construct_time("DATE-END")

    @property
    def observer_coordinate(self):
        lon_unit = u.deg
        lat_unit = u.deg
        radius_unit = u.m
        lon_key = "HGLN_OBS"
        lat_key = "HGLT_OBS"
        kwargs = {'lon': u.Quantity(self.get(lon_key), unit=self._get_unit(lon_key)).to_value(lon_unit),
                  'lat': u.Quantity(self.get(lat_key), unit=self._get_unit(lat_key)).to_value(lat_unit),
                  'radius': self.distance_to_sun.to_value(radius_unit),
                  'unit': (lon_unit, lat_unit, radius_unit),
                  'frame': HeliographicStonyhurst}
        return SkyCoord(obstime=self.date_reference, **kwargs)

    @property
    def version(self):
        return self.get("VERSION")

    # ---------- SPICE-specific metadata properties ----------
    @property
    def observing_mode_id_solar_orbiter(self):
        return self.get("OBS_ID")

    @property
    def darkmap_subtracted_onboard(self):
        return bool(self.get("DARKMAP"))

    @property
    def bias_frame_subtracted_onboard(self):
        return bool(self.get("BLACKLEV"))

    @property
    def window_type(self):
        return self.get("WIN_TYPE")

    @property
    def window_table_id(self):
        return self.get("WINTABID")

    @property
    def slit_id(self):
        return self.get("SLIT_ID")

    @property
    def slit_width(self):
        return self._construct_quantity("SLIT_WID")

    @property
    def dumbbell(self):
        dumbell_types = ["none", "lower", "upper"]
        dumbell_idx = self.get("DUMBBELL")
        return dumbell_types[dumbell_idx]

    @property
    def solar_B0(self):
        """Tilt angle of solar north toward spacecraft."""
        return self._construct_quantity("SOLAR_B0")

    @property
    def solar_P0(self):
        """Angle from spacecraft celestial north to solar north."""
        return self._construct_quantity("SOLAR_P0")

    @property
    def solar_ep(self):
        """Angle from spacecraft ecliptic north to solar north angle."""
        return self._construct_quantity("SOLAR_EP")

    @property
    def carrington_rotation(self):
        """Carrington Rotation number of observation."""
        return self.get("CAR_ROT")

    @property
    def date_start_earth(self):
        """Time at which photons reaching SPICE at start time would have reach Earth."""
        return self._construct_time("DATE_EAR")

    @property
    def date_start_sun(self):
        """
        Time at which photons reaching SPICE at start time would have left Sun.

        The Sun is defined as the center of the Sun assuming photon was not impeded.
        """
        return self._construct_time("DATE_SUN")
