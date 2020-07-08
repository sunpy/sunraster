
from astropy.io import fits
import astropy.units as u


__all__ = ["SPICEMeta"]


class SPICEMeta(fits.header.Header, SlitSpectrographMetaABC):
    def _get_unit(self, key):
        try:
            return [s.split("]") for s in self.comments[key].split("[")[1:]][0][:-1][0]
        except IndexError:
            return None
        
    def _construct_quantity(self, key):
        val = self.get(key, None)
        if val:
            val *= u.Unit(self._get_unit(key))
        return val
    
    def _construct_time(self, key):
        val = self.get(key, None)
        scale = self._get_unit(key).lower()
        if val:
            val = Time(val, format="fits", scale=scale)
        return val
        
    @property
    def spectral_window(self):
        spectral_window = self["EXTNAME"]
        redundant_txt = "WINDOW"
        if redundant_txt in spectral_window:
            spectral_window = np.asanyarray(spectral_window.split("_"))
            idx = np.array([redundant_txt not in window_chunk for window_chunk in spectral_window])
            spectral_window = "_".join(spectral_window[idx])
        return spectral_window
    
    @property
    def detector(self):
        return self.get("DETECTOR", None)

    @property
    def instrument(self):
        return self.get("INSTRUME", None)
    
    @property
    def observatory(self):
        return self.get("OBSRVTRY", None)
    
    @property
    def processing_level(self):
        return self.get("LEVEL", None)
    
    @property
    def rsun_meters(self):
        return self._construct_quantity("RSUN_REF")
    
    @property
    def rsun_obs(self):
        return self._construct_quantity("RSUN_ARC")
    
    @property
    def obsid(self):
        return self.get("OBS_ID", None)
    
    @property
    def obsvr(self):
        return self.get("OBS_VR", None)

    @property
    def dsun(self):
        return self._construct_quantity("DSUN_OBS")
    
    @property
    def date(self):
        return self._construct_time("DATE-OBS")

    @property
    def observer_coordinate(self):
        lon_unit = u.deg
        lat_unit = u.deg
        radius_unit = u.m
        lon_key = "HGLN_OBS"
        lat_key = "HGLT_OBS"
        kwargs = {'lon': u.Quantity(self.get(lon_key), unit=self._get_unit(lon_key)).to_value(lon_unit),
                  'lat': u.Quantity(self.get(lat_key), unit=self._get_unit(lat_key)).to_value(lat_unit),
                  'radius': self.dsun.to_value(radius_unit),
                  'unit': (lon_unit, lat_unit, radius_unit),
                  'frame': "heliographic_stonyhurst"}
        return SkyCoord(obstime=self.date, **kwargs)
    
    def __str__(self):
        return textwrap.dedent(f"""\
                SPICEMeta
                ---------
                Observatory:\t\t{self.observatory}
                Instrument:\t\t{self.instrument}
                Detector:\t\t{self.detector}
                Spectral Window:\t{self.spectral_window}
                Date:\t\t\t{self.date}
                OBS ID:\t\t\t{self.obsid}
                """)

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"
