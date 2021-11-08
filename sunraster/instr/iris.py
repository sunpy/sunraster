import textwrap

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from ndcube import NDCollection
from sunpy.coordinates import Helioprojective

from sunraster import RasterSequence, SpectrogramCube
from sunraster.meta import Meta, SlitSpectrographMetaABC

__all__ = ["read_iris_spectrograph_level2_fits"]

# Define some properties of IRIS detectors. Source: IRIS instrument paper.
DETECTOR_GAIN = {"NUV": 18.0, "FUV": 6.0}
DETECTOR_YIELD = {"NUV": 1.0, "FUV": 1.5}
DN_UNIT = {
    "NUV": u.def_unit("DN_IRIS_NUV", DETECTOR_GAIN["NUV"] / DETECTOR_YIELD["NUV"] * u.photon),
    "FUV": u.def_unit("DN_IRIS_FUV", DETECTOR_GAIN["FUV"] / DETECTOR_YIELD["FUV"] * u.photon),
}
READOUT_NOISE = {"NUV": 1.2 * DN_UNIT["NUV"], "FUV": 3.1 * DN_UNIT["FUV"]}


def read_iris_spectrograph_level2_fits(filenames, spectral_windows=None, uncertainty=True, memmap=False):
    """
    Reads IRIS level 2 spectrograph FITS from an OBS into an IRISSpectrograph
    instance.

    Parameters
    ----------
    filenames: `list` of `str` or `str`
        Filename of filenames to be read. They must all be associated with the same
        OBS number.
    spectral_windows: iterable of `str` or `str`
        Spectral windows to extract from files. Default=None, implies, extract all
        spectral windows.

    Returns
    -------
    result: `ndcube.NDCollection`
    """
    if isinstance(filenames, str):
        filenames = [filenames]
    for f, filename in enumerate(filenames):
        hdulist = fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap)
        hdulist.verify("fix")
        if f == 0:
            # Collecting the window observations.
            windows_in_obs = np.array(
                [hdulist[0].header["TDESC{0}".format(i)] for i in range(1, hdulist[0].header["NWIN"] + 1)]
            )
            # If spectral_window is not set then get every window.
            # Else take the appropriate windows
            if not spectral_windows:
                spectral_windows_req = windows_in_obs
                window_fits_indices = range(1, len(hdulist) - 2)
            else:
                if isinstance(spectral_windows, str):
                    spectral_windows_req = [spectral_windows]
                else:
                    spectral_windows_req = spectral_windows
                spectral_windows_req = np.asarray(spectral_windows_req, dtype="U")
                window_is_in_obs = np.asarray([window in windows_in_obs for window in spectral_windows_req])
                if not all(window_is_in_obs):
                    missing_windows = window_is_in_obs == False
                    raise ValueError(
                        "Spectral windows {0} not in file {1}".format(spectral_windows[missing_windows], filenames[0])
                    )
                window_fits_indices = np.nonzero(np.in1d(windows_in_obs, spectral_windows))[0] + 1
            # Create a empty list for every spectral window and each
            # spectral window is a key for the dictionary.
            data_dict = dict([(window_name, list()) for window_name in spectral_windows_req])
        # Extract axis-aligned metadata.
        times = Time(hdulist[0].header["STARTOBS"]) + TimeDelta(
            hdulist[-2].data[:, hdulist[-2].header["TIME"]], format="sec"
        )
        fov_center = SkyCoord(
            Tx=hdulist[-2].data[:, hdulist[-2].header["XCENIX"]],
            Ty=hdulist[-2].data[:, hdulist[-2].header["YCENIX"]],
            unit=u.arcsec,
            frame=Helioprojective,
        )
        obs_vrix = hdulist[-2].data[:, hdulist[-2].header["OBS_VRIX"]] * u.m / u.s
        ophaseix = hdulist[-2].data[:, hdulist[-2].header["OPHASEIX"]]
        exposure_times_fuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEF"]] * u.s
        exposure_times_nuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEN"]] * u.s
        for i, window_name in enumerate(spectral_windows_req):
            # Define metadata object for window.
            meta = IRISSGMeta(
                hdulist[0].header,
                window_name,
                data_shape=hdulist[window_fits_indices[i]].data.shape,
            )
            # Determine values of properties dependent on detector type.
            if "FUV" in meta.detector:
                exposure_times = exposure_times_fuv
                DN_unit = DN_UNIT["FUV"]
                readout_noise = READOUT_NOISE["FUV"]
            else:
                exposure_times = exposure_times_nuv
                DN_unit = DN_UNIT["NUV"]
                readout_noise = READOUT_NOISE["NUV"]
            meta.add("exposure time", exposure_times, None, 0)
            meta.add("exposure FOV center", fov_center, None, 0)
            meta.add("observer radial velocity", obs_vrix, None, 0)
            meta.add("orbital phase", ophaseix, None, 0)
            # Derive WCS, data and mask for NDCube from file.
            # Sit-and-stare have a CDELT of 0 which causes issues in astropy WCS.
            # In this case, set CDELT to a tiny non-zero number.
            if hdulist[window_fits_indices[i]].header["CDELT3"] == 0:
                hdulist[window_fits_indices[i]].header["CDELT3"] = 1e-10
            wcs_ = WCS(hdulist[window_fits_indices[i]].header)
            if not memmap:
                data_mask = hdulist[window_fits_indices[i]].data == -200.0
            else:
                data_mask = None
            # Derive uncertainty of data
            if uncertainty:
                out_uncertainty = (
                    u.Quantity(
                        np.sqrt(
                            (hdulist[window_fits_indices[i]].data * DN_unit).to(u.photon).value
                            + readout_noise.to(u.photon).value ** 2
                        ),
                        unit=u.photon,
                    )
                    .to(DN_unit)
                    .value
                )
            else:
                out_uncertainty = None
            # Appending NDCube instance to the corresponding window key in dictionary's list.
            cube = SpectrogramCube(
                hdulist[window_fits_indices[i]].data,
                wcs=wcs_,
                uncertainty=out_uncertainty,
                unit=DN_unit,
                meta=meta,
                mask=data_mask,
            )
            cube.extra_coords.add("time", 0, times, physical_types="time")
            data_dict[window_name].append(cube)
        hdulist.close()
    # Construct dictionary of SpectrogramSequences for spectral windows
    window_data_pairs = [
        (window_name, RasterSequence(data_dict[window_name], common_axis=0)) for window_name in spectral_windows_req
    ]
    # Initialize an NDCollection object.
    return NDCollection(window_data_pairs, aligned_axes=(0, 1, 2))


class IRISSGMeta(Meta, metaclass=SlitSpectrographMetaABC):
    def __init__(self, header, spectral_window, **kwargs):
        super().__init__(header, **kwargs)
        spectral_windows = np.array([self["TDESC{0}".format(i)] for i in range(1, self["NWIN"] + 1)])
        window_mask = np.array([spectral_window in window for window in spectral_windows])
        if window_mask.sum() < 1:
            raise ValueError(
                "Spectral window not found. "
                f"Input spectral window: {spectral_window}; "
                f"Spectral windows in header: {spectral_windows}"
            )
        elif window_mask.sum() > 1:
            raise ValueError(
                "Spectral window must be unique. "
                f"Input spectral window: {spectral_window}; "
                f"Ambiguous spectral windows in header: {spectral_windows[window_mask]}"
            )
        self._iwin = np.arange(len(spectral_windows))[window_mask][0] + 1

    def __str__(self):
        return textwrap.dedent(
            f"""\
                IRISMeta
                --------
                Observatory:\t\t{self.observatory}
                Instrument:\t\t{self.instrument}
                Detector:\t\t{self.detector}
                Spectral Window:\t{self.spectral_window}
                Spectral Range:\t\t{self.spectral_range}
                Date:\t\t\t{self.date_reference}
                OBS ID:\t\t\t{self.observing_mode_id}
                OBS Description:\t{self.observing_mode_description}
                """
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    # ---------- IRIS-specific convenience methods ----------
    def _construct_time(self, key):
        val = self.get(key, None)
        if val is not None:
            val = Time(val, format="fits", scale="utc")
        return val

    # ---------- Inherited ABC properties ----------
    @property
    def spectral_window(self):
        return self.get(f"TDESC{self._iwin}")

    @property
    def detector(self):
        return self.get(f"TDET{self._iwin}")

    @property
    def instrument(self):
        return self.get("INSTRUME")

    @property
    def observatory(self):
        return self.get("TELESCOP")

    @property
    def processing_level(self):
        return self.get("DATA_LEV")

    @property
    def distance_to_sun(self):
        return self.get("DSUN_OBS") * u.m

    @property
    def date_reference(self):
        return self._construct_time("DATE_OBS")

    @property
    def date_start(self):
        return self.date_reference

    @property
    def date_end(self):
        return self._construct_time("DATE_END")

    @property
    def observing_mode_id(self):
        return int(self.get("OBSID"))

    # ---------- IRIS-specific metadata properties ----------

    @property
    def observing_mode_description(self):
        return self.get("OBS_DESC")

    @property
    def observing_campaign_start(self):
        """
        Start time of observing campaign.
        """
        return self._construct_time("STARTOBS")

    @property
    def observing_campaign_end(self):
        """
        End time of observing mode.
        """
        return self._construct_time("ENDOBS")

    @property
    def observation_includes_SAA(self):
        """
        Whether IRIS passed through SAA during observations.
        """
        return bool(self.get("SAA"))

    @property
    def satellite_rotation(self):
        """
        Satellite roll from solar north.
        """
        return self.get("SAT_ROT") * u.deg

    @property
    def exposure_control_triggers_in_observation(self):
        """
        Number of times automatic exposure control triggered during observing
        campaign.
        """
        return self.get("AECNOBS")

    @property
    def exposure_control_triggers_in_raster(self):
        """
        Number of times automatic exposure control was triggered during this
        raster.
        """
        return self.get("AECNRAS")

    @property
    def number_raster_positions(self):
        """
        Number of positions in raster.
        """
        self.get("NRASTERP")

    @property
    def spectral_range(self):
        """
        The spectral range of the spectral window.
        """
        return [self.get(f"TWMIN{self._iwin}"), self.get(f"TWMAX{self._iwin}")] * u.AA

    @property
    def raster_FOV_width_y(self):
        """
        Width of the field of view of the raster in the Y (slit) direction.
        """
        return self.get("FOVY") * u.arcsec

    @property
    def raster_FOV_width_x(self):
        """
        Width of the field of view of the raster in the X (rastering)
        direction.
        """
        return self.get("FOVX") * u.arcsec

    @property
    def FOV_center(self):
        """
        Location of the center of the field of view.
        """
        return SkyCoord(
            Tx=self.get("XCEN"),
            Ty=self.get("YCEN"),
            unit=u.arcsec,
            frame=Helioprojective,
        )

    @property
    def automatic_exposure_control_enabled(self):
        return bool(self.get("IAECFLAG"))

    @property
    def tracking_mode_enabled(self):
        return bool(self.get("TR_MODE"))

    @property
    def observatory_at_high_latitude(self):
        """
        Whether IRIS passed through high Earth latitude during observations.
        """
        return bool(self.get("HLZ"))

    @property
    def spatial_summing_factor(self):
        """
        Number of pixels summed together in the spatial (Y/slit) direction.
        """
        return self.get("SUMSPAT")

    @property
    def spectral_summing_factor(self):
        """
        Number of pixels summed together in the spectral direction.
        """
        if "fuv" in self.detector.lower():
            return self.get("SUMSPTRF")
        else:
            return self.get("SUMSPTRN")
