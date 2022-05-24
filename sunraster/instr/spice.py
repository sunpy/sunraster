import copy
import numbers
import textwrap

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from ndcube import NDCollection

from sunraster import RasterSequence, SpectrogramCube, SpectrogramSequence
from sunraster.meta import Meta, SlitSpectrographMetaABC

__all__ = ["read_spice_l2_fits", "SPICEMeta"]


INCORRECT_OBSID_MESSAGE = "File has incorrect SPIOBSID."


def read_spice_l2_fits(filenames, windows=None, memmap=True, read_dumbbells=False):
    """
    Read SPICE level 2 FITS file.

    Parameters
    ----------
    filenames: iterable of `str`
        The name(s), including path, of the SPICE FITS file(s) to read.
    windows: iterable of `str`
        The names of the windows to read.
        All windows must of the same type: dumbbell and regular.
        Default=None implies all narrow-slit or dumbbell windows read out
        depending on value of read_dumbells kwarg. See below.
    memmap: `bool`
        If True, FITS file is reading with memory mapping.
    read_dumbbells: `bool`
        Defines whether dumbbell or regular windows are returned.
        If True, returns the dumbbell windows.
        If False, returns regular windows.
        Default=False
        Ignored if windows kwarg is set.

    Returns
    -------
    output: `ndcube.NDCollection` or `sunraster.SpectrogramCube`, `sunraster.RasterSequence`,
        `sunraster.SpectrogramSequence`
        A collection of spectrogram/raster cubes/sequences, one for each window.
        If only one window present or requested, a single spectrogram cube
        or sequence is returned.
    """
    # Sanitize inputs.
    if isinstance(filenames, str):
        filenames = [filenames]
    # Read first file.
    first_cubes = _read_single_spice_l2_fits(
        filenames[0], windows=windows, memmap=memmap, read_dumbbells=read_dumbbells
    )
    # Derive information for consistency checks between files and read subsequent files.
    if len(filenames) > 1:
        # Wrap windows from first file in lists
        # so windows from other files can be appended.
        cube_lists = dict([(key, [value]) for key, value in first_cubes.items()])
        # Get info from first file for consistency checks between files.
        first_meta = _get_meta_from_last_added(cube_lists)
        first_obs_id = _get_obsid(first_meta)
        if windows is None:
            windows = list(cube_lists.keys())
        # Read subsequent files and append output to relevant window in cube_lists.
        for i, filename in enumerate(filenames[1:]):
            try:
                cube_lists = _read_single_spice_l2_fits(
                    filename,
                    windows=windows,
                    memmap=memmap,
                    read_dumbbells=read_dumbbells,
                    output=cube_lists,
                    spice_id=first_obs_id,
                )
            except ValueError as err:
                err_message = err.args[0]
                if INCORRECT_OBSID_MESSAGE in err_message:
                    this_obs_id = err_message.split()[-1]
                    raise ValueError(
                        "All files must correspond to same observing campaign/SPICE OBS ID. "
                        f"First file SPICE OBS ID: {first_obs_id}; "
                        f"{i+1}th file SPICE OBS ID: {this_obs_id}"
                    )
        # Depending on type of file, combine data from different files into
        # SpectrogramSequences and RasterSequences.
        is_raster = "ras" in first_meta.get("FILENAME") and not any(
            [window[1].meta.contains_dumbbell for window in cube_lists.values()]
        )
        if is_raster:
            sequence_class = RasterSequence
        else:
            sequence_class = SpectrogramSequence
        window_sequences = [
            (key, sequence_class([v[0] for v in value], common_axis=-1)) for key, value in cube_lists.items()
        ]
    else:
        # If only one file being read, leave data in SpectrogramCube objects.
        window_sequences = list(first_cubes.items())
    if len(window_sequences) > 1:
        # Data should be aligned along all axes except the spectral axis.
        # But they should be aligned along all axes if they come from the
        # same spectral window, e.g. because they are dumbbell windows.
        first_sequence = window_sequences[0][1]
        first_spectral_window = first_sequence[0].meta.spectral_window
        if all([window[1][0].meta.spectral_window == first_spectral_window for window in window_sequences]):
            aligned_axes = tuple(range(len(first_sequence.dimensions)))
        else:
            aligned_axes = tuple(
                i for i, phys_type in enumerate(first_sequence.array_axis_physical_types) if phys_type != ("em.wl",)
            )
    else:
        aligned_axes = None
    return NDCollection(window_sequences, aligned_axes=aligned_axes)


def _get_meta_from_last_added(obj):
    return list(obj.values())[0][-1].meta


def _get_obsid(spice_meta):
    return spice_meta.spice_observation_id


def _read_single_spice_l2_fits(
    filename,
    windows=None,
    memmap=True,
    read_dumbbells=False,
    output=None,
    spice_id=None,
):
    """
    Read SPICE level 2 FITS file(s).

    Parameters
    ----------
    filename: `str`
        The name, including path, of the SPICE FITS file to read.
    windows: iterable of `str`
        The names of the windows to read.
        All windows must of the same type: dumbbell and regular.
        Default=None implies all narrow-slit or dumbbell windows read out
        depending on value of read_dumbells kwarg. See below.
    memmap: `bool`
        If True, FITS file is reading with memory mapping.
    read_dumbbells: `bool`
        Defines whether dumbbell or regular windows are returned.
        If True, returns the dumbbell windows.
        If False, returns regular windows.
        Default=False
        Ignored if windows kwarg is set.
    output: `dict` of `list`s (optional)
        A dictionary of lists with the same keys are the windows kwarg.
        The output for each window will be appended to the list corresponding
        the window's name.
    spice_id: `int` (optional)
        If not None, file must have a SPIOBSID equal to this value.
        Otherwise an error is raised

    Returns
    -------
    output: `dict` of `sunraster.SpectrogramCube`
        A collection of spectrogram cubes, one for each window.
    """
    window_cubes = []
    dumbbell_label = "DUMBBELL"
    with fits.open(filename, memmap=memmap) as hdulist:
        if isinstance(spice_id, numbers.Integral) and hdulist[0].header["SPIOBSID"] != spice_id:
            raise ValueError(f"{INCORRECT_OBSID_MESSAGE} Expected {spice_id}. Got {hdulist[0].header['SPIOBSID']}.")
        # Derive names of windows to be read.
        if windows is None:
            if read_dumbbells:
                windows = [hdu.header["EXTNAME"] for hdu in hdulist if dumbbell_label in hdu.header["EXTNAME"]]
            else:
                windows = [
                    hdu.header["EXTNAME"]
                    for hdu in hdulist
                    if (hdu.header["EXTNAME"] != "VARIABLE_KEYWORDS" and dumbbell_label not in hdu.header["EXTNAME"])
                ]
        dumbbells_requested = [dumbbell_label in window for window in windows]
        if any(dumbbells_requested) and not all(dumbbells_requested):
            raise ValueError("Cannot read dumbbell and other window types simultaneously.")
        # Retrieve window names from FITS file.
        for i, hdu in enumerate(hdulist):
            if hdu.header["EXTNAME"] in windows:
                # Define metadata object.
                meta = SPICEMeta(
                    hdu.header,
                    comments=_convert_fits_comments_to_key_value_pairs(hdu.header),
                    data_shape=hdu.data.shape,
                )
                # Rename WCS time axis to time.
                meta.update([("CTYPE4", "TIME")])
                new_header = copy.deepcopy(hdu.header)
                new_header["CTYPE4"] = "TIME"
                # Define WCS from new header
                wcs = WCS(new_header)
                # Define exposure times from metadata.
                exp_times = u.Quantity(np.zeros(hdu.data.shape[-1]) + meta.get("XPOSURE"), unit=u.s)
                # Define data cube.
                data = hdu.data
                spectrogram = SpectrogramCube(
                    data=data,
                    wcs=wcs,
                    mask=np.isnan(data),
                    unit=u.adu,
                    meta=meta,
                    instrument_axes=("raster scan", "spectral", "slit", "slit step"),
                )
                spectrogram.meta.add("exposure time", exp_times, None, 3)
                window_name = meta.get("EXTNAME")
                if output is None:
                    window_cubes.append((window_name, spectrogram))
                else:
                    output[window_name].append(spectrogram)
    if output is None:
        return dict(window_cubes)
    else:
        return output


def _convert_fits_comments_to_key_value_pairs(fits_header):
    keys = np.unique(np.array(list(fits_header.keys())))
    keys = keys[keys != ""]
    return [(key, fits_header.comments[key]) for key in keys]


class SPICEMeta(Meta, metaclass=SlitSpectrographMetaABC):
    # ---------- SPICE-specific convenience methods ----------
    def _get_unit(self, key):
        comment = self.comments.get(key)
        if comment:
            try:
                return [s.split("]") for s in comment.split("[")[1:]][0][:-1][0]
            except IndexError:
                pass
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
        return textwrap.dedent(
            f"""\
                SPICEMeta
                ---------
                Observatory:\t\t\t\t{self.observatory}
                Instrument:\t\t\t\t{self.instrument}
                Detector:\t\t\t\t{self.detector}
                Spectral Window:\t\t\t{self.spectral_window}
                Date:\t\t\t\t\t{self.date_reference}
                OBS_ID (SOC Observation ID):\t\t{self.observing_mode_id_solar_orbiter}
                SPIOBSID (SPICE Observation ID):\t{self.spice_observation_id}
                """
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    # ---------- Inherited ABC properties ----------
    @property
    def spectral_window(self):
        spectral_window = self.get("EXTNAME")
        # Remove redundant text associated with dumbbells.
        joiner = "_"
        if self.contains_dumbbell:
            dummy_txt = ""
            spectral_window = spectral_window.replace("DUMBBELL", dummy_txt)
            spectral_window = spectral_window.replace("UPPER", dummy_txt)
            spectral_window = spectral_window.replace("LOWER", dummy_txt)
            spectral_window = joiner.join(list(filter((dummy_txt).__ne__, spectral_window.split(joiner))))
        # Remove other redundant text from window name.
        redundant_txt = "WINDOW"
        if redundant_txt in spectral_window:
            spectral_window = joiner.join([comp for comp in spectral_window.split(joiner) if "WINDOW" not in comp])
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
    def spice_observation_id(self):
        return self.get("SPIOBSID")

    @property
    def observer_radial_velocity(self):
        return self._construct_quantity("OBS_VR")

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
    def observer_location(self):
        from sunpy.coordinates import HeliographicStonyhurst

        lon_unit = u.deg
        lat_unit = u.deg
        radius_unit = u.m
        lon_key = "HGLN_OBS"
        lat_key = "HGLT_OBS"
        kwargs = {
            "lon": u.Quantity(self.get(lon_key), unit=self._get_unit(lon_key)).to_value(lon_unit),
            "lat": u.Quantity(self.get(lat_key), unit=self._get_unit(lat_key)).to_value(lat_unit),
            "radius": self.distance_to_sun.to_value(radius_unit),
            "unit": (lon_unit, lat_unit, radius_unit),
            "frame": HeliographicStonyhurst,
        }
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
    def slit_id(self):
        return self.get("SLIT_ID")

    @property
    def slit_width(self):
        return self._construct_quantity("SLIT_WID")

    @property
    def contains_dumbbell(self):
        return self.get("DUMBBELL") in [1, 2]

    @property
    def dumbbell_type(self):
        dumbbell_types = [None, "lower", "upper"]
        dumbbell_idx = self.get("DUMBBELL")
        return dumbbell_types[dumbbell_idx]

    @property
    def solar_B0(self):
        """
        Tilt angle of solar north toward spacecraft.
        """
        return self._construct_quantity("SOLAR_B0")

    @property
    def solar_P0(self):
        """
        Angle from spacecraft celestial north to solar north.
        """
        return self._construct_quantity("SOLAR_P0")

    @property
    def solar_ep(self):
        """
        Angle from spacecraft ecliptic north to solar north angle.
        """
        return self._construct_quantity("SOLAR_EP")

    @property
    def carrington_rotation(self):
        """
        Carrington Rotation number of observation.
        """
        return self.get("CAR_ROT")

    @property
    def date_start_earth(self):
        """
        Time at which photons reaching SPICE at start time would have reach
        Earth.
        """
        return self._construct_time("DATE_EAR")

    @property
    def date_start_sun(self):
        """
        Time at which photons reaching SPICE at start time would have left Sun.

        The Sun is defined as the center of the Sun assuming photon was
        not impeded.
        """
        return self._construct_time("DATE_SUN")
