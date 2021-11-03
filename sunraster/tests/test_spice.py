import os.path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from ndcube import NDCollection
from sunpy.coordinates import HeliographicStonyhurst

from sunraster import RasterSequence, SpectrogramCube, SpectrogramSequence
from sunraster.instr.spice import SPICEMeta, read_spice_l2_fits
from sunraster.tests import test_data_dir

READ_SPICE_L2_FITS_RETURN_TYPE = NDCollection
SPECTRAL_WINDOW = ("WINDOW0_74.73", "Extension name")
DETECTOR = ("SW", "Detector array name")
INSTRUMENT = ("SPICE", "Instrument name")
OBSERVATORY = ("Solar Orbiter", "Observatory Name")
PROCESSING_LEVEL = ("L2", "Data processing level")
RSUN_METERS = (695700000.0, "[m]      Assumed  photospheric Solar radius")
RSUN_ANGULAR = (1764.0728936, "[arcsec] Apparent photospheric Solar radius")
OBSERVING_MODE_ID = (10, "")
OBSERVATORY_RADIAL_VELOCITY = (
    -7036.06122832,
    "[m/s] Radial velocity of S/C away from the Sun",
)
DISTANCE_TO_SUN = (81342963151.0, "[m]  S/C distance from Sun")
DATE_REFERENCE = ("2020-06-02T07:47:58.017", "[UTC] Equals DATE-BEG")
DATE_START = ("2020-06-02T07:47:58.017", "[UTC] Beginning of data acquisition")
DATE_END = ("2020-06-02T07:47:58.117", "[UTC] End of data acquisition")
HGLN_OBS = (35.8382263864, "[deg] S/C Heliographic longitude")
HGLT_OBS = (4.83881036748, "[deg] S/C Heliographic latitude (B0 angle)")
SPICE_OBSERVING_MODE_ID = (12583744, "SPICE Observation ID")
DARKMAP = (0, "If set, a dark map was subtracted on-board")
BLACKLEV = (0, "If set, a bias frame was subtracted on-board")
WINDOW_TYPE = ("Full Detector Narrow-slit", "Description of window type")
WINDOW_TABLE_ID = (255, "Index in on-board window data table (0-255)")
SLIT_ID = (2, "Slit ID (0-3)")
SLIT_WIDTH = (4, "[arcsec] Slit width")
DUMBBELL = (0, "0/1/2: not a dumbbell/lower dumbbel/upper dumbb")
SOLAR_B0 = (4.83881036748, "[deg] Tilt angle of Solar North toward S/C")
SOLAR_P0 = (1.49702480927, "[deg] S/C Celestial North to Solar North angle")
SOLAR_EP = (-6.14143491727, "[deg] S/C Ecliptic  North to Solar North angle")
CARRINGTON_ROTATION_NUMBER = (2231, "Carrington rotation number")
DATE_START_EARTH = ("2020-06-02T07:51:52.799", "[UTC] DATE-BEG + EAR_TDEL")
DATE_START_SUN = ("2020-06-02T07:43:26.686", "[UTC] DATE-BEG - SUN_TIME")


@pytest.fixture
def spice_fits_header():
    hdr = fits.Header()
    hdr.append(tuple(["EXTNAME"] + list(SPECTRAL_WINDOW)))
    hdr.append(tuple(["DETECTOR"] + list(DETECTOR)))
    hdr.append(tuple(["INSTRUME"] + list(INSTRUMENT)))
    hdr.append(tuple(["OBSRVTRY"] + list(OBSERVATORY)))
    hdr.append(tuple(["LEVEL"] + list(PROCESSING_LEVEL)))
    hdr.append(tuple(["RSUN_REF"] + list(RSUN_METERS)))
    hdr.append(tuple(["RSUN_ARC"] + list(RSUN_ANGULAR)))
    hdr.append(tuple(["OBS_ID"] + list(OBSERVING_MODE_ID)))
    hdr.append(tuple(["OBS_VR"] + list(OBSERVATORY_RADIAL_VELOCITY)))
    hdr.append(tuple(["DSUN_OBS"] + list(DISTANCE_TO_SUN)))
    hdr.append(tuple(["DATE-OBS"] + list(DATE_REFERENCE)))
    hdr.append(tuple(["DATE-BEG"] + list(DATE_START)))
    hdr.append(tuple(["DATE-END"] + list(DATE_END)))
    hdr.append(tuple(["HGLN_OBS"] + list(HGLN_OBS)))
    hdr.append(tuple(["HGLT_OBS"] + list(HGLT_OBS)))
    hdr.append(tuple(["SPIOBSID"] + list(SPICE_OBSERVING_MODE_ID)))
    hdr.append(tuple(["DARKMAP"] + list(DARKMAP)))
    hdr.append(tuple(["BLACKLEV"] + list(BLACKLEV)))
    hdr.append(tuple(["WIN_TYPE"] + list(WINDOW_TYPE)))
    hdr.append(tuple(["WINTABID"] + list(WINDOW_TABLE_ID)))
    hdr.append(tuple(["SLIT_ID"] + list(SLIT_ID)))
    hdr.append(tuple(["SLIT_WID"] + list(SLIT_WIDTH)))
    hdr.append(tuple(["DUMBBELL"] + list(DUMBBELL)))
    hdr.append(tuple(["SOLAR_B0"] + list(SOLAR_B0)))
    hdr.append(tuple(["SOLAR_P0"] + list(SOLAR_P0)))
    hdr.append(tuple(["SOLAR_EP"] + list(SOLAR_EP)))
    hdr.append(tuple(["CAR_ROT"] + list(CARRINGTON_ROTATION_NUMBER)))
    hdr.append(tuple(["DATE_EAR"] + list(DATE_START_EARTH)))
    hdr.append(tuple(["DATE_SUN"] + list(DATE_START_SUN)))
    return hdr


@pytest.fixture
def spice_meta(spice_fits_header):
    return SPICEMeta(
        spice_fits_header,
        comments=zip(spice_fits_header.keys(), spice_fits_header.comments),
    )


@pytest.fixture
def spice_rasdb_filename(tmp_path):
    """
    Inserts data into a raster SPICE FITS file with dumbbells and returns new
    filename.

    A new FITS file is saved in a tmp file path.
    """
    filename = "solo_L2_spice-n-ras-db_20200602T081733_V01_12583760-000.fits"
    with fits.open(os.path.join(test_data_dir, filename)) as hdulist:
        new_hdulist = fits.HDUList()
        new_hdulist.append(fits.PrimaryHDU(np.random.rand(1, 48, 832, 30), header=hdulist[0].header))
        new_hdulist.append(fits.ImageHDU(np.random.rand(1, 48, 832, 30), header=hdulist[1].header))
        new_hdulist.append(fits.ImageHDU(np.random.rand(1, 56, 64, 30), header=hdulist[2].header))
        new_hdulist.append(fits.ImageHDU(np.random.rand(1, 56, 64, 30), header=hdulist[3].header))
        new_hdulist.append(hdulist[-1])
        tmp_spice_path = tmp_path / "spice"
        if not os.path.exists(tmp_spice_path):
            tmp_spice_path.mkdir()
        new_filename = os.path.join(tmp_spice_path, filename)
        new_hdulist.writeto(new_filename, overwrite=True)
    return new_filename


@pytest.fixture
def spice_sns_filename(tmp_path):
    """
    Inserts data into a sit-and-stare SPICE FITS file and returns new filename.

    A new FITS file is saved in a tmp file path.
    """
    filename = "solo_L2_spice-n-sit_20200620T235901_V01_16777431-000.fits"
    with fits.open(os.path.join(test_data_dir, filename)) as hdulist:
        new_hdulist = fits.HDUList()
        new_hdulist.append(fits.PrimaryHDU(np.random.rand(32, 48, 1024, 1), header=hdulist[0].header))
        new_hdulist.append(fits.ImageHDU(np.random.rand(32, 48, 1024, 1), header=hdulist[1].header))
        new_hdulist.append(hdulist[-1])
        tmp_spice_path = tmp_path / "spice"
        if not os.path.exists(tmp_spice_path):
            tmp_spice_path.mkdir()
        new_filename = os.path.join(tmp_spice_path, filename)
        new_hdulist.writeto(new_filename, output_verify="fix+ignore", overwrite=True)
    return new_filename


def _construct_expected_time(date_info):
    return Time(date_info[0], format="fits", scale=date_info[1][1:4].lower())


def test_meta_spectral_window(spice_meta):
    assert spice_meta.spectral_window == SPECTRAL_WINDOW[0][8:]


def test_meta_detector(spice_meta):
    assert spice_meta.detector == DETECTOR[0]


def test_meta_instrument(spice_meta):
    assert spice_meta.instrument == INSTRUMENT[0]


def test_meta_observatory(spice_meta):
    assert spice_meta.observatory == OBSERVATORY[0]


def test_meta_processing_level(spice_meta):
    assert spice_meta.processing_level == PROCESSING_LEVEL[0]


def test_meta_rsun_meters(spice_meta):
    assert spice_meta.rsun_meters == RSUN_METERS[0] * u.m


def test_meta_rsun_angular(spice_meta):
    assert spice_meta.rsun_angular == RSUN_ANGULAR[0] * u.arcsec


def test_meta_spice_observation_id(spice_meta):
    assert spice_meta.spice_observation_id == SPICE_OBSERVING_MODE_ID[0]


def test_meta_observer_radial_velocity(spice_meta):
    assert spice_meta.observer_radial_velocity == OBSERVATORY_RADIAL_VELOCITY[0] * u.m / u.s


def test_meta_distance_to_sun(spice_meta):
    assert spice_meta.distance_to_sun == DISTANCE_TO_SUN[0] * u.m


def test_meta_date_reference(spice_meta):
    assert spice_meta.date_reference == _construct_expected_time(DATE_REFERENCE)


def test_meta_date_start(spice_meta):
    assert spice_meta.date_start == _construct_expected_time(DATE_START)


def test_meta_date_end(spice_meta):
    assert spice_meta.date_end == _construct_expected_time(DATE_END)


def test_meta_observer_location(spice_meta):
    obstime = _construct_expected_time(DATE_REFERENCE)
    observer_location = SkyCoord(
        lon=HGLN_OBS[0],
        lat=HGLT_OBS[0],
        radius=DISTANCE_TO_SUN[0],
        unit=(u.deg, u.deg, u.m),
        obstime=obstime,
        frame=HeliographicStonyhurst,
    )
    assert spice_meta.observer_location.lon == observer_location.lon
    assert spice_meta.observer_location.lat == observer_location.lat
    assert spice_meta.observer_location.radius == observer_location.radius
    assert spice_meta.observer_location.obstime == observer_location.obstime
    assert spice_meta.observer_location.frame.name == observer_location.frame.name


def test_meta_observing_mode_id_solar_orbiter(spice_meta):
    assert spice_meta.observing_mode_id_solar_orbiter == OBSERVING_MODE_ID[0]


def test_meta_darkmap_subtracted_onboard(spice_meta):
    assert spice_meta.darkmap_subtracted_onboard == False


def test_meta_bias_frame_subtracted_onboard(spice_meta):
    assert spice_meta.bias_frame_subtracted_onboard == False


def test_meta_window_type(spice_meta):
    assert spice_meta.window_type == WINDOW_TYPE[0]


def test_meta_slit_id(spice_meta):
    assert spice_meta.slit_id == SLIT_ID[0]


def test_meta_slit_width(spice_meta):
    assert spice_meta.slit_width == SLIT_WIDTH[0] * u.arcsec


def test_meta_contains_dumbbell(spice_meta):
    assert spice_meta.contains_dumbbell is False


def test_meta_dumbbell_type(spice_meta):
    assert spice_meta.dumbbell_type is None


def test_meta_solar_B0(spice_meta):
    assert spice_meta.solar_B0 == SOLAR_B0[0] * u.deg


def test_meta_solar_P0(spice_meta):
    assert spice_meta.solar_P0 == SOLAR_P0[0] * u.deg


def test_meta_solar_ep(spice_meta):
    assert spice_meta.solar_ep == SOLAR_EP[0] * u.deg


def test_meta_carrington_rotation(spice_meta):
    assert spice_meta.carrington_rotation == CARRINGTON_ROTATION_NUMBER[0]


def test_meta_date_start_earth(spice_meta):
    date_start_earth = _construct_expected_time(DATE_START_EARTH)
    assert spice_meta.date_start_earth == date_start_earth


def test_meta_date_start_sun(spice_meta):
    date_start_sun = _construct_expected_time(DATE_START_SUN)
    assert spice_meta.date_start_sun == date_start_sun


def test_read_spice_l2_fits_single_file_multiple_windows(spice_rasdb_filename):
    filename = spice_rasdb_filename
    result = read_spice_l2_fits(filename)
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert set(result.aligned_axes.values()) == {(0, 2, 3)}
    assert len(result) == 2
    assert all(isinstance(window, SpectrogramCube) for window in result.values())


def test_read_spice_l2_fits_single_file_window(spice_rasdb_filename):
    filename = spice_rasdb_filename
    result = read_spice_l2_fits(filename, windows=["WINDOW0_70.51"])
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert result.aligned_axes is None
    assert len(result) == 1
    assert all(isinstance(window, SpectrogramCube) for window in result.values())


def test_read_spice_l2_fits_single_file_dumbbells(spice_rasdb_filename):
    filename = spice_rasdb_filename
    result = read_spice_l2_fits(filename, read_dumbbells=True)
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert all(window.meta.contains_dumbbell for window in result.values())
    assert set(result.aligned_axes.values()) == {tuple(range(4))}
    assert all(isinstance(window, SpectrogramCube) for window in result.values())


def test_read_spice_l2_fits_multiple_rasters_multiple_windows(spice_rasdb_filename):
    filenames = [spice_rasdb_filename] * 2
    result = read_spice_l2_fits(filenames)
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert set(result.aligned_axes.values()) == {(0, 2, 3)}
    assert len(result) == 2
    assert all([window.dimensions[0].value == len(filenames) for window in result.values()])
    assert all([isinstance(window, RasterSequence) for window in result.values()])


def test_read_spice_l2_fits_multiple_rasters_single_window(spice_rasdb_filename):
    filenames = [spice_rasdb_filename] * 2
    result = read_spice_l2_fits(filenames, windows=["WINDOW0_70.51"])
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert result.aligned_axes is None
    assert len(result) == 1
    assert all(window.dimensions[0].value == len(filenames) for window in result.values())
    assert all(isinstance(window, RasterSequence) for window in result.values())


def test_read_spice_l2_fits_multiple_sns_multiple_windows(spice_sns_filename):
    filenames = [spice_sns_filename] * 2
    result = read_spice_l2_fits(filenames)
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert set(result.aligned_axes.values()) == {(0, 2, 3)}
    assert len(result) == 2
    assert all(window.dimensions[0].value == len(filenames) for window in result.values())
    assert all(isinstance(window, SpectrogramSequence) for window in result.values())


def test_read_spice_l2_fits_multiple_files_dumbbells(spice_rasdb_filename):
    filenames = [spice_rasdb_filename] * 2
    result = read_spice_l2_fits(filenames, read_dumbbells=True)
    assert isinstance(result, READ_SPICE_L2_FITS_RETURN_TYPE)
    assert all(window[0].meta.contains_dumbbell for window in result.values())
    assert set(result.aligned_axes.values()) == {tuple(range(4))}
    assert all(window.dimensions[0].value == len(filenames) for window in result.values())
    assert all(isinstance(window, SpectrogramSequence) for window in result.values())


def test_read_spice_l2_fits_incompatible_files(spice_rasdb_filename, spice_sns_filename):
    with pytest.raises(ValueError):
        filenames = [spice_rasdb_filename, spice_sns_filename]
        read_spice_l2_fits(filenames)
