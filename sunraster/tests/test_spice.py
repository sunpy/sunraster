import pytest
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from sunpy.coordinates import HeliographicStonyhurst

from sunraster.instr.spice import SPICEMeta


SPECTRAL_WINDOW = ('WINDOW0_74.73', 'Extension name')
DETECTOR = ('SW', 'Detector array name')
INSTRUMENT = ('SPICE', 'Instrument name')
OBSERVATORY = ('Solar Orbiter', 'Observatory Name')
PROCESSING_LEVEL = ('L2', 'Data processing level')
RSUN_METERS = (695700000.0, '[m]      Assumed  photospheric Solar radius')
RSUN_ANGULAR = (1764.0728936, '[arcsec] Apparent photospheric Solar radius')
OBSERVING_MODE_ID = (10, '')
OBSERVATORY_RADIAL_VELOCITY = (-7036.06122832, '[m/s] Radial velocity of S/C away from the Sun')
DISTANCE_TO_SUN = (81342963151.0, '[m]  S/C distance from Sun')
DATE_REFERENCE = ('2020-06-02T07:47:58.017', '[UTC] Equals DATE-BEG')
DATE_START = ('2020-06-02T07:47:58.017', '[UTC] Beginning of data acquisition')
DATE_END = ('2020-06-02T07:47:58.117', '[UTC] End of data acquisition')
HGLN_OBS = (35.8382263864, '[deg] S/C Heliographic longitude')
HGLT_OBS = (4.83881036748, '[deg] S/C Heliographic latitude (B0 angle)')


@pytest.fixture
def spice_fits_header():
    hdr = fits.Header()
    hdr["EXTNAME"] = SPECTRAL_WINDOW
    hdr["DETECTOR"] = DETECTOR
    hdr["INSTRUME"] = INSTRUMENT
    hdr["OBSRVTRY"] = OBSERVATORY
    hdr["LEVEL"] = PROCESSING_LEVEL
    hdr["RSUN_REF"] = RSUN_METERS
    hdr["RSUN_ARC"] = RSUN_ANGULAR
    hdr["OBS_ID"] = OBSERVING_MODE_ID
    hdr["OBS_VR"] = OBSERVATORY_RADIAL_VELOCITY
    hdr["DSUN_OBS"] = DISTANCE_TO_SUN
    hdr["DATE-OBS"] = DATE_REFERENCE
    hdr["DATE-BEG"] = DATE_START
    hdr["DATE-END"] = DATE_END
    hdr["HGLN_OBS"] = HGLN_OBS
    hdr["HGLT_OBS"] = HGLT_OBS
    return hdr


@pytest.fixture
def spice_meta(spice_fits_header):
    return SPICEMeta(spice_fits_header)


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


def test_meta_observing_mode_id(spice_meta):
    assert spice_meta.observing_mode_id == OBSERVING_MODE_ID[0]


def test_meta_observatory_radial_velocity(spice_meta):
    assert spice_meta.observatory_radial_velocity == OBSERVATORY_RADIAL_VELOCITY[0]


def test_meta_distance_to_sun(spice_meta):
    assert spice_meta.distance_to_sun == DISTANCE_TO_SUN[0] * u.m


def test_meta_date_reference(spice_meta):
    assert spice_meta.date_reference == DATE_REFERENCE[0]


def test_meta_date_start(spice_meta):
    assert spice_meta.date_start == DATE_START[0]


def test_meta_date_end(spice_meta):
    assert spice_meta.date_end == DATE_END[0]


def test_meta_observer_coordinate(spice_meta):
    obstime = Time(DATE_REFERENCE[0], format="fits", scale="utc")
    observer_coordinate = SkyCoord(
            lon=HGLN_OBS[0]*u.deg, lat=HGLT_OBS[0]*u.deg, radius=DISTANCE_TO_SUN[0],
            unit=(u.deg, u.deg, u.m), obstime=obstime, frame=HeliographicStonyhurst)
    assert spice_meta.observer_coordinate.lon == observer_coordinate.lon
    assert spice_meta.observer_coordinate.lat == observer_coordinate.lat
    assert spice_meta.observer_coordinate.radius == observer_coordinate.radius
    assert spice_meta.observer_coordinate.obstime == observer_coordinate.obstime
    assert spice_meta.observer_coordinate.frame.name == observer_coordinate.frame.name
