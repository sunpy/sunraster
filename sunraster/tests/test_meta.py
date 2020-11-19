
import pytest
from astropy.io import fits

from sunraster.meta import Meta


KEY = "KEY"
VALUE = 0
COMMENT = "This is a comment."
MISSING_KEY = "NO KEY"
ADDITIVE_KEY = "HISTORY"
ADDITIVE_ENTRY = "1st entry."

@pytest.fixture
def fits_header():
    hdr = fits.Header()
    hdr[KEY] = (VALUE, COMMENT)
    hdr[ADDITIVE_KEY] = ADDITIVE_ENTRY
    return hdr


@pytest.fixture
def meta_fits(fits_header):
    return Meta(fits_header)


def test_meta_get_fits(meta_fits):
    assert meta_fits.get(KEY) == VALUE
