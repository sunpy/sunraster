
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


def test_meta_get_fits_default(meta_fits):
    assert meta_fits.get(MISSING_KEY, "None") == "None"


def test_meta_get_fits_comment(meta_fits):
    assert meta_fits.get_comment(KEY) == COMMENT


@pytest.mark.parametrize("key,new_value", [(KEY, 1),
                                           (MISSING_KEY, 1)])
def test_meta_update_fits(meta_fits, key, new_value):
    meta_fits.update([(key, new_value)])
    assert meta_fits.get(key) == new_value


def test_meta_update_fits_history(meta_fits):
    second_entry = "2nd entry."
    meta_fits.update([(ADDITIVE_KEY, second_entry)])
    assert list(meta_fits.get(ADDITIVE_KEY)) == [ADDITIVE_ENTRY, second_entry]


def test_meta_update_comments(meta_fits):
    new_comment = "This is a new comment."
    meta_fits.update_comments([(KEY, new_comment)])
    assert meta_fits.get_comment(KEY) == new_comment
