import abc
import copy

from astropy.io import fits

__all__ = ["Meta", "RemoteSensorMetaABC", "SlitSpectrographMetaABC"]

class MetaABC(abc.ABCMeta):
    @abc.abstractproperty
    def detector(self):
        pass

    @abc.abstractproperty
    def instrument(self):
        pass

    @abc.abstractproperty
    def observatory(self):
        pass

    @abc.abstractproperty
    def processing_level(self):
        """The level to which the data has been processed."""
        pass

    @abc.abstractproperty
    def observer_coordinate(self):
        """Coordinate of observatory location based on header info."""
        pass

    @abc.abstractproperty
    def date_reference(self):
        """The base time from which time axis values are measured.

        Often the same or very similar to date_start.
        """
        pass

    @abc.abstractproperty
    def date_start(self):
        pass

    @abc.abstractproperty
    def date_end(self):
        pass

    @abc.abstractproperty
    def version(self):
        """The data version."""
        pass


class RemoteSensorMetaABC(MetaABC):
    @abc.abstractproperty
    def rsun_meters(self):
        """Solar radius in units of length."""
        pass

    @abc.abstractproperty
    def rsun_angular(self):
        """Solar radius in angular units as seen from observatory."""
        pass

    @abc.abstractproperty
    def distance_to_sun(self):
        """Distance to Sun center from observatory."""
        pass


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @abc.abstractproperty
    def spectral_window(self):
        pass

    @abc.abstractproperty
    def observing_mode_id(self):
        """Unique identifier for the observing mode. Often referred to as OBS ID."""
        pass

    @abc.abstractproperty
    def observatory_radial_velocity(self):
        """Velocity of observatory in direction of source."""
        pass


class Meta(dict):
    def __init__(self, header, comments=None):
        super().__init__(header)
        self.original_header = header
        if comments is None:
            comments = {}
        else:
            comments = dict(comments)
        self.comments = comments
