import abc

from sunraster.extern.meta import Meta

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
        """
        The level to which the data has been processed.
        """

    @abc.abstractproperty
    def observer_location(self):
        """
        Coordinate of observatory location based on header info.
        """

    @abc.abstractproperty
    def date_reference(self):
        """
        The base time from which time axis values are measured.

        Often the same or very similar to date_start.
        """

    @abc.abstractproperty
    def date_start(self):
        pass

    @abc.abstractproperty
    def date_end(self):
        pass

    @abc.abstractproperty
    def version(self):
        """
        The data version.
        """


class RemoteSensorMetaABC(MetaABC):
    @abc.abstractproperty
    def rsun_meters(self):
        """
        Solar radius in units of length.
        """

    @abc.abstractproperty
    def rsun_angular(self):
        """
        Solar radius in angular units as seen from observatory.
        """

    @abc.abstractproperty
    def distance_to_sun(self):
        """
        Distance to Sun center from observatory.
        """


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @abc.abstractproperty
    def spectral_window(self):
        pass

    @abc.abstractproperty
    def observing_mode_id(self):
        """
        Unique identifier for the observing mode.

        Often referred to as OBS ID.
        """

    @abc.abstractproperty
    def observer_radial_velocity(self):
        """
        Velocity of observatory in direction of source.
        """
