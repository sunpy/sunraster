import abc

from ndcube.meta import NDMetaABC

__all__ = ["MetaABC", "RemoteSensorMetaABC", "SlitSpectrographMetaABC"]


class MetaABC(NDMetaABC):
    @property
    @abc.abstractmethod
    def detector(self):
        pass

    @property
    @abc.abstractmethod
    def instrument(self):
        pass

    @property
    @abc.abstractmethod
    def observatory(self):
        pass

    @property
    @abc.abstractmethod
    def processing_level(self):
        """
        The level to which the data has been processed.
        """

    @property
    @abc.abstractmethod
    def observer_location(self):
        """
        Coordinate of observatory location based on header info.
        """

    @property
    @abc.abstractmethod
    def date_reference(self):
        """
        The base time from which time axis values are measured.

        Often the same or very similar to date_start.
        """

    @property
    @abc.abstractmethod
    def date_start(self):
        pass

    @property
    @abc.abstractmethod
    def date_end(self):
        pass

    @property
    @abc.abstractmethod
    def version(self):
        """
        The data version.
        """


class RemoteSensorMetaABC(MetaABC):
    @property
    @abc.abstractmethod
    def rsun_meters(self):
        """
        Solar radius in units of length.
        """

    @property
    @abc.abstractmethod
    def rsun_angular(self):
        """
        Solar radius in angular units as seen from observatory.
        """

    @property
    @abc.abstractmethod
    def distance_to_sun(self):
        """
        Distance to Sun center from observatory.
        """


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @property
    @abc.abstractmethod
    def spectral_window(self):
        pass

    @property
    @abc.abstractmethod
    def observing_mode_id(self):
        """
        Unique identifier for the observing mode.

        Often referred to as OBS ID.
        """

    @property
    @abc.abstractmethod
    def observer_radial_velocity(self):
        """
        Velocity of observatory in direction of source.
        """
