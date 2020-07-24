import abc
import copy

__all__ = ["Meta", "RemoteSensorMetaABC", "SlitSpectrographMetaABC"]


class RemoteSensorMetaABC(abc.ABCMeta):
    @abc.abstractmethod
    def detector(self):
        pass

    @abc.abstractmethod
    def instrument(self):
        pass

    @abc.abstractmethod
    def observatory(self):
        pass

    @abc.abstractmethod
    def processing_level(self):
        """The level to which the data has been processed."""
        pass

    @abc.abstractmethod
    def rsun_meters(self):
        """Solar radius in units of length."""
        pass

    @abc.abstractmethod
    def rsun_angular(self):
        """Solar radius in angular units as seen from observatory."""
        pass

    @abc.abstractmethod
    def distance_to_sun(self):
        """Distance to Sun center from observatory."""
        pass

    @abc.abstractmethod
    def observer_coordinate(self):
        """Coordinate of observatory location based on header info."""
        pass

    @abc.abstractmethod
    def date_reference(self):
        """The base time from which time axis values are measured.

        Often the same or very similar to date_start.
        """
        pass

    @abc.abstractmethod
    def date_start(self):
        pass

    @abc.abstractmethod
    def date_end(self):
        pass


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @abc.abstractmethod
    def spectral_window(self):
        pass

    @abc.abstractmethod
    def observing_mode_id(self):
        """Unique identifier for the observing mode. Often referred to as OBS ID."""
        pass

    @abc.abstractmethod
    def observatory_radial_velocity(self):
        """Velocity of observatory in direction of source."""
        pass


class Meta():
    def __init__(self, header):
        self.raw_meta = header
        self.original_header = copy.deepcopy(header)

    def get(self, key, default=None):
        return self._get_fits(key, default)

    def _get_fits(self, key, default=None):
        return self.raw_meta.get(key, default)

    def get_comment(self, key):
        return self._get_comment_fits(key)

    def _get_comment_fits(self, key):
        return self.raw_meta.comments[key]

    def update(self, key_value_pairs):
        """Update or add an entry in the metadata.

        Parameters
        ----------
        key_value_pairs: iterable
            An iterable of (key, value) pair tuples giving the metadata key/name
            and the value with which to update it.
            If key doesn't exist, a new entry is added.
        """
        self._update_fits(key_value_pairs)

    def _update_fits(self, key_value_pairs):
        for key, value in key_value_pairs:
            self.raw_meta[key] = value

    def update_comments(self, key_comment_pairs):
        """Update comment associated with existing metadata entry.

        Parameters
        ----------
        key_comment_pairs: iterable
            An iterable of (key, comment) pair tuples giving the metadata key/name
            and the new comment to assign to it.
        """
        self._update_comments_fits(key_comment_pairs)

    def _update_comments_fits(self, key_comment_pairs):
        for key, comment in key_comment_pairs:
            self.raw_meta.comments[key] = comment
