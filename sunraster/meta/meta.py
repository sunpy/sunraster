import abc

__all__ = ["RemoteSensorMetaABC", "SlitSpectrographMetaABC"]

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
        pass
    
    @abc.abstractmethod
    def rsun_meters(self):
        pass
    
    @abc.abstractmethod
    def rsun_obs(self):
        pass
    
    @abc.abstractmethod
    def dsun(self):
        pass
    
    @abc.abstractmethod
    def observer_coordinate(self):
        pass
    
    @abc.abstractmethod
    def date(self):
        pass


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @abc.abstractmethod
    def spectral_window(self):
        pass
    
    @abc.abstractmethod
    def obsid(self):
        pass
    
    @abc.abstractmethod
    def obsvr(self):
        pass
