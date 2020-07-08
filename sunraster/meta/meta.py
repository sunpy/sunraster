__all__ = ["RemoteSensorMetaABC", "SlitSpectrographMetaABC"]

class RemoteSensorMetaABC():
    @property
    def detector(self):
        pass
    
    @property
    def instrument(self):
        pass
    
    @property
    def observatory(self):
        pass
    
    @property
    def processing_level(self):
        pass
    
    @property
    def rsun_meters(self):
        pass
    
    @property
    def rsun_obs(self):
        pass
    
    @property
    def dsun(self):
        pass
    
    @property
    def observer_coordinate(self):
        pass
    
    @property
    def date(self):
        pass


class SlitSpectrographMetaABC(RemoteSensorMetaABC):
    @property
    def spectral_window(self):
        pass
    
    @property
    def obsid(self):
        pass
    
    @property
    def obsvr(self):
        pass
