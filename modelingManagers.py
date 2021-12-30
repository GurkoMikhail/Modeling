from multiprocessing.managers import BaseManager
import sources
from manualProxy import SourceProxy
from inspect import isclass


class SourceManager(BaseManager):

    def __init__(self):
        super().__init__()
        for sourceName in dir(sources):
            SourceClass = getattr(sources, sourceName)
            if isclass(SourceClass):
                self.register(sourceName, SourceClass, SourceProxy)
        self.start()
        


