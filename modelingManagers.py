from multiprocessing.managers import BaseManager
import multiprocessing.process
import itertools
import sources
from manualProxy import SourceProxy
from inspect import isclass


class SourceManager(BaseManager):

    def __init__(self):
        super().__init__()
        count = next(multiprocessing.process._process_counter) - 1
        multiprocessing.process._process_counter = itertools.count(count)
        for sourceName in dir(sources):
            SourceClass = getattr(sources, sourceName)
            if isclass(SourceClass):
                self.register(sourceName, SourceClass, SourceProxy)
        self.start()
        


