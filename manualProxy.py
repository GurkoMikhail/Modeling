from multiprocessing.managers import BaseProxy
from sources import Source


def MakeProxy(target):
    """ Create a derived Proxy class for `target`. """
    
    def __getattr__(self, key):
        result = self._callmethod('__getattribute__', (key,))
        if callable(result):
            def wrapper(*args, **kwargs):
                return self._callmethod(key, args, **kwargs)
            return wrapper
        return result

    def __setattr__(self, key, value):
        self._callmethod('__setattr__', (key, value))

    def __delattr__(self, key):
        self._callmethod('__delattr__', (key,))

    dic = {'__getattr__': __getattr__}
    proxy_name = target.__name__ + "Proxy"
    ProxyType = type(proxy_name, (BaseProxy,), dic)
    exposed = dir(target)
    # exposed += ('__setattr__', '__delattr__')
    ProxyType._exposed_ = tuple(exposed)

    return ProxyType


SourceProxy = MakeProxy(Source)
