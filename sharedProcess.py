import multiprocessing as mp
from multiprocessing import Process, Queue


class FakePopen:

    def poll(self, *args):
        pass
    
    def terminate(self, *args):
        pass
    
    def wait(self, *args):
        pass


class SharedMethod:

    def __init__(self, name, input_queue, output_queue):
        self.name = name
        self.input_queue = input_queue
        self.output_queue = output_queue

    def __call__(self, *args):
        self.input_queue.put([self.name, *args])
        return self.output_queue.get()


class SharedProperty:

    def __init__(self, name):
        self.name = name

    def __get__(self, instance, owner):
        return instance.__dict__[self.name]()

    def __set__(self, instance, value):
        pass
        # instance.__dict__[self.name] = value


class SharedProcess(Process):

    def __init__(self, daemon=None):
        super().__init__(daemon=daemon)
        self.input_queue = Queue(maxsize=1)
        self.output_queue = Queue(maxsize=1)
        self._exceptions = ['input_queue', 'output_queue', *dir(Process)]
    
    def start(self):
        super().start()
        if mp.get_start_method() == 'spawn':
            self._popen = FakePopen()
        for name in dir(self):
            if name not in self._exceptions and name[0] != '_':
                isCallable = callable(getattr(self, name))
                setattr(self, name, SharedMethod(name, self.input_queue, self.output_queue))
                if not isCallable:
                    setattr(self.__class__, name, SharedProperty(name))
                    
    def run(self):
        print(f'{self.name} started')
        for name, *args in iter(self.input_queue.get, None):
            output = getattr(self, name)
            if callable(output):
                output = output(*args)
            self.output_queue.put(output)

