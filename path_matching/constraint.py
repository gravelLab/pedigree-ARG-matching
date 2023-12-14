from abc import abstractmethod


class Constraint:

    @abstractmethod
    def verify(self):
        pass
