from abc import ABC, abstractmethod
from typing import Union

from dnaio import SequenceRecord


class Filter(ABC):
    threshold: Union[int, float]
    passed: int
    total: int

    @abstractmethod
    def __init__(self): ...

    def __call__(self, __record: SequenceRecord) -> bool: ...

