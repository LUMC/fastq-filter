from abc import ABC, abstractmethod
from typing import Union

class Filter(ABC):
    threshold: Union[int, float]
    passed: int
    total: int

    @abstractmethod
    def __init__(self): ...

    @abstractmethod
    def passes_filter(self, __values: str) -> bool: ...

