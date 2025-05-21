from importlib import metadata

__all__ = ["IsotonicPEP", "__version__"]

try:
    __version__ = metadata.version(__name__)
except metadata.PackageNotFoundError:
    __version__ = "0.0.0"