from Aries.outputs import LoggingConfig, PackageLogFilter
from .main import main


if __name__ == '__main__':
    with LoggingConfig(filters=[PackageLogFilter(packages="Cancer")]):
        main()
