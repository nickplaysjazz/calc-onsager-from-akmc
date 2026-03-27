from pathlib import Path

from constants import expected_options, SinkGeometry
from error_handler import raise_error
from io_handler import read_settings

def main():
    # Read in values from settings.ini
    settings = read_settings("settings.ini", expected_options)

    if settings["sink_type"] == SinkGeometry.BULK or settings["sink_type"] == SinkGeometry.GRAIN_BOUNDARY:
        raise_error(Path(__file__).name, f"Sink geometry {settings["sink_type"]} not implemented!")


    print(settings)

if __name__ == "__main__":
    main()
