import configparser
import os

from pathlib import Path

from error_handler import raise_error


def parse_ini(filename: str) -> configparser.ConfigParser:
    """Attempts to read in a .ini file in the same directory.

    Returns configparser object."""

    BASE_DIR = Path(__file__).resolve().parent
    config_path = BASE_DIR / filename

    if not os.path.isfile(config_path):
        raise_error(
            Path(__file__).name, f"Configuration file {filename} does not exist."
        )

    config = configparser.ConfigParser()
    try:
        config.read(config_path)
    except IOError as e:
        raise_error(
            Path(__file__).name, f"Unable to read configuration file {filename}!", e
        )
    except FileNotFoundError as e:
        raise_error(
            Path(__file__).name, f"Configuration file {filename} does not exist", e
        )
    except PermissionError as e:
        raise_error(
            Path(__file__).name,
            f"Adequate permission to access {filename} not granted",
            e,
        )

    return config


def read_settings(filename: str, schema: dict) -> dict:
    """Attempts to read in settings.ini in the same directory.

    Returns settings dictionary with key-value pairs obtained from settings.ini."""

    try:
        # First, obtain configparser object
        settings = parse_ini(filename)

        # Then translate Options section into dictionary
        ret_dict = settings["Options"]

        # Throw an error for unexpected options or types
        ret_dict = validate_config(ret_dict, schema)

    except configparser.ParsingError as e:
        raise_error(
            Path(__file__).name, "Configuration file does not follow legal syntax", e
        )
    except configparser.NoSectionError as e:
        raise_error(
            Path(__file__).name, "Configuration file does not follow expected syntax", e
        )
    except KeyError as e:
        raise_error(
            Path(__file__).name, "Configuration file does not follow expected syntax", e
        )

    return ret_dict


def validate_config(raw_data, schema):
    """Takes in a raw_data dictionary and a schema dictionary with variable types as values.

    Checks that raw_data has all variable names given in the schema (and no extras). Also checks that each variable type is correct. Otherwise throws an error.
    """

    expected_keys = set(schema.keys())
    provided_keys = set(raw_data.keys())

    if expected_keys != provided_keys:
        set_diff = provided_keys - expected_keys
        if set_diff:
            raise_error(
                Path(__file__).name,
                f"Found option(s) {set_diff} and did not recognize!",
            )
        set_diff = expected_keys - provided_keys
        raise_error(
            Path(__file__).name, f"Expected option(s) {set_diff} and could not locate!"
        )

    validated = {}
    for key, type_func in schema.items():
        try:
            # Try to cast the string from the config file
            validated[key] = type_func(raw_data[key])
        except (ValueError, KeyError) as e:
            raise_error(
                Path(__file__).name,
                f"Setting '{key}' must be of type {type_func.__name__}",
                e,
            )

    return validated
