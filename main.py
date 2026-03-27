from pathlib import Path

from constants import expected_options, SinkGeometry
from error_handler import raise_error
from io_handler import read_settings

def get_pairwise_interactions(settings):
    pairwise = {}

    pairwise["eps_aa"] = settings["cohesive_aa"] / 6
    pairwise["eps_bb"] = settings["cohesive_bb"] / 6
    pairwise["eps_ab"] = (pairwise["eps_aa"] + pairwise["eps_bb"] + settings["ordering_ab"]) / 2

    pairwise["eps_av"] = (settings["vac_for_a"] + 6*pairwise["eps_aa"])/12
    pairwise["eps_bv"] = (settings["vac_for_b"] + 6*pairwise["eps_bb"])/12
    pairwise["eps_a,aa"] = (settings["int_aa_for_a"] + 6*pairwise["eps_aa"])/12
    pairwise["eps_a,ab"] = (settings["int_ab_for_a"] + 6*pairwise["eps_aa"])/12
    pairwise["eps_a,bb"] = (settings["int_bb_for_a"] + 6*pairwise["eps_aa"])/12
    pairwise["eps_b,aa"] = (settings["int_aa_for_b"] + 6*pairwise["eps_bb"])/12
    pairwise["eps_b,ab"] = (settings["int_ab_for_b"] + 6*pairwise["eps_bb"])/12
    pairwise["eps_b,bb"] = (settings["int_bb_for_b"] + 6*pairwise["eps_bb"])/12

    return pairwise

def main():
    # Read in values from settings.ini
    settings = read_settings("settings.ini", expected_options)

    if (
        settings["sink_type"] == SinkGeometry.BULK
        or settings["sink_type"] == SinkGeometry.GRAIN_BOUNDARY
    ):
        raise_error(
            Path(__file__).name,
            f"Sink geometry {settings["sink_type"]} not implemented!",
        )

    pairwise = get_pairwise_interactions(settings)

    for key in pairwise: print(key, pairwise[key])


if __name__ == "__main__":
    main()
