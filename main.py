import math

from pathlib import Path

from constants import expected_options, SinkGeometry
from error_handler import raise_error
from io_handler import read_settings


def get_pairwise_interactions(settings):
    pairwise = {}

    pairwise["eps_aa"] = settings["cohesive_aa"] / 6
    pairwise["eps_bb"] = settings["cohesive_bb"] / 6
    pairwise["eps_ab"] = (
        pairwise["eps_aa"] + pairwise["eps_bb"] + settings["ordering_ab"]
    ) / 2

    pairwise["eps_av"] = (settings["vac_for_a"] + 6 * pairwise["eps_aa"]) / 12
    pairwise["eps_bv"] = (settings["vac_for_b"] + 6 * pairwise["eps_bb"]) / 12
    pairwise["eps_a,aa"] = (settings["int_aa_for_a"] + 6 * pairwise["eps_aa"]) / 12
    pairwise["eps_a,ab"] = (settings["int_ab_for_a"] + 6 * pairwise["eps_aa"]) / 12
    pairwise["eps_a,bb"] = (settings["int_bb_for_a"] + 6 * pairwise["eps_aa"]) / 12
    pairwise["eps_b,aa"] = (settings["int_aa_for_b"] + 6 * pairwise["eps_bb"]) / 12
    pairwise["eps_b,ab"] = (settings["int_ab_for_b"] + 6 * pairwise["eps_bb"]) / 12
    pairwise["eps_b,bb"] = (settings["int_bb_for_b"] + 6 * pairwise["eps_bb"]) / 12

    return pairwise


def get_five_jump_frequencies(settings, pairwise):
    five_jump_energies = {}

    five_jump_energies["e0"] = (
        settings["saddle_av"] - 11 * pairwise["eps_aa"] - 12 * pairwise["eps_av"]
    )
    five_jump_energies["e1"] = (
        settings["saddle_av"]
        - 11 * pairwise["eps_av"]
        - 1 * pairwise["eps_bv"]
        - 10 * pairwise["eps_aa"]
        - 1 * pairwise["eps_ab"]
    )
    five_jump_energies["e2"] = (
        settings["saddle_bv"]
        - 1 * pairwise["eps_bv"]
        - 11 * pairwise["eps_ab"]
        - 11 * pairwise["eps_av"]
    )
    five_jump_energies["e3"] = (
        settings["saddle_av"]
        - 11 * pairwise["eps_aa"]
        - 11 * pairwise["eps_av"]
        - 1 * pairwise["eps_bv"]
    )
    five_jump_energies["e4"] = (
        settings["saddle_av"]
        - 12 * pairwise["eps_av"]
        - 10 * pairwise["eps_aa"]
        - 1 * pairwise["eps_ab"]
    )

    five_jump_freq = {}
    for k in five_jump_energies:
        new_k = "w" + k.strip("e")
        five_jump_freq[new_k] = settings["att_freq"] * math.exp(
            -five_jump_energies[k] / settings["kt"]
        )

    # Koiwa & Ishioka coefficients for five-jump model
    c = five_jump_freq["w4"] / five_jump_freq["w0"]
    B = [180.3, 924.3, 1338.1, 40.1, 253.3, 596, 435.3]
    F = 1 - ((10*c**4 + B[0]*c**3 + B[1]*c**2 +B[2]*c)/(2*c**4+B[3]*c**3+B[4]*c**2+B[5]*c+B[6])/7)

    return five_jump_freq, F


def calc_vacancy_onsager_coeff(settings, five_jump_model, vac_conc):
    """Tucker, Najafabadi, Allen, Morgan, JNM, 2010"""
    n = 1 # number of lattice sites
    nn = 12 # number of nearest neighbors
    s = settings["lattice_param"] * math.sqrt(2)/2 # nn jump distance for fcc

    # Assume that the input solute concentration and equil vac conc
    # are the isolated fractions of solute & vacancies
    # In reality it should be that the input is the TOTAL number, but
    # the difference is small at equilibrium anyway
    bv_pair_conc = nn * settings["solute_composition"] * vac_conc * (five_jump_model["w3"]/five_jump_model["w4"])







def main():
    # Read in values from settings.ini
    settings = read_settings("settings.ini", expected_options)

    # Throw error for things not yet implemented
    if (
        settings["sink_type"] == SinkGeometry.BULK
        or settings["sink_type"] == SinkGeometry.GRAIN_BOUNDARY
    ):
        raise_error(
            Path(__file__).name,
            f"Sink geometry {settings["sink_type"]} not implemented!",
        )

    pairwise = get_pairwise_interactions(settings)

    eq_vac_conc = math.exp(-settings["vac_for_a"] / settings["kt"])
    eq_int_conc = math.exp(-settings["int_aa_for_a"] / settings["kt"])

    five_jump_model, F = get_five_jump_frequencies(settings, pairwise)

    for key in five_jump_model:
        print(key, five_jump_model[key])


if __name__ == "__main__":
    main()
