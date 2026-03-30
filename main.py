import math

from pathlib import Path

from constants import expected_options, SinkGeometry
from error_handler import raise_error
from io_handler import read_settings, read_csv


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
    F = 1 - (
        (10 * c**4 + B[0] * c**3 + B[1] * c**2 + B[2] * c)
        / (2 * c**4 + B[3] * c**3 + B[4] * c**2 + B[5] * c + B[6])
        / 7
    )

    return five_jump_freq, F


def get_interstitial_dilute_jump_frequencies(settings, pairwise):
    """Tucker, Najafabadi, Allen, Morgan, JNM, 2010"""
    interstitial_dilute_energy = {}

    # AA jump to NN
    interstitial_dilute_energy["e0"] = (
        settings["saddle_a_aa"] - 12 * pairwise["eps_a,aa"] - 11 * pairwise["eps_aa"]
    )
    # AB jump to NN
    interstitial_dilute_energy["eI"] = (
        settings["saddle_a_ab"] - 12 * pairwise["eps_a,ab"] - 11 * pairwise["eps_aa"]
    )
    # AA rotation
    interstitial_dilute_energy["eR0"] = 0.05
    # AB rotation
    interstitial_dilute_energy["eR"] = 0.05
    # AA+B->A+AB
    interstitial_dilute_energy["e2p"] = (
        settings["saddle_b_aa"]
        - 11 * pairwise["eps_a,aa"]
        - 1 * pairwise["eps_b,aa"]
        - 11 * pairwise["eps_ab"]
    )
    # A+AB->B+AA
    interstitial_dilute_energy["e2"] = (
        settings["saddle_a_ab"] - 12 * pairwise["eps_a,ab"] - 11 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to AA,A
    interstitial_dilute_energy["e1"] = (
        settings["saddle_a_aa"]
        - 1 * pairwise["eps_b,aa"]
        - 11 * pairwise["eps_a,aa"]
        - 1 * pairwise["eps_ab"]
        - 10 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to AA initially
    interstitial_dilute_energy["e3"] = (
        settings["saddle_a_aa"]
        - 1 * pairwise["eps_b,aa"]
        - 11 * pairwise["eps_a,aa"]
        - 11 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to A initially
    interstitial_dilute_energy["e4"] = (
        settings["saddle_a_aa"]
        - 12 * pairwise["eps_a,aa"]
        - 1 * pairwise["eps_ab"]
        - 10 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to AA,A; initial/final rotated by 90 deg
    interstitial_dilute_energy["e1p"] = (
        settings["saddle_a_aa"]
        - 1 * pairwise["eps_b,aa"]
        - 11 * pairwise["eps_a,aa"]
        - 1 * pairwise["eps_ab"]
        - 10 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to AA inetially; initial/final rotated by 90 deg
    interstitial_dilute_energy["e3p"] = (
        settings["saddle_a_aa"]
        - 1 * pairwise["eps_b,aa"]
        - 11 * pairwise["eps_a,aa"]
        - 11 * pairwise["eps_aa"]
    )
    # AA+A->A+AA, B 1NN to A inieially; initial/final rotated by 90 deg
    interstitial_dilute_energy["e4p"] = (
        settings["saddle_a_aa"]
        - 12 * pairwise["eps_a,aa"]
        - 1 * pairwise["eps_ab"]
        - 10 * pairwise["eps_aa"]
    )

    interstitial_dilute_freq = {}
    for k in interstitial_dilute_energy:
        new_k = "w" + k.strip("e")
        interstitial_dilute_freq[new_k] = settings["att_freq"] * math.exp(
            -interstitial_dilute_energy[k] / settings["kt"]
        )

    return interstitial_dilute_freq


def calc_vacancy_onsager_coeff(settings, five_jump_model, five_jump_F_factor, vac_conc):
    """Tucker, Najafabadi, Allen, Morgan, JNM, 2010"""
    n = 1  # number of lattice sites
    s = settings["lattice_param"] * math.sqrt(2) / 2  # nn jump distance for fcc

    w0 = five_jump_model["w0"]
    w1 = five_jump_model["w1"]
    w2 = five_jump_model["w2"]
    w3 = five_jump_model["w3"]
    w4 = five_jump_model["w4"]

    # Assume that the input solute concentration and equil vac conc
    # are the isolated fractions of solute & vacancies
    # In reality it should be that the input is the TOTAL number, but
    # the difference is small at equilibrium anyway
    bv_pair_conc = 12 * settings["solute_composition"] * vac_conc * (w4 / w3)

    omega = 2 * w1 + 2 * w2 + 7 * w3 * five_jump_F_factor
    a_aa_0 = 4 * w1 + 14 * w3
    a_aa_1 = (
        -2 * (3 * w3 - 2 * w1) ** 2
        + 28 * w3 * (1 - five_jump_F_factor) * ((w0 - w4) / w4)
        - 14
        * w3
        * (1 - five_jump_F_factor)
        * (2 * w1 + 2 * w2 + 7 * w3)
        * ((w0 - w4) / w4) ** 2
    ) / omega
    a_ab_1 = (
        w2
        * (
            2 * (3 * w3 - 2 * w1)
            + 14 * w3 * (1 - five_jump_F_factor) * ((w0 - w4) / w4)
        )
    ) / omega

    onsager_coeff_vac = {}
    onsager_coeff_vac["l_aa_v"] = (((n) * s**2) / 6) * (
        12 * vac_conc * (1 - 7 * settings["solute_composition"]) * w0
        + bv_pair_conc * a_aa_0
        + bv_pair_conc * a_aa_1
    )
    onsager_coeff_vac["l_ab_v"] = (((n) * s**2) / 6) * bv_pair_conc * a_ab_1
    onsager_coeff_vac["l_bb_v"] = (((n) * s**2) / 6) * (
        bv_pair_conc * w2 - ((2 * bv_pair_conc * w2**2) / omega)
    )
    onsager_coeff_vac["l_av_v"] = (
        -onsager_coeff_vac["l_aa_v"] - onsager_coeff_vac["l_ab_v"]
    )
    onsager_coeff_vac["l_bv_v"] = (
        -onsager_coeff_vac["l_bb_v"] - onsager_coeff_vac["l_ab_v"]
    )
    onsager_coeff_vac["l_vv_v"] = (
        onsager_coeff_vac["l_aa_v"]
        + 2 * onsager_coeff_vac["l_ab_v"]
        + onsager_coeff_vac["l_bb_v"]
    )

    norm_onsager_coeff_vac = {k: (v / vac_conc) for k, v in onsager_coeff_vac.items()}

    return norm_onsager_coeff_vac


def calc_interstitial_onsager_coeff(
    settings, pairwise, interstitial_dilute_frequencies, int_conc
):
    """Tucker, Najafabadi, Allen, Morgan, JNM, 2010"""
    n = 1  # number of lattice sites; this is 1 for FCC materials even for interstitials
    s = settings["lattice_param"] * math.sqrt(2) / 2  # nn jump distance for fcc

    aa_b_binding = (
        pairwise["eps_a,aa"]
        - pairwise["eps_aa"]
        - pairwise["eps_b,aa"]
        + pairwise["eps_ab"]
    )
    bi_pair_conc = (
        8
        * int_conc
        * settings["solute_composition"]
        * math.exp(aa_b_binding / settings["kt"])
    )
    a_type_bi_pair_conc = (
        bi_pair_conc  # Assume all orientations are equally available for nn jumps
    )

    w0 = interstitial_dilute_frequencies["w0"]
    wI = interstitial_dilute_frequencies["wI"]
    wR0 = interstitial_dilute_frequencies["wR0"]
    wR = interstitial_dilute_frequencies["wR"]
    w2p = interstitial_dilute_frequencies["w2p"]
    w2 = interstitial_dilute_frequencies["w2"]
    w1 = interstitial_dilute_frequencies["w1"]
    w3 = interstitial_dilute_frequencies["w3"]
    w4 = interstitial_dilute_frequencies["w4"]
    w1p = interstitial_dilute_frequencies["w1p"]
    w3p = interstitial_dilute_frequencies["w3p"]
    w4p = interstitial_dilute_frequencies["w4p"]

    A = (wR + wI) * (5 * w3 + w2p) + 5 * w3 * w2
    L_AA_f = (4 * n * s**2 * int_conc * w0) / 3
    L_AA_p = ((n * s**2 * a_type_bi_pair_conc * w3) / 6) * (
        (-7 * w0 / w4)
        + (16 * (w3 + 2 * w1 + w2p) / (5 * w3 + 2 * w1 + w2p))
        + (
            (12 * (wR + wI) * (2 * w3 + w2p) + 2 * w3 * w2)
            / ((wR + wI) * (5 * w3 + w2p) + 5 * w3 * w2)
        )
        + (6 * w4p * (w3p + 2 * w1p) / (w4 * (2 * w3p + w1p)))
    )

    onsager_coeff_int = {}
    onsager_coeff_int["l_aa_i"] = L_AA_f + L_AA_p
    onsager_coeff_int["l_ab_i"] = (n * s**2 * a_type_bi_pair_conc / 6) * (
        6 * wI * w2p * w3 / A
    )
    onsager_coeff_int["l_bb_i"] = (
        (n * s**2 * a_type_bi_pair_conc / 6)
        * ((w2p * wI) / (w2 * A))
        * (wR * (5 * w3 + w2p) + 5 * w3 * w2)
    )
    onsager_coeff_int["l_ai_i"] = (
        onsager_coeff_int["l_aa_i"] + onsager_coeff_int["l_ab_i"]
    )
    onsager_coeff_int["l_bi_i"] = (
        onsager_coeff_int["l_bb_i"] + onsager_coeff_int["l_ab_i"]
    )
    onsager_coeff_int["l_ii_i"] = (
        onsager_coeff_int["l_aa_i"]
        + 2 * onsager_coeff_int["l_ab_i"]
        + onsager_coeff_int["l_bb_i"]
    )

    norm_onsager_coeff_vac = {k: (v / int_conc) for k, v in onsager_coeff_int.items()}

    return norm_onsager_coeff_vac


def get_point_defect_conc(settings, eq_vac_conc, eq_int_conc, D_v, D_i):
    s = settings["lattice_param"] * math.sqrt(2) / 2  # nn jump distance for fcc

    vac_conc = eq_vac_conc + (settings["dim1"] * s) ** 2 * settings["dpa_rate"] / (
        12 * D_v
    )
    int_conc = eq_int_conc + (settings["dim1"] * s) ** 2 * settings["dpa_rate"] / (
        12 * D_i
    )

    return vac_conc, int_conc


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


    #----------------------------------------------------------------------------------
    # EDIT BELOW

    # EXAMPLE OVERRIDING VALUES WITH CSV VALUES
    # Override solute concentration and temperature with values read in from a .csv file
    temp, B_conc, pd_conc = read_csv("example.csv")

    pairwise = get_pairwise_interactions(settings)

    for t in range(len(temp)):
        settings["kt"] = temp[t]
        settings["solute_composition"] = B_conc[t]

        eq_vac_conc = math.exp(-settings["vac_for_a"] / settings["kt"])
        eq_int_conc = math.exp(-settings["int_aa_for_a"] / settings["kt"])

        five_jump_model, five_jump_F_factor = get_five_jump_frequencies(
            settings, pairwise
        )

        norm_onsager_coeff_vac = calc_vacancy_onsager_coeff(
            settings, five_jump_model, five_jump_F_factor, eq_vac_conc
        )

        interstitial_jumps = get_interstitial_dilute_jump_frequencies(
            settings, pairwise
        )

        norm_onsager_coeff_int = calc_interstitial_onsager_coeff(
            settings, pairwise, interstitial_jumps, eq_int_conc
        )

        vac_conc, int_conc = get_point_defect_conc(
            settings,
            eq_vac_conc,
            eq_int_conc,
            norm_onsager_coeff_vac["l_vv_v"],
            norm_onsager_coeff_int["l_ii_i"],
        )

        onsager_vac = {k: (v * vac_conc) for k, v in norm_onsager_coeff_vac.items()}
        onsager_int = {k: (v * int_conc) for k, v in norm_onsager_coeff_int.items()}

        ris_alpha_sign_factor = (onsager_vac["l_bv_v"] / onsager_vac["l_av_v"]) - (
            onsager_int["l_bi_i"] / onsager_int["l_ai_i"]
        )


        print(B_conc[t], ris_alpha_sign_factor)


if __name__ == "__main__":
    main()
