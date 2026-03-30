from enum import Enum


class SinkGeometry(Enum):
    DISLOCATION = "Dislocation"
    BULK = "Bulk"
    GRAIN_BOUNDARY = "Grain Boundary"


expected_options = {
    "solute_composition": float,
    "dim1": int,
    "dim2": int,
    "dim3": int,
    "sink_type": SinkGeometry,
    "dpa_rate": float,
    "ordering_ab": float,
    "cohesive_aa": float,
    "cohesive_bb": float,
    "vac_for_a": float,
    "vac_for_b": float,
    "int_aa_for_a": float,
    "int_ab_for_a": float,
    "int_bb_for_a": float,
    "int_aa_for_b": float,
    "int_ab_for_b": float,
    "int_bb_for_b": float,
    "saddle_av": float,
    "saddle_bv": float,
    "saddle_a_aa": float,
    "saddle_a_ab": float,
    "saddle_a_bb": float,
    "saddle_b_aa": float,
    "saddle_b_ab": float,
    "saddle_b_bb": float,
    "att_freq": float,
    "kt": float,
    "lattice_param": float,
}
