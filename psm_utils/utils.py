"""Various utility functions."""

from typing import Optional

from pyteomics.mass import nist_mass


def mass_to_mz(mass: float, charge: int, adduct_mass: Optional[float] = None) -> float:
    """
    Convert mass to m/z.

    Parameters
    ----------
    mass
        Mass of the uncharged ion without adducts.
    charge
        Charge of the ion.
    adduct_mass
        Mass of the charge-carrying adduct. Defaults to the mass of a proton.

    """
    if adduct_mass is None:
        adduct_mass = nist_mass["H"][1][0]
    return (mass + charge * adduct_mass) / charge


def mz_to_mass(mz: float, charge: int, adduct_mass: Optional[float] = None) -> float:
    """
    Convert m/z to mass.

    Parameters
    ----------
    mz
        m/z of the charged ion and adducts.
    charge
        Charge of the ion.
    adduct_mass
        Mass of the charge-carrying adduct. Defaults to the mass of a proton.

    """
    if adduct_mass is None:
        adduct_mass = nist_mass["H"][1][0]
    return mz * charge - charge * adduct_mass
