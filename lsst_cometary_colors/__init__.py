import logging
from collections import defaultdict
import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.io import ascii
from sbpy.activity import Afrho, Haser
from sbpy.data import Ephem
import sbpy.units as sbu
from sbpy.spectroscopy import SpectralGradient, Reddening
from sbpy.calib import solar_fluxd

# Filter reference data, from Willmer (2018)
bandwidth = {
    "u": 547 * u.AA,
    "g": 1333 * u.AA,
    "r": 1338 * u.AA,
    "i": 1209 * u.AA,
    "z": 994 * u.AA,
    "y": 814 * u.AA,
}

lambda_0 = {  # mean energy wavelength
    "u": 3681 * u.AA,
    "g": 4864 * u.AA,
    "r": 6250 * u.AA,
    "i": 7565 * u.AA,
    "z": 8702 * u.AA,
    "y": 9722 * u.AA,
}

# effective wavelength for LSST filters already built-into sbpy
lambda_eff = {
    filt: solar_fluxd.get()[f"LSST {filt}(lambda eff)"]
    for filt in "ugrizy"
}


class Coma:
    """Dust and gas super class."""

    pass


class Gas(Coma):
    """Gas coma, based on Haser models.


    Parameters
    ----------
    rh : astropy Quantity
        Heliocentric distance.

    Q_CN : astropy Quantity
        CN production rate of the comet.

    NH_CN, C3_CN, C2_CN, NH2 _CN: float
        Molecular production rates normalized by Q(CN).  Default parameters are
        "typical" comets from A'Hearn et al. (1995) (based on CN/OH, NH/OH,
        C3/OH, and C2/CN), and Fink (2009) (NH2/H2O, CN/H2O).

    """

    def __init__(self, eph, QCN, NH_CN=1.3, C3_CN=0.08, C2_CN=0.87, NH2_CN=1.3):
        self.eph = eph

        self.Q = {
            "CN": QCN,
            "NH": QCN * NH_CN,
            "C3": QCN * C3_CN,
            "C2": QCN * C2_CN,
            "NH2": QCN * NH2_CN,
        }
        logging.info("Initialized gas coma with production rates:")
        logging.info(str(self.Q))

        # gamma and g-factor rh scale factor
        rh_scale = (self.eph["rh"][0] / u.au) ** -2

        # key by molecule
        self.haser = {}
        for row in ascii.read("scalelengths.csv", delimiter=","):
            self.haser[row["molecule"]] = Haser(
                self.Q[row["molecule"]],
                0.85 * u.km / u.s,
                row["gamma parent"] * u.km * rh_scale,
                row["gamma daughter"] * u.km * rh_scale,
            )

        # key by molecule+band name, only save filters with throughput > 0
        self.bands = {}
        for row in ascii.read("emission-bands.csv", delimiter=","):
            if not np.isfinite(row["gfactor"]):
                continue
            k = f'{row["molecule"]} {row["band"]}'
            self.bands[k] = {
                "name": row["band"],
                "molecule": row["molecule"],
                "gfactor": row["gfactor"] / u.s * rh_scale,
                "wave": row["wave"] * u.AA,
                "throughput": {filt: row[filt] for filt in "ugrizy" if row[filt] > 0},
                "reference": row["reference"],
            }

    def bands_in_filter(self, filt, molecule=None, band_name=None):
        """List of bands in a filter, optionally limited to a specific molecule or band."""
        bands = []
        for band in self.bands.values():
            if molecule not in [None, band["molecule"]]:
                continue

            if band_name not in [None, band["name"]]:
                continue

            if filt in band["throughput"]:
                bands.append(band)

        return bands

    def fluxd(self, filt, aper, molecule=None, band_name=None):
        """Calculate effective flux density in requested LSST filter.


        Parameters
        ----------
        filt : string
            u, g, r, i, z, or y.

        aper : sbpy Aperture or astropy Quantity (length)
            Photometric aperture.

        molecule : str, optional
            Limit calculations to just this molecule.

        band_name : str, optional
            Limit calculations to just this band.

        """

        fluxd_total = 0 * u.W / u.m**2 / u.um
        A = 4 * np.pi * self.eph["Delta"][0] ** 2
        for band in self.bands_in_filter(filt, molecule=molecule, band_name=band_name):
            haser = self.haser[band["molecule"]]
            N = haser.total_number(aper)
            flux = (N * band["gfactor"] * const.h *
                    const.c / band["wave"]).to("W")
            fluxd = (flux * band["wave"] / lambda_0[filt] / bandwidth[filt] / A).to(
                "W/(m2 um)"
            )

            fluxd_total = fluxd_total + fluxd * band["throughput"][filt]
            logging.debug(
                "%s %s %f %s %s %s",
                filt,
                band["molecule"],
                np.log10(N),
                str(flux),
                str(fluxd),
                str(fluxd_total),
            )

        return fluxd_total

    def m(self, filt, aper):
        """Same as fluxd, but for magnitudes."""

        return self.fluxd(filt, aper).to(u.ABmag, u.spectral_density(lambda_eff[filt]))


class Dust(Coma):
    """Dust coma, based on Afρ model.


    Parameters
    ----------
    eph : sbpy Ephem
        Target ephemeris.

    a0frho : astropy Quantity or sbpy Afrho, length
        Dust production rate proxy at 0 deg phase angle at 0.55 μm (V-band),
        units of length.

    S_V : astropy Quantity, inverse length
        Spectral gradient at 0.55 μm (V-band).

    """

    def __init__(self, eph, a0frho, S_V=10 * u.percent / sbu.hundred_nm):
        self.eph = eph
        self.a0frho = Afrho(a0frho)
        self.S_V = SpectralGradient(S_V, wave0=0.55 * u.um)

    def fluxd(self, filt, aper):
        """Calculate effective flux density in requested LSST filter.


        Parameters
        ----------
        filt : string
            u, g, r, i, z, or y.

        aper : sbpy Aperture or astropy Quantity (length)
            Photometric aperture.


        Notes
        -----
        Reddening, as calculated here, is approximate.

        """

        fluxd = self.a0frho.to_fluxd(
            f"LSST {filt}", aper, self.eph, unit="W/(m2 um)", phasecor=True
        )
        r = Reddening(self.S_V)(lambda_eff[filt])
        return fluxd * r

    def m(self, filt, aper):
        """Same as fluxd, but for magnitudes."""

        return self.fluxd(filt, aper).to(u.ABmag, u.spectral_density(lambda_eff[filt]))


class Comet:
    """Comet = gas + dust.


    Parameters
    ----------
    eph : sbpy Ephem
        rh, delta, and phase.

    a0frho : astropy Quantity
        Dust production rate proxy at 0 deg phase angle at 0.55 μm, units of
        length.

    S_V : astropy Quantity, inverse length
        Spectral gradient at 0.55 μm (V-band).

    QCN_afrho : astropy Quantity
        Ratio of CN production rate to Afρ value.  Default is from A'Hearn et al.
        (1995):  CN / OH * OH / Afρ  = 10**(-2.50 + 25.82) = 2.1e23 / cm s.

    **gas_ratios :
        Additional keyword arguments are passed to ``Gas()``, e.g., NH_CN, C3_CN,
        C2_CN, and/or NH2_CN.

    """

    def __init__(
        self,
        eph,
        a0frho,
        S_V=0.1 / sbu.hundred_nm,
        QCN_afrho=2.1e23 / u.s / u.cm,
        **gas_ratios,
    ):
        self._QCN_afrho = QCN_afrho
        self._eph = eph

        self.dust = Dust(eph, a0frho, S_V=S_V)
        QCN = (QCN_afrho * a0frho.to("cm")).to("1/s")
        self.gas = Gas(eph, QCN, **gas_ratios)

    @property
    def Q(self):
        return self.gas.Q

    @property
    def eph(self):
        return self._eph

    def fluxd(self, filt, aper):
        """Calculate effective flux density in requested LSST filter.


        Parameters
        ----------
        filt : string
            u, g, r, i, z, or y.

        aper : sbpy Aperture or astropy Quantity (length)
            Photometric aperture.

        """

        dust = self.dust.fluxd(filt, aper)
        gas = self.gas.fluxd(filt, aper)
        return dust + gas

    def m(self, filt, aper):
        """Same as fluxd, but for magnitudes."""

        return self.fluxd(filt, aper).to(u.ABmag, u.spectral_density(lambda_eff[filt]))

    def fraction(self, filt, aper):
        """Same as fluxd, but for relative fraction.


        Returns
        -------
        frac : dict
            Relative fraction keyed by source (dust, gas, gas band, or
            molecule).

        """

        frac = {}

        dust = self.dust.fluxd(filt, aper)
        gas = self.gas.fluxd(filt, aper)
        total = dust + gas

        frac["dust"] = float(dust / total)
        frac["gas"] = float(gas / total)

        for band in self.gas.bands_in_filter(filt):
            molecule = band["molecule"]
            band_name = band["name"]
            fluxd = self.gas.fluxd(
                filt, aper, molecule=molecule, band_name=band_name)
            f = float(fluxd / total)
            frac[molecule] = frac.get(molecule, {})
            frac[molecule]["total"] = f + frac[molecule].get("total", 0)
            frac[molecule][band_name] = f

        return frac
