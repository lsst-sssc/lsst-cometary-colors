import numpy as np
import astropy.units as u
from sbpy.activity import Afrho, CircularAperture
from sbpy.data import Ephem
from sbpy.calib import solar_fluxd
from sbpy.spectroscopy import SpectralGradient
import sbpy.units as sbu
from . import Gas, Dust, Comet, lambda_eff


def test_no_gas():
    eph = Ephem.from_dict({"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
    comet = Comet(eph, 1000 * u.cm, QCN_afrho=0 / u.s / u.cm)
    rho = 1e4 * u.km

    dust_only = comet.dust.fluxd("r", rho)
    gas_only = comet.gas.fluxd("r", rho)
    total = comet.fluxd("r", rho)
    print(dust_only, gas_only, total)
    assert u.isclose(dust_only, total)
    assert u.isclose(0, gas_only.value)


def test_gas_emits_photons():
    eph = Ephem.from_dict({"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
    comet = Comet(eph, 1000 * u.cm, QCN_afrho=0 / u.s / u.cm)
    rho = 1e4 * u.km

    comet = Comet(eph, 1000 * u.cm, QCN_afrho=1e23 / u.s / u.cm)
    dust_only = comet.dust.fluxd("r", rho)
    total = comet.fluxd("r", rho)
    assert dust_only < total


def test_spectral_reddening():
    eph = Ephem.from_dict({"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
    comet = Comet(eph, 1000 * u.cm, QCN_afrho=0 / u.s / u.cm)
    rho = 1e4 * u.km
    u_r = (
        comet.dust.fluxd("r", rho)
        / comet.dust.fluxd("u", rho)
        * solar_fluxd.get()["LSST u"]
        / solar_fluxd.get()["LSST r"]
    )

    # default reddening is 10% / 100 nm at 0.55 μm
    dwave_u_v = (0.55 * u.um - lambda_eff["u"]).to("um") / (0.1 * u.um)
    dwave_v_r = (lambda_eff["r"] - 0.55 * u.um).to("um") / (0.1 * u.um)

    assert np.isclose((1 + dwave_v_r * 0.1) / (1 - dwave_u_v * 0.1), u_r)


def test_252P():
    # Li et al. 2017
    eph = Ephem.from_dict(
        {"rh": 1.098 * u.au, "delta": 0.164 * u.au, "phase": 51.3 * u.deg}
    )
    eph["phase"] = 0 * u.deg  # avoid phase corrections
    rho = CircularAperture(19.2 * u.arcsec).as_length(eph["delta"][0])

    # flux density in CN filter, assume dominated by CN
    m_CN = 10.3 * u.ABmag
    fluxd_CN = m_CN.to("W/(m2 um)", u.spectral_density(lambda_eff["u"]))

    S_V = SpectralGradient(6.5 * u.percent / sbu.hundred_nm, wave0=0.55 * u.um)
    # BC, RC afrho = 24.3 and 28.9
    afrho = 26.5 * u.cm  # reddened from BC to V

    comet = Comet(eph, afrho, QCN_afrho=120e23 / u.s / afrho, S_V=S_V)
    fluxd = comet.gas.fluxd("u", rho)
    assert u.isclose(fluxd, fluxd_CN)
