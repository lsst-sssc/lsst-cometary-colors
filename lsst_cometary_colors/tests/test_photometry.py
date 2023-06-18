import numpy as np
import astropy.units as u
import astropy.constants as const
from sbpy.activity import Afrho, CircularAperture
from sbpy.data import Ephem
from sbpy.calib import solar_fluxd
from sbpy.spectroscopy import SpectralGradient
import sbpy.units as sbu
from lsst_cometary_colors import Gas, Dust, Comet, bandwidth, lambda_eff


def test_no_gas():
    eph = Ephem.from_dict(
        {"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
    comet = Comet(eph, 1000 * u.cm, QCN_afrho=0 / u.s / u.cm)
    rho = 1e4 * u.km

    dust_only = comet.dust.fluxd("r", rho)
    gas_only = comet.gas.fluxd("r", rho)
    total = comet.fluxd("r", rho)
    print(dust_only, gas_only, total)
    assert u.isclose(dust_only, total)
    assert u.isclose(0, gas_only.value)


def test_gas_emits_photons():
    eph = Ephem.from_dict(
        {"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
    comet = Comet(eph, 1000 * u.cm, QCN_afrho=0 / u.s / u.cm)
    rho = 1e4 * u.km

    comet = Comet(eph, 1000 * u.cm, QCN_afrho=1e23 / u.s / u.cm)
    dust_only = comet.dust.fluxd("r", rho)
    total = comet.fluxd("r", rho)
    assert dust_only < total


def test_spectral_reddening():
    eph = Ephem.from_dict(
        {"rh": 2 * u.au, "delta": 1 * u.au, "phase": 0 * u.deg})
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
    """

    Calculations for the Li et al. 2017 paper from Michael S. P. Kelley's notes.
    The results were from an online tool at Lowell Obs. by Dave Schleicher.

    252P, 11:30, 40 px

          R (AU) =  1.098
      Delta (AU) =  0.16
        V (km/s) =   9.94
    Aper rad (") =   9.6

        Rho (km) =   1142.   log(Rho) =  3.06
    Area (cmxcm) = 4.10E+16

                  OH(0-0)  NH(delv=0) CN(delv=0)  C3     C2(delv=0)
          Flux:   2.17E-14  0.00E+00  1.14E-14  0.00E+00  0.00E+00
      log Flux:   -13.66      0.00    -13.94      0.00      0.00

      g-factor:   2.62E-15  7.10E-14  3.94E-13  8.29E-13  3.73E-13
        M(rho):   6.29E+26  0.00E+00  2.20E+24  0.00E+00  0.00E+00
    log M(rho):    26.80      0.00     24.34      0.00      0.00

        N(rho):   1.54E+10  0.00E+00  5.36E+07  0.00E+00  0.00E+00
    log N(rho):    10.19      0.00      7.73      0.00      0.00

            Xp:   2.89E+04  6.03E+04  1.57E+04  3.38E+03  2.65E+04
            Xd:   1.93E+05  1.81E+05  2.53E+05  3.26E+04  7.96E+04
    H Fraction:   5.27E-04  3.01E-04  6.66E-04  1.46E-02  1.30E-03

        M(tot):   1.19E+30  0.00E+00  3.30E+27  0.00E+00  0.00E+00
    log M(tot):    30.08      0.00     27.52      0.00      0.00

          TAUd:   1.93E+05  1.81E+05  2.53E+05  3.26E+04  7.96E+04

             Q:   6.18E+24  0.00E+00  1.30E+22  0.00E+00  0.00E+00
         log Q:    24.79      0.00     22.11      0.00      0.00

       Q/Q(OH):   1.00E+00  0.00E+00  2.11E-03  0.00E+00  0.00E+00

    """
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

    comet = Comet(
        eph,
        afrho,
        QCN_afrho=120e23 / u.s / afrho,
        S_V=S_V,
        C3_CN=0,
        C2_CN=0,
        NH_CN=0,
        NH2_CN=0,
    )

    band = comet.gas.bands_in_filter("u", "CN")[0]
    hnu = (const.h * const.c / 387 / u.nm).to("erg")
    g = 3.94e-13 * u.erg / u.s / hnu
    # Agreement to 20%?  Great!  (We're not accounting for heliocentric velocity.)
    assert np.isclose(band["gfactor"], g, rtol=0.2)

    fluxd = comet.gas.fluxd("u", rho)
    # CN filter width is about 56 AA, but u-band is wider
    scaled_fluxd_CN = fluxd_CN * 56 * u.AA / bandwidth["u"]
    # This is a crude estimate, so just check for agreement within 50%.
    assert u.isclose(fluxd, scaled_fluxd_CN, rtol=0.5)
