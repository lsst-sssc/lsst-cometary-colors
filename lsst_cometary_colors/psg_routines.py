import os
import requests
from urllib.parse import quote

import astropy.units as u
from synphot import SourceSpectrum

@u.quantity_input(fwhm=u.arcsec, resolution=u.nm)
def generate_psg_config_file(output_location, helio_dist=2.87796, geo_dist=2.9737, fwhm=4*u.arcsec, resolution=1*u.nm, afrho=100, activity=7.5e27):
    """Generates a comet configuration file in <output_location> for the
    NASA Planetary Spectrum Generator (PSG) which can be uploaded or submitted
    via the API (see `generate_psg_spectrum()`).
    The optional configurable parameters are:
    * helio_dist: heliocentric distance of comet from the Sun (au),
    * geo_dist: geocentric distance of comet from the Earth (au),
    * fwhm: the FWHM of the instrument FOV (Astropy Quantity, transformable to arcsec),
    * resolution: the resolution of the generated spectrum (covers 300-1200nm; Astropy Quantity, transformable to nm),
    * afrho: Afrho (cm; used as unit of dust abundance for PSG),
    * activity: Atmospheric activity (molecules/sec)

    The pathname of the config file is returned.
    """

    config_lines = ['<OBJECT>Comet',
                    '<OBJECT-DIAMETER>7.00',
                    '<OBJECT-GRAVITY>0.600',
                    '<OBJECT-GRAVITY-UNIT>rho',
                   f'<OBJECT-STAR-DISTANCE>{helio_dist:}',
                    '<OBJECT-STAR-TYPE>G',
                    '<OBJECT-STAR-TEMPERATURE>5777',
                    '<OBJECT-STAR-RADIUS>1.0',
                    '<OBJECT-PERIAPSIS>0.00',
                    '<OBJECT-ECCENTRICITY>0.00000',
                    '<OBJECT-INCLINATION>90.00',
                    '<OBJECT-OBS-LONGITUDE>0.00',
                    '<OBJECT-OBS-LATITUDE>0.00',
                    '<OBJECT-STAR-LONGITUDE>0.00',
                    '<OBJECT-STAR-LATITUDE>0.00',
                    '<OBJECT-SEASON>0.00',
                    '<GEOMETRY>Observatory',
                    '<GEOMETRY-OFFSET-NS>0.0',
                    '<GEOMETRY-OFFSET-EW>0.0',
                    '<GEOMETRY-OFFSET-UNIT>arcsec',
                   f'<GEOMETRY-OBS-ALTITUDE>{geo_dist:}',
                    '<GEOMETRY-ALTITUDE-UNIT>AU',
                    '<GEOMETRY-STAR-DISTANCE>-1.000000e+00',
                    '<GEOMETRY-SOLAR-ANGLE>48.121',
                    '<GEOMETRY-OBS-ANGLE>48.121',
                    '<GEOMETRY-PLANET-FRACTION>1.000e+00',
                    '<GEOMETRY-BRDFSCALER>0.999',
                    '<GEOMETRY-AZIMUTH>0.000',
                    '<GEOMETRY-STAR-FRACTION>0.000000e+00',
                    '<GEOMETRY-ROTATION>0.00,0.00',
                    '<ATMOSPHERE-STRUCTURE>Coma',
                   f'<ATMOSPHERE-PRESSURE>{activity:}',
                    '<ATMOSPHERE-PUNIT>gas',
                    '<ATMOSPHERE-TEMPERATURE>38.1',
                    '<ATMOSPHERE-WEIGHT>471.57',
                    '<ATMOSPHERE-NGAS>4',
                    '<ATMOSPHERE-GAS>OH,CN,C2,NH',
                    '<ATMOSPHERE-TYPE>GSFC[OH],GSFC[CN],GSFC[C2],GSFC[NH]',
                    '<ATMOSPHERE-ABUN>100,0.35,0.2,0.3',
                    '<ATMOSPHERE-UNIT>pct,pct,pct,pct',
                    '<ATMOSPHERE-TAU>2.4e4 1.6e5,1.3e4 2.1e5,2.2e4 6.6e4,5.0e4 1.5e5',
                    '<SURFACE-TEMPERATURE>162.3',
                    '<SURFACE-ALBEDO>0.04',
                    '<SURFACE-EMISSIVITY>0.2',
                    '<SURFACE-NSURF>0',
                    '<SURFACE-SURF>',
                    '<SURFACE-TYPE>',
                    '<SURFACE-ABUN>',
                    '<SURFACE-UNIT>',
                    '<SURFACE-THICK>',
                   f'<SURFACE-GAS-RATIO>{afrho:}',
                    '<SURFACE-GAS-UNIT>afrho',
                    '<SURFACE-MODEL>Lommel-Seeliger,HG1,0.000',
                    '<ATMOSPHERE-NMAX>0',
                    '<ATMOSPHERE-LMAX>0',
                    '<ATMOSPHERE-NAERO>0',
                    '<ATMOSPHERE-AEROS>',
                    '<ATMOSPHERE-ATYPE>',
                    '<ATMOSPHERE-AABUN>',
                    '<ATMOSPHERE-AUNIT>',
                    '<ATMOSPHERE-ASIZE>',
                    '<ATMOSPHERE-ASUNI>',
                    '<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all',
                    '<SURFACE-PHASEG>',
                    '<GENERATOR-INSTRUMENT>user',
                    '<GENERATOR-RANGE1>300',
                    '<GENERATOR-RANGE2>1200',
                    '<GENERATOR-RANGEUNIT>nm',
                   f'<GENERATOR-RESOLUTION>{resolution.to(u.nm).value:}',
                   f'<GENERATOR-RESOLUTIONUNIT>{resolution.to(u.nm).unit:}',
                    '<GENERATOR-RADUNITS>Wm2um',
                   f'<GENERATOR-BEAM>{fwhm.to(u.arcsec).value:}',
                   f'<GENERATOR-BEAM-UNIT>{fwhm.to(u.arcsec).unit:}',
                    '<GENERATOR-RESOLUTIONKERNEL>N',
                    '<GENERATOR-TELESCOPE>SINGLE',
                    '<GENERATOR-NOISE>NO',
                    '<GENERATOR-TRANS-APPLY>N',
                    '<GENERATOR-TRANS-SHOW>N',
                    '<GENERATOR-TRANS>02-01',
                    '<GENERATOR-LOGRAD>N',
                    '<GENERATOR-GAS-MODEL>Y',
                    '<GENERATOR-CONT-MODEL>Y',
                    '<GENERATOR-CONT-STELLAR>Y',
                   ]

    filename= os.path.join(output_location, f"psg_config.txt")

    with open(filename, 'w') as fh:
        fh.writelines([line + '\n' for line in config_lines])

    return filename

def generate_psg_spectrum(config_file):
    """Calls the Planetary Spectrum Generator via the API service to compute
    a radiance spectrum with the configuration from <config_file> and saves it to disk
    as 'psg_spectrum.txt' in the same directory as the config_file (overwrites any
    existing version).
    The path to the spectrum is returned"""

    psg_url = "https://psg.gsfc.nasa.gov/api.php"
    with open(config_file, 'r') as fh:
        psg_data = fh.read()

    # URLencode the contents of the PSG config file and then pass the whole thing
    # as the value for the `file` parameter in the POST
    payload = quote(psg_data)
    payload = 'file=' + payload

    # Create output file in the same directory as the config file (Could maybe
    # make a unique name or pull some of the params out of the config file
    filename = os.path.join(os.path.dirname(os.path.abspath(config_file)), 'psg_spectrum.txt')

    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    resp = requests.post(psg_url, stream=True, headers=headers, data=payload)
    if resp.status_code in [200, 201]:
        with open(filename, 'wb') as fd:
            for chunk in resp.iter_content(chunk_size=128):
                fd.write(chunk)
    else:
        print(f"Error retrieving spectrum from PSG. (Status code={resp.status_code:})")
        filename = None

    return filename

def read_psg_spectrum(spectrum):
    """Reads a spectrum file produced by PSG from <spectrum> and returns a
    `synphot.SourceSpectrum`
    """

    # Units for PSG config files from generate_psg_config_file()
    # (Could read the units in the PSG header if we wanted to)
    psg_units = u.W/u.m**2/u.um

    sourcespec = SourceSpectrum.from_file(spectrum, wave_unit=u.nm, flux_unit=psg_units)

    return sourcespec
