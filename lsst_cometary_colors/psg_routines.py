import os

import astropy.units as u

@u.quantity_input(fwhm=u.arcsec)
def generate_psg_config_file(output_location, helio_dist=2.87796, geo_dist=2.9737, fwhm=4*u.arcsec):

    config_lines = ['<OBJECT>Comet',
                    '<OBJECT-DIAMETER>7.00',
                    '<OBJECT-GRAVITY>0.600',
                    '<OBJECT-GRAVITY-UNIT>rho',
                   f'<OBJECT-STAR-DISTANCE>{helio_dist:}',
                    '<OBJECT-STAR-TYPE>G',
                    '<OBJECT-STAR-TEMPERATURE>5777',
                    '<OBJECT-STAR-RADIUS>1.0',
                    '<GEOMETRY>Observatory',
                    '<GEOMETRY-OFFSET-NS>0.0',
                    '<GEOMETRY-OFFSET-EW>0.0',
                    '<GEOMETRY-OFFSET-UNIT>arcsec',
                   f'<GEOMETRY-OBS-ALTITUDE>{geo_dist:}',
                    '<GEOMETRY-ALTITUDE-UNIT>AU',
                    '<ATMOSPHERE-STRUCTURE>Coma',
                    '<ATMOSPHERE-PRESSURE>7.5e+27',
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
                    '<SURFACE-GAS-RATIO>100',
                    '<SURFACE-GAS-UNIT>afrho',
                    '<SURFACE-MODEL>Lommel-Seeliger,HG1,0.000',
                    '<GENERATOR-RESOLUTIONKERNEL>N',
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
                    '<GENERATOR-RESOLUTION>1',
                    '<GENERATOR-RESOLUTIONUNIT>nm',
                    '<GENERATOR-RADUNITS>Wm2um',
                   f'<GENERATOR-BEAM>{fwhm.to(u.arcsec).value:}',
                   f'<GENERATOR-BEAM-UNIT>{fwhm.to(u.arcsec).unit:}',
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