import os
import tempfile, shutil
import pytest
import warnings

import astropy.units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SourceSpectrum

from lsst_cometary_colors.psg_routines import *


class TestGenerateConfig(object):

    @pytest.fixture
    def test_dir(self, tmp_path_factory):
        return tmp_path_factory.mktemp('tmp_lsst_', numbered=True)

    def test_exists(self, test_dir):

        config_file = generate_psg_config_file(test_dir)

        assert os.path.exists(config_file) is True

    def test_change_rh(self, test_dir):

        r_h = 2.1234

        config_file = generate_psg_config_file(test_dir, helio_dist=r_h)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        line = [line.rstrip() for line in lines if 'OBJECT-STAR-DISTANCE' in line]
        assert f'<OBJECT-STAR-DISTANCE>{r_h:}' == line[0]

    def test_change_delta(self, test_dir):

        delta = 1.001

        config_file = generate_psg_config_file(test_dir, geo_dist=delta)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        line = [line.rstrip() for line in lines if 'GEOMETRY-OBS-ALTITUDE' in line]
        assert f'<GEOMETRY-OBS-ALTITUDE>{delta:}' == line[0]

    def test_change_both(self, test_dir):

        r_h = 2.1234
        delta = 1.001

        config_file = generate_psg_config_file(test_dir, r_h, delta)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        line = [line.rstrip() for line in lines if 'GEOMETRY-OBS-ALTITUDE' in line]
        assert f'<GEOMETRY-OBS-ALTITUDE>{delta:}' == line[0]
        line = [line.rstrip() for line in lines if 'OBJECT-STAR-DISTANCE' in line]
        assert f'<OBJECT-STAR-DISTANCE>{r_h:}' == line[0]

    def test_fwhm_default(self, test_dir):

        config_file = generate_psg_config_file(test_dir)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        checks = { 'GENERATOR-BEAM' : 4.0,
                   'GENERATOR-BEAM-UNIT' : 'arcsec'
                 }
        for key, value in checks.items():
            line = [line.rstrip() for line in lines if key in line]
            assert f'<{key:}>{value:}' == line[0]

    def test_fwhm_newval(self, test_dir):

        fwhm = 5.0*u.arcsec
        config_file = generate_psg_config_file(test_dir, fwhm=fwhm)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        checks = { 'GENERATOR-BEAM' : fwhm.value,
                   'GENERATOR-BEAM-UNIT' : 'arcsec'
                 }
        for key, value in checks.items():
            line = [line.rstrip() for line in lines if key in line]
            assert f'<{key:}>{value:}' == line[0]

    def test_fwhm_newunits(self, test_dir):

        fwhm = 0.1*u.arcmin
        config_file = generate_psg_config_file(test_dir, fwhm=fwhm)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        checks = { 'GENERATOR-BEAM' : fwhm.to(u.arcsec).value,
                   'GENERATOR-BEAM-UNIT' : fwhm.to(u.arcsec).unit
                 }
        for key, value in checks.items():
            line = [line.rstrip() for line in lines if key in line]
            assert f'<{key:}>{value:}' == line[0]

    def test_fwhm_badunits(self, test_dir):

        fwhm = 5
        with pytest.raises(TypeError) as execinfo:
            config_file = generate_psg_config_file(test_dir, fwhm=fwhm)

    def test_resolution_newunits(self, test_dir):

        resolution = 1000*u.AA
        config_file = generate_psg_config_file(test_dir, resolution=resolution)

        with open(config_file, 'r') as fh:
            lines = fh.readlines()

        checks = { 'GENERATOR-RESOLUTION' : resolution.to(u.nm).value,
                   'GENERATOR-RESOLUTIONUNIT' : resolution.to(u.nm).unit
                 }
        for key, value in checks.items():
            line = [line.rstrip() for line in lines if key in line]
            assert f'<{key:}>{value:}' == line[0]

    def test_resolution_badunits(self, test_dir):

        resolution = 5
        with pytest.raises(TypeError) as execinfo:
            config_file = generate_psg_config_file(test_dir, resolution=resolution)


class TestGeneratePSGSpectrum:

    @pytest.fixture(scope='session')
    def config_file(self, tmp_path_factory):
        config_path = tmp_path_factory.mktemp('tmp_psg_')
        config_file = generate_psg_config_file(config_path, resolution=100*u.nm)
        return config_file

    def test_exists(self, config_file):

        spectrum_file = generate_psg_spectrum(config_file)

        assert os.path.exists(config_file) is True
        assert os.path.exists(spectrum_file) is True
        assert os.path.getsize(spectrum_file) == 2339

class TestReadPSGSpectrum:

    funit = u.W/u.m**2/u.micron

    test_spectrum_fp = os.path.abspath(os.path.join(__package__, 'lsst_cometary_colors', "tests", "data", "psg_spectrum.txt"))

    def test_spectrum_exists(self):
        assert os.path.exists(self.test_spectrum_fp) is True

    def test_spectrum(self):

        source_spec = read_psg_spectrum(self.test_spectrum_fp)

        assert type(source_spec) == SourceSpectrum
        assert_quantity_allclose(source_spec.waveset[0], 300 * u.nm)
        assert_quantity_allclose(source_spec(source_spec.waveset[0], flux_unit='FLAM'), 4.4279706e-16*self.funit)
        assert_quantity_allclose(source_spec.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(source_spec(source_spec.waveset[-1], flux_unit='FLAM'), 4.7801348e-16*self.funit)
