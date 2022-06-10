import os
import tempfile, shutil
import pytest

import astropy.units as u

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

