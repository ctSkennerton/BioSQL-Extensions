#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `biosqlx` package."""

import os
import unittest
from io import StringIO
from click.testing import CliRunner

from biosqlx import biosqlx
from biosqlx import cli
from . import connection_parameters

class TestExportSequence(unittest.TestCase):
    """Tests for `biosqlx` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        testdb, dbdriver, dbuser, dbpassword, dbhost = connection_parameters()
        self.database_connection_params = ['-d', testdb, '-r', dbdriver]
        self.common_params = ['export', 'sequence']

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_taxonomy_feat_prot(self):
        """Export from taxonomy as protein features."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['--taxonomy', 'Geobacter', '-o', 'feat-prot'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual('>', result.output[0])

    def test_taxonomy_fasta(self):
        """Export from taxonomy as fasta."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['--taxonomy', 'Geobacter'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual('>', result.output[0])

    def test_export_sequence_feat_nucl(self):
        """Export from taxonomy as nucleotide features."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['--taxonomy', 'Geobacter', '-o', 'feat-nucl'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual('>', result.output[0])

    def test_export_sequence_genbank(self):
        """Export from taxonomy as genbank."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['--taxonomy', 'Geobacter', '-o', 'gb'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual('L', result.output[0])

    def test_export_sequence_csv(self):
        """Export from taxonomy as csv."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['--taxonomy', 'Geobacter', '-o', 'csv'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        line = next(f)
        f.close()
        self.assertEqual(True, ',' in line)

