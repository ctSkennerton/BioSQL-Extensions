#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `biosqlx` package."""

import os
import unittest
from io import StringIO
from click.testing import CliRunner

from Bio import SeqIO

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

    def test_export_sequence_qv_feat_prot(self):
        """Export from qualifier and value as translated proteins."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-q', 'gene', '-v', 'dnaA', '-o', 'feat-prot'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(3, count)

    def test_export_sequence_qv_feat_nucl(self):
        """Export from qualifier and value as untranslated proteins."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-q', 'gene', '-v', 'dnaA', '-o', 'feat-nucl'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(3, count)

    def test_export_sequence_qv_fasta(self):
        """Export from qualifier and value as fasta."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-q', 'gene', '-v', 'dnaA', '-o', 'fasta'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            self.assertEqual(True, rec.id in ['NC_004347.2', 'NC_002939.5', 'NC_000913.3'])
            count += 1

        f.close()
        self.assertEqual(3, count)

    def test_export_sequence_qv_genbank(self):
        """Export from qualifier and value as fasta."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-q', 'gene', '-v', 'dnaA', '-o', 'gb'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'genbank'):
            self.assertEqual(True, rec.id in ['NC_004347.2', 'NC_002939.5', 'NC_000913.3'])
            count += 1

        f.close()
        self.assertEqual(3, count)

    def test_export_sequence_qv_taxonomy_feat_prot(self):
        """Export from taxonomy, qualifier and value as translated proteins."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-t', 'Geobacter', '-q', 'gene', '-v', 'dnaA', '-o', 'feat-prot'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(1, count)

    def test_export_sequence_qv_taxonomy_feat_nucl(self):
        """Export from taxonomy, qualifier and value as untranslated proteins."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-t', 'Geobacter', '-q', 'gene', '-v', 'dnaA', '-o', 'feat-nucl'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(1, count)

    def test_export_sequence_qv_taxonomy_fasta(self):
        """Export from taxonomy, qalifier and value as fasta."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-t', 'Geobacter', '-q', 'gene', '-v', 'dnaA', '-o', 'fasta'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            self.assertEqual(rec.id, 'NC_002939.5')
            count += 1

        f.close()
        self.assertEqual(1, count)

    def test_export_sequence_qv_taxonomy_genbank(self):
        """Export from taxonomy, qualifier and value as genbank."""
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + ['-t', 'Geobacter', '-q', 'gene', '-v', 'dnaA', '-o', 'gb'])
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'genbank'):
            self.assertEqual(rec.id, 'NC_002939.5')
            count += 1

        f.close()
        self.assertEqual(1, count)

    def test_export_sequence_qv_namespace_feat_prot(self):
        """Export from namespace, qualifier and value as translated proteins."""
        fixture_params = ['-D', 'gammas', '-q', 'gene', '-v', 'dnaA', '-o', 'feat-prot']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)

        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(2, count)

    def test_export_sequence_qv_namespace_feat_nucl(self):
        """Export from namespace, qualifier and value as untranslated proteins."""
        fixture_params = ['-D', 'gammas', '-q', 'gene', '-v', 'dnaA', '-o', 'feat-nucl']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(2, count)

    def test_export_sequence_qv_namespace_fasta(self):
        """Export from namespace, qalifier and value as fasta."""
        fixture_params = ['-D', 'gammas', '-q', 'gene', '-v', 'dnaA', '-o', 'fasta']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            self.assertEqual(True, rec.id in ['NC_004347.2', 'NC_000913.3'])
            count += 1

        f.close()
        self.assertEqual(2, count)

    def test_export_sequence_qv_namespace_genbank(self):
        """Export from taxonomy, qualifier and value as genbank."""
        fixture_params = ['-D', 'gammas', '-q', 'gene', '-v', 'dnaA', '-o', 'gb']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'genbank'):
            self.assertEqual(True, rec.id in ['NC_004347.2', 'NC_000913.3'])
            count += 1

        f.close()
        self.assertEqual(2, count)

    def test_export_sequence_namespace_feat_prot(self):
        """Export from namespace as translated proteins."""
        fixture_params = ['-D', 'deltas', '-o', 'feat-prot']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)

        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(3430, count)

    def test_export_sequence_namespace_feat_nucl(self):
        """Export from namespace as untranslated proteins."""
        fixture_params = ['-D', 'deltas', '-o', 'feat-nucl']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            count += 1

        f.close()
        self.assertEqual(3485, count)

    def test_export_sequence_namespace_fasta(self):
        """Export from namespace as fasta."""
        fixture_params = ['-D', 'deltas', '-o', 'fasta']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'fasta'):
            self.assertEqual('NC_002939.5', rec.id)
            count += 1

        f.close()
        self.assertEqual(1, count)

    def test_export_sequence_namespace_genbank(self):
        """Export from taxonomy, qualifier and value as genbank."""
        fixture_params = ['-D', 'deltas', '-o', 'gb']
        runner = CliRunner()
        result = runner.invoke(cli.main, self.database_connection_params + self.common_params + fixture_params)
        self.assertEqual(result.exit_code, 0)
        f = StringIO(result.output)
        count = 0
        for rec in SeqIO.parse(f, 'genbank'):
            self.assertEqual('NC_002939.5', rec.id)
            count += 1

        f.close()
        self.assertEqual(1, count)
