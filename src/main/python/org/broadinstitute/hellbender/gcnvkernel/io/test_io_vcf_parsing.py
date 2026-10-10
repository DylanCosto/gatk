from unittest import TestCase
from unittest.mock import patch
from types import SimpleNamespace
import os
import gcnvkernel.io.io_vcf_parsing as io


class test_io_vcf_parsing(TestCase):
    def test_read_sample_segments_and_calls(self):

        current_dir = os.getcwd()
        #for GATK PythonUnitTestRunner/Java tests
        test_sub_dir = current_dir + "/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        # for running in IntelliJ/Python tests
        # test_sub_dir = current_dir + "/../../../../../../../../src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"

        clustered_vcf = test_sub_dir + 'clustered.1000G.vcf.gz'
        pesky_intervals_vcf = test_sub_dir + 'genotyped-intervals-HG00099.mapped.ILLUMINA.bwa.GBR.exome.20130415.bam.cram.vcf.gz'
        pesky_sample_name = 'HG00099'
        contig = '14'
        debug_path = io.read_sample_segments_and_calls(pesky_intervals_vcf, clustered_vcf, pesky_sample_name, contig)
        self.assertTrue(len(debug_path) == 12)  #should match number of chr 14 lines in clustered VCF (12)

        clustered_vcf = current_dir + '/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/clustering/threeSamples.vcf.gz'
        intervals_vcf = test_sub_dir + 'intervals_output_SAMPLE_000.vcf.gz'
        sample_name = 'SAMPLE_000'

        contig = "1"
        path = io.read_sample_segments_and_calls(intervals_vcf, clustered_vcf, sample_name, contig)
        # no segments on contig 1
        self.assertTrue(len(path) == 0)

        contig = "2"
        path = io.read_sample_segments_and_calls(intervals_vcf, clustered_vcf, sample_name, contig)
        self.assertTrue(len(path) == 2)

        sample1_intervals_vcf = test_sub_dir + 'intervals_output_SAMPLE_001.vcf.gz'
        sample1_name = 'SAMPLE_001'
        contig = "2"
        path1 = io.read_sample_segments_and_calls(sample1_intervals_vcf, clustered_vcf, sample1_name, contig)
        # all samples should have the same number of intervals
        self.assertTrue(len(path) == len(path1))


class TestIndexedSegmentLookup(TestCase):
    @staticmethod
    def record(start, end, copy_number=2):
        return SimpleNamespace(POS=start, INFO={'END': end},
                               genotype=lambda sample: {'CN': copy_number})

    def read(self, intervals, segments):
        self.records_read = 0

        def reader(filename):
            records = intervals if filename == 'intervals' else segments

            def fetch(contig):
                if records is None:
                    raise ValueError('Missing contig')

                def iterator():
                    for record in records:
                        self.records_read += 1
                        yield record
                return iterator()
            return SimpleNamespace(fetch=fetch)

        with patch.object(io.vcf, 'Reader', side_effect=reader):
            return io.read_sample_segments_and_calls('intervals', 'segments', 'sample', '1')

    def test_overlapping_and_nested_segments(self):
        intervals = [self.record(100 * i + 1, 100 * i + 100) for i in range(5)]
        segments = [self.record(1, 500, 1), self.record(1, 100, 2),
                    self.record(101, 300, 3), self.record(301, 500, 4)]
        self.assertEqual(self.read(intervals, segments),
                         [(1, 0, 4), (2, 0, 0), (3, 1, 2), (4, 3, 4)])

    def test_gaps_and_boundaries(self):
        intervals = [self.record(1, 100), self.record(201, 300), self.record(401, 500)]
        segments = [self.record(1, 100), self.record(50, 250),
                    self.record(201, 300), self.record(401, 600)]
        self.assertEqual(self.read(intervals, segments),
                         [(2, 0, 0), (2, 1, 1), (2, 1, 1), (2, 2, 2)])

    def test_missing_or_empty_segments(self):
        intervals = [self.record(1, 100)]
        self.assertEqual(self.read(intervals, None), [])
        self.assertEqual(self.read(intervals, []), [])
        self.assertEqual(self.records_read, 0)

    def test_missing_intervals_and_out_of_range_start(self):
        with self.assertRaises(StopIteration):
            self.read([], [self.record(1, 100)])
        with self.assertRaises(StopIteration):
            self.read([self.record(1, 100)], [self.record(101, 200)])

    def test_reject_nonmonotonic_modeling_intervals(self):
        for intervals in ([self.record(101, 200), self.record(1, 300)],
                          [self.record(1, 300), self.record(101, 200)]):
            with self.subTest(intervals=intervals), self.assertRaises(ValueError):
                self.read(intervals, [self.record(1, 100)])

    def test_each_interval_and_segment_is_parsed_once(self):
        n = 200
        intervals = [self.record(100 * i + 1, 100 * i + 100) for i in range(n)]
        segments = list(intervals)
        self.assertEqual(self.read(intervals, segments), [(2, i, i) for i in range(n)])
        self.assertEqual(self.records_read, len(intervals) + len(segments))
