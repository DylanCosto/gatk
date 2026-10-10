import logging
from bisect import bisect_left
from itertools import chain
import vcf
from typing import List, Tuple

_logger = logging.getLogger(__name__)


def read_sample_segments_and_calls(intervals_vcf: str,
                                   clustered_vcf: str,
                                   sample_name: str,
                                   contig: str) -> List[Tuple[int, int, int]]:
    """Get (copy number, start index, inclusive end index) for clustered segments.

    Modeling intervals must have nondecreasing starts and ends. Clustered
    segments may overlap or be nested, so their end indices need not increase.
    """
    intervals = vcf.Reader(filename=intervals_vcf)
    segments = vcf.Reader(filename=clustered_vcf)

    try:
        segments_iter = iter(segments.fetch(contig))
    except ValueError:
        return []
    first_segment = next(segments_iter, None)
    if first_segment is None:
        return []

    # Keep only coordinates: retaining parsed VCF records also retains every
    # sample's FORMAT data. Each interval is parsed once per contig.
    interval_starts: List[int] = []
    interval_ends: List[int] = []
    for interval in intervals.fetch(contig):
        start = interval.POS
        end = try_getting_info_attribute(interval, 'END')
        if interval_starts and (start < interval_starts[-1] or end < interval_ends[-1]):
            raise ValueError("Modeling interval starts and ends must be sorted on contig " + contig)
        interval_starts.append(start)
        interval_ends.append(end)
    if not interval_starts:
        raise StopIteration("No modeling intervals on contig " + contig)

    path: List[Tuple[int, int, int]] = []
    for segment in chain((first_segment,), segments_iter):
        segment_start_index = bisect_left(interval_starts, segment.POS)
        if segment_start_index == len(interval_starts):
            raise StopIteration("Ran out of intervals with unmatched segments remaining")
        segment_end_index = bisect_left(interval_ends, try_getting_info_attribute(segment, 'END'))
        if segment_end_index == len(interval_ends):
            # Preserve the previous behavior for a segment extending past the
            # last modeled interval: use the last available interval.
            _logger.warning("Ran out of intervals with segment end unmatched")
            segment_end_index -= 1
        if segment_end_index < segment_start_index:
            _logger.warning("Sample %s contains segment at %s:%s with end index less than start index",
                            sample_name, contig, segment.POS)
        segment_copy_number = try_getting_format_attribute(segment, sample_name, 'CN')
        path.append((segment_copy_number, segment_start_index, segment_end_index))
    return path


def try_getting_info_attribute(record,
                               attribute: str) -> int:
    try:
        value = record.INFO[attribute]
    except AttributeError:
        print('No {} field for record at position:{}'.format(attribute, record.POS))
    else:
        return value


def try_getting_format_attribute(record,
                                 sample_name: str,
                                 attribute: str) -> int:
    try:
        value = record.genotype(sample_name)[attribute]
    except AttributeError:
        print('No {} field for {} intervals at position:{}'.format(attribute, sample_name, record.POS))
    else:
        return value
