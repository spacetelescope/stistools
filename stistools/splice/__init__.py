from .splice import (
    splice,
    read_spectrum,
    find_overlap,
    merge_overlap,
    concatenate_sections,
    nearest_index,
    __doc__ as __doc_splice__)
from .calculate_sensitivity import (
    calculate_sensitivity,
    eval_tds,
    extraction_height_correction,
    aperture_throughput,
    __doc__ as __doc_calculate_sensitivity__)

__doc__ = __doc_splice__ + '\n\n' + '-'*50 + '\n\n' + __doc_calculate_sensitivity__
