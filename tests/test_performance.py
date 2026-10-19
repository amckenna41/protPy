###########################################################################
###############      protPy - Performance Test (Phase 8)        ###########
###########################################################################

import os
import random
import time
import unittest

import protpy


@unittest.skipUnless(os.environ.get("PROTPY_PERF"), "set PROTPY_PERF=1 to run the performance test")
class ProtpyPerformanceTests(unittest.TestCase):

    def test_all_descriptors_500_residues(self):
        random.seed(0)
        seq = "".join(random.choice(protpy.amino_acids) for _ in range(500))
        start = time.perf_counter()
        protpy.all_descriptors(seq, descriptor_set="all", _show_descriptor_list=False)
        self.assertLess(time.perf_counter() - start, 2.0)


if __name__ == "__main__":
    unittest.main()
