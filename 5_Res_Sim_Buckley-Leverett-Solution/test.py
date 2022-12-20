#!/usr/bin/env python

# Copyright 2018-2020 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import unittest
import nbconvert

with open("assignment19.ipynb") as f:
    exporter = nbconvert.PythonExporter()
    python_file, _ = exporter.from_file(f)


with open("assignment19.py", "w") as f:
    f.write(python_file)

from assignment19 import BuckleyLeverett

class TestSolution(unittest.TestCase):

    def setUp(self):

        self.inputs = {
                'reservoir': { 
                    'oil': {
                        'residual saturation': 0.2,
                        'corey-brooks exponent': 3.0,
                        'max relative permeability': 0.2,
                        },
                    'water': {
                        'critical saturation': 0.2,
                        'corey-brooks exponent': 3.0,
                        'max relative permeability': 1.0,
                        },
                    },
                'fluid': {
                    'oil': {
                        'viscosity': 1.0,
                    },
                    'water': {
                        'viscosity': 1.0,
                        },
                    },
                'initial conditions': {
                    'water saturation': 0.2,
                    },
                }

    def test_water_front_saturation(self):

        bl = BuckleyLeverett(self.inputs)

        Swf = bl.compute_saturation_front()

        self.assertAlmostEqual(Swf, 0.49999999999970984, 7)


if __name__ == '__main__':
    unittest.main()
