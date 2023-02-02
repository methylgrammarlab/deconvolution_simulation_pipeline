# ###################################################
#
# MIT License
#
# Copyright (c) 2022 irene unterman and ben berman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###################################################

import pytest
from scripts.epistate_simulator import generate_atlas

def test_epistate_atlas():
    config = {'theta_high': [0.8, 0.9, 0.7, 0.8, 0.8], 'theta_low': 0.2, 'lambda_high': 1, 'lambda_low': 0,
              'm_per_region': 5, 'regions_per_t': 1, 't': 4}
    h, l, t = generate_atlas(config)
    assert len(h) == config["regions_per_t"]*config['t']
    assert len(l) == config["regions_per_t"]*config['t']
    assert len(t) == config["regions_per_t"]*config['t']
    assert len(h[0]) == len(l[0])
    assert len(h[0]) == config["m_per_region"]
