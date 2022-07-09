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
