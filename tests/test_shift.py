import numpy as np

from donorconf.shift import donor_permutation_test, shift_screen


def test_exact_test_minimum_p_and_flags():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(6, 2))
    B = rng.normal(size=(3, 2))
    r = donor_permutation_test(A, B)
    assert r["exact"] and abs(r["min_attainable_p"] - 1 / 84) < 1e-12  # C(9,3)=84
    assert 0 < r["p_value"] <= 1


def test_screen_warns_when_it_cannot_reject():
    rng = np.random.default_rng(0)
    r = shift_screen(rng.normal(size=(3, 2)), rng.normal(size=(2, 2)), level=0.05)
    assert r["can_reject"] is False  # C(5,2)=10 -> min p = 0.1


def test_type_one_error_controlled_under_exchangeability():
    rng = np.random.default_rng(3)
    R, rej = 400, 0
    for _ in range(R):
        A = rng.normal(size=(12, 3))
        B = rng.normal(size=(8, 3))
        rej += donor_permutation_test(A, B, n_perm=499, max_exact=0, seed=int(rng.integers(1e9)))["p_value"] <= 0.05
    assert rej / R <= 0.05 + 3 * np.sqrt(0.05 * 0.95 / R)


def test_power_against_a_large_shift():
    rng = np.random.default_rng(4)
    R, rej = 100, 0
    for _ in range(R):
        A = rng.normal(size=(12, 3))
        B = rng.normal(size=(8, 3)) + 2.0
        rej += donor_permutation_test(A, B, n_perm=499, max_exact=0, seed=int(rng.integers(1e9)))["p_value"] <= 0.05
    assert rej / R > 0.9
