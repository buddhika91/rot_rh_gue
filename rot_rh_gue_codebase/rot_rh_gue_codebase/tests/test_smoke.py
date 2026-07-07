from rot_rh_gue.arithmetic import mangoldt_array
from rot_rh_gue.audit import run_frozen_audit


def test_mangoldt_small_values():
    lam = mangoldt_array(10)
    assert lam[2] > 0
    assert lam[4] == lam[2]
    assert lam[6] == 0


def test_tiny_audit_runs(tmp_path):
    prefix = tmp_path / "tiny"
    out = run_frozen_audit(
        N_list=[200, 300],
        jacobi_dims=[16],
        controls=1,
        seed=1,
        out_prefix=str(prefix),
        control_modes=["permuted"],
        oversample=4,
        verbose=False,
    )
    assert out["summary"]
    assert (tmp_path / "tiny_summary.csv").exists()
