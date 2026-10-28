"""Read-only audit of original HCT116 Hi-C contact-count values."""

import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-integer-audit")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = (
    "hct116_chr10_20000000_21000000_"
    "500bp_midpoint_32kb"
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=1024,
    timeout=180,
)
def audit_original_counts():
    import numpy as np

    root = Path("/research")

    source_path = root / "matrices" / f"{SOURCE}.npz"
    metadata_path = root / "metadata" / f"{SOURCE}.json"

    source_bytes = source_path.read_bytes()
    metadata = json.loads(metadata_path.read_text())

    digest = hashlib.sha256(source_bytes).hexdigest()

    if digest != metadata["sha256"]:
        raise RuntimeError("Source SHA-256 mismatch")

    with np.load(
        io.BytesIO(source_bytes),
        allow_pickle=False,
    ) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if not (x.shape == y.shape == counts.shape):
        raise RuntimeError("Source array shape mismatch")

    if not np.all(np.isfinite(counts)):
        raise RuntimeError("Nonfinite contact values")

    if np.any(counts < 0):
        raise RuntimeError("Negative contact values")

    integral = counts == np.rint(counts)
    fractional_indices = np.flatnonzero(~integral)

    positive = counts > 0
    off_diagonal = x < y

    if counts.size != metadata["selected_records"]:
        raise RuntimeError("Record-count mismatch")

    result = {
        "source_sha256_verified": True,
        "dtype": str(counts.dtype),
        "total_records": int(counts.size),
        "off_diagonal_records": int(
            np.count_nonzero(off_diagonal)
        ),
        "diagonal_records": int(
            np.count_nonzero(~off_diagonal)
        ),
        "positive_records": int(
            np.count_nonzero(positive)
        ),
        "zero_valued_records": int(
            np.count_nonzero(counts == 0)
        ),
        "integer_valued_records": int(
            np.count_nonzero(integral)
        ),
        "fractional_records": int(
            fractional_indices.size
        ),
        "minimum_count": float(counts.min()),
        "maximum_count": float(counts.max()),
        "total_contact_signal": float(counts.sum()),
        "all_integer_valued": bool(np.all(integral)),
        "within_exact_float_integer_range": bool(
            np.all(counts <= 2**53 - 1)
        ),
        "fractional_examples": [
            {
                "bin_x_bp": int(x[i]),
                "bin_y_bp": int(y[i]),
                "count": float(counts[i]),
            }
            for i in fractional_indices[:10]
        ],
    }

    return result


@app.local_entrypoint()
def main():
    result = audit_original_counts.remote()

    print("\nHCT116 ORIGINAL CONTACT-COUNT AUDIT")
    print("=" * 65)
    print(json.dumps(result, indent=2))

    if result["all_integer_valued"]:
        print("\nPASS: All stored contact values are integer-valued")
    else:
        print(
            "\nWARNING: Fractional contact values detected; "
            "Poisson count assumptions require further investigation"
        )

    print("\nAUDIT COMPLETE — SOURCE DATA UNCHANGED")
