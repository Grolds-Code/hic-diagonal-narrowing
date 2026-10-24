"""Inspect available .hic normalization metadata and reader capabilities."""

import modal

app = modal.App("hic-diagonal-narrowing-norm-capability")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install(
        "libcurl4-openssl-dev",
        "zlib1g-dev",
    )
    .pip_install("hic-straw==1.3.1")
)

URL = (
    "https://huggingface.co/buckets/FatimaInstitute/"
    "3d-genomics-prediction/resolve/HCT116/"
    "ENCFF573OPJ.hic?download=true"
)


@app.function(image=image, cpu=1, memory=1024, timeout=240)
def inspect_normalization_capabilities():
    import hicstraw

    hic = hicstraw.HiCFile(URL)

    resolutions = list(hic.getResolutions())

    if 500 not in resolutions:
        raise RuntimeError("500 bp resolution is unavailable")

    public_methods = [
        name for name in dir(hic)
        if not name.startswith("_")
    ]

    norm_methods = [
        name for name in public_methods
        if "norm" in name.lower()
    ]

    module_norm_methods = [
        name for name in dir(hicstraw)
        if "norm" in name.lower()
    ]

    result = {
        "genome": hic.getGenomeID(),
        "resolution_500bp": "AVAILABLE",
        "hicfile_methods": public_methods,
        "normalization_related_methods": norm_methods,
        "module_normalization_methods": module_norm_methods,
        "normalization_types": "NOT_EXPOSED",
        "direct_vector_access": "NOT_CONFIRMED",
        "method_documentation": {},
    }

    method = getattr(
        hic, "getNormalizationTypes", None
    )

    if callable(method):
        try:
            result["normalization_types"] = [
                str(value) for value in method()
            ]
        except Exception as exc:
            result["normalization_types"] = (
                f"QUERY_FAILED: {type(exc).__name__}: {exc}"
            )

    for name in norm_methods:
        member = getattr(hic, name)
        result["method_documentation"][name] = (
            str(getattr(member, "__doc__", None))[:1500]
        )

    if any(
        name in norm_methods
        for name in (
            "getNormVector",
            "getNormalizationVector",
        )
    ):
        result["direct_vector_access"] = (
            "POTENTIALLY_AVAILABLE_NEEDS_SIGNATURE_CHECK"
        )

    return result


@app.local_entrypoint()
def main():
    result = inspect_normalization_capabilities.remote()

    print("\nHCT116 .HIC NORMALIZATION CAPABILITY PROBE")
    print("=" * 65)

    for key, value in result.items():
        print(f"\n{key}: {value}")

    print("\nRead-only metadata inspection: PASS")
    print("Bin validity: NOT YET DETERMINED")
