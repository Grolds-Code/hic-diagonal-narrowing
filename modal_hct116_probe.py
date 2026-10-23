"""Test remote HCT116 Hi-C file access from Modal."""

import urllib.request

import modal


app = modal.App("hic-diagonal-narrowing-probe")

image = modal.Image.debian_slim(
    python_version="3.12"
)

URL = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)


@app.function(
    image=image,
    cpu=1,
    memory=512,
    timeout=90,
)
def check_hct116_access():
    """Retrieve only the first few bytes of the remote file."""

    request = urllib.request.Request(
        URL,
        headers={"Range": "bytes=0-15"},
    )

    with urllib.request.urlopen(
        request,
        timeout=30,
    ) as response:

        status = response.status
        content_range = response.headers.get(
            "Content-Range", "missing"
        )

        if status != 206:
            raise RuntimeError(
                f"HTTP {status}: byte-range request "
                "was not honored. Refusing full download."
            )

        header = response.read(16)

    if not header.startswith(b"HIC\x00"):
        raise RuntimeError(
            f"Unexpected Hi-C file signature: {header[:8]!r}"
        )

    return {
        "http_status": status,
        "content_range": content_range,
        "hic_signature": "VALID",
        "bytes_read": len(header),
    }


@app.local_entrypoint()
def main():
    result = check_hct116_access.remote()

    print("\nHCT116 REMOTE DATA ACCESS")
    print("=" * 45)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nREMOTE DATA PROBE: PASS")
