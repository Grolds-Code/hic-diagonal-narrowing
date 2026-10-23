"""Modal infrastructure for Hi-C Diagonal Narrowing research."""

from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = modal.Image.debian_slim(
    python_version="3.12"
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=512,
    timeout=120,
)
def initialize_storage():
    """Initialize persistent cloud research directories."""

    root = Path("/research")

    directories = [
        "raw",
        "matrices",
        "experiments",
        "results",
        "figures",
        "metadata",
    ]

    for directory in directories:
        (root / directory).mkdir(
            parents=True,
            exist_ok=True,
        )

    # Make writes visible to subsequent Modal functions.
    volume.commit()

    return {
        "volume": "hic-diagonal-narrowing-data",
        "directories": sorted(
            p.name for p in root.iterdir()
            if p.is_dir()
        ),
    }


@app.local_entrypoint()
def main():
    result = initialize_storage.remote()

    print("\nCLOUD RESEARCH STORAGE")
    print("=" * 40)
    print("Volume:", result["volume"])

    for directory in result["directories"]:
        print("  /research/" + directory)

    print("\nRemote initialization: PASS")
