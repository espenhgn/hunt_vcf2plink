"""Fetch UCSC liftOver and the hg19ToHg38 chain file (Snakemake script, rule setup_liftover).

liftOver is the official UCSC build for the platform running the workflow
(LIFTOVER_PLATFORMS), unless LIFTOVER_BINARY gives a URL or local path. The chain
file (CHAIN_FILE) is a URL or local path. Both are checked against md5 sums when
given. UCSC liftOver and chain files are free for non-commercial use
(https://genome.ucsc.edu/license/).
"""

import hashlib
import platform
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def fetch(src, dest, expected_md5=""):
    """Download src (URL) or copy it (local path) to dest; check md5 if given."""
    tmp = dest.with_name(dest.name + ".part")
    if "://" in src:
        print(f"downloading {src}")
        urllib.request.urlretrieve(src, tmp)
    else:
        print(f"copying {src}")
        shutil.copyfile(src, tmp)
    if expected_md5:
        h = hashlib.md5()
        with open(tmp, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        if h.hexdigest() != expected_md5:
            tmp.unlink()
            sys.exit(f"md5 mismatch for {src}: {h.hexdigest()}, expected {expected_md5}")
        print(f"md5 ok ({expected_md5})")
    tmp.replace(dest)


def main(params, output):
    liftover, liftover_md5 = params.liftover, params.liftover_md5
    if not liftover:
        key = f"{platform.system()}-{platform.machine()}"
        if key not in params.platforms:
            sys.exit(f"no UCSC liftOver build configured for {key} (LIFTOVER_PLATFORMS); "
                     "set LIFTOVER_BINARY")
        build = params.platforms[key]
        liftover = f"{params.url_base}/{build['dir']}/liftOver"
        liftover_md5 = build.get("md5") or ""

    liftover_bin = Path(output.liftover)
    liftover_bin.parent.mkdir(parents=True, exist_ok=True)
    fetch(liftover, liftover_bin, liftover_md5)
    liftover_bin.chmod(0o755)

    # liftOver without arguments prints its usage and exits non-zero
    res = subprocess.run([str(liftover_bin)], capture_output=True, text=True, check=False)
    if "liftOver" not in res.stdout + res.stderr:
        sys.exit(f"{liftover_bin} does not run on this system:\n{res.stdout}{res.stderr}")

    fetch(params.chain, Path(output.chain), params.chain_md5)


if __name__ == "__main__":
    with open(snakemake.log[0], "w") as log:  # noqa: F821
        sys.stdout = sys.stderr = log
        main(snakemake.params, snakemake.output)  # noqa: F821
