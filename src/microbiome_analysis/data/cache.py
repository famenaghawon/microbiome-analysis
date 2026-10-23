"""Cache + fetch for real HMP data. Honest: no fake URLs."""
import hashlib
import os
import sys
from pathlib import Path

import requests

CACHE_DIR = Path(os.getenv("MB_CACHE", "data/cache"))
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, expected_sha256: str | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if expected_sha256 and sha256(dest) != expected_sha256:
            print(f"[cache] checksum mismatch on {dest}; refetching", file=sys.stderr)
        else:
            return dest
    print(f"[cache] downloading {url} -> {dest}", file=sys.stderr)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            f.writelines(r.iter_content(chunk_size=1 << 20))
        tmp.replace(dest)
    if expected_sha256 and sha256(dest) != expected_sha256:
        raise ValueError(f"Checksum mismatch after fetch: {dest}")
    return dest
