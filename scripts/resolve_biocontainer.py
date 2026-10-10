#!/usr/bin/env python3
"""Print docker://quay.io/biocontainers/<tool>:<version>--<build> for the newest build of a version.
Usage: resolve_biocontainer.py <tool> <version>
(QUAY_API env var overrides the API root; used for testing.)"""
import json, os, sys, urllib.parse, urllib.request

def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    tool, ver = sys.argv[1], sys.argv[2]
    api = os.environ.get("QUAY_API", "https://quay.io/api/v1").rstrip("/")
    tags, page = [], 1
    while page <= 15:
        url = (f"{api}/repository/biocontainers/{tool}/tag/?limit=100&page={page}"
               f"&onlyActiveTags=true&filter_tag_name=like:{urllib.parse.quote(ver)}")
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
        tags += d.get("tags", [])
        if not d.get("has_additional"):
            break
        page += 1
    cand = [t for t in tags if t.get("name", "").startswith(ver + "--")]
    if not cand:
        names = ", ".join(sorted({t.get("name", "") for t in tags})[:15])
        sys.exit(f"no tag starting with '{ver}--' for {tool}. Examples seen: {names}")
    cand.sort(key=lambda t: t.get("start_ts", 0), reverse=True)
    print(f"docker://quay.io/biocontainers/{tool}:{cand[0]['name']}")

if __name__ == "__main__":
    main()
