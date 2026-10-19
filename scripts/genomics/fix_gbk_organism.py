#!/usr/bin/env python3

import sys
import re
from pathlib import Path

antismash_dir = Path(__file__).resolve().parents[2] / "derivedData/antismash"

for gbk in antismash_dir.rglob("*.region*.gbk"):
    sample = gbk.parent.name
    text = gbk.read_text()
    text = re.sub(r'^(SOURCE\s+)\.$', rf'\g<1>{sample}', text, flags=re.MULTILINE)
    text = re.sub(r'^(  ORGANISM\s+)\.$', rf'\g<1>{sample}', text, flags=re.MULTILINE)
    gbk.write_text(text)

print(f"Patched {sum(1 for _ in antismash_dir.rglob('*.region*.gbk'))} GBK files")
