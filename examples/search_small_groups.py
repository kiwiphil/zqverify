#!/usr/bin/env python3
import json
from pathlib import Path
from zqverify.search import run_catalog
out=run_catalog()
path=Path(__file__).resolve().parents[1]/'reports'/'small_group_search.json'
path.write_text(json.dumps(out,indent=2,sort_keys=True))
print(json.dumps(out,indent=2,sort_keys=True))
