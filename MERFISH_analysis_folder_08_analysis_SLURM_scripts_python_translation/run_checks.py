from __future__ import annotations

import compileall
import pathlib
import subprocess
import sys

root = pathlib.Path(__file__).resolve().parent
translation = root / "python_translation"
compile_ok = compileall.compile_dir(str(translation), force=True, quiet=1)
print(f"compile_ok={compile_ok}")
if not compile_ok:
    sys.exit(1)
result = subprocess.run([sys.executable, "-m", "pytest", str(root / "tests"), "-q"], text=True, capture_output=True)
print(result.stdout)
print(result.stderr)
if result.returncode != 0:
    sys.exit(result.returncode)
