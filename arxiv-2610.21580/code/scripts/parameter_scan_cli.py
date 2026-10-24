#!/usr/bin/env python3
"""在缺少 matplotlib 时运行参数扫描的非绘图阶段。"""

from __future__ import annotations

import runpy
import sys
import types
from pathlib import Path


UPSTREAM = Path(r"D:\Project\harnessing-quantum-fork\submission\quantum.harness\scripts\parameter_scan.py")


def main() -> int:
    if not UPSTREAM.is_file():
        raise SystemExit("未找到参数扫描脚本")
    try:
        __import__("matplotlib")
    except ModuleNotFoundError:
        if len(sys.argv) > 1 and sys.argv[1] == "plot":
            raise SystemExit("plot 阶段需要 matplotlib；其余阶段可正常运行")
        matplotlib = types.ModuleType("matplotlib")
        matplotlib.use = lambda *_args, **_kwargs: None
        pyplot = types.ModuleType("matplotlib.pyplot")
        sys.modules["matplotlib"] = matplotlib
        sys.modules["matplotlib.pyplot"] = pyplot
    runpy.run_path(str(UPSTREAM), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
