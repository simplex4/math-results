#!/usr/bin/env python3
"""Run the certified d=5 interval verifier in independent parallel shards."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
from pathlib import Path
import subprocess


def run_one(binary: Path, log_dir: Path, index: int, count: int,
            max_boxes: int, max_depth: int) -> tuple[int, int]:
    log = log_dir / f"shard-{index:04d}-of-{count:04d}.log"
    command = [str(binary), "--shard-index", str(index),
               "--shard-count", str(count), "--max-boxes", str(max_boxes),
               "--max-depth", str(max_depth)]
    with log.open("w", encoding="utf-8") as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT,
                                check=False)
    return index, result.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("./makeev_d5_interval"))
    parser.add_argument("--shards", type=int, default=256)
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--max-boxes", type=int, default=0,
                        help="per-shard limit; zero means unlimited")
    parser.add_argument("--max-depth", type=int, default=400)
    parser.add_argument("--log-dir", type=Path, default=Path("interval-logs"))
    args = parser.parse_args()
    if args.shards < 1 or args.shards & (args.shards - 1):
        parser.error("--shards must be a power of two")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary not found: {binary}")
    args.log_dir.mkdir(parents=True, exist_ok=True)

    failed: list[int] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, binary, args.log_dir, i, args.shards,
                               args.max_boxes, args.max_depth)
                   for i in range(args.shards)]
        for future in concurrent.futures.as_completed(futures):
            index, code = future.result()
            print(f"shard {index}/{args.shards}: exit {code}", flush=True)
            if code:
                failed.append(index)

    if failed:
        failed.sort()
        print("INCOMPLETE shards:", " ".join(map(str, failed)))
        raise SystemExit(2)
    print("Certificate verified: every shard completed.")


if __name__ == "__main__":
    main()
