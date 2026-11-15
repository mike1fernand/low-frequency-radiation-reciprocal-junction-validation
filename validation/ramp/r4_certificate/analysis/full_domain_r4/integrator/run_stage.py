#!/usr/bin/env python3
"""One bounded stage. Repeat the identical command to resume."""
import argparse
import json
from pathlib import Path
from compact_controller import integrate_stage,summary

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--settings',required=True,help='JSON settings file')
    p.add_argument('--output',required=True)
    p.add_argument('--max-seconds',type=float,default=20)
    a=p.parse_args()
    print(json.dumps(summary(integrate_stage(json.loads(Path(a.settings).read_text()),a.output,a.max_seconds))))
