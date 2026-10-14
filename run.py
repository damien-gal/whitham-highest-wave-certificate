#!/usr/bin/env python3
"""run.py: compile the programs, run the stages and record the certificate.

Python standard library only. Every proof inequality, quadrature, lattice
refinement and scalar closure runs in the C programs. This script checks the
exact inputs against the sources, the exit code and success marker of every
stage, the coverage records in the logs, and the provenance: hashes of the
sources, inputs, binaries, logs and linked libraries. It never performs proof
arithmetic and never accepts a numerical sample as evidence.

Usage:
    python3 run.py                   compile, run all stages, write artifacts/
    python3 run.py --check-existing  re-validate the artifacts/ of a previous run
"""
import argparse
import fcntl
import hashlib
import json
import platform
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parent
ARTIFACTS = ROOT / 'artifacts'
CERTIFICATE = ARTIFACTS / 'certificate.json'
PAPER = ROOT.parent / 'paper/appendices/c-uniqueness-kernel-estimates.tex'

PRECISION_BITS = 128
WORKERS = 2                 # stages run two at a time, as does the build
STAGE_TIMEOUT = 7200        # seconds, far above the longest stage
GRID_SHA256 = 'b0d0c1783ca0700ff2cf65e520d2e59fe7513718103e00671181113bbe6f26d2'

# The files whose hashes the certificate records: the package root files and
# everything under the source and input directories.
HASHED_FILES = ['Makefile', 'precision.h', 'PROOF.md', 'import-manifest.json', 'run.py']
HASHED_DIRS = ['core', 'existence', 'convexity', 'uniqueness', 'inputs']

SCHEMA = 'whitham-unified-c-v1'
SCOPE = ('Existence, derivative/convexity and uniqueness, conditional on the documented '
         'analytic arguments and FLINT/Arb arithmetic base.')
ARITHMETIC = ('All proof arithmetic in C/Arb at 128 working bits and exact signed 128-bit '
              'lattice refinement. Python standard library orchestrates only.')
REMAINING = ['Independent human proof/code review', 'Approved manuscript exposition integration']


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

class Stage(NamedTuple):
    name: str
    command: list       # run from the package root
    exit_code: int      # expected exit code: 0 for a proof, 2 for a rejection test
    marker: str         # must appear in the log


STAGES = [
    # existence
    Stage('existence', ['build/defect', '0.05', 'compact'], 0, 'u0/F0 compact coverage: PASS'),
    Stage('T0', ['build/norm_t0'], 0, 'T0 compact coverage PASS'),
    Stage('supplemental_scalars', ['build/supplemental_scalars'], 0,
          'PASS supplemental scalar comparisons'),
    Stage('scalars', ['build/scalars'], 0, 'SCALAR CONDITIONS PASS'),
    Stage('evaluators', ['build/test_evaluators'], 0, 'regression checks: PASS'),
    # derivative and convexity
    Stage('derivative_cusp', ['build/derivative_cusp', '0.05'], 0, 'DERIVATIVE CUSP PASS'),
    Stage('derivative_compact', ['build/derivative_compact', '0.05'], 0,
          'DERIVATIVE COMPACT PASS'),
    Stage('local1', ['build/local_t1', '0.002'], 0, 'LOCAL T1 PASS'),
    Stage('local2', ['build/local_t2', '0.002'], 0, 'LOCAL T2 PASS'),
    Stage('operator_regression', ['build/operators', 'crosscheck'], 0,
          'global kernel crosscheck PASS'),
    Stage('T1', ['build/operators', '1', '0.002', 'cover', '0.998'], 0, 'OP1 COMPACT PASS'),
    Stage('T2', ['build/operators', '2', '0.002', 'cover', '0.998'], 0, 'OP2 COMPACT PASS'),
    # uniqueness
    Stage('sign', ['build/sign'], 0, 'PASS kernel-sign'),
    Stage('refinement_tests', ['build/test_refine'], 0, 'PASS 300'),
    Stage('kernel_tests', ['build/test_kernel'], 0, 'PASS ball conversion tests'),
    Stage('closure', ['build/closure'], 0, 'PASS scalar closure CONDITIONAL'),
    Stage('matrix', ['build/matrix', 'inputs/grid720.bin', '720', 'artifacts/moments.bin'], 0,
          'PASS: kernel moment'),
    # rejection tests: each program must refuse an impossible target
    Stage('reject_defect', ['build/defect', '0.05', 'compact', '--reject-test'], 2, 'UNRESOLVED'),
    Stage('reject_derivatives', ['build/derivative_compact', '0.05', '--reject-test'], 2,
          'UNRESOLVED'),
    Stage('reject_T1', ['build/operators', '1', '.1', 'reject'], 2, 'UNRESOLVED'),
    Stage('reject_T2', ['build/operators', '2', '.1', 'reject'], 2, 'UNRESOLVED'),
    Stage('reject_local1', ['build/local_t1', '0.5'], 2, 'UNRESOLVED'),
    Stage('reject_local2', ['build/local_t2', '0.003'], 2, 'UNRESOLVED'),
    Stage('reject_sign', ['build/sign', 'reject'], 2, 'UNRESOLVED'),
]

# The comparison reads the moment matrices that the matrix stage writes, so it
# runs after all the stages above have finished.
LATE_STAGES = [
    Stage('uniqueness', ['build/comparison', 'inputs/grid720.bin', 'artifacts/moments.bin'], 0,
          'PASS uniqueness full cusp and 719 compact cells'),
    Stage('reject_uniqueness',
          ['build/comparison', 'inputs/grid720.bin', 'artifacts/moments.bin', 'reject'], 2,
          'UNRESOLVED: cusp comparison target'),
]

STAGE_BY_NAME = {stage.name: stage for stage in STAGES + LATE_STAGES}

# Stages whose log must record a complete recursive cover: a binary tree with
# nodes = 2 leaves - 1, printed after the success marker.
RECURSIVE_COVERS = ['existence', 'T0', 'derivative_compact', 'T1', 'T2']


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    """Write the file atomically, so that a crash never leaves a half record."""
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def count_lines(pattern, text):
    return len(re.findall(pattern, text, re.M))


# ---------------------------------------------------------------------------
# Provenance hashes
# ---------------------------------------------------------------------------

def source_hashes():
    """Hash of every file the certificate depends on, keyed by relative path."""
    paths = [ROOT / name for name in HASHED_FILES]
    for directory in HASHED_DIRS:
        paths += [path for path in (ROOT / directory).rglob('*') if path.is_file()]
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(paths)}


def library_hashes():
    """The shared libraries the programs link against, read from ldd."""
    text = subprocess.check_output(['ldd', str(ROOT / 'build/comparison')], text=True)
    return {path: sha256(path) for path in re.findall(r'(?:=> )?(/\S+)\s+\(', text)}


def header_hashes():
    """The system headers the build included, read from the dependency files."""
    paths = set()
    for dependency_file in (ROOT / 'build').glob('*.d'):
        paths.update(re.findall(r'(?<!\S)(/[^\s\\]+\.h)(?=\s|:|$)', dependency_file.read_text()))
    return {path: sha256(path) for path in sorted(paths)}


# ---------------------------------------------------------------------------
# Input checks: the exact inputs of the proof must agree with the sources
# ---------------------------------------------------------------------------

def check_precision():
    text = (ROOT / 'precision.h').read_text()
    require(re.search(rf'^#define CAP_PREC {PRECISION_BITS}$', text, re.M),
            'Common precision must be 128')


def check_approximation():
    """core/basic_functions.c assigns exactly the decimals of inputs/approximation.json."""
    expected = json.loads((ROOT / 'inputs/approximation.json').read_text())['assignments']
    code = (ROOT / 'core/basic_functions.c').read_text()
    assigned = dict(re.findall(r'arb_set_str\(up->([^,]+),"([^"]+)",prec\);', code))
    require(assigned == expected,
            'Exact approximation manifest differs from core/basic_functions.c')


# Name of each scalar in convexity/closure.c and its key in inputs/closure.json.
CLOSURE_SCALARS = {'A': 'A', 'r': 'r', 'q0': 'q0', 'd0': 'delta0', 'd1': 'delta1',
                   'd2': 'delta2', 'a1': 'alpha1', 'a2': 'alpha2', 'af': 'alphaf'}


def check_closure():
    """convexity/closure.c declares the scalars of inputs/closure.json, under the fixed contract."""
    closure = json.loads((ROOT / 'inputs/closure.json').read_text())
    code = (ROOT / 'convexity/closure.c').read_text()
    declared = dict(re.findall(r'(\w+)\s*=\s*decimal\("([^" ]+)"\)', code))
    for c_name, key in CLOSURE_SCALARS.items():
        require(Decimal(declared[c_name]) == Decimal(closure[key]), 'Scalar input mismatch ' + key)
    contractions = all(Decimal(closure[key]) == Decimal('.998') for key in ['q0', 'q1', 'q2'])
    require(closure['operator_mode'] == 'raw' and contractions
            and Decimal(closure['curvature_budget']) == Decimal('.089'),
            'Unsupported scalar contract')


def check_kernel_sign():
    """uniqueness/kernel.c holds the thirteen coefficients of each root polynomial.

    Returns the manifest inputs/kernel-sign.json, whose coefficients and targets
    are then bound to the paper.
    """
    manifest = json.loads((ROOT / 'inputs/kernel-sign.json').read_text())
    code = (ROOT / 'uniqueness/kernel.c').read_text()
    for side in ['minus', 'plus']:
        pattern = r'arb_set_str\(' + side + r'\s*\+\s*(\d+),\s*"([^"]+)",\s*CAP_PREC\)'
        found = re.findall(pattern, code)
        indices = [int(index) for index, _ in found]
        values = [value for _, value in found]
        require(indices == list(range(13)) and values == manifest['coeffs'][side],
                'C root coefficients differ')
    return manifest


def check_grid():
    require(sha256(ROOT / 'inputs/grid720.bin') == GRID_SHA256, 'Exact grid input')


# One term "+ 1.234 x^{5}" of a root polynomial as printed in the paper.
TERM = r'([+-]?)\s*(\d+\.\d+)\s*x(?:\^\{?(\d+)\}?)?'

# The seven targets of the kernel-sign lemma, as written in its proof once the
# whitespace is removed; a second group is a power of ten.
TARGETS = {
    'epsilon': r'\\ep:=([0-9.]+)',
    'epsilon_star': r'\\ep\^\*:=([0-9.]+)\\cdot10\^\{(-?\d+)\}',
    'epsilon_K': r'\\ep_\\cK:=([0-9.]+)\\cdot10\^\{(-?\d+)\}',
    'r_minus': r'r_-:=([0-9.]+)',
    'R_minus': r'R_-:=([0-9.]+)',
    'r_plus': r'r_\+:=([0-9.]+)',
    'R_plus': r'R_\+:=([0-9.]+)',
}


def bind_paper(manifest):
    """Check the kernel-sign inputs against the paper and describe the binding.

    The coefficients of the two root polynomials and the seven targets of the
    lemma are read from the uniqueness appendix and compared, as exact decimals,
    with inputs/kernel-sign.json. Nothing else in the manuscript is a premise of
    the computation. The ancillary package distributed with the paper has no
    paper source; the manifest is then the only reference.
    """
    if not PAPER.exists():
        return {'path': str(PAPER),
                'status': 'absent (standalone package); exact kernel-sign coefficients and '
                          'targets taken from inputs/kernel-sign.json'}
    text = '\n'.join(line.split('%', 1)[0] for line in PAPER.read_text().splitlines())
    polynomials = [('minus', 'eq:gamma-', 'y^*_-(x)'), ('plus', 'eq:gamma+', 'y^*_+(x)')]
    for side, label, symbol in polynomials:
        block = text.split('\\label{' + label + '}', 1)[1].split('\\end{multline}', 1)[0]
        body = block.split(symbol + ':=', 1)[1].split('\\,', 1)[0].replace('\\\\', '')
        terms = re.findall(TERM, body)
        require(not re.sub(TERM, '', body).strip(), 'Unparsed paper polynomial')
        require([int(power or 1) for _, _, power in terms] == list(range(1, 14)),
                'Polynomial powers')
        coefficients = [(sign if sign == '-' else '') + value for sign, value, _ in terms]
        require(coefficients == manifest['coeffs'][side], 'Paper root coefficients differ')
    proof = text.split('\\label{L.signK}', 1)[1].split('\\end{proof}', 1)[0]
    proof = re.sub(r'\s+', '', proof)
    for name, pattern in TARGETS.items():
        match = re.search(pattern, proof)
        require(match is not None, 'Missing target ' + name)
        value = Decimal(match[1])
        if len(match.groups()) == 2:
            value *= Decimal(10) ** int(match[2])
        require(value == Decimal(manifest['targets'][name]), 'Changed target ' + name)
    return {'path': str(PAPER), 'sha256_at_run': sha256(PAPER),
            'binding': 'Parsed coefficient powers, exact decimal strings and seven exact targets. '
                       'Other manuscript text is not a computational premise.'}


def check_inputs():
    """Run every input check and return the record of the paper binding."""
    check_precision()
    check_approximation()
    check_closure()
    manifest = check_kernel_sign()
    check_grid()
    return bind_paper(manifest)


# ---------------------------------------------------------------------------
# Coverage: a log must record a complete cover, not only a success marker
# ---------------------------------------------------------------------------

def check_coverage(stage, text):
    if stage.name in RECURSIVE_COVERS:
        match = re.search(re.escape(stage.marker) + r' leaves=(\d+) nodes=(\d+)', text)
        require(match is not None, 'Missing coverage counts ' + stage.name)
        leaves, nodes = (int(group) for group in match.groups())
        require(leaves > 0 and nodes == 2 * leaves - 1, 'Incomplete recursive cover ' + stage.name)
        if stage.name in ['T1', 'T2']:
            require(count_lines(r'^OP[12] LEAF ', text) == leaves, 'Missing operator leaves')
    elif stage.name == 'uniqueness':
        cells = [int(index) for index in re.findall(r'^COMPARISON_CELL (\d+) ', text, re.M)]
        require(cells == list(range(719)), 'Missing comparison interval')
        require(count_lines(r'^BOUND ', text) == 720, 'Missing profile bounds')
    elif stage.name == 'sign':
        # A forest of 448 binary trees: 64 cusp cells and 384 compact cells.
        match = re.search(r'PASS kernel-sign.*leaves=(\d+) nodes=(\d+)', text)
        require(match is not None, 'Missing sign coverage')
        leaves, nodes = (int(group) for group in match.groups())
        require(nodes == 2 * leaves - 448 and count_lines(r'^SIGN_LEAF ', text) == leaves,
                'Incomplete sign forest')


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def run_stage(stage, artifacts=ARTIFACTS):
    """Run one stage, check its outcome and write its log and record."""
    binary = ROOT / stage.command[0]
    log_path = artifacts / (stage.name + '.log')
    print('START', stage.name, flush=True)
    start = time.monotonic()
    binary_hash = sha256(binary)
    with log_path.open('w') as log:
        process = subprocess.run(stage.command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                 timeout=STAGE_TIMEOUT)
    text = log_path.read_text()
    require(process.returncode == stage.exit_code and stage.marker in text,
            'Stage failed: ' + stage.name)
    require(binary_hash == sha256(binary), 'Binary changed during ' + stage.name)
    check_coverage(stage, text)
    result = {
        'status': 'PASS',
        'precision_bits': PRECISION_BITS,
        'command': stage.command,
        'exit_code': process.returncode,
        'success_marker': stage.marker,
        'elapsed_seconds': time.monotonic() - start,
        'binary_sha256': binary_hash,
        'log_sha256': sha256(log_path),
    }
    write_json(artifacts / (stage.name + '.json'), result)
    print('PASS', stage.name, round(result['elapsed_seconds'], 2), 's', flush=True)
    return result


def build():
    with (ARTIFACTS / 'build.log').open('w') as log:
        subprocess.run(['make', '-B', f'-j{WORKERS}', 'all'], cwd=ROOT, stdout=log,
                       stderr=subprocess.STDOUT, check=True)


def verify():
    """Compile, run every stage and write the certificate."""
    ARTIFACTS.mkdir(exist_ok=True)
    lock = (ARTIFACTS / '.run.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    start = time.monotonic()
    record = {
        'schema': SCHEMA,
        'status': 'RUNNING',
        'started_utc': utc_now(),
        'working_precision_bits': PRECISION_BITS,
        'human_review': 'pending',
        'scope': SCOPE,
    }
    write_json(CERTIFICATE, record)
    try:
        record['manuscript_inputs'] = check_inputs()
        record['source_sha256'] = source_hashes()
        build()
        compiler = subprocess.check_output(['gcc', '--version'], text=True)
        record['compiler'] = compiler.splitlines()[0]
        record['environment'] = {'platform': platform.platform(),
                                 'python_orchestrator': sys.version}
        libraries = library_hashes()
        headers = header_hashes()
        results = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for stage, result in zip(STAGES, pool.map(run_stage, STAGES)):
                results[stage.name] = result
        for stage in LATE_STAGES:
            results[stage.name] = run_stage(stage)
        # Nothing may have changed underneath the run.
        require(record['source_sha256'] == source_hashes(), 'Source/input change during full run')
        require(libraries == library_hashes() and headers == header_hashes(),
                'Arithmetic dependency change during run')
        check_inputs()
        record.update(
            status='PASS',
            finished_utc=utc_now(),
            elapsed_seconds=time.monotonic() - start,
            jobs=results,
            libraries_sha256=libraries,
            headers_sha256=headers,
            artifacts_sha256={'artifacts/moments.bin': sha256(ARTIFACTS / 'moments.bin')},
            arithmetic=ARITHMETIC,
            remaining=REMAINING,
        )
        write_json(CERTIFICATE, record)
        print(f'PASS complete all-C CAP at {PRECISION_BITS} bits', flush=True)
    except BaseException as error:
        record.update(status='UNRESOLVED', finished_utc=utc_now(), error=repr(error))
        write_json(CERTIFICATE, record)
        raise


def check_certificate(certificate):
    """Re-validate a certificate record against the current sources, binaries and logs."""
    require(certificate['status'] == 'PASS', 'No aggregate PASS')
    require(certificate['source_sha256'] == source_hashes(), 'Source/input mismatch')
    check_inputs()
    for name, result in certificate['jobs'].items():
        stage = STAGE_BY_NAME[name]
        require(result['command'] == stage.command and result['exit_code'] == stage.exit_code,
                'Changed command/contract')
        require(sha256(ROOT / stage.command[0]) == result['binary_sha256'],
                'Changed binary ' + name)
        log_path = ARTIFACTS / (name + '.log')
        require(sha256(log_path) == result['log_sha256'], 'Changed log ' + name)
        check_coverage(stage, log_path.read_text())
    require(set(certificate['jobs']) == set(STAGE_BY_NAME), 'Missing stages')
    for name, digest in certificate['artifacts_sha256'].items():
        require(sha256(ROOT / name) == digest, 'Changed artifact ' + name)
    for path, digest in (certificate['libraries_sha256'] | certificate['headers_sha256']).items():
        require(sha256(path) == digest, 'Changed arithmetic dependency ' + path)


def check_existing():
    check_certificate(json.loads(CERTIFICATE.read_text()))
    print('PASS existing all-C certificate: hashes, inputs and coverage contracts')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check-existing', action='store_true',
                        help='re-validate artifacts/ without compiling or running anything')
    if parser.parse_args().check_existing:
        check_existing()
    else:
        verify()


if __name__ == '__main__':
    main()
