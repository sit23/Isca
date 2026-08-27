#!/usr/bin/env python3
"""Select which trip_test test cases are relevant to the current changes,
then run trip_test_command_line on exactly that set.

Selection and execution happen in this one script/process by design: the
set of test cases that gets *run* can never drift out of sync with the set
that was *selected*, because there is only one place that decides and only
one place that acts on it.

Relevance is decided in three tiers:
  1. A "global trigger" file changed (the trip_test harness itself, the
     shared isca python framework, compiler/env config, this script, or the
     CI workflow) -> run everything, no cleverness.
  2. A specific test case's own definition (exp/test_cases/<dir>/*) changed
     -> that test case runs, regardless of any codebase-file overlap.
  3. Otherwise, a test case runs only if a changed src/ file appears in the
     path_names file of the codebase that test case is built with
     (src/extra/model/<codebase>/path_names), e.g. touching
     atmos_spectral_barotropic/ will never select shallow_water_stirring,
     and vice versa.

The test-case -> (directory, codebase) table is not hand-maintained here:
it's parsed straight out of trip_test_functions.py's get_nml_diag, which is
the one real definition of that mapping, so this script can't silently
drift from what trip_test itself actually does.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TRIP_TEST_FUNCTIONS = REPO_ROOT / 'exp' / 'test_cases' / 'trip_test' / 'trip_test_functions.py'
TRIP_TEST_COMMAND_LINE = REPO_ROOT / 'exp' / 'test_cases' / 'trip_test' / 'trip_test_command_line'
MODEL_DIR = REPO_ROOT / 'src' / 'extra' / 'model'

# CodeBase class name (as used in trip_test_functions.py) -> CodeBase.name
# (as used for the src/extra/model/<name>/path_names file). See
# src/extra/python/isca/codebase.py.
CODEBASE_NAME = {
    'IscaCodeBase': 'isca',
    'DryCodeBase': 'dry',
    'GreyCodeBase': 'grey',
    'ColumnCodeBase': 'column',
    'SocratesCodeBase': 'socrates',
    'ShallowCodeBase': 'shallow',
    'BarotropicCodeBase': 'barotropic',
}

# Changing any of these means the selection logic itself can't be trusted
# to reason about the change -- fall back to running everything.
GLOBAL_TRIGGER_PREFIXES = [
    'exp/test_cases/trip_test/',
    'src/extra/python/isca/',
    'src/extra/env/',
    'ci/',
    '.github/workflows/',
]


def parse_test_case_table():
    """test_case_name -> (test_case_dir, codebase_name), parsed out of
    trip_test_functions.py's get_nml_diag rather than duplicated here."""
    text = TRIP_TEST_FUNCTIONS.read_text()
    blocks = re.split(r"\n(?=    if '.*?' in test_case_name:)", text)
    table = {}
    for block in blocks:
        name_match = re.match(r"    if '(.*?)' in test_case_name:", block)
        if not name_match:
            continue
        name = name_match.group(1)
        dir_match = re.search(r"exp/test_cases/([^/']+)/", block)
        codebase_match = re.search(r"codebase_to_use\s*=\s*(\w+)", block)
        if not (dir_match and codebase_match):
            continue
        table[name] = (dir_match.group(1), CODEBASE_NAME.get(codebase_match.group(1)))
    return table


def default_test_case_list():
    """The 'all' list, parsed from list_all_test_cases_implemented_in_trip_test
    so this script only ever offers test cases trip_test itself knows about.
    Entries commented out of exps_implemented (opt-in-only test cases, e.g.
    ones needing data files not in the repo) are deliberately excluded."""
    text = TRIP_TEST_FUNCTIONS.read_text()
    m = re.search(r"exps_implemented = \[(.*?)\]\n\n\s*return", text, re.S)
    active_lines = [
        line for line in m.group(1).splitlines()
        if not line.strip().startswith('#')
    ]
    return re.findall(r"'(\w+)'", '\n'.join(active_lines))


def load_path_names(codebase_name):
    path_names_file = MODEL_DIR / codebase_name / 'path_names'
    if not path_names_file.exists():
        return set()
    return {line.strip() for line in path_names_file.read_text().splitlines() if line.strip()}


def changed_files(base_sha, head_ref):
    result = subprocess.run(
        ['git', 'diff', '--name-only', f'{base_sha}...{head_ref}'],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def select_test_cases(base_sha, head_ref):
    files = changed_files(base_sha, head_ref)
    all_cases = default_test_case_list()

    if any(f.startswith(prefix) for f in files for prefix in GLOBAL_TRIGGER_PREFIXES):
        return all_cases, files, 'a global-trigger file changed'

    table = parse_test_case_table()
    exp_changed = [f for f in files if f.startswith('exp/test_cases/')]
    src_changed = {f[len('src/'):] for f in files if f.startswith('src/')}

    path_names_cache = {}
    selected = set()

    for name in all_cases:
        if name not in table:
            continue
        test_dir, codebase_name = table[name]

        if any(f.startswith(f'exp/test_cases/{test_dir}/') for f in exp_changed):
            selected.add(name)
            continue

        if codebase_name is None:
            continue
        if codebase_name not in path_names_cache:
            path_names_cache[codebase_name] = load_path_names(codebase_name)
        if src_changed & path_names_cache[codebase_name]:
            selected.add(name)

    return sorted(selected), files, None


def drop_unavailable_codebases(test_cases, excluded_codebases):
    """Remove test cases built with a codebase this environment can't
    actually build. Currently used for 'socrates': SocratesCodeBase needs
    externally-licensed source at GFDL_SOC that public CI doesn't have, and
    trip_test_functions.py doesn't catch that failure per-test-case -- it's
    an uncaught OSError that kills the whole comparison run, so anything
    still queued behind a Socrates test case never gets a chance to run.
    See the socrates_opensource_auto branch (auto-fetches Socrates instead
    of requiring a pre-supplied GFDL_SOC) -- once that lands, CI may be able
    to build Socrates itself and this exclusion can likely be dropped."""
    if not excluded_codebases:
        return test_cases, []
    table = parse_test_case_table()
    excluded = set(excluded_codebases)
    kept, dropped = [], []
    for name in test_cases:
        codebase_name = table.get(name, (None, None))[1]
        (dropped if codebase_name in excluded else kept).append(name)
    return kept, dropped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base_ref', help="e.g. 'origin/master' -- resolved to a SHA before use")
    parser.add_argument('head_ref', help='commit to compare against base_ref, e.g. github.sha')
    parser.add_argument('--force-full', action='store_true',
                         help='escape hatch: skip selection, run every test case trip_test knows about')
    parser.add_argument('-n', '--num-cores', type=int, default=2)
    parser.add_argument('-r', '--repo', required=True, help='repo URL for trip_test_command_line -r')
    parser.add_argument('--exclude-codebases', nargs='*', default=[],
                         help="codebase name(s) (e.g. 'socrates') this environment cannot build, "
                              'so any test case using them is dropped regardless of selection')
    parser.add_argument('--dry-run', action='store_true',
                         help='print the selection and exit without invoking trip_test')
    args = parser.parse_args()

    base_sha = subprocess.run(
        ['git', 'rev-parse', args.base_ref],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    ).stdout.strip()

    if args.force_full:
        selected, files, reason = default_test_case_list(), None, 'forced via escape hatch'
    else:
        selected, files, reason = select_test_cases(base_sha, args.head_ref)

    selected, dropped = drop_unavailable_codebases(selected, args.exclude_codebases)

    print('=' * 70)
    print(f'trip_test selection: base={args.base_ref} ({base_sha}) head={args.head_ref}')
    if reason:
        print(f'running the FULL suite: {reason}')
    else:
        print(f'{len(files)} changed file(s) under consideration')
    if dropped:
        print(f'excluding {len(dropped)} test case(s) needing unavailable codebase(s) '
              f'{args.exclude_codebases}: {dropped}')
    print(f'selected {len(selected)} test case(s): {selected}')
    print('=' * 70)

    if not selected:
        print('No relevant test cases selected -- nothing to run.')
        return

    if args.dry_run:
        return

    cmd = [
        sys.executable, str(TRIP_TEST_COMMAND_LINE),
        base_sha, args.head_ref,
        '-e', *selected,
        '-n', str(args.num_cores),
        '-r', args.repo,
    ]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


if __name__ == '__main__':
    main()
