---
name: external-api-change
description: >-
  Workflow for fixing or changing RNAlysis code that talks to an EXTERNAL WEB SERVICE — UniProt,
  Ensembl, PANTHER, PhylomeDB, OrthoInspector, KEGG, or GO. Use before touching
  rnalysis/utils/io.py or its callers (enrichment.py, enrichment_runner.py, ontology.py) to fix a
  broken/changed API, handle a format change, or add a new call to one of these services. Also use
  when asked to "capture a fixture for an external API", "the GO/UniProt/Ensembl/PANTHER/KEGG/
  PhylomeDB/OrthoInspector integration broke", or "mock an external service in a test". Skip for
  CLI-tool integrations (kallisto, bowtie2, R) — those are a different fragile area, not this one.
---

# Changing code that talks to an external web service

This is, by the project's own account, **the most fragile part of the codebase** (see
`CLAUDE.md` gotchas and `.claude/context.md` fragile spots #1). UniProt, Ensembl, PANTHER,
PhylomeDB, OrthoInspector, KEGG, and GO change response formats, rename fields, and go down
without notice — all outside RNAlysis' control. The relevant code lives in
**`rnalysis/utils/io.py`** (async `aiohttp` + `tenacity` retries + `aiolimiter` rate limiting +
response caches) and its callers: **`enrichment.py`**, **`utils/enrichment_runner.py`**,
**`utils/ontology.py`**.

Real precedents worth knowing before you start: OrthoInspector relocated its whole API to
`api.bigest-icube.fr`; Ensembl silently dropped cross-division ortholog data (worm→human now
returns nothing); PANTHER can return an empty HTTP 200 body instead of an error. None of these
were visible from reading old code — each needed a live check.

This is TDD like everything else in the repo, with one twist: **step 1 happens before you write
any test**, because you cannot write an honest test against a service whose current behavior you
haven't confirmed.

---

## When to use this vs. skip it

Use this skill for any change to code that issues a request to one of the six services above, or
adds a new one. Skip it for CLI-tool bridges (kallisto/bowtie2/cutadapt/R) — those have their own
fragility (see `CLAUDE.md`'s R/permissions gotcha and `.claude/workflows.md`'s R-bridge workflow)
but aren't a web-API problem.

---

## Step 1 — Confirm the service's REAL current behavior before writing any code

Do not trust the existing code, your training data, or memory of "how this API works" — these
services change often enough that stale assumptions are the #1 way this class of fix goes wrong.
Before touching `io.py` or its callers:

- Hit the endpoint live (curl / browser / a throwaway Python snippet) and read the actual
  response — status code, headers, body shape, field names.
- Or use the `research` skill to investigate and record what you found.
- If network access isn't available in your environment, say so explicitly and get a human (or
  an agent with network access) to confirm before coding against assumptions.

Write down what you found (in the PR description at minimum) — the next person debugging this
service needs that evidence trail, not just the resulting diff.

## Step 2 — Preserve the `io.py` scaffolding; degrade gracefully

Keep the existing resilience layers intact while you fix the specific breakage:

- **`tenacity`** retry decorators around the request.
- **`aiolimiter`** rate limiting (services like UniProt idmapping throttle hard; CI already hits
  this — see the CI note below).
- The **response caches** (don't bypass them while debugging and then forget to restore them).

A changed or dead service must **never crash the app**. If the new response shape can't be
parsed, fail toward an empty/partial result with a clear, user-facing message — not an unhandled
exception. This is a non-programmer-facing GUI (per `CLAUDE.md`); a stack trace is not an
acceptable failure mode for "PANTHER is down today."

## Step 3 — Capture a real payload as a test fixture, then mock it

Tests must never depend on the live service (**flaky by construction**, and CI's network tier is
narrow — see below). Capture one:

```bash
python packaging/snapshot_api_payload.py https://rest.uniprot.org/idmapping/status/abc123
python packaging/snapshot_api_payload.py https://api.example.org/search \
    --param query=BRCA1 --header "Accept: application/json" \
    --out tests/test_files/example_search_response.json
```

This provider-neutral helper does a plain GET and writes the raw response body byte-for-byte to
`tests/test_files/` (default) or `--out`, deriving a filename from the URL when you don't name
one. Run `python packaging/snapshot_api_payload.py --help` for the full option list. **It makes a
real live request** — run it locally/by hand where you have network, not in a sandboxed or
offline environment. It does not itself write or run any test; it only produces the fixture file.

Then write the test in the matching `tests/test_*.py` (`test_io.py` for most of these) mocking
the HTTP layer against that fixture. The repo's established convention is `requests_mock` for
synchronous `requests` calls (see the many `with requests_mock.Mocker() as m:` blocks already in
`tests/test_io.py`) — match whichever HTTP client (`requests` vs `aiohttp`) the code under test
actually uses, and mock at that layer.

Red → green → refactor from here as usual: write the failing test against the captured fixture,
make it pass, clean up, run the module's test file.

## Step 4 — Know your CI network tier

Read `tests/conftest.py` and `.github/workflows/build_ci.yml` if you need the exact current
rules; as of this writing:

- `tests/conftest.py` auto-assigns `test_io.py` to the **`integration_net`** tier (network/web-API
  tests), a leaf tier under the `integration` umbrella marker. A test elsewhere gets the same tier
  automatically if it carries a `skipif` whose reason mentions "available" (the live-network tests
  scattered through `test_enrichment.py` etc.).
- `build_ci.yml` runs the `integration_net` tier on **exactly one** matrix cell — `ubuntu-latest`
  / **Python 3.13** — out of the full 3 OS × 3 Python matrix; the other 8 cells skip that step
  outright (the comment in the workflow explains why: avoiding 9x redundant load on rate-limited
  services). **A network-tier red is a 1-of-1 failure, not a 1-of-9 flake** — don't wave it off
  as "just one matrix cell."
- Run it locally with `pytest -m integration_net tests/test_io.py` if you have network access;
  otherwise say plainly in the PR that you could not run it (per `CLAUDE.md` rule 1).

---

## Definition of done

- You confirmed the service's current real behavior (Step 1) and said how (live hit, or the
  `research` skill), not just re-read the old code.
- `io.py`'s retry/rate-limit/cache scaffolding is intact; a broken/changed/down service degrades
  to an empty or partial result with a clear message, never an uncaught crash.
- A real captured payload lives under `tests/test_files/` (via `snapshot_api_payload.py` or
  equivalent) and the test mocks against it — no test in the suite depends on the live service to
  pass.
- Relevant tests pass locally; you ran `integration_net` if you had network, and said so plainly
  if you didn't.
- If this also changes a public API signature or the GUI reflection contract, that's a *second*,
  separate risk — see `CLAUDE.md` rule 3 (plan-first) before combining the two.
