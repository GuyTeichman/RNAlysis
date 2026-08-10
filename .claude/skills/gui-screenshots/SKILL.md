---
name: gui-screenshots
description: >-
  Capture and attach RNAlysis GUI screenshots to a PR whenever a change makes a VISIBLE
  difference to a GUI dialog. Use when adding/renaming/removing a parameter of a public
  Filter/FeatureSet/fastq/enrichment function, changing a parameter's type annotation, adding or
  editing @readable_name, or editing :param: help text — before opening or updating the PR. Skip
  for behind-the-scenes changes with no visible effect on any dialog. Also use when asked to
  "screenshot a dialog", "grab the GUI for <function>", or "show what the dialog looks like".
---

# GUI screenshots for PRs

RNAlysis' GUI is generated from the API by reflection (see `CLAUDE.md`). That cuts both ways:
a signature or annotation edit on a public function **silently reshapes a GUI dialog**, and — the
useful direction — that same dialog can be **rebuilt and grabbed to a PNG straight from the API**,
headlessly, without clicking through the running app. This skill captures those PNGs and attaches
them to the PR so a reviewer can see the visible change.

The engine is a plain, provider-neutral script: **`packaging/capture_gui_dialog.py`**. This skill
adds the *when*, the collision-safe conventions, and the assets-branch/PR-attach workflow.

---

## Step 1 — Is this change GUI-visible? (skip the skill if not)

Attach screenshots only when a **user would see a difference in a dialog**. It is GUI-visible if
the change is to a **public** (non-`_`) function on `Filter`/`CountFilter`/`DESeqFilter`/
`FoldChangeFilter`/`FeatureSet`/`RankedSet` or in `fastq`/`enrichment`, **and** it does any of:

- adds, renames, removes, or reorders a **parameter**;
- changes a parameter's **type annotation** (a different `param_typing` type → a different widget);
- adds/changes **`@readable_name('…')`** (the dialog/button title);
- edits a `:param x:` **docstring line** (per-parameter help/tooltip text);
- changes a **default** that the widget shows (e.g. a combo's preselected value).

**Skip** (no screenshot needed) for: internal helpers (leading `_`, or in a widget's
`EXCLUDED_FUNCS`/`EXCLUDED_PARAMS`), pure implementation/refactor with an unchanged signature,
performance work, docstring prose outside `:param:` lines, tests, CI, and non-GUI modules
(`utils/`, `io.py`, R templates, etc.). When unsure, capture it — a redundant screenshot is
cheaper than a silent GUI regression slipping through.

For a **renamed/retyped** parameter, capture **before and after** so the diff is legible: grab
once from `development` (or the merge-base) and once from your branch.

---

## Step 2 — Capture

Run in the project's dev environment (the one with **PyQt6** installed — locally that's the
`rnalysis` conda env). Pass one or more dotted targets (`<module>.<Class>.<method>` or
`<module>.<function>`), resolved from the `rnalysis` package:

```bash
python packaging/capture_gui_dialog.py \
    filtering.CountFilter.filter_low_reads \
    enrichment.FeatureSet.go_enrichment \
    --out .gui_shots
```

One PNG per target lands in `--out`, named `<dotted.target>.png`. The script:

- uses the **native** window platform so text renders with real fonts (Qt's bare `offscreen`
  platform renders text as missing-glyph boxes — only pass `--offscreen` for a quick layout check);
- hides the same infra parameters the real app hides (`self`, `backend`, `gui_mode`, …), so the
  shot matches what a user sees. Use `--show-all` to include them, or `--exclude a,b` to hide more;
- grabs the **whole form** (every parameter, no scrollbar), which is what makes it documentation-grade.

**Headless CI note:** on a runner with no display, wrap the command in `xvfb-run` on Linux so a
real framebuffer resolves fonts. macOS/Windows CI runners have a window session and work directly.

Sanity-check the PNGs before publishing (open them / view them). A dialog full of tofu boxes means
you're on the `offscreen` platform without a framebuffer — fix that, don't publish it.

---

## Step 3 — Publish to the `assets` branch (collision-safe)

Screenshots live on a dedicated **`assets`** branch, not in the PR's code diff — keeps binaries out
of mainline history, and any branch's blobs are linkable. **Namespace every PR's images under
`pr-<N>/`** so parallel agents never clobber each other's files or race on the same paths.

Do the publish in a **throwaway git worktree** so it never disturbs your feature-branch worktree
(or Guy's working tree — see the repo's worktree convention). From the repo root:

```bash
N=<pr-number-or-branch-slug>
git fetch origin
# check out the assets branch in a scratch worktree (create it on first use)
git worktree add ../rnalysis-assets assets 2>/dev/null \
  || git worktree add -b assets ../rnalysis-assets origin/assets 2>/dev/null \
  || git worktree add -b assets ../rnalysis-assets    # first time ever: fresh branch

mkdir -p "../rnalysis-assets/pr-$N"
cp .gui_shots/*.png "../rnalysis-assets/pr-$N/"

git -C ../rnalysis-assets add "pr-$N"
git -C ../rnalysis-assets commit -m "GUI screenshots for PR #$N"
git -C ../rnalysis-assets pull --rebase origin assets 2>/dev/null || true   # absorb concurrent pushes
git -C ../rnalysis-assets push origin assets

git worktree remove ../rnalysis-assets
```

If the push is rejected (a parallel agent pushed first), re-run the `pull --rebase` + `push`. Because
each PR writes only under its own `pr-$N/`, rebases apply cleanly — the collision is a push race, not
a content conflict.

> First-time-ever bootstrap: if no `assets` branch exists yet, an **orphan** branch keeps it out of
> code history — `git worktree add --detach ../rnalysis-assets && git -C ../rnalysis-assets switch
> --orphan assets && …` — but the plain fresh-branch fallback above is fine too. Ask Guy which he
> prefers before creating it.

---

## Step 4 — Attach to the PR

Embed the raw URLs in the PR body or a comment. Get `owner/repo` from `gh` rather than hardcoding:

```bash
REPO=$(gh repo view --json nameWithOwner -q .nameWithOwner)
# per image: https://raw.githubusercontent.com/$REPO/assets/pr-$N/<file>.png
```

Write the markdown to a file and post with `--body-file` (never inline a body containing `—`/`→` on
Windows — a Python/text round-trip corrupts UTF-8; use a file):

```markdown
### GUI changes

**`filter_low_reads` dialog**

![filter_low_reads dialog](https://raw.githubusercontent.com/OWNER/RNAlysis/assets/pr-123/filtering.CountFilter.filter_low_reads.png)
```

```bash
gh pr comment "$N" --body-file gui_shots_comment.md
```

For a before/after, put the two images side by side (an HTML `<table>` with one `<img>` per cell
renders in GitHub markdown) and label which is `development` and which is the branch.

---

## Fidelity & limits (be honest in the PR)

- The capture is **faithful for the parameter form** — widgets, labels, defaults, help buttons —
  because that form *is* the reflection output. It does **not** capture the OS window title bar
  (Qt `.grab()` never does), nor live end-to-end states (loaded data, result plots, hover tooltips,
  open menus, multi-window flows). Those "display-level" tutorial shots are a human task; don't fake
  them with this tool.
- Rendering varies by **OS theme (light/dark), DPI, and font availability**, so an agent's PNG may
  not pixel-match a given user's screen. That's fine for "here's the visible change" — do **not**
  repurpose these PNGs as a visual-regression baseline.
- The `ToggleSwitch` may clip a long label ("False"→"alse"); that's the fixed-size custom widget
  itself, same as in the live app — not a capture artifact.

---

## Definition of done

A PR that makes a GUI-visible change is not finished until the relevant dialog screenshot(s) are on
the `assets` branch under `pr-<N>/` and linked from the PR. Note in the PR body anything you could
**not** render (e.g. a display-level flow) and why.
