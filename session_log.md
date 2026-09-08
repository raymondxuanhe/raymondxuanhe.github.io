# Session Log

A running log of what occurred in each Claude Code session. Newest entries at the top.

Each entry records: date, a summary of what was done, and files touched.

---

## 2026-09-08 — Research page redesign, image cleanup, site font fix

**Summary**
- **Deleted 7 unused template demo images** from `images/` (verified unreferenced across
  the whole repo). Kept the 3 in use: `bio-photo.jpg`, `bio-photo-2.jpg`
  (`_data/authors.yml`), `raymondhe1.jpg` (`_config.yml`). Note: images live in `images/`,
  not `files/images/`.
- **Redesigned the Research page** (`_pages/research.md`, the live page — `research.html`
  is empty template leftover since no `research` collection is defined) to mimic
  Ryohei Oishi's site (ryohei-oishi.github.io). Each Working Paper is now a card: bold
  title, gray coauthor/status line, an always-visible control row with a "▶ Abstract"
  toggle plus outlined link buttons (red `Paper`, blue others). **Hovering a paper**
  expands its abstract (which sits above the controls), rotates the triangle, and
  highlights the card. Works-in-Progress entries are plain titles (nothing to expand).
  Abstracts are justified. Fully dark-mode aware via local CSS variables overridden under
  `html[data-theme="dark"]` (fixes the previously washed-out "Works in Progress" text in
  night mode). Titles are plain text (not links); the red `Paper` button is the PDF link.
- **Fixed the site-wide font** (`_sass/_themes.scss`): changed the `$sans-serif` stack from
  the stale Safari-only `".SFNSText-Regular"`/`"San Francisco"` names to
  `-apple-system, BlinkMacSystemFont, "Roboto", …` — matching the Minimal Mistakes default
  that Oishi's site uses. `BlinkMacSystemFont` makes Chrome render the San Francisco system
  font too (previously Chrome fell through to Helvetica Neue). Applies to the whole site
  since everything inherits `$sans-serif`.

**Local preview environment (workarounds, no committed project files changed for these)**
- System Ruby was 2.6 (too old). Installed Ruby via Homebrew; settled on **`ruby@3.3`**
  (keg-only) for this old GitHub-Pages-pinned Jekyll 3.9 stack.
- Served with `--safe` (skips a bogus `jekyll-theme` plugin entry in `_config.yml` that
  GitHub Pages also ignores) plus a `RUBYOPT` shim in the scratchpad that re-adds the
  no-op `tainted?`/`taint`/`untaint` methods Ruby 3.2+ removed but Liquid 4.0.3 still calls.
- Run command:
  `PATH="/opt/homebrew/opt/ruby@3.3/bin:$PATH" RUBYOPT="-r<scratchpad>/taint_shim.rb" bundle exec jekyll serve --safe -l -H localhost`

**Files touched**
- `images/` (deleted 7 unused `.jpg` files)
- `_pages/research.md` (rewritten as interactive hover-expand cards)
- `_sass/_themes.scss` (updated `$sans-serif` font stack)

---

## 2026-09-08

**Summary**
- Created `README.md` documenting the repository's directory structure (replaced the
  upstream academicpages template README).
- Added a **Ground rules (non-negotiable)** section to `CLAUDE.md`: (1) never delete or
  destructively change data/code without express permission, (2) never operate outside
  this repository directory without express permission.
- Created this `session_log.md`.

**Files touched**
- `README.md` (rewritten)
- `CLAUDE.md` (added ground rules section)
- `session_log.md` (new)
