# Session Log

A running log of what occurred in each Claude Code session. Newest entries at the top.

Each entry records: date, a summary of what was done, and files touched.

---

## 2026-09-08 — Homepage redesign (centered nav + two-column hero)

**Summary**
- **Centered the top navigation** site-wide: added a small `<style>` to
  `_includes/head/custom.html` (`.greedy-nav .visible-links { margin: 0 auto }`). The nav
  is a `display:table`, so auto side-margins center the whole group (site title + links +
  theme toggle). Global (all pages), theme-aware.
- **Created a dedicated homepage layout** `_layouts/home.html` and pointed `about.md` at it
  (`layout: home`). Iterated on the design in three steps:
  1. First a *centered single-column* hero (profile photo + name + links centered on top),
     scoped to the homepage only (other pages keep the left-sidebar layout).
  2. Then, per request, switched to a **two-column hero**: rectangular photo (rounded
     corners, ~250px, natural portrait aspect) on the **left**; on the **right** the name,
     a subtitle, field pills, the biography, and contact links. Responsive: stacks and
     centers on screens ≤600px.
- **Right-column content** (driven by `about.md` front matter so it stays editable):
  - `subtitle: "Ph.D. in Economics Candidate, University of Texas at Austin"`
  - `fields:` list rendered as rounded outline **pills** — International Macroeconomics,
    Labor Economics, Sovereign Debt.
  - Biography, then contact links.
- **Contact links**: reuse the `author-profile.html` include inside the hero with its
  avatar/name/Follow-button hidden via CSS, showing only the links as a left-aligned row.
- **Removed `location` and `linkedin`** from the `author:` block in `_config.yml`
  (commented out, not deleted) — removes "Austin, Texas" and "LinkedIn" everywhere (homepage
  hero + sidebars on other pages). NOTE: `_config.yml` is not hot-reloaded by `jekyll serve`,
  so the server must be restarted after this change.
- **Rewrote the homepage bio** in `about.md` as a mock, research-focused blurb (placeholder
  for Raymond to edit) to remove repetition with the new subtitle + field pills; kept the
  "on the academic job market in Fall 2026" line.
- All homepage CSS is scoped inside `_layouts/home.html` and uses theme variables, so dark
  mode works and no other page is affected. The top-nav centering is the only global change.

**Files touched**
- `_includes/head/custom.html` (global nav-centering `<style>`)
- `_layouts/home.html` (new homepage layout: two-column hero + subtitle + field pills)
- `_pages/about.md` (`layout: home`, `subtitle`/`fields` front matter, new research blurb)
- `_config.yml` (commented out `author.location` and `author.linkedin`)

**Open items / next time**
- Everything is uncommitted (working tree only), alongside earlier uncommitted work.
- Raymond is going to sit with the homepage design and revisit later — the bio blurb, field
  list, subtitle, and photo crop are all placeholders/easy to tweak.
- Possible follow-ups discussed but not done: forcing a portrait crop on the photo; centering
  the masthead nav on the homepage only (currently global); enlarging the name.

---

## 2026-09-08 — Research page redesign, image cleanup, site font fix, unused-page purge

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

- **Purged unused pages and template cruft** (all verified unreferenced across nav,
  permalinks, links, config, and a clean `jekyll build --safe` afterward):
  - Deleted 10 unused `.html` pages from `_pages/`: `research.html` and
    `computational_benchmarks.html` (empty duplicates of the live `.md` pages — removing
    them also cleared the permalink conflicts), plus the template demo pages
    `collection-archive`, `page-archive`, `year-archive`, `talkmap`, `category-archive`,
    `tag-archive`, `portfolio`, `publications`.
  - Deleted 4 unused demo `.md` pages from `_pages/`: `markdown.md`,
    `archive-layout-with-content.md`, `non-menu-page.md`, `terms.md`.
  - Deleted the leftover academicpages **sample** collections: all 5 `_publications/`
    entries ("Paper Title Number 1"…), both `_portfolio/` entries, and the now-empty
    `_publications/` and `_portfolio/` directories.
  - `_config.yml`: removed `publications:`/`portfolio:` from `collections:`, removed their
    `defaults:` scopes, and removed the `category_archive:`/`tag_archive:` settings (site
    uses no tags/categories taxonomy and those archive pages are gone).
  - Remaining `_pages/` (all real): `404.md`, `about.md`, `research.md`,
    `computational_benchmarks.md`, `sitemap.md`, `teaching.html` (live Teaching page, in
    nav). `teaching.html` was deliberately kept — it is the display page for the real
    `_teaching` collection.

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
- `_pages/` (deleted 10 unused `.html` + 4 demo `.md` pages)
- `_publications/`, `_portfolio/` (deleted all sample entries + the empty dirs)
- `_config.yml` (removed publications/portfolio collections + defaults, and
  category/tag archive settings)

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
