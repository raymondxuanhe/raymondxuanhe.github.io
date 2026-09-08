# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Ground rules (non-negotiable)

1. **Never delete data or code** — do not remove, overwrite, or destructively change any
   file, directory, or content unless the user gives express permission for that specific
   deletion.
2. **Stay inside this directory** — do not read, write, or operate on anything outside this
   repository directory unless the user gives express permission.

## What this is

Personal academic website for Raymond He (economics, UT Austin), hosted at
https://raymondxuanhe.github.io via GitHub Pages. It is a customized fork of the
**academicpages** Jekyll template (itself derived from Minimal Mistakes). Most of the
directory structure and template plumbing is generic academicpages; the sections below
focus on what is specific to this site or otherwise non-obvious.

## Commands

```bash
bundle install                              # install Ruby/Jekyll deps (delete Gemfile.lock if it errors)
bundle exec jekyll serve -l -H localhost    # serve + live-reload at localhost:4000
docker compose up                           # alternative: serve via the provided Dockerfile
npm run build:js                            # rebuild assets/js/main.min.js after editing assets/js/**
```

Note: `_config.yml` is **not** reloaded by `jekyll serve` — restart the server after editing it.

There is no test suite. GitHub Pages builds and deploys automatically on push to `master`.

## Content model

Content lives in front-matter markdown/HTML, not code. Pages are in `_pages/`; the four
Jekyll **collections** (configured in `_config.yml`) are `_publications/`, `_talks/`,
`_teaching/`, and `_portfolio/`. Site-wide config and the top navigation are in
`_config.yml` and `_data/navigation.yml`. Layouts (`_layouts/`), partials (`_includes/`),
and styles (`_sass/`) are the generic template internals — change these only for
structural/theming work.

Publications and talks can be bulk-generated from TSV via the scripts in
`markdown_generator/` (`publications.py`/`talks.py` read `*.tsv` and emit collection
markdown). This is optional; entries are also written by hand.

## Files and downloads

`files/` holds all downloadable assets served at `/files/...`, organized by purpose
(`CV/`, `research/`, `teaching/`, `benchmarks/`, `julia_files/`, `Bio/`). The CV is
authored in **LaTeX** at `files/CV/raymond_he_CV.tex` and committed alongside its compiled
`.pdf` (plus latexmk build artifacts). The job-market paper PDF lives under
`files/research/jmp/`.

`scripts/update_cv_json.sh` + `scripts/cv_markdown_to_json.py` convert a markdown CV into
`_data/cv.json`. Heads up: the shell script expects `_pages/cv.md`, which does not exist in
this repo — this pipeline is inherited from the template and is not the live CV workflow
(the CV is the LaTeX file above).

## Julia / computational benchmarks

Unusually for a Jekyll site, this repo also carries a Julia project (`Project.toml`,
`Manifest.toml`) used to produce the economics computation benchmarks and charts shown on
the site. The Julia sources live in `files/julia_files/` (e.g. `opt_growth.jl`,
`generate_benchmark_charts.jl`); results are surfaced via `_pages/computational_benchmarks.*`.
This is separate from the website build — you do not need Julia to serve or deploy the site.
